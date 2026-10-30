#include <math.h>
#include <string.h>
#include "Florence_f.h"
#include "Florence_c.h"

/* One Florence evaluation (one Fortran run) for rotational time TR, split as
 *   r[0] = P/T1M          (paramagnetic rate scaled by P, no exchange)
 *   r[1..3] = fraction of the inner sphere that is dipolar, contact, cross
 *   r[4] = outer sphere
 * RUNO(i) = P*R_i/(1+R1M*TAUM), so R_i/R1M = RUNO(i)/inner, and P/T1M
 * follows from inner = P/(T1M+TAUM). IREL = 1 (T1), S4M = 0, ACQ = 1. */
static void florence_split(double r[5], double f, double P, double SI, double GAMMAI, double SPIN, double DELTA2, double TR, double TAUVM, double TAUMM, double TAUDELTA, double DPARAM, double EPARAM, double GXM, double GYM, double GZM, double AXM, double AYM, double AZM, double DM, double DDM, double CONCM, double AMOLFRAM, double RKM, double ACONTM, double THETAM, double PHIM, double FLAG)
{
  const double parameters[28] = { SI, GAMMAI, SPIN, 1.0, DELTA2, TR, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, 0.0, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, 1.0, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG };
  double c[5], inner;
  int i;
  FlorenceRUNO(f, parameters, c);
  /* c[0] is the exact total; c[4] the outer sphere */
  inner = c[0] - c[4];
  r[4] = c[4];
  if (inner == 0.0) {
    r[0] = r[1] = r[2] = r[3] = 0.0;
    return;
  }
  r[0] = P/(P/inner - TAUMM);
  for (i = 1; i <= 3; i++)
    r[i] = c[i]/inner;
}

/* Lipari-Szabo Florence: R1M = S^2 R1M(TAURM) + (1-S^2) R1M(TAURM*TAULM/(TAURM+TAULM)),
 * with SLS = S, the order parameter (not S^2),
 * then exchange P/(T1M+TAUMM), then the outer sphere.
 * index 1,2,3 = dipolar, contact, cross share of the exchange-limited inner sphere,
 * index 4 = outer sphere, anything else = total. Indices 1..4 add up to the total,
 * so the model can be written as FlorenceLS(1,...)\+FlorenceLS(2,...)\+...
 * The engine (funcsn) calls the \+ pieces one after another for the same point,
 * so all five results are computed once and kept: the other pieces cost no
 * Florence evaluation. Two Fortran runs per point (TAURM and the fast time),
 * one when SLS^2 is 0 or 1. Calls must stay serialized, as for Florence_c.c. */
double FlorenceLS(double index, double f, double SI, double GAMMAI, double SPIN, double DELTA2, double TAURM, double TAUVM, double TAUMM, double TAUDELTA, double DPARAM, double EPARAM, double GXM, double GYM, double GZM, double AXM, double AYM, double AZM, double CONCM, double RKM, double DM, double DDM, double ACONTM, double AMOLFRAM, double THETAM, double PHIM, double SLS, double TAULM, double FLAG)
{
  static int valid = 0;
  static double last_in[28], last_out[5];
  const double in[28] = { f, SI, GAMMAI, SPIN, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, GXM, GYM, GZM, AXM, AYM, AZM, CONCM, RKM, DM, DDM, ACONTM, AMOLFRAM, THETAM, PHIM, SLS, TAULM, FLAG };
  double P = AMOLFRAM*1e-6/(111.0*CONCM), S2 = SLS*SLS, s[5], q[5], R, x;
  int i, k = index >= 1.0 && index < 5.0 ? (int)index : 0;

  if (valid && memcmp(in, last_in, sizeof(in)) == 0)
    return last_out[k];

  if (S2 != 0.0)
    florence_split(s, f, P, SI, GAMMAI, SPIN, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG);
  if (S2 != 1.0)
    florence_split(q, f, P, SI, GAMMAI, SPIN, DELTA2, TAURM*TAULM/(TAURM+TAULM), TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG);
  if (S2 == 0.0)
    memcpy(s, q, sizeof(s));
  else if (S2 == 1.0)
    memcpy(q, s, sizeof(q));

  /* P*R1M of the Lipari-Szabo mixture; the exchange factor is shared by all mechanisms */
  R = S2*s[0] + (1.0-S2)*q[0];
  x = 1.0/(1.0 + R*TAUMM/P);
  for (i = 1; i <= 3; i++)
    last_out[i] = (S2*s[0]*s[i] + (1.0-S2)*q[0]*q[i])*x;
  last_out[4] = s[4];   /* outer sphere does not depend on TAURM */
  last_out[0] = R*x + s[4];

  memcpy(last_in, in, sizeof(in));
  valid = 1;
  for (i = 0; i <= 4; i++)
    if (!isfinite(last_out[i])) valid = 0;
  return last_out[k];
}
