#include <stdio.h>
#include <math.h>
#include "Florence_c.h"

/* FlorenceLS: the \+ pieces 1..4 add up to the total, the total equals the
 * Lipari-Szabo combination built independently from FlorenceN, and SLS = 1
 * gives plain Florence. Slow and fast exchange, outer sphere on and off. */
static double SI = 0.5, GAMMAI = 2.675222e8, SPIN = 3.5, DELTA2 = 0.0236;
static double TAURM = 2.88e-11, TAUVM = 2.86e-12, TAUDELTA = 2.0, G = 2.003;
static double CONCM = 0.001, AMOLFRAM = 12, RKM = 3.1, ACONTM = 0.748, TAULM = 5e-12;

static double N(int i, double f, double TR, double TM, double DM, double DDM)
{
    return FlorenceN(i, 4, f, SI, GAMMAI, SPIN, 1, DELTA2, TR, TAUVM, TM, TAUDELTA,
                     0, 0, 0, G, G, G, 0, 0, 0, DM, DDM, CONCM, 1, AMOLFRAM, RKM,
                     ACONTM, 0, 0, 1);
}

static double LS(int i, double f, double TM, double DM, double DDM, double SLS)
{
    return FlorenceLS(i, f, SI, GAMMAI, SPIN, DELTA2, TAURM, TAUVM, TM, TAUDELTA,
                      0, 0, G, G, G, 0, 0, 0, CONCM, RKM, DM, DDM, ACONTM,
                      AMOLFRAM, 0, 0, SLS, TAULM, 1);
}

int main(void)
{
    double cases[][4] = { /* TAUMM, DM, DDM, SLS^2 */
        {1e-7, 0, 0, 0.5}, {1e-6, 3.6, 2.3e-5, 0.5}, {1e-5, 3.6, 2.3e-5, 0.5},
        {3e-5, 0, 0, 0.5}, {1e-5, 3.6, 2.3e-5, 1.0}, {1e-5, 3.6, 2.3e-5, 0.0},
    };
    double P = AMOLFRAM * 1e-6 / (111.0 * CONCM), TRf = TAURM * TAULM / (TAURM + TAULM);
    int c, k, i, fail = 0;

    for (c = 0; c < 6; ++c) {
        double TM = cases[c][0], DM = cases[c][1], DDM = cases[c][2], S2 = cases[c][3];
        for (k = 0; k <= 40; ++k) {
            double f = pow(10.0, 4.0 + k * 0.1), sum = 0, total, os, T1s, T1f, ref;
            for (i = 1; i <= 4; ++i) sum += LS(i, f, TM, DM, DDM, sqrt(S2));
            total = LS(0, f, TM, DM, DDM, sqrt(S2));
            os = N(4, f, TAURM, TM, DM, DDM);
            T1s = P / (N(0, f, TAURM, TM, DM, DDM) - os) - TM;
            T1f = P / (N(0, f, TRf, TM, DM, DDM) - N(4, f, TRf, TM, DM, DDM)) - TM;
            ref = P / (1.0 / (S2 / T1s + (1.0 - S2) / T1f) + TM) + os;
            if (fabs(sum / total - 1) > 1e-12 || fabs(total / ref - 1) > 1e-12 ||
                (S2 == 1.0 && fabs(total / N(0, f, TAURM, TM, DM, DDM) - 1) > 1e-12)) {
                printf("FAIL case %d f=%g: pieces %.15e total %.15e reference %.15e\n",
                       c, f, sum, total, ref);
                fail = 1;
            }
        }
    }
    if (!fail) printf("FlorenceLS: pieces = total = Lipari-Szabo reference\n");
    return fail;
}
