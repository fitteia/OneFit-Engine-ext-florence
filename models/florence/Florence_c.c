#include <math.h>
#include <string.h>
#include "Florence_f.h"
#include "Florence_c.h"

/* A small exact-input cache shares the five outputs across all public wrappers.
 * Unlike the historical call-count cache, changing ANY input forces evaluation.
 * Eight slots retain interleaved species/components without unbounded growth.
 * No approximate matching: optimizers must see every parameter perturbation.
 * Calls must remain serialized: the Fortran backend uses shared COMMON blocks.
 * This cache does not make that backend thread-safe.
 */
#define FLORENCE_CACHE_SIZE 8
struct florence_cache_entry {
    int valid;
    double frequency;
    double parameters[28];
    double result[10];
};
static struct florence_cache_entry florence_cache[FLORENCE_CACHE_SIZE];
static unsigned int florence_next;

/* All five outputs RUNO(1..5) = total, dipolar, contact, cross, outer sphere
 * for one input set: one Fortran run, or none when the cache has it. */
void FlorenceRUNO(double frequency, const double parameters[28], double runo[5])
{
    unsigned int i;
    double input[28], log_frequency, result[10] = {0};
    struct florence_cache_entry *entry;

    /* Check before log10. */
    if (!isfinite(frequency) || frequency <= 0.0) {
        for (i = 0; i < 5; ++i) runo[i] = NAN;
        return;
    }
    for (i = 0; i < 28; ++i)
        if (!isfinite(parameters[i])) {
            for (i = 0; i < 5; ++i) runo[i] = NAN;
            return;
        }
    for (i = 0; i < FLORENCE_CACHE_SIZE; ++i) {
        entry = &florence_cache[i];
        if (entry->valid && entry->frequency == frequency &&
            memcmp(entry->parameters, parameters, sizeof(entry->parameters)) == 0) {
            memcpy(runo, entry->result, 5 * sizeof(double));
            return;
        }
    }
    memcpy(input, parameters, sizeof(input));
    log_frequency = log10(frequency);
    if (parameters[27] >= 2.0 && parameters[27] < 3.0)
        modflor_(input, &log_frequency, result);
    else
        florencef77_(input, &log_frequency, result);
    memcpy(runo, result, 5 * sizeof(double));
    /* Do not retain numerical failures. Only the first five slots are outputs
     * of the legacy backend, despite its ten-element array declaration. */
    for (i = 0; i < 5; ++i)
        if (!isfinite(result[i])) return;
    entry = &florence_cache[florence_next];
    entry->frequency = frequency;
    memcpy(entry->parameters, parameters, sizeof(entry->parameters));
    memcpy(entry->result, result, sizeof(result));
    entry->valid = 1;
    florence_next = (florence_next + 1) % FLORENCE_CACHE_SIZE;
}

static double florence_component(double frequency, const double parameters[28],
                                 double index)
{
    double runo[5];

    /* Keep the historical finite selector convention: truncate 1..4,
     * otherwise select total. */
    if (!isfinite(index))
        return NAN;
    FlorenceRUNO(frequency, parameters, runo);
    return runo[index >= 1.0 && index < 5.0 ? (int)index : 0];
}

double FlorenceN(
	  double index,
	  double n,
      double FREQ,
      double SI,
      double GAMMAI,
      double SPIN,
      double IREL,
      double DELTA2,
      double TAURM,
      double TAUVM,
      double TAUMM,
      double TAUDELTA,
      double DPARAM,
      double EPARAM,
      double S4M,
      double GXM,
      double GYM,
      double GZM,
      double AXM,
      double AYM,
      double AZM,
      double DM,
      double DDM,
      double CONCM,
      double ACQ,
      double AMOLFRAM,
      double RKM,
      double ACONTM,
      double THETAM,
      double PHIM,
      double FLAG
)
{
    const double parameters[28] = { SI, GAMMAI, SPIN, IREL, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, S4M, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, ACQ, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG };
    (void)n; /* Legacy grouping hint; cache identity now uses actual inputs. */
    return florence_component(FREQ, parameters, index);
}

double FlorenceN4LS(
	  double index,
	  double n,
      double FREQ,
      double SI,
      double GAMMAI,
      double SPIN,
      double IREL,
      double DELTA2,
      double TAURM,
      double TAUVM,
      double TAUMM,
      double TAUDELTA,
      double DPARAM,
      double EPARAM,
      double S4M,
      double GXM,
      double GYM,
      double GZM,
      double AXM,
      double AYM,
      double AZM,
      double DM,
      double DDM,
      double CONCM,
      double ACQ,
      double AMOLFRAM,
      double RKM,
      double ACONTM,
      double THETAM,
      double PHIM,
      double FLAG
)
{
    const double parameters[28] = { SI, GAMMAI, SPIN, IREL, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, S4M, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, ACQ, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG };
    (void)n; /* Legacy grouping hint; cache identity now uses actual inputs. */
    return florence_component(FREQ, parameters, index);
}

double Florence(
      double FREQ,
      double SI,
      double GAMMAI,
      double SPIN,
      double IREL,
      double DELTA2,
      double TAURM,
      double TAUVM,
      double TAUMM,
      double TAUDELTA,
      double DPARAM,
      double EPARAM,
      double S4M,
      double GXM,
      double GYM,
      double GZM,
      double AXM,
      double AYM,
      double AZM,
      double DM,
      double DDM,
      double CONCM,
      double ACQ,
      double AMOLFRAM,
      double RKM,
      double ACONTM,
      double THETAM,
      double PHIM,
      double FLAG
)
{
    const double parameters[28] = { SI, GAMMAI, SPIN, IREL, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, S4M, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, ACQ, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG };
    return florence_component(FREQ, parameters, 0.0);
}

double Florence4(
	  double index,
      double FREQ,
      double SI,
      double GAMMAI,
      double SPIN,
      double IREL,
      double DELTA2,
      double TAURM,
      double TAUVM,
      double TAUMM,
      double TAUDELTA,
      double DPARAM,
      double EPARAM,
      double S4M,
      double GXM,
      double GYM,
      double GZM,
      double AXM,
      double AYM,
      double AZM,
      double DM,
      double DDM,
      double CONCM,
      double ACQ,
      double AMOLFRAM,
      double RKM,
      double ACONTM,
      double THETAM,
      double PHIM,
      double FLAG
)
{
    const double parameters[28] = { SI, GAMMAI, SPIN, IREL, DELTA2, TAURM, TAUVM, TAUMM, TAUDELTA, DPARAM, EPARAM, S4M, GXM, GYM, GZM, AXM, AYM, AZM, DM, DDM, CONCM, ACQ, AMOLFRAM, RKM, ACONTM, THETAM, PHIM, FLAG };
    return florence_component(FREQ, parameters, index);
}
