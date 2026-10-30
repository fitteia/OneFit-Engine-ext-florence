/* Synthetic FlorenceLS R1 dispersion with known parameters and a fixed
 * pseudo-random noise sequence (1% Gaussian, LCG + Box-Muller), written as
 * "frequency R1 error" with error = 1% of the noiseless value. An argument
 * N > 0 starts the noise sequence from seed N instead (pulls.pl). */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Florence_c.h"

static unsigned long long lcg = 20261030ULL;
static double uniform(void)
{
    lcg = lcg * 6364136223846793005ULL + 1442695040888963407ULL;
    return ((lcg >> 11) + 0.5) / 9007199254740992.0;
}

int main(int argc, char **argv)
{
    const double TAURM = 2.88e-11, DELTA2 = 0.0236, TAUVM = 2.86e-12, SLS = 0.7;
    int k;
    if (argc > 1 && atoll(argv[1]) > 0)
        lcg = (unsigned long long)atoll(argv[1]);
    for (k = 0; k < 31; ++k) {
        double f = pow(10.0, 4.0 + k * 0.15);
        double r = FlorenceLS(0, f, 0.5, 2.675222e8, 3.5, DELTA2, TAURM, TAUVM, 1e-7, 2,
                              0, 0, 2.003, 2.003, 2.003, 0, 0, 0, 0.001, 3.1, 0, 0, 0.748,
                              12, 0, 0, SLS, 5e-12, 1);
        double g = sqrt(-2 * log(uniform())) * cos(2 * M_PI * uniform());
        printf("%.6e %.6e %.6e\n", f, r * (1 + 0.01 * g), 0.01 * r);
    }
    return 0;
}
