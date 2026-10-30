#include <stdio.h>
#include "Florence_c.h"

/* A case the modified model (FLAG = 2) does not implement: zero-field
 * splitting (DPARAM != 0) with SI = 0.5 and a contact term. Florence must
 * fail visibly (message, exit status 1) - it used to stop silently with
 * status 0, which ended a fit as if it had succeeded. run_tests.sh checks
 * both; reaching the printf below is a failure. */
int main(void) {
    double v = FlorenceN(0, 5, 1e4, 0.5, 2.675222e8, 2.5, 1, 0.01, 1e-11,
        2.86e-12, 1e-7, 2, 0.05, 0.01, 0, 2.003, 2.003, 2.003, 0, 0, 0,
        2.2e-9, 3.6e-10, 0.001, 1, 12, 2.767, 0.748, 0, 0, 2);
    printf("returned %.6e instead of failing\n", v);
    return 0;
}
