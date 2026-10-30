# A converged Florence fit with known answers

`FlorenceLS-R1.dat` is R1 against the proton Larmor frequency from
FlorenceLS with known parameters (`make_data.c`: TAURM = 2.88e-11 s,
DELTA2 = 0.0236 x 1e20 s^-2, SLS = 0.7, TAUVM = 2.86e-12 s, Gd(III)-like
S = 7/2), 31 frequencies from 1e4 to 3.2e8 Hz, with a fixed 1% Gaussian
noise sequence and 1% errors. `fit_reference.t` fits TAURM, DELTA2 and
SLS twice, from below and from above the true values (20-30% off), and
checks that each fit:

- converges (MIGRAD) with an accurate error matrix, and MINOS succeeds;
- each fitted value lies within 3 sigma of the true one, with MINOS
  errors within 20% of the parabolic ones;
- ends with the reference chi2 (38.331 for 28 degrees of freedom),
  values and errors - so both reach the same minimum.

MINUIT prints the error-matrix state either as "ERROR MATRIX ACCURATE" or,
once it has an estimate, as "ERROR MATRIX UNCERTAINTY= x%", for accurate
and approximate matrices alike; `FitRef.pm` then reads the state from the
parameter table MINUIT prints after MIGRAD (headed APPROXIMATE when it is
not accurate) and from MINUIT's NOT POS-DEF/APPROXIMATE notes.

Reference (OFE-Go 0.5.5 and OFE-Raku 1.1.0, identical):

| | true | fitted | error |
| --- | --- | --- | --- |
| TAURM | 2.88e-11 | 2.8804e-11 | 0.0505e-11 |
| DELTA2 | 0.0236 | 0.023611 | 0.000065 |
| SLS | 0.7 | 0.70143 | 0.0068 |

TAUVM is fixed: these data depend on DELTA2 x TAUVM only (with both
free, their correlation is -1.000). TAURM and SLS stay correlated
(-0.98), as the Lipari-Szabo form makes them, but both are well
determined.

It needs an installed OneFit Engine with this extension, so the
repository's CI, which builds the extension alone, does not run it:

```sh
ONEFITE=/opt/onefite-go/bin/onefite OFE_PATH=/opt/onefite-go prove -v fit_reference.t
```

Without `ONEFITE` it uses `onefite-go`, else `onefite`, from PATH, and
skips when there is none. With `OFE_PATH` it also rebuilds the data with
`make_data.c` and compares.

## Are the errors right?

One fit shows the errors are consistent with the truth, not that they are
right. `pulls.pl` refits data sets that differ only in their noise
(`make_data.c SEED`, seeds 1..N) and reports, per parameter, the pulls
(fitted - true) / error. Run on .113 with OFE-Go 0.5.5 (2026-10-30,
about 2 minutes):

```
50 fits, 0 not converged cleanly
TR  pull mean +0.046  spread 1.064  within 1 sigma 64%
D2  pull mean -0.037  spread 0.962  within 1 sigma 68%
S   pull mean -0.081  spread 1.046  within 1 sigma 68%
```

Means near 0 and spreads near 1 (for 50 fits, within about 0.14 and
0.10), with about 68% within one sigma: the errors MINUIT reports for
this fit are the right size.

```sh
OFE_PATH=/opt/onefite-go ONEFITE=/opt/onefite-go/bin/onefite perl pulls.pl 50
```

Unlike the timing fixtures in `../perf-comparison/`, which do not
converge, this fit is a reference for parameters and uncertainties.
