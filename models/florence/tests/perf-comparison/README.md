# Florence eigensolver replacement: reproducing the comparison

This directory reproduces the performance and accuracy comparison between
the **original NAG-derived** Florence complex-Hermitian eigensolver and its
**replacement** (`OFHEIG`, formerly `F02AXF`, backed by a hand-rolled Householder-
tridiagonalization + implicit-QL eigensolver — no NAG or LAPACK code at
all), fitting the same three real datasets against both.

## Why this replacement exists

The Florence model bundle (`OneFit-Engine-ext-florence`) previously
contained NAG-copyrighted Fortran routines (`F01BCF`, `F02AYF`, plus the
NAG-derived body of `F02AXF`) whose redistribution terms required review
before the bundle could be published. Those routines have been fully
removed and replaced with an independent, from-scratch implementation of
the same algorithm (Householder reduction + implicit QL with Wilkinson
shifts — standard, decades-old published numerical linear algebra, not
NAG's, LAPACK's, or EISPACK's code). See
`../NOTICE` for the full history and verification details.

## What to expect

Fitting is CPU-bound and this eigensolver is called tens of millions of
times per fit (once per MINUIT chi² evaluation), on small matrices
(N=6-8 for these three test files). The replacement is **honestly
slower** than the original NAG code — about **1.6x** on the fit measured
in most detail (`florenceN4LS_fH_R1OSabh_test.json`), 20.2s → 33.0s.
This isn't a bug: it's the real, measured cost of removing a dependency
that could not be redistributed. Three earlier replacement attempts
(LAPACK's `ZHEEV`, a Jacobi-method solver, and a `COMPLEX*16` version of
this same Householder+QL algorithm) were all measured and are documented
in `../NOTICE` for context — this is the fastest of everything tried.

Do not read chi² differences between the two trees as a parity test.
None of the three fits converges cleanly (audit 2026-10-30, OFE-Go
0.5.0): `florenceN_fH_test.json` (chi² about 15.0) and
`florenceN_fH_SS_test.json` (about 54) end MIGRAD at its call limit,
and `florenceN4LS_fH_R1OSabh_test.json` (about 1.44) without a
positive-definite error matrix, although each run exits 0. From such a
stopping point, round-off differences in the model can lead MINUIT to
a different chi², so the 0.3-6% differences recorded below say nothing
about either eigensolver, and the fitted parameters and errors of these
fits are not reference values. These files are kept as timing
benchmarks; `../fit-reference/` has a converged fit with known answers,
checked values and uncertainties. Parity is checked on the model itself: the golden test
(`../run_tests.sh`), the random-matrix sweep against `ZHEEV`
(`../test_cheevhqr.f`, linked with the production solver), and a
comparison of the Florence outputs of this bundle and the NAG one.

## A note on what's reproducible from this public repo

This package builds and lets you re-time the **NAG-free** tree only.
The original NAG-derived code is not in this repository (deliberately
- see the top-level README) and lives in a separate, private repo. If
you separately hold a NAG Library license and have access to that
repo, you can build a second tree from it yourself (same
`doctor -install -enable-extensions -extensions-ref=<their ref>`
pattern) and set `OLD_OFE_PATH` to it when running
`run-comparison.sh`, to reproduce the full old-vs-new comparison.
Otherwise, the recorded numbers below are the reference result.

## Prerequisites

- `git`, `go` (1.21+), `gfortran`, `gcc`, `make`, `ar` - the same
  prerequisites `onefite-go doctor` itself checks for.
- No special credentials needed for this package - everything it
  builds is public.

## Reproducing

```bash
./setup-comparison.sh                 # builds onefite-go + the engine tree
                                       # (~a couple minutes)

ONEFITE_GO="$PWD/comparison-workdir/onefite-go" \
OLD_OFE_PATH="$PWD/comparison-workdir/ofe-new" \
NEW_OFE_PATH="$PWD/comparison-workdir/ofe-new" \
./run-comparison.sh \
  florenceN_fH_SS_test.json \
  florenceN4LS_fH_R1OSabh_test.json \
  florenceN_fH_test.json
```

(`OLD_OFE_PATH`/`NEW_OFE_PATH` point at the same tree here since only
one is built; swap `OLD_OFE_PATH` for your own NAG-derived tree if you
have one, per the note above, to get a real old-vs-new run.)

`run-comparison.sh` fits each file 5 times against each tree (set
`NTRIALS=N` to change that) and writes a `report.txt` plus full
`fit.out`/`fit.log` for every run under a directory it prints at the end
(`OUTDIR`, or a fresh `mktemp -d` if you don't set it).

## Pinned versions (this comparison was originally run against)

| Repo | Ref |
|---|---|
| Florence extensions (this replacement) | `4e8acf36c174b0e4778dd022c02a3e755787ff4d` |
| `OneFit-Engine-go` (`onefite-go` itself) | `fd6d3de18f6e8d3f7ff2b82f96caa55ab3de06d4` |
| `onefite-c-code` (recorded, not independently pinnable by `doctor`) | `6633911acb2c376267c202e43dc64ed9e2b709df` |
| `minuit` (recorded, not independently pinnable by `doctor`) | `b842e5b02052f598f3f14606bcc01f78e95d8fd5` |

`onefite-c-code` and `minuit` always track their own default branch in
`doctor -install` (no per-tree pinning flag).

## Real numbers from the original run

5 trials each, `florenceN4LS_fH_R1OSabh_test.json` (N=8):

| Version | Time | vs. NAG | chi² |
|---|---|---|---|
| Original NAG | 20.2s | baseline | 1.44858 |
| This replacement (Householder+QL) | 33.0s | +63% | 1.534496 |

See `../NOTICE` for the full table including the three intermediate
attempts (`ZHEEV`, Jacobi, `COMPLEX*16` Householder+QL) and the two
performance angles investigated and ruled out (QL-iteration count on
real matrices; exact-input-matrix caching, which doesn't apply since
the eigensolver's input matrix changes on every single call within a
fit).
