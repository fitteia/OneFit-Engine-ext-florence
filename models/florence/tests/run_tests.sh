#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL="$(cd "$HERE/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "== building Florence extension =="
gfortran -c -ffixed-form -std=legacy "$MODEL/Florence_f.f" -o "$WORK/f.o"
gcc -c -I"$MODEL" "$MODEL/Florence_c.c" -o "$WORK/c.o"
gcc -c -I"$MODEL" "$MODEL/FlorenceLS.c" -o "$WORK/ls.o"
gcc -I"$MODEL" "$HERE/test_ls.c" "$WORK/f.o" "$WORK/c.o" "$WORK/ls.o" -lgfortran -lm -o "$WORK/test_ls"
gcc -I"$MODEL" "$HERE/test_driver.c" "$WORK/f.o" "$WORK/c.o" -lgfortran -lm -o "$WORK/driver"
gcc -I"$MODEL" "$HERE/staleness_demo.c" "$WORK/f.o" "$WORK/c.o" -lgfortran -lm -o "$WORK/staleness_demo"
gcc -I"$MODEL" "$HERE/test_unsupported.c" "$WORK/f.o" "$WORK/c.o" -lgfortran -lm -o "$WORK/test_unsupported"
gfortran -c -ffixed-form -std=legacy "$HERE/test_eigensolver.f" -o "$WORK/test_eigensolver.o"
gfortran "$WORK/test_eigensolver.o" "$WORK/f.o" -o "$WORK/test_eigensolver"
gfortran -ffixed-form -std=legacy "$HERE/test_solver_fatal.f" "$WORK/f.o" -o "$WORK/test_solver_fatal"

echo "== eigensolver known-answer regression =="
"$WORK/test_eigensolver"

echo "== Householder+QL eigensolver internals regression =="
gfortran -O2 -ffixed-form -std=legacy "$HERE/test_cheevhqr.f" "$WORK/f.o" -llapack -lblas -o "$WORK/test_cheevhqr"
"$WORK/test_cheevhqr"

echo "== golden regression =="
: > "$WORK/actual.txt"
for freq in 1e3 1e4 1e5 1e6 2.5e6 1e7 5e7 1e8 1e9 5e9; do
  for flag in 1 2; do
    { echo "=== freq=$freq flag=$flag ==="; "$WORK/driver" "$freq" "$flag"; } >> "$WORK/actual.txt"
  done
done
diff -u "$HERE/golden_output.txt" "$WORK/actual.txt"
echo "PASS: Florence golden output"

echo "== TAUDELTA pass-through regression =="
"$WORK/driver" 1e6 1 2 > "$WORK/taudelta2.txt"
"$WORK/driver" 1e6 1 1 > "$WORK/taudelta1.txt"
if cmp -s "$WORK/taudelta2.txt" "$WORK/taudelta1.txt"; then
  echo "FAIL: TAUDELTA does not reach the Fortran backend"
  exit 1
fi
echo "PASS: Florence TAUDELTA pass-through"

echo "== stale-cache regression =="
"$WORK/staleness_demo" | tee "$WORK/stale.txt"
grep -q 'MATCH: yes' "$WORK/stale.txt"
echo "PASS: Florence cache regression"

echo "== FlorenceLS regression =="
"$WORK/test_ls"
echo "PASS: FlorenceLS pieces and Lipari-Szabo total"

echo "== unsupported case fails visibly =="
if "$WORK/test_unsupported" > "$WORK/unsup.out" 2> "$WORK/unsup.err"; then
  echo "FAIL: unsupported FLAG=2 case exited with status 0"; cat "$WORK/unsup.out"; exit 1
fi
grep -q 'TUNOstock' "$WORK/unsup.err" || { echo "FAIL: no message for the unsupported case"; exit 1; }
echo "PASS: Florence unsupported case reports an error"

echo "== eigensolver failure stops with status 1 =="
if "$WORK/test_solver_fatal" > "$WORK/fatal.out" 2> "$WORK/fatal.err"; then
  echo "FAIL: an OFHEIG failure exited with status 0"; cat "$WORK/fatal.out"; exit 1
fi
grep -q 'eigensolver failed, code 2' "$WORK/fatal.err" || { echo "FAIL: no message for the OFHEIG failure"; cat "$WORK/fatal.err"; exit 1; }
echo "PASS: Florence eigensolver failure reports an error"
