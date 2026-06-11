#!/bin/sh
# Regenerate reference tables under tests/data/ from the bundled
# Fortran cdflib.f90 source. Run from the repository root:
# `tests/regenerate/regenerate.sh`.
#
# Re-run when the parameter grid or the cdflib.f90 reference changes;
# the generated CSVs are committed.

set -eu

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

BUILD_DIR="tests/regenerate/build"
mkdir -p "$BUILD_DIR" tests/data

REFS_DIR="tests/regenerate/refs"
# One generator per logical group. Each writes one or more CSVs under
# tests/data/. The Rust test files in tests/ read those CSVs.
GENERATORS="
    erf_normal_kernels
    normal_distribution
    gamma_kernels
    beta_kernels
    beta_distributions
    discrete_distributions
    noncentral_distributions
    dispatchers
    psi_dt1_kernels
"
for name in $GENERATORS; do
    SRC="tests/regenerate/gen_${name}.f90"
    if [ ! -f "$SRC" ]; then
        continue
    fi
    BIN="$BUILD_DIR/gen_${name}"
    echo "compiling $SRC -> $BIN"
    # -fdefault-real-8 promotes unsuffixed literals (e.g. the bare 0.15
    # and 0.1 in gamma_inc_inv) to double precision; -fdefault-double-8
    # keeps D+00 literals at 8 bytes instead of letting the first flag
    # push them to 16. Together every constant is IEEE binary64, matching
    # the Rust port.
    gfortran -O2 -fdefault-real-8 -fdefault-double-8 \
        -Wall -Wno-unused-variable -Wno-unused-dummy-argument \
        -J "$BUILD_DIR" \
        -o "$BIN" "$SRC" "$REFS_DIR/cdflib.f90"
    echo "running $BIN"
    "$BIN"
done

echo "done"
