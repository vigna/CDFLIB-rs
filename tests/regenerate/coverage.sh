#!/bin/sh
# Checks that the reference-table generators execute every line of
# tests/regenerate/refs/cdflib.f90, except the lines listed with a reason
# in tests/regenerate/unreachable.txt. Run from the repository root:
# `tests/regenerate/coverage.sh`. Requires gfortran and gcov.
#
# The generators run in a temporary directory, so the committed CSVs are
# not touched. The script prints every executable line that no generator
# reaches and that unreachable.txt does not list, and every listed line
# that a generator does reach, and exits with status 1 if there is any.
# Lines that only print messages are not counted.

set -eu

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
REGEN="$ROOT/tests/regenerate"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

mkdir -p "$WORK/tests/data"
cd "$WORK"

FLAGS="-O0 --coverage -fdefault-real-8 -fdefault-double-8 -ffp-contract=off"
# The same generators, in the same order, as regenerate.sh.
GENERATORS=$(sed -n '/^GENERATORS="/,/^"/p' "$REGEN/regenerate.sh" | grep -v '"')

gfortran $FLAGS -c "$REGEN/refs/cdflib.f90" -o cdflib.o
for name in $GENERATORS; do
    gfortran $FLAGS -o "gen_$name" "$REGEN/gen_$name.f90" cdflib.o
    "./gen_$name" > /dev/null
done
gcov -o . cdflib.o > /dev/null

# Executed-line report: unexecuted lines of the library routines (the
# *_values tables excluded), message-printing lines excluded.
awk -F: '
    { cnt = $1; gsub(/ /, "", cnt); ln = $2 + 0; src = $0; sub(/^[^:]*:[^:]*:/, "", src) }
    src ~ /^(recursive )?(subroutine|function) / {
        fn = src; sub(/^(recursive )?(subroutine|function) +/, "", fn); sub(/ *\(.*/, "", fn)
    }
    cnt == "#####" && fn !~ /_values$/ && src !~ /write \(/ && !seen[ln]++ { print ln, fn, src }
' cdflib.f90.gcov > uncovered.txt

# Line numbers listed in unreachable.txt, ranges expanded.
grep -v '^#' "$REGEN/unreachable.txt" | awk '{
    split($1, r, "-"); lo = r[1]; hi = (r[2] == "" ? r[1] : r[2])
    for (i = lo; i <= hi; i++) print i
}' | sort -n > listed.txt

awk 'FILENAME == ARGV[1] { listed[$1] = 1; next } !($1 in listed)' listed.txt uncovered.txt > missing.txt

# Lines that some generator executes: a positive count in any listing
# (gcov appends * to the count of a line some of whose blocks did not run).
awk -F: '
    { cnt = $1; gsub(/ /, "", cnt); sub(/\*$/, "", cnt) }
    cnt ~ /^[0-9]+$/ && cnt + 0 > 0 { print $2 + 0 }
' cdflib.f90.gcov | sort -un > executed.txt

# Listed lines that a generator executes after all.
awk 'FILENAME == ARGV[1] { executed[$1] = 1; next } ($1 in executed)' executed.txt listed.txt > stale.txt

status=0
if [ -s missing.txt ]; then
    echo "Lines of cdflib.f90 that no generator executes:"
    cat missing.txt
    status=1
fi
if [ -s stale.txt ]; then
    echo "Lines listed in unreachable.txt that a generator executes:"
    cat stale.txt
    status=1
fi
if [ $status -ne 0 ]; then
    exit 1
fi
echo "Every executable line of cdflib.f90 is reached by a generator or listed in"
echo "unreachable.txt, and no listed line is reached."
