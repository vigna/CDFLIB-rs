#!/bin/sh
# Checks that the tests comparing the Rust port with the reference tables
# under tests/data/ take both sides of every branch of the ported code,
# except the branch sides listed with a reason in
# tests/regenerate/branch_allowlist.txt. Run from the repository root:
# tests/regenerate/rust_branch_coverage.sh. Requires a nightly toolchain
# (TOOLCHAIN, default nightly), the llvm-tools component and python3.
#
# coverage.sh checks the F90 side line by line. This script checks the
# Rust side branch by branch, operands of && and || included: gfortran
# evaluates .and. and .or. without branching, so gcov cannot see which
# operand decided a condition.
#
# Tests: the integration tests in TESTS, which compare with tests/data/,
# except those in SKIP, which check Rust-only behavior, and the unit
# tests in UNIT_TESTS, which replay F90 traces and tables. The script
# fails if an integration test reads tests/data/ and is not in TESTS.
#
# Scope: the functions in src/special, src/search and src/dist, outside
# the test modules, except the Rust-only convenience functions, which are
# excluded by name (pdf, ln_pdf, pmf, ln_pmf, mean, variance, entropy and
# fmt; the accessors have no branches) or by a documentation comment
# starting with "Rust only:". A closure belongs to the function that
# contains it.
#
# Allowlist: one branch side per line, in the form
#   <file> <function> [<condition>] <true|false>: <reason>
# where the condition is its source text with whitespace collapsed; an
# entry stands for every branch with the same file, function, condition
# and side. Lines starting with # are comments. The script fails on every
# untaken branch side the allowlist does not list, and on every entry for
# a side that the tests take or that does not exist.
#
# The build goes to target/branch-coverage, so the usual target
# directory is not touched.

set -eu

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"

TOOLCHAIN=${TOOLCHAIN:-nightly}
TARGET="$ROOT/target/branch-coverage"
PROFILES="$TARGET/profiles"
ALLOWLIST="$ROOT/tests/regenerate/branch_allowlist.txt"

TESTS="
    beta_distributions
    beta_routines
    discrete_distributions
    dispatcher_calls
    dispatchers
    erf_normal_routines
    gamma_routines
    kernel_coverage
    noncentral_distributions
    normal_distribution
    psi_dt1_routines
"
SKIP="
    cdfbin_start_outside_search_range_is_an_error
    error_f_signed_zero_and_nan_follow_f90
    gamma_inc_accuracy_envelopes
"
UNIT_TESTS="
    solver_traces_match_f90
    exparg_matches_reference
"

# The lists must match the test files.
status=0
for f in tests/*.rs; do
    name=$(basename "$f" .rs)
    reads=no
    if grep -q 'read_csv(.*tests/data' "$f"; then
        reads=yes
    fi
    case " $(echo $TESTS) " in
        *" $name "*) listed=yes ;;
        *) listed=no ;;
    esac
    if [ $reads != $listed ]; then
        echo "tests/$name.rs: reads tests/data/: $reads, listed in TESTS: $listed" >&2
        status=1
    fi
done
for t in $SKIP; do
    if ! grep -q "fn $t()" tests/*.rs; then
        echo "no integration test named $t (SKIP)" >&2
        status=1
    fi
done
for t in $UNIT_TESTS; do
    if ! grep -rq "fn $t()" src; then
        echo "no unit test named $t (UNIT_TESTS)" >&2
        status=1
    fi
done
if [ $status -ne 0 ]; then
    exit 1
fi

# llvm-profdata and llvm-cov must come from the same LLVM major version as
# the compiler.
if ! rustc +"$TOOLCHAIN" -vV > /dev/null 2>&1; then
    echo "toolchain $TOOLCHAIN not found; install it with: rustup toolchain install $TOOLCHAIN" >&2
    exit 1
fi
HOST=$(rustc +"$TOOLCHAIN" -vV | sed -n 's/^host: //p')
LLVM_MAJOR=$(rustc +"$TOOLCHAIN" -vV | sed -n 's/^LLVM version: \([0-9]*\).*/\1/p')
LLVM_BIN=""
for dir in "$(rustc +"$TOOLCHAIN" --print sysroot)/lib/rustlib/$HOST/bin" \
    "$(rustc --print sysroot)/lib/rustlib/$HOST/bin"; do
    if [ -x "$dir/llvm-profdata" ] && [ -x "$dir/llvm-cov" ]; then
        major=$("$dir/llvm-profdata" --version | sed -n 's/.*LLVM version \([0-9]*\).*/\1/p' | head -n 1)
        if [ "$major" = "$LLVM_MAJOR" ]; then
            LLVM_BIN=$dir
            break
        fi
    fi
done
if [ -z "$LLVM_BIN" ]; then
    echo "llvm-profdata and llvm-cov for LLVM $LLVM_MAJOR not found; install them with:" >&2
    echo "    rustup component add llvm-tools --toolchain $TOOLCHAIN" >&2
    exit 1
fi

TEST_ARGS=""
for t in $TESTS; do
    TEST_ARGS="$TEST_ARGS --test $t"
done
SKIP_ARGS=""
for t in $SKIP; do
    SKIP_ARGS="$SKIP_ARGS --skip $t"
done

# The .cargo/config.toml flags are replaced by RUSTFLAGS, so target-cpu is
# repeated here.
RUSTFLAGS="-C target-cpu=native -C instrument-coverage -Z coverage-options=branch"
CARGO_TARGET_DIR="$TARGET"
LLVM_PROFILE_FILE="$PROFILES/%p-%m.profraw"
export RUSTFLAGS CARGO_TARGET_DIR LLVM_PROFILE_FILE

rm -rf "$PROFILES"
mkdir -p "$PROFILES"
EXECUTABLES=$(cargo +"$TOOLCHAIN" test --quiet --no-run --message-format=json --lib $TEST_ARGS |
    sed -n 's/.*"executable":"\([^"]*\)".*/\1/p')
cargo +"$TOOLCHAIN" test --quiet $TEST_ARGS -- $SKIP_ARGS
cargo +"$TOOLCHAIN" test --quiet --lib -- $UNIT_TESTS

"$LLVM_BIN/llvm-profdata" merge -sparse "$PROFILES"/*.profraw -o "$TARGET/tests.profdata"
OBJECTS=""
for e in $EXECUTABLES; do
    OBJECTS="$OBJECTS -object $e"
done
"$LLVM_BIN/llvm-cov" export -format=text -instr-profile="$TARGET/tests.profdata" \
    $OBJECTS > "$TARGET/coverage.json"

python3 -I - "$ROOT" "$TARGET/coverage.json" "$ALLOWLIST" << 'EOF'
import json
import os
import re
import sys

root, export, allowlist = sys.argv[1:4]
SCOPE = ("src/special/", "src/search/", "src/dist/")
EXCLUDED = {"pdf", "ln_pdf", "pmf", "ln_pmf", "mean", "variance", "entropy", "fmt"}
FN = re.compile(r"\bfn\s+([A-Za-z_][A-Za-z0-9_]*)")
ENTRY = re.compile(r"^(\S+) (\S+) \[(.*)\] (true|false): (\S.*)$")

sources = {}


def source(rel):
    if rel not in sources:
        with open(os.path.join(root, rel), "rb") as fh:
            sources[rel] = fh.read().decode().split("\n")
    return sources[rel]


def relative(path):
    path = os.path.realpath(path)
    prefix = os.path.realpath(root) + os.sep
    if not path.startswith(prefix):
        return None
    rel = path[len(prefix):]
    return rel if rel.startswith(SCOPE) else None


test_modules = {}


def test_module(rel):
    # The first line of the test module, or one past the end of the file.
    if rel not in test_modules:
        lines = source(rel)
        start = len(lines) + 1
        for i, line in enumerate(lines):
            if line.strip() == "#[cfg(test)]":
                rest = [l for l in lines[i + 1:] if l.strip()]
                if rest and rest[0].lstrip().startswith("mod "):
                    start = i + 1
                    break
        test_modules[rel] = start
    return test_modules[rel]


def function(rel, line):
    # The function defined at or above a line, and the line of its fn.
    lines = source(rel)
    for i in range(line - 1, -1, -1):
        m = FN.search(lines[i])
        if m:
            return m.group(1), i + 1
    return "?", 0


def rust_only(rel, fn_line):
    # Whether the documentation of the function starts with "Rust only:".
    lines = source(rel)
    i = fn_line - 2
    while i >= 0 and lines[i].strip().startswith("#["):
        i -= 1
    first = None
    while i >= 0 and lines[i].strip().startswith("///"):
        first = lines[i].strip()
        i -= 1
    return first is not None and first.startswith("/// Rust only:")


def condition(rel, l1, c1, l2, c2):
    # Columns are 1-based byte offsets; the end is exclusive.
    lines = [l.encode() for l in source(rel)]
    if l1 == l2:
        text = lines[l1 - 1][c1 - 1:c2 - 1]
    else:
        text = b" ".join([lines[l1 - 1][c1 - 1:]] + lines[l1:l2 - 1] + [lines[l2 - 1][:c2 - 1]])
    return " ".join(text.decode().split())


with open(export) as fh:
    data = json.load(fh)

# Branch counts summed over every instantiation in every test binary.
sites = {}
for d in data["data"]:
    for func in d["functions"]:
        files = func["filenames"]
        first = func["regions"][0]
        rel = relative(files[first[5]])
        if rel is None or first[0] >= test_module(rel):
            continue
        name, fn_line = function(rel, first[0])
        if name in EXCLUDED or rust_only(rel, fn_line):
            continue
        for b in func.get("branches", []):
            l1, c1, l2, c2, taken, not_taken, file_id = b[:7]
            if relative(files[file_id]) != rel:
                continue
            site = sites.setdefault((rel, l1, c1, l2, c2), [name, 0, 0])
            site[1] += taken
            site[2] += not_taken

untaken = {}
for (rel, l1, c1, l2, c2), (name, taken, not_taken) in sites.items():
    cond = condition(rel, l1, c1, l2, c2)
    for side, count in (("true", taken), ("false", not_taken)):
        if count == 0:
            untaken.setdefault((rel, name, cond, side), []).append(l1)

status = 0
listed = {}
with open(allowlist) as fh:
    for n, line in enumerate(fh, 1):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"):
            continue
        m = ENTRY.match(line)
        if m is None:
            print(f"branch_allowlist.txt:{n}: malformed entry: {line}")
            status = 1
            continue
        key = m.group(1, 2, 3, 4)
        if key in listed:
            print(f"branch_allowlist.txt:{n}: duplicate of line {listed[key]}")
            status = 1
        listed[key] = n

missing = sorted(k for k in untaken if k not in listed)
stale = sorted((n, k) for k, n in listed.items() if k not in untaken)

if missing:
    print("Branch sides that the tests never take and branch_allowlist.txt does not list:")
    for rel, name, cond, side in missing:
        lines = ", ".join(str(l) for l in sorted(untaken[(rel, name, cond, side)]))
        print(f"  {rel} {name} [{cond}] {side} (line {lines})")
    status = 1
if stale:
    print("Entries of branch_allowlist.txt for branch sides that the tests take or that do not exist:")
    for n, (rel, name, cond, side) in stale:
        print(f"  line {n}: {rel} {name} [{cond}] {side}")
    status = 1
if status:
    sys.exit(1)
print(f"{len(sites)} branches, {len(untaken)} untaken branch sides, all listed in")
print("branch_allowlist.txt, and every listed side is untaken.")
EOF
