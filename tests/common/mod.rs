//! Shared helpers for integration tests.
//!
//! Three things live here:
//! - [`assert_exact`], the comparison used by every reference-table test
//!   against the F90 fixtures: bit for bit on the platform that generates
//!   them, within the last-bit differences of the system libm elsewhere.
//! - [`assert_close`] / [`assert_close_eps`], the tolerance-based comparison
//!   used by the round-trip tests.
//! - [`read_csv`], a tiny line-based reader for the fixture CSVs under
//!   `tests/data/`.
//!
//! The CSV format used by `tests/data/*.csv` is intentionally minimal: an
//! optional header line starting with `#`, followed by one record per line,
//! comma-separated, with every field a `f64` parsable by `f64::from_str`.

#![allow(dead_code)] // helpers are used by some test files, not all

use std::path::Path;

// ---------------------------------------------------------------------
// Tolerance constants
//
// The comparisons against the F90 fixtures under tests/data/ are exact on
// the platform that generates them (see assert_exact). The tolerances
// below serve the round-trip tests, which check the Rust port against
// itself.
// ---------------------------------------------------------------------

/// Default relative tolerance: one digit shy of `f64::EPSILON`.
pub const DEFAULT_REL_TOL: f64 = 1e-14;

/// Default absolute tolerance, used at boundary points where the
/// relative criterion is meaningless (e.g. `cdf(x) = 0.0`).
pub const DEFAULT_ABS_TOL: f64 = 1e-300;

/// `dinvr`-driven inverses and round-trip tests where the forward CDF
/// is computed by a direct or iterative routine. The search matches
/// CDFLIB's `dstinv` configuration with reltol = 1e-8; round-trip
/// residuals are bounded by that search tolerance plus the CDF's
/// Lipschitz factor near the queried quantile. 5e-8 leaves ~5x margin.
pub const INVERSE_REL_TOL: f64 = 5e-8;

/// `dinvr`-driven inverses where the forward CDF chains through an
/// iterative routine (`StudentsT::inverse_cdf`, `Beta::inverse_cdf`,
/// `ChiSquared::inverse_cdf`, …). With the search matching CDFLIB's
/// reltol = 1e-8, the worst-case projection through `1/|f'(x)|` near
/// low-pdf quantiles (e.g., t(df=4) at 0.975) reaches ~5e-7.
pub const CHAINED_INVERSE_REL_TOL: f64 = 5e-7;

/// Whether the tests run on the platform that generates the fixtures:
/// gfortran on macOS on Apple silicon, whose `log`, `exp` and `pow` come
/// from the same system libm that Rust calls there.
pub const REFERENCE_PLATFORM: bool = cfg!(all(target_os = "macos", target_arch = "aarch64"));

/// Relative tolerance of [`assert_exact`] on other platforms. Their libm
/// may differ from Apple's in the last bit, and conditioning amplifies
/// that: `exp` near underflow has condition number up to about 745. The
/// largest relative difference observed on x86_64 Linux (glibc 2.41) is
/// 2.3e-13.
pub const FOREIGN_LIBM_REL_TOL: f64 = 1e-12;

/// Assert that `got` and `expected` are the same IEEE binary64 value, bit
/// for bit, on the [`REFERENCE_PLATFORM`], and close within
/// [`FOREIGN_LIBM_REL_TOL`] elsewhere (any NaN matches any NaN). `ctx`,
/// typically the CSV row, is printed on failure.
///
/// The fixtures are generated with `-ffp-contract=off`, so the F90 results
/// contain no fused multiply-adds, and on the [`REFERENCE_PLATFORM`] the
/// Rust port reproduces them exactly. Elsewhere [`DEFAULT_ABS_TOL`] also
/// applies near zero.
#[track_caller]
pub fn assert_exact(got: f64, expected: f64, ctx: &[f64]) {
    if got.is_nan() && expected.is_nan() {
        return;
    }
    if REFERENCE_PLATFORM || got.to_bits() == expected.to_bits() {
        assert_eq!(
            got.to_bits(),
            expected.to_bits(),
            "got {got:e}, expected {expected:e} for {ctx:?}",
        );
        return;
    }
    let diff = (got - expected).abs();
    let rel = diff / got.abs().max(expected.abs());
    assert!(
        diff <= DEFAULT_ABS_TOL || rel <= FOREIGN_LIBM_REL_TOL,
        "got {got:e}, expected {expected:e} (relative difference {rel:e}) for {ctx:?}",
    );
}

/// Assert that `got` is close to `expected` using the default tolerances.
#[track_caller]
pub fn assert_close(got: f64, expected: f64) {
    assert_close_eps(got, expected, DEFAULT_REL_TOL, DEFAULT_ABS_TOL);
}

/// Assert that `got` is close to `expected` under a mixed relative/absolute
/// tolerance scheme.
///
/// The test passes if `|got - expected| ≤ abs_tol` OR `|got - expected| /
/// max(|got|, |expected|) ≤ rel_tol`. The denominator uses the maximum of
/// the two magnitudes so that the comparison is symmetric and well-defined
/// when one of them is exactly zero (in which case `abs_tol` carries the
/// test).
///
/// `NaN`/`Inf` mismatches always fail.
#[track_caller]
pub fn assert_close_eps(got: f64, expected: f64, rel_tol: f64, abs_tol: f64) {
    if got.is_nan() || expected.is_nan() {
        assert!(
            got.is_nan() && expected.is_nan(),
            "NaN mismatch: got {got}, expected {expected}",
        );
        return;
    }
    if got.is_infinite() || expected.is_infinite() {
        assert_eq!(
            got, expected,
            "infinity mismatch: got {got}, expected {expected}",
        );
        return;
    }

    let diff = (got - expected).abs();
    if diff <= abs_tol {
        return;
    }
    let scale = got.abs().max(expected.abs());
    let rel = if scale == 0.0 { 0.0 } else { diff / scale };
    assert!(
        rel <= rel_tol,
        "got {got}, expected {expected}, abs diff {diff}, rel diff {rel} (rel_tol {rel_tol}, abs_tol {abs_tol})",
    );
}

/// Read a reference-table CSV from disk.
///
/// Lines beginning with `#` are treated as comments/header and skipped.
/// Empty lines are skipped. Every remaining line must be a list of
/// comma-separated `f64` values, all rows the same length. Returns one
/// `Vec<f64>` per data row, and panics if there is none.
///
/// Paths are resolved relative to `CARGO_MANIFEST_DIR` so tests can refer
/// to fixtures by repository-relative paths like `"tests/data/erf.csv"`.
pub fn read_csv(rel_path: &str) -> Vec<Vec<f64>> {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let path = Path::new(manifest_dir).join(rel_path);
    let raw = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("could not read fixture {}: {e}", path.display()));

    let mut rows = Vec::new();
    let mut width: Option<usize> = None;
    for (lineno, line) in raw.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let row: Vec<f64> = trimmed
            .split(',')
            .map(|cell| {
                cell.trim().parse::<f64>().unwrap_or_else(|e| {
                    panic!(
                        "{}:{}: could not parse {cell:?} as f64: {e}",
                        path.display(),
                        lineno + 1,
                    )
                })
            })
            .collect();
        if let Some(w) = width {
            assert_eq!(
                row.len(),
                w,
                "{}:{}: row width {} differs from earlier rows ({w})",
                path.display(),
                lineno + 1,
                row.len(),
            );
        } else {
            width = Some(row.len());
        }
        rows.push(row);
    }
    assert!(
        !rows.is_empty(),
        "fixture {} has no data rows",
        path.display()
    );
    rows
}
