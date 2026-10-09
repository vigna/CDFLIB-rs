#![cfg(not(miri))]

//! Sanity tests for the crate scaffolding: error type, traits, and the
//! shared test helpers in `tests/common/mod.rs`.

mod common;

use cdflib::{ContinuousCdf, DiscreteCdf, SearchError};
use common::{assert_close, assert_close_eps, read_csv};
use std::f64::consts::PI;

// --- assert_close ---------------------------------------------------------

#[test]
fn assert_close_passes_when_within_default_tol() {
    assert_close(1.0 + 1e-15, 1.0);
    assert_close(-PI, -PI);
}

#[test]
fn assert_close_handles_boundary_values() {
    // At expected = 0.0 only the absolute tolerance applies.
    assert_close_eps(1e-310, 0.0, 1e-14, 1e-300);
    // At expected = 1.0 likewise.
    assert_close_eps(1.0, 1.0, 1e-14, 1e-300);
}

#[test]
#[should_panic(expected = "rel diff")]
fn assert_close_fails_outside_tol() {
    assert_close_eps(1.0, 2.0, 1e-14, 1e-300);
}

#[test]
#[should_panic(expected = "NaN mismatch")]
fn assert_close_fails_on_nan_mismatch() {
    assert_close_eps(f64::NAN, 0.0, 1e-14, 1e-300);
}

#[test]
fn assert_close_accepts_both_nan() {
    assert_close_eps(f64::NAN, f64::NAN, 1e-14, 1e-300);
}

#[test]
#[should_panic(expected = "infinity mismatch")]
fn assert_close_fails_on_infinity_mismatch() {
    assert_close_eps(f64::INFINITY, 1e308, 1e-14, 1e-300);
}

// --- read_csv -------------------------------------------------------------

#[test]
fn read_csv_parses_a_simple_table() {
    // Write a temporary fixture in the target directory, then read it back
    // by its absolute path, which read_csv does not join to the manifest
    // directory.
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("read_csv_simple.csv");
    std::fs::write(
        &path,
        "# header line\n1.0, 2.0, 3.0\n# comment in the middle\n4.5,5.5,6.5\n\n",
    )
    .unwrap();

    let rows = read_csv(path.to_str().unwrap());
    assert_eq!(rows, vec![vec![1.0, 2.0, 3.0], vec![4.5, 5.5, 6.5]]);
}

#[test]
#[should_panic(expected = "no data rows")]
fn read_csv_rejects_a_table_without_rows() {
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("read_csv_empty.csv");
    std::fs::write(&path, "# header line\n\n").unwrap();
    read_csv(path.to_str().unwrap());
}

// --- error / trait shape sanity ------------------------------------------

#[test]
fn search_error_displays_useful_messages() {
    let e = SearchError::AnswerBelowLowerBound { bound: -1.0 };
    assert!(e.to_string().contains("-1"), "got: {e}");
    let e = SearchError::AnswerAboveUpperBound { bound: 1.0 };
    assert!(e.to_string().contains('1'), "got: {e}");
}

// A tiny stub distribution to prove the trait shape compiles.

#[derive(Debug)]
struct StubContinuous;

impl ContinuousCdf for StubContinuous {
    type Error = SearchError;
    fn cdf(&self, x: f64) -> f64 {
        // Uniform on [0..1].
        x.clamp(0.0, 1.0)
    }
    fn inverse_cdf(&self, p: f64) -> Result<f64, SearchError> {
        Ok(p.clamp(0.0, 1.0))
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }
}

#[test]
fn continuous_cdf_trait_shape_compiles() {
    let d = StubContinuous;
    assert_close(d.ccdf(0.3), 0.7);
    assert_close(d.inverse_cdf(0.3).unwrap(), 0.3);
}

#[derive(Debug)]
struct StubDiscrete;

impl DiscreteCdf for StubDiscrete {
    type Error = SearchError;
    fn cdf(&self, x: u64) -> f64 {
        if x >= 1 {
            1.0
        } else {
            0.0
        }
    }
    fn ccdf(&self, x: u64) -> f64 {
        if x >= 1 {
            0.0
        } else {
            1.0
        }
    }
    fn inverse_cdf(&self, p: f64) -> Result<u64, SearchError> {
        Ok(if p > 0.0 { 1 } else { 0 })
    }
}

#[test]
fn discrete_cdf_trait_shape_compiles() {
    let d = StubDiscrete;
    assert_close(d.ccdf(0), 1.0);
    assert_close(d.ccdf(1), 0.0);
    assert_eq!(d.inverse_cdf(1.0).unwrap(), 1);
    assert_eq!(d.inverse_cdf(0.0).unwrap(), 0);
}
