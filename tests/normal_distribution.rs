#![cfg(not(miri))]

//! Reference-table tests for the Normal distribution against `cdfnor`.

mod common;

use cdflib::{ContinuousCdf, Normal};
use common::{assert_exact, read_csv};

#[test]
fn cdf_and_ccdf_match_cdfnor_reference() {
    for row in read_csv("tests/data/normal_cdf.csv") {
        let [mean, sd, x, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let n = Normal::new(mean, sd);
        // Normal::cdf is a direct cumnor wrapper; no iterative routine.
        assert_exact(n.cdf(x), expected_cdf, &row);
        assert_exact(n.ccdf(x), expected_sf, &row);
    }
}

#[test]
fn inverse_cdf_matches_cdfnor_reference() {
    // Each row supplies (p, q) with full precision in both tails (one
    // generated as cum, the other as ccum, neither via 1 - other). The
    // trait single-arg API can only carry one of them with full
    // precision, so we route each row to its accurate branch:
    //   - p ≤ 0.5 → inverse_cdf(p) is the accurate call
    //   - q ≤ 0.5 → inverse_ccdf(q) is the accurate call
    for row in read_csv("tests/data/normal_inverse_cdf.csv") {
        let [mean, sd, p, q, expected_x] = row[..] else {
            panic!("width");
        };
        let n = Normal::new(mean, sd);
        if p <= 0.5 {
            assert_exact(n.inverse_cdf(p).unwrap(), expected_x, &row);
        }
        if q <= 0.5 {
            assert_exact(n.inverse_ccdf(q).unwrap(), expected_x, &row);
        }
    }
}
