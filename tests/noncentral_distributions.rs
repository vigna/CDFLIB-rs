#![cfg(not(miri))]

//! Reference-table tests for the noncentral chi-squared and noncentral
//! F distributions (CDFLIB's `cumchn` and `cumfnc`).
//!
//! Both functions are Poisson-mixture series whose internal convergence
//! tolerances are configured to `1e-5` (chi²) and `1e-4` (F) inside
//! CDFLIB. The assertions here are tuned to the measured error on the
//! committed fixture grid, not to machine epsilon.

mod common;

use cdflib::{ChiSquaredNoncentral, ContinuousCdf, FisherSnedecorNoncentral};
use common::{assert_exact, read_csv};

#[test]
fn chi_squared_noncentral_matches_cumchn_reference() {
    for row in read_csv("tests/data/chi_squared_noncentral_cdf.csv") {
        let [df, ncp, x, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = ChiSquaredNoncentral::new(df, ncp);
        assert_exact(d.cdf(x), expected_cdf, &row);
        assert_exact(d.ccdf(x), expected_sf, &row);
    }
}

#[test]
fn f_noncentral_matches_cumfnc_reference() {
    for row in read_csv("tests/data/fisher_snedecor_noncentral_cdf.csv") {
        let [dfn, dfd, ncp, fx, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = FisherSnedecorNoncentral::new(dfn, dfd, ncp);
        assert_exact(d.cdf(fx), expected_cdf, &row);
        assert_exact(d.ccdf(fx), expected_sf, &row);
    }
}
