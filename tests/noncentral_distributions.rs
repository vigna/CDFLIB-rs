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
use common::{
    assert_close_eps, read_csv, DEFAULT_ABS_TOL, NONCENTRAL_CHI_REL_TOL, NONCENTRAL_F_ABS_TOL,
    NONCENTRAL_F_REL_TOL,
};

#[test]
fn chi_squared_noncentral_matches_cumchn_reference() {
    for row in read_csv("tests/data/chi_squared_noncentral_cdf.csv") {
        let [df, ncp, x, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = ChiSquaredNoncentral::new(df, ncp);
        assert_close_eps(
            d.cdf(x),
            expected_cdf,
            NONCENTRAL_CHI_REL_TOL,
            DEFAULT_ABS_TOL,
        );
        assert_close_eps(
            d.ccdf(x),
            expected_sf,
            NONCENTRAL_CHI_REL_TOL,
            DEFAULT_ABS_TOL,
        );
    }
}

#[test]
fn f_noncentral_matches_cumfnc_reference() {
    for row in read_csv("tests/data/fisher_snedecor_noncentral_cdf.csv") {
        let [dfn, dfd, ncp, fx, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = FisherSnedecorNoncentral::new(dfn, dfd, ncp);
        assert_close_eps(
            d.cdf(fx),
            expected_cdf,
            NONCENTRAL_F_REL_TOL,
            NONCENTRAL_F_ABS_TOL,
        );
        assert_close_eps(
            d.ccdf(fx),
            expected_sf,
            NONCENTRAL_F_REL_TOL,
            NONCENTRAL_F_ABS_TOL,
        );
    }
}
