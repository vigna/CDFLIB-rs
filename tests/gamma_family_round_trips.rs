#![cfg(not(miri))]

//! Round-trip integration tests for ChiSquared and Gamma: invoke
//! `inverse_cdf` / `inverse_ccdf` / search, then verify by re-evaluating
//! `cdf` at the answer. All assertions go through `INVERSE_REL_TOL`.

mod common;

use cdflib::{ChiSquared, ContinuousCdf, Gamma};
use common::{assert_close_eps, INVERSE_REL_TOL};

#[test]
fn chi_squared_inverse_cdf_round_trip() {
    for &df in &[1.0, 2.0, 5.0, 10.0, 30.0, 100.0] {
        let c = ChiSquared::new(df);
        for &p in &[0.01, 0.05, 0.5, 0.9, 0.95, 0.99] {
            let x = c.inverse_cdf(p).unwrap();
            assert_close_eps(c.cdf(x), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
        }
    }
}

#[test]
fn chi_squared_inverse_ccdf_round_trip() {
    for &df in &[1.0, 2.0, 5.0, 10.0, 30.0] {
        let c = ChiSquared::new(df);
        for &q in &[0.01, 0.05, 0.1, 0.5] {
            let x = c.inverse_ccdf(q).unwrap();
            assert_close_eps(c.ccdf(x), q, INVERSE_REL_TOL, INVERSE_REL_TOL);
        }
    }
}

#[test]
fn chi_squared_search_df_round_trip() {
    for &(p_target, x) in &[(0.95, 3.84), (0.99, 6.63), (0.5, 2.0)] {
        let df = ChiSquared::search_df(p_target, 1.0 - p_target, x).unwrap();
        let p_back = ChiSquared::new(df).cdf(x);
        assert_close_eps(p_back, p_target, INVERSE_REL_TOL, INVERSE_REL_TOL);
    }
}

#[test]
fn gamma_inverse_cdf_round_trip() {
    for &(shape, rate) in &[(1.0, 1.0), (2.0, 3.0), (0.5, 2.0), (10.0, 0.5)] {
        let g = Gamma::new(shape, rate);
        for &p in &[0.01, 0.1, 0.5, 0.9, 0.99] {
            let x = g.inverse_cdf(p).unwrap();
            assert_close_eps(g.cdf(x), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
        }
    }
}

#[test]
fn gamma_search_round_trip() {
    let shape = Gamma::search_shape(0.95, 0.05, 5.0, 2.0).unwrap();
    let back = Gamma::new(shape, 2.0).cdf(5.0);
    assert_close_eps(back, 0.95, INVERSE_REL_TOL, INVERSE_REL_TOL);

    let rate = Gamma::search_rate(0.5, 0.5, 4.0, 2.0).unwrap();
    let back = Gamma::new(2.0, rate).cdf(4.0);
    assert_close_eps(back, 0.5, INVERSE_REL_TOL, INVERSE_REL_TOL);
}
