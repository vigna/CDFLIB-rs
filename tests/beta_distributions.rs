#![cfg(not(miri))]

//! Reference-table tests for the Beta family of distributions: Beta,
//! Student's *t*, Fisher–Snedecor (F). Each table was produced by the
//! Fortran `cumbet`/`cumt`/`cumf` routines from
//! `tests/regenerate/refs/cdflib.f90`.

mod common;

use cdflib::{Beta, ContinuousCdf, FisherSnedecor, StudentsT};
use common::{assert_exact, read_csv};

#[test]
fn beta_cdf_matches_cumbet_reference() {
    for row in read_csv("tests/data/beta_cdf.csv") {
        let [a, b, x, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = Beta::new(a, b);
        assert_exact(d.cdf(x), expected_cdf, &row);
        assert_exact(d.ccdf(x), expected_sf, &row);
    }
}

#[test]
fn students_t_cdf_matches_cumt_reference() {
    for row in read_csv("tests/data/students_t_cdf.csv") {
        let [df, t, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = StudentsT::new(df);
        assert_exact(d.cdf(t), expected_cdf, &row);
        assert_exact(d.ccdf(t), expected_sf, &row);
    }
}

#[test]
fn f_cdf_matches_cumf_reference() {
    for row in read_csv("tests/data/f_cdf.csv") {
        let [dfn, dfd, fx, expected_cdf, expected_sf] = row[..] else {
            panic!("width");
        };
        let d = FisherSnedecor::new(dfn, dfd);
        assert_exact(d.cdf(fx), expected_cdf, &row);
        assert_exact(d.ccdf(fx), expected_sf, &row);
    }
}
