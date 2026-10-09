#![cfg(not(miri))]

//! Reference-table tests for the beta routines.

mod common;

use cdflib::special::internal::dbetrm;
use cdflib::special::{beta_inc, beta_log};
use common::{assert_exact, read_csv};

#[test]
fn beta_log_matches_reference() {
    for row in read_csv("tests/data/beta_log.csv") {
        let [a, b, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(beta_log(a, b), expected, &row);
    }
}

#[test]
fn beta_inc_matches_reference() {
    for row in read_csv("tests/data/beta_inc.csv") {
        let [a, b, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = beta_inc(a, b, x, 1.0 - x);
        assert_exact(p, expected_p, &row);
        assert_exact(q, expected_q, &row);
    }
}

#[test]
fn dbetrm_matches_reference() {
    for row in read_csv("tests/data/dbetrm.csv") {
        let [a, b, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(dbetrm(a, b), expected, &row);
    }
}
