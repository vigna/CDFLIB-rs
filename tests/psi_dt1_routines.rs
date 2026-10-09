#![cfg(not(miri))]

//! Reference-table tests for psi, dlanor, dt1, and stvaln.
//! Each row is compared against the Fortran reference output produced
//! by `tests/regenerate/gen_psi_dt1_kernels.f90`.

mod common;

use cdflib::special::internal::stvaln;
use cdflib::special::{dlanor, dt1, psi};
use common::{assert_exact, read_csv};

#[test]
fn psi_matches_reference() {
    for row in read_csv("tests/data/psi.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(psi(x), expected, &row);
    }
}

#[test]
fn dlanor_matches_reference() {
    for row in read_csv("tests/data/dlanor.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(dlanor(x), expected, &row);
    }
}

#[test]
fn dt1_matches_reference() {
    for row in read_csv("tests/data/dt1.csv") {
        let [p, q, df, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(dt1(p, q, df), expected, &row);
    }
}

#[test]
fn stvaln_matches_reference() {
    for row in read_csv("tests/data/stvaln.csv") {
        let [p, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(stvaln(p), expected, &row);
    }
}
