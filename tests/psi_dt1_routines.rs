#![cfg(not(miri))]

//! Reference-table tests for psi, dlanor, dt1, and stvaln.
//! Each row is compared against the Fortran reference output produced
//! by `tests/regenerate/gen_psi_dt1_kernels.f90`.

mod common;

use cdflib::special::internal::stvaln;
use cdflib::special::{dlanor, dt1, psi};
use common::{assert_close_eps, read_csv, DEFAULT_ABS_TOL, KERNEL_REL_TOL, STVALN_ABS_TOL};

#[test]
fn psi_matches_reference() {
    for row in read_csv("tests/data/psi.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(psi(x), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

#[test]
fn dlanor_matches_reference() {
    for row in read_csv("tests/data/dlanor.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(dlanor(x), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

#[test]
fn dt1_matches_reference() {
    for row in read_csv("tests/data/dt1.csv") {
        let [p, q, df, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(dt1(p, q, df), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

#[test]
fn stvaln_matches_reference() {
    for row in read_csv("tests/data/stvaln.csv") {
        let [p, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(stvaln(p), expected, KERNEL_REL_TOL, STVALN_ABS_TOL);
    }
}
