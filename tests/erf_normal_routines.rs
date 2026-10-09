#![cfg(not(miri))]

//! Reference-table tests for the erf and standard-normal routines.
//! Each row is compared against the Fortran reference output produced
//! by `tests/regenerate/gen_erf_normal_kernels.f90`.

mod common;

use cdflib::special::{cumnor, dinvnr, error_f, error_fc, error_fc_scaled};
use common::{assert_exact, read_csv};

#[test]
fn error_f_matches_reference() {
    for row in read_csv("tests/data/error_f.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(error_f(x), expected, &row);
    }
}

#[test]
fn error_fc_matches_reference() {
    for row in read_csv("tests/data/error_fc.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(error_fc(x), expected, &row);
    }
}

#[test]
fn error_fc_scaled_matches_reference() {
    for row in read_csv("tests/data/error_fc_scaled.csv") {
        let [x, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(error_fc_scaled(x), expected, &row);
    }
}

#[test]
fn cumnor_matches_reference() {
    for row in read_csv("tests/data/cumnor.csv") {
        let [x, expected_cum, expected_ccum] = row[..] else {
            panic!("width");
        };
        let (cum, ccum) = cumnor(x);
        assert_exact(cum, expected_cum, &row);
        assert_exact(ccum, expected_ccum, &row);
    }
}

#[test]
fn dinvnr_matches_reference() {
    for row in read_csv("tests/data/dinvnr.csv") {
        let [p, q, expected_x] = row[..] else {
            panic!("width");
        };
        assert_exact(dinvnr(p, q), expected_x, &row);
    }
}

// F90 error_f computes ax * (top / bot) in the |x| <= 0.5 branch and
// negates only when x < 0.0, so -0.0 yields +0.0; a NaN argument falls
// through every branch comparison into the saturated 1.0 tail and the
// final x < 0.0 test is false, so NaN yields exactly 1.0.
#[test]
fn error_f_signed_zero_and_nan_follow_f90() {
    assert_eq!(error_f(0.0).to_bits(), 0.0f64.to_bits());
    assert_eq!(error_f(-0.0).to_bits(), 0.0f64.to_bits());
    assert_eq!(error_f(f64::NAN), 1.0);
}
