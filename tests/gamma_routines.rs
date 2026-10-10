#![cfg(not(miri))]

//! Reference-table tests for the gamma routines.

mod common;

use cdflib::special::internal::dstrem;
use cdflib::special::{
    gamma, gamma_inc, gamma_inc_inv, gamma_inc_with_acc, gamma_log, GammaIncAcc,
};
use common::{assert_exact, read_csv};

#[test]
fn gamma_log_matches_reference() {
    for row in read_csv("tests/data/gamma_log.csv") {
        let [a, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(gamma_log(a), expected, &row);
    }
}

#[test]
fn gamma_matches_reference() {
    for row in read_csv("tests/data/gamma.csv") {
        let [a, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(gamma(a), expected, &row);
    }
}

#[test]
fn gamma_inc_matches_reference() {
    for row in read_csv("tests/data/gamma_inc.csv") {
        let [a, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = gamma_inc(a, x);
        assert_exact(p, expected_p, &row);
        assert_exact(q, expected_q, &row);
    }
}

// Truncation-depth fidelity at Digits6 and Digits3: Rust and the F90
// reference use the same shallower expansions at these accuracy levels, so
// the results agree exactly although the regimes are nominally accurate to
// 6 and 3 digits only.
#[test]
fn gamma_inc_digits6_matches_reference() {
    for row in read_csv("tests/data/gamma_inc_d6.csv") {
        let [a, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = gamma_inc_with_acc(a, x, GammaIncAcc::Digits6);
        assert_exact(p, expected_p, &row);
        assert_exact(q, expected_q, &row);
    }
}

#[test]
fn gamma_inc_digits3_matches_reference() {
    for row in read_csv("tests/data/gamma_inc_d3.csv") {
        let [a, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = gamma_inc_with_acc(a, x, GammaIncAcc::Digits3);
        assert_exact(p, expected_p, &row);
        assert_exact(q, expected_q, &row);
    }
}

// The declared envelope of each accuracy regime: Digits6 stays within ~1e-6
// of Max, Digits3 within ~1e-3 of Max. Across the full reference grid.
#[test]
fn gamma_inc_accuracy_envelopes() {
    for row in read_csv("tests/data/gamma_inc.csv") {
        let [a, x, _, _] = row[..] else {
            panic!("width");
        };
        let (p_max, _) = gamma_inc(a, x);
        let (p6, _) = gamma_inc_with_acc(a, x, GammaIncAcc::Digits6);
        let (p3, _) = gamma_inc_with_acc(a, x, GammaIncAcc::Digits3);
        assert!(
            (p_max - p6).abs() <= 1.0e-6,
            "Digits6 envelope busted at a={a}, x={x}: |{p_max} - {p6}| > 1e-6",
        );
        assert!(
            (p_max - p3).abs() <= 1.0e-3,
            "Digits3 envelope busted at a={a}, x={x}: |{p_max} - {p3}| > 1e-3",
        );
    }
}

#[test]
fn gamma_inc_inv_matches_reference() {
    for row in read_csv("tests/data/gamma_inc_inv.csv") {
        let [a, p, q, expected_x, expected_ierr] = row[..] else {
            panic!("width");
        };
        let (x, ierr) = gamma_inc_inv(a, -1.0, p, q);
        assert_exact(x, expected_x, &row);
        assert_eq!(f64::from(ierr), expected_ierr, "{row:?}");
    }
}

#[test]
fn dstrem_matches_reference() {
    for row in read_csv("tests/data/dstrem.csv") {
        let [z, expected] = row[..] else {
            panic!("width");
        };
        assert_exact(dstrem(z), expected, &row);
    }
}

// The a <= 0, |a| < 15 branch of gamma must keep its negative rows in the
// fixture: the -1.0e-8 row in particular distinguishes the two-step
// (x + 0.5) + 0.5 rounding from a single x + 1.0 rounding.
#[test]
fn gamma_negative_small_matches_reference_exactly() {
    let mut rows = 0;
    for row in read_csv("tests/data/gamma.csv") {
        let [a, expected] = row[..] else {
            panic!("width");
        };
        if a >= 0.0 || a <= -15.0 {
            continue;
        }
        let got = gamma(a);
        assert_eq!(
            got.to_bits(),
            expected.to_bits(),
            "a = {a}: got {got:e}, expected {expected:e}"
        );
        rows += 1;
    }
    assert!(rows >= 12, "fixture lost its negative rows ({rows})");
}
