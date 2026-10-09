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

// Truncation-depth fidelity at Digits6 and Digits3: both Rust and the F90
// reference use the same shallower expansion at these accuracy levels, so the
// per-row agreement should still be at iterative-routine precision (the
// regimes are nominally "6 digits" and "3 digits", but each implementation
// rounds the same way on the same expansion).
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
    // At the shallowest truncation, Rust and F90 share the same arithmetic
    // but the loose convergence tolerance (acc = 5e-4 → tol ≈ 2.5e-4 in the
    // Taylor and continued-fraction loops) lets one or two extra iterations
    // separate the two implementations on rows where the final iterate
    // straddles the tolerance. Empirically the worst row sits at
    // ~1.2e-5 absolute, which is one order tighter than the regime spec
    // (1 unit of the 3rd significant digit, ~1e-3). We assert the spec.
    for row in read_csv("tests/data/gamma_inc_d3.csv") {
        let [a, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = gamma_inc_with_acc(a, x, GammaIncAcc::Digits3);
        let dp = (p - expected_p).abs();
        let dq = (q - expected_q).abs();
        assert!(
            dp <= 1.0e-3,
            "Digits3 vs F90 ind=2 disagrees at a={a}, x={x}: |{p} - {expected_p}| = {dp}",
        );
        assert!(
            dq <= 1.0e-3,
            "Digits3 q disagrees at a={a}, x={x}: |{q} - {expected_q}| = {dq}",
        );
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
        let [a, p, q, expected_x, _ierr] = row[..] else {
            panic!("width");
        };
        let (x, _) = gamma_inc_inv(a, -1.0, p, q);
        assert_exact(x, expected_x, &row);
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
