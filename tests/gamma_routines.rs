#![cfg(not(miri))]

//! Reference-table tests for the gamma routines.

mod common;

use cdflib::special::internal::dstrem;
use cdflib::special::{gamma, gamma_inc, gamma_inc_inv, gamma_inc_with_acc, gamma_log, GammaIncAcc};
use common::{
    assert_close_eps, read_csv, DEFAULT_ABS_TOL, ITERATIVE_KERNEL_ABS_TOL,
    ITERATIVE_KERNEL_REL_TOL, KERNEL_REL_TOL,
};

#[test]
fn gamma_log_matches_reference() {
    for row in read_csv("tests/data/gamma_log.csv") {
        let [a, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(gamma_log(a), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

#[test]
fn gamma_matches_reference() {
    for row in read_csv("tests/data/gamma.csv") {
        let [a, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(gamma(a), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

#[test]
fn gamma_inc_matches_reference() {
    for row in read_csv("tests/data/gamma_inc.csv") {
        let [a, x, expected_p, expected_q] = row[..] else {
            panic!("width");
        };
        let (p, q) = gamma_inc(a, x);
        assert_close_eps(
            p,
            expected_p,
            ITERATIVE_KERNEL_REL_TOL,
            ITERATIVE_KERNEL_ABS_TOL,
        );
        assert_close_eps(
            q,
            expected_q,
            ITERATIVE_KERNEL_REL_TOL,
            ITERATIVE_KERNEL_ABS_TOL,
        );
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
        assert_close_eps(
            p,
            expected_p,
            ITERATIVE_KERNEL_REL_TOL,
            ITERATIVE_KERNEL_ABS_TOL,
        );
        assert_close_eps(
            q,
            expected_q,
            ITERATIVE_KERNEL_REL_TOL,
            ITERATIVE_KERNEL_ABS_TOL,
        );
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
        assert_close_eps(
            x,
            expected_x,
            ITERATIVE_KERNEL_REL_TOL,
            ITERATIVE_KERNEL_ABS_TOL,
        );
    }
}

#[test]
fn dstrem_matches_reference() {
    for row in read_csv("tests/data/dstrem.csv") {
        let [z, expected] = row[..] else {
            panic!("width");
        };
        assert_close_eps(dstrem(z), expected, KERNEL_REL_TOL, DEFAULT_ABS_TOL);
    }
}

// The a <= 0, |a| < 15 branch of gamma is libm-free (peel-loop products
// and one rational evaluation), so the Rust port matches the F90 bit for
// bit on the committed fixture. The -1.0e-8 row in particular
// distinguishes the two-step (x + 0.5) + 0.5 rounding from a single
// x + 1.0 rounding. The |a| >= 15 reflection rows go through libm and
// are covered at KERNEL_REL_TOL by gamma_matches_reference instead.
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
