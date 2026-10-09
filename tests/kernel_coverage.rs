#![cfg(not(miri))]

//! Bit-exact reference tests for the kernels of cdflib.f90 that
//! `tests/regenerate/gen_kernel_coverage.f90` drives through every branch,
//! and for the error exits of `gamma_user`, `psi`, `gamma_inc`,
//! `gamma_inc_inv` and `beta_inc`.
//!
//! The fixtures are produced with `-ffp-contract=off`, so the Rust port,
//! which never fuses multiply-adds, must reproduce every value exactly.
//! Where CDFLIB signals an error with a sentinel (`ierr`, `ans = 2`, a zero
//! result), the test checks that Rust returns the matching enum variant.

mod common;

use cdflib::special::internal::{
    algdiv, alnrel, apser, bcorr, beta_asym, beta_frac, beta_grat, beta_pser, beta_rcomp,
    beta_rcomp1, beta_up, dexpm1, esum, fpser, gam1, gamma_ln1, gamma_rat1, gsumln, rcomp, rexp,
    rlog, rlog1,
};
use cdflib::special::{
    beta, try_beta_inc, try_gamma, try_gamma_inc_inv, try_gamma_inc_with_acc, try_psi,
    BetaIncError, GammaIncAcc, GammaIncInvError,
};
use common::{assert_exact, read_csv};

#[test]
fn algdiv_matches_reference() {
    for r in read_csv("tests/data/algdiv.csv") {
        assert_exact(algdiv(r[0], r[1]), r[2], &r);
    }
}

#[test]
fn alnrel_matches_reference() {
    for r in read_csv("tests/data/alnrel.csv") {
        assert_exact(alnrel(r[0]), r[1], &r);
    }
}

#[test]
fn apser_matches_reference() {
    for r in read_csv("tests/data/apser.csv") {
        assert_exact(apser(r[0], r[1], r[2], r[3]), r[4], &r);
    }
}

#[test]
fn bcorr_matches_reference() {
    for r in read_csv("tests/data/bcorr.csv") {
        assert_exact(bcorr(r[0], r[1]), r[2], &r);
    }
}

#[test]
fn beta_matches_reference() {
    for r in read_csv("tests/data/beta.csv") {
        assert_exact(beta(r[0], r[1]), r[2], &r);
    }
}

#[test]
fn beta_asym_matches_reference() {
    for r in read_csv("tests/data/beta_asym.csv") {
        assert_exact(beta_asym(r[0], r[1], r[2], r[3]), r[4], &r);
    }
}

#[test]
fn beta_frac_matches_reference() {
    for r in read_csv("tests/data/beta_frac.csv") {
        assert_exact(beta_frac(r[0], r[1], r[2], r[3], r[4], r[5]), r[6], &r);
    }
}

#[test]
fn beta_grat_matches_reference() {
    // a, b, x, y, w, eps, w_out, ierr
    for r in read_csv("tests/data/beta_grat.csv") {
        match beta_grat(r[0], r[1], r[2], r[3], r[4], r[5]) {
            Ok(w) => {
                assert_eq!(r[7], 0.0, "{r:?}");
                assert_exact(w, r[6], &r);
            }
            Err(_) => {
                // CDFLIB ierr = 1 leaves w unchanged.
                assert_eq!(r[7], 1.0, "{r:?}");
                assert_exact(r[4], r[6], &r);
            }
        }
    }
}

#[test]
fn beta_pser_matches_reference() {
    for r in read_csv("tests/data/beta_pser.csv") {
        assert_exact(beta_pser(r[0], r[1], r[2], r[3]), r[4], &r);
    }
}

#[test]
fn beta_rcomp_matches_reference() {
    for r in read_csv("tests/data/beta_rcomp.csv") {
        assert_exact(beta_rcomp(r[0], r[1], r[2], r[3]), r[4], &r);
    }
}

#[test]
fn beta_rcomp1_matches_reference() {
    for r in read_csv("tests/data/beta_rcomp1.csv") {
        assert_exact(beta_rcomp1(r[0] as i32, r[1], r[2], r[3], r[4]), r[5], &r);
    }
}

#[test]
fn beta_up_matches_reference() {
    for r in read_csv("tests/data/beta_up.csv") {
        assert_exact(beta_up(r[0], r[1], r[2], r[3], r[4] as i32, r[5]), r[6], &r);
    }
}

#[test]
fn dexpm1_matches_reference() {
    for r in read_csv("tests/data/dexpm1.csv") {
        assert_exact(dexpm1(r[0]), r[1], &r);
    }
}

#[test]
fn esum_matches_reference() {
    for r in read_csv("tests/data/esum.csv") {
        assert_exact(esum(r[0] as i32, r[1]), r[2], &r);
    }
}

#[test]
fn fpser_matches_reference() {
    for r in read_csv("tests/data/fpser.csv") {
        assert_exact(fpser(r[0], r[1], r[2], r[3]), r[4], &r);
    }
}

#[test]
fn gam1_matches_reference() {
    for r in read_csv("tests/data/gam1.csv") {
        assert_exact(gam1(r[0]), r[1], &r);
    }
}

#[test]
fn gamma_ln1_matches_reference() {
    for r in read_csv("tests/data/gamma_ln1.csv") {
        assert_exact(gamma_ln1(r[0]), r[1], &r);
    }
}

#[test]
fn gamma_rat1_matches_reference() {
    // a, x, r, eps, p, q
    for r in read_csv("tests/data/gamma_rat1.csv") {
        let (p, q) = gamma_rat1(r[0], r[1], r[2], r[3]);
        assert_exact(p, r[4], &r);
        assert_exact(q, r[5], &r);
    }
}

#[test]
fn gsumln_matches_reference() {
    for r in read_csv("tests/data/gsumln.csv") {
        assert_exact(gsumln(r[0], r[1]), r[2], &r);
    }
}

#[test]
fn rcomp_matches_reference() {
    for r in read_csv("tests/data/rcomp.csv") {
        assert_exact(rcomp(r[0], r[1]), r[2], &r);
    }
}

#[test]
fn rexp_matches_reference() {
    for r in read_csv("tests/data/rexp.csv") {
        assert_exact(rexp(r[0]), r[1], &r);
    }
}

#[test]
fn rlog_matches_reference() {
    for r in read_csv("tests/data/rlog.csv") {
        assert_exact(rlog(r[0]), r[1], &r);
    }
}

#[test]
fn rlog1_matches_reference() {
    for r in read_csv("tests/data/rlog1.csv") {
        assert_exact(rlog1(r[0]), r[1], &r);
    }
}

#[test]
fn gamma_edge_cases_match_reference() {
    // gamma_user returns 0 when the Γ function cannot be computed; Rust
    // returns a GammaDomainError instead.
    for r in read_csv("tests/data/gamma_edge.csv") {
        match try_gamma(r[0]) {
            Ok(g) => assert_exact(g, r[1], &r),
            Err(_) => assert_eq!(r[1], 0.0, "{r:?}"),
        }
        if r[1] == 0.0 {
            assert!(try_gamma(r[0]).is_err(), "{r:?}");
        }
    }
}

#[test]
fn psi_edge_cases_match_reference() {
    // psi returns 0 on its error exits; Rust returns a PsiError instead.
    for r in read_csv("tests/data/psi_edge.csv") {
        match try_psi(r[0]) {
            Ok(v) => assert_exact(v, r[1], &r),
            Err(_) => assert_eq!(r[1], 0.0, "{r:?}"),
        }
        if r[1] == 0.0 {
            assert!(try_psi(r[0]).is_err(), "{r:?}");
        }
    }
}

#[test]
fn gamma_inc_edge_cases_match_reference() {
    // a, x, ind, ans, qans, error; error is 1 where gamma_inc sets ans = 2.
    for r in read_csv("tests/data/gamma_inc_edge.csv") {
        let acc = match r[2] as i32 {
            0 => GammaIncAcc::Max,
            1 => GammaIncAcc::Digits6,
            _ => GammaIncAcc::Digits3,
        };
        match try_gamma_inc_with_acc(r[0], r[1], acc) {
            Ok((p, q)) => {
                assert_eq!(r[5], 0.0, "{r:?}");
                assert_exact(p, r[3], &r);
                assert_exact(q, r[4], &r);
            }
            Err(_) => assert_eq!(r[5], 1.0, "{r:?}"),
        }
    }
}

#[test]
fn gamma_inc_inv_edge_cases_match_reference() {
    // a, x0, p, q, x, ierr
    for r in read_csv("tests/data/gamma_inc_inv_edge.csv") {
        let (x, ierr) = (r[4], r[5] as i32);
        let got = try_gamma_inc_inv(r[0], r[1], r[2], r[3]);
        match ierr {
            // q == 0: F90 returns x = huge(x) with ierr = 0.
            0 if x == f64::MAX => assert_eq!(got, Err(GammaIncInvError::AtInfinity), "{r:?}"),
            k if k >= 0 => {
                let (gx, iters) = got.unwrap_or_else(|e| panic!("{r:?}: {e}"));
                assert_exact(gx, x, &r);
                assert_eq!(iters as i32, k, "{r:?}");
            }
            -2 => assert!(
                matches!(got, Err(GammaIncInvError::ANotPositive(_))),
                "{r:?}"
            ),
            -3 => assert_eq!(got, Err(GammaIncInvError::NoSolution), "{r:?}"),
            -4 => assert_eq!(got, Err(GammaIncInvError::InconsistentPq), "{r:?}"),
            -6 => match got {
                Err(GammaIncInvError::NotConverged { partial }) => assert_exact(partial, x, &r),
                _ => panic!("{r:?}: {got:?}"),
            },
            -7 => assert_eq!(got, Err(GammaIncInvError::IterationFailed), "{r:?}"),
            -8 => match got {
                Err(GammaIncInvError::UncertainAccuracy { value }) => assert_exact(value, x, &r),
                _ => panic!("{r:?}: {got:?}"),
            },
            _ => panic!("unexpected ierr {ierr}"),
        }
    }
}

#[test]
fn beta_inc_regimes_match_reference() {
    // a, b, x, y, w, w1, ierr
    for r in read_csv("tests/data/beta_inc_regimes.csv") {
        let got = try_beta_inc(r[0], r[1], r[2], r[3]);
        match r[6] as i32 {
            0 => {
                let (w, w1) = got.unwrap_or_else(|e| panic!("{r:?}: {e}"));
                assert_exact(w, r[4], &r);
                assert_exact(w1, r[5], &r);
            }
            1 => assert!(
                matches!(got, Err(BetaIncError::NegativeParameter { .. })),
                "{r:?}"
            ),
            2 => assert_eq!(got, Err(BetaIncError::BothZero), "{r:?}"),
            3 => assert!(matches!(got, Err(BetaIncError::XOutOfRange(_))), "{r:?}"),
            4 => assert!(matches!(got, Err(BetaIncError::YOutOfRange(_))), "{r:?}"),
            5 => assert!(
                matches!(got, Err(BetaIncError::InconsistentSum { .. })),
                "{r:?}"
            ),
            6 => assert_eq!(got, Err(BetaIncError::XZeroAndAZero), "{r:?}"),
            7 => assert_eq!(got, Err(BetaIncError::YZeroAndBZero), "{r:?}"),
            k => panic!("unexpected ierr {k}"),
        }
    }
}
