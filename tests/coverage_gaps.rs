#![cfg(not(miri))]

//! Targeted tests that drive control flow into specific branches of the
//! special functions and distributions. Each test names the branch by its
//! triggering condition rather than by source location. Numerical
//! correctness against the F90 is the job of the reference tables, in
//! particular tests/kernel_coverage.rs and tests/dispatcher_calls.rs,
//! whose fixtures reach every reachable line of cdflib.f90 (see
//! tests/regenerate/coverage.sh).

mod common;

use cdflib::special::internal::{
    apser, beta_grat, beta_rcomp, beta_rcomp1, beta_up, fpser, gam1, gamma_rat1, rcomp, rexp,
};
use cdflib::special::{
    beta_inc, gamma, gamma_inc, gamma_inc_inv, gamma_log, psi, try_gamma_inc, try_gamma_inc_inv,
    GammaIncError, GammaIncInvError,
};
use cdflib::{ContinuousCdf, DiscreteCdf, FisherSnedecorNoncentral, NegativeBinomial, Poisson};

// ---------- gamma.rs ----------

#[test]
fn rexp_large_positive_and_negative() {
    // Drives rexp's x.exp() fallback for |x| > 0.15, both the
    // positive and negative branches.
    let pos = rexp(0.5); // exp(0.5) - 1
    assert!((pos - (0.5_f64.exp() - 1.0)).abs() < 1e-13);

    let neg = rexp(-0.5); // exp(-0.5) - 1
    assert!((neg - ((-0.5_f64).exp() - 1.0)).abs() < 1e-13);
}

#[test]
fn gam1_above_one_and_negative() {
    // gam1(1.2): the a > 1 branch (d > 0 and t > 0).
    // 1/Γ(2.2) - 1 ≈ 1/1.10180 - 1 ≈ -0.0924.
    let g_big = gam1(1.2);
    let expected_big = 1.0 / gamma(2.2) - 1.0;
    assert!((g_big - expected_big).abs() < 1e-12, "gam1(1.2) = {g_big}");

    // gam1(-0.3): the a < 0 branch (t < 0 AND d ≤ 0).
    // 1/Γ(0.7) - 1.
    let g_neg = gam1(-0.3);
    let expected_neg = 1.0 / gamma(0.7) - 1.0;
    assert!(
        (g_neg - expected_neg).abs() < 1e-12,
        "gam1(-0.3) = {g_neg}, want {expected_neg}"
    );
}

#[test]
fn gamma_negative_argument_and_overflow_paths() {
    use cdflib::special::{try_gamma, GammaDomainError};

    // |a| ≥ 15 reflection branch with t > 0.9. Γ(-15.95) routes through
    // t = 0.95 → t = 1 - 0.95 = 0.05.
    let g_reflected = gamma(-15.95);
    assert!(
        g_reflected.is_finite() && g_reflected != 0.0,
        "g_reflected = {g_reflected}"
    );

    // Reflection at a negative integer with |a| ≥ 15: sin(πt) = 0 → pole.
    assert_eq!(try_gamma(-20.0), Err(GammaDomainError::Pole(-20.0)));
    assert_eq!(try_gamma(-25.0), Err(GammaDomainError::Pole(-25.0)));

    // Same pole via the |a| < 15 branch.
    assert_eq!(try_gamma(-2.0), Err(GammaDomainError::Pole(-2.0)));
    assert_eq!(try_gamma(-5.0), Err(GammaDomainError::Pole(-5.0)));

    // Reflection-branch underflow: w, the logarithm of the magnitude of
    // Γ(-a), exceeds 0.99999 * exparg(0) for a far enough negative
    // non-integer, so Γ(a) underflows.
    assert_eq!(try_gamma(-200.7), Err(GammaDomainError::Underflow(-200.7)));

    // Large positive overflow: a ≥ 15 and 0.99999 * exparg(0) < w.
    assert_eq!(try_gamma(200.0), Err(GammaDomainError::Overflow(200.0)));
}

#[test]
fn psi_reflection_and_overflow_branches() {
    use cdflib::special::{try_psi, PsiError};

    // psi(0) is a pole.
    assert_eq!(try_psi(0.0), Err(PsiError::Pole(0.0)));

    // psi's small-|x| reflection path: aug = -1/x. Use a tiny positive x
    // below xsmall = 1e-9. CDFLIB: psi(x) ≈ -1/x - γ + O(x).
    let small = 1e-12;
    let p_small = psi(small);
    assert!(
        (p_small + 1.0 / small).abs() < 1.0,
        "psi(1e-12) = {p_small}"
    );

    // Sign-negation branch in psi's cotan reduction. psi(-0.6) routes
    // through it.
    let p462 = psi(-0.6);
    assert!(p462.is_finite());

    // Singularity at a negative integer: z = 0 inside the m+m == n
    // branch. psi(-2.0) hits this exactly (Pole variant).
    assert_eq!(try_psi(-2.0), Err(PsiError::Pole(-2.0)));
    assert_eq!(try_psi(-3.0), Err(PsiError::Pole(-3.0)));

    // w ≥ xmax1 branch in psi's reflection: huge negative argument →
    // Overflow variant.
    assert_eq!(try_psi(-1e20), Err(PsiError::Overflow(-1e20)));

    // x ≥ xmax1 in psi's asymptotic block: huge positive argument
    // bypasses the rational and returns aug + ln(x).
    let big = 1e20_f64;
    let p_big = psi(big);
    assert!((p_big - big.ln()).abs() < 1.0, "psi(1e20) = {p_big}");
}

#[test]
#[should_panic(expected = "psi(0.0): ψ has a pole at 0.0")]
fn psi_panics_on_pole() {
    let _ = psi(0.0);
}

#[test]
fn rcomp_branches() {
    // Exercise both rcomp's a < 20 and a ≥ 20 paths. Values are
    // positive; just sanity-check finite.
    //
    // a < 1 branch: a · exp(t) · (1 + gam1(a))
    let r1 = rcomp(0.5, 1.0);
    assert!(r1.is_finite() && r1 > 0.0);

    // 1 ≤ a < 20: exp(t) / Γ(a)
    let r2 = rcomp(5.0, 3.0);
    assert!(r2.is_finite() && r2 > 0.0);

    // a ≥ 20: asymptotic. Use moderate x near a so u ≠ 0.
    let r3 = rcomp(50.0, 50.0);
    assert!(r3.is_finite() && r3 > 0.0);

    // a ≥ 20 with x so tiny that x/a underflows to zero.
    // x = 0 hits the a·x == 0 short-circuit upstream, so we need
    // x > 0 but x/a == 0. f64::MIN_POSITIVE / 1e20 underflows.
    let r4 = rcomp(1e20, f64::MIN_POSITIVE);
    // Result must be exactly 0 by the early-return.
    assert_eq!(r4, 0.0);
}

#[test]
fn gamma_inc_taylor_series_through_label_200() {
    // The Taylor series for P(A,X)/X^A (label 160) for a < 1 and x < 1.1,
    // through label 200, whose clamp of a negative qans to 0 no input
    // reaches (it is listed in tests/regenerate/unreachable.txt).
    let (p, q) = gamma_inc(0.99, 1.0e-8);
    assert!(p.is_finite() && q.is_finite());
    assert!((p + q - 1.0).abs() < 1e-10);
}

#[test]
fn gamma_inc_temme_indeterminate_sentinel() {
    // The error value of the Temme expansion for L = 1 (label 330). It is
    // returned when s ≈ 0 and 3.28e-3 < a·ε², that is, when the expansion
    // cannot resolve P from Q; this needs a > 3.28e-3/ε² ≈ 6.6e28.
    // For a = x = 1e30, gamma_inc reaches label 30 with l = x/a = 1, so
    // s = 0 and label 330 returns the error.
    assert!(matches!(
        try_gamma_inc(1e30, 1e30),
        Err(GammaIncError::Indeterminate { .. })
    ));
}

// The other error values of gamma_inc (labels 270 and 410) are reached by
// tests/data/gamma_inc_edge.csv; the lines that no input reaches are listed
// with their reason in tests/regenerate/unreachable.txt.

// ---------- beta.rs ----------

#[test]
fn fpser_full_body() {
    // fpser body past the early t < exparg(1) exit. With a = 2.0 and
    // x = 0.1, t = 2·ln(0.1) ≈ -4.6 is above exparg(1), so we fall through
    // to the series.
    let eps = f64::EPSILON.max(1e-15);
    let r = fpser(2.0, 0.1, 0.1, eps);
    // I_{0.1}(2, 0.1) is small but positive; sanity-check.
    assert!(r.is_finite() && r > 0.0, "fpser(2, 0.1, 0.1) = {r}");

    // Also exercise the a ≤ 1e-3 · eps skip, with a so small the
    // exp(t) prefactor is left at 1.0.
    let r2 = fpser(1e-20, 1.0, 0.5, eps);
    assert!(r2.is_finite());
}

#[test]
fn apser_large_b_eps_else_branch() {
    // apser's else branch: b · eps > 0.02. The eps argument here is
    // treated as a tolerance, not the machine epsilon, so we can pass a
    // "large" eps directly to force the else.
    let r = apser(1e-15, 5.0, 0.05, 0.1);
    assert!(r.is_finite());
}

#[test]
fn beta_rcomp_degenerate_and_small_branches() {
    // x == 0 or y == 0 short-circuit.
    assert_eq!(beta_rcomp(2.0, 3.0, 0.0, 1.0), 0.0);
    assert_eq!(beta_rcomp(2.0, 3.0, 1.0, 0.0), 0.0);

    // x > 0.375 AND y ≤ 0.375 → alnrel(-y), y.ln() path.
    let r1 = beta_rcomp(2.0, 5.0, 0.8, 0.2);
    assert!(r1.is_finite() && r1 > 0.0);

    // a0 ≥ 1 path (both a, b ≥ 1 but min(a, b) < 8).
    let r2 = beta_rcomp(2.0, 3.0, 0.4, 0.6);
    assert!(r2.is_finite() && r2 > 0.0);

    // b0 > 1 path inside small a0.
    let r3 = beta_rcomp(0.5, 3.5, 0.3, 0.7);
    assert!(r3.is_finite() && r3 > 0.0);

    // b0 ≤ 1 path (small a0 AND small b0).
    let r4 = beta_rcomp(0.4, 0.8, 0.3, 0.7);
    assert!(r4.is_finite() && r4 > 0.0);
}

#[test]
fn gamma_rat1_branches() {
    // a · x == 0: trivial short-circuit (caller would normally avoid
    // this, but the guard exists).
    let eps = f64::EPSILON.max(1e-15);
    // a == 0 and x == 1: x > a → (1.0, 0.0).
    let (p, q) = gamma_rat1(0.0, 1.0, 1.0, eps);
    assert_eq!((p, q), (1.0, 0.0));
    // a == 1 and x == 0: x ≤ a → (0.0, 1.0).
    let (p, q) = gamma_rat1(1.0, 0.0, 1.0, eps);
    assert_eq!((p, q), (0.0, 1.0));

    // a == 1/2 AND x < 0.25 → erf branch.
    let (p, q) = gamma_rat1(0.5, 0.1, 0.0, eps);
    assert!(p.is_finite() && q.is_finite() && (p + q - 1.0).abs() < 1e-12);

    // a == 1/2 AND x ≥ 0.25 → erfc branch.
    let (p, q) = gamma_rat1(0.5, 0.5, 0.0, eps);
    assert!((p + q - 1.0).abs() < 1e-12);

    // x < 1.1 through label 50: small a and small x, so that
    // -0.13394 < z.
    let (p, q) = gamma_rat1(0.05, 0.3, 0.0, eps);
    assert!((p + q - 1.0).abs() < 1e-10);

    // Continued-fraction branch (x ≥ 1.1).
    let (p, q) = gamma_rat1(0.7, 2.0, 1.0, eps);
    assert!((p + q - 1.0).abs() < 1e-10);
}

#[test]
fn beta_rcomp1_branches() {
    // a0 ≥ 8 path (a, b ≥ 8). Use mu = 0 so the result matches
    // beta_rcomp · exp(0); keeps the assertion simple.
    let r1 = beta_rcomp1(0, 10.0, 12.0, 0.45, 0.55);
    let r1_ref = beta_rcomp(10.0, 12.0, 0.45, 0.55);
    assert!(r1.is_finite() && (r1 - r1_ref).abs() / r1_ref.abs() < 1e-12);

    // a0 ≥ 8 with the b ≤ a sub-branch.
    let r2 = beta_rcomp1(0, 15.0, 10.0, 0.6, 0.4);
    let r2_ref = beta_rcomp(15.0, 10.0, 0.6, 0.4);
    assert!(r2.is_finite() && (r2 - r2_ref).abs() / r2_ref.abs() < 1e-12);

    // |e| > 0.6 branch in u. Needs |lambda/a| > 0.6 with a0 ≥ 8.
    // lambda = a - (a+b)·x = 10 - 22·0.1 = 7.8, |7.8/10| = 0.78.
    let r_large_e_in_u = beta_rcomp1(0, 10.0, 12.0, 0.1, 0.9);
    assert!(r_large_e_in_u.is_finite() && r_large_e_in_u > 0.0);

    // |e| > 0.6 branch in v. Symmetric: a > b sub-branch makes
    // lambda = (a+b)y - b; |lambda/b| large for skewed y.
    let r_large_e_in_v = beta_rcomp1(0, 15.0, 8.0, 0.95, 0.05);
    assert!(r_large_e_in_v.is_finite() && r_large_e_in_v > 0.0);

    // a0 < 8 paths.
    let r3 = beta_rcomp1(0, 0.5, 12.0, 0.3, 0.7);
    let r3_ref = beta_rcomp(0.5, 12.0, 0.3, 0.7);
    assert!(r3.is_finite() && (r3 - r3_ref).abs() / r3_ref.abs() < 1e-12);

    // b0 ≤ 1 path with apb > 1. a = 0.5, b = 0.7 → both ≤ 1,
    // apb = 1.2 > 1, so the (1 + gam1(u))/apb branch fires.
    let r_apb_above_1 = beta_rcomp1(0, 0.5, 0.7, 0.4, 0.6);
    assert!(r_apb_above_1.is_finite() && r_apb_above_1 > 0.0);

    // b0 ≤ 1 path with esum(mu, z) underflowing to 0. mu sufficiently
    // negative drives exp(mu + z) to 0. mu = -800 puts the exponent
    // below exparg(1) ≈ -708.
    let r_underflow = beta_rcomp1(-800, 0.5, 0.5, 0.4, 0.6);
    assert_eq!(r_underflow, 0.0);
}

#[test]
fn beta_up_b_gt_1_branches() {
    // beta_up's b > 1 path with n large enough to drive the
    // "decreasing terms" inner loop.
    let eps = f64::EPSILON.max(1e-15);
    let r = beta_up(2.0, 3.0, 0.4, 0.6, 10, eps);
    assert!(r.is_finite() && r > 0.0);

    // Also a path where y ≤ 1e-4 (k = nm1 branch).
    let r2 = beta_up(2.0, 3.0, 0.99999, 1e-5, 5, eps);
    assert!(r2.is_finite());
}

#[test]
fn beta_grat_overflow_sentinel() {
    use cdflib::special::internal::BetaGratError;
    // b · z == 0.0 early-return. With b ≈ 0 and z finite, b·z is
    // exactly 0 and beta_grat returns BzZero.
    let eps = f64::EPSILON.max(1e-15);
    assert_eq!(
        beta_grat(20.0, 0.0, 0.5, 0.5, 0.0, eps),
        Err(BetaGratError::BzZero)
    );

    // Happy-path smoke test: deep-tail call doesn't panic.
    let w2 = beta_grat(100.0, 0.5, 0.5, 0.5, 0.0, eps).unwrap();
    assert!(w2.is_finite());
}

#[test]
fn beta_inc_fpser_apser_dispatch() {
    // fpser branch of beta_inc (label 90). Needs b0 < min(eps, eps·a0),
    // that is, b strictly less than about 1e-15 with a moderate.
    let (w, w1) = beta_inc(5.0, 1e-17, 0.5, 0.5);
    assert!((w + w1 - 1.0).abs() < 1e-10);

    // apser branch of beta_inc (label 100). Needs a0 < min(eps, eps·b0)
    // and b0·x0 ≤ 1.
    let (w, w1) = beta_inc(1e-17, 5.0, 0.05, 0.95);
    assert!((w + w1 - 1.0).abs() < 1e-10);
}

// ---------- distribution single-line gaps ----------

#[test]
fn poisson_inverse_cdf_high_quantile() {
    // The doubling search for the upper end of the bracket in
    // Poisson::inverse_cdf (a Rust-only integer quantile), checked for
    // consistency with cdf.
    let d = Poisson::new(4.0);
    let p = 1.0 - 1e-12;
    let s = d.inverse_cdf(p).unwrap();
    assert!(d.cdf(s) >= p);
    assert!(s > 0 && s < 100);
}

#[test]
fn negative_binomial_inverse_cdf_high_quantile() {
    // The same doubling search in NegativeBinomial::inverse_cdf.
    let d = NegativeBinomial::new(5, 0.05);
    let p = 1.0 - 1e-10;
    let s = d.inverse_cdf(p).unwrap();
    assert!(d.cdf(s) >= p);
}

#[test]
fn fisher_snedecor_noncentral_cdf_basic() {
    // The forward sum of cumfnc. Its aup - 1 + b == 0 branch is
    // unreachable, since dfn and dfd are at least 1 (see
    // tests/regenerate/unreachable.txt); exercise the surrounding code
    // with a representative input.
    let d = FisherSnedecorNoncentral::new(4.0, 8.0, 2.5);
    let x = 1.0;
    let c = d.cdf(x);
    assert!(c.is_finite() && (0.0..=1.0).contains(&c));
}

// ---------- gamma_log smoke check ----------

#[test]
fn gamma_log_does_not_regress() {
    // Canary against accidental constant changes in the asymptotic branch.
    let v = gamma_log(50.0);
    // ln Γ(50) ≈ 144.5658...
    assert!((v - 144.56574394634488).abs() < 1e-10);
}

// ---------- gamma_inc_inv extreme-value paths ----------

// The reference-table fixtures for gamma_inc_inv sample a and p on a
// moderate grid; the tests below drive branches that need genuinely extreme
// inputs (subnormals, a ≥ 10²², caller-supplied bad initial approximations).
// They are coverage-driven only: numerical correctness is tested elsewhere.

// The paths of gamma_inc_inv that no input reaches are listed with their
// reason in tests/regenerate/unreachable.txt: qg == 0 (q·Γ(a + 1) with
// q > 0 does not round to 0, since Γ(a + 1) ≥ 0.88 for a in (0..1)),
// b == 0 (b = qg/a with a < 1 only magnifies qg), r == 0 in the Schroder
// iterations (rcomp underflows only where gamma_inc already returns P or Q
// equal to 0), and x ≤ 0 after the second-order step (entry requires
// |t| ≤ 0.1 and |w·t| ≤ 0.1, so |h| ≤ 0.11 and x = xn·(1 − h) > 0).

#[test]
fn gamma_inc_inv_small_a_label_30_early_return() {
    // Label 30 with b ≤ bmin(iop) returns the c1..c5 approximation
    // directly. f64::EPSILON ≈ 2.22e-16 is not above 1e-10, so iop = 1
    // and bmin(1) = 1e-28. With a = 0.5 and q = 1e-29,
    // b = q·Γ(1.5)/0.5 ≈ 1.8e-29 < 1e-28.
    let q = 1.0e-29;
    let p = 1.0 - q;
    let (r, _) = gamma_inc_inv(0.5, -1.0, p, q);
    assert!(r.is_finite() && r > 0.0);
}

#[test]
fn gamma_inc_inv_amin_early_return() {
    // amin(iop) = 500 ≤ a (iop = 1) with |d| = |1 - xn/a| ≤ dmin(iop) =
    // 1e-6. With p = 0.5 the rational s ≈ 0, so xn ≈ a + (s² − 1)/3 ≈
    // a − 1/3 and d ≈ 1/(3a); a = 1e7 gives d ≈ 3.3e-8 < 1e-6, so the
    // routine returns xn.
    let (r, _) = gamma_inc_inv(1.0e7, -1.0, 0.5, 0.5);
    assert!(r.is_finite() && r > 0.0);
    assert!((r - 1.0e7).abs() < 1.0);
}

#[test]
fn gamma_inc_inv_label_40_bq_branch() {
    // Drives the b·q ≤ 1e-8 branch of label 40, reached when 0.6·a < qg;
    // there b·q = q·qg/a. With a = q = 1e-9: qg ≈ 1e-9, 0.6·a = 6e-10 < qg,
    // b = 1, b·q = 1e-9 ≤ 1e-8.
    let r = try_gamma_inc_inv(1.0e-9, -1.0, 1.0 - 1.0e-9, 1.0e-9);
    // The routine may legitimately return Ok or a soft-failure error;
    // the only goal here is to drive the code path.
    match r {
        Ok((x, _iters)) => assert!(x.is_finite() && x > 0.0),
        Err(GammaIncInvError::NoSolution)
        | Err(GammaIncInvError::IterationFailed)
        | Err(GammaIncInvError::UncertainAccuracy { .. }) => {}
        Err(e) => panic!("unexpected error {e:?}"),
    }
}

// Schroder iteration exits, reached with a caller-supplied x0. Each
// corresponds to an F90 ierr of -6, -7 or -8.

#[test]
fn gamma_inc_inv_schroder_p_subnormal_p() {
    // The p ≤ 1e10·tiny exit of the Schroder iteration using P (label
    // 170), reached with x0 > 0 and p ≤ 0.5.
    let p = 1.0e-300;
    let q = 1.0; // 1 - 1e-300 == 1.0 in f64
    let r = try_gamma_inc_inv(2.0, 1.0, p, q);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_q_subnormal_q() {
    // The q ≤ 1e10·tiny exit of the Schroder iteration using Q (label
    // 220), reached with x0 > 0 and 0.5 < p.
    let q = 1.0e-300;
    let p = 1.0;
    let r = try_gamma_inc_inv(2.0, 1.0, p, q);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_p_amax_certify_fail() {
    // The amax < a exit with |d| ≤ e2 of the Schroder iteration using P.
    // amax = 0.4e-10/ε² ≈ 8.1e20; a = 1e25 and x0 = a give d = 0.
    let a = 1.0e25;
    let r = try_gamma_inc_inv(a, a, 0.5, 0.5);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_q_amax_certify_fail() {
    // The same exit of the Schroder iteration using Q (0.5 < p).
    let a = 1.0e25;
    let r = try_gamma_inc_inv(a, a, 0.7, 0.3);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_p_saturates_to_zero() {
    // gamma_inc(a, x) returns (0, 1) for x far below the mode; in the
    // Schroder iteration using P (p ≤ 0.5), pn = 0 gives ierr = -8.
    let r = try_gamma_inc_inv(10.0, 1.0e-100, 0.1, 0.9);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_q_saturates_to_zero() {
    // The same on the Schroder iteration using Q: x far above the mode
    // gives qn = 0.
    let r = try_gamma_inc_inv(10.0, 1.0e10, 0.9, 0.1);
    assert!(matches!(r, Err(GammaIncInvError::UncertainAccuracy { .. })));
}

#[test]
fn gamma_inc_inv_schroder_p_first_order_negative() {
    // First-order Schröder step with t = (pn − p)/r ≥ 1 ⇒ x = xn·(1−t) ≤ 0.
    // Choose x0 well above the true x: gamma_inc(2, 10) ≈ (0.9995, 5e-4),
    // r ≈ 4.5e-3. For p = 0.01: t ≈ 0.99/4.5e-3 ≈ 220 ≫ 1.
    let r = try_gamma_inc_inv(2.0, 10.0, 0.01, 0.99);
    assert!(matches!(r, Err(GammaIncInvError::IterationFailed)));
}

#[test]
fn gamma_inc_inv_schroder_q_first_order_negative() {
    // The same on the q branch, used for p > 0.5. From x0 = 0.01, where
    // qn ≈ 1 and r ≈ 1e-4, the first step t = (q − qn)/r ≈ −1e4 overshoots
    // to x ≈ 100, where qn ≈ 4e-42 is far below q = 0.01, so the next
    // first-order step has t ≫ 1 and x = xn·(1−t) ≤ 0.
    let r = try_gamma_inc_inv(2.0, 0.01, 0.99, 0.01);
    assert!(matches!(r, Err(GammaIncInvError::IterationFailed)));
}
