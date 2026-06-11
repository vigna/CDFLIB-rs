// Regression tests for the boundary-input contract: inverse_cdf/inverse_sf at
// p ∈ {0, 1} return the support endpoints, search_* reject NaN/Inf with typed
// errors instead of panicking or hanging, and cdf/sf propagate NaN without
// panicking through beta_inc / gamma_inc.

#![cfg(not(miri))]

use cdflib::traits::{ContinuousCdf, DiscreteCdf};
use cdflib::{
    Beta, Binomial, ChiSquared, ChiSquaredNoncentral, FisherSnedecor, FisherSnedecorNoncentral,
    Gamma, NegativeBinomial, Normal, Poisson, StudentsT,
};

// ---- Continuous endpoint contract ----

#[test]
fn normal_endpoints() {
    let n = Normal::new(0.0, 1.0);
    assert_eq!(n.inverse_cdf(0.0).unwrap(), f64::NEG_INFINITY);
    assert_eq!(n.inverse_cdf(1.0).unwrap(), f64::INFINITY);
    assert_eq!(n.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(n.inverse_ccdf(1.0).unwrap(), f64::NEG_INFINITY);
}

#[test]
fn gamma_endpoints() {
    let g = Gamma::new(2.0, 1.5);
    assert_eq!(g.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(g.inverse_cdf(1.0).unwrap(), f64::INFINITY);
    assert_eq!(g.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(g.inverse_ccdf(1.0).unwrap(), 0.0);
}

#[test]
fn chi_squared_endpoints() {
    let c = ChiSquared::new(5.0);
    assert_eq!(c.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(c.inverse_cdf(1.0).unwrap(), f64::INFINITY);
    assert_eq!(c.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(c.inverse_ccdf(1.0).unwrap(), 0.0);
}

#[test]
fn chi_squared_noncentral_endpoints() {
    let c = ChiSquaredNoncentral::new(5.0, 2.0);
    assert_eq!(c.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(c.inverse_cdf(1.0).unwrap(), f64::INFINITY);
}

#[test]
fn beta_endpoints() {
    let b = Beta::new(2.0, 5.0);
    assert_eq!(b.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(b.inverse_cdf(1.0).unwrap(), 1.0);
    assert_eq!(b.inverse_ccdf(0.0).unwrap(), 1.0);
    assert_eq!(b.inverse_ccdf(1.0).unwrap(), 0.0);
}

#[test]
fn fisher_snedecor_endpoints() {
    let f = FisherSnedecor::new(5.0, 10.0);
    assert_eq!(f.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(f.inverse_cdf(1.0).unwrap(), f64::INFINITY);
    assert_eq!(f.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(f.inverse_ccdf(1.0).unwrap(), 0.0);
}

#[test]
fn fisher_snedecor_noncentral_endpoints() {
    let f = FisherSnedecorNoncentral::new(5.0, 10.0, 2.0);
    assert_eq!(f.inverse_cdf(0.0).unwrap(), 0.0);
    assert_eq!(f.inverse_cdf(1.0).unwrap(), f64::INFINITY);
}

#[test]
fn students_t_endpoints() {
    let t = StudentsT::new(10.0);
    assert_eq!(t.inverse_cdf(0.0).unwrap(), f64::NEG_INFINITY);
    assert_eq!(t.inverse_cdf(1.0).unwrap(), f64::INFINITY);
    assert_eq!(t.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(t.inverse_ccdf(1.0).unwrap(), f64::NEG_INFINITY);
}

// ---- Discrete endpoint contract ----

#[test]
fn binomial_endpoints() {
    let b = Binomial::new(10, 0.3);
    assert_eq!(b.inverse_cdf(0.0).unwrap(), 0);
    assert_eq!(b.inverse_cdf(1.0).unwrap(), 10);
    // inverse_sf returns the real-valued F90 cdfbin which=2 quantile.
    // At q=0 (p=1) the search converges at s=n; at q=1 (p=0) it walks
    // to the lower bound and fails per F90's status=1.
    let s = b.inverse_ccdf(0.0).unwrap();
    assert!((s - 10.0).abs() < 1e-6, "got s={s}");
    assert!(matches!(
        b.inverse_ccdf(1.0),
        Err(cdflib::BinomialError::Search(_))
    ));
}

#[test]
fn poisson_endpoints() {
    let p = Poisson::new(3.0);
    assert_eq!(p.inverse_cdf(0.0).unwrap(), 0);
    assert_eq!(p.inverse_cdf(1.0).unwrap(), u64::MAX);
    // inverse_sf returns the real-valued F90 cdfpoi which=2 quantile.
    // At q=0 the search walks to a large s where sf < abs_tol (F90 dstinv
    // converges by absolute tolerance, not by sign change); at q=1 it
    // hits the lower search bound and reports F90 status=1.
    let s_zero = p.inverse_ccdf(0.0).unwrap();
    assert!(s_zero > 10.0 && s_zero.is_finite(), "got {s_zero}");
    assert!(matches!(
        p.inverse_ccdf(1.0),
        Err(cdflib::PoissonError::Search(_))
    ));
}

#[test]
fn negative_binomial_endpoints() {
    let nb = NegativeBinomial::new(5, 0.5);
    assert_eq!(nb.inverse_cdf(0.0).unwrap(), 0);
    assert_eq!(nb.inverse_cdf(1.0).unwrap(), u64::MAX);
    // Same F90 cdfnbn which=2 behavior: q=0 converges by abs_tol at large s;
    // q=1 hits the lower search bound.
    let s_zero = nb.inverse_ccdf(0.0).unwrap();
    assert!(s_zero > 10.0 && s_zero.is_finite(), "got {s_zero}");
    assert!(matches!(
        nb.inverse_ccdf(1.0),
        Err(cdflib::NegativeBinomialError::Search(_))
    ));
}

// ---- search_* NaN rejection (must produce typed errors, not hang or panic) ----

#[test]
fn normal_solve_rejects_nan_x() {
    use cdflib::NormalError;
    assert!(matches!(
        Normal::search_mean(0.5, 0.5, f64::NAN, 1.0),
        Err(NormalError::XNotFinite(_))
    ));
    assert!(matches!(
        Normal::search_sd(0.5, 0.5, f64::NAN, 0.0),
        Err(NormalError::XNotFinite(_))
    ));
}

#[test]
fn gamma_solve_rejects_nan_x() {
    use cdflib::GammaError;
    assert!(matches!(
        Gamma::search_shape(0.5, 0.5, f64::NAN, 2.0),
        Err(GammaError::XNotFinite(_))
    ));
    assert!(matches!(
        Gamma::search_rate(0.5, 0.5, f64::NAN, 2.0),
        Err(GammaError::XNotFinite(_))
    ));
    assert!(matches!(
        Gamma::search_shape(0.5, 0.5, 1.0, f64::NAN),
        Err(GammaError::RateNotFinite(_))
    ));
}

#[test]
fn chi_squared_solve_rejects_nan_x() {
    use cdflib::ChiSquaredError;
    assert!(matches!(
        ChiSquared::search_df(0.5, 0.5, f64::NAN),
        Err(ChiSquaredError::XNotFinite(_))
    ));
}

#[test]
fn chi_squared_noncentral_solve_rejects_nan() {
    use cdflib::ChiSquaredNoncentralError;
    assert!(matches!(
        ChiSquaredNoncentral::search_df(0.5, f64::NAN, 2.0),
        Err(ChiSquaredNoncentralError::XNotFinite(_))
    ));
    assert!(matches!(
        ChiSquaredNoncentral::search_ncp(0.5, f64::NAN, 5.0),
        Err(ChiSquaredNoncentralError::XNotFinite(_))
    ));
    assert!(matches!(
        ChiSquaredNoncentral::search_df(0.5, 1.0, f64::NAN),
        Err(ChiSquaredNoncentralError::NcpNotFinite(_))
    ));
}

#[test]
fn students_t_solve_rejects_nan_t() {
    use cdflib::StudentsTError;
    assert!(matches!(
        StudentsT::search_df(0.5, 0.5, f64::NAN),
        Err(StudentsTError::TNotFinite(_))
    ));
}

#[test]
fn fisher_snedecor_noncentral_solve_rejects_nan() {
    use cdflib::FisherSnedecorNoncentralError;
    assert!(matches!(
        FisherSnedecorNoncentral::search_dfn(0.5, f64::NAN, 5.0, 1.0),
        Err(FisherSnedecorNoncentralError::FNotFinite(_))
    ));
    assert!(matches!(
        FisherSnedecorNoncentral::search_dfd(0.5, 1.0, f64::NAN, 1.0),
        Err(FisherSnedecorNoncentralError::DfnNotFinite(_))
    ));
    assert!(matches!(
        FisherSnedecorNoncentral::search_ncp(0.5, 1.0, 5.0, f64::NAN),
        Err(FisherSnedecorNoncentralError::DfdNotFinite(_))
    ));
}

// ---- cdf/sf propagate NaN (do not panic through beta_inc / gamma_inc) ----

#[test]
fn continuous_cdf_nan_returns_nan() {
    assert!(Normal::new(0.0, 1.0).cdf(f64::NAN).is_nan());
    assert!(Gamma::new(2.0, 1.0).cdf(f64::NAN).is_nan());
    assert!(ChiSquared::new(5.0).cdf(f64::NAN).is_nan());
    assert!(ChiSquaredNoncentral::new(5.0, 2.0).cdf(f64::NAN).is_nan());
    assert!(Beta::new(2.0, 5.0).cdf(f64::NAN).is_nan());
    assert!(FisherSnedecor::new(5.0, 10.0).cdf(f64::NAN).is_nan());
    assert!(FisherSnedecorNoncentral::new(5.0, 10.0, 2.0)
        .cdf(f64::NAN)
        .is_nan());
    assert!(StudentsT::new(10.0).cdf(f64::NAN).is_nan());
}

#[test]
fn continuous_ccdf_nan_returns_nan() {
    assert!(Normal::new(0.0, 1.0).ccdf(f64::NAN).is_nan());
    assert!(Gamma::new(2.0, 1.0).ccdf(f64::NAN).is_nan());
    assert!(ChiSquared::new(5.0).ccdf(f64::NAN).is_nan());
    assert!(Beta::new(2.0, 5.0).ccdf(f64::NAN).is_nan());
    assert!(StudentsT::new(10.0).ccdf(f64::NAN).is_nan());
}

// ---- F90 parity at degenerate search_pr inputs ----

#[test]
fn binomial_search_pr_all_successes_errors_instead_of_panicking() {
    use cdflib::BinomialError;
    // s == n pins cumbin to (1, 0) for every pr (cdflib.f90:6636-6645),
    // so dzror sees no sign change and reports a search failure, as the
    // F90 does with status -1 mapped through qleft/qhi.
    assert!(matches!(
        Binomial::search_pr(0.5, 0.5, 7, 7),
        Err(BinomialError::Search(_))
    ));
}

#[test]
fn negative_binomial_search_pr_r_zero_converges_like_f90() {
    // With r = 0 the cumulative is a step from 0 (pr = 0, via cumbet's
    // endpoint guard, cdflib.f90:6563-6571) to 1 (pr > 0), so the search
    // converges to the discontinuity at 0 instead of panicking inside
    // beta_inc.
    let pr = NegativeBinomial::search_pr(0.5, 0.5, 0, 5).unwrap();
    assert!(pr.abs() < 1e-7, "pr = {pr}");
}

#[test]
fn binomial_zero_trials_rejected() {
    use cdflib::BinomialError;
    // cdfbin rejects xn <= 0 with status -5 for every which except 3.
    assert!(matches!(
        Binomial::try_new(0, 0.5),
        Err(BinomialError::TrialsZero)
    ));
    assert!(matches!(
        Binomial::search_pr(0.5, 0.5, 0, 0),
        Err(BinomialError::TrialsZero)
    ));
}

#[test]
fn students_t_tail_saturates_when_t_squared_overflows() {
    // For |t| large enough that t*t overflows, cumt's beta_inc arguments
    // become (0, NaN); the x == 0 short-circuit then yields the exact
    // 0/1 tails as in the F90 cumbet path.
    let t10 = StudentsT::new(10.0);
    assert_eq!(t10.cdf(1e160), 1.0);
    assert_eq!(t10.ccdf(1e160), 0.0);
    assert_eq!(t10.cdf(-1e160), 0.0);
    assert_eq!(t10.ccdf(-1e160), 1.0);
}

#[test]
fn binomial_search_pr_all_successes_upper_tail_branch() {
    use cdflib::BinomialError;
    // Same degenerate s == n input through the p > q side: the search
    // runs on ompr and cumbin's guard pins the survival function to 0.
    assert!(matches!(
        Binomial::search_pr(0.7, 0.3, 7, 7),
        Err(BinomialError::Search(_))
    ));
}

#[test]
fn negative_binomial_search_pr_r_zero_upper_tail_branch() {
    // r = 0 through the p > q side: cumbet's endpoint guards keep the
    // ompr search away from beta_inc at both endpoints, as in the F90.
    let pr = NegativeBinomial::search_pr(0.7, 0.3, 0, 5).unwrap();
    assert!(pr.abs() < 1e-7, "pr = {pr}");
}
