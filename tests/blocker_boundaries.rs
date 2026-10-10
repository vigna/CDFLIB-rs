#![cfg(not(miri))]

// Tests of the boundary-input contract: inverse_cdf and inverse_ccdf at
// p ∈ {0, 1} return the support endpoints, search_* reject NaN and infinite
// arguments with typed errors instead of panicking or hanging, and cdf and
// ccdf propagate NaN without panicking through beta_inc and gamma_inc.

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
    // inverse_ccdf returns the real-valued F90 cdfbin which=2 quantile.
    // At q=0 (p=1) it returns n, the value the F90 search converges to
    // (a row of tests/data/cdfbin_calls.csv); at q=1 (p=0) the search
    // walks to the lower bound and fails per F90's status=1.
    assert_eq!(b.inverse_ccdf(0.0).unwrap(), 10.0);
    // Rust only: n also where the F90 search stops at its start (pr = 0)
    // or cannot start (n < 5).
    assert_eq!(Binomial::new(10, 0.0).inverse_ccdf(0.0).unwrap(), 10.0);
    assert_eq!(Binomial::new(3, 0.5).inverse_ccdf(0.0).unwrap(), 3.0);
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
    // inverse_ccdf returns the real-valued F90 cdfpoi which=2 quantile.
    // At q=0 it returns +inf (Rust only: the F90 search stops at a finite
    // s where ccdf is below its absolute tolerance), also for lambda = 0;
    // at q=1 it hits the lower search bound and reports F90 status=1.
    assert_eq!(p.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(Poisson::new(0.0).inverse_ccdf(0.0).unwrap(), f64::INFINITY);
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
    // As for the Poisson distribution: +inf at q=0, also for pr = 1; q=1
    // hits the lower search bound.
    assert_eq!(nb.inverse_ccdf(0.0).unwrap(), f64::INFINITY);
    assert_eq!(
        NegativeBinomial::new(5, 1.0).inverse_ccdf(0.0).unwrap(),
        f64::INFINITY
    );
    assert!(matches!(
        nb.inverse_ccdf(1.0),
        Err(cdflib::NegativeBinomialError::Search(_))
    ));
}

// ---- search_* NaN rejection (must produce typed errors, not hang or panic) ----

#[test]
fn normal_search_rejects_nan_x() {
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
fn gamma_search_rejects_nan_x() {
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
fn chi_squared_search_rejects_nan_x() {
    use cdflib::ChiSquaredError;
    assert!(matches!(
        ChiSquared::search_df(0.5, 0.5, f64::NAN),
        Err(ChiSquaredError::XNotFinite(_))
    ));
}

#[test]
fn chi_squared_noncentral_search_rejects_nan() {
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
fn students_t_search_rejects_nan_t() {
    use cdflib::StudentsTError;
    assert!(matches!(
        StudentsT::search_df(0.5, 0.5, f64::NAN),
        Err(StudentsTError::TNotFinite(_))
    ));
}

#[test]
fn fisher_snedecor_noncentral_search_rejects_nan() {
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

// ---- cdf and ccdf propagate NaN (no panic through beta_inc or gamma_inc) ----

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
    assert!(ChiSquaredNoncentral::new(5.0, 2.0).ccdf(f64::NAN).is_nan());
    assert!(Beta::new(2.0, 5.0).ccdf(f64::NAN).is_nan());
    assert!(FisherSnedecor::new(5.0, 10.0).ccdf(f64::NAN).is_nan());
    assert!(FisherSnedecorNoncentral::new(5.0, 10.0, 2.0)
        .ccdf(f64::NAN)
        .is_nan());
    assert!(StudentsT::new(10.0).ccdf(f64::NAN).is_nan());
}

// ---- F90 parity at degenerate search_pr inputs ----

#[test]
fn binomial_search_pr_all_successes_errors_instead_of_panicking() {
    use cdflib::BinomialError;
    // s == n pins cumbin to (1, 0) for every pr (cdflib.f90:6847-6856),
    // so dzror sees no sign change and reports a search failure, as the
    // F90 does with status -1 mapped through qleft/qhi.
    assert!(matches!(
        Binomial::search_pr(0.5, 0.5, 7, 7),
        Err(BinomialError::Search(_))
    ));
}

#[test]
fn negative_binomial_search_pr_r_zero_is_rejected() {
    use cdflib::NegativeBinomialError;
    // cdfnbn accepts r = 0 and converges to the step of the cumulative at
    // pr = 0; the Rust rejects r = 0, as NegativeBinomial::try_new does.
    assert_eq!(
        NegativeBinomial::search_pr(0.5, 0.5, 0, 5),
        Err(NegativeBinomialError::RNotPositive)
    );
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
fn negative_binomial_search_pr_r_zero_upper_tail_branch_is_rejected() {
    use cdflib::NegativeBinomialError;
    // The same rejection on the p > q side, where cdfnbn searches on ompr.
    assert_eq!(
        NegativeBinomial::search_pr(0.7, 0.3, 0, 5),
        Err(NegativeBinomialError::RNotPositive)
    );
}

// ---- cdf and ccdf are exact at the infinite ends of the support ----

#[test]
fn continuous_cdf_at_infinity_is_exact() {
    let check = |name: &str, cdf: &dyn Fn(f64) -> f64, ccdf: &dyn Fn(f64) -> f64| {
        assert_eq!(cdf(f64::INFINITY), 1.0, "{name} cdf(+inf)");
        assert_eq!(ccdf(f64::INFINITY), 0.0, "{name} ccdf(+inf)");
        assert_eq!(cdf(f64::NEG_INFINITY), 0.0, "{name} cdf(-inf)");
        assert_eq!(ccdf(f64::NEG_INFINITY), 1.0, "{name} ccdf(-inf)");
    };
    let d = Normal::new(1.0, 2.0);
    check("Normal", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = Gamma::new(2.0, 3.0);
    check("Gamma", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = ChiSquared::new(3.0);
    check("ChiSquared", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = ChiSquaredNoncentral::new(3.0, 2.0);
    check("ChiSquaredNoncentral", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = FisherSnedecor::new(3.0, 5.0);
    check("FisherSnedecor", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = FisherSnedecorNoncentral::new(3.0, 5.0, 2.0);
    check("FisherSnedecorNoncentral", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = StudentsT::new(4.0);
    check("StudentsT", &|x| d.cdf(x), &|x| d.ccdf(x));
    let d = Beta::new(2.0, 3.0);
    check("Beta", &|x| d.cdf(x), &|x| d.ccdf(x));
}
