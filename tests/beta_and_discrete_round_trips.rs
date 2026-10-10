#![cfg(not(miri))]

//! Round-trip tests for Beta, StudentsT, FisherSnedecor and the discrete
//! distributions (Binomial, Poisson, NegativeBinomial).

mod common;

use cdflib::traits::{Mean, Variance};
use cdflib::{
    Beta, Binomial, ContinuousCdf, Discrete, DiscreteCdf, FisherSnedecor, NegativeBinomial,
    Poisson, StudentsT,
};
use common::{assert_close_eps, CHAINED_INVERSE_REL_TOL, DEFAULT_REL_TOL, INVERSE_REL_TOL};

#[test]
fn beta_round_trip() {
    for &(a, b) in &[(2.0, 5.0), (0.5, 0.5), (10.0, 3.0)] {
        let d = Beta::new(a, b);
        for &p in &[0.01, 0.1, 0.5, 0.9, 0.99] {
            let x = d.inverse_cdf(p).unwrap();
            assert_close_eps(d.cdf(x), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
        }
    }
}

#[test]
fn students_t_known_quantiles() {
    // T(10) two-sided 95% critical value, computed by qt(0.975, 10)
    // in R to 15 digits. The Student's t has a particularly slow CDF
    // at this quantile (f' ≈ 0.05), so the inverse precision is
    // function-noise-limited; see CHAINED_INVERSE_REL_TOL.
    const T10_INV_975: f64 = 2.2281388519649425;
    let d = StudentsT::new(10.0);
    let x = d.inverse_cdf(0.975).unwrap();
    assert_close_eps(
        x,
        T10_INV_975,
        CHAINED_INVERSE_REL_TOL,
        CHAINED_INVERSE_REL_TOL,
    );
    // Symmetry: Pr[T ≤ 0] = 0.5 exactly.
    for &df in &[1.0, 3.0, 30.0] {
        let d = StudentsT::new(df);
        assert_eq!(d.cdf(0.0), 0.5);
    }
}

#[test]
fn students_t_round_trip() {
    for &df in &[2.0, 5.0, 30.0] {
        let d = StudentsT::new(df);
        for &p in &[0.05, 0.5, 0.95] {
            let t = d.inverse_cdf(p).unwrap();
            assert_close_eps(d.cdf(t), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
        }
    }
}

#[test]
fn students_t_extreme_quantiles_stay_finite() {
    let d = StudentsT::new(10.0);
    for &p in &[1.0e-12, 1.0 - 1.0e-12] {
        let t = d.inverse_cdf(p).unwrap();
        assert!(t.is_finite(), "p={p}, t={t}");
        assert_close_eps(d.cdf(t), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
    }
}

#[test]
fn students_t_moments_switch_at_documented_df_thresholds() {
    assert!(StudentsT::new(1.0).mean().is_nan());
    assert_eq!(StudentsT::new(1.0 + 1.0e-12).mean(), 0.0);

    assert!(StudentsT::new(1.0).variance().is_nan());
    assert!(StudentsT::new(1.5).variance().is_infinite());
    assert_eq!(StudentsT::new(2.0).variance(), f64::INFINITY);
    assert!(StudentsT::new(2.0 + 1.0e-12).variance().is_finite());
}

#[test]
fn f_distribution_round_trip() {
    let d = FisherSnedecor::new(5.0, 10.0);
    for &p in &[0.1, 0.5, 0.95] {
        let x = d.inverse_cdf(p).unwrap();
        assert_close_eps(d.cdf(x), p, INVERSE_REL_TOL, INVERSE_REL_TOL);
    }
}

#[test]
fn binomial_cdf_sums_to_one() {
    let d = Binomial::new(20, 0.3);
    let mut sum = 0.0;
    for s in 0..=20 {
        sum += d.pmf(s);
    }
    assert!((sum - 1.0).abs() < DEFAULT_REL_TOL, "sum = {sum}");
}

#[test]
fn binomial_cdf_matches_cumulative_pmf() {
    // Cumulative ∑pmf accumulates rounding linearly; cdf computes the
    // identity in a single beta_inc call. The two paths differ by O(n·ε).
    let d = Binomial::new(15, 0.4);
    let mut running = 0.0;
    for s in 0..=15 {
        running += d.pmf(s);
        assert!((d.cdf(s) - running).abs() < DEFAULT_REL_TOL, "s={s}");
    }
}

#[test]
fn poisson_cdf_matches_cumulative_pmf() {
    let d = Poisson::new(5.0);
    let mut running = 0.0;
    for s in 0..50 {
        running += d.pmf(s);
        assert!((d.cdf(s) - running).abs() < DEFAULT_REL_TOL, "s={s}");
    }
}

#[test]
fn negative_binomial_cdf_matches_cumulative_pmf() {
    let d = NegativeBinomial::new(5, 0.4);
    let mut running = 0.0;
    for s in 0..100 {
        running += d.pmf(s);
        assert!((d.cdf(s) - running).abs() < DEFAULT_REL_TOL, "s={s}");
    }
}

#[test]
fn discrete_inverse_cdf_contract() {
    let p = Poisson::new(12.0);
    for &target in &[0.1, 0.5, 0.95] {
        let s = p.inverse_cdf(target).unwrap();
        // Discrete-inverse contract: smallest s with cdf(s) >= target.
        assert!(p.cdf(s) >= target);
        if s > 0 {
            assert!(p.cdf(s - 1) < target);
        }
    }
}

#[test]
fn discrete_inverse_ccdf_contract() {
    // inverse_ccdf returns the real-valued quantile of CDFLIB's which = 2.
    // Round-trip contract: at the returned real s, the integer-floor s
    // satisfies the discrete bound ccdf(floor(s)) >= q (when interior).
    let binomial = Binomial::new(20, 0.3);
    for &target in &[0.5, 0.1] {
        let s = binomial.inverse_ccdf(target).unwrap();
        let s_floor = s.floor() as u64;
        assert!(
            binomial.ccdf(s_floor) >= target,
            "binomial: s={s}, floor={s_floor}, q={target}"
        );
    }

    let poisson = Poisson::new(2.0);
    for &target in &[0.5, 0.1, 0.05] {
        let s = poisson.inverse_ccdf(target).unwrap();
        let s_floor = s.floor() as u64;
        assert!(
            poisson.ccdf(s_floor) >= target,
            "poisson: s={s}, floor={s_floor}, q={target}"
        );
    }

    let negbin = NegativeBinomial::new(5, 0.4);
    for &target in &[0.5, 0.1] {
        let s = negbin.inverse_ccdf(target).unwrap();
        let s_floor = s.floor() as u64;
        assert!(
            negbin.ccdf(s_floor) >= target,
            "negbin: s={s}, floor={s_floor}, q={target}"
        );
    }
}
