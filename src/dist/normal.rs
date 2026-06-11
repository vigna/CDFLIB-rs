use crate::special::{cumnor, dinvnr};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use std::f64::consts::{E, PI};
use thiserror::Error;

/// Normal (Gaussian) distribution *N*(*μ*, *σ*²) with mean *μ* and standard
/// deviation *σ*.
///
/// # Example
///
/// ```
/// use cdflib::Normal;
/// use cdflib::traits::ContinuousCdf;
///
/// let n = Normal::new(0.0, 1.0);
///
/// // Pr[X ≤ 1.96] ≈ 0.975
/// let p = n.cdf(1.96);
///
/// // Standard normal quantile for 0.95
/// let x = n.inverse_cdf(0.95).unwrap();
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Normal {
    mean: f64,
    sd: f64,
}

/// Errors that can arise constructing a [`Normal`] or evaluating its
/// inverse routines.
///
/// [`Normal`]: crate::Normal
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum NormalError {
    /// The standard deviation *σ* was not strictly positive.
    #[error("standard deviation must be positive, got {0}")]
    SdNotPositive(f64),
    /// The mean *μ* was not finite.
    #[error("mean must be finite, got {0}")]
    MeanNotFinite(f64),
    /// The standard deviation *σ* was not finite.
    #[error("standard deviation must be finite, got {0}")]
    SdNotFinite(f64),
    /// The argument *x* was not finite.
    #[error("argument x must be finite, got {0}")]
    XNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary (|*p* + *q* − 1| > 3ε).
    /// Mirrors CDFLIB's `cdfnor` status 3 (cdflib.f90:5659).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
}

impl Normal {
    /// Construct a normal distribution with mean *μ* and standard deviation
    /// *σ* > 0.
    ///
    /// # Panics
    ///
    /// Panics if either argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(mean: f64, sd: f64) -> Self {
        Self::try_new(mean, sd).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a [`NormalError`]
    /// instead of panicking.
    ///
    /// Returns [`MeanNotFinite`], [`SdNotFinite`], or [`SdNotPositive`] if
    /// either argument fails its respective test.
    ///
    /// [`MeanNotFinite`]: NormalError::MeanNotFinite
    /// [`SdNotFinite`]: NormalError::SdNotFinite
    /// [`SdNotPositive`]: NormalError::SdNotPositive
    #[inline]
    pub fn try_new(mean: f64, sd: f64) -> Result<Self, NormalError> {
        if !mean.is_finite() {
            return Err(NormalError::MeanNotFinite(mean));
        }
        if !sd.is_finite() {
            return Err(NormalError::SdNotFinite(sd));
        }
        if sd <= 0.0 {
            return Err(NormalError::SdNotPositive(sd));
        }
        Ok(Self { mean, sd })
    }

    /// Constructs a standard normal distribution *N*(0, 1).
    #[inline]
    pub const fn standard() -> Self {
        Self { mean: 0.0, sd: 1.0 }
    }

    /// Returns the mean *μ*.
    #[inline]
    pub const fn mean(&self) -> f64 {
        self.mean
    }

    /// Returns the standard deviation *σ*.
    #[inline]
    pub const fn sd(&self) -> f64 {
        self.sd
    }

    /// Returns the mean *μ* satisfying *p* = Pr[*X* ≤ *x*] given *σ*.
    ///
    /// CDFLIB's `cdfnor` with `which = 3` (cdflib.f90:5695). Caller passes
    /// both *p* and *q* = 1 − *p*; consistency is enforced within
    /// 3ε via [`PQSumNotOne`]. Passing the pair preserves
    /// tail precision when one tail is much smaller than the other.
    ///
    /// [`PQSumNotOne`]: NormalError::PQSumNotOne
    #[inline]
    pub fn search_mean(p: f64, q: f64, x: f64, sd: f64) -> Result<f64, NormalError> {
        check_pq(p, q)?;
        if !x.is_finite() {
            return Err(NormalError::XNotFinite(x));
        }
        if !sd.is_finite() {
            return Err(NormalError::SdNotFinite(sd));
        }
        if sd <= 0.0 {
            return Err(NormalError::SdNotPositive(sd));
        }
        let z = dinvnr(p, q);
        Ok(x - sd * z)
    }

    /// Returns the standard deviation *σ* satisfying *p* = Pr[*X* ≤ *x*] given *μ*.
    ///
    /// CDFLIB's `cdfnor` with `which = 4` (cdflib.f90:5702). Caller passes
    /// both *p* and *q*; see [`search_mean`] for the (*p*, *q*) convention.
    ///
    /// The case *p* = 1/2 with *x* = *μ* is underdetermined (every *σ* > 0
    /// satisfies the equation); the formula returns a meaningless value
    /// (typically 0 since the numerator is 0 and *dinvnr* converges to a
    /// tiny non-zero denominator). F90 produces the same value.
    ///
    /// [`search_mean`]: Self::search_mean
    #[inline]
    pub fn search_sd(p: f64, q: f64, x: f64, mean: f64) -> Result<f64, NormalError> {
        check_pq(p, q)?;
        if !x.is_finite() {
            return Err(NormalError::XNotFinite(x));
        }
        if !mean.is_finite() {
            return Err(NormalError::MeanNotFinite(mean));
        }
        let z = dinvnr(p, q);
        Ok((x - mean) / z)
    }
}

#[inline]
fn check_p(p: f64) -> Result<(), NormalError> {
    if !(0.0..=1.0).contains(&p) || !p.is_finite() {
        Err(NormalError::PNotInRange(p))
    } else {
        Ok(())
    }
}

#[inline]
fn check_q(q: f64) -> Result<(), NormalError> {
    if !(0.0..=1.0).contains(&q) || !q.is_finite() {
        Err(NormalError::QNotInRange(q))
    } else {
        Ok(())
    }
}

#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), NormalError> {
    check_p(p)?;
    check_q(q)?;
    // F90 cdflib.f90:5659 uses 3 * epsilon as the consistency tolerance.
    if (p + q - 1.0).abs() > 3.0 * f64::EPSILON {
        return Err(NormalError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl ContinuousCdf for Normal {
    type Error = NormalError;

    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        let (cum, _ccum) = cumnor((x - self.mean) / self.sd);
        cum
    }

    /// Direct complementary-CDF computation, not 1 − cdf(*x*). Crucial for
    /// preserving precision in the right tail (where cdf(*x*) saturates to
    /// 1.0 well before the true value reaches it).
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        let (_cum, ccum) = cumnor((x - self.mean) / self.sd);
        ccum
    }

    /// Quantile: *x* such that Pr[*X* ≤ *x*] = *p*.
    ///
    /// Maximum precision is achieved when *p* ≤ 1/2. For *p* > 1/2, the
    /// internal *q* = 1 − *p* loses precision near *p* = 1; users with a
    /// known small right-tail probability *q* should call [`inverse_ccdf`]
    /// directly. (A single-argument API cannot carry both *p* and *q*
    /// with full precision; CDFLIB's (*p*, *q*) pair convention
    /// exists for exactly this reason.)
    ///
    /// [`inverse_ccdf`]: Self::inverse_ccdf
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, NormalError> {
        check_p(p)?;
        if p == 0.0 {
            return Ok(f64::NEG_INFINITY);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let q = 1.0 - p;
        let z = dinvnr(p, q);
        Ok(self.mean + self.sd * z)
    }
}

impl Normal {
    /// Returns the quantile *x* such that [ccdf]\(*x*\) = *q*.
    ///
    /// Mirrors CDFLIB's `cdfnor` with `which = 2`, routed through the
    /// upper-tail input so a small right-tail probability *q* keeps its
    /// precision.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, NormalError> {
        check_q(q)?;
        if q == 0.0 {
            return Ok(f64::INFINITY);
        }
        if q == 1.0 {
            return Ok(f64::NEG_INFINITY);
        }
        let p = 1.0 - q;
        let z = dinvnr(p, q);
        Ok(self.mean + self.sd * z)
    }
}

impl Continuous for Normal {
    #[inline]
    fn pdf(&self, x: f64) -> f64 {
        self.ln_pdf(x).exp()
    }

    #[inline]
    fn ln_pdf(&self, x: f64) -> f64 {
        let z = (x - self.mean) / self.sd;
        -0.5 * z * z - self.sd.ln() - 0.5 * (2.0 * PI).ln()
    }
}

impl Mean for Normal {
    #[inline]
    fn mean(&self) -> f64 {
        self.mean
    }
}

impl Variance for Normal {
    #[inline]
    fn variance(&self) -> f64 {
        self.sd * self.sd
    }
}

impl Entropy for Normal {
    /// Differential entropy: ½ ln(2π e *σ*²).
    #[inline]
    fn entropy(&self) -> f64 {
        0.5 * (2.0 * PI * E * self.sd * self.sd).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_rejects_bad_sd() {
        assert!(matches!(
            Normal::try_new(0.0, -1.0),
            Err(NormalError::SdNotPositive(_))
        ));
        assert!(matches!(
            Normal::try_new(0.0, 0.0),
            Err(NormalError::SdNotPositive(_))
        ));
        assert!(matches!(
            Normal::try_new(0.0, f64::NAN),
            Err(NormalError::SdNotFinite(_))
        ));
        assert!(matches!(
            Normal::try_new(f64::INFINITY, 1.0),
            Err(NormalError::MeanNotFinite(_))
        ));
    }

    #[test]
    fn standard_normal_at_one_sigma() {
        let n = Normal::standard();
        let p = n.cdf(1.0);
        assert!((p - 0.8413447460685429).abs() < 1e-14, "p = {p}");
    }

    #[test]
    fn ccdf_matches_1_minus_cdf_at_moderate_x() {
        let n = Normal::new(2.0, 3.0);
        for &x in &[-1.0, 0.0, 2.0, 4.0] {
            let s = (n.ccdf(x) + n.cdf(x) - 1.0).abs();
            assert!(s < 1e-14, "x = {x}: sum - 1 = {s}");
        }
    }

    #[test]
    fn ccdf_stays_accurate_in_deep_right_tail() {
        // For x = mean + 10*sd the CDF saturates to 1.0; the SF should
        // not be 0. CDFLIB-grade tail accuracy is the whole point.
        let n = Normal::new(0.0, 1.0);
        let s = n.ccdf(10.0);
        assert!(s > 0.0 && s < 1e-22, "sf(10) = {s}");
    }

    #[test]
    fn inverse_cdf_round_trip() {
        let n = Normal::new(-1.0, 2.5);
        for &x in &[-5.0, -1.0, 0.0, 3.0] {
            let p = n.cdf(x);
            let back = n.inverse_cdf(p).unwrap();
            assert!((back - x).abs() < 1e-10, "x={x}, back={back}");
        }
    }

    #[test]
    fn inverse_ccdf_handles_tiny_tails() {
        let n = Normal::standard();
        let q = n.ccdf(7.0); // ~1.28e-12
        let back = n.inverse_ccdf(q).unwrap();
        assert!((back - 7.0).abs() < 1e-9, "back = {back}");
    }

    #[test]
    fn search_mean_inverts_cdf_relation() {
        // If Pr[X ≤ 2] = 0.975 with sd = 1, mean should be 2 - 1.96 ≈ 0.04.
        let mean = Normal::search_mean(0.975, 0.025, 2.0, 1.0).unwrap();
        let expected = 2.0 - 1.9599639845400545;
        assert!((mean - expected).abs() < 1e-10, "mean = {mean}");
    }

    #[test]
    fn search_sd_inverts_cdf_relation() {
        // If Pr[X ≤ 2] = 0.975 with mean = 0, sd should be 2/1.96 ≈ 1.02.
        let sd = Normal::search_sd(0.975, 0.025, 2.0, 0.0).unwrap();
        let expected = 2.0 / 1.9599639845400545;
        assert!((sd - expected).abs() < 1e-10, "sd = {sd}");
    }

    #[test]
    fn search_sd_underdetermined_matches_f90() {
        // p = 1/2 makes z ≈ 0 and x = mean makes the numerator zero, so
        // every sd > 0 satisfies the equation. F90 returns the meaningless
        // value (x - mean) / dinvnr(0.5, 0.5) ≈ 0/tiny ≈ 0; we let that
        // propagate rather than catching it with a typed error.
        let r = Normal::search_sd(0.5, 0.5, 3.0, 3.0).unwrap();
        assert_eq!(r, 0.0, "expected the F90 underdetermined value 0; got {r}");
    }

    #[test]
    fn search_mean_rejects_bad_inputs() {
        // p out of range
        assert!(matches!(
            Normal::search_mean(-0.1, 1.1, 0.0, 1.0),
            Err(NormalError::PNotInRange(_))
        ));
        assert!(matches!(
            Normal::search_mean(1.1, -0.1, 0.0, 1.0),
            Err(NormalError::PNotInRange(_))
        ));
        // p + q != 1
        assert!(matches!(
            Normal::search_mean(0.3, 0.3, 0.0, 1.0),
            Err(NormalError::PQSumNotOne { .. })
        ));
        // sd not finite
        assert!(matches!(
            Normal::search_mean(0.5, 0.5, 0.0, f64::NAN),
            Err(NormalError::SdNotFinite(_))
        ));
        // sd <= 0
        assert!(matches!(
            Normal::search_mean(0.5, 0.5, 0.0, -1.0),
            Err(NormalError::SdNotPositive(_))
        ));
        assert!(matches!(
            Normal::search_mean(0.5, 0.5, 0.0, 0.0),
            Err(NormalError::SdNotPositive(_))
        ));
    }

    #[test]
    fn search_sd_rejects_bad_inputs() {
        assert!(matches!(
            Normal::search_sd(-0.1, 1.1, 0.0, 0.0),
            Err(NormalError::PNotInRange(_))
        ));
        assert!(matches!(
            Normal::search_sd(1.5, -0.5, 0.0, 0.0),
            Err(NormalError::PNotInRange(_))
        ));
        assert!(matches!(
            Normal::search_sd(0.5, 0.5, 0.0, f64::NAN),
            Err(NormalError::MeanNotFinite(_))
        ));
        assert!(matches!(
            Normal::search_sd(0.5, 0.5, 0.0, f64::INFINITY),
            Err(NormalError::MeanNotFinite(_))
        ));
    }

    #[test]
    fn search_mean_tail_precision_with_independent_q() {
        // The F90 (p, q) pair convention: when q is tiny and known
        // precisely, deriving q' = 1 - p loses it. search_mean should
        // use the precise q.
        let q = 1.0e-15;
        let p = 1.0 - q;
        let mean_independent = Normal::search_mean(p, q, 2.0, 1.0).unwrap();
        // For p ~ 1 - 1e-15, z ≈ +7.94, so mean ≈ 2 - 7.94 ≈ -5.94.
        assert!(mean_independent < 0.0);
    }

    #[test]
    fn pdf_at_mean_is_1_over_sd_sqrt_2pi() {
        for sd in [0.5, 1.0, 3.7] {
            let n = Normal::new(0.0, sd);
            let expected = 1.0 / (sd * (2.0 * PI).sqrt());
            let got = n.pdf(0.0);
            assert!((got - expected).abs() < 1e-15, "sd = {sd}");
        }
    }

    // Entropy compares two ln values to a 1e-15 (sub-ULP) tolerance,
    // which miri's soft-float libm cannot match. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn moments() {
        let n = Normal::new(-2.0, 3.0);
        assert_eq!(n.mean(), -2.0);
        assert_eq!(n.variance(), 9.0);
        assert_eq!(n.std_dev(), 3.0);
        let expected_entropy = 0.5 * (2.0 * PI * E * 9.0).ln();
        assert!((n.entropy() - expected_entropy).abs() < 1e-15);
    }
}
