use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::gamma_log;
use crate::special::GammaIncError;
use crate::traits::{Discrete, DiscreteCdf, Mean, Variance};
use thiserror::Error;

use super::chi_squared::cumchi;
use super::integer_quantile;

// Parameters of cdfpoi (cdflib.f90:5971-5983).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Poisson distribution with rate parameter *λ*.
///
/// Models the number of events occurring in a fixed interval of time or
/// space, given a known constant mean rate *λ* and independent occurrences.
/// The CDF reduces to the regularized upper incomplete Γ:
/// Pr[*X* ≤ *s*] = *Q*(*s* + 1, *λ*).
///
/// The methods correspond to CDFLIB's `cdfpoi` (cdflib.f90:5880), whose
/// XLAM is *λ*: `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is
/// [`inverse_ccdf`], `which = 3` is [`search_lambda`].
///
/// # Notes
///
/// [`Entropy`] is not implemented.
///
/// # Example
///
/// ```
/// use cdflib::Poisson;
/// use cdflib::traits::{Discrete, DiscreteCdf, Mean};
///
/// let p = Poisson::new(3.0);
/// assert_eq!(p.mean(), 3.0);
///
/// // Probability of observing exactly 2 events, 9/2 e⁻³
/// let pmf = p.pmf(2);
/// assert!((pmf - 0.22404180765538775).abs() < 1e-12);
///
/// // Probability of observing 2 or fewer events
/// let cdf = p.cdf(2);
/// assert!((cdf - 0.42319008112684353).abs() < 1e-12);
///
/// // Compute lambda given Pr[X ≤ 3] = 0.5
/// let lambda = Poisson::search_lambda(0.5, 0.5, 3).unwrap();
/// assert!((lambda - 3.672060749).abs() < 1e-6);
/// ```
///
/// [`Entropy`]: crate::traits::Entropy
/// [`cdf`]: DiscreteCdf::cdf
/// [`ccdf`]: DiscreteCdf::ccdf
/// [`inverse_ccdf`]: Poisson::inverse_ccdf
/// [`search_lambda`]: Poisson::search_lambda
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Poisson {
    lambda: f64,
}

/// Errors arising from constructing a [`Poisson`] or from its parameter search.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfpoi`.
///
/// [`Poisson`]: crate::Poisson
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum PoissonError {
    /// The rate parameter *λ* was negative (`cdfpoi` status −5). *λ* = 0,
    /// a degenerate distribution concentrated at 0, is accepted.
    #[error("lambda must be nonnegative, got {0}")]
    LambdaNegative(f64),
    /// The rate parameter *λ* was not finite, or so large that 2*λ*, the χ²
    /// argument of `cumpoi`, is not finite (checked only in Rust).
    #[error("lambda and twice lambda must be finite, got {0}")]
    LambdaNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdfpoi` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfpoi` status −3); NaN is
    /// also rejected.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfpoi` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfpoi` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative Poisson distribution (*cum*, *ccum*) of *s* or
/// fewer events with mean *xlam* (`cumpoi`, cdflib.f90:7796).
///
/// The error value of `gamma_inc`, which CDFLIB passes on unchecked, is
/// returned as a [`GammaIncError`].
///
/// [`GammaIncError`]: crate::special::GammaIncError
#[inline]
pub(crate) fn cumpoi(s: f64, xlam: f64) -> Result<(f64, f64), GammaIncError> {
    let df = 2.0 * (s + 1.0);
    let chi = 2.0 * xlam;
    let (ccum, cum) = cumchi(chi, df)?;
    Ok((cum, ccum))
}

// Rust only: cdfpoi has no status for an error value of gamma_inc, which it
// passes on as a probability. It needs s + 1 beyond 6.6e28 with xlam within
// a few ulps of it, so only the search over a real s in inverse_ccdf can
// reach it; for a u64 s, s + 1 is at most about 1.8e19.
#[inline]
fn cumpoi_or_panic(s: f64, xlam: f64) -> (f64, f64) {
    match cumpoi(s, xlam) {
        Ok(r) => r,
        Err(e) => panic!("cumpoi({s}, {xlam}): {e}"),
    }
}

// cdflib.f90:6012-6030 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), PoissonError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(PoissonError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:6032-6050 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), PoissonError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(PoissonError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:6065-6076 (status -5). Rust also rejects an xlam for which
// chi = 2 * xlam in cumpoi is not finite, which gives NaN there.
#[inline]
fn check_xlam(xlam: f64) -> Result<(), PoissonError> {
    if xlam < 0.0 {
        return Err(PoissonError::LambdaNegative(xlam));
    }
    if !(2.0 * xlam).is_finite() {
        return Err(PoissonError::LambdaNotFinite(xlam));
    }
    Ok(())
}

// cdflib.f90:6078-6088 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), PoissonError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(PoissonError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl Poisson {
    /// Construct a Poisson(*λ*) distribution with *λ* ≥ 0.
    ///
    /// # Panics
    ///
    /// Panics if *λ* is invalid; use [`try_new`] for a fallible variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(lambda: f64) -> Self {
        Self::try_new(lambda).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`PoissonError`] instead of panicking.
    ///
    /// Returns [`LambdaNegative`] if *λ* < 0, and [`LambdaNotFinite`] if *λ*
    /// is NaN or so large that 2*λ* is not finite.
    ///
    /// [`LambdaNegative`]: PoissonError::LambdaNegative
    /// [`LambdaNotFinite`]: PoissonError::LambdaNotFinite
    #[inline]
    pub fn try_new(lambda: f64) -> Result<Self, PoissonError> {
        check_xlam(lambda)?;
        Ok(Self { lambda })
    }

    /// Returns the rate parameter *λ*.
    #[inline]
    pub const fn lambda(&self) -> f64 {
        self.lambda
    }

    /// Returns the rate parameter *λ* satisfying Pr[*X* ≤ *s*] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfpoi` with `which = 3`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε.
    ///
    /// At *p* = 0, where the answer is +∞, the search stops, as the F90 does,
    /// where the computed probability becomes 0, and returns that finite value:
    /// `search_lambda(0.0, 1.0, 3)` returns about 1957.5.
    #[inline]
    pub fn search_lambda(p: f64, q: f64, s: u64) -> Result<f64, PoissonError> {
        check_p(p)?;
        check_q(q)?;
        // Rust only: the test s < 0 (cdflib.f90:6054-6063, status -4) is
        // dropped, since s is a u64.
        check_pq(p, q)?;
        let s = s as f64;

        // cdflib.f90:6145-6183
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let xlam = d.x();
            let (cum, ccum) = cumpoi_or_panic(s, xlam);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: INF }
            }
            .into());
        }
        Ok(d.x())
    }

    /// Returns the real-valued *s* such that [cdf]\(*s*\) = 1 − *q* on the
    /// continuous extension of the CDF, searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfpoi` with `which = 2`, with *p* = 1 − *q*. At *q* = 0
    /// returns +∞, also for *λ* = 0, as [`inverse_cdf`] returns
    /// [`u64::MAX`] at *p* = 1.
    ///
    /// # Panics
    ///
    /// Panics if the search evaluates `gamma_inc` where it cannot compute
    /// its result, which needs *λ* beyond 6.6 · 10²⁸ and *s* + 1 within a
    /// few ulps of it.
    ///
    /// [cdf]: crate::traits::DiscreteCdf::cdf
    /// [`inverse_cdf`]: crate::traits::DiscreteCdf::inverse_cdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, PoissonError> {
        check_q(q)?;
        // Rust only: exact endpoint. The F90 search stops at a finite s
        // where ccdf is below its absolute tolerance (or, for lambda = 0,
        // where ccdf is constant, at its start).
        if q == 0.0 {
            return Ok(f64::INFINITY);
        }
        let p = 1.0 - q;
        let xlam = self.lambda;

        // cdflib.f90:6101-6139
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let s = d.x();
            let (cum, ccum) = cumpoi_or_panic(s, xlam);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: INF }
            }
            .into());
        }
        Ok(d.x())
    }
}

impl DiscreteCdf for Poisson {
    type Error = PoissonError;

    /// CDFLIB's `cdfpoi` with `which = 1`.
    #[inline]
    fn cdf(&self, s: u64) -> f64 {
        // Rust only: the test s < 0 (cdflib.f90:6054-6063) is vacuous for a
        // u64.
        // cdflib.f90:6094
        cumpoi_or_panic(s as f64, self.lambda).0
    }

    /// CDFLIB's `cdfpoi` with `which = 1`.
    #[inline]
    fn ccdf(&self, s: u64) -> f64 {
        // Rust only: the test s < 0 (cdflib.f90:6054-6063) is vacuous for a
        // u64.
        // cdflib.f90:6094
        cumpoi_or_panic(s as f64, self.lambda).1
    }

    /// Rust only: the smallest integer *s* with [`cdf`](Self::cdf)(*s*) ≥
    /// *p* ([`u64::MAX`] if none). CDFLIB has no counterpart; its
    /// `which = 2` solves for a real *s* (see [`inverse_ccdf`]).
    ///
    /// [`inverse_ccdf`]: Poisson::inverse_ccdf
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<u64, PoissonError> {
        check_p(p)?;
        if p == 0.0 {
            return Ok(0);
        }
        if p == 1.0 {
            return Ok(u64::MAX);
        }
        Ok(integer_quantile(p, u64::MAX, |s| self.cdf(s)))
    }
}

impl Discrete for Poisson {
    #[inline]
    fn pmf(&self, s: u64) -> f64 {
        self.ln_pmf(s).exp()
    }
    #[inline]
    fn ln_pmf(&self, s: u64) -> f64 {
        if self.lambda == 0.0 {
            return if s == 0 { 0.0 } else { f64::NEG_INFINITY };
        }
        let sf = s as f64;
        sf * self.lambda.ln() - self.lambda - gamma_log(sf + 1.0)
    }
}

impl Mean for Poisson {
    #[inline]
    fn mean(&self) -> f64 {
        self.lambda
    }
}

impl Variance for Poisson {
    #[inline]
    fn variance(&self) -> f64 {
        self.lambda
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_rejects_bad_lambda() {
        assert!(matches!(
            Poisson::try_new(-1.0),
            Err(PoissonError::LambdaNegative(_))
        ));
        // λ = 0 is the degenerate point mass at 0; CDFLIB accepts it
        // (cdflib.f90:6065-6076 reject only xlam < 0), so the Rust port does too.
        assert!(Poisson::try_new(0.0).is_ok());
        assert!(matches!(
            Poisson::try_new(f64::NAN),
            Err(PoissonError::LambdaNotFinite(_))
        ));
        assert!(matches!(
            Poisson::try_new(f64::INFINITY),
            Err(PoissonError::LambdaNotFinite(_))
        ));
    }

    #[test]
    fn inverse_ccdf_matches_integer_quantile_at_integer_boundary() {
        // At the integer quantile boundary, the continuous search should
        // return the exact integer (because cumpoi(s_int, λ) for integer s
        // is the discrete CDF value at s_int).
        let p = Poisson::new(3.0);
        let q_target = p.ccdf(2); // = Pr[X > 2] for Poisson(3)
        let s = p.inverse_ccdf(q_target).unwrap();
        assert!((s - 2.0).abs() < 1e-6, "got s = {s}");
    }

    #[test]
    fn inverse_ccdf_between_integers() {
        // For a target q strictly between two integer SF values, the result
        // should be a real value strictly between the two integers.
        let p = Poisson::new(3.0);
        let hi_sf = p.ccdf(2);
        let lo_sf = p.ccdf(3);
        let q_target = 0.5 * (lo_sf + hi_sf);
        let s = p.inverse_ccdf(q_target).unwrap();
        assert!(s > 2.0 && s < 3.0, "got s = {s}");
    }

    #[test]
    fn search_lambda_rejects_bad_p() {
        assert!(matches!(
            Poisson::search_lambda(-0.1, 1.1, 3),
            Err(PoissonError::PNotInRange(_))
        ));
        assert!(matches!(
            Poisson::search_lambda(1.5, -0.5, 3),
            Err(PoissonError::PNotInRange(_))
        ));
        assert!(matches!(
            Poisson::search_lambda(f64::NAN, 0.5, 3),
            Err(PoissonError::PNotInRange(_))
        ));
    }

    #[test]
    fn inverse_cdf_p_zero_returns_zero() {
        let p = Poisson::new(5.0);
        assert_eq!(p.inverse_cdf(0.0).unwrap(), 0);
    }

    #[test]
    fn inverse_cdf_rejects_bad_p() {
        let p = Poisson::new(5.0);
        assert!(matches!(
            p.inverse_cdf(-0.1),
            Err(PoissonError::PNotInRange(_))
        ));
        assert!(matches!(
            p.inverse_cdf(1.1),
            Err(PoissonError::PNotInRange(_))
        ));
    }

    // Search convergence in this regime depends on the host FPU's exact
    // ln/exp results; miri's soft-float libm shims accumulate enough drift
    // through gamma_inc that the dinvr bracketing cannot certify a sign
    // change. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn search_lambda_uses_precision_pivot_at_both_tails() {
        // Compute lambda when p is near 1 (i.e., q near 0). The cdfpoi
        // precision pivot uses ccum-q in this regime so the residual
        // stays small. Round-trip should still recover the original.
        let lambda = 5.0_f64;
        let s = 10u64; // mean+5σ-ish, cdf will be very close to 1
        let dist = Poisson::new(lambda);
        let p_target = dist.cdf(s);
        let q_target = dist.ccdf(s);
        let recovered = Poisson::search_lambda(p_target, q_target, s).unwrap();
        assert!(
            (recovered - lambda).abs() < 1e-5,
            "p_target={p_target}, recovered={recovered}"
        );
    }

    #[test]
    fn extreme_right_tail_matches_high_precision_reference() {
        let p = Poisson::new(200.0);
        let expected_cdf = 0.999_999_993_591_493_9;
        let expected_sf = 6.408_506_071_899_014e-9;
        assert!((p.cdf(285) - expected_cdf).abs() < 1e-15);
        assert!((p.ccdf(285) - expected_sf).abs() < 1e-22);
    }

    #[test]
    fn moments_match_lambda() {
        // Verify exact-value mean/variance assertions
        let p = Poisson::new(4.0);
        assert_eq!(p.mean(), 4.0);
        assert_eq!(p.variance(), 4.0);
        assert!(p.ln_pmf(3).is_finite());
    }

    #[test]
    fn pmf_at_zero_rate() {
        let p = Poisson::new(0.0);
        assert_eq!(p.pmf(0), 1.0);
        assert_eq!(p.pmf(1), 0.0);
        assert_eq!(p.ln_pmf(1), f64::NEG_INFINITY);
    }
}
