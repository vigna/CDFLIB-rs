use crate::error::SearchError;
use crate::search::{dstinv, dstzr};
use crate::special::gamma_log;
use crate::traits::{Discrete, DiscreteCdf, Mean, Variance};
use thiserror::Error;

use super::beta::cumbet;
use super::integer_quantile;

// Parameters of cdfnbn (cdflib.f90:5308-5323).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Negative binomial distribution with target successes *r* and success
/// probability *pr*.
///
/// Models the “number of failures before the *r*-th success” in a sequence
/// of independent Bernoulli trials. The CDF reduces to the incomplete Β
/// (Abramowitz–Stegun 26.5.26): Pr[*F* ≤ *s*] = *I*ₚᵣ(*r*, *s* + 1).
///
/// The methods correspond to CDFLIB's `cdfnbn` (cdflib.f90:5197), whose F
/// (the number of failures) is the argument *s* of [`cdf`] and whose S
/// (the number of successes) is *r*: `which = 1` is [`cdf`] / [`ccdf`],
/// `which = 2` is [`inverse_ccdf`], `which = 3` is [`search_r`],
/// `which = 4` is [`search_pr`].
///
/// # Notes
///
/// [`Entropy`] is not implemented.
///
/// # Example
///
/// ```
/// use cdflib::NegativeBinomial;
/// use cdflib::traits::{Discrete, DiscreteCdf};
///
/// let nb = NegativeBinomial::new(5, 0.5);
///
/// // Probability of 3 or fewer failures before 5th success
/// let cdf = nb.cdf(3);
///
/// // Compute success probability given Pr[F ≤ 5] = 0.9 and r = 10
/// let pr = NegativeBinomial::search_pr(0.9, 0.1, 10, 5).unwrap();
/// ```
///
/// [`Entropy`]: crate::traits::Entropy
/// [`cdf`]: DiscreteCdf::cdf
/// [`ccdf`]: DiscreteCdf::ccdf
/// [`inverse_ccdf`]: NegativeBinomial::inverse_ccdf
/// [`search_r`]: NegativeBinomial::search_r
/// [`search_pr`]: NegativeBinomial::search_pr
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NegativeBinomial {
    r: u64,
    pr: f64,
}

/// Errors arising from constructing a [`NegativeBinomial`] or from its
/// parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfnbn`.
///
/// [`NegativeBinomial`]: crate::NegativeBinomial
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum NegativeBinomialError {
    /// The success probability *pr* fell outside (0 . . 1] (`cdfnbn` status
    /// −6); NaN is also rejected. Rust also excludes *pr* = 0, which would
    /// never produce a success.
    #[error("success probability {0} outside (0..1]")]
    PrOutOfRange(f64),
    /// The target number of successes *r* was zero (checked only in Rust;
    /// `cdfnbn` status −5 rejects only *s* < 0).
    #[error("`r` must be positive")]
    RNotPositive,
    /// The probability *p* fell outside [0 . . 1] (`cdfnbn` status −2); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfnbn` status −3); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfnbn` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfnbn` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative negative binomial distribution (*cum*, *ccum*):
/// the probability of *f* or fewer failures before the *s*-th success,
/// each trial having success probability *pr* (`cumnbn`, cdflib.f90:7518).
///
/// *ompr* = 1 − *pr* is passed separately to preserve its precision.
#[inline]
pub(crate) fn cumnbn(f: f64, s: f64, pr: f64, ompr: f64) -> (f64, f64) {
    cumbet(pr, ompr, s, f + 1.0)
}

// cdflib.f90:5352-5371 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), NegativeBinomialError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(NegativeBinomialError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:5372-5391 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), NegativeBinomialError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(NegativeBinomialError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:5418-5437 (status -6). Rust also rejects NaN and pr = 0. OMPR =
// 1 - PR is derived from pr, so the F90 checks on ompr (status -7,
// cdflib.f90:5438-5457) and on pr + ompr (status 4, cdflib.f90:5470-5481)
// cannot fail.
#[inline]
fn check_pr(pr: f64) -> Result<(), NegativeBinomialError> {
    if pr <= 0.0 || 1.0 < pr || pr.is_nan() {
        return Err(NegativeBinomialError::PrOutOfRange(pr));
    }
    Ok(())
}

// cdflib.f90:5458-5469 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), NegativeBinomialError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(NegativeBinomialError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl NegativeBinomial {
    /// Construct a NegBin(*r*, *pr*) distribution with target successes
    /// *r* ≥ 1 and success probability *pr* ∈ (0 . . 1].
    ///
    /// # Panics
    ///
    /// Panics if either argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(r: u64, pr: f64) -> Self {
        Self::try_new(r, pr).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`NegativeBinomialError`] instead of panicking.
    #[inline]
    pub fn try_new(r: u64, pr: f64) -> Result<Self, NegativeBinomialError> {
        // Rust only: CDFLIB accepts s = 0 successes.
        if r == 0 {
            return Err(NegativeBinomialError::RNotPositive);
        }
        check_pr(pr)?;
        Ok(Self { r, pr })
    }

    /// Returns the target number of successes *r*.
    #[inline]
    pub const fn r(&self) -> u64 {
        self.r
    }

    /// Returns the success probability *pr*.
    #[inline]
    pub const fn pr(&self) -> f64 {
        self.pr
    }

    /// Returns the (continuous) target number of successes *r* satisfying
    /// Pr[*F* ≤ *s*] = *p* given the success probability, searched for in
    /// [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfnbn` with `which = 3`, with *ompr* = 1 − *pr*. The
    /// caller passes both *p* and *q* = 1 − *p*; they must sum to 1 within
    /// 3ε.
    #[inline]
    pub fn search_r(p: f64, q: f64, pr: f64, s: u64) -> Result<f64, NegativeBinomialError> {
        check_p(p)?;
        check_q(q)?;
        // Rust only: the test f < 0 (cdflib.f90:5395-5404, status -4) is
        // dropped, since s is a u64.
        check_pr(pr)?;
        let ompr = 1.0 - pr;
        check_pq(p, q)?;
        // F90 F, the number of failures.
        let f = s as f64;

        // cdflib.f90:5538-5576
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let s = d.x();
            let (cum, ccum) = cumnbn(f, s, pr, ompr);
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

    /// Returns the success probability *pr* satisfying Pr[*F* ≤ *s*] = *p*
    /// given *r*.
    ///
    /// CDFLIB's `cdfnbn` with `which = 4`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε. When *p* > *q* the
    /// search runs on *ompr* = 1 − *pr* and returns *pr* = 1 − *ompr*.
    #[inline]
    pub fn search_pr(p: f64, q: f64, r: u64, s: u64) -> Result<f64, NegativeBinomialError> {
        check_p(p)?;
        check_q(q)?;
        // Rust only: the tests f < 0 (cdflib.f90:5395-5404, status -4) and
        // s < 0 (cdflib.f90:5408-5417, status -5) are dropped, since s and r
        // are u64.
        check_pq(p, q)?;
        // F90 F (failures) and S (successes).
        let f = s as f64;
        let s = r as f64;

        // cdflib.f90:5582-5634
        let search = dstzr(0.0, 1.0, ATOL, TOL);
        let (pr, z) = if p <= q {
            let mut z = search.dzror();
            let mut ompr = 1.0 - z.x();
            while z.status() == 1 {
                let (cum, _ccum) = cumnbn(f, s, z.x(), ompr);
                let fx = cum - p;
                z.dzror(fx);
                ompr = 1.0 - z.x();
            }
            (z.x(), z)
        } else {
            let mut z = search.dzror();
            let mut pr = 1.0 - z.x();
            while z.status() == 1 {
                let (_cum, ccum) = cumnbn(f, s, pr, z.x());
                let fx = ccum - q;
                z.dzror(fx);
                pr = 1.0 - z.x();
            }
            (pr, z)
        };
        if z.status() == -1 {
            return Err(if z.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: 1.0 }
            }
            .into());
        }
        Ok(pr)
    }

    /// Returns the real-valued *s* such that [cdf]\(*s*\) = 1 − *q* on the
    /// continuous extension of the CDF, searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfnbn` with `which = 2`, with *p* = 1 − *q*.
    ///
    /// [cdf]: crate::traits::DiscreteCdf::cdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, NegativeBinomialError> {
        check_q(q)?;
        let p = 1.0 - q;
        // F90 S, the number of successes.
        let s = self.r as f64;
        let pr = self.pr;
        let ompr = 1.0 - pr;

        // cdflib.f90:5494-5532
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let f = d.x();
            let (cum, ccum) = cumnbn(f, s, pr, ompr);
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

impl DiscreteCdf for NegativeBinomial {
    type Error = NegativeBinomialError;

    /// CDFLIB's `cdfnbn` with `which = 1`, with *ompr* = 1 − *pr*.
    #[inline]
    fn cdf(&self, s: u64) -> f64 {
        // Rust only: the tests f < 0 and s < 0 (cdflib.f90:5395-5417) are
        // vacuous for u64 arguments.
        // cdflib.f90:5487
        cumnbn(s as f64, self.r as f64, self.pr, 1.0 - self.pr).0
    }

    /// CDFLIB's `cdfnbn` with `which = 1`, with *ompr* = 1 − *pr*.
    #[inline]
    fn ccdf(&self, s: u64) -> f64 {
        // Rust only: the tests f < 0 and s < 0 (cdflib.f90:5395-5417) are
        // vacuous for u64 arguments.
        // cdflib.f90:5487
        cumnbn(s as f64, self.r as f64, self.pr, 1.0 - self.pr).1
    }

    /// Rust only: the smallest integer *s* with [`cdf`](Self::cdf)(*s*) ≥
    /// *p* ([`u64::MAX`] if none). CDFLIB has no counterpart; its
    /// `which = 2` solves for a real *s* (see [`inverse_ccdf`]).
    ///
    /// [`inverse_ccdf`]: NegativeBinomial::inverse_ccdf
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<u64, NegativeBinomialError> {
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

impl Discrete for NegativeBinomial {
    #[inline]
    fn pmf(&self, s: u64) -> f64 {
        self.ln_pmf(s).exp()
    }
    #[inline]
    fn ln_pmf(&self, s: u64) -> f64 {
        if self.pr == 1.0 {
            return if s == 0 { 0.0 } else { f64::NEG_INFINITY };
        }
        let rf = self.r as f64;
        let sf = s as f64;
        // ln C(s+r-1, s) + r ln pr + s ln(1-pr)
        let log_c = gamma_log(sf + rf) - gamma_log(sf + 1.0) - gamma_log(rf);
        log_c + rf * self.pr.ln() + sf * (1.0 - self.pr).ln()
    }
}

impl Mean for NegativeBinomial {
    #[inline]
    fn mean(&self) -> f64 {
        self.r as f64 * (1.0 - self.pr) / self.pr
    }
}

impl Variance for NegativeBinomial {
    #[inline]
    fn variance(&self) -> f64 {
        let m = self.mean();
        m / self.pr
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_parameters() {
        assert!(matches!(
            NegativeBinomial::try_new(0, 0.5),
            Err(NegativeBinomialError::RNotPositive)
        ));
        assert!(matches!(
            NegativeBinomial::try_new(1, 0.0),
            Err(NegativeBinomialError::PrOutOfRange(0.0))
        ));
    }

    #[test]
    fn inverse_zero_and_moments() {
        let d = NegativeBinomial::new(5, 0.4);
        assert_eq!(d.inverse_cdf(0.0).unwrap(), 0);
        assert!(d.ln_pmf(3).is_finite());
        assert!(d.mean().is_finite());
        assert!(d.variance().is_finite());
    }

    #[test]
    fn search_helpers_reject_invalid_inputs() {
        assert!(matches!(
            NegativeBinomial::search_r(-0.1, 1.1, 0.5, 3),
            Err(NegativeBinomialError::PNotInRange(-0.1))
        ));
        assert!(matches!(
            NegativeBinomial::search_r(0.5, 0.5, 0.0, 3),
            Err(NegativeBinomialError::PrOutOfRange(0.0))
        ));
    }

    #[test]
    fn pmf_at_unit_success_probability() {
        let d = NegativeBinomial::new(3, 1.0);
        assert_eq!(d.pmf(0), 1.0);
        assert_eq!(d.pmf(1), 0.0);
        assert_eq!(d.ln_pmf(1), f64::NEG_INFINITY);
    }
}
