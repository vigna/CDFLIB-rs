use crate::error::SearchError;
use crate::search::{search_bounded_zero, search_monotone, SEARCH_BOUND};
use crate::special::beta_inc;
use crate::special::gamma_log;
use crate::traits::{Discrete, DiscreteCdf, Mean, Variance};
use thiserror::Error;

/// Negative binomial distribution with target successes *r* and success
/// probability *p*.
///
/// Models the “number of failures before the *r*-th success” in a sequence
/// of independent Bernoulli trials. The CDF reduces to the incomplete Β:
/// Pr[*F* ≤ *s*] = *Iₚ*(*r*, *s* + 1).
///
/// # Notes
///
/// [`Entropy`] is not implemented.
///
/// [`Entropy`]: crate::traits::Entropy
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
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NegativeBinomial {
    r: u64,
    pr: f64,
}

/// Errors arising from constructing a [`NegativeBinomial`] or from its
/// parameter searches.
///
/// [`NegativeBinomial`]: crate::NegativeBinomial
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum NegativeBinomialError {
    /// The success probability *pr* fell outside (0 . . 1] (or was non-finite).
    /// The lower endpoint is excluded because *pr* = 0 would never produce a
    /// success.
    #[error("success probability {0} outside (0..1]")]
    PrOutOfRange(f64),
    /// The target number of successes *r* was zero.
    #[error("`r` must be positive")]
    RNotPositive,
    /// The probability *p* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary (|*p* + *q* − 1| > 3ε).
    /// Mirrors CDFLIB's `cdfnbn` status 3.
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The internal root-finder failed; see [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
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
        if r == 0 {
            return Err(NegativeBinomialError::RNotPositive);
        }
        if !(pr > 0.0 && pr <= 1.0 && pr.is_finite()) {
            return Err(NegativeBinomialError::PrOutOfRange(pr));
        }
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

    /// Returns the target number of successes *r* satisfying
    /// Pr[*F* ≤ *s*] = *p* given the success probability.
    ///
    /// Mirrors CDFLIB's `cdfnbn` with `which = 3`. Caller passes both
    /// *p* and *q* = 1 − *p*; consistency is enforced within 3ε.
    #[inline]
    pub fn search_r(p: f64, q: f64, pr: f64, s: u64) -> Result<f64, NegativeBinomialError> {
        check_pq(p, q)?;
        if !(pr > 0.0 && pr <= 1.0) {
            return Err(NegativeBinomialError::PrOutOfRange(pr));
        }
        let sf = s as f64;
        // I_pr(r, s+1) is the negative-binomial CDF (A-S 26.5.26), so
        // beta_inc's cum here is the CDF and ccum is the SF.
        let f = |r: f64| {
            let (cum, ccum) = beta_inc(r, sf + 1.0, pr, 1.0 - pr);
            if p <= q {
                cum - p
            } else {
                ccum - q
            }
        };
        // Match cdfnbn's which=3: range (0, inf), start = 5.0.
        Ok(search_monotone(
            0.0,
            SEARCH_BOUND,
            5.0,
            0.0,
            SEARCH_BOUND,
            f,
        )?)
    }

    /// Returns the success probability *pr* satisfying Pr[*F* ≤ *s*] = *p* given *r*.
    ///
    /// Mirrors CDFLIB's `cdfnbn` with `which = 4` (cdflib.f90:5400-5430).
    /// Caller passes both *p* and *q* = 1 − *p*; consistency is enforced
    /// within 3ε. When *p* > *q* the search runs on *ompr* = 1 − *pr*
    /// (F90's variable-switch precision strategy) and returns *pr* = 1 − *ompr*.
    #[inline]
    pub fn search_pr(p: f64, q: f64, r: u64, s: u64) -> Result<f64, NegativeBinomialError> {
        check_pq(p, q)?;
        let rf = r as f64;
        let sf = s as f64;
        // Match cdfnbn's which=4 exactly: drive dzror directly on pr when
        // p<=q, else on ompr = 1-pr with the upper-tail residual.
        // The F90 loop evaluates cumnbn through cumbet, whose endpoint
        // guards (cdflib.f90:6563-6571) fire before beta_inc, so pr = 0
        // yields cum = 0 even when r = 0.
        if p <= q {
            let f = |pr: f64| {
                let cum = if pr <= 0.0 {
                    0.0
                } else if 1.0 - pr <= 0.0 {
                    1.0
                } else {
                    let (cum, _ccum) = beta_inc(rf, sf + 1.0, pr, 1.0 - pr);
                    cum
                };
                cum - p
            };
            Ok(search_bounded_zero(0.0, 1.0, f)?)
        } else {
            let f = |ompr: f64| {
                let ccum = if 1.0 - ompr <= 0.0 {
                    1.0
                } else if ompr <= 0.0 {
                    0.0
                } else {
                    let (_cum, ccum) = beta_inc(rf, sf + 1.0, 1.0 - ompr, ompr);
                    ccum
                };
                ccum - q
            };
            let ompr = search_bounded_zero(0.0, 1.0, f)?;
            Ok(1.0 - ompr)
        }
    }
}

#[inline]
fn check_p(p: f64) -> Result<(), NegativeBinomialError> {
    if !(0.0..=1.0).contains(&p) || !p.is_finite() {
        Err(NegativeBinomialError::PNotInRange(p))
    } else {
        Ok(())
    }
}

#[inline]
fn check_q(q: f64) -> Result<(), NegativeBinomialError> {
    if !(0.0..=1.0).contains(&q) || !q.is_finite() {
        Err(NegativeBinomialError::QNotInRange(q))
    } else {
        Ok(())
    }
}

#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), NegativeBinomialError> {
    check_p(p)?;
    check_q(q)?;
    if (p + q - 1.0).abs() > 3.0 * f64::EPSILON {
        return Err(NegativeBinomialError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl DiscreteCdf for NegativeBinomial {
    type Error = NegativeBinomialError;

    #[inline]
    fn cdf(&self, s: u64) -> f64 {
        let (cum, _) = beta_inc(self.r as f64, s as f64 + 1.0, self.pr, 1.0 - self.pr);
        cum
    }

    #[inline]
    fn ccdf(&self, s: u64) -> f64 {
        let (_, ccum) = beta_inc(self.r as f64, s as f64 + 1.0, self.pr, 1.0 - self.pr);
        ccum
    }

    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<u64, NegativeBinomialError> {
        check_p(p)?;
        if p == 0.0 {
            return Ok(0);
        }
        if p == 1.0 {
            return Ok(u64::MAX);
        }
        let pr = self.pr;
        let r = self.r as f64;
        // Smallest s with cdf(s) >= p; sample then halve the integer range.
        let mean = r * (1.0 - pr) / pr;
        let sd = (mean / pr).sqrt();
        let mut hi = (mean + 10.0 * sd + 10.0).ceil() as u64;
        while self.cdf(hi) < p && hi < u64::MAX / 2 {
            hi *= 2;
        }
        // Saturate at u64::MAX rather than silently return a wrong answer if
        // expansion exits without finding a sign change. Unreachable for realistic
        // parameters, but preserves the "smallest s with cdf(s) ≥ p" contract.
        if self.cdf(hi) < p {
            return Ok(u64::MAX);
        }
        let mut lo = 0u64;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if self.cdf(mid) < p {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        Ok(lo)
    }
}

impl NegativeBinomial {
    /// Returns the real-valued *s* such that [cdf]\(*s*\) = 1 − *q* on the
    /// smooth continuous extension via *I*ₚᵣ(*r*, *s*+1).
    ///
    /// Mirrors CDFLIB's `cdfnbn` with `which = 2` (cdflib.f90:5335-5380):
    /// single dinvr loop, residual cum-p if p≤q else ccum-q.
    ///
    /// [cdf]: crate::traits::DiscreteCdf::cdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, NegativeBinomialError> {
        check_q(q)?;
        let rf = self.r as f64;
        let pr = self.pr;
        let p = 1.0 - q;
        // F90 cumnbn(f, s_r, pr, ompr, cum, ccum) reduces to cumbet(pr, ompr,
        // s_r, f+1.0, cum, ccum), which is beta_inc(s_r, f+1, pr, ompr) with
        // outputs in (cum, ccum) order matching beta_inc's (P, Q).
        let f = |s: f64| {
            let (cum, ccum) = beta_inc(rf, s + 1.0, pr, 1.0 - pr);
            if p <= q {
                cum - p
            } else {
                ccum - q
            }
        };
        // F90 dstinv(0.0, inf, 0.5, 0.5, 5.0, atol, tol); s = 5.0.
        Ok(search_monotone(
            0.0,
            SEARCH_BOUND,
            5.0,
            0.0,
            SEARCH_BOUND,
            f,
        )?)
    }
}

impl Discrete for NegativeBinomial {
    #[inline]
    fn pmf(&self, s: u64) -> f64 {
        self.ln_pmf(s).exp()
    }
    #[inline]
    fn ln_pmf(&self, s: u64) -> f64 {
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
}
