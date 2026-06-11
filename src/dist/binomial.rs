use crate::error::SearchError;
use crate::search::{search_bounded_zero, search_monotone, SEARCH_BOUND};
use crate::special::beta_inc;
use crate::special::gamma_log;
use crate::traits::{Discrete, DiscreteCdf, Mean, Variance};
use thiserror::Error;

/// Binomial distribution with *n* trials and success probability *p*.
///
/// Models the number of successes in a sequence of *n* independent
/// Bernoulli trials. The CDF reduces to the incomplete Β
/// (Abramowitz–Stegun 26.5.24):
/// Pr[*S* ≤ *s*] = *I*₁ ₋ *ₚ*(*n* − *s*, *s* + 1).
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
/// use cdflib::Binomial;
/// use cdflib::traits::{Discrete, DiscreteCdf};
///
/// let b = Binomial::new(10, 0.3);
///
/// // Probability of 3 or fewer successes in 10 trials
/// let cdf = b.cdf(3);
///
/// // Compute success probability given Pr[S ≤ 2] = 0.5 and n = 10
/// let pr = Binomial::search_pr(0.5, 0.5, 10, 2).unwrap();
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Binomial {
    n: u64,
    pr: f64,
}

/// Errors arising from constructing a [`Binomial`] or from its parameter searches.
///
/// [`Binomial`]: crate::Binomial
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum BinomialError {
    /// The success probability *pr* fell outside [0 . . 1] (or was non-finite).
    #[error("success probability {0} outside [0..1]")]
    PrOutOfRange(f64),
    /// The number of trials *n* is zero. Mirrors CDFLIB's `cdfbin`
    /// status −5 (*xn* ≤ 0).
    #[error("number of trials is zero")]
    TrialsZero,
    /// The observed number of successes *s* exceeds the number of trials *n*.
    #[error("number of successes {s} exceeds the number of trials {n}")]
    SuccessesExceedTrials { s: u64, n: u64 },
    /// The probability *p* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary (|*p* + *q* − 1| > 3ε).
    /// Mirrors CDFLIB's `cdfbin` status 3.
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The internal root-finder failed; see [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

impl Binomial {
    /// Construct a Binomial(*n*, *pr*) distribution with *n* trials and success
    /// probability *pr* ∈ [0 . . 1].
    ///
    /// # Panics
    ///
    /// Panics if *pr* is invalid; use [`try_new`] for a fallible variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(n: u64, pr: f64) -> Self {
        Self::try_new(n, pr).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`BinomialError`] instead of panicking.
    ///
    /// Returns [`PrOutOfRange`] if *pr* falls outside [0 . . 1] or is
    /// non-finite.
    ///
    /// Returns [`TrialsZero`] if *n* is zero, mirroring CDFLIB's
    /// `cdfbin` status −5.
    ///
    /// [`TrialsZero`]: BinomialError::TrialsZero
    /// [`PrOutOfRange`]: BinomialError::PrOutOfRange
    #[inline]
    pub fn try_new(n: u64, pr: f64) -> Result<Self, BinomialError> {
        if n == 0 {
            return Err(BinomialError::TrialsZero);
        }
        if !(0.0..=1.0).contains(&pr) || !pr.is_finite() {
            return Err(BinomialError::PrOutOfRange(pr));
        }
        Ok(Self { n, pr })
    }

    /// Returns the number of trials *n*.
    #[inline]
    pub const fn n(&self) -> u64 {
        self.n
    }

    /// Returns the success probability *pr*.
    #[inline]
    pub const fn pr(&self) -> f64 {
        self.pr
    }

    /// Returns the (continuous) number of trials *n* satisfying
    /// Pr[*S* ≤ *s*] = *p* given the success probability. The search
    /// works on the continuous extension of the CDF.
    ///
    /// Mirrors CDFLIB's `cdfbin` with `which = 3`. Caller passes both
    /// *p* and *q* = 1 − *p*; consistency is enforced within 3ε.
    #[inline]
    pub fn search_trials(p: f64, q: f64, pr: f64, s: u64) -> Result<f64, BinomialError> {
        check_pq(p, q)?;
        if !(0.0..=1.0).contains(&pr) || !pr.is_finite() {
            return Err(BinomialError::PrOutOfRange(pr));
        }
        let sf = s as f64;
        // For pr fixed, Pr[S ≤ s] is decreasing in n.
        let f = |n: f64| {
            if sf >= n {
                // Degenerate left side: CDF == 1, SF == 0.
                return if p <= q { 1.0 - p } else { -q };
            }
            let (sf_bin, cdf_bin) = beta_inc(sf + 1.0, n - sf, pr, 1.0 - pr);
            if p <= q {
                cdf_bin - p
            } else {
                sf_bin - q
            }
        };
        // Match cdfbin's which=3: range (zero, inf), start = 5.0.
        Ok(search_monotone(
            0.0,
            SEARCH_BOUND,
            5.0,
            0.0,
            SEARCH_BOUND,
            f,
        )?)
    }

    /// Returns the success probability *pr* satisfying Pr[*S* ≤ *s*] = *p* given *n*.
    ///
    /// Mirrors CDFLIB's `cdfbin` with `which = 4` (cdflib.f90:3232-3265).
    /// Caller passes both *p* and *q* = 1 − *p*; consistency is enforced
    /// within 3ε. When *p* > *q* the search runs on *ompr* = 1 − *pr*
    /// (F90's variable-switch precision strategy) and returns *pr* = 1 − *ompr*.
    #[inline]
    pub fn search_pr(p: f64, q: f64, n: u64, s: u64) -> Result<f64, BinomialError> {
        check_pq(p, q)?;
        if n == 0 {
            return Err(BinomialError::TrialsZero);
        }
        if s > n {
            return Err(BinomialError::SuccessesExceedTrials { s, n });
        }
        let nf = n as f64;
        let sf = s as f64;
        // Match cdfbin's which=4 exactly: drive dzror directly on pr when
        // p<=q, else on ompr = 1-pr with the upper-tail residual.
        // The F90 loop evaluates cumbin, whose s < xn guard
        // (cdflib.f90:6636-6645) pins the cumulative to (1, 0) when
        // s >= n, keeping beta_inc away from b = 0 at the pr = 1 endpoint.
        if p <= q {
            let f = |pr: f64| {
                let cdf_bin = if s < n {
                    let (_sf_bin, cdf_bin) = beta_inc(sf + 1.0, nf - sf, pr, 1.0 - pr);
                    cdf_bin
                } else {
                    1.0
                };
                cdf_bin - p
            };
            Ok(search_bounded_zero(0.0, 1.0, f)?)
        } else {
            let f = |ompr: f64| {
                let sf_bin = if s < n {
                    let (sf_bin, _cdf_bin) = beta_inc(sf + 1.0, nf - sf, 1.0 - ompr, ompr);
                    sf_bin
                } else {
                    0.0
                };
                sf_bin - q
            };
            let ompr = search_bounded_zero(0.0, 1.0, f)?;
            Ok(1.0 - ompr)
        }
    }
}

#[inline]
fn check_p(p: f64) -> Result<(), BinomialError> {
    if !(0.0..=1.0).contains(&p) || !p.is_finite() {
        Err(BinomialError::PNotInRange(p))
    } else {
        Ok(())
    }
}

#[inline]
fn check_q(q: f64) -> Result<(), BinomialError> {
    if !(0.0..=1.0).contains(&q) || !q.is_finite() {
        Err(BinomialError::QNotInRange(q))
    } else {
        Ok(())
    }
}

#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), BinomialError> {
    check_p(p)?;
    check_q(q)?;
    if (p + q - 1.0).abs() > 3.0 * f64::EPSILON {
        return Err(BinomialError::PQSumNotOne { p, q });
    }
    Ok(())
}

/// `cumbin`: cumulative binomial via beta reduction.
fn cumbin(s: u64, n: u64, pr: f64) -> (f64, f64) {
    if s >= n {
        return (1.0, 0.0);
    }
    let sf = s as f64;
    let nf = n as f64;
    // cumbet(pr, ompr, s+1, n-s) returns (P, Q); cumbin swaps them.
    let (p, q) = beta_inc(sf + 1.0, nf - sf, pr, 1.0 - pr);
    (q, p)
}

impl DiscreteCdf for Binomial {
    type Error = BinomialError;

    #[inline]
    fn cdf(&self, s: u64) -> f64 {
        cumbin(s, self.n, self.pr).0
    }

    #[inline]
    fn ccdf(&self, s: u64) -> f64 {
        cumbin(s, self.n, self.pr).1
    }

    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<u64, BinomialError> {
        check_p(p)?;
        // Smallest s with cdf(s) >= p.
        if p == 0.0 {
            return Ok(0);
        }
        if p == 1.0 {
            return Ok(self.n);
        }
        let mut lo = 0u64;
        let mut hi = self.n;
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

impl Binomial {
    /// Returns the real-valued *s* such that [cdf]\(*s*\) = 1 − *q* on the
    /// smooth continuous extension via *I*₁₋ₚᵣ(*n*−*s*, *s*+1).
    ///
    /// Mirrors CDFLIB's `cdfbin` with `which = 2` (cdflib.f90:3138-3185):
    /// single dinvr loop, residual cum-p if p≤q else ccum-q.
    ///
    /// [cdf]: crate::traits::DiscreteCdf::cdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, BinomialError> {
        check_q(q)?;
        let nf = self.n as f64;
        let pr = self.pr;
        let p = 1.0 - q;
        // F90 cumbin(s, xn, pr, ompr, cum, ccum) handles s >= xn natively by
        // setting cum=1, ccum=0 (cdflib.f90:6636-6645); otherwise reduces to
        // beta_inc, which returns (P, Q) where Rust's binding is (ccum, cum)
        // matching F90 cumbet's argument order.
        let f = |s: f64| {
            let (cum, ccum) = if s >= nf {
                (1.0, 0.0)
            } else {
                let (cb_ccum, cb_cum) = beta_inc(s + 1.0, nf - s, pr, 1.0 - pr);
                (cb_cum, cb_ccum)
            };
            if p <= q {
                cum - p
            } else {
                ccum - q
            }
        };
        // F90 dstinv(0.0, xn, 0.5, 0.5, 5.0, atol, tol); s = 5.0.
        Ok(search_monotone(0.0, nf, 5.0, 0.0, nf, f)?)
    }
}

impl Discrete for Binomial {
    #[inline]
    fn pmf(&self, s: u64) -> f64 {
        if s > self.n {
            return 0.0;
        }
        self.ln_pmf(s).exp()
    }
    #[inline]
    fn ln_pmf(&self, s: u64) -> f64 {
        if s > self.n {
            return f64::NEG_INFINITY;
        }
        let n = self.n as f64;
        let sf = s as f64;
        let pr = self.pr;
        // ln C(n,s) + s ln pr + (n-s) ln(1-pr)
        let log_c = gamma_log(n + 1.0) - gamma_log(sf + 1.0) - gamma_log(n - sf + 1.0);
        let log_pr = if pr == 0.0 {
            if s == 0 {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        } else {
            sf * pr.ln()
        };
        let log_q = if pr == 1.0 {
            if s == self.n {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        } else {
            (n - sf) * (1.0 - pr).ln()
        };
        log_c + log_pr + log_q
    }
}

impl Mean for Binomial {
    #[inline]
    fn mean(&self) -> f64 {
        self.n as f64 * self.pr
    }
}

impl Variance for Binomial {
    #[inline]
    fn variance(&self) -> f64 {
        self.n as f64 * self.pr * (1.0 - self.pr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_inputs() {
        assert!(matches!(
            Binomial::try_new(10, -0.1),
            Err(BinomialError::PrOutOfRange(-0.1))
        ));
        assert!(matches!(
            Binomial::search_trials(-0.1, 1.1, 0.5, 3),
            Err(BinomialError::PNotInRange(-0.1))
        ));
        assert!(matches!(
            Binomial::search_trials(0.5, 0.5, f64::NAN, 3),
            Err(BinomialError::PrOutOfRange(x)) if x.is_nan()
        ));
        assert!(matches!(
            Binomial::search_pr(0.5, 0.5, 3, 4),
            Err(BinomialError::SuccessesExceedTrials { s: 4, n: 3 })
        ));
    }

    // ln_pmf(0) on a degenerate Bernoulli (pr = 0) is exactly 0.0 only
    // when the FPU evaluates the two gamma_log(11) calls bit-identically.
    // Miri's soft-float libm shims drift by ~1 ULP; skip under miri.
    #[cfg(not(miri))]
    #[test]
    fn edge_and_moment_cases() {
        let b = Binomial::new(10, 0.3);
        assert_eq!(b.inverse_cdf(0.0).unwrap(), 0);
        assert_eq!(b.pmf(11), 0.0);
        assert_eq!(b.ln_pmf(11), f64::NEG_INFINITY);
        assert_eq!(Binomial::new(10, 0.0).ln_pmf(0), 0.0);
        assert_eq!(Binomial::new(10, 1.0).ln_pmf(10), 0.0);
        assert_eq!(b.mean(), 3.0);
        assert!((b.variance() - 2.1).abs() < 1e-15);
    }
}
