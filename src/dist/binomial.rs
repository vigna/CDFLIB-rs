use crate::error::SearchError;
use crate::search::{dstinv, dstzr};
use crate::special::gamma_log;
use crate::traits::{Discrete, DiscreteCdf, Mean, Variance};
use thiserror::Error;

use super::beta::cumbet;
use super::integer_quantile;

// Parameters of cdfbin (cdflib.f90:3000-3014).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Binomial distribution with *n* trials and success probability *pr*.
///
/// Models the number of successes in a sequence of *n* independent
/// Bernoulli trials. The CDF reduces to the incomplete Β
/// (Abramowitz–Stegun 26.5.24):
/// Pr[*S* ≤ *s*] = *I*₁ ₋ *ₚ*(*n* − *s*, *s* + 1).
///
/// The methods correspond to CDFLIB's `cdfbin` (cdflib.f90:2890), whose
/// XN is *n*: `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is
/// [`inverse_ccdf`], `which = 3` is [`search_trials`], `which = 4` is
/// [`search_pr`].
///
/// # Notes
///
/// [`Entropy`] is not implemented.
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
///
/// [`Entropy`]: crate::traits::Entropy
/// [`cdf`]: DiscreteCdf::cdf
/// [`ccdf`]: DiscreteCdf::ccdf
/// [`inverse_ccdf`]: Binomial::inverse_ccdf
/// [`search_trials`]: Binomial::search_trials
/// [`search_pr`]: Binomial::search_pr
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Binomial {
    n: u64,
    pr: f64,
}

/// Errors arising from constructing a [`Binomial`] or from its parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfbin`.
///
/// [`Binomial`]: crate::Binomial
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum BinomialError {
    /// The success probability *pr* fell outside [0 . . 1] (`cdfbin` status −6);
    /// NaN is also rejected.
    #[error("success probability {0} outside [0..1]")]
    PrOutOfRange(f64),
    /// The number of trials *n* is zero (`cdfbin` status −5, *xn* ≤ 0).
    /// [`search_trials`] also returns it, checked only in Rust, when the
    /// *n* it computes is 0.
    ///
    /// [`search_trials`]: crate::Binomial::search_trials
    #[error("number of trials is zero")]
    TrialsZero,
    /// The number of successes *s* exceeds the number of trials *n*
    /// (`cdfbin` status −4).
    #[error("number of successes {s} exceeds the number of trials {n}")]
    SuccessesExceedTrials { s: u64, n: u64 },
    /// The probability *p* fell outside [0 . . 1] (`cdfbin` status −2); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfbin` status −3); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfbin` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfbin` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative binomial distribution (*cum*, *ccum*) of *s* or
/// fewer successes in *xn* trials with success probability *pr*
/// (`cumbin`, cdflib.f90:6790).
///
/// *ompr* = 1 − *pr* is passed separately to preserve its precision.
#[inline]
pub(crate) fn cumbin(s: f64, xn: f64, pr: f64, ompr: f64) -> (f64, f64) {
    if s < xn {
        let (ccum, cum) = cumbet(pr, ompr, s + 1.0, xn - s);
        (cum, ccum)
    } else {
        (1.0, 0.0)
    }
}

// cdflib.f90:3044-3063 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), BinomialError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(BinomialError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:3064-3083 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), BinomialError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(BinomialError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:3084-3096 (status -5).
#[inline]
fn check_xn(n: u64) -> Result<(), BinomialError> {
    if n == 0 {
        return Err(BinomialError::TrialsZero);
    }
    Ok(())
}

// cdflib.f90:3117-3136 (status -6). Rust also rejects NaN. OMPR = 1 - PR is
// derived from pr, so the F90 checks on ompr (status -7, cdflib.f90:3137-3156)
// and on pr + ompr (status 4, cdflib.f90:3169-3181) cannot fail.
#[inline]
fn check_pr(pr: f64) -> Result<(), BinomialError> {
    if pr < 0.0 || 1.0 < pr || pr.is_nan() {
        return Err(BinomialError::PrOutOfRange(pr));
    }
    Ok(())
}

// cdflib.f90:3157-3168 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), BinomialError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(BinomialError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl Binomial {
    /// Construct a Binomial(*n*, *pr*) distribution with *n* trials and success
    /// probability *pr* ∈ [0 . . 1].
    ///
    /// # Panics
    ///
    /// Panics if *n* is 0 or *pr* is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(n: u64, pr: f64) -> Self {
        Self::try_new(n, pr).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`BinomialError`] instead of panicking.
    ///
    /// Returns [`TrialsZero`] if *n* is zero (`cdfbin` status −5), and
    /// [`PrOutOfRange`] if *pr* falls outside [0 . . 1] or is NaN (`cdfbin`
    /// status −6).
    ///
    /// [`TrialsZero`]: BinomialError::TrialsZero
    /// [`PrOutOfRange`]: BinomialError::PrOutOfRange
    #[inline]
    pub fn try_new(n: u64, pr: f64) -> Result<Self, BinomialError> {
        check_xn(n)?;
        check_pr(pr)?;
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
    /// Pr[*S* ≤ *s*] = *p* given the success probability, searched for in
    /// [0 . . 10³⁰⁰]. The search works on the continuous extension of the CDF.
    ///
    /// CDFLIB's `cdfbin` with `which = 3`, with *ompr* = 1 − *pr*. The
    /// caller passes both *p* and *q* = 1 − *p*; they must sum to 1 within
    /// 3ε.
    ///
    /// A computed *n* of 0, at the lower end of the search interval, is
    /// reported as [`TrialsZero`].
    ///
    /// [`TrialsZero`]: BinomialError::TrialsZero
    #[inline]
    pub fn search_trials(p: f64, q: f64, pr: f64, s: u64) -> Result<f64, BinomialError> {
        check_p(p)?;
        check_q(q)?;
        // Rust only: the test s < 0 (cdflib.f90:3100-3107, status -4) is
        // dropped, since s is a u64.
        check_pr(pr)?;
        let ompr = 1.0 - pr;
        check_pq(p, q)?;
        let s = s as f64;

        // cdflib.f90:3238-3278
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let xn = d.x();
            let (cum, ccum) = cumbin(s, xn, pr, ompr);
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
        // Rust only: the F90 returns the answer whatever its value; at the
        // lower end of the search interval there are no trials.
        let xn = d.x();
        if xn == 0.0 {
            return Err(BinomialError::TrialsZero);
        }
        Ok(xn)
    }

    /// Returns the success probability *pr* satisfying Pr[*S* ≤ *s*] = *p*
    /// given *n*.
    ///
    /// CDFLIB's `cdfbin` with `which = 4`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε. When *p* > *q* the
    /// search runs on *ompr* = 1 − *pr* and returns *pr* = 1 − *ompr*.
    #[inline]
    pub fn search_pr(p: f64, q: f64, n: u64, s: u64) -> Result<f64, BinomialError> {
        check_p(p)?;
        check_q(q)?;
        check_xn(n)?;
        // cdflib.f90:3097-3116 (status -4). Rust only: the test s < 0 is
        // dropped, since s is a u64.
        if n < s {
            return Err(BinomialError::SuccessesExceedTrials { s, n });
        }
        check_pq(p, q)?;
        let xn = n as f64;
        let s = s as f64;

        // cdflib.f90:3284-3334
        let search = dstzr(0.0, 1.0, ATOL, TOL);
        let (pr, z) = if p <= q {
            let mut z = search.dzror();
            let mut ompr = 1.0 - z.x();
            while z.status() == 1 {
                let (cum, _ccum) = cumbin(s, xn, z.x(), ompr);
                let fx = cum - p;
                z.dzror(fx);
                ompr = 1.0 - z.x();
            }
            (z.x(), z)
        } else {
            let mut z = search.dzror();
            let mut pr = 1.0 - z.x();
            while z.status() == 1 {
                let (_cum, ccum) = cumbin(s, xn, pr, z.x());
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
    /// continuous extension of the CDF, searched for in [0 . . *n*].
    ///
    /// CDFLIB's `cdfbin` with `which = 2`, with *p* = 1 − *q*. CDFLIB
    /// starts the search at *s* = 5, so it fails with
    /// [`SearchError::StartOutOfRange`] when *n* < 5.
    ///
    /// [cdf]: crate::traits::DiscreteCdf::cdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, BinomialError> {
        check_q(q)?;
        let p = 1.0 - q;
        let xn = self.n as f64;
        let pr = self.pr;
        let ompr = 1.0 - pr;

        // cdflib.f90:3194-3232
        let mut d = dstinv(0.0, xn, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let s = d.x();
            let (cum, ccum) = cumbin(s, xn, pr, ompr);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: xn }
            }
            .into());
        }
        Ok(d.x())
    }
}

impl DiscreteCdf for Binomial {
    type Error = BinomialError;

    /// CDFLIB's `cdfbin` with `which = 1`, with *ompr* = 1 − *pr*.
    #[inline]
    fn cdf(&self, s: u64) -> f64 {
        // Rust only: no status -4 for n < s (cdflib.f90:3108-3114); cumbin
        // returns (1, 0) there. The test s < 0 is vacuous for a u64.
        // cdflib.f90:3187
        cumbin(s as f64, self.n as f64, self.pr, 1.0 - self.pr).0
    }

    /// CDFLIB's `cdfbin` with `which = 1`, with *ompr* = 1 − *pr*.
    #[inline]
    fn ccdf(&self, s: u64) -> f64 {
        // Rust only: no status -4 for n < s (cdflib.f90:3108-3114); cumbin
        // returns (1, 0) there. The test s < 0 is vacuous for a u64.
        // cdflib.f90:3187
        cumbin(s as f64, self.n as f64, self.pr, 1.0 - self.pr).1
    }

    /// Rust only: the smallest integer *s* ≤ *n* with
    /// [`cdf`](Self::cdf)(*s*) ≥ *p*. CDFLIB has no counterpart; its
    /// `which = 2` solves for a real *s* (see [`inverse_ccdf`]).
    ///
    /// [`inverse_ccdf`]: Binomial::inverse_ccdf
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<u64, BinomialError> {
        check_p(p)?;
        if p == 0.0 {
            return Ok(0);
        }
        if p == 1.0 {
            return Ok(self.n);
        }
        Ok(integer_quantile(p, self.n, |s| self.cdf(s)))
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
