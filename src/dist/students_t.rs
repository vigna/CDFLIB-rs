use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::{beta_log, dt1, pow2, psi};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use thiserror::Error;

use super::beta::cumbet;

// Parameters of cdft (cdflib.f90:6286-6301).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e30;
const MAXDF: f64 = 1.0e10;
const TOL: f64 = 1.0e-8;

/// Student's *t* distribution with *df* > 0 degrees of freedom.
///
/// The methods correspond to CDFLIB's `cdft` (cdflib.f90:6189):
/// `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is [`inverse_cdf`] /
/// [`inverse_ccdf`], `which = 3` is [`search_df`].
///
/// # Example
///
/// ```
/// use cdflib::StudentsT;
/// use cdflib::traits::ContinuousCdf;
///
/// let d = StudentsT::new(10.0);
///
/// // Two-sided 95% critical value
/// let t = d.inverse_cdf(0.975).unwrap();
/// assert!((t - 2.228138852).abs() < 1e-6);
///
/// // Pr[T ≤ 2.228] ≈ 0.975
/// let p = d.cdf(2.228);
/// assert!((p - 0.9749941140914443).abs() < 1e-12);
/// ```
///
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`inverse_ccdf`]: StudentsT::inverse_ccdf
/// [`search_df`]: StudentsT::search_df
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StudentsT {
    df: f64,
}

/// Errors arising from constructing a [`StudentsT`] or from its parameter search.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdft`.
///
/// [`StudentsT`]: crate::StudentsT
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum StudentsTError {
    /// The degrees of freedom *df* was not strictly positive (`cdft`
    /// status −5).
    #[error("degrees of freedom must be positive, got {0}")]
    DfNotPositive(f64),
    /// The degrees of freedom *df* was not finite (checked only in Rust).
    #[error("degrees of freedom must be finite, got {0}")]
    DfNotFinite(f64),
    /// The argument *t* was not finite (checked only in Rust).
    #[error("argument t must be finite, got {0}")]
    TNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdft` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdft` status −3); NaN is
    /// also rejected.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdft` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdft` status 1 or 2); see
    /// [`SearchError`]. The quantile methods also return
    /// [`StartOutOfRange`] when the starting value given by `dt1` falls
    /// outside the search interval, or is NaN, as for a tiny *df*, where
    /// CDFLIB stops with a fatal error.
    ///
    /// [`SearchError`]: crate::error::SearchError
    /// [`StartOutOfRange`]: crate::error::SearchError::StartOutOfRange
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative *t* distribution (*cum*, *ccum*) at *t* with
/// *df* degrees of freedom (`cumt`, cdflib.f90:7855).
#[inline]
pub(crate) fn cumt(t: f64, df: f64) -> (f64, f64) {
    let xx = df / (df + pow2(t));
    let yy = pow2(t) / (df + pow2(t));
    let (a, oma) = cumbet(xx, yy, 0.5 * df, 0.5);
    if t <= 0.0 {
        let cum = 0.5 * a;
        let ccum = oma + cum;
        (cum, ccum)
    } else {
        let ccum = 0.5 * a;
        let cum = oma + ccum;
        (cum, ccum)
    }
}

// cdflib.f90:6328-6347 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), StudentsTError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(StudentsTError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:6348-6367 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), StudentsTError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(StudentsTError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:6368-6380 (status -5). Rust also rejects a non-finite df.
#[inline]
fn check_df(df: f64) -> Result<(), StudentsTError> {
    if df <= 0.0 {
        return Err(StudentsTError::DfNotPositive(df));
    }
    if !df.is_finite() {
        return Err(StudentsTError::DfNotFinite(df));
    }
    Ok(())
}

// cdflib.f90:6381-6393 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), StudentsTError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(StudentsTError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl StudentsT {
    /// Construct a Student's *t* distribution with *df* > 0 degrees of
    /// freedom.
    ///
    /// # Panics
    ///
    /// Panics if *df* is invalid; use [`try_new`] for a fallible variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(df: f64) -> Self {
        Self::try_new(df).unwrap()
    }

    /// Fallible counterpart of [`new`] returning a
    /// [`StudentsTError`] instead of panicking.
    ///
    /// Returns [`DfNotPositive`] or [`DfNotFinite`] otherwise.
    ///
    /// [`DfNotPositive`]: StudentsTError::DfNotPositive
    /// [`DfNotFinite`]: StudentsTError::DfNotFinite
    /// [`new`]: Self::new
    #[inline]
    pub fn try_new(df: f64) -> Result<Self, StudentsTError> {
        check_df(df)?;
        Ok(Self { df })
    }

    /// Returns the degrees of freedom *df*.
    #[inline]
    pub const fn df(&self) -> f64 {
        self.df
    }

    /// Returns the degrees of freedom *df* satisfying Pr[*T* ≤ *t*] = *p*,
    /// searched for in [1 . . 10¹⁰].
    ///
    /// CDFLIB's `cdft` with `which = 3`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε. As in CDFLIB, the
    /// lower bound reported on failure is 0, not 1.
    #[inline]
    pub fn search_df(p: f64, q: f64, t: f64) -> Result<f64, StudentsTError> {
        check_p(p)?;
        check_q(q)?;
        check_pq(p, q)?;
        // Rust only: a finite t.
        if !t.is_finite() {
            return Err(StudentsTError::TNotFinite(t));
        }

        // cdflib.f90:6450-6488
        let mut d = dstinv(1.0, MAXDF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let df = d.x();
            let (cum, ccum) = cumt(t, df);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: MAXDF }
            }
            .into());
        }
        Ok(d.x())
    }

    /// CDFLIB's `cdft` with `which = 2`: returns *t* given (*p*, *q*),
    /// already checked.
    fn search_t(&self, p: f64, q: f64) -> Result<f64, StudentsTError> {
        let df = self.df;
        // cdflib.f90:6406-6444
        let mut d = dstinv(-INF, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(dt1(p, q, df))?;
        while d.status() == 1 {
            let t = d.x();
            let (cum, ccum) = cumt(t, df);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: -INF }
            } else {
                SearchError::AnswerAboveUpperBound { bound: INF }
            }
            .into());
        }
        Ok(d.x())
    }

    /// Returns the quantile *t* such that [ccdf]\(*t*\) = *q*, searched for
    /// in [−10³⁰ . . 10³⁰] starting from [`dt1`].
    ///
    /// CDFLIB's `cdft` with `which = 2`, with *p* = 1 − *q*.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    /// [`dt1`]: crate::special::dt1
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, StudentsTError> {
        check_q(q)?;
        // Rust only: exact endpoints.
        if q == 0.0 {
            return Ok(f64::INFINITY);
        }
        if q == 1.0 {
            return Ok(f64::NEG_INFINITY);
        }
        let p = 1.0 - q;
        self.search_t(p, q)
    }
}

impl ContinuousCdf for StudentsT {
    type Error = StudentsTError;

    /// CDFLIB's `cdft` with `which = 1`.
    #[inline]
    fn cdf(&self, t: f64) -> f64 {
        // Rust only: NaN for a NaN t, which cumt passes to beta_inc; as in
        // the F90, beta_inc can then return a value computed from the other
        // arguments.
        if t.is_nan() {
            return f64::NAN;
        }
        // cdflib.f90:6399
        cumt(t, self.df).0
    }

    /// CDFLIB's `cdft` with `which = 1`.
    #[inline]
    fn ccdf(&self, t: f64) -> f64 {
        // Rust only: NaN for a NaN t, as in cdf.
        if t.is_nan() {
            return f64::NAN;
        }
        // cdflib.f90:6399
        cumt(t, self.df).1
    }

    /// CDFLIB's `cdft` with `which = 2`, with *q* = 1 − *p*.
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, StudentsTError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(f64::NEG_INFINITY);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let q = 1.0 - p;
        self.search_t(p, q)
    }
}

impl Continuous for StudentsT {
    #[inline]
    fn pdf(&self, t: f64) -> f64 {
        self.ln_pdf(t).exp()
    }
    #[inline]
    fn ln_pdf(&self, t: f64) -> f64 {
        let df = self.df;
        // ln f(t) = -ln(√df · Β(df/2, 1/2)) - (df + 1)/2 · ln(1 + t²/df).
        // For large df, beta_log and ln_1p keep the precision that
        // ln Γ((df + 1)/2) - ln Γ(df/2) and ln(1 + t²/df) would lose. Where
        // t * t overflows, ln(1 + t²/df) = ln((df + t²) / df) is computed as
        // 2 ln hypot(√df, t) - ln df, which cannot overflow; there t²/df > 1,
        // so the difference does not cancel. This includes t = ±inf.
        let r = t * t / df;
        let ln_term = if r.is_finite() {
            r.ln_1p()
        } else {
            2.0 * df.sqrt().hypot(t).ln() - df.ln()
        };
        -0.5 * df.ln() - beta_log(0.5 * df, 0.5) - 0.5 * (df + 1.0) * ln_term
    }
}

impl Mean for StudentsT {
    /// Defined only for *df* > 1; we return 0 for *df* > 1 and NaN
    /// for *df* ≤ 1.
    #[inline]
    fn mean(&self) -> f64 {
        if self.df > 1.0 {
            0.0
        } else {
            f64::NAN
        }
    }
}

impl Variance for StudentsT {
    /// Defined as *df*/(*df* − 2) for *df* > 2, ∞ for 1 < *df* ≤ 2, NaN
    /// otherwise.
    #[inline]
    fn variance(&self) -> f64 {
        if self.df > 2.0 {
            self.df / (self.df - 2.0)
        } else if self.df > 1.0 {
            f64::INFINITY
        } else {
            f64::NAN
        }
    }
}

impl Entropy for StudentsT {
    #[inline]
    fn entropy(&self) -> f64 {
        let df = self.df;
        // For the smallest subnormal df, df/2 is 0, where ψ has a pole; the
        // entropy tends to +inf as df tends to 0.
        if df / 2.0 == 0.0 {
            return f64::INFINITY;
        }
        // H = (df+1)/2 · [ψ((df+1)/2) - ψ(df/2)] + ln(√df · Β(df/2, 1/2))
        // = (df+1)/2 · [ψ((df+1)/2) - ψ(df/2)] + 0.5·ln(df) + ln Β(df/2, 1/2)
        0.5 * (df + 1.0) * (psi((df + 1.0) / 2.0) - psi(df / 2.0))
            + 0.5 * df.ln()
            + beta_log(df / 2.0, 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_rejects_nonpositive_or_nonfinite_df() {
        assert!(matches!(
            StudentsT::try_new(0.0),
            Err(StudentsTError::DfNotPositive(0.0))
        ));
        assert!(matches!(
            StudentsT::try_new(-1.0),
            Err(StudentsTError::DfNotPositive(-1.0))
        ));
        assert!(matches!(
            StudentsT::try_new(f64::INFINITY),
            Err(StudentsTError::DfNotFinite(x)) if x.is_infinite()
        ));
        assert!(matches!(
            StudentsT::try_new(f64::NAN),
            Err(StudentsTError::DfNotFinite(_))
        ));
    }

    #[test]
    fn rejects_probability_out_of_range() {
        let d = StudentsT::new(10.0);
        assert!(matches!(
            d.inverse_cdf(-1.0),
            Err(StudentsTError::PNotInRange(-1.0))
        ));
        assert!(matches!(
            d.inverse_ccdf(2.0),
            Err(StudentsTError::QNotInRange(2.0))
        ));
    }

    #[test]
    fn inverse_ccdf_is_zero_at_median() {
        let d = StudentsT::new(7.0);
        assert!(d.inverse_ccdf(0.5).unwrap().abs() < 1e-10);
        let t = d.inverse_ccdf(0.25).unwrap();
        assert!(t.is_finite());
        assert!((d.ccdf(t) - 0.25).abs() < 1e-8);
    }

    #[test]
    fn extreme_left_tail_matches_high_precision_reference() {
        let d = StudentsT::new(100.0);
        let t = -6.5;
        let expected_cdf = 1.589_507_013_117_725_5e-9;
        let expected_sf = 0.999_999_998_410_493;
        assert!((d.cdf(t) - expected_cdf).abs() < 1e-23);
        assert!((d.ccdf(t) - expected_sf).abs() < 1e-15);
    }

    #[test]
    fn pdf_ln_pdf_and_entropy_are_finite() {
        let d = StudentsT::new(5.0);
        let x = 1.25;
        let ln_pdf = d.ln_pdf(x);
        assert!(ln_pdf.is_finite());
        assert!((d.pdf(x) - ln_pdf.exp()).abs() < 1e-15);
        assert!(d.entropy().is_finite());
    }

    #[test]
    fn nan_argument_gives_nan() {
        let d = StudentsT::new(5e-324);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
        assert!(StudentsT::new(10.0).pdf(f64::NAN).is_nan());
    }

    #[test]
    fn log_density_where_t_squared_overflows() {
        let ln_pdf = StudentsT::new(0.1).ln_pdf(-1e160);
        assert!((ln_pdf + 408.43131890946705).abs() < 1e-10);
        assert_eq!(
            StudentsT::new(10.0).ln_pdf(f64::INFINITY),
            f64::NEG_INFINITY
        );
    }
}
