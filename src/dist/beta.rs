use crate::error::SearchError;
use crate::search::{dstinv, dstzr};
use crate::special::beta_inc;
use crate::special::{beta_log, psi};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use thiserror::Error;

// Parameters of cdfbet (cdflib.f90:2559-2571).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Β distribution with shape parameters *a* > 0 and *b* > 0.
///
/// Defined over the interval [0 . . 1], with density
/// *f*(*x*; *a*, *b*) = *xᵃ* ⁻ ¹ (1 − *x*)*ᵇ* ⁻ ¹ / Β(*a*, *b*).
///
/// The methods correspond to CDFLIB's `cdfbet` (cdflib.f90:2451): `which = 1`
/// is [`cdf`] / [`ccdf`], `which = 2` is [`inverse_cdf`] /
/// [`inverse_ccdf`], `which = 3` is [`search_a`], `which = 4` is
/// [`search_b`].
///
/// # Example
///
/// ```
/// use cdflib::Beta;
/// use cdflib::traits::ContinuousCdf;
///
/// let b = Beta::new(2.0, 5.0);
///
/// // Pr[X ≤ 0.3]
/// let p = b.cdf(0.3);
///
/// // Compute parameter a given Pr[X ≤ 0.5] = 0.9 and b = 2.0
/// let a = Beta::search_a(0.9, 0.1, 0.5, 2.0).unwrap();
/// ```
///
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`inverse_ccdf`]: Beta::inverse_ccdf
/// [`search_a`]: Beta::search_a
/// [`search_b`]: Beta::search_b
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beta {
    a: f64,
    b: f64,
}

/// Errors arising from constructing a [`Beta`] or from its parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfbet`.
///
/// [`Beta`]: crate::Beta
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum BetaError {
    /// The shape parameter *a* was not strictly positive (`cdfbet` status −6).
    /// [`search_a`] also returns it, checked only in Rust, when the *a* it
    /// computes is 0.
    ///
    /// [`search_a`]: crate::Beta::search_a
    #[error("shape parameter `a` must be positive, got {0}")]
    ANotPositive(f64),
    /// The shape parameter *a* was not finite (checked only in Rust).
    #[error("shape parameter `a` must be finite, got {0}")]
    ANotFinite(f64),
    /// The shape parameter *b* was not strictly positive (`cdfbet` status −7).
    /// [`search_b`] also returns it, checked only in Rust, when the *b* it
    /// computes is 0.
    ///
    /// [`search_b`]: crate::Beta::search_b
    #[error("shape parameter `b` must be positive, got {0}")]
    BNotPositive(f64),
    /// The shape parameter *b* was not finite (checked only in Rust).
    #[error("shape parameter `b` must be finite, got {0}")]
    BNotFinite(f64),
    /// The argument *x* fell outside [0 . . 1] (`cdfbet` status −4); NaN is also
    /// rejected.
    #[error("argument x must be in [0..1], got {0}")]
    XOutOfRange(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdfbet` status −2); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfbet` status −3); NaN is
    /// also rejected.
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfbet` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfbet` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative Β distribution (*cum*, *ccum*) at (*x*, *y*)
/// with parameters *a* and *b* (`cumbet`, cdflib.f90:6723).
///
/// *y* = 1 − *x* is passed separately to preserve its precision.
#[inline]
pub(crate) fn cumbet(x: f64, y: f64, a: f64, b: f64) -> (f64, f64) {
    if x <= 0.0 {
        (0.0, 1.0)
    } else if y <= 0.0 {
        (1.0, 0.0)
    } else {
        // F90 ignores the ierr of beta_inc. Every caller passes x + y = 1
        // and a, b >= 0 not both 0, so with 0 < x and 0 < y beta_inc has no
        // error exit.
        beta_inc(a, b, x, y)
    }
}

// cdflib.f90:2601-2619 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), BetaError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(BetaError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:2621-2639 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), BetaError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(BetaError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:2641-2659 (status -4). Rust also rejects NaN.
#[inline]
fn check_x(x: f64) -> Result<(), BetaError> {
    if x < 0.0 || 1.0 < x || x.is_nan() {
        return Err(BetaError::XOutOfRange(x));
    }
    Ok(())
}

// cdflib.f90:2681-2692 (status -6). Rust also rejects a non-finite a.
#[inline]
fn check_a(a: f64) -> Result<(), BetaError> {
    if a <= 0.0 {
        return Err(BetaError::ANotPositive(a));
    }
    if !a.is_finite() {
        return Err(BetaError::ANotFinite(a));
    }
    Ok(())
}

// cdflib.f90:2694-2705 (status -7). Rust also rejects a non-finite b.
#[inline]
fn check_b(b: f64) -> Result<(), BetaError> {
    if b <= 0.0 {
        return Err(BetaError::BNotPositive(b));
    }
    if !b.is_finite() {
        return Err(BetaError::BNotFinite(b));
    }
    Ok(())
}

// cdflib.f90:2707-2717 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), BetaError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(BetaError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl Beta {
    /// Construct a Β(*a*, *b*) distribution with the given shape parameters.
    ///
    /// # Panics
    ///
    /// Panics if either parameter is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(a: f64, b: f64) -> Self {
        Self::try_new(a, b).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a [`BetaError`]
    /// instead of panicking.
    ///
    /// Returns [`ANotPositive`], [`ANotFinite`], [`BNotPositive`], or
    /// [`BNotFinite`] if either parameter fails its validity check.
    ///
    /// [`ANotPositive`]: BetaError::ANotPositive
    /// [`ANotFinite`]: BetaError::ANotFinite
    /// [`BNotPositive`]: BetaError::BNotPositive
    /// [`BNotFinite`]: BetaError::BNotFinite
    #[inline]
    pub fn try_new(a: f64, b: f64) -> Result<Self, BetaError> {
        check_a(a)?;
        check_b(b)?;
        Ok(Self { a, b })
    }

    /// Returns the shape parameter *a*.
    #[inline]
    pub const fn a(&self) -> f64 {
        self.a
    }

    /// Returns the shape parameter *b*.
    #[inline]
    pub const fn b(&self) -> f64 {
        self.b
    }

    /// Returns the shape parameter *a* satisfying Pr[*X* ≤ *x*] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfbet` with `which = 3`, with *y* = 1 − *x*. The caller
    /// passes both *p* and *q* = 1 − *p*, so that a small value of either
    /// keeps its precision; they must sum to 1 within 3ε.
    ///
    /// A computed *a* of 0, at the lower end of the search interval, is
    /// reported as [`ANotPositive`].
    ///
    /// [`ANotPositive`]: BetaError::ANotPositive
    #[inline]
    pub fn search_a(p: f64, q: f64, x: f64, b: f64) -> Result<f64, BetaError> {
        check_p(p)?;
        check_q(q)?;
        check_x(x)?;
        // Rust only: the API takes x alone, so y = 1 - x and the F90 checks
        // on y (status -5, cdflib.f90:2661-2679) and on x + y (status 4,
        // cdflib.f90:2719-2729) cannot fail.
        let y = 1.0 - x;
        check_b(b)?;
        check_pq(p, q)?;

        // cdflib.f90:2800-2840
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let a = d.x();
            let (cum, ccum) = cumbet(x, y, a, b);
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
        // lower end of the search interval it is not a valid a.
        let a = d.x();
        check_a(a)?;
        Ok(a)
    }

    /// Returns the shape parameter *b* satisfying Pr[*X* ≤ *x*] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfbet` with `which = 4`, with *y* = 1 − *x*. The caller
    /// passes both *p* and *q* = 1 − *p*, so that a small value of either
    /// keeps its precision; they must sum to 1 within 3ε.
    ///
    /// A computed *b* of 0, at the lower end of the search interval, is
    /// reported as [`BNotPositive`].
    ///
    /// [`BNotPositive`]: BetaError::BNotPositive
    #[inline]
    pub fn search_b(p: f64, q: f64, x: f64, a: f64) -> Result<f64, BetaError> {
        check_p(p)?;
        check_q(q)?;
        check_x(x)?;
        // Rust only: the API takes x alone, so y = 1 - x and the F90 checks
        // on y (status -5, cdflib.f90:2661-2679) and on x + y (status 4,
        // cdflib.f90:2719-2729) cannot fail.
        let y = 1.0 - x;
        check_a(a)?;
        check_pq(p, q)?;

        // cdflib.f90:2846-2884
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let b = d.x();
            let (cum, ccum) = cumbet(x, y, a, b);
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
        // lower end of the search interval it is not a valid b.
        let b = d.x();
        check_b(b)?;
        Ok(b)
    }

    /// CDFLIB's `cdfbet` with `which = 2`: returns (*x*, *y*) given (*p*,
    /// *q*), already checked.
    fn search_x_y(&self, p: f64, q: f64) -> Result<(f64, f64), BetaError> {
        let a = self.a;
        let b = self.b;
        // cdflib.f90:2742-2794
        let search = dstzr(0.0, 1.0, ATOL, TOL);
        let (x, y, z) = if p <= q {
            let mut z = search.dzror();
            let mut y = 1.0 - z.x();
            while z.status() == 1 {
                let (cum, _ccum) = cumbet(z.x(), y, a, b);
                let fx = cum - p;
                z.dzror(fx);
                y = 1.0 - z.x();
            }
            (z.x(), y, z)
        } else {
            let mut z = search.dzror();
            let mut x = 1.0 - z.x();
            while z.status() == 1 {
                let (_cum, ccum) = cumbet(x, z.x(), a, b);
                let fx = ccum - q;
                z.dzror(fx);
                x = 1.0 - z.x();
            }
            (x, z.x(), z)
        };
        if z.status() == -1 {
            return Err(if z.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: 1.0 }
            }
            .into());
        }
        Ok((x, y))
    }

    /// Returns the quantile *x* such that [ccdf]\(*x*\) = *q*.
    ///
    /// CDFLIB's `cdfbet` with `which = 2`, with *p* = 1 − *q*.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, BetaError> {
        check_q(q)?;
        // Rust only: exact endpoints.
        if q == 1.0 {
            return Ok(0.0);
        }
        if q == 0.0 {
            return Ok(1.0);
        }
        let p = 1.0 - q;
        Ok(self.search_x_y(p, q)?.0)
    }
}

impl ContinuousCdf for Beta {
    type Error = BetaError;

    /// CDFLIB's `cdfbet` with `which = 1`, with *y* = 1 − *x*.
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: no status -4 for x outside [0..1] (cdflib.f90:2641-2659);
        // cumbet returns (0, 1) for x <= 0 and (1, 0) for 1 < x.
        // cdflib.f90:2735
        cumbet(x, 1.0 - x, self.a, self.b).0
    }

    /// CDFLIB's `cdfbet` with `which = 1`, with *y* = 1 − *x*.
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        // Rust only: no status -4 for x outside [0..1] (cdflib.f90:2641-2659);
        // cumbet returns (0, 1) for x <= 0 and (1, 0) for 1 < x.
        // cdflib.f90:2735
        cumbet(x, 1.0 - x, self.a, self.b).1
    }

    /// CDFLIB's `cdfbet` with `which = 2`, with *q* = 1 − *p*.
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, BetaError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(1.0);
        }
        let q = 1.0 - p;
        Ok(self.search_x_y(p, q)?.0)
    }
}

impl Continuous for Beta {
    #[inline]
    fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 || 1.0 < x {
            return 0.0;
        }
        self.ln_pdf(x).exp()
    }
    #[inline]
    fn ln_pdf(&self, x: f64) -> f64 {
        if x < 0.0 || 1.0 < x {
            return f64::NEG_INFINITY;
        }
        // For a = 1 the term (a - 1) ln x is 0 at x = 0, and for b = 1 the
        // term (b - 1) ln(1 - x) is 0 at x = 1, where they would be
        // 0 · (-inf).
        let ln_x_term = if self.a == 1.0 && x == 0.0 {
            0.0
        } else {
            (self.a - 1.0) * x.ln()
        };
        let ln_y_term = if self.b == 1.0 && x == 1.0 {
            0.0
        } else {
            (self.b - 1.0) * (1.0 - x).ln()
        };
        ln_x_term + ln_y_term - beta_log(self.a, self.b)
    }
}

impl Mean for Beta {
    #[inline]
    fn mean(&self) -> f64 {
        self.a / (self.a + self.b)
    }
}

impl Variance for Beta {
    #[inline]
    fn variance(&self) -> f64 {
        let s = self.a + self.b;
        self.a * self.b / (s * s * (s + 1.0))
    }
}

impl Entropy for Beta {
    #[inline]
    fn entropy(&self) -> f64 {
        // H = ln Β(a,b) - (a-1)ψ(a) - (b-1)ψ(b) + (a+b-2)ψ(a+b)
        beta_log(self.a, self.b) - (self.a - 1.0) * psi(self.a) - (self.b - 1.0) * psi(self.b)
            + (self.a + self.b - 2.0) * psi(self.a + self.b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_parameters() {
        assert!(matches!(
            Beta::try_new(0.0, 1.0),
            Err(BetaError::ANotPositive(0.0))
        ));
        assert!(matches!(
            Beta::try_new(1.0, 0.0),
            Err(BetaError::BNotPositive(0.0))
        ));
        assert!(matches!(
            Beta::try_new(f64::NAN, 1.0),
            Err(BetaError::ANotFinite(_))
        ));
        assert!(matches!(
            Beta::try_new(1.0, f64::INFINITY),
            Err(BetaError::BNotFinite(_))
        ));
    }

    #[test]
    fn inverse_boundaries_and_density_edges() {
        let d = Beta::new(2.0, 3.0);
        assert_eq!(d.cdf(0.0), 0.0);
        assert_eq!(d.cdf(1.0), 1.0);
        assert_eq!(d.ccdf(0.0), 1.0);
        assert_eq!(d.ccdf(1.0), 0.0);
        assert_eq!(d.inverse_cdf(0.0).unwrap(), 0.0);
        assert_eq!(d.inverse_cdf(1.0).unwrap(), 1.0);
        assert_eq!(d.inverse_ccdf(1.0).unwrap(), 0.0);
        assert_eq!(d.inverse_ccdf(0.0).unwrap(), 1.0);
        assert_eq!(d.pdf(0.0), 0.0);
        assert_eq!(d.pdf(1.0), 0.0);
        assert_eq!(d.ln_pdf(0.0), f64::NEG_INFINITY);
        assert_eq!(d.ln_pdf(1.0), f64::NEG_INFINITY);
        assert!(d.pdf(0.4).is_finite());
        assert!(d.ln_pdf(0.4).is_finite());
        assert!(d.inverse_ccdf(0.4).unwrap().is_finite());
        assert!(d.mean().is_finite());
        assert!(d.variance().is_finite());
        assert!(d.entropy().is_finite());
    }

    #[test]
    fn search_parameter_rejects_invalid_inputs() {
        assert!(matches!(
            Beta::search_a(-0.1, 1.1, 0.5, 2.0),
            Err(BetaError::PNotInRange(-0.1))
        ));
        assert!(matches!(
            Beta::search_a(0.5, 0.5, 0.5, 0.0),
            Err(BetaError::BNotPositive(0.0))
        ));
        assert!(matches!(
            Beta::search_b(0.5, 0.5, 0.5, 0.0),
            Err(BetaError::ANotPositive(0.0))
        ));
        assert!(matches!(
            Beta::search_a(0.5, 0.5, 1.5, 2.0),
            Err(BetaError::XOutOfRange(1.5))
        ));
        assert!(matches!(
            Beta::search_b(0.5, 0.5, -0.1, 2.0),
            Err(BetaError::XOutOfRange(x)) if x == -0.1
        ));
    }

    #[test]
    fn density_at_the_ends_is_the_limit() {
        assert!((Beta::new(1.0, 3.0).pdf(0.0) - 3.0).abs() < 1e-12);
        assert!((Beta::new(2.0, 1.0).pdf(1.0) - 2.0).abs() < 1e-12);
        assert_eq!(Beta::new(0.5, 2.0).pdf(0.0), f64::INFINITY);
        assert_eq!(Beta::new(2.0, 0.5).pdf(1.0), f64::INFINITY);
        assert_eq!(Beta::new(2.0, 1.0).pdf(0.0), 0.0);
        assert_eq!(Beta::new(2.0, 1.0).pdf(1.5), 0.0);
        assert!(Beta::new(1.0, 1.0).pdf(f64::NAN).is_nan());
        assert!(Beta::new(1.0, 1.0).ln_pdf(f64::NAN).is_nan());
    }
}
