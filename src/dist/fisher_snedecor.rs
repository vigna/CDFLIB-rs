use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::beta_inc;
use crate::special::{beta_log, psi};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use thiserror::Error;

// Parameters of cdff (cdflib.f90:4132-4148).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Fisher–Snedecor (*F*) distribution with *dfn* numerator and *dfd*
/// denominator degrees of freedom.
///
/// The CDF reduces to the incomplete Β (Abramowitz–Stegun 26.6.2).
///
/// The methods correspond to CDFLIB's `cdff` (cdflib.f90:4028):
/// `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is [`inverse_cdf`] /
/// [`inverse_ccdf`], `which = 3` is [`search_dfn`], `which = 4` is
/// [`search_dfd`]. As CDFLIB warns, the CDF is not necessarily monotone in
/// either degrees of freedom, so the searches assume monotonicity and find
/// one of possibly two solutions.
///
/// # Example
///
/// ```
/// use cdflib::FisherSnedecor;
/// use cdflib::traits::ContinuousCdf;
///
/// let f = FisherSnedecor::new(5.0, 10.0);
///
/// // Pr[X ≤ 3.33]
/// let p = f.cdf(3.33);
/// assert!((p - 0.9501687242027787).abs() < 1e-12);
///
/// // Compute numerator df given Pr[X ≤ 3.33] = 0.95 and dfd = 10
/// let dfn = FisherSnedecor::search_dfn(0.95, 0.05, 3.33, 10.0).unwrap();
/// assert!((dfn - 4.967304933).abs() < 1e-6);
/// ```
///
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`inverse_ccdf`]: FisherSnedecor::inverse_ccdf
/// [`search_dfn`]: FisherSnedecor::search_dfn
/// [`search_dfd`]: FisherSnedecor::search_dfd
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FisherSnedecor {
    dfn: f64,
    dfd: f64,
}

/// Errors arising from constructing a [`FisherSnedecor`] or from its
/// parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdff`.
///
/// [`FisherSnedecor`]: crate::FisherSnedecor
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum FisherSnedecorError {
    /// The numerator degrees of freedom *dfn* was not strictly positive
    /// (`cdff` status −5). Rust also rejects the smallest subnormal number,
    /// whose half, which `cumf` passes to `beta_inc`, is 0.
    #[error("numerator df must be positive and not the smallest subnormal, got {0}")]
    DfnNotPositive(f64),
    /// The numerator degrees of freedom *dfn* was not finite (checked only
    /// in Rust).
    #[error("numerator df must be finite, got {0}")]
    DfnNotFinite(f64),
    /// The denominator degrees of freedom *dfd* was not strictly positive
    /// (`cdff` status −6). Rust also rejects the smallest subnormal number,
    /// whose half, which `cumf` passes to `beta_inc`, is 0.
    #[error("denominator df must be positive and not the smallest subnormal, got {0}")]
    DfdNotPositive(f64),
    /// The denominator degrees of freedom *dfd* was not finite (checked
    /// only in Rust).
    #[error("denominator df must be finite, got {0}")]
    DfdNotFinite(f64),
    /// The value *f* (the point at which the CDF is evaluated) was not
    /// strictly positive (`cdff` status −4 for *f* < 0; Rust also rejects
    /// *f* = 0, where Pr[*X* ≤ *f*] = 0 whatever the parameters).
    #[error("f must be positive, got {0}")]
    FNotPositive(f64),
    /// The value *f* was not finite (checked only in Rust).
    #[error("f must be finite, got {0}")]
    FNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdff` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdff` status −3); NaN is
    /// also rejected.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdff` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdff` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative *F* distribution (*cum*, *ccum*) at *f* with
/// *dfn* and *dfd* degrees of freedom (`cumf`, cdflib.f90:7134).
///
/// The complement of *xx* = *dfd* / (*dfd* + *dfn*·*f*) is computed
/// directly when *xx* > 0.5 to avoid cancellation.
#[inline]
pub(crate) fn cumf(f: f64, dfn: f64, dfd: f64) -> (f64, f64) {
    if f <= 0.0 {
        return (0.0, 1.0);
    }
    let prod = dfn * f;
    let dsum = dfd + prod;
    let mut xx = dfd / dsum;
    let yy;
    if 0.5 < xx {
        yy = prod / dsum;
        xx = 1.0 - yy;
    } else {
        yy = 1.0 - xx;
    }
    // F90 ignores the ierr of beta_inc. With dfn and dfd positive and
    // halving to nonzero values, which the Rust checks ensure, beta_inc
    // has no error exit.
    let (ccum, cum) = beta_inc(0.5 * dfd, 0.5 * dfn, xx, yy);
    (cum, ccum)
}

// cdflib.f90:4175-4194 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), FisherSnedecorError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(FisherSnedecorError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:4195-4214 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), FisherSnedecorError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(FisherSnedecorError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:4215-4227 (status -4). Rust also rejects f = 0 and a
// non-finite f.
#[inline]
fn check_f(f: f64) -> Result<(), FisherSnedecorError> {
    if f <= 0.0 {
        return Err(FisherSnedecorError::FNotPositive(f));
    }
    if !f.is_finite() {
        return Err(FisherSnedecorError::FNotFinite(f));
    }
    Ok(())
}

// cdflib.f90:4228-4240 (status -5). Rust also rejects a non-finite dfn.
#[inline]
fn check_dfn(dfn: f64) -> Result<(), FisherSnedecorError> {
    if dfn <= 0.0 {
        return Err(FisherSnedecorError::DfnNotPositive(dfn));
    }
    // Rust only: cumf then passes b = 0.5 * dfn = 0 to beta_inc, which
    // fails when yy = dfn * f / (dfd + dfn * f) is 0 (ierr 7) or when dfd
    // halves to 0 too (ierr 2); the F90 ignores the error.
    if 0.5 * dfn == 0.0 {
        return Err(FisherSnedecorError::DfnNotPositive(dfn));
    }
    if !dfn.is_finite() {
        return Err(FisherSnedecorError::DfnNotFinite(dfn));
    }
    Ok(())
}

// cdflib.f90:4241-4253 (status -6). Rust also rejects a non-finite dfd.
#[inline]
fn check_dfd(dfd: f64) -> Result<(), FisherSnedecorError> {
    if dfd <= 0.0 {
        return Err(FisherSnedecorError::DfdNotPositive(dfd));
    }
    // Rust only: cumf then passes a = 0.5 * dfd = 0 to beta_inc, which
    // fails when xx = dfd / (dfd + dfn * f) is 0, as soon as dfn * f >= 2
    // (ierr 6), or when dfn halves to 0 too (ierr 2); the F90 ignores the
    // error.
    if 0.5 * dfd == 0.0 {
        return Err(FisherSnedecorError::DfdNotPositive(dfd));
    }
    if !dfd.is_finite() {
        return Err(FisherSnedecorError::DfdNotFinite(dfd));
    }
    Ok(())
}

// cdflib.f90:4254-4265 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), FisherSnedecorError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(FisherSnedecorError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl FisherSnedecor {
    /// Construct an *F*(*dfn*, *dfd*) distribution with strictly positive
    /// numerator and denominator degrees of freedom.
    ///
    /// # Panics
    ///
    /// Panics if either argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(dfn: f64, dfd: f64) -> Self {
        Self::try_new(dfn, dfd).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`FisherSnedecorError`] instead of panicking.
    ///
    /// Returns [`DfnNotPositive`], [`DfnNotFinite`], [`DfdNotPositive`], or
    /// [`DfdNotFinite`] if either argument fails its validity check.
    ///
    /// [`DfnNotPositive`]: FisherSnedecorError::DfnNotPositive
    /// [`DfnNotFinite`]: FisherSnedecorError::DfnNotFinite
    /// [`DfdNotPositive`]: FisherSnedecorError::DfdNotPositive
    /// [`DfdNotFinite`]: FisherSnedecorError::DfdNotFinite
    #[inline]
    pub fn try_new(dfn: f64, dfd: f64) -> Result<Self, FisherSnedecorError> {
        check_dfn(dfn)?;
        check_dfd(dfd)?;
        Ok(Self { dfn, dfd })
    }

    /// Returns the numerator degrees of freedom *dfn*.
    #[inline]
    pub const fn dfn(&self) -> f64 {
        self.dfn
    }

    /// Returns the denominator degrees of freedom *dfd*.
    #[inline]
    pub const fn dfd(&self) -> f64 {
        self.dfd
    }

    /// Returns the numerator degrees of freedom *dfn* satisfying
    /// Pr[*X* ≤ *f*] = *p*, searched for in [1 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdff` with `which = 3`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε.
    #[inline]
    pub fn search_dfn(p: f64, q: f64, f: f64, dfd: f64) -> Result<f64, FisherSnedecorError> {
        check_p(p)?;
        check_q(q)?;
        check_f(f)?;
        check_dfd(dfd)?;
        check_pq(p, q)?;

        // cdflib.f90:4329-4374. The lower bound is 1, since dfn = 0 makes
        // beta_inc fail inside cumf (cdflib.f90:4322-4325).
        let bound_lo = 1.0;
        let bound_hi = INF;
        let mut d = dstinv(bound_lo, bound_hi, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let dfn = d.x();
            let (cum, ccum) = cumf(f, dfn, dfd);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: bound_lo }
            } else {
                SearchError::AnswerAboveUpperBound { bound: bound_hi }
            }
            .into());
        }
        Ok(d.x())
    }

    /// Returns the denominator degrees of freedom *dfd* satisfying
    /// Pr[*X* ≤ *f*] = *p*, searched for in [1 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdff` with `which = 4`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε.
    #[inline]
    pub fn search_dfd(p: f64, q: f64, f: f64, dfn: f64) -> Result<f64, FisherSnedecorError> {
        check_p(p)?;
        check_q(q)?;
        check_f(f)?;
        check_dfn(dfn)?;
        check_pq(p, q)?;

        // cdflib.f90:4385-4426. The lower bound is 1, since dfd = 0 makes
        // beta_inc fail inside cumf (cdflib.f90:4378-4381).
        let bound_lo = 1.0;
        let bound_hi = INF;
        let mut d = dstinv(bound_lo, bound_hi, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let dfd = d.x();
            let (cum, ccum) = cumf(f, dfn, dfd);
            let fx = if p <= q { cum - p } else { ccum - q };
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: bound_lo }
            } else {
                SearchError::AnswerAboveUpperBound { bound: bound_hi }
            }
            .into());
        }
        Ok(d.x())
    }

    /// CDFLIB's `cdff` with `which = 2`: returns *f* given (*p*, *q*),
    /// already checked.
    fn search_f(&self, p: f64, q: f64) -> Result<f64, FisherSnedecorError> {
        let dfn = self.dfn;
        let dfd = self.dfd;
        // cdflib.f90:4278-4318
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let f = d.x();
            let (cum, ccum) = cumf(f, dfn, dfd);
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

    /// Returns the quantile *f* such that [ccdf]\(*f*\) = *q*, searched for
    /// in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdff` with `which = 2`, with *p* = 1 − *q*.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, FisherSnedecorError> {
        check_q(q)?;
        // Rust only: exact endpoints.
        if q == 1.0 {
            return Ok(0.0);
        }
        if q == 0.0 {
            return Ok(f64::INFINITY);
        }
        let p = 1.0 - q;
        self.search_f(p, q)
    }
}

impl ContinuousCdf for FisherSnedecor {
    type Error = FisherSnedecorError;

    /// CDFLIB's `cdff` with `which = 1`.
    ///
    /// The result can be NaN when *dfn* or *dfd* is above about
    /// 2 · 10³⁰⁷, where the F90 `beta_inc` never returns; the inverses and
    /// the parameter searches then return meaningless values (see
    /// [`SearchError`]).
    ///
    /// [`SearchError`]: crate::SearchError
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x, which cumf passes to beta_inc; as in
        // the F90, beta_inc can then return a value computed from the other
        // arguments.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: no status -4 for f < 0 (cdflib.f90:4215-4227); cumf
        // returns (0, 1) there.
        // cdflib.f90:4271
        cumf(x, self.dfn, self.dfd).0
    }

    /// CDFLIB's `cdff` with `which = 1`.
    ///
    /// The result can be NaN for huge parameters, as for [`cdf`].
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x, as in cdf.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: no status -4 for f < 0 (cdflib.f90:4215-4227); cumf
        // returns (0, 1) there.
        // cdflib.f90:4271
        cumf(x, self.dfn, self.dfd).1
    }

    /// CDFLIB's `cdff` with `which = 2`, with *q* = 1 − *p*.
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, FisherSnedecorError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let q = 1.0 - p;
        self.search_f(p, q)
    }
}

impl Continuous for FisherSnedecor {
    #[inline]
    fn pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            return 0.0;
        }
        self.ln_pdf(x).exp()
    }
    #[inline]
    fn ln_pdf(&self, x: f64) -> f64 {
        if x < 0.0 {
            return f64::NEG_INFINITY;
        }
        // The density tends to 0 as x tends to +inf, where the expression
        // below would be inf - inf, or 0 * inf for dfn = 2.
        if x == f64::INFINITY {
            return f64::NEG_INFINITY;
        }
        let dfn = self.dfn;
        let dfd = self.dfd;
        // f(x) = (dfn/dfd)^(dfn/2) · x^(dfn/2-1) · (1 + dfn·x/dfd)^(-(dfn+dfd)/2) / Β(dfn/2, dfd/2)
        let half_dfn = dfn / 2.0;
        let half_dfd = dfd / 2.0;
        // For dfn = 2 the term (dfn/2 - 1) ln x is 0 at x = 0, where it
        // would be 0 · (-inf).
        let ln_x_term = if half_dfn == 1.0 && x == 0.0 {
            0.0
        } else {
            (half_dfn - 1.0) * x.ln()
        };
        // ln(dfn/dfd) is a difference of logarithms where the ratio
        // overflows, underflows or is subnormal.
        let ratio = dfn / dfd;
        let ln_ratio = if ratio != 0.0 && ratio.is_finite() {
            ratio.ln()
        } else {
            dfn.ln() - dfd.ln()
        };
        // y = dfn·x/dfd. For y ≥ 2^53, ln(1 + y) = ln y in double precision;
        // substituting it, the density is computed in the form
        // ln f(x) = -(dfd/2) ln(dfn/dfd) - (1 + dfd/2) ln x - ln Β(dfn/2, dfd/2),
        // which cannot overflow and has no cancellation between large terms.
        let ln_y = ln_ratio + x.ln();
        if ln_y >= 53.0 * std::f64::consts::LN_2 {
            return -half_dfd * ln_ratio - (1.0 + half_dfd) * x.ln() - beta_log(half_dfn, half_dfd);
        }
        let y = if ratio != 0.0 && ratio.is_finite() {
            ratio * x
        } else {
            ln_y.exp()
        };
        // ln_1p keeps the precision of ln(1 + dfn·x/dfd) when dfd is large.
        half_dfn * ln_ratio + ln_x_term
            - (half_dfn + half_dfd) * y.ln_1p()
            - beta_log(half_dfn, half_dfd)
    }
}

impl Mean for FisherSnedecor {
    /// Defined for *dfd* > 2.
    #[inline]
    fn mean(&self) -> f64 {
        if self.dfd > 2.0 {
            self.dfd / (self.dfd - 2.0)
        } else {
            f64::NAN
        }
    }
}

impl Variance for FisherSnedecor {
    /// Defined for *dfd* > 4.
    #[inline]
    fn variance(&self) -> f64 {
        let dfn = self.dfn;
        let dfd = self.dfd;
        if dfd > 4.0 {
            // 2 dfd² (dfn + dfd - 2) / (dfn (dfd - 2)² (dfd - 4)), written
            // with m = dfd / (dfd - 2), the mean, and (dfn + dfd - 2) /
            // (dfn (dfd - 4)) = 1/dfn + (1 + 2/dfn) / (dfd - 4), so that an
            // intermediate overflows only where the variance does.
            let m = dfd / (dfd - 2.0);
            2.0 * m * m * (1.0 / dfn + (1.0 + 2.0 / dfn) / (dfd - 4.0))
        } else {
            f64::NAN
        }
    }
}

impl Entropy for FisherSnedecor {
    /// The result is NaN when both *dfn* and *dfd* are below about
    /// 1.1 · 10⁻³⁰⁸, where *ψ*(*dfn*/2) and *ψ*(*dfd*/2) overflow to −∞
    /// and their terms cancel.
    #[inline]
    fn entropy(&self) -> f64 {
        // Closed-form: H = ln(dfd/dfn · Β(dfn/2, dfd/2))
        //                + (1 - dfn/2) ψ(dfn/2) - (1 + dfd/2) ψ(dfd/2)
        //                + (dfn+dfd)/2 · ψ((dfn+dfd)/2)
        // The logarithm of dfd/dfn is a difference so that the ratio cannot
        // overflow or underflow.
        let dfn = self.dfn;
        let dfd = self.dfd;
        dfd.ln() - dfn.ln() + beta_log(dfn / 2.0, dfd / 2.0) + (1.0 - dfn / 2.0) * psi(dfn / 2.0)
            - (1.0 + dfd / 2.0) * psi(dfd / 2.0)
            + 0.5 * (dfn + dfd) * psi((dfn + dfd) / 2.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_parameters() {
        assert!(matches!(
            FisherSnedecor::try_new(0.0, 1.0),
            Err(FisherSnedecorError::DfnNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecor::try_new(1.0, 0.0),
            Err(FisherSnedecorError::DfdNotPositive(0.0))
        ));
    }

    #[test]
    fn inverse_and_density_edges() {
        let d = FisherSnedecor::new(5.0, 10.0);
        assert_eq!(d.inverse_cdf(0.0).unwrap(), 0.0);
        assert_eq!(d.inverse_ccdf(1.0).unwrap(), 0.0);
        assert_eq!(d.pdf(0.0), 0.0);
        assert_eq!(d.ln_pdf(0.0), f64::NEG_INFINITY);
        assert!(d.inverse_ccdf(0.25).unwrap().is_finite());
        assert!(d.pdf(1.5).is_finite());
        assert!(d.ln_pdf(1.5).is_finite());
        assert!(d.entropy().is_finite());
    }

    #[test]
    fn moment_thresholds_and_invalid_solves() {
        assert!(FisherSnedecor::new(5.0, 2.0).mean().is_nan());
        assert!(FisherSnedecor::new(5.0, 4.0).variance().is_nan());
        assert!(FisherSnedecor::new(5.0, 10.0).mean().is_finite());
        assert!(FisherSnedecor::new(5.0, 10.0).variance().is_finite());
        assert!(matches!(
            FisherSnedecor::search_dfn(-0.1, 1.1, 1.0, 5.0),
            Err(FisherSnedecorError::PNotInRange(-0.1))
        ));
        assert!(matches!(
            FisherSnedecor::search_dfn(0.5, 0.5, 0.0, 5.0),
            Err(FisherSnedecorError::FNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecor::search_dfn(0.5, 0.5, 1.0, 0.0),
            Err(FisherSnedecorError::DfdNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecor::search_dfd(0.5, 0.5, 0.0, 5.0),
            Err(FisherSnedecorError::FNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecor::search_dfd(0.5, 0.5, 1.0, 0.0),
            Err(FisherSnedecorError::DfnNotPositive(0.0))
        ));
    }

    #[test]
    fn df_halving_to_zero_is_rejected() {
        // cumf passes dfn / 2 and dfd / 2 to beta_inc, which fails on 0.
        let tiny = f64::from_bits(1);
        assert_eq!(
            FisherSnedecor::try_new(tiny, 1.0),
            Err(FisherSnedecorError::DfnNotPositive(tiny))
        );
        assert_eq!(
            FisherSnedecor::try_new(1.0, tiny),
            Err(FisherSnedecorError::DfdNotPositive(tiny))
        );
        let d = FisherSnedecor::new(2.0 * tiny, 1.0);
        assert!(d.cdf(0.1).is_finite());
    }

    #[test]
    fn density_at_zero_is_the_limit() {
        assert!((FisherSnedecor::new(2.0, 5.0).pdf(0.0) - 1.0).abs() < 1e-12);
        assert_eq!(FisherSnedecor::new(1.0, 5.0).pdf(0.0), f64::INFINITY);
        assert_eq!(FisherSnedecor::new(3.0, 5.0).pdf(0.0), 0.0);
    }

    #[test]
    fn nan_argument_gives_nan() {
        let d = FisherSnedecor::new(1e-300, 1e-300);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
    }

    #[test]
    fn density_keeps_precision_for_large_dfd() {
        // The F(1, dfd) density at 1 tends to the χ²(1) density at 1.
        let pdf = FisherSnedecor::new(1.0, 1e17).pdf(1.0);
        assert!((pdf / 0.24197072451914335 - 1.0).abs() < 1e-12);
        let pdf = FisherSnedecor::new(2.0, 1e20).pdf(0.5);
        assert!((pdf / 0.6065306597126334 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn variance_does_not_overflow() {
        // For large dfd the variance tends to 2 / dfn.
        assert!((FisherSnedecor::new(5.0, 1e103).variance() - 0.4).abs() < 1e-12);
        assert!((FisherSnedecor::new(5.0, f64::MAX).variance() - 0.4).abs() < 1e-12);
    }

    #[test]
    fn entropy_with_a_subnormal_parameter_is_the_limit() {
        assert_eq!(
            FisherSnedecor::new(1e-323, 5.0).entropy(),
            f64::NEG_INFINITY
        );
        assert_eq!(FisherSnedecor::new(5.0, 1e-323).entropy(), f64::INFINITY);
    }
}
