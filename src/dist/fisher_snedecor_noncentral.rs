use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::beta_inc;
use crate::special::gamma_log;
use crate::traits::{ContinuousCdf, Mean, Variance};
use thiserror::Error;

use super::fisher_snedecor::cumf;

// Parameters of cdffnc (cdflib.f90:4571-4587).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e30;
const TENT4: f64 = 1.0e4;
const TOL: f64 = 1.0e-8;

/// Noncentral *F* distribution with numerator df *dfn*, denominator df *dfd*,
/// and noncentrality *λ* ≥ 0.
///
/// The methods correspond to CDFLIB's `cdffnc` (cdflib.f90:4432), whose
/// PNONC is *λ*: `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is
/// [`inverse_cdf`], `which = 3` is [`search_dfn`], `which = 4` is
/// [`search_dfd`], `which = 5` is [`search_ncp`]. CDFLIB's searches use
/// *p* only, and `cumfnc` requires *dfn* ≥ 1 and *dfd* ≥ 1.
///
/// # Notes
///
/// Neither [`Continuous`] nor [`Entropy`] is implemented.
///
/// # Example
///
/// ```
/// use cdflib::FisherSnedecorNoncentral;
/// use cdflib::traits::ContinuousCdf;
///
/// let d = FisherSnedecorNoncentral::new(5.0, 10.0, 2.0);
///
/// // Pr[X ≤ 4.0]
/// let p = d.cdf(4.0);
/// assert!((p - 0.9281722121).abs() < 1e-5);
///
/// // Compute noncentrality λ given Pr[X ≤ 4.0] = 0.5, dfn = 5, dfd = 10
/// let ncp = FisherSnedecorNoncentral::search_ncp(0.5, 4.0, 5.0, 10.0).unwrap();
/// assert!((ncp - 14.766).abs() < 1e-2);
/// ```
///
/// [`Continuous`]: crate::traits::Continuous
/// [`Entropy`]: crate::traits::Entropy
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`search_dfn`]: FisherSnedecorNoncentral::search_dfn
/// [`search_dfd`]: FisherSnedecorNoncentral::search_dfd
/// [`search_ncp`]: FisherSnedecorNoncentral::search_ncp
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FisherSnedecorNoncentral {
    dfn: f64,
    dfd: f64,
    ncp: f64,
}

/// Errors arising from constructing a [`FisherSnedecorNoncentral`] or from
/// its parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdffnc` and
/// to the fatal errors of `cumfnc`.
///
/// [`FisherSnedecorNoncentral`]: crate::FisherSnedecorNoncentral
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum FisherSnedecorNoncentralError {
    /// The numerator degrees of freedom *dfn* was not strictly positive
    /// (`cdffnc` status −5).
    #[error("numerator df must be positive, got {0}")]
    DfnNotPositive(f64),
    /// The numerator degrees of freedom *dfn* was less than 1, where
    /// `cumfnc` stops with a fatal error (cdflib.f90:7317-7322). It is
    /// checked after the `cdffnc` status checks.
    #[error("numerator df must be at least 1, got {0}")]
    DfnTooSmall(f64),
    /// The numerator degrees of freedom *dfn* was not finite (checked only
    /// in Rust).
    #[error("numerator df must be finite, got {0}")]
    DfnNotFinite(f64),
    /// The denominator degrees of freedom *dfd* was not strictly positive
    /// (`cdffnc` status −6).
    #[error("denominator df must be positive, got {0}")]
    DfdNotPositive(f64),
    /// The denominator degrees of freedom *dfd* was less than 1, where
    /// `cumfnc` stops with a fatal error (cdflib.f90:7324-7329). It is
    /// checked after the `cdffnc` status checks.
    #[error("denominator df must be at least 1, got {0}")]
    DfdTooSmall(f64),
    /// The denominator degrees of freedom *dfd* was not finite (checked
    /// only in Rust).
    #[error("denominator df must be finite, got {0}")]
    DfdNotFinite(f64),
    /// The noncentrality parameter *λ* was negative (`cdffnc` status −7).
    #[error("noncentrality parameter must be nonnegative, got {0}")]
    NcpNegative(f64),
    /// The noncentrality parameter *λ* was not finite (checked only in
    /// Rust).
    #[error("noncentrality parameter must be finite, got {0}")]
    NcpNotFinite(f64),
    /// The value *f* was not strictly positive (`cdffnc` status −4 for
    /// *f* < 0; Rust also rejects *f* = 0, where Pr[*X* ≤ *f*] = 0 whatever
    /// the parameters).
    #[error("f must be positive, got {0}")]
    FNotPositive(f64),
    /// The value *f* was not finite (checked only in Rust).
    #[error("f must be finite, got {0}")]
    FNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdffnc` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1]; NaN is also rejected.
    /// No method of [`FisherSnedecorNoncentral`] returns it, since CDFLIB's
    /// `cdffnc` does not use *q*.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The search for the answer failed (`cdffnc` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

/// Returns the cumulative noncentral *F* distribution (*cum*, *ccum*) at
/// *f* with *dfn* and *dfd* degrees of freedom and noncentrality parameter
/// *pnonc* (`cumfnc`, cdflib.f90:7219).
///
/// The series of incomplete Β functions weighted by Poisson terms is
/// summed backwards and forwards from the central term.
///
/// # Panics
///
/// Panics if *dfn* < 1 or *dfd* < 1, where CDFLIB stops with a fatal
/// error; the public API rejects such values first.
#[allow(clippy::assign_op_pattern)]
pub(crate) fn cumfnc(f: f64, dfn: f64, dfd: f64, pnonc: f64) -> (f64, f64) {
    const EPS: f64 = 0.0001;

    if f <= 0.0 {
        return (0.0, 1.0);
    }

    if dfn < 1.0 {
        panic!("cumfnc: dfn < 1");
    }

    if dfd < 1.0 {
        panic!("cumfnc: dfd < 1");
    }
    // Handle the case in which the noncentrality parameter is essentially
    // zero.
    if pnonc < 1.0e-10 {
        return cumf(f, dfn, dfd);
    }
    // Rust only: a NaN f, dfn or dfd passes the tests above, and the F90
    // then never returns, inside beta_inc or in the forward sum below.
    if f.is_nan() || dfn.is_nan() || dfd.is_nan() {
        return (f64::NAN, f64::NAN);
    }

    let xnonc = pnonc / 2.0;
    // Rust only: int(xnonc) and icent + 1 below overflow the default
    // integer, which is undefined in Fortran; a NaN pnonc is caught here
    // too.
    if xnonc.is_nan() || f64::from(i32::MAX) <= xnonc {
        panic!("cumfnc: integer overflow for pnonc = {pnonc}");
    }
    // Calculate the central term of the Poisson weighting factor.
    let mut icent = xnonc as i32;

    if icent == 0 {
        icent = 1;
    }
    // Compute the central weight term.
    let centwt = (-xnonc + icent as f64 * xnonc.ln() - gamma_log((icent + 1) as f64)).exp();
    // Compute the central incomplete Β term. Ensure that the minimum of
    // the argument to Β and 1 minus the argument is computed accurately.
    let prod = dfn * f;
    let dsum = dfd + prod;
    let mut yy = dfd / dsum;
    let xx;

    if 0.5 < yy {
        xx = prod / dsum;
        yy = 1.0 - xx;
    } else {
        xx = 1.0 - yy;
    }

    let arg1 = 0.5 * dfn + icent as f64;
    // F90 ignores the ierr of beta_inc; with dfn, dfd >= 1 and 0 < f,
    // beta_inc has no error exit.
    let (mut betdn, _dummy) = beta_inc(arg1, 0.5 * dfd, xx, yy);

    let mut adn = dfn / 2.0 + icent as f64;
    let mut aup = adn;
    let b = dfd / 2.0;
    let mut betup = betdn;
    let mut sum1 = centwt * betdn;
    // Now sum terms backward from icent until convergence or all done.
    let mut xmult = centwt;
    let mut i = icent;
    let mut dnterm =
        (gamma_log(adn + b) - gamma_log(adn + 1.0) - gamma_log(b) + adn * xx.ln() + b * yy.ln())
            .exp();

    loop {
        if i <= 0 {
            break;
        }

        if sum1 < f64::EPSILON || xmult * betdn < EPS * sum1 {
            break;
        }

        xmult = xmult * (i as f64 / xnonc);
        i = i - 1;
        adn = adn - 1.0;
        dnterm = (adn + 1.0) / ((adn + b) * xx) * dnterm;
        betdn = betdn + dnterm;
        sum1 = sum1 + xmult * betdn;
    }

    i = icent + 1;
    // Now sum forward until convergence.
    xmult = centwt;

    let expon = if (aup - 1.0 + b) == 0.0 {
        -gamma_log(aup) - gamma_log(b) + (aup - 1.0) * xx.ln() + b * yy.ln()
    } else {
        gamma_log(aup - 1.0 + b) - gamma_log(aup) - gamma_log(b)
            + (aup - 1.0) * xx.ln()
            + b * yy.ln()
    };
    // CDFLIB guards against exp producing values so small that combining
    // them with ordinary quantities raises floating-point errors.
    let mut upterm = if expon <= f64::EPSILON.ln() {
        0.0
    } else {
        expon.exp()
    };

    loop {
        // Rust only: i + 1 overflows the default integer, which is
        // undefined in Fortran.
        if i == i32::MAX {
            panic!("cumfnc: integer overflow for pnonc = {pnonc}");
        }
        // Rust only: once sum1 is NaN neither exit test below can succeed,
        // and the F90 loop never ends; this happens when dfn or dfd is so
        // large that the terms lose all their digits.
        if sum1.is_nan() {
            panic!("cumfnc: the sum is NaN for dfn = {dfn}, dfd = {dfd}, pnonc = {pnonc}");
        }
        xmult = xmult * (xnonc / i as f64);
        i = i + 1;
        aup = aup + 1.0;
        upterm = (aup + b - 2.0) * xx / (aup - 1.0) * upterm;
        betup = betup - upterm;
        sum1 = sum1 + xmult * betup;

        if sum1 < f64::EPSILON || xmult * betup < EPS * sum1 {
            break;
        }
    }

    let cum = sum1;
    let ccum = 0.5 + (0.5 - cum);
    (cum, ccum)
}

// cdflib.f90:4614-4633 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(FisherSnedecorNoncentralError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:4634-4646 (status -4). Rust also rejects f = 0 and a
// non-finite f.
#[inline]
fn check_f(f: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if f <= 0.0 {
        return Err(FisherSnedecorNoncentralError::FNotPositive(f));
    }
    if !f.is_finite() {
        return Err(FisherSnedecorNoncentralError::FNotFinite(f));
    }
    Ok(())
}

// cdflib.f90:4647-4659 (status -5). Rust also rejects a non-finite dfn.
#[inline]
fn check_dfn(dfn: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if dfn <= 0.0 {
        return Err(FisherSnedecorNoncentralError::DfnNotPositive(dfn));
    }
    if !dfn.is_finite() {
        return Err(FisherSnedecorNoncentralError::DfnNotFinite(dfn));
    }
    Ok(())
}

// cdflib.f90:4660-4672 (status -6). Rust also rejects a non-finite dfd.
#[inline]
fn check_dfd(dfd: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if dfd <= 0.0 {
        return Err(FisherSnedecorNoncentralError::DfdNotPositive(dfd));
    }
    if !dfd.is_finite() {
        return Err(FisherSnedecorNoncentralError::DfdNotFinite(dfd));
    }
    Ok(())
}

// cdflib.f90:4673-4685 (status -7). Rust also rejects a non-finite pnonc.
#[inline]
fn check_pnonc(pnonc: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if pnonc < 0.0 {
        return Err(FisherSnedecorNoncentralError::NcpNegative(pnonc));
    }
    if !pnonc.is_finite() {
        return Err(FisherSnedecorNoncentralError::NcpNotFinite(pnonc));
    }
    Ok(())
}

// cdflib.f90:7317-7322, where cumfnc stops. Rust checks it after the
// status checks of cdffnc, before cumfnc is called.
#[inline]
fn check_cumfnc_dfn(dfn: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if dfn < 1.0 {
        return Err(FisherSnedecorNoncentralError::DfnTooSmall(dfn));
    }
    Ok(())
}

// cdflib.f90:7324-7329, where cumfnc stops. Rust checks it after the
// status checks of cdffnc, before cumfnc is called.
#[inline]
fn check_cumfnc_dfd(dfd: f64) -> Result<(), FisherSnedecorNoncentralError> {
    if dfd < 1.0 {
        return Err(FisherSnedecorNoncentralError::DfdTooSmall(dfd));
    }
    Ok(())
}

impl FisherSnedecorNoncentral {
    /// Construct a noncentral *F*(*dfn*, *dfd*, *λ*) distribution.
    ///
    /// # Panics
    ///
    /// Panics if any argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(dfn: f64, dfd: f64, ncp: f64) -> Self {
        Self::try_new(dfn, dfd, ncp).unwrap()
    }

    /// Fallible counterpart of [`new`] returning a
    /// [`FisherSnedecorNoncentralError`] instead of panicking.
    ///
    /// Returns [`DfnNotPositive`], [`DfnNotFinite`], [`DfdNotPositive`],
    /// [`DfdNotFinite`], [`NcpNegative`], or [`NcpNotFinite`] if an argument
    /// fails its validity check, and then [`DfnTooSmall`] or [`DfdTooSmall`]
    /// if *dfn* < 1 or *dfd* < 1, where `cumfnc` stops with a fatal error.
    ///
    /// [`DfnNotPositive`]: FisherSnedecorNoncentralError::DfnNotPositive
    /// [`DfnNotFinite`]: FisherSnedecorNoncentralError::DfnNotFinite
    /// [`DfdNotPositive`]: FisherSnedecorNoncentralError::DfdNotPositive
    /// [`DfdNotFinite`]: FisherSnedecorNoncentralError::DfdNotFinite
    /// [`NcpNegative`]: FisherSnedecorNoncentralError::NcpNegative
    /// [`NcpNotFinite`]: FisherSnedecorNoncentralError::NcpNotFinite
    /// [`DfnTooSmall`]: FisherSnedecorNoncentralError::DfnTooSmall
    /// [`DfdTooSmall`]: FisherSnedecorNoncentralError::DfdTooSmall
    /// [`new`]: Self::new
    #[inline]
    pub fn try_new(dfn: f64, dfd: f64, ncp: f64) -> Result<Self, FisherSnedecorNoncentralError> {
        check_dfn(dfn)?;
        check_dfd(dfd)?;
        check_pnonc(ncp)?;
        check_cumfnc_dfn(dfn)?;
        check_cumfnc_dfd(dfd)?;
        Ok(Self { dfn, dfd, ncp })
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

    /// Returns the noncentrality parameter *λ*.
    #[inline]
    pub const fn ncp(&self) -> f64 {
        self.ncp
    }

    /// Returns the numerator degrees of freedom *dfn* satisfying
    /// Pr[*X* ≤ *f*] = *p*, searched for in [1 . . 10³⁰].
    ///
    /// CDFLIB's `cdffnc` with `which = 3`. As in CDFLIB, the lower bound
    /// reported on failure is 0, not 1. The precision limits of [`cdf`]
    /// apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    pub fn search_dfn(
        p: f64,
        f: f64,
        dfd: f64,
        ncp: f64,
    ) -> Result<f64, FisherSnedecorNoncentralError> {
        check_p(p)?;
        check_f(f)?;
        check_dfd(dfd)?;
        check_pnonc(ncp)?;
        check_cumfnc_dfd(dfd)?;
        let pnonc = ncp;

        // cdflib.f90:4738-4771
        let mut d = dstinv(1.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let dfn = d.x();
            let (cum, _ccum) = cumfnc(f, dfn, dfd, pnonc);
            let fx = cum - p;
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

    /// Returns the denominator degrees of freedom *dfd* satisfying
    /// Pr[*X* ≤ *f*] = *p*, searched for in [1 . . 10³⁰].
    ///
    /// CDFLIB's `cdffnc` with `which = 4`. As in CDFLIB, the lower bound
    /// reported on failure is 0, not 1. The precision limits of [`cdf`]
    /// apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    pub fn search_dfd(
        p: f64,
        f: f64,
        dfn: f64,
        ncp: f64,
    ) -> Result<f64, FisherSnedecorNoncentralError> {
        check_p(p)?;
        check_f(f)?;
        check_dfn(dfn)?;
        check_pnonc(ncp)?;
        check_cumfnc_dfn(dfn)?;
        let pnonc = ncp;

        // cdflib.f90:4777-4809
        let mut d = dstinv(1.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let dfd = d.x();
            let (cum, _ccum) = cumfnc(f, dfn, dfd, pnonc);
            let fx = cum - p;
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

    /// Returns the noncentrality parameter *λ* satisfying
    /// Pr[*X* ≤ *f*] = *p*, searched for in [0 . . 10⁴].
    ///
    /// CDFLIB's `cdffnc` with `which = 5`. At *p* = 0, where the answer is
    /// +∞, the search stops, as the F90 does, where the computed probability
    /// becomes 0, and returns that finite value:
    /// `search_ncp(0.0, 2.0, 5.0, 10.0)` returns about 9770. The precision
    /// limits of [`cdf`] apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    pub fn search_ncp(
        p: f64,
        f: f64,
        dfn: f64,
        dfd: f64,
    ) -> Result<f64, FisherSnedecorNoncentralError> {
        check_p(p)?;
        check_f(f)?;
        check_dfn(dfn)?;
        check_dfd(dfd)?;
        check_cumfnc_dfn(dfn)?;
        check_cumfnc_dfd(dfd)?;

        // cdflib.f90:4815-4848
        let mut d = dstinv(0.0, TENT4, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let pnonc = d.x();
            let (cum, _ccum) = cumfnc(f, dfn, dfd, pnonc);
            let fx = cum - p;
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(if d.qleft() {
                SearchError::AnswerBelowLowerBound { bound: 0.0 }
            } else {
                SearchError::AnswerAboveUpperBound { bound: TENT4 }
            }
            .into());
        }
        Ok(d.x())
    }
}

impl ContinuousCdf for FisherSnedecorNoncentral {
    type Error = FisherSnedecorNoncentralError;

    /// CDFLIB's `cdffnc` with `which = 1`.
    ///
    /// The series stops when a term is less than 10⁻⁴ times the sum. For a
    /// large *λ*, whose Poisson weights spread over many terms, this happens
    /// before the sum is complete, and, as in the F90, the result is badly
    /// wrong: with *dfn* = 10 and *dfd* = 1000 the true cdf at twice the mean
    /// is close to 1, but the computed one is 0.9992 for *λ* = 10³, 0.988 for
    /// *λ* = 10⁵ and 0.51 for *λ* = 10⁸. For very large degrees of freedom
    /// (*dfd* beyond about 10¹⁴, for example) the terms lose their digits, so
    /// that the result can fall outside [0 . . 1]. For a small *λ* the sum can
    /// also exceed 1 by a few ulps in the right tail, so that [`ccdf`] is
    /// slightly negative: with *dfn* = 3.7, *dfd* = 30 and *λ* = 10⁻⁹, the
    /// cdf at 100 is 1 + 4.4 · 10⁻¹⁶.
    ///
    /// # Panics
    ///
    /// Panics if *λ*/2 ≥ 2³¹ − 1, or if the forward sum of `cumfnc` runs
    /// past index 2³¹ − 1, where CDFLIB's default integers overflow; the
    /// latter needs *λ*/2 within a few hundred thousand of 2³¹. Panics
    /// also when the degrees of freedom are so large that the sum of the
    /// series is NaN, where the F90 never returns.
    ///
    /// [`ccdf`]: ContinuousCdf::ccdf
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: exact endpoint at +inf, where cumfnc truncates its sum short of 1.
        if x == f64::INFINITY {
            return 1.0;
        }
        // Rust only: no status -4 for f < 0 (cdflib.f90:4634-4646); cumfnc
        // returns (0, 1) there.
        // cdflib.f90:4691
        cumfnc(x, self.dfn, self.dfd, self.ncp).0
    }

    /// CDFLIB's `cdffnc` with `which = 1`.
    ///
    /// Unlike the central distributions, CDFLIB computes this as
    /// 1 − [`cdf`] (except for *λ* < 10⁻¹⁰, where it uses the central *F*),
    /// from a series that stops when a term is less than 10⁻⁴ times the
    /// sum, so the result has no relative precision in the right tail. The
    /// precision limits of [`cdf`] apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: exact endpoint at +inf, where cumfnc truncates its sum short of 1.
        if x == f64::INFINITY {
            return 0.0;
        }
        // Rust only: no status -4 for f < 0 (cdflib.f90:4634-4646); cumfnc
        // returns (0, 1) there.
        // cdflib.f90:4691
        cumfnc(x, self.dfn, self.dfd, self.ncp).1
    }

    /// CDFLIB's `cdffnc` with `which = 2`, searched for in [0 . . 10³⁰].
    ///
    /// The precision limits of [`cdf`] apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, FisherSnedecorNoncentralError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let dfn = self.dfn;
        let dfd = self.dfd;
        let pnonc = self.ncp;

        // cdflib.f90:4698-4732
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let f = d.x();
            let (cum, _ccum) = cumfnc(f, dfn, dfd, pnonc);
            let fx = cum - p;
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

impl Mean for FisherSnedecorNoncentral {
    /// Defined for *dfd* > 2.
    #[inline]
    fn mean(&self) -> f64 {
        if self.dfd > 2.0 {
            // dfd (dfn + λ) / (dfn (dfd - 2)), written as a product of ratios
            // so that no intermediate overflows.
            self.dfd / (self.dfd - 2.0) * ((self.dfn + self.ncp) / self.dfn)
        } else {
            f64::NAN
        }
    }
}

impl Variance for FisherSnedecorNoncentral {
    /// Defined for *dfd* > 4.
    #[inline]
    fn variance(&self) -> f64 {
        let dfn = self.dfn;
        let dfd = self.dfd;
        let ncp = self.ncp;
        if dfd > 4.0 {
            // 2 dfd² ((dfn + λ)² + (dfd - 2)(dfn + 2λ))
            //     / (dfn² (dfd - 2)² (dfd - 4)),
            // written with m = dfd / (dfd - 2), r = (dfn + λ) / dfn and
            // t = (dfn + 2λ) / dfn² so that no intermediate overflows.
            let m = dfd / (dfd - 2.0);
            let r = (dfn + ncp) / dfn;
            let t = (dfn + 2.0 * ncp) / dfn / dfn;
            2.0 * m * m * (r * r / (dfd - 4.0) + (dfd - 2.0) / (dfd - 4.0) * t)
        } else {
            f64::NAN
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_inputs() {
        assert!(matches!(
            FisherSnedecorNoncentral::try_new(0.0, 5.0, 1.0),
            Err(FisherSnedecorNoncentralError::DfnNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::try_new(0.5, 5.0, 1.0),
            Err(FisherSnedecorNoncentralError::DfnTooSmall(0.5))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::try_new(5.0, 0.0, 1.0),
            Err(FisherSnedecorNoncentralError::DfdNotPositive(0.0))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::try_new(5.0, 0.5, 1.0),
            Err(FisherSnedecorNoncentralError::DfdTooSmall(0.5))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::try_new(5.0, 5.0, -1.0),
            Err(FisherSnedecorNoncentralError::NcpNegative(-1.0))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::search_ncp(-0.1, 1.0, 5.0, 10.0),
            Err(FisherSnedecorNoncentralError::PNotInRange(-0.1))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::search_dfd(0.5, 1.0, 0.5, 1.0),
            Err(FisherSnedecorNoncentralError::DfnTooSmall(0.5))
        ));
        assert!(matches!(
            FisherSnedecorNoncentral::search_ncp(0.5, 1.0, 5.0, 0.5),
            Err(FisherSnedecorNoncentralError::DfdTooSmall(0.5))
        ));
    }

    #[test]
    fn inverse_and_moment_edges() {
        let d = FisherSnedecorNoncentral::new(5.0, 10.0, 2.0);
        assert_eq!(d.inverse_cdf(0.0).unwrap(), 0.0);
        assert!(d.inverse_cdf(0.25).unwrap().is_finite());
        assert!(d.mean().is_finite());
        assert!(d.variance().is_finite());
        assert!(FisherSnedecorNoncentral::new(5.0, 2.0, 2.0).mean().is_nan());
        assert!(FisherSnedecorNoncentral::new(5.0, 4.0, 2.0)
            .variance()
            .is_nan());
    }

    #[test]
    fn central_reduction_path_is_consistent() {
        let d = FisherSnedecorNoncentral::new(5.0, 10.0, 0.0);
        let x = 1.5;
        let cdf = d.cdf(x);
        let ccdf = d.ccdf(x);
        assert!((cdf + ccdf - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dfn_below_1_is_checked_after_the_status_checks() {
        // F90 returns status -6 for dfd = -1 before cumfnc can stop on
        // dfn = 0.5.
        assert_eq!(
            FisherSnedecorNoncentral::search_ncp(0.5, 1.0, 0.5, -1.0),
            Err(FisherSnedecorNoncentralError::DfdNotPositive(-1.0))
        );
        assert_eq!(
            FisherSnedecorNoncentral::search_ncp(0.5, 1.0, 0.5, 2.0),
            Err(FisherSnedecorNoncentralError::DfnTooSmall(0.5))
        );
        assert_eq!(
            FisherSnedecorNoncentral::try_new(0.5, 2.0, -1.0),
            Err(FisherSnedecorNoncentralError::NcpNegative(-1.0))
        );
    }

    #[test]
    fn nan_f_gives_nan() {
        let d = FisherSnedecorNoncentral::new(3.0, 4.0, 2.0);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
    }

    #[test]
    #[should_panic(expected = "integer overflow")]
    fn huge_ncp_overflows_the_f90_integers() {
        // pnonc / 2 >= 2^31 - 1: int(xnonc) overflows in F90.
        FisherSnedecorNoncentral::new(3.0, 4.0, 5.0e9).cdf(1.0);
    }

    #[test]
    fn moments_do_not_overflow() {
        // For large dfd the mean tends to (dfn + λ) / dfn and the variance
        // to 2 (dfn + 2λ) / dfn² + ((dfn + λ) / dfn)² · 2 / dfd.
        let d = FisherSnedecorNoncentral::new(1.0, f64::MAX, 1.0);
        assert_eq!(d.mean(), 2.0);
        assert_eq!(d.variance(), 6.0);
    }
}
