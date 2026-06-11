use crate::error::SearchError;
use crate::search::{search_monotone_with_atol, SEARCH_BOUND};
use crate::special::gamma_inc;
use crate::special::gamma_log;
use crate::traits::{ContinuousCdf, Mean, Variance};
use thiserror::Error;

/// Noncentral χ² distribution with *df* > 0 degrees of freedom and
/// noncentrality parameter *λ* ≥ 0.
///
/// # Notes
///
/// Neither [`Continuous`] nor [`Entropy`] is implemented.
///
/// [`Continuous`]: crate::traits::Continuous
/// [`Entropy`]: crate::traits::Entropy
///
/// # Example
///
/// ```
/// use cdflib::ChiSquaredNoncentral;
/// use cdflib::traits::ContinuousCdf;
///
/// let d = ChiSquaredNoncentral::new(5.0, 10.0);
///
/// // Probability of observing a value ≤ 15.0
/// let p = d.cdf(15.0);
///
/// // Compute noncentrality *λ* given Pr[X ≤ 15] = 0.5 and df = 5
/// let ncp = ChiSquaredNoncentral::search_ncp(0.5, 15.0, 5.0).unwrap();
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChiSquaredNoncentral {
    df: f64,
    ncp: f64,
}

/// Errors arising from constructing a [`ChiSquaredNoncentral`] or from its
/// parameter searches.
///
/// [`ChiSquaredNoncentral`]: crate::ChiSquaredNoncentral
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum ChiSquaredNoncentralError {
    /// The degrees of freedom *df* was not strictly positive.
    #[error("degrees of freedom must be positive, got {0}")]
    DfNotPositive(f64),
    /// The degrees of freedom *df* was not finite.
    #[error("degrees of freedom must be finite, got {0}")]
    DfNotFinite(f64),
    /// The noncentrality parameter *λ* was negative.
    #[error("noncentrality parameter must be ≥ 0, got {0}")]
    NcpNegative(f64),
    /// The noncentrality parameter *λ* was not finite.
    #[error("noncentrality parameter must be finite, got {0}")]
    NcpNotFinite(f64),
    /// The argument *x* was not strictly positive.
    #[error("argument x must be positive, got {0}")]
    XNotPositive(f64),
    /// The argument *x* was not finite.
    #[error("argument x must be finite, got {0}")]
    XNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (or was non-finite).
    #[error("probability {0} outside [0..1]")]
    QNotInRange(f64),
    /// The internal root-finder failed; see [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

impl ChiSquaredNoncentral {
    /// Construct a noncentral χ²(*df*, *λ*) distribution with *df* > 0
    /// degrees of freedom and noncentrality *λ* ≥ 0.
    ///
    /// # Panics
    ///
    /// Panics if either argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(df: f64, ncp: f64) -> Self {
        Self::try_new(df, ncp).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a
    /// [`ChiSquaredNoncentralError`] instead of panicking.
    ///
    /// Returns [`DfNotPositive`], [`DfNotFinite`], [`NcpNegative`], or
    /// [`NcpNotFinite`] if either argument fails its validity check.
    ///
    /// [`DfNotPositive`]: ChiSquaredNoncentralError::DfNotPositive
    /// [`DfNotFinite`]: ChiSquaredNoncentralError::DfNotFinite
    /// [`NcpNegative`]: ChiSquaredNoncentralError::NcpNegative
    /// [`NcpNotFinite`]: ChiSquaredNoncentralError::NcpNotFinite
    #[inline]
    pub fn try_new(df: f64, ncp: f64) -> Result<Self, ChiSquaredNoncentralError> {
        if !df.is_finite() {
            return Err(ChiSquaredNoncentralError::DfNotFinite(df));
        }
        if df <= 0.0 {
            return Err(ChiSquaredNoncentralError::DfNotPositive(df));
        }
        if !ncp.is_finite() {
            return Err(ChiSquaredNoncentralError::NcpNotFinite(ncp));
        }
        if ncp < 0.0 {
            return Err(ChiSquaredNoncentralError::NcpNegative(ncp));
        }
        Ok(Self { df, ncp })
    }

    /// Returns the degrees of freedom *df*.
    #[inline]
    pub const fn df(&self) -> f64 {
        self.df
    }

    /// Returns the noncentrality parameter *λ*.
    #[inline]
    pub const fn ncp(&self) -> f64 {
        self.ncp
    }

    /// Returns the degrees of freedom *df* satisfying Pr[*X* ≤ *x*] = *p*
    /// given *λ*. Mirrors CDFLIB's `cdfchn` with `which = 3`.
    ///
    /// Unlike most `cdf*` searches, this one does not take *q*: CDFLIB
    /// (cdflib.f90:3685) documents *q* as "generally not used by this
    /// subroutine and is only included for similarity with other routines",
    /// so it is dropped from the Rust surface.
    #[inline]
    pub fn search_df(p: f64, x: f64, ncp: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        if !x.is_finite() {
            return Err(ChiSquaredNoncentralError::XNotFinite(x));
        }
        if x <= 0.0 {
            return Err(ChiSquaredNoncentralError::XNotPositive(x));
        }
        if !ncp.is_finite() {
            return Err(ChiSquaredNoncentralError::NcpNotFinite(ncp));
        }
        if ncp < 0.0 {
            return Err(ChiSquaredNoncentralError::NcpNegative(ncp));
        }
        let f = |df: f64| cumchn(x, df, ncp).0 - p;
        // Match cdfchn's which=3: range (zero, inf), start = 5.0, atol = 1e-50.
        Ok(search_monotone_with_atol(
            0.0,
            SEARCH_BOUND,
            5.0,
            0.0,
            SEARCH_BOUND,
            1.0e-50,
            f,
        )?)
    }

    /// Returns the noncentrality *λ* satisfying Pr[*X* ≤ *x*] = *p* given *df*.
    ///
    /// Mirrors CDFLIB's `cdfchn` with `which = 4`. The search runs over
    /// (0, 10⁴] because `cumchn`'s iteration cost grows with *λ*. As in
    /// [`search_df`](Self::search_df), *q* is dropped from the Rust surface.
    #[inline]
    pub fn search_ncp(p: f64, x: f64, df: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        if !x.is_finite() {
            return Err(ChiSquaredNoncentralError::XNotFinite(x));
        }
        if x <= 0.0 {
            return Err(ChiSquaredNoncentralError::XNotPositive(x));
        }
        if !df.is_finite() {
            return Err(ChiSquaredNoncentralError::DfNotFinite(df));
        }
        if df <= 0.0 {
            return Err(ChiSquaredNoncentralError::DfNotPositive(df));
        }
        let f = |ncp: f64| cumchn(x, df, ncp).0 - p;
        // CDF is decreasing in ncp (shift right). Upper bound 1e4
        // matches CDFLIB's hard cap; cumchn's iteration cost grows
        // with ncp so unbounded searches are intentionally avoided.
        // atol = 1e-50 per cdfchn (cdflib.f90:3719).
        Ok(search_monotone_with_atol(
            0.0, 1.0e4, 5.0, 0.0, 1.0e4, 1.0e-50, f,
        )?)
    }
}

#[inline]
fn check_p(p: f64) -> Result<(), ChiSquaredNoncentralError> {
    if !(0.0..=1.0).contains(&p) || !p.is_finite() {
        Err(ChiSquaredNoncentralError::PNotInRange(p))
    } else {
        Ok(())
    }
}

/// `cumchn`: noncentral χ² CDF and SF.
fn cumchn(x: f64, df: f64, pnonc: f64) -> (f64, f64) {
    if x.is_nan() || df.is_nan() || pnonc.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    if x <= 0.0 {
        return (0.0, 1.0);
    }
    if pnonc <= 1e-10 {
        let (p, q) = gamma_inc(df / 2.0, x / 2.0);
        return (p, q);
    }

    let eps = 1e-5;
    let ntired: i32 = 1000;
    let xnonc = pnonc / 2.0;
    let mut icent = xnonc as i32;
    if icent == 0 {
        icent = 1;
    }
    let chid2 = x / 2.0;

    let lfact = gamma_log((icent + 1) as f64);
    let lcntwt = -xnonc + (icent as f64) * xnonc.ln() - lfact;
    let centwt = lcntwt.exp();

    let dg = |i: i32| df + 2.0 * (i as f64);
    let (pcent, _) = gamma_inc(dg(icent) / 2.0, chid2);

    let dfd2 = dg(icent) / 2.0;
    let lfact = gamma_log(1.0 + dfd2);
    let lcntaj = dfd2 * chid2.ln() - chid2 - lfact;
    let centaj = lcntaj.exp();
    let mut sum = centwt * pcent;

    // Sum backwards.
    let mut iterb: i32 = 0;
    let mut sumadj = 0.0;
    let mut adj = centaj;
    let mut wt = centwt;
    let mut i = icent;
    loop {
        let dfd2 = dg(i) / 2.0;
        adj = adj * dfd2 / chid2;
        sumadj += adj;
        let pterm = pcent + sumadj;
        wt *= i as f64 / xnonc;
        let term = wt * pterm;
        sum += term;
        i -= 1;
        iterb += 1;
        let small = sum < 1e-20 || term < eps * sum;
        if iterb > ntired || small || i == 0 {
            break;
        }
    }

    // Sum forwards.
    let mut iterf: i32 = 0;
    let mut adj = centaj;
    let mut sumadj = centaj;
    let mut wt = centwt;
    let mut i = icent;
    loop {
        wt *= xnonc / (i + 1) as f64;
        let pterm = pcent - sumadj;
        let term = wt * pterm;
        sum += term;
        i += 1;
        let dfd2 = dg(i) / 2.0;
        adj = adj * chid2 / dfd2;
        sumadj += adj;
        iterf += 1;
        let small = sum < 1e-20 || term < eps * sum;
        if iterf > ntired || small {
            break;
        }
    }

    let cum = sum;
    (cum, 0.5 + (0.5 - cum))
}

impl ContinuousCdf for ChiSquaredNoncentral {
    type Error = ChiSquaredNoncentralError;

    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        cumchn(x, self.df, self.ncp).0
    }
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        cumchn(x, self.df, self.ncp).1
    }
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let df = self.df;
        let ncp = self.ncp;
        let f = |x: f64| cumchn(x, df, ncp).0 - p;
        // Match cdfchn's which=2: range (0, inf), start = 5.0, atol = 1e-50.
        Ok(search_monotone_with_atol(
            0.0,
            SEARCH_BOUND,
            5.0,
            0.0,
            SEARCH_BOUND,
            1.0e-50,
            f,
        )?)
    }
}

impl Mean for ChiSquaredNoncentral {
    #[inline]
    fn mean(&self) -> f64 {
        self.df + self.ncp
    }
}

impl Variance for ChiSquaredNoncentral {
    #[inline]
    fn variance(&self) -> f64 {
        2.0 * (self.df + 2.0 * self.ncp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_invalid_inputs() {
        assert!(matches!(
            ChiSquaredNoncentral::try_new(0.0, 1.0),
            Err(ChiSquaredNoncentralError::DfNotPositive(0.0))
        ));
        assert!(matches!(
            ChiSquaredNoncentral::try_new(1.0, -1.0),
            Err(ChiSquaredNoncentralError::NcpNegative(-1.0))
        ));
        assert!(matches!(
            ChiSquaredNoncentral::search_df(-0.1, 1.0, 2.0),
            Err(ChiSquaredNoncentralError::PNotInRange(-0.1))
        ));
    }

    #[test]
    fn inverse_and_moment_edges() {
        let d = ChiSquaredNoncentral::new(5.0, 2.0);
        assert_eq!(d.inverse_cdf(0.0).unwrap(), 0.0);
        assert!(d.inverse_cdf(0.25).unwrap().is_finite());
        assert!(d.mean().is_finite());
        assert!(d.variance().is_finite());
    }

    #[test]
    fn central_limit_path_is_consistent() {
        let d = ChiSquaredNoncentral::new(4.0, 0.0);
        let x = 3.0;
        let cdf = d.cdf(x);
        let ccdf = d.ccdf(x);
        assert!((cdf + ccdf - 1.0).abs() < 1e-12);
    }
}
