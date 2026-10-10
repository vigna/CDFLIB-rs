use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::gamma_log;
use crate::traits::{ContinuousCdf, Mean, Variance};
use thiserror::Error;

use super::chi_squared::cumchi;

// Parameters of cdfchn (cdflib.f90:3801-3815).
const ATOL: f64 = 1.0e-50;
const INF: f64 = 1.0e300;
const TENT4: f64 = 1.0e4;
const TOL: f64 = 1.0e-8;

/// Noncentral χ² distribution with *df* > 0 degrees of freedom and
/// noncentrality parameter *λ* ≥ 0.
///
/// The methods correspond to CDFLIB's `cdfchn` (cdflib.f90:3683), whose
/// PNONC is *λ*: `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is
/// [`inverse_cdf`], `which = 3` is [`search_df`], `which = 4` is
/// [`search_ncp`]. CDFLIB's searches use *p* only.
///
/// # Notes
///
/// Neither [`Continuous`] nor [`Entropy`] is implemented.
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
/// assert!((p - 0.5531405533).abs() < 1e-6);
///
/// // Compute noncentrality λ given Pr[X ≤ 15] = 0.5 and df = 5
/// let ncp = ChiSquaredNoncentral::search_ncp(0.5, 15.0, 5.0).unwrap();
/// assert!((ncp - 10.946).abs() < 1e-3);
/// ```
///
/// [`Continuous`]: crate::traits::Continuous
/// [`Entropy`]: crate::traits::Entropy
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`search_df`]: ChiSquaredNoncentral::search_df
/// [`search_ncp`]: ChiSquaredNoncentral::search_ncp
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChiSquaredNoncentral {
    df: f64,
    ncp: f64,
}

/// Errors arising from constructing a [`ChiSquaredNoncentral`] or from its
/// parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfchn`.
///
/// [`ChiSquaredNoncentral`]: crate::ChiSquaredNoncentral
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum ChiSquaredNoncentralError {
    /// The degrees of freedom *df* was not strictly positive (`cdfchn`
    /// status −5). [`search_df`] also returns it, checked only in Rust,
    /// when the *df* it computes is 0.
    ///
    /// [`search_df`]: crate::ChiSquaredNoncentral::search_df
    #[error("degrees of freedom must be positive, got {0}")]
    DfNotPositive(f64),
    /// The degrees of freedom *df* was not finite (checked only in Rust).
    #[error("degrees of freedom must be finite, got {0}")]
    DfNotFinite(f64),
    /// The noncentrality parameter *λ* was negative (`cdfchn` status −6).
    #[error("noncentrality parameter must be nonnegative, got {0}")]
    NcpNegative(f64),
    /// The noncentrality parameter *λ* was not finite (checked only in
    /// Rust).
    #[error("noncentrality parameter must be finite, got {0}")]
    NcpNotFinite(f64),
    /// The argument *x* was not strictly positive (`cdfchn` status −4 for
    /// *x* < 0; Rust also rejects *x* = 0, where Pr[*X* ≤ *x*] = 0 whatever
    /// the parameter searched for).
    #[error("argument x must be positive, got {0}")]
    XNotPositive(f64),
    /// The argument *x* was not finite (checked only in Rust).
    #[error("argument x must be finite, got {0}")]
    XNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdfchn` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1]; NaN is also rejected.
    /// No method of [`ChiSquaredNoncentral`] returns it, since CDFLIB's
    /// `cdfchn` does not use *q*.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The search for the answer failed (`cdfchn` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
}

// Rust only: cdfchn has no status for an error value of gamma_inc, which
// cumchn passes on as a probability. It occurs only when the degrees of
// freedom passed to cumchi (df, or df + 2 icent) exceed about 1.3e29 and x
// is within a few ulps of them.
#[inline]
fn cumchi_or_panic(x: f64, df: f64) -> (f64, f64) {
    match cumchi(x, df) {
        Ok(r) => r,
        Err(e) => panic!("cumchi({x}, {df}): {e}"),
    }
}

/// Returns the cumulative noncentral χ² distribution (*cum*, *ccum*) at *x*
/// with *df* degrees of freedom and noncentrality parameter *pnonc*
/// (`cumchn`, cdflib.f90:6910).
///
/// The series is summed backwards and forwards from the central term, the
/// one with the greatest Poisson weight; each sum stops when a term is
/// less than `EPS` times the sum or after `NTIRED` terms.
#[allow(clippy::assign_op_pattern)]
pub(crate) fn cumchn(x: f64, df: f64, pnonc: f64) -> (f64, f64) {
    const EPS: f64 = 0.00001;
    const NTIRED: i32 = 1000;
    // F90 statement functions (cdflib.f90:7002-7003). The F90 qsmall reads
    // sum1 from the enclosing scope; the closure takes it as an argument.
    let qsmall = |sum1: f64, xx: f64| sum1 < 1.0e-20 || xx < EPS * sum1;
    let dg = |i: i32| df + 2.0 * i as f64;

    if x <= 0.0 {
        return (0.0, 1.0);
    }
    // When the noncentrality parameter is (essentially) zero, use the
    // cumulative χ² distribution.
    if pnonc <= 1.0e-10 {
        return cumchi_or_panic(x, df);
    }

    let xnonc = pnonc / 2.0;
    // Rust only: int(xnonc) and icent + 1 below overflow the default
    // integer, which is undefined in Fortran; a NaN pnonc is caught here
    // too.
    if xnonc.is_nan() || f64::from(i32::MAX) <= xnonc {
        panic!("cumchn: integer overflow for pnonc = {pnonc}");
    }
    // Weight, χ² and adjustment term of the central term of the series,
    // the one in which the Poisson weight is greatest. The adjustment term
    // is the amount that must be subtracted from the χ² to move up two
    // degrees of freedom.
    let mut icent = xnonc as i32;
    if icent == 0 {
        icent = 1;
    }

    // Rust only: for the smallest subnormal x, chid2 = x/2 is 0, and the
    // first backward step below computes adj = 0 * dfd2 / 0 = NaN, as in
    // the F90. Sum the series at the next floating-point number instead,
    // whose half is the smallest subnormal: the terms are powers
    // (x/2)^(df/2 + j), which change by a factor of at most 2^(df/2 + j)
    // and are subnormal unless df is tiny.
    let x = if x / 2.0 == 0.0 {
        2.0 * f64::from_bits(1)
    } else {
        x
    };
    let chid2 = x / 2.0;
    // Central weight term.
    let lfact = gamma_log((icent + 1) as f64);
    let lcntwt = -xnonc + icent as f64 * xnonc.ln() - lfact;
    let centwt = lcntwt.exp();
    // Central χ².
    let (pcent, _ccum) = cumchi_or_panic(x, dg(icent));
    // Central adjustment term.
    let dfd2 = dg(icent) / 2.0;
    let lfact = gamma_log(1.0 + dfd2);
    let lcntaj = dfd2 * chid2.ln() - chid2 - lfact;
    let centaj = lcntaj.exp();
    let mut sum1 = centwt * pcent;

    // Sum backwards from the central term towards zero. Quit whenever
    // either the zero term is reached, or the term gets small relative to
    // the sum, or more than NTIRED terms are totaled.
    let mut iterb = 0;
    let mut sumadj = 0.0;
    let mut adj = centaj;
    let mut wt = centwt;
    let mut i = icent;
    #[allow(unused_assignments)]
    let mut term = 0.0;

    loop {
        let dfd2 = dg(i) / 2.0;
        // Adjust the χ² for two fewer degrees of freedom; the adjusted
        // value ends up in pterm.
        adj = adj * dfd2 / chid2;
        sumadj = sumadj + adj;
        let pterm = pcent + sumadj;
        // Adjust the Poisson weight for j decreased by one.
        wt = wt * (i as f64 / xnonc);
        term = wt * pterm;
        sum1 = sum1 + term;
        i = i - 1;
        iterb = iterb + 1;

        if NTIRED < iterb || qsmall(sum1, term) || i == 0 {
            break;
        }
    }

    let mut iterf = 0;
    // Now sum forward from the central term towards infinity. Quit when
    // either the term gets small relative to the sum, or more than NTIRED
    // terms are totaled.
    let mut sumadj = centaj;
    let mut adj = centaj;
    let mut wt = centwt;
    let mut i = icent;
    // Update weights for the next higher j.
    loop {
        // Rust only: i + 1 overflows the default integer, which is
        // undefined in Fortran.
        if i == i32::MAX {
            panic!("cumchn: integer overflow for pnonc = {pnonc}");
        }
        wt = wt * (xnonc / (i + 1) as f64);
        // Calculate pterm and add the term to the sum.
        let pterm = pcent - sumadj;
        term = wt * pterm;
        sum1 = sum1 + term;
        // Update the adjustment term for df for the next iteration.
        i = i + 1;
        let dfd2 = dg(i) / 2.0;
        adj = adj * chid2 / dfd2;
        sumadj = sumadj + adj;
        iterf = iterf + 1;

        if NTIRED < iterf || qsmall(sum1, term) {
            break;
        }
    }

    let cum = sum1;
    let ccum = 0.5 + (0.5 - cum);
    (cum, ccum)
}

// cdflib.f90:3843-3862 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), ChiSquaredNoncentralError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(ChiSquaredNoncentralError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:3863-3875 (status -4). Rust also rejects x = 0 and a
// non-finite x.
#[inline]
fn check_x(x: f64) -> Result<(), ChiSquaredNoncentralError> {
    if x <= 0.0 {
        return Err(ChiSquaredNoncentralError::XNotPositive(x));
    }
    if !x.is_finite() {
        return Err(ChiSquaredNoncentralError::XNotFinite(x));
    }
    Ok(())
}

// cdflib.f90:3876-3888 (status -5). Rust also rejects a non-finite df.
#[inline]
fn check_df(df: f64) -> Result<(), ChiSquaredNoncentralError> {
    if df <= 0.0 {
        return Err(ChiSquaredNoncentralError::DfNotPositive(df));
    }
    if !df.is_finite() {
        return Err(ChiSquaredNoncentralError::DfNotFinite(df));
    }
    Ok(())
}

// cdflib.f90:3889-3901 (status -6). Rust also rejects a non-finite pnonc.
#[inline]
fn check_pnonc(pnonc: f64) -> Result<(), ChiSquaredNoncentralError> {
    if pnonc < 0.0 {
        return Err(ChiSquaredNoncentralError::NcpNegative(pnonc));
    }
    if !pnonc.is_finite() {
        return Err(ChiSquaredNoncentralError::NcpNotFinite(pnonc));
    }
    Ok(())
}

impl ChiSquaredNoncentral {
    /// Construct a noncentral χ²(*df*, *λ*) distribution.
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
        check_df(df)?;
        check_pnonc(ncp)?;
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

    /// Returns the degrees of freedom *df* satisfying Pr[*X* ≤ *x*] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfchn` with `which = 3`. A computed *df* of 0, at the
    /// lower end of the search interval, is reported as [`DfNotPositive`].
    /// At *p* = 0, where the answer is +∞, the search stops, as the F90 does,
    /// where the computed probability becomes 0, and returns that finite
    /// value: `search_df(0.0, 2.0, 1.0)` returns about 395. The precision
    /// limits of [`cdf`] apply.
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does. Since the search evaluates its upper bound
    /// 10³⁰⁰, this happens in particular when *x* is within a few ulps of
    /// 10³⁰⁰.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    /// [`DfNotPositive`]: ChiSquaredNoncentralError::DfNotPositive
    #[inline]
    pub fn search_df(p: f64, x: f64, ncp: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        check_x(x)?;
        check_pnonc(ncp)?;
        let pnonc = ncp;

        // cdflib.f90:3952-3984
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let df = d.x();
            let (cum, _ccum) = cumchn(x, df, pnonc);
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
        // Rust only: the F90 returns the answer whatever its value; at the
        // lower end of the search interval it is not a valid df.
        let df = d.x();
        check_df(df)?;
        Ok(df)
    }

    /// Returns the noncentrality parameter *λ* satisfying
    /// Pr[*X* ≤ *x*] = *p*, searched for in [0 . . 10⁴]: the cost of
    /// `cumchn` grows with *λ*, which is why CDFLIB bounds the search.
    ///
    /// CDFLIB's `cdfchn` with `which = 4`.
    ///
    /// At *p* = 0, where the answer is +∞, the search stops, as the F90 does,
    /// where the computed probability becomes 0, and returns that finite
    /// value: `search_ncp(0.0, 2.0, 1.0)` returns about 395.
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` fails inside `cumchn`, which needs *df*
    /// beyond about 1.3 · 10²⁹ and *x* within a few ulps of it.
    #[inline]
    pub fn search_ncp(p: f64, x: f64, df: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        check_x(x)?;
        check_df(df)?;

        // cdflib.f90:3990-4022
        let mut d = dstinv(0.0, TENT4, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let pnonc = d.x();
            let (cum, _ccum) = cumchn(x, df, pnonc);
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

impl ContinuousCdf for ChiSquaredNoncentral {
    type Error = ChiSquaredNoncentralError;

    /// CDFLIB's `cdfchn` with `which = 1`.
    ///
    /// The series sums at most about 1000 terms on each side of its central
    /// term, and stops earlier when a term is less than 10⁻⁵ times the sum.
    /// For a large *λ* these terms do not cover the Poisson weights, and, as
    /// in the F90, the result is badly wrong: with *df* = 1 the true
    /// cdf(*df* + *λ*) is close to 0.5, but the computed one is 0.42 for
    /// *λ* = 10⁶ and 0.056 for *λ* = 10⁸. For very large *df* (beyond about
    /// 10¹⁵) the terms lose their digits, so that the result can fall
    /// outside [0 . . 1], or be NaN, as for *df* = 10²² and *λ* = 1 at
    /// *x* = 9.999999999996 · 10²¹. For a small *λ* the sum can also exceed
    /// 1 by a few ulps in the right tail, so that [`ccdf`] is slightly
    /// negative: with *df* = 3 and *λ* = 10⁻⁹, the cdf at 100 is
    /// 1 + 4.4 · 10⁻¹⁶.
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` fails inside `cumchn`, which needs *df*
    /// beyond about 1.3 · 10²⁹ and *x* within a few ulps of it. It may also
    /// panic for *λ*/2 ≥ 2³¹ − 1001, where CDFLIB's default integers can
    /// overflow.
    ///
    /// [`ccdf`]: ContinuousCdf::ccdf
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: exact endpoint at +inf, where cumchn gives NaN.
        if x == f64::INFINITY {
            return 1.0;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:3866-3875); cumchn
        // returns (0, 1) there.
        // cdflib.f90:3907
        cumchn(x, self.df, self.ncp).0
    }

    /// CDFLIB's `cdfchn` with `which = 1`.
    ///
    /// Unlike the central distributions, CDFLIB computes this as
    /// 1 − [`cdf`] (except for *λ* ≤ 10⁻¹⁰, where it uses the central χ²),
    /// from a series that stops when a term is less than 10⁻⁵ times the
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
        // Rust only: exact endpoint at +inf, where cumchn gives NaN.
        if x == f64::INFINITY {
            return 0.0;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:3866-3875); cumchn
        // returns (0, 1) there.
        // cdflib.f90:3907
        cumchn(x, self.df, self.ncp).1
    }

    /// CDFLIB's `cdfchn` with `which = 2`, searched for in [0 . . 10³⁰⁰].
    ///
    /// The precision limits of [`cdf`] apply: with *df* = 1 and *λ* ≥ 10⁷,
    /// for example, the search for the median fails with
    /// [`AnswerAboveUpperBound`].
    ///
    /// # Panics
    ///
    /// Panics as [`cdf`] does. Since the search evaluates its upper bound
    /// 10³⁰⁰, this happens in particular when *df* is within a few ulps of
    /// 10³⁰⁰.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    /// [`AnswerAboveUpperBound`]: crate::error::SearchError::AnswerAboveUpperBound
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, ChiSquaredNoncentralError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let df = self.df;
        let pnonc = self.ncp;

        // cdflib.f90:3914-3946
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let x = d.x();
            let (cum, _ccum) = cumchn(x, df, pnonc);
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

    #[test]
    fn nan_x_gives_nan() {
        // As in F90, a NaN x makes every term NaN, so the sums stop only at
        // i = 0 or after NTIRED terms, and the result is NaN.
        let d = ChiSquaredNoncentral::new(3.0, 20.0);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
    }

    #[test]
    #[should_panic(expected = "integer overflow")]
    fn huge_ncp_overflows_the_f90_integers() {
        // pnonc / 2 >= 2^31 - 1: int(xnonc) overflows in F90.
        ChiSquaredNoncentral::new(1.0, 5.0e9).cdf(1.0);
    }
}
