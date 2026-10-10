use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::GammaIncError;
use crate::special::{gamma_log, psi};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use thiserror::Error;

use super::gamma::cumgam;

// Parameters of cdfchi (cdflib.f90:3442-3455).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// χ² distribution with *df* degrees of freedom.
///
/// χ²(*df*) is Γ(*df*/2, 2) in shape-scale parameterization. The
/// CDF reduces to the regularized incomplete Γ function:
/// *F*(*x*; *df*) = *P*(*df*/2, *x*/2).
///
/// The methods correspond to CDFLIB's `cdfchi` (cdflib.f90:3340):
/// `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is [`inverse_cdf`] /
/// [`inverse_ccdf`], `which = 3` is [`search_df`].
///
/// # Example
///
/// ```
/// use cdflib::ChiSquared;
/// use cdflib::traits::ContinuousCdf;
///
/// let c = ChiSquared::new(5.0);
///
/// // Pr[X ≤ 11.07] ≈ 0.95
/// let p = c.cdf(11.07);
/// assert!((p - 0.9499903813775945).abs() < 1e-12);
///
/// // Compute df given Pr[X ≤ 3.84] = 0.95
/// let df = ChiSquared::search_df(0.95, 0.05, 3.84).unwrap();
/// assert!((df - 0.9994101496).abs() < 1e-6);
/// ```
///
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`inverse_ccdf`]: ChiSquared::inverse_ccdf
/// [`search_df`]: ChiSquared::search_df
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChiSquared {
    df: f64,
}

/// Errors arising from constructing a [`ChiSquared`] or from its parameter search.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfchi`.
///
/// [`ChiSquared`]: crate::ChiSquared
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum ChiSquaredError {
    /// The degrees of freedom *df* was not strictly positive (`cdfchi`
    /// status −5). [`search_df`] also returns it, checked only in Rust,
    /// when the *df* it computes is 0.
    ///
    /// [`search_df`]: crate::ChiSquared::search_df
    #[error("degrees of freedom must be positive, got {0:?}")]
    DfNotPositive(f64),
    /// The degrees of freedom *df* was not finite (checked only in Rust).
    #[error("degrees of freedom must be finite, got {0:?}")]
    DfNotFinite(f64),
    /// The argument *x* was not strictly positive (`cdfchi` status −4 for
    /// *x* < 0; Rust also rejects *x* = 0).
    #[error("argument x must be positive, got {0:?}")]
    XNotPositive(f64),
    /// The argument *x* was not finite (checked only in Rust).
    #[error("argument x must be finite, got {0:?}")]
    XNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdfchi` status −2); NaN is
    /// also rejected.
    #[error("probability p {0:?} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfchi` status −3); NaN is
    /// also rejected.
    #[error("probability q {0:?} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfchi` status 3).
    #[error("p ({p:?}) and q ({q:?}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfchi` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
    /// The incomplete Γ function failed during the search (`cdfchi` status
    /// 10, cdflib.f90:3602-3605, :3652-3655); see
    /// [`GammaIncError`].
    ///
    /// [`GammaIncError`]: crate::special::GammaIncError
    #[error(transparent)]
    GammaInc(#[from] GammaIncError),
}

/// Returns the cumulative χ² distribution (*cum*, *ccum*) at *x* with *df*
/// degrees of freedom (`cumchi`, cdflib.f90:6860).
///
/// The error value of `gamma_inc`, which CDFLIB passes on unchecked, is
/// returned as a [`GammaIncError`].
///
/// [`GammaIncError`]: crate::special::GammaIncError
#[inline]
pub(crate) fn cumchi(x: f64, df: f64) -> Result<(f64, f64), GammaIncError> {
    let a = df * 0.5;
    let xx = x * 0.5;
    cumgam(xx, a)
}

// cdflib.f90:3484-3502 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), ChiSquaredError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(ChiSquaredError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:3504-3522 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), ChiSquaredError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(ChiSquaredError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:3524-3535 (status -4). Rust also rejects x = 0 and a
// non-finite x.
#[inline]
fn check_x(x: f64) -> Result<(), ChiSquaredError> {
    if x <= 0.0 {
        return Err(ChiSquaredError::XNotPositive(x));
    }
    if !x.is_finite() {
        return Err(ChiSquaredError::XNotFinite(x));
    }
    Ok(())
}

// cdflib.f90:3537-3548 (status -5). Rust also rejects a non-finite df.
#[inline]
fn check_df(df: f64) -> Result<(), ChiSquaredError> {
    if df <= 0.0 {
        return Err(ChiSquaredError::DfNotPositive(df));
    }
    if !df.is_finite() {
        return Err(ChiSquaredError::DfNotFinite(df));
    }
    Ok(())
}

// cdflib.f90:3550-3560 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), ChiSquaredError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(ChiSquaredError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl ChiSquared {
    /// Construct a χ²(*df*) distribution with *df* > 0 degrees of freedom.
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
    /// [`ChiSquaredError`] instead of panicking.
    ///
    /// Returns [`DfNotPositive`] or [`DfNotFinite`] otherwise.
    ///
    /// [`DfNotFinite`]: ChiSquaredError::DfNotFinite
    /// [`DfNotPositive`]: ChiSquaredError::DfNotPositive
    /// [`new`]: Self::new
    /// [`ChiSquaredError`]: crate::ChiSquaredError
    #[inline]
    pub fn try_new(df: f64) -> Result<Self, ChiSquaredError> {
        check_df(df)?;
        Ok(Self { df })
    }

    /// Returns the degrees of freedom *df*.
    #[inline]
    pub const fn df(&self) -> f64 {
        self.df
    }

    /// Returns the degrees of freedom *df* satisfying Pr\[*X* ≤ *x*\] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfchi` with `which = 3`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε. A computed *df* of 0 is
    /// reported as [`DfNotPositive`]. At *p* = 0, where the answer is +∞, the
    /// search returns the finite point where the computed probability
    /// becomes 0 (about 395 for `search_df(0.0, 1.0, 1.0)`).
    ///
    /// [`DfNotPositive`]: ChiSquaredError::DfNotPositive
    #[inline]
    pub fn search_df(p: f64, q: f64, x: f64) -> Result<f64, ChiSquaredError> {
        check_p(p)?;
        check_q(q)?;
        check_x(x)?;
        check_pq(p, q)?;

        // cdflib.f90:3634-3677
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let df = d.x();
            // F90 tests 1.5 < fx + porq, that is, 1.5 < cum or 1.5 < ccum,
            // for the error value of gamma_inc (status 10,
            // cdflib.f90:3652-3655), but gamma_inc leaves qans unset on
            // error, so the test on ccum cannot see it. cumchi returns the
            // error itself, whichever of p and q is smaller.
            let (cum, ccum) = cumchi(x, df)?;
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
        // lower end of the search interval it is not a valid df.
        let df = d.x();
        check_df(df)?;
        Ok(df)
    }

    /// CDFLIB's `cdfchi` with `which = 2`: returns *x* given (*p*, *q*),
    /// already checked.
    fn search_x(&self, p: f64, q: f64) -> Result<f64, ChiSquaredError> {
        let df = self.df;
        // cdflib.f90:3584-3628
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let x = d.x();
            // F90 tests 1.5 < fx + porq, that is, 1.5 < cum or 1.5 < ccum,
            // for the error value of gamma_inc (status 10,
            // cdflib.f90:3602-3605), but gamma_inc leaves qans unset on
            // error, so the test on ccum cannot see it. cumchi returns the
            // error itself, whichever of p and q is smaller.
            let (cum, ccum) = cumchi(x, df)?;
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

    /// Returns the quantile *x* such that [ccdf]\(*x*\) = *q*.
    ///
    /// CDFLIB's `cdfchi` with `which = 2`, with *p* = 1 − *q*.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, ChiSquaredError> {
        check_q(q)?;
        // Rust only: exact endpoints.
        if q == 1.0 {
            return Ok(0.0);
        }
        if q == 0.0 {
            return Ok(f64::INFINITY);
        }
        let p = 1.0 - q;
        self.search_x(p, q)
    }
}

impl ContinuousCdf for ChiSquared {
    type Error = ChiSquaredError;

    /// CDFLIB's `cdfchi` with `which = 1`.
    ///
    /// The result is 0 where the product of *df*/2 and *x*/2 underflows to
    /// 0, even if the true value is far from 0: with *df* = 10⁻⁵,
    /// cdf(10⁻³²⁰) is 0 instead of 0.996.
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` cannot compute its result, which needs
    /// *df* > 1.3 · 10²⁹ and *x* within a few ulps of *df*.
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: exact endpoint at +inf, where cumchi gives NaN.
        if x == f64::INFINITY {
            return 1.0;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:3524-3535); cumchi
        // returns (0, 1) there.
        // cdflib.f90:3570-3578. F90 sets status 10 for the error value of
        // gamma_inc by testing porq, which is not set when which = 1.
        // Rust only: panic on the GammaIncError of cumchi.
        match cumchi(x, self.df) {
            Ok((cum, _ccum)) => cum,
            Err(e) => panic!("cumchi({x:?}, {:?}): {e}", self.df),
        }
    }

    /// CDFLIB's `cdfchi` with `which = 1`.
    ///
    /// The result is 1 where the product of *df*/2 and *x*/2 underflows to
    /// 0, as described for [`cdf`].
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` cannot compute its result, which needs
    /// *df* > 1.3 · 10²⁹ and *x* within a few ulps of *df*.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: exact endpoint at +inf, where cumchi gives NaN.
        if x == f64::INFINITY {
            return 0.0;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:3524-3535); cumchi
        // returns (0, 1) there.
        // cdflib.f90:3570-3578. F90 sets status 10 for the error value of
        // gamma_inc by testing porq, which is not set when which = 1.
        // Rust only: panic on the GammaIncError of cumchi.
        match cumchi(x, self.df) {
            Ok((_cum, ccum)) => ccum,
            Err(e) => panic!("cumchi({x:?}, {:?}): {e}", self.df),
        }
    }

    /// CDFLIB's `cdfchi` with `which = 2`, with *q* = 1 − *p*.
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, ChiSquaredError> {
        check_p(p)?;
        // Rust only: exact endpoints.
        if p == 0.0 {
            return Ok(0.0);
        }
        if p == 1.0 {
            return Ok(f64::INFINITY);
        }
        let q = 1.0 - p;
        self.search_x(p, q)
    }
}

impl Continuous for ChiSquared {
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
        // below would be inf - inf, or 0 * inf for df = 2.
        if x == f64::INFINITY {
            return f64::NEG_INFINITY;
        }
        let k = self.df / 2.0;
        // ln f(x) = -(k ln 2 + ln Γ(k)) + (k - 1) ln x - x/2. Where k is
        // subnormal or 0, so that df/2 loses digits or underflows, ln Γ(k)
        // is computed as -ln k = ln 2 - ln df, with an absolute error below
        // 1e-307.
        let ln_gamma_k = if k < f64::MIN_POSITIVE {
            2.0_f64.ln() - self.df.ln()
        } else {
            gamma_log(k)
        };
        // For k = 1 the term (k - 1) ln x is 0 at x = 0, where it would be
        // 0 · (-inf).
        let ln_x_term = if k == 1.0 && x == 0.0 {
            0.0
        } else {
            (k - 1.0) * x.ln()
        };
        -(k * 2.0_f64.ln() + ln_gamma_k) + ln_x_term - x / 2.0
    }
}

impl Mean for ChiSquared {
    #[inline]
    fn mean(&self) -> f64 {
        self.df
    }
}

impl Variance for ChiSquared {
    #[inline]
    fn variance(&self) -> f64 {
        2.0 * self.df
    }
}

impl Entropy for ChiSquared {
    /// *H* = *k* + ln 2 + ln Γ(*k*) + (1 − *k*) *ψ*(*k*) with *k* = *df*/2.
    #[inline]
    fn entropy(&self) -> f64 {
        let k = self.df / 2.0;
        // For the smallest subnormal df, k is 0, where ψ has a pole; the
        // entropy tends to -inf as df tends to 0.
        if k == 0.0 {
            return f64::NEG_INFINITY;
        }
        k + 2.0_f64.ln() + gamma_log(k) + (1.0 - k) * psi(k)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cdf_at_simple_points() {
        let c = ChiSquared::new(2.0);
        // For df=2, χ² ≡ Exp(1/2); Pr[X ≤ x] = 1 - exp(-x/2).
        for &x in &[0.5_f64, 1.0, 3.84, 10.0] {
            let expected = 1.0 - (-x / 2.0).exp();
            assert!((c.cdf(x) - expected).abs() < 1e-13, "x={x}");
        }
    }

    #[test]
    fn cdf_at_3_84_with_df_1() {
        // χ²₁ at 3.841 ≈ 0.95 (classic statistics-textbook value).
        let c = ChiSquared::new(1.0);
        let p = c.cdf(3.841458820694124);
        assert!((p - 0.95).abs() < 1e-10, "p = {p}");
    }

    #[test]
    fn moments() {
        let c = ChiSquared::new(7.0);
        assert_eq!(c.mean(), 7.0);
        assert_eq!(c.variance(), 14.0);
    }

    #[test]
    fn pdf_nonzero_in_body() {
        let c = ChiSquared::new(4.0);
        for &x in &[1.0, 2.0, 4.0, 8.0] {
            let p = c.pdf(x);
            assert!(p > 0.0 && p < 1.0, "x={x}: pdf={p}");
        }
        // At the mode (df-2 for df>=2): mode of χ²₄ is at 2.
        let m = c.pdf(2.0);
        assert!(m > c.pdf(0.5));
        assert!(m > c.pdf(10.0));
    }

    #[test]
    fn new_rejects_bad_df() {
        assert!(matches!(
            ChiSquared::try_new(f64::NAN),
            Err(ChiSquaredError::DfNotFinite(_))
        ));
        assert!(matches!(
            ChiSquared::try_new(f64::INFINITY),
            Err(ChiSquaredError::DfNotFinite(_))
        ));
        assert!(matches!(
            ChiSquared::try_new(-1.0),
            Err(ChiSquaredError::DfNotPositive(_))
        ));
        assert!(matches!(
            ChiSquared::try_new(0.0),
            Err(ChiSquaredError::DfNotPositive(_))
        ));
    }

    #[test]
    fn search_df_rejects_bad_inputs() {
        assert!(matches!(
            ChiSquared::search_df(-0.1, 1.1, 3.0),
            Err(ChiSquaredError::PNotInRange(_))
        ));
        assert!(matches!(
            ChiSquared::search_df(1.5, -0.5, 3.0),
            Err(ChiSquaredError::PNotInRange(_))
        ));
        assert!(matches!(
            ChiSquared::search_df(0.3, 0.3, 3.0),
            Err(ChiSquaredError::PQSumNotOne { .. })
        ));
        assert!(matches!(
            ChiSquared::search_df(0.5, 0.5, 0.0),
            Err(ChiSquaredError::XNotPositive(0.0))
        ));
        assert!(matches!(
            ChiSquared::search_df(0.5, 0.5, -1.0),
            Err(ChiSquaredError::XNotPositive(-1.0))
        ));
    }

    #[test]
    fn search_df_precision_pivot_at_upper_tail() {
        // For x near the upper tail (p close to 1), the cum-p residual is
        // dominated by 1-cum-eps; the ccum-q form is numerically better.
        // Verify round-trip works in both halves.
        for (p_target, x) in [(0.99, 6.63), (0.999, 10.83), (0.95, 3.84), (0.5, 0.455)] {
            let df = ChiSquared::search_df(p_target, 1.0 - p_target, x).unwrap();
            let cdf_back = ChiSquared::new(df).cdf(x);
            assert!(
                (cdf_back - p_target).abs() < 1e-6,
                "p={p_target}, x={x}, df={df}, cdf_back={cdf_back}"
            );
        }
    }

    #[test]
    fn cdf_at_x_zero_is_zero() {
        let c = ChiSquared::new(5.0);
        assert_eq!(c.cdf(0.0), 0.0);
        assert_eq!(c.cdf(-1.0), 0.0);
    }

    #[test]
    fn ccdf_at_x_zero_is_one() {
        let c = ChiSquared::new(5.0);
        assert_eq!(c.ccdf(0.0), 1.0);
        assert_eq!(c.ccdf(-1.0), 1.0);
    }

    #[test]
    fn inverse_cdf_p_zero_returns_zero() {
        let c = ChiSquared::new(5.0);
        assert_eq!(c.inverse_cdf(0.0).unwrap(), 0.0);
    }

    #[test]
    fn inverse_cdf_rejects_bad_p() {
        let c = ChiSquared::new(5.0);
        assert!(matches!(
            c.inverse_cdf(-0.1),
            Err(ChiSquaredError::PNotInRange(_))
        ));
        assert!(matches!(
            c.inverse_cdf(1.5),
            Err(ChiSquaredError::PNotInRange(_))
        ));
    }

    #[test]
    fn inverse_ccdf_q_one_returns_zero() {
        let c = ChiSquared::new(5.0);
        assert_eq!(c.inverse_ccdf(1.0).unwrap(), 0.0);
    }

    #[test]
    fn inverse_ccdf_rejects_bad_q() {
        let c = ChiSquared::new(5.0);
        assert!(matches!(
            c.inverse_ccdf(-0.1),
            Err(ChiSquaredError::QNotInRange(_))
        ));
        assert!(matches!(
            c.inverse_ccdf(1.5),
            Err(ChiSquaredError::QNotInRange(_))
        ));
    }

    #[test]
    fn pdf_at_x_zero_for_df_le_2_handled() {
        let c = ChiSquared::new(3.0);
        assert_eq!(c.pdf(0.0), 0.0);
        assert_eq!(c.pdf(-1.0), 0.0);
        assert_eq!(c.ln_pdf(0.0), f64::NEG_INFINITY);
        assert_eq!(c.ln_pdf(-1.0), f64::NEG_INFINITY);
    }

    #[test]
    fn entropy_finite_for_df_ge_1() {
        for df in [1.0_f64, 2.0, 5.0, 10.0, 30.0] {
            let h = ChiSquared::new(df).entropy();
            assert!(h.is_finite(), "df={df}: entropy={h}");
        }
    }

    #[test]
    fn density_at_zero_is_the_limit() {
        assert!((ChiSquared::new(2.0).pdf(0.0) - 0.5).abs() < 1e-12);
        assert_eq!(ChiSquared::new(1.0).pdf(0.0), f64::INFINITY);
        assert_eq!(ChiSquared::new(3.0).pdf(0.0), 0.0);
    }

    #[test]
    fn nan_argument_gives_nan() {
        let d = ChiSquared::new(3.0);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
    }

    #[test]
    fn density_at_zero_for_a_subnormal_df() {
        let d = ChiSquared::new(f64::from_bits(1));
        assert_eq!(d.pdf(0.0), f64::INFINITY);
        assert_eq!(d.ln_pdf(0.0), f64::INFINITY);
    }

    // ln df and ln x, about -744, cancel in ln_pdf, so the result depends on
    // the last bits of ln, which Miri's float shims do not guarantee; skip
    // under miri.
    #[cfg(not(miri))]
    #[test]
    fn density_for_a_subnormal_df() {
        // df/2 underflows to 0 for the smallest subnormal df, and rounds to
        // 4/3 of the true value for df = 1.5e-323.
        let d = ChiSquared::new(f64::from_bits(1));
        assert!((d.ln_pdf(0.01) / -740.5330489159531 - 1.0).abs() < 1e-13);
        assert!((d.pdf(5e-324) - 0.5).abs() < 1e-13);
        let d = ChiSquared::new(1.5e-323);
        assert!((d.ln_pdf(1.0) / -744.5346068132731 - 1.0).abs() < 1e-13);
    }
}
