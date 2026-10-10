use crate::error::SearchError;
use crate::search::dstinv;
use crate::special::{
    gamma_log, psi, try_gamma_inc, try_gamma_inc_inv, GammaIncError, GammaIncInvError,
};
use crate::traits::{Continuous, ContinuousCdf, Entropy, Mean, Variance};
use thiserror::Error;

// Parameters of cdfgam (cdflib.f90:4958-4972).
const ATOL: f64 = 1.0e-10;
const INF: f64 = 1.0e300;
const TOL: f64 = 1.0e-8;

/// Γ distribution with *α* > 0 (shape) and *β* > 0 (rate). Mean = *α*/*β*.
///
/// Density *f*(*x*; *α*, *β*) = (*βᵅ* / Γ(*α*)) · *xᵅ* ⁻ ¹ · exp(−*β*·*x*) for
/// *x* > 0. The CDF reduces to the regularized incomplete Γ function:
/// *F*(*x*; *α*, *β*) = *P*(*α*, *β*·*x*).
///
/// The methods correspond to CDFLIB's `cdfgam` (cdflib.f90:4854):
/// `which = 1` is [`cdf`] / [`ccdf`], `which = 2` is [`inverse_cdf`] /
/// [`inverse_ccdf`], `which = 3` is [`search_shape`], `which = 4` is
/// [`search_rate`].
///
/// # Note on naming
///
/// CDFLIB's `cdfgam` calls its second parameter SCALE, but defines the
/// density as proportional to *t*^(SHAPE − 1) · exp(−SCALE · *t*) and
/// computes `cumgam(x * scale, shape, …)`: the parameter is mathematically
/// the **rate** *β* (mean = *α*/*β*), not the conventional scale *θ*
/// (mean = *α*·*θ*, with CDF *P*(*α*, *x*/*θ*)). Users with shape-scale
/// parameters should pass `rate = 1.0 / scale`.
///
/// # Example
///
/// ```
/// use cdflib::Gamma;
/// use cdflib::traits::ContinuousCdf;
///
/// let g = Gamma::new(2.0, 1.0);
///
/// // Pr[X ≤ 2.0] = 1 - 3 e⁻²
/// let p = g.cdf(2.0);
/// assert!((p - 0.5939941502901619).abs() < 1e-12);
///
/// // Compute shape parameter given Pr[X ≤ 5.0] = 0.9 and rate = 2.0
/// let shape = Gamma::search_shape(0.9, 0.1, 5.0, 2.0).unwrap();
/// assert!((shape - 6.574843866).abs() < 1e-6);
/// ```
///
/// [`cdf`]: ContinuousCdf::cdf
/// [`ccdf`]: ContinuousCdf::ccdf
/// [`inverse_cdf`]: ContinuousCdf::inverse_cdf
/// [`inverse_ccdf`]: Gamma::inverse_ccdf
/// [`search_shape`]: Gamma::search_shape
/// [`search_rate`]: Gamma::search_rate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gamma {
    shape: f64,
    rate: f64,
}

/// Errors arising from constructing a [`Gamma`] or from its parameter searches.
///
/// The variants correspond to the `status` codes of CDFLIB's `cdfgam`.
///
/// [`Gamma`]: crate::Gamma
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum GammaError {
    /// The shape parameter *α* was not strictly positive (`cdfgam` status −5).
    /// [`search_shape`] also returns it, checked only in Rust, when the
    /// shape it computes is 0.
    ///
    /// [`search_shape`]: crate::Gamma::search_shape
    #[error("shape must be positive, got {0}")]
    ShapeNotPositive(f64),
    /// The rate parameter *β* (CDFLIB's SCALE) was not strictly positive
    /// (`cdfgam` status −6). [`search_rate`] also returns it, checked only
    /// in Rust, when the rate it computes is not strictly positive, as for
    /// *p* = 0.
    ///
    /// [`search_rate`]: crate::Gamma::search_rate
    #[error("rate must be positive, got {0}")]
    RateNotPositive(f64),
    /// The shape parameter *α* was not finite (checked only in Rust).
    #[error("shape must be finite, got {0}")]
    ShapeNotFinite(f64),
    /// The rate parameter *β* was not finite (checked only in Rust).
    /// [`search_rate`] also returns it when the rate it computes is not
    /// finite, as for *q* = 0.
    ///
    /// [`search_rate`]: crate::Gamma::search_rate
    #[error("rate must be finite, got {0}")]
    RateNotFinite(f64),
    /// The argument *x* was not strictly positive (`cdfgam` status −4 for
    /// *x* < 0; Rust also rejects *x* = 0, where Pr[*X* ≤ *x*] = 0 whatever
    /// the parameters).
    #[error("argument x must be positive, got {0}")]
    XNotPositive(f64),
    /// The argument *x* was not finite (checked only in Rust).
    #[error("argument x must be finite, got {0}")]
    XNotFinite(f64),
    /// The probability *p* fell outside [0 . . 1] (`cdfgam` status −2); NaN is
    /// also rejected.
    #[error("probability p {0} outside [0..1]")]
    PNotInRange(f64),
    /// The probability *q* fell outside [0 . . 1] (`cdfgam` status −3); NaN is
    /// also rejected.
    #[error("probability q {0} outside [0..1]")]
    QNotInRange(f64),
    /// The pair (*p*, *q*) is not complementary: 3ε < |*p* + *q* − 1|
    /// (`cdfgam` status 3).
    #[error("p ({p}) and q ({q}) are not complementary: |p + q - 1| > 3ε")]
    PQSumNotOne { p: f64, q: f64 },
    /// The search for the answer failed (`cdfgam` status 1 or 2); see
    /// [`SearchError`].
    ///
    /// [`SearchError`]: crate::error::SearchError
    #[error(transparent)]
    Search(#[from] SearchError),
    /// The inverse incomplete Γ function failed (`cdfgam` status 10 after
    /// `gamma_inc_inv` returns a negative `ierr`); see [`GammaIncInvError`].
    /// [`AtInfinity`] does not occur: *q* = 0 is handled as an endpoint
    /// in the quantile methods and reported as [`RateNotFinite`] by
    /// [`search_rate`].
    ///
    /// [`RateNotFinite`]: GammaError::RateNotFinite
    /// [`search_rate`]: crate::Gamma::search_rate
    ///
    /// [`AtInfinity`]: crate::special::GammaIncInvError::AtInfinity
    ///
    /// [`GammaIncInvError`]: crate::special::GammaIncInvError
    #[error(transparent)]
    GammaIncInv(#[from] GammaIncInvError),
    /// The incomplete Γ function failed during the search (`cdfgam` status
    /// 10, cdflib.f90:5150-5156); see [`GammaIncError`].
    ///
    /// [`GammaIncError`]: crate::special::GammaIncError
    #[error(transparent)]
    GammaInc(#[from] GammaIncError),
}

/// Returns the cumulative Γ distribution (*cum*, *ccum*) at *x* with shape
/// parameter *a* and unit rate (`cumgam`, cdflib.f90:7455).
///
/// CDFLIB passes on the error value of `gamma_inc` unchecked; here it is
/// returned as the [`GammaIncError`] of [`try_gamma_inc`].
///
/// [`GammaIncError`]: crate::special::GammaIncError
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
#[inline]
pub(crate) fn cumgam(x: f64, a: f64) -> Result<(f64, f64), GammaIncError> {
    if x <= 0.0 {
        Ok((0.0, 1.0))
    } else {
        try_gamma_inc(a, x)
    }
}

// cdflib.f90:5003-5021 (status -2). Rust also rejects NaN.
#[inline]
fn check_p(p: f64) -> Result<(), GammaError> {
    if p < 0.0 || 1.0 < p || p.is_nan() {
        return Err(GammaError::PNotInRange(p));
    }
    Ok(())
}

// cdflib.f90:5023-5041 (status -3). Rust also rejects NaN.
#[inline]
fn check_q(q: f64) -> Result<(), GammaError> {
    if q < 0.0 || 1.0 < q || q.is_nan() {
        return Err(GammaError::QNotInRange(q));
    }
    Ok(())
}

// cdflib.f90:5043-5054 (status -4). Rust also rejects x = 0 and a
// non-finite x.
#[inline]
fn check_x(x: f64) -> Result<(), GammaError> {
    if x <= 0.0 {
        return Err(GammaError::XNotPositive(x));
    }
    if !x.is_finite() {
        return Err(GammaError::XNotFinite(x));
    }
    Ok(())
}

// cdflib.f90:5056-5067 (status -5). Rust also rejects a non-finite shape.
#[inline]
fn check_shape(shape: f64) -> Result<(), GammaError> {
    if shape <= 0.0 {
        return Err(GammaError::ShapeNotPositive(shape));
    }
    if !shape.is_finite() {
        return Err(GammaError::ShapeNotFinite(shape));
    }
    Ok(())
}

// cdflib.f90:5069-5080 (status -6). Rust also rejects a non-finite rate.
#[inline]
fn check_rate(rate: f64) -> Result<(), GammaError> {
    if rate <= 0.0 {
        return Err(GammaError::RateNotPositive(rate));
    }
    if !rate.is_finite() {
        return Err(GammaError::RateNotFinite(rate));
    }
    Ok(())
}

// cdflib.f90:5082-5092 (status 3).
#[inline]
fn check_pq(p: f64, q: f64) -> Result<(), GammaError> {
    if 3.0 * f64::EPSILON < ((p + q) - 1.0).abs() {
        return Err(GammaError::PQSumNotOne { p, q });
    }
    Ok(())
}

impl Gamma {
    /// Construct a Γ(*α*, *β*) distribution with shape *α* > 0 and rate
    /// *β* > 0.
    ///
    /// # Panics
    ///
    /// Panics if either argument is invalid; use [`try_new`] for a fallible
    /// variant.
    ///
    /// [`try_new`]: Self::try_new
    #[inline]
    pub fn new(shape: f64, rate: f64) -> Self {
        Self::try_new(shape, rate).unwrap()
    }

    /// Fallible counterpart of [`new`](Self::new) returning a [`GammaError`]
    /// instead of panicking.
    ///
    /// Returns [`ShapeNotPositive`], [`ShapeNotFinite`], [`RateNotPositive`],
    /// or [`RateNotFinite`] if either argument fails its respective test.
    ///
    /// [`ShapeNotFinite`]: GammaError::ShapeNotFinite
    /// [`RateNotFinite`]: GammaError::RateNotFinite
    /// [`ShapeNotPositive`]: GammaError::ShapeNotPositive
    /// [`RateNotPositive`]: GammaError::RateNotPositive
    #[inline]
    pub fn try_new(shape: f64, rate: f64) -> Result<Self, GammaError> {
        check_shape(shape)?;
        check_rate(rate)?;
        Ok(Self { shape, rate })
    }

    /// Returns the shape parameter *α*.
    #[inline]
    pub const fn shape(&self) -> f64 {
        self.shape
    }

    /// Returns the rate parameter *β*.
    #[inline]
    pub const fn rate(&self) -> f64 {
        self.rate
    }

    /// Returns the shape parameter *α* satisfying Pr[*X* ≤ *x*] = *p*,
    /// searched for in [0 . . 10³⁰⁰].
    ///
    /// CDFLIB's `cdfgam` with `which = 3`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε.
    ///
    /// A computed shape of 0, at the lower end of the search interval, is
    /// reported as [`ShapeNotPositive`]. At *p* = 0, where the answer is +∞,
    /// the search stops, as the F90 does, where the computed probability
    /// becomes 0, and returns that finite value:
    /// `search_shape(0.0, 1.0, 1.0, 1.0)` returns about 395.
    ///
    /// [`ShapeNotPositive`]: GammaError::ShapeNotPositive
    #[inline]
    pub fn search_shape(p: f64, q: f64, x: f64, rate: f64) -> Result<f64, GammaError> {
        check_p(p)?;
        check_q(q)?;
        check_x(x)?;
        check_rate(rate)?;
        check_pq(p, q)?;
        // F90 SCALE.
        let scale = rate;

        // cdflib.f90:5128-5178
        let xscale = x * scale;
        let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
        while d.status() == 1 {
            let shape = d.x();
            // F90 tests 1.5 < cum or 1.5 < ccum for the error value of
            // gamma_inc (status 10, cdflib.f90:5150-5156), but gamma_inc
            // leaves qans unset on error, so the test on ccum cannot see
            // it. cumgam returns the error itself, whichever of p and q is
            // smaller.
            let (cum, ccum) = cumgam(xscale, shape)?;
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
        // lower end of the search interval it is not a valid shape.
        let shape = d.x();
        check_shape(shape)?;
        Ok(shape)
    }

    /// Returns the rate parameter *β* (CDFLIB's SCALE) satisfying
    /// Pr[*X* ≤ *x*] = *p*.
    ///
    /// CDFLIB's `cdfgam` with `which = 4`. The caller passes both *p* and
    /// *q* = 1 − *p*; they must sum to 1 within 3ε. A computed rate that is
    /// not positive or not finite, as at *p* = 0 or *p* = 1, is reported as
    /// [`RateNotPositive`] or [`RateNotFinite`].
    ///
    /// [`RateNotPositive`]: GammaError::RateNotPositive
    /// [`RateNotFinite`]: GammaError::RateNotFinite
    #[inline]
    pub fn search_rate(p: f64, q: f64, x: f64, shape: f64) -> Result<f64, GammaError> {
        check_p(p)?;
        check_q(q)?;
        check_x(x)?;
        check_shape(shape)?;
        check_pq(p, q)?;

        // cdflib.f90:5182-5191. A negative ierr (status 10) is the
        // GammaIncInvError of try_gamma_inc_inv.
        let xx = match try_gamma_inc_inv(shape, -1.0, p, q) {
            Ok((xx, _ierr)) => xx,
            // Rust only: for q = 0, gamma_inc_inv gives xx = huge(xx) with
            // ierr = 0 (cdflib.f90:11460-11463), and the F90 returns
            // scale = huge(xx) / x; the rate is +inf.
            Err(GammaIncInvError::AtInfinity) => {
                return Err(GammaError::RateNotFinite(f64::INFINITY));
            }
            Err(e) => return Err(e.into()),
        };
        let scale = xx / x;
        // Rust only: the F90 returns scale whatever its value; for p = 0,
        // gamma_inc_inv gives xx = 0.
        if scale <= 0.0 {
            return Err(GammaError::RateNotPositive(scale));
        }
        if !scale.is_finite() {
            return Err(GammaError::RateNotFinite(scale));
        }
        Ok(scale)
    }

    /// CDFLIB's `cdfgam` with `which = 2`: returns *x* given (*p*, *q*),
    /// already checked, with *q* > 0.
    fn search_x(&self, p: f64, q: f64) -> Result<f64, GammaError> {
        // F90 SCALE.
        let scale = self.rate;
        // cdflib.f90:5114-5124. A negative ierr (status 10) is the
        // GammaIncInvError of try_gamma_inc_inv.
        let (xx, _ierr) = try_gamma_inc_inv(self.shape, -1.0, p, q)?;
        let x = xx / scale;
        Ok(x)
    }

    /// Returns the quantile *x* such that [ccdf]\(*x*\) = *q*.
    ///
    /// CDFLIB's `cdfgam` with `which = 2`, with *p* = 1 − *q*, so that a
    /// tiny *q* keeps its precision.
    ///
    /// [ccdf]: crate::traits::ContinuousCdf::ccdf
    #[inline]
    pub fn inverse_ccdf(&self, q: f64) -> Result<f64, GammaError> {
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

impl ContinuousCdf for Gamma {
    type Error = GammaError;

    /// CDFLIB's `cdfgam` with `which = 1`.
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` cannot compute its result, which needs
    /// shape > 6.6 · 10²⁸ and *β*·*x* within a few ulps of the shape.
    #[inline]
    fn cdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:5043-5054); cumgam
        // returns (0, 1) there.
        // cdflib.f90:5102-5110. F90 sets status 10 for the error value of
        // gamma_inc by testing porq, which is not set when which = 1.
        let xscale = x * self.rate;
        // Rust only: exact endpoint where xscale is +inf, for which cumgam
        // gives NaN; this includes x = +inf.
        if xscale == f64::INFINITY {
            return 1.0;
        }
        // Rust only: panic on the GammaIncError of cumgam.
        match cumgam(xscale, self.shape) {
            Ok((cum, _ccum)) => cum,
            Err(e) => panic!("cumgam({xscale}, {}): {e}", self.shape),
        }
    }

    /// CDFLIB's `cdfgam` with `which = 1`.
    ///
    /// # Panics
    ///
    /// Panics where `gamma_inc` cannot compute its result, which needs
    /// shape > 6.6 · 10²⁸ and *β*·*x* within a few ulps of the shape.
    #[inline]
    fn ccdf(&self, x: f64) -> f64 {
        // Rust only: NaN for a NaN x.
        if x.is_nan() {
            return f64::NAN;
        }
        // Rust only: no status -4 for x < 0 (cdflib.f90:5043-5054); cumgam
        // returns (0, 1) there.
        // cdflib.f90:5102-5110. F90 sets status 10 for the error value of
        // gamma_inc by testing porq, which is not set when which = 1.
        let xscale = x * self.rate;
        // Rust only: exact endpoint where xscale is +inf, for which cumgam
        // gives NaN; this includes x = +inf.
        if xscale == f64::INFINITY {
            return 0.0;
        }
        // Rust only: panic on the GammaIncError of cumgam.
        match cumgam(xscale, self.shape) {
            Ok((_cum, ccum)) => ccum,
            Err(e) => panic!("cumgam({xscale}, {}): {e}", self.shape),
        }
    }

    /// CDFLIB's `cdfgam` with `which = 2`, with *q* = 1 − *p*.
    #[inline]
    fn inverse_cdf(&self, p: f64) -> Result<f64, GammaError> {
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

impl Continuous for Gamma {
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
        // below would be inf - inf, or 0 * inf for shape = 1.
        if x == f64::INFINITY {
            return f64::NEG_INFINITY;
        }
        // ln f = shape·ln(rate) - ln Γ(shape) + (shape-1) ln x - rate·x; for
        // shape = 1 the term (shape-1) ln x is 0 at x = 0, where it would be
        // 0 · (-inf).
        let ln_x_term = if self.shape == 1.0 && x == 0.0 {
            0.0
        } else {
            (self.shape - 1.0) * x.ln()
        };
        self.shape * self.rate.ln() - gamma_log(self.shape) + ln_x_term - self.rate * x
    }
}

impl Mean for Gamma {
    #[inline]
    fn mean(&self) -> f64 {
        self.shape / self.rate
    }
}

impl Variance for Gamma {
    #[inline]
    fn variance(&self) -> f64 {
        // shape / rate², divided in two steps so that rate² cannot overflow
        // or underflow.
        self.shape / self.rate / self.rate
    }
}

impl Entropy for Gamma {
    /// *H* = *α* − ln *β* + ln Γ(*α*) + (1 − *α*) *ψ*(*α*).
    #[inline]
    fn entropy(&self) -> f64 {
        self.shape - self.rate.ln() + gamma_log(self.shape) + (1.0 - self.shape) * psi(self.shape)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cdf_reduces_to_exponential_for_shape_1() {
        // Gamma(1, β) ≡ Exp(β): CDF = 1 - exp(-β·x).
        let g = Gamma::new(1.0, 2.0);
        for &x in &[0.5_f64, 1.0, 4.0, 10.0] {
            let expected = 1.0 - (-x * 2.0).exp();
            assert!((g.cdf(x) - expected).abs() < 1e-13, "x={x}");
        }
    }

    #[test]
    fn moments() {
        // Gamma(shape=3, rate=2): mean = 3/2, variance = 3/4.
        let g = Gamma::new(3.0, 2.0);
        assert_eq!(g.mean(), 1.5);
        assert_eq!(g.variance(), 0.75);
    }

    #[test]
    fn pdf_at_mode() {
        // For shape > 1, the mode of Gamma(α, β) is at (α-1)/β.
        let g = Gamma::new(3.0, 2.0);
        let mode = (3.0 - 1.0) / 2.0;
        let pm = g.pdf(mode);
        assert!(pm > g.pdf(mode * 0.5));
        assert!(pm > g.pdf(mode * 2.0));
    }

    #[test]
    fn rejects_invalid_parameters_and_probabilities() {
        assert!(matches!(
            Gamma::try_new(0.0, 1.0),
            Err(GammaError::ShapeNotPositive(0.0))
        ));
        assert!(matches!(
            Gamma::try_new(1.0, 0.0),
            Err(GammaError::RateNotPositive(0.0))
        ));
        assert!(matches!(
            Gamma::try_new(f64::INFINITY, 1.0),
            Err(GammaError::ShapeNotFinite(x)) if x.is_infinite()
        ));
        assert!(matches!(
            Gamma::try_new(1.0, f64::INFINITY),
            Err(GammaError::RateNotFinite(x)) if x.is_infinite()
        ));
        assert!(matches!(
            Gamma::search_shape(-0.1, 1.1, 1.0, 1.0),
            Err(GammaError::PNotInRange(-0.1))
        ));
    }

    #[test]
    fn inverse_and_density_edges() {
        let g = Gamma::new(2.0, 3.0);
        assert_eq!(g.inverse_cdf(0.0).unwrap(), 0.0);
        assert_eq!(g.inverse_ccdf(1.0).unwrap(), 0.0);
        assert_eq!(g.pdf(0.0), 0.0);
        assert_eq!(g.ln_pdf(0.0), f64::NEG_INFINITY);
        assert_eq!(g.cdf(-1.0), 0.0);
        assert_eq!(g.ccdf(-1.0), 1.0);
        assert!(g.ccdf(1.0).is_finite());
        assert!(g.inverse_ccdf(0.25).unwrap().is_finite());
        assert!(g.entropy().is_finite());
    }

    #[test]
    fn search_parameter_rejects_nonpositive_inputs() {
        assert!(matches!(
            Gamma::search_shape(0.5, 0.5, 0.0, 1.0),
            Err(GammaError::XNotPositive(0.0))
        ));
        assert!(matches!(
            Gamma::search_shape(0.5, 0.5, 1.0, 0.0),
            Err(GammaError::RateNotPositive(0.0))
        ));
        assert!(matches!(
            Gamma::search_rate(0.5, 0.5, 1.0, 0.0),
            Err(GammaError::ShapeNotPositive(0.0))
        ));
        assert!(matches!(
            Gamma::search_rate(0.5, 0.5, -0.1, 2.0),
            Err(GammaError::XNotPositive(x)) if x == -0.1
        ));
    }

    #[test]
    fn search_rate_endpoints() {
        // q = 0: gamma_inc_inv gives huge(x) with ierr = 0 in F90.
        assert_eq!(
            Gamma::search_rate(1.0, 0.0, 2.0, 3.0),
            Err(GammaError::RateNotFinite(f64::INFINITY))
        );
        // p = 0: gamma_inc_inv gives x = 0, so the F90 rate is 0.
        assert_eq!(
            Gamma::search_rate(0.0, 1.0, 2.0, 3.0),
            Err(GammaError::RateNotPositive(0.0))
        );
        // q = 0 with p + q - 1 beyond epsilon but within 3 epsilon: cdfgam
        // passes, and gamma_inc_inv returns ierr = -4 before testing q = 0.
        let p = 1.0 - 4.0 * f64::EPSILON / 2.0;
        assert_eq!(
            Gamma::search_rate(p, 0.0, 2.0, 3.0),
            Err(GammaError::GammaIncInv(GammaIncInvError::InconsistentPq))
        );
    }

    #[test]
    fn density_at_zero_is_the_limit() {
        assert!((Gamma::new(1.0, 2.0).pdf(0.0) - 2.0).abs() < 1e-12);
        assert_eq!(Gamma::new(0.5, 1.0).pdf(0.0), f64::INFINITY);
        assert_eq!(Gamma::new(2.0, 1.0).pdf(0.0), 0.0);
        assert_eq!(Gamma::new(1.0, 2.0).pdf(-1.0), 0.0);
    }

    #[test]
    fn nan_argument_gives_nan() {
        let d = Gamma::new(2.0, 1.0);
        assert!(d.cdf(f64::NAN).is_nan());
        assert!(d.ccdf(f64::NAN).is_nan());
    }

    #[test]
    fn variance_does_not_overflow() {
        assert!((Gamma::new(1e300, 1e200).variance() / 1e-100 - 1.0).abs() < 1e-12);
    }
}
