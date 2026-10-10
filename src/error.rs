//! Errors shared across distributions.
//!
//! Each distribution has its own error enum. Those whose inverses or
//! searches use the root finder have a [`SearchError`] variant, and a few
//! also pass through the errors of the special functions they call (see
//! [`GammaError`]).
//!
//! [`SearchError`]: crate::error::SearchError
//! [`GammaError`]: crate::GammaError

use thiserror::Error;

/// Errors of the root finder used by the parameter searches and by the
/// inverses that are not in closed form.
///
/// The out-of-bounds variants are CDFLIB's `status = 1` and `status = 2`
/// (as in `cdfbet`, cdflib.f90:2545-2546); `bound` is the violated end of
/// the search interval.
///
/// As in CDFLIB, a search stops once it has located the answer within the
/// larger of an absolute tolerance of 10⁻¹⁰ (10⁻⁵⁰ for the noncentral χ²)
/// and a relative tolerance of 10⁻⁸, so answers below about 10⁻¹⁰ are not
/// resolved: [`ChiSquared`]`::new(1.0).inverse_cdf(1e-6)` returns 0
/// instead of 1.6 · 10⁻¹². The inverses of the normal and Γ distributions
/// do not use these searches.
///
/// Also as in CDFLIB, a search ends without an error at a meaningless
/// point where the function it searches is NaN (cdflib.f90:8469-8474), as
/// for Β and *F* distributions with a parameter above about 2 · 10³⁰⁷, and
/// at a wrong answer where that function is wrong, as in the left tails of
/// the noncentral χ² and *F* distributions:
/// [`ChiSquaredNoncentral`]`::search_ncp(1e-10, 0.1, 0.5)` returns about
/// 22 instead of about 46.3.
///
/// Finally, a search whose target probability (the smaller of *p* and *q*)
/// is below about 10⁻¹⁵⁸ can end without an error away from the answer, or
/// fail naming the wrong bound, because `dzror` compares signs through a
/// product that underflows (cdflib.f90:9100):
/// [`ChiSquared`]`::new(10.0).inverse_ccdf(1e-200)` returns about 1176.25,
/// where the upper tail is about 2 · 10⁻²⁴⁶.
///
/// [`ChiSquared`]: crate::ChiSquared
/// [`ChiSquaredNoncentral`]: crate::ChiSquaredNoncentral
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum SearchError {
    /// The solution lay below the lower search bound (CDFLIB `status = 1`).
    #[error("answer fell below lower search bound {bound:?}")]
    AnswerBelowLowerBound { bound: f64 },
    /// The solution lay above the upper search bound (CDFLIB `status = 2`).
    #[error("answer fell above upper search bound {bound:?}")]
    AnswerAboveUpperBound { bound: f64 },
    /// The initial guess *start* fell outside the range
    /// [*small* . . *big*]. Mirrors CDFLIB's `dinvr` fatal-error abort
    /// at cdflib.f90:8258-8263.
    #[error("start {start:?} fell outside [{small:?}..{big:?}]")]
    StartOutOfRange { start: f64, small: f64, big: f64 },
}
