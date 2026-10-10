//! Errors shared across distributions.
//!
//! Each distribution module declares its own narrow error enum (so `match`
//! arms stay meaningful). The enums for distributions whose inverse routines
//! go through the reverse-communication root-finder carry a [`SearchError`]
//! variant; distributions that are closed-form everywhere (e.g. [`Normal`])
//! do not. A few distributions whose routines bubble up their own structured
//! errors carry additional pass-through variants (see [`GammaError`] for an
//! example).
//!
//! [`SearchError`]: crate::error::SearchError
//! [`Normal`]: crate::Normal
//! [`GammaError`]: crate::GammaError

use thiserror::Error;

/// Errors of the internal root-finder used by parameter searches and
/// non-closed-form inverse CDFs.
///
/// The two out-of-bounds variants mirror CDFLIB's `status = 1` and `status = 2`
/// (as in `cdfbet`, cdflib.f90:2545-2546): the answer fell below the lowest
/// search bound or above the highest, respectively. `bound` carries the
/// violated endpoint (CDFLIB's `bound` output).
///
/// As in CDFLIB, a search stops once it has located the answer within the
/// larger of an absolute tolerance of 10⁻¹⁰ (10⁻⁵⁰ for the noncentral χ²)
/// and a relative tolerance of 10⁻⁸, so answers below about 10⁻¹⁰ are not
/// resolved: [`ChiSquared`]`::new(1.0).inverse_cdf(1e-6)` returns 0, while
/// the true quantile is 1.6 · 10⁻¹². The inverses of the normal and Γ
/// distributions do not use these searches.
///
/// Also as in CDFLIB, a search can end at a meaningless point without an
/// error when the function it searches evaluates to NaN: CDFLIB's `dinvr`
/// reports a failure of its final phase as a success (cdflib.f90:8469-8474).
/// This happens, for example, for Β and *F* distributions with a parameter
/// above about 2 · 10³⁰⁷, whose CDF is then NaN.
///
/// [`ChiSquared`]: crate::ChiSquared
#[derive(Debug, Clone, Copy, PartialEq, Error)]
pub enum SearchError {
    /// The solution lay below the lower search bound (CDFLIB `status = 1`).
    #[error("answer fell below lower search bound {bound}")]
    AnswerBelowLowerBound { bound: f64 },
    /// The solution lay above the upper search bound (CDFLIB `status = 2`).
    #[error("answer fell above upper search bound {bound}")]
    AnswerAboveUpperBound { bound: f64 },
    /// The initial guess *start* fell outside the range
    /// [*small* . . *big*]. Mirrors CDFLIB's `dinvr` fatal-error abort
    /// at cdflib.f90:8258-8263.
    #[error("start {start} fell outside [{small}..{big}]")]
    StartOutOfRange { start: f64, small: f64, big: f64 },
}
