//! Distribution implementations.
//!
//! Each distribution lives in its own private submodule; the value type
//! and its error enum are re-exported here (and from the crate root).
//! Trait impls ([`ContinuousCdf`], [`Continuous`], etc. from
//! [`crate::traits`]) live alongside the struct definition.
//!
//! [`ContinuousCdf`]: crate::traits::ContinuousCdf
//! [`Continuous`]: crate::traits::Continuous
//! [`crate::traits`]: crate::traits

pub(crate) mod beta;
pub(crate) mod binomial;
pub(crate) mod chi_squared;
pub(crate) mod chi_squared_noncentral;
pub(crate) mod fisher_snedecor;
pub(crate) mod fisher_snedecor_noncentral;
pub(crate) mod gamma;
pub(crate) mod negative_binomial;
pub(crate) mod normal;
pub(crate) mod poisson;
pub(crate) mod students_t;

pub use beta::{Beta, BetaError};
pub use binomial::{Binomial, BinomialError};
pub use chi_squared::{ChiSquared, ChiSquaredError};
pub use chi_squared_noncentral::{ChiSquaredNoncentral, ChiSquaredNoncentralError};
pub use fisher_snedecor::{FisherSnedecor, FisherSnedecorError};
pub use fisher_snedecor_noncentral::{FisherSnedecorNoncentral, FisherSnedecorNoncentralError};
pub use gamma::{Gamma, GammaError};
pub use negative_binomial::{NegativeBinomial, NegativeBinomialError};
pub use normal::{Normal, NormalError};
pub use poisson::{Poisson, PoissonError};
pub use students_t::{StudentsT, StudentsTError};

/// Rust only: returns the smallest integer *s* ≤ `max` with `cdf(s)` ≥ *p*,
/// or `max` if there is none, for 0 < *p* < 1 and a nondecreasing `cdf`.
///
/// Used by [`DiscreteCdf::inverse_cdf`]. The bracket is found by doubling
/// and narrowed by bisection.
///
/// [`DiscreteCdf::inverse_cdf`]: crate::traits::DiscreteCdf::inverse_cdf
pub(crate) fn integer_quantile(p: f64, max: u64, cdf: impl Fn(u64) -> f64) -> u64 {
    let mut lo = 0u64;
    let mut hi = 1u64.min(max);
    while cdf(hi) < p {
        if hi == max {
            return max;
        }
        lo = hi + 1;
        hi = hi.saturating_mul(2).min(max);
    }
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        if cdf(mid) < p {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
