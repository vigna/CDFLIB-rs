//! Root finders for the distribution-parameter inverters.
//!
//! CDFLIB's `dinvr` / `dzror` reverse-communication state machines plus
//! a closure-driven convenience wrapper ([`search_monotone`]) that drives
//! them with an internal loop.
//!
//! The state machines live in [`dinvr`] and [`dzror`]; see those modules
//! for the algorithm description. Variable names and the iteration trace
//! match the CDFLIB source line-for-line so a debugger comparison is
//! straightforward.
//!
//! # Constants
//!
//! The F90 cdf* dispatchers each declare their own local search-setup
//! constants: `abs_step = rel_step = 0.5`, `stp_mul = 5.0`, `tol = 1e-8`,
//! `atol = 1e-10`. The values are identical across the nine solver-driven
//! routines that declare them; `cdfchn` uses a tighter `atol = 1e-50`,
//! and `cdfnor` needs no solver at all.
//! Rather than restate the duplicated values in eleven Rust callsites,
//! this module centralizes them in the [`ABS_STEP`], [`REL_STEP`],
//! [`STP_MUL`], [`ABS_TOL`], [`REL_TOL`] constants below, and exposes
//! [`search_monotone_with_atol`] for the one case that needs the tighter
//! tolerance. Any future audit against the F90 needs to verify only
//! these two locations.
//!
//! [`ABS_STEP`]: self::ABS_STEP
//! [`REL_STEP`]: self::REL_STEP
//! [`STP_MUL`]: self::STP_MUL
//! [`ABS_TOL`]: self::ABS_TOL
//! [`REL_TOL`]: self::REL_TOL
//! [`search_monotone_with_atol`]: self::search_monotone_with_atol

mod dinvr;
mod dzror;

use crate::error::SearchError;
use dinvr::{InvrAction, InvrConfig, InvrState};
use dzror::{ZrorAction, ZrorConfig, ZrorState};

pub(crate) const SEARCH_BOUND: f64 = 1.0e300;

// CDFLIB's cdf* dispatchers all set up dstinv with the same K-block
// constants: abs_step = rel_step = 0.5 (its K3/K4/K8), stp_mul = 5.0
// (K4/K5/K9), tol = 1e-8. Match them so that dinvr's iteration trace
// is bit-identical to CDFLIB at the dispatcher level. Callers that want
// a tighter converged value can drive InvrState directly with their
// own config.
//
// The default abs_tol = 1e-10 matches every F90 cdf* dispatcher
// EXCEPT cdfchn, which uses 1.0D-50 (cdflib.f90:3719). Callers in
// that regime use [search_monotone_with_atol] to override the default.
const ABS_STEP: f64 = 0.5;
const REL_STEP: f64 = 0.5;
const STP_MUL: f64 = 5.0;
pub(crate) const ABS_TOL: f64 = 1.0e-10;
const REL_TOL: f64 = 1.0e-8;

/// Returns *x* such that *f*(*x*) = 0 on a monotone function, driving
/// CDFLIB's `dinvr` state machine internally. `small` and `big` map onto
/// F90 `dstinv`'s first two arguments (`zsmall`, `zbig`, cdflib.f90:8256);
/// `start` is the initial `x` passed to `dinvr` itself. Monotonicity is
/// inferred by `dinvr` from the endpoint evaluations.
///
/// `qleft_bound` and `qhi_bound` are the literal values each F90 `cdf*`
/// dispatcher writes into the `bound` output for `status = 1` and
/// `status = 2`, respectively (the values reported in
/// [`SearchError::AnswerBelowLowerBound`] / [`AnswerAboveUpperBound`]).
/// Most dispatchers pass `small` and `big`; a few (`cdft which=3`,
/// `cdffnc which=3`, `cdffnc which=4`) write `0.0` for qleft even when
/// `small > 0`, so the bound is explicit here per F90 site.
///
/// Uses the default absolute tolerance [`ABS_TOL`] = 1e-10, matching every
/// CDFLIB `cdf*` dispatcher except `cdfchn`; that one wants a tighter
/// `atol` and uses [`search_monotone_with_atol`] instead.
///
/// [`AnswerAboveUpperBound`]: SearchError::AnswerAboveUpperBound
#[inline]
pub(crate) fn search_monotone(
    small: f64,
    big: f64,
    start: f64,
    qleft_bound: f64,
    qhi_bound: f64,
    f: impl FnMut(f64) -> f64,
) -> Result<f64, SearchError> {
    search_monotone_with_atol(small, big, start, qleft_bound, qhi_bound, ABS_TOL, f)
}

/// Returns *x* such that *f*(*x*) = 0, like [`search_monotone`] but with a
/// caller-supplied `abs_tol`. Used by the noncentral-χ² dispatchers to
/// match `cdfchn`'s `atol = 1e-50`.
#[inline]
pub(crate) fn search_monotone_with_atol(
    small: f64,
    big: f64,
    start: f64,
    qleft_bound: f64,
    qhi_bound: f64,
    abs_tol: f64,
    mut f: impl FnMut(f64) -> f64,
) -> Result<f64, SearchError> {
    // Cap the upper bound at 1e300 the way CDFLIB's cdf* callers do:
    // f64::MAX causes many cdflib::special::* evaluators (e.g.
    // gamma_inc(a, MAX)) to NaN due to Inf-Inf cancellation in their
    // tail formulas. 1e300 is several orders of magnitude beyond any
    // realistic distribution argument.
    let big = big.min(SEARCH_BOUND);
    let small = small.max(-SEARCH_BOUND);

    // CDFLIB's dinvr aborts (ftnstop) if start ∉ [small..big]
    // (cdflib.f90:8020-8024). Return a typed error instead.
    if !(small <= start && start <= big) {
        return Err(SearchError::StartOutOfRange { start, small, big });
    }

    let cfg = InvrConfig {
        small,
        big,
        abs_step: ABS_STEP,
        rel_step: REL_STEP,
        stp_mul: STP_MUL,
        abs_tol,
        rel_tol: REL_TOL,
    };

    let mut state = InvrState::new(cfg, start);
    let mut fx = 0.0;
    // F90's dinvr has no eval cap; iteration runs until the state machine
    // reports Converged or Failed (cdflib.f90:E0000 reverse-communication).
    loop {
        match state.step(fx) {
            InvrAction::NeedEval(x) => {
                fx = f(x);
            }
            InvrAction::Converged(x) => return Ok(x),
            InvrAction::Failed { qleft, .. } => {
                return Err(if qleft {
                    SearchError::AnswerBelowLowerBound { bound: qleft_bound }
                } else {
                    SearchError::AnswerAboveUpperBound { bound: qhi_bound }
                });
            }
        }
    }
}

/// Returns *x* such that *f*(*x*) = 0 on a a range [xlo, xhi], driving
/// CDFLIB's `dzror` state machine directly.
#[inline]
pub(crate) fn search_bounded_zero(
    xlo: f64,
    xhi: f64,
    mut f: impl FnMut(f64) -> f64,
) -> Result<f64, SearchError> {
    search_bounded_zero_with_tol(xlo, xhi, ABS_TOL, REL_TOL, &mut f)
}

#[inline]
pub(crate) fn search_bounded_zero_with_tol(
    xlo: f64,
    xhi: f64,
    abs_tol: f64,
    rel_tol: f64,
    f: &mut impl FnMut(f64) -> f64,
) -> Result<f64, SearchError> {
    let cfg = ZrorConfig {
        xlo,
        xhi,
        abstol: abs_tol,
        reltol: rel_tol,
    };
    let mut state = ZrorState::new(cfg);
    let mut fx = 0.0;
    loop {
        match state.step(fx) {
            ZrorAction::NeedEval(x) => fx = f(x),
            // The F90 direct-dzror dispatchers read the answer from the
            // reverse-communication variable x, which still holds the
            // last evaluation point; the swap at label 80 updates xlo
            // but not x.
            ZrorAction::Converged { x, .. } => return Ok(x),
            ZrorAction::Failed { qleft, .. } => {
                return Err(if qleft {
                    SearchError::AnswerBelowLowerBound { bound: cfg.xlo }
                } else {
                    SearchError::AnswerAboveUpperBound { bound: cfg.xhi }
                })
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_increasing_function() {
        // f(x) = x³ - 8; root at x = 2.
        let r = search_monotone(0.0, 100.0, 1.0, 0.0, 100.0, |x| x.powi(3) - 8.0).unwrap();
        assert!((r - 2.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn solves_decreasing_function() {
        // f(x) = 1/x - 0.25; root at x = 4.
        let r = search_monotone(0.01, 1000.0, 10.0, 0.01, 1000.0, |x| 1.0 / x - 0.25).unwrap();
        assert!((r - 4.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn solves_root_at_moderate_value() {
        // f(x) = ln(x) - 1 → root at x = e.
        let r = search_monotone(1e-10, 1000.0, 1.0, 1e-10, 1000.0, |x| x.ln() - 1.0).unwrap();
        let e = std::f64::consts::E;
        assert!((r - e).abs() / e < 1e-8, "r = {r}, e = {e}");
    }

    // ============================ Failure paths in dinvr ============================
    //
    // These cover the four range-validity branches and the qlim overshoot
    // failures. Each one constructs a function where the [small..big] range
    // does NOT enclose a root, or the root lies outside even after expansion.

    #[test]
    fn increasing_fsmall_positive_fails_at_small() {
        // f is monotone increasing but already positive at small.
        let r = search_monotone(1.0, 10.0, 5.0, 1.0, 10.0, |x| x + 1.0);
        assert!(matches!(
            r,
            Err(SearchError::AnswerBelowLowerBound { bound }) if bound == 1.0
        ));
    }

    #[test]
    fn increasing_fbig_negative_fails_at_big() {
        // Monotone increasing but f(big) still negative.
        let r = search_monotone(1.0, 10.0, 5.0, 1.0, 10.0, |x| x - 100.0);
        assert!(matches!(
            r,
            Err(SearchError::AnswerAboveUpperBound { bound }) if bound == 10.0
        ));
    }

    #[test]
    fn decreasing_fsmall_negative_fails_at_small() {
        // Monotone decreasing but f(small) already negative.
        let r = search_monotone(1.0, 10.0, 5.0, 1.0, 10.0, |x| -x - 1.0);
        assert!(matches!(
            r,
            Err(SearchError::AnswerBelowLowerBound { bound }) if bound == 1.0
        ));
    }

    #[test]
    fn decreasing_fbig_positive_fails_at_big() {
        // Monotone decreasing but f(big) still positive.
        let r = search_monotone(1.0, 10.0, 5.0, 1.0, 10.0, |x| 100.0 - x);
        assert!(matches!(
            r,
            Err(SearchError::AnswerAboveUpperBound { bound }) if bound == 10.0
        ));
    }

    #[test]
    fn converges_immediately_when_fstart_zero() {
        // start happens to be the root → AwaitInitial returns Converged.
        let r = search_monotone(0.0, 10.0, 3.0, 0.0, 10.0, |x| x - 3.0).unwrap();
        assert!((r - 3.0).abs() < 1e-15);
    }

    #[test]
    fn qleft_qhi_bounds_are_passed_through_distinctly() {
        // The qleft_bound / qhi_bound parameters must surface in the
        // failure variant verbatim, independent of small / big. F90
        // cdf* dispatchers exploit this (e.g. cdft which=3 writes
        // bound = 0.0D+00 even though small = 1.0).
        // qleft path: increasing function positive everywhere → small fails.
        let err = search_monotone(1.0, 10.0, 5.0, 99.0, 999.0, |x| x + 1.0).unwrap_err();
        assert!(
            matches!(err, SearchError::AnswerBelowLowerBound { bound } if bound == 99.0),
            "expected AnswerBelowLowerBound {{ bound: 99.0 }}, got {err:?}"
        );
        // qhi path: increasing function negative everywhere → big fails.
        let err = search_monotone(1.0, 10.0, 5.0, 99.0, 999.0, |x| x - 100.0).unwrap_err();
        assert!(
            matches!(err, SearchError::AnswerAboveUpperBound { bound } if bound == 999.0),
            "expected AnswerAboveUpperBound {{ bound: 999.0 }}, got {err:?}"
        );
    }

    #[test]
    fn nan_objective_surfaces_as_search_failure() {
        // A NaN-returning objective surfaces as AnswerBelowLowerBound
        // (the initial NaN at start is neither < 0 nor > 0, so the
        // range-expansion logic falls through).
        let err = search_monotone(0.0, 1.0, 0.5, 0.0, 1.0, |_| f64::NAN).unwrap_err();
        assert!(matches!(
            err,
            SearchError::AnswerBelowLowerBound { bound: 0.0 }
        ));
    }
}
