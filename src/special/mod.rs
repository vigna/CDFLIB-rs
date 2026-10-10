//! Special functions underlying the distributions.
//!
//! This module is split into two surfaces.
//!
//! The top-level [`cdflib::special`](crate::special) namespace exposes the
//! user-facing special functions a statistical user is likely to call directly:
//! [`beta`], [`beta_log`], [`beta_inc`], [`gamma`], [`gamma_log`],
//! [`gamma_inc`], [`gamma_inc_inv`], [`psi`], [`error_f`], [`error_fc`],
//! [`error_fc_scaled`], [`cumnor`], [`dinvnr`], [`dlanor`], and [`dt1`].
//!
//! The functions that can fail have a fallible `try_*` form
//! ([`try_beta_inc`], [`try_gamma`], [`try_gamma_inc`], [`try_gamma_inc_inv`],
//! [`try_psi`]), except [`dlanor`], which panics when its argument is outside
//! the range of its asymptotic expansion; [`gamma_inc_with_acc`] exposes the
//! accuracy selector of CDFLIB's `gamma_inc`.
//!
//! The companion [`internal`] submodule exposes the CDFLIB-style helper
//! routines used inside the routines above ([`algdiv`], [`bcorr`], [`gam1`],
//! [`rlog`], etc.). They are public so users porting C/Fortran code that calls these
//! directly can find each routine under its CDFLIB name, but they are not part
//! of the user-facing statistical API.
//!
//! The two-output (*cum*, *ccum*) convention from CDFLIB is preserved on the
//! routines that drive distribution tail accuracy: [`cumnor`], [`gamma_inc`],
//! and [`beta_inc`]. Returning both tail probabilities directly is essential to
//! the library's tail accuracy.
//!
//! CDFLIB's other `cum*` helpers (`cumbet`, `cumbin`, `cumchi`, `cumchn`,
//! `cumf`, `cumfnc`, `cumgam`, `cumnbn`, `cumpoi`, `cumt`) are folded into the
//! corresponding distribution modules and are not exposed here; if you want
//! their behavior, use the distribution's [`ContinuousCdf::cdf`] /
//! [`ContinuousCdf::ccdf`] or [`DiscreteCdf::cdf`] / [`DiscreteCdf::ccdf`]
//! methods.
//!
//! [`beta`]: crate::special::beta()
//! [`beta_log`]: crate::special::beta_log
//! [`beta_inc`]: crate::special::beta_inc
//! [`gamma`]: crate::special::gamma()
//! [`gamma_log`]: crate::special::gamma_log
//! [`gamma_inc`]: crate::special::gamma_inc
//! [`gamma_inc_inv`]: crate::special::gamma_inc_inv
//! [`psi`]: crate::special::psi
//! [`error_f`]: crate::special::error_f
//! [`error_fc`]: crate::special::error_fc
//! [`error_fc_scaled`]: crate::special::error_fc_scaled
//! [`cumnor`]: crate::special::cumnor
//! [`dinvnr`]: crate::special::dinvnr
//! [`dlanor`]: crate::special::dlanor
//! [`dt1`]: crate::special::dt1
//! [`try_beta_inc`]: crate::special::try_beta_inc
//! [`try_gamma`]: crate::special::try_gamma
//! [`try_gamma_inc`]: crate::special::try_gamma_inc
//! [`try_gamma_inc_inv`]: crate::special::try_gamma_inc_inv
//! [`try_psi`]: crate::special::try_psi
//! [`gamma_inc_with_acc`]: crate::special::gamma_inc_with_acc
//! [`internal`]: crate::special::internal
//! [`algdiv`]: crate::special::internal::algdiv
//! [`bcorr`]: crate::special::internal::bcorr
//! [`gam1`]: crate::special::internal::gam1
//! [`rlog`]: crate::special::internal::rlog
//! [`ContinuousCdf::cdf`]: crate::traits::ContinuousCdf::cdf
//! [`ContinuousCdf::ccdf`]: crate::traits::ContinuousCdf::ccdf
//! [`DiscreteCdf::cdf`]: crate::traits::DiscreteCdf::cdf
//! [`DiscreteCdf::ccdf`]: crate::traits::DiscreteCdf::ccdf

pub(crate) mod beta;
pub(crate) mod erf;
pub(crate) mod gamma;
pub mod internal;
pub(crate) mod normal;
pub(crate) mod students_t;
pub use beta::{beta, beta_inc, beta_log, try_beta_inc, BetaIncError};
pub use erf::{error_f, error_fc, error_fc_scaled};
pub use gamma::{
    gamma, gamma_inc, gamma_inc_inv, gamma_inc_with_acc, gamma_log, psi, try_gamma, try_gamma_inc,
    try_gamma_inc_inv, try_gamma_inc_with_acc, try_psi, GammaDomainError, GammaIncAcc,
    GammaIncError, GammaIncInvError, PsiError,
};
pub use normal::{cumnor, dinvnr, dlanor};
pub use students_t::dt1;

/// Returns *x*², the F90 `x**2`, which gfortran computes as the product
/// `x * x`.
#[inline]
pub(crate) fn pow2(x: f64) -> f64 {
    x * x
}

/// Evaluates the polynomial *a*₀ + *a*₁·*x* + … + *aₙ*·*xⁿ* by Horner's
/// method (cdflib.f90:9709). The degree *n* is `a.len() - 1`.
#[inline]
#[allow(clippy::needless_range_loop)]
pub(crate) fn eval_pol(a: &[f64], x: f64) -> f64 {
    let n = a.len() - 1;
    let mut term = a[n];
    for i in (0..n).rev() {
        term = term * x + a[i];
    }
    term
}

/// Returns the integer machine constant number *i* (cdflib.f90:12551).
///
/// The table is the IEEE 754 binary64 configuration of CDFLIB's `imach`
/// array: 1 = base of integer arithmetic, 2 = number of base-2 digits of an
/// integer, 3 = largest integer, 4 = base of floating-point arithmetic,
/// 5-7 = digits, minimum and maximum exponent of single precision,
/// 8-10 = digits, minimum and maximum exponent of double precision.
///
/// # Panics
///
/// Panics if *i* is not in [1 . . 10].
#[inline]
pub(crate) fn ipmpar(i: usize) -> i32 {
    const IMACH: [i32; 10] = [2, 31, 2147483647, 2, 24, -125, 128, 53, -1021, 1024];
    IMACH[i - 1]
}

/// Returns the largest positive *w* for which exp(*w*) can be computed
/// (*l* = 0), or the largest negative *w* for which the computed value of
/// exp(*w*) is nonzero (*l* ≠ 0); only an approximate value is returned
/// (cdflib.f90:9761).
#[inline]
#[allow(clippy::approx_constant)]
pub(crate) fn exparg(l: i32) -> f64 {
    // Get the arithmetic base.
    let b = ipmpar(4);
    // Compute the logarithm of the arithmetic base.
    let lnb = if b == 2 {
        0.69314718055995
    } else if b == 8 {
        2.0794415416798
    } else if b == 16 {
        2.7725887222398
    } else {
        (b as f64).ln()
    };

    if l != 0 {
        let m = ipmpar(9) - 1;
        0.99999 * (m as f64 * lnb)
    } else {
        let m = ipmpar(10);
        0.99999 * (m as f64 * lnb)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // tests/data/exparg.csv is written by tests/regenerate/gen_kernel_coverage.f90.
    #[test]
    #[cfg(not(miri))]
    fn exparg_matches_reference() {
        let text = std::fs::read_to_string("tests/data/exparg.csv").unwrap();
        let mut rows = 0;
        for line in text.lines().filter(|l| !l.starts_with('#')) {
            let v: Vec<f64> = line.split(',').map(|t| t.trim().parse().unwrap()).collect();
            assert_eq!(exparg(v[0] as i32).to_bits(), v[1].to_bits(), "{line}");
            rows += 1;
        }
        assert_eq!(rows, 3);
    }
}
