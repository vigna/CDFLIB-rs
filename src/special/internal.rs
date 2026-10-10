//! CDFLIB helper routines that the user-facing routines build on.
//!
//! Every function here is part of the public API and is documented with its
//! CDFLIB role; they live in a separate module so the user-facing
//! [`cdflib::special`] namespace stays focused on the
//! routines a statistical user is likely to call directly ([`beta_inc`],
//! [`gamma_inc`], [`error_f`], [`cumnor`], etc.).
//!
//! Users porting code that calls these helpers by name (for example
//! [`algdiv`], [`bcorr`], [`gam1`], [`rlog`], [`stvaln`]) can find
//! each routine here under its CDFLIB name. Each one is a port of the routine
//! of the same name in the Fortran 90 CDFLIB, and its documentation follows
//! the Fortran header; `ierr` out-parameters are surfaced through the
//! matching Rust `Result` types (see for example [`BetaGratError`]).
//!
//! [`beta_inc`]: crate::special::beta_inc
//! [`gamma_inc`]: crate::special::gamma_inc
//! [`error_f`]: crate::special::error_f
//! [`cumnor`]: crate::special::cumnor
//! [`algdiv`]: crate::special::internal::algdiv
//! [`bcorr`]: crate::special::internal::bcorr
//! [`gam1`]: crate::special::internal::gam1
//! [`rlog`]: crate::special::internal::rlog
//! [`stvaln`]: crate::special::internal::stvaln
//! [`BetaGratError`]: crate::special::internal::BetaGratError
//! [`cdflib::special`]: crate::special

pub use super::beta::{
    algdiv, apser, bcorr, beta_asym, beta_frac, beta_grat, beta_pser, beta_rcomp, beta_rcomp1,
    beta_up, dbetrm, esum, fpser, gamma_rat1, BetaGratError,
};
pub use super::gamma::{alnrel, dexpm1, dstrem, gam1, gamma_ln1, gsumln, rcomp, rexp, rlog, rlog1};
pub use super::normal::stvaln;
