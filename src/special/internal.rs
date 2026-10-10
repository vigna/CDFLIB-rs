//! CDFLIB helper routines that the user-facing routines build on.
//!
//! Each function here ([`algdiv`], [`bcorr`], [`gam1`], [`rlog`],
//! [`stvaln`], etc.) is a port of the routine of the same name in the
//! Fortran 90 CDFLIB, documented after its Fortran header; `ierr`
//! out-parameters become `Result` types (see for example
//! [`BetaGratError`]). They live here so that [`cdflib::special`] keeps to
//! the routines a statistical user is likely to call directly.
//!
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
