//! Γ function family: ln Γ, Γ, ψ, and the incomplete Γ ratio functions
//! *P*(*a*, *x*) and *Q*(*a*, *x*).
//!
//! [`gamma_inc`] selects among the algorithms of CDFLIB's `gamma_inc`:
//!
//! - **Taylor series** for *P*/*R* and for *P*(*a*, *x*)/*x*ᵃ.
//!
//! - **Continued fraction expansion**.
//!
//! - **Asymptotic expansion** and the **Temme expansion** for large *a*.
//!
//! - **Finite sums for *Q*** when 1 ≤ *a* and 2*a* is an integer.
//!
//! - **Special cases** for *a* = 1/2 (reduces to [`error_f`] / [`error_fc`])
//!   and for *a*·*x* = 0.
//!
//! [`gamma_inc`]: crate::special::gamma_inc
//! [`error_f`]: crate::special::error_f
//! [`error_fc`]: crate::special::error_fc

#![allow(clippy::approx_constant, clippy::excessive_precision)]

use super::erf::{error_f, error_fc, error_fc_scaled};
use super::{eval_pol, pow2};

// =====================================================================
// Low-level helpers
// =====================================================================

/// Evaluates the function ln(1 + *a*) (cdflib.f90:118).
///
/// For |*a*| ≤ 0.375 a rational approximation in *t* = *a*/(*a* + 2) is
/// used; otherwise the logarithm of 1 + *a* is computed directly.
#[inline]
pub fn alnrel(a: f64) -> f64 {
    const P1: f64 = -0.129418923021993e1;
    const P2: f64 = 0.405303492862024;
    const P3: f64 = -0.178874546012214e-1;
    const Q1: f64 = -0.162752256355323e1;
    const Q2: f64 = 0.747811014037616;
    const Q3: f64 = -0.845104217945565e-1;

    if a.abs() <= 0.375 {
        let t = a / (a + 2.0);
        let t2 = t * t;

        let w = (((P3 * t2 + P2) * t2 + P1) * t2 + 1.0) / (((Q3 * t2 + Q2) * t2 + Q1) * t2 + 1.0);

        2.0 * t * w
    } else {
        let x = 1.0 + a;
        x.ln()
    }
}

/// Evaluates the function exp(*x*) − 1 (cdflib.f90:13887).
///
/// For |*x*| ≤ 0.15 a rational approximation is used, avoiding the
/// cancellation of exp(*x*) − 1 near zero.
#[inline]
pub fn rexp(x: f64) -> f64 {
    const P1: f64 = 0.914041914819518e-9;
    const P2: f64 = 0.238082361044469e-1;
    const Q1: f64 = -0.499999999085958;
    const Q2: f64 = 0.107141568980644;
    const Q3: f64 = -0.119041179760821e-1;
    const Q4: f64 = 0.595130811860248e-3;

    if x.abs() <= 0.15 {
        x * (((P2 * x + P1) * x + 1.0) / ((((Q4 * x + Q3) * x + Q2) * x + Q1) * x + 1.0))
    } else {
        let w = x.exp();

        if x <= 0.0 {
            (w - 0.5) - 0.5
        } else {
            w * (0.5 + (0.5 - 1.0 / w))
        }
    }
}

/// Evaluates the function exp(*x*) − 1 (cdflib.f90:7970).
///
/// Same algorithm as [`rexp`]: a rational approximation for |*x*| ≤ 0.15,
/// exp(*x*) otherwise.
///
/// # Example
///
/// ```
/// use cdflib::special::internal::dexpm1;
///
/// // dexpm1(0.0) = 0 exactly.
/// assert_eq!(dexpm1(0.0), 0.0);
/// // dexpm1 matches f64::exp(x) - 1 to ~1 ULP for moderate x.
/// let y = dexpm1(0.5);
/// # #[cfg(not(miri))]
/// assert!((y - (0.5_f64.exp() - 1.0)).abs() < 1e-15);
/// ```
///
/// [`rexp`]: crate::special::internal::rexp
#[inline]
pub fn dexpm1(x: f64) -> f64 {
    const P1: f64 = 0.914041914819518e-9;
    const P2: f64 = 0.238082361044469e-1;
    const Q1: f64 = -0.499999999085958;
    const Q2: f64 = 0.107141568980644;
    const Q3: f64 = -0.119041179760821e-1;
    const Q4: f64 = 0.595130811860248e-3;

    if x.abs() <= 0.15 {
        let top = (P2 * x + P1) * x + 1.0;
        let bot = (((Q4 * x + Q3) * x + Q2) * x + Q1) * x + 1.0;
        x * (top / bot)
    } else {
        let w = x.exp();

        if x <= 0.0 {
            (w - 0.5) - 0.5
        } else {
            w * (0.5 + (0.5 - 1.0 / w))
        }
    }
}

/// Computes *x* − 1 − ln(*x*) (cdflib.f90:13954).
#[inline]
#[allow(clippy::assign_op_pattern, clippy::needless_late_init)]
pub fn rlog(x: f64) -> f64 {
    const A: f64 = 0.566749439387324e-1;
    const B: f64 = 0.456512608815524e-1;
    const HALF: f64 = 0.5;
    const P0: f64 = 0.333333333333333;
    const P1: f64 = -0.224696413112536;
    const P2: f64 = 0.620886815375787e-2;
    const Q1: f64 = -0.127408923933623e1;
    const Q2: f64 = 0.354508718369557;
    const TWO: f64 = 2.0;

    if x < 0.61 {
        let r = (x - 0.5) - 0.5;
        r - x.ln()
    } else if x < 1.57 {
        let mut u;
        let w1;

        if x < 0.82 {
            u = x - 0.7;
            u = u / 0.7;
            w1 = A - u * 0.3;
        } else if x < 1.18 {
            u = (x - HALF) - HALF;
            w1 = 0.0;
        } else {
            // cdflib.f90 tests x < 1.57 again here; it always holds.
            u = 0.75 * x - 1.0;
            w1 = B + u / 3.0;
        }

        let r = u / (u + TWO);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        TWO * t * (1.0 / (1.0 - r) - r * w) + w1
    } else if 1.57 <= x {
        let r = (x - HALF) - HALF;
        r - x.ln()
    } else {
        // Rust only: x is NaN, for which cdflib.f90 leaves rlog unassigned.
        f64::NAN
    }
}

/// Evaluates the function *x* − ln(1 + *x*) (cdflib.f90:14048).
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn rlog1(x: f64) -> f64 {
    const A: f64 = 0.566749439387324e-1;
    const B: f64 = 0.456512608815524e-1;
    const HALF: f64 = 0.5;
    const P0: f64 = 0.333333333333333;
    const P1: f64 = -0.224696413112536;
    const P2: f64 = 0.620886815375787e-2;
    const Q1: f64 = -0.127408923933623e1;
    const Q2: f64 = 0.354508718369557;
    const TWO: f64 = 2.0;

    if x < -0.39 {
        let w = (x + HALF) + HALF;
        x - w.ln()
    } else if x < -0.18 {
        let mut h = x + 0.3;
        h = h / 0.7;
        let w1 = A - h * 0.3;

        let r = h / (h + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        TWO * t * (1.0 / (1.0 - r) - r * w) + w1
    } else if x <= 0.18 {
        let h = x;
        let w1 = 0.0;

        let r = h / (h + TWO);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        TWO * t * (1.0 / (1.0 - r) - r * w) + w1
    } else if x <= 0.57 {
        let h = 0.75 * x - 0.25;
        let w1 = B + h / 3.0;

        let r = h / (h + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        TWO * t * (1.0 / (1.0 - r) - r * w) + w1
    } else {
        let w = (x + HALF) + HALF;
        x - w.ln()
    }
}

/// Computes 1/Γ(*a* + 1) − 1 for −0.5 ≤ *a* ≤ 1.5 (cdflib.f90:10176).
#[inline]
pub fn gam1(a: f64) -> f64 {
    const P: [f64; 7] = [
        0.577215664901533,
        -0.409078193005776,
        -0.230975380857675,
        0.597275330452234e-1,
        0.766968181649490e-2,
        -0.514889771323592e-2,
        0.589597428611429e-3,
    ];
    const Q: [f64; 5] = [
        0.100000000000000e1,
        0.427569613095214,
        0.158451672430138,
        0.261132021441447e-1,
        0.423244297896961e-2,
    ];
    const R: [f64; 9] = [
        -0.422784335098468,
        -0.771330383816272,
        -0.244757765222226,
        0.118378989872749,
        0.930357293360349e-3,
        -0.118290993445146e-1,
        0.223047661158249e-2,
        0.266505979058923e-3,
        -0.132674909766242e-3,
    ];
    const S1: f64 = 0.273076135303957;
    const S2: f64 = 0.559398236957378e-1;

    let d = a - 0.5;

    let t = if 0.0 < d { d - 0.5 } else { a };

    if t == 0.0 {
        0.0
    } else if 0.0 < t {
        let top = (((((P[6] * t + P[5]) * t + P[4]) * t + P[3]) * t + P[2]) * t + P[1]) * t + P[0];

        let bot = (((Q[4] * t + Q[3]) * t + Q[2]) * t + Q[1]) * t + 1.0;

        let w = top / bot;

        if d <= 0.0 {
            a * w
        } else {
            (t / a) * ((w - 0.5) - 0.5)
        }
    } else if t < 0.0 {
        let top = (((((((R[8] * t + R[7]) * t + R[6]) * t + R[5]) * t + R[4]) * t + R[3]) * t
            + R[2])
            * t
            + R[1])
            * t
            + R[0];

        let bot = (S2 * t + S1) * t + 1.0;
        let w = top / bot;

        if d <= 0.0 {
            a * ((w + 0.5) + 0.5)
        } else {
            t * w / a
        }
    } else {
        // Rust only: a is NaN, for which cdflib.f90 leaves gam1 unassigned.
        f64::NAN
    }
}

// =====================================================================
// Γ(a): the Γ function itself
// =====================================================================

/// Errors of [`try_gamma`].
///
/// CDFLIB's `gamma_user` (cdflib.f90:10300) sets its result to 0 on entry
/// (cdflib.f90:10379) and returns that sentinel on each path where Γ(*a*)
/// cannot be computed; each such path maps onto one variant.
///
/// [`try_gamma`]: crate::special::try_gamma
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaDomainError {
    /// Argument is zero or a negative integer; Γ has a pole there.
    #[error("Γ has a pole at {0}")]
    Pole(f64),
    /// Result would overflow f64: 1/*t* can overflow for tiny |*a*|,
    /// *a* ≥ 1000, or exp(*w*) in the final assembly would overflow for
    /// positive *a* (beyond about 171.6).
    #[error("Γ({0}) overflows f64")]
    Overflow(f64),
    /// Result would underflow f64: the argument is too negative. CDFLIB's
    /// `gamma_user` returns its sentinel 0 for every *a* ≤ −1000
    /// (cdflib.f90:10457-10459), and for negative *a* beyond about −171.6
    /// when exp(*w*), the magnitude of Γ(−*a*), would overflow in the
    /// final assembly (cdflib.f90:10499-10501). Γ(*a*) underflows there,
    /// except at the negative integers, which are poles and get the same
    /// sentinel for *a* ≤ −1000.
    #[error("Γ({0}) underflows f64")]
    Underflow(f64),
}

/// Evaluates the Γ function.
///
/// # Panics
///
/// Panics on a [`GammaDomainError`] (pole at zero or a negative integer
/// greater than −1000, overflow, or underflow for negative *a* beyond
/// about −171.6). Use [`try_gamma`] for the fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::gamma;
///
/// let y = gamma(3.0);
/// assert!((y - 2.0).abs() < 1e-14);
/// ```
///
/// [`try_gamma`]: crate::special::try_gamma
#[inline]
pub fn gamma(a: f64) -> f64 {
    try_gamma(a).unwrap_or_else(|e| panic!("gamma({a}): {e}"))
}

/// Evaluates the Γ function, port of CDFLIB's `gamma_user`
/// (cdflib.f90:10300).
///
/// Returns [`GammaDomainError`] where `gamma_user` returns its error
/// sentinel 0: at a pole, or when the result would overflow or underflow.
///
/// # Example
///
/// ```
/// use cdflib::special::{try_gamma, GammaDomainError};
///
/// assert!(matches!(try_gamma(-3.0), Err(GammaDomainError::Pole(_))));
/// assert!(matches!(try_gamma(2000.0), Err(GammaDomainError::Overflow(_))));
/// assert!((try_gamma(3.0).unwrap() - 2.0).abs() < 1e-14);
/// ```
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn try_gamma(a: f64) -> Result<f64, GammaDomainError> {
    const D: f64 = 0.41893853320467274178;
    const P: [f64; 7] = [
        0.539637273585445e-3,
        0.261939260042690e-2,
        0.204493667594920e-1,
        0.730981088720487e-1,
        0.279648642639792,
        0.553413866010467,
        1.0,
    ];
    const PI: f64 = 3.1415926535898;
    const Q: [f64; 7] = [
        -0.832979206704073e-3,
        0.470059485860584e-2,
        0.225211131035340e-1,
        -0.170458969313360,
        -0.567902761974940e-1,
        0.113062953091122e1,
        1.0,
    ];
    const R1: f64 = 0.820756370353826e-3;
    const R2: f64 = -0.595156336428591e-3;
    const R3: f64 = 0.793650663183693e-3;
    const R4: f64 = -0.277777777770481e-2;
    const R5: f64 = 0.833333333333333e-1;

    // cdflib.f90:10379 presets gamma_user to the error sentinel 0; here each
    // path that returns it returns a GammaDomainError instead.
    let mut x = a;

    if a.abs() < 15.0 {
        // Evaluation of Γ(a) for |a| < 15.
        let mut t = 1.0;
        let mut m = a as i32 - 1;

        // Let t be the product of a - j when 2 <= a.
        if 0 <= m {
            for _j in 1..=m {
                x = x - 1.0;
                t = x * t;
            }

            x = x - 1.0;
        } else {
            // Let t be the product of a + j when a < 1.
            t = a;

            if a <= 0.0 {
                m = -m - 1;

                for _j in 1..=m {
                    x = x + 1.0;
                    t = x * t;
                }

                x = (x + 0.5) + 0.5;
                t = x * t;
                if t == 0.0 {
                    return Err(GammaDomainError::Pole(a));
                }
            }

            // Check if 1/t can overflow.
            if t.abs() < 1.0e-30 {
                if 1.0001 < t.abs() * f64::MAX {
                    return Ok(1.0 / t);
                }
                return Err(GammaDomainError::Overflow(a));
            }
        }

        // Compute Γ(1 + x) for 0 <= x < 1.
        let mut top = P[0];
        let mut bot = Q[0];
        for i in 1..7 {
            top = top * x + P[i];
            bot = bot * x + Q[i];
        }

        let mut gamma_user = top / bot;

        // Termination.
        if 1.0 <= a {
            gamma_user = gamma_user * t;
        } else {
            gamma_user = gamma_user / t;
        }

        Ok(gamma_user)
    } else {
        // Evaluation of Γ(a) for 15 <= |a|.
        if 1000.0 <= a.abs() {
            // Rust only: cdflib.f90 returns the same sentinel for both signs;
            // the sign of a selects the error variant.
            if a < 0.0 {
                return Err(GammaDomainError::Underflow(a));
            }
            return Err(GammaDomainError::Overflow(a));
        }

        let mut t;
        // Rust only: the F90 leaves s unset when 0 < a; it is read only
        // when a < 0, after the reflection below sets it.
        let mut s = 0.0;

        if a <= 0.0 {
            x = -a;
            let n = x as i32;
            t = x - n as f64;

            if 0.9 < t {
                t = 1.0 - t;
            }

            s = (PI * t).sin() / PI;

            if n % 2 == 0 {
                s = -s;
            }

            if s == 0.0 {
                return Err(GammaDomainError::Pole(a));
            }
        }

        // Compute the modified asymptotic sum.
        t = 1.0 / (x * x);

        let mut g = ((((R1 * t + R2) * t + R3) * t + R4) * t + R5) / x;

        let lnx = x.ln();

        // Final assembly.
        let z = x;
        g = (D + g) + (z - 0.5) * (lnx - 1.0);
        let w = g;
        t = g - w;

        if 0.99999 * super::exparg(0) < w {
            // Rust only: exp(w) is the magnitude of Γ(-a) when a < 0, and
            // then Γ(a) underflows; the sign of a selects the error variant.
            if a < 0.0 {
                return Err(GammaDomainError::Underflow(a));
            }
            return Err(GammaDomainError::Overflow(a));
        }

        let mut gamma_user = w.exp() * (1.0 + t);

        if a < 0.0 {
            gamma_user = (1.0 / (gamma_user * s)) / x;
        }

        Ok(gamma_user)
    }
}

/// Evaluates ln(Γ(1 + *a*)), for −0.2 ≤ *a* ≤ 1.25 (cdflib.f90:12016).
#[inline]
pub fn gamma_ln1(a: f64) -> f64 {
    const P0: f64 = 0.577215664901533;
    const P1: f64 = 0.844203922187225;
    const P2: f64 = -0.168860593646662;
    const P3: f64 = -0.780427615533591;
    const P4: f64 = -0.402055799310489;
    const P5: f64 = -0.673562214325671e-1;
    const P6: f64 = -0.271935708322958e-2;
    const Q1: f64 = 0.288743195473681e1;
    const Q2: f64 = 0.312755088914843e1;
    const Q3: f64 = 0.156875193295039e1;
    const Q4: f64 = 0.361951990101499;
    const Q5: f64 = 0.325038868253937e-1;
    const Q6: f64 = 0.667465618796164e-3;
    const R0: f64 = 0.422784335098467;
    const R1: f64 = 0.848044614534529;
    const R2: f64 = 0.565221050691933;
    const R3: f64 = 0.156513060486551;
    const R4: f64 = 0.170502484022650e-1;
    const R5: f64 = 0.497958207639485e-3;
    const S1: f64 = 0.124313399877507e1;
    const S2: f64 = 0.548042109832463;
    const S3: f64 = 0.101552187439830;
    const S4: f64 = 0.713309612391000e-2;
    const S5: f64 = 0.116165475989616e-3;

    if a < 0.6 {
        let top = (((((P6 * a + P5) * a + P4) * a + P3) * a + P2) * a + P1) * a + P0;

        let bot = (((((Q6 * a + Q5) * a + Q4) * a + Q3) * a + Q2) * a + Q1) * a + 1.0;

        -(a * (top / bot))
    } else {
        let x = (a - 0.5) - 0.5;

        let top = ((((R5 * x + R4) * x + R3) * x + R2) * x + R1) * x + R0;

        let bot = ((((S5 * x + S4) * x + S3) * x + S2) * x + S1) * x + 1.0;

        x * (top / bot)
    }
}

/// Errors of [`try_psi`].
///
/// CDFLIB's `psi` (cdflib.f90:13446) documents that “PSI is assigned the
/// value 0 when the psi function is undefined”; each path that returns
/// that sentinel maps onto one variant.
///
/// [`try_psi`]: crate::special::try_psi
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum PsiError {
    /// Argument is zero or a negative integer; ψ has a pole there.
    #[error("ψ has a pole at {0}")]
    Pole(f64),
    /// Argument is at or below −*xmax1*, where *xmax1* = 2³¹ − 1 is the
    /// largest integer (CDFLIB's `ipmpar(3)`), the lower bound on
    /// acceptable negative arguments.
    #[error("ψ argument {0} is too large in magnitude (overflow)")]
    Overflow(f64),
}

/// Evaluates the ψ or digamma function, d/d*x* ln Γ(*x*).
///
/// # Panics
///
/// Panics on a [`PsiError`] (pole or overflow). Use [`try_psi`] for the
/// fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::psi;
///
/// // ψ(1) = −γ (Euler–Mascheroni constant)
/// let y = psi(1.0);
/// assert!((y + 0.57721566).abs() < 1e-8);
/// ```
///
/// [`try_psi`]: crate::special::try_psi
#[inline]
pub fn psi(xx: f64) -> f64 {
    try_psi(xx).unwrap_or_else(|e| panic!("psi({xx}): {e}"))
}

/// Evaluates the ψ or digamma function, d/d*x* ln Γ(*x*), port of CDFLIB's
/// `psi` (cdflib.f90:13446).
///
/// The main computation involves evaluation of rational Chebyshev
/// approximations. `psi` was written at Argonne National Laboratory for
/// FUNPACK, and subsequently modified by A. H. Morris of NSWC.
///
/// Returns [`PsiError`] where CDFLIB returns its sentinel 0 because ψ is
/// undefined: at a pole, or when *xx* ≤ −*xmax1*.
///
/// # References
///
/// William Cody, Anthony Strecok, Henry Thacher, “Chebyshev Approximations
/// for the Psi Function”, *Mathematics of Computation*, 27:123–127, 1973.
///
/// # Example
///
/// ```
/// use cdflib::special::{try_psi, PsiError};
///
/// assert!(matches!(try_psi(0.0), Err(PsiError::Pole(_))));
/// assert!((try_psi(1.0).unwrap() + 0.57721566).abs() < 1e-8);
/// ```
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn try_psi(xx: f64) -> Result<f64, PsiError> {
    const DX0: f64 = 1.461632144968362341262659542325721325;
    const P1: [f64; 7] = [
        0.895385022981970e-2,
        0.477762828042627e1,
        0.142441585084029e3,
        0.118645200713425e4,
        0.363351846806499e4,
        0.413810161269013e4,
        0.130560269827897e4,
    ];
    const P2: [f64; 4] = [
        -0.212940445131011e1,
        -0.701677227766759e1,
        -0.448616543918019e1,
        -0.648157123766197,
    ];
    const PIOV4: f64 = 0.785398163397448;
    // Coefficients for rational approximation of psi(x) / (x - x0),
    // 0.5 <= x <= 3.0.
    const Q1: [f64; 6] = [
        0.448452573429826e2,
        0.520752771467162e3,
        0.221000799247830e4,
        0.364127349079381e4,
        0.190831076596300e4,
        0.691091682714533e-5,
    ];
    const Q2: [f64; 4] = [
        0.322703493791143e2,
        0.892920700481861e2,
        0.546117738103215e2,
        0.777788548522962e1,
    ];

    // xmax1 is the largest positive floating point constant with entirely
    // integer representation. It is also used as negative of lower bound on
    // acceptable negative arguments and as the positive argument beyond
    // which psi may be represented as ln(x).
    let mut xmax1 = super::ipmpar(3) as f64;
    xmax1 = xmax1.min(1.0 / f64::EPSILON);

    // xsmall is the absolute argument below which pi * cotan(pi * x) may be
    // represented by 1/x.
    let xsmall = 1.0e-9;

    let mut x = xx;
    let mut aug = 0.0;

    if x == 0.0 {
        return Err(PsiError::Pole(xx));
    }

    // x < 0.5: use reflection formula psi(1 - x) = psi(x) + pi * cotan(pi * x).
    if x < 0.5 {
        'l40: {
            // 0 < |x| <= xsmall: use 1/x as a substitute for pi * cotan(pi * x).
            if x.abs() <= xsmall {
                aug = -1.0 / x;
                break 'l40;
            }

            // Reduction of argument for cotangent.
            let mut w = -x;
            let mut sgn = PIOV4;

            if w <= 0.0 {
                w = -w;
                sgn = -sgn;
            }

            // Make an error exit if x <= -xmax1.
            if xmax1 <= w {
                return Err(PsiError::Overflow(xx));
            }

            let mut nq = w as i32;
            w = w - nq as f64;
            nq = (w * 4.0) as i32;
            w = 4.0 * (w - nq as f64 * 0.25);

            // w is now related to the fractional part of 4 * x. Adjust
            // argument to correspond to values in first quadrant and
            // determine sign.
            let mut n = nq / 2;
            if n + n != nq {
                w = 1.0 - w;
            }

            let z = PIOV4 * w;
            let mut m = n / 2;

            if m + m != n {
                sgn = -sgn;
            }

            // Determine final value for -pi * cotan(pi * x).
            n = (nq + 1) / 2;
            m = n / 2;
            m = m + m;

            if m == n {
                if z == 0.0 {
                    return Err(PsiError::Pole(xx));
                }

                // black_box keeps the optimizer from merging sin and cos
                // into one sincos call, whose last bit can differ from the
                // separate calls that the F90 makes.
                aug = 4.0 * sgn * (z.cos() / std::hint::black_box(z).sin());
            } else {
                aug = 4.0 * sgn * (z.sin() / std::hint::black_box(z).cos());
            }
        }

        x = 1.0 - x;
    }

    if x <= 3.0 {
        // 0.5 <= x <= 3.
        let mut den = x;
        let mut upper = P1[0] * x;

        for i in 1..=5 {
            den = (den + Q1[i - 1]) * x;
            upper = (upper + P1[i]) * x;
        }

        den = (upper + P1[6]) / (den + Q1[5]);
        let xmx0 = x - DX0;
        Ok(den * xmx0 + aug)
    } else if x < xmax1 {
        // 3 < x < xmax1.
        let w = (1.0 / x) / x;
        let mut den = w;
        let mut upper = P2[0] * w;

        for i in 1..=3 {
            den = (den + Q2[i - 1]) * w;
            upper = (upper + P2[i]) * w;
        }

        aug = upper / (den + Q2[3]) - 0.5 / x + aug;
        Ok(aug + x.ln())
    } else {
        // xmax1 <= x.
        Ok(aug + x.ln())
    }
}

/// Evaluates ln(Γ(*a*)) for positive *a* (cdflib.f90:12120).
///
/// # Example
///
/// ```
/// use cdflib::special::gamma_log;
///
/// let y = gamma_log(3.0);
/// assert!((y - 2.0_f64.ln()).abs() < 1e-14);
/// ```
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn gamma_log(a: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-1;
    const C1: f64 = -0.277777777760991e-2;
    const C2: f64 = 0.793650666825390e-3;
    const C3: f64 = -0.595202931351870e-3;
    const C4: f64 = 0.837308034031215e-3;
    const C5: f64 = -0.165322962780713e-2;
    const D: f64 = 0.418938533204673;

    if a <= 0.8 {
        gamma_ln1(a) - a.ln()
    } else if a <= 2.25 {
        let t = (a - 0.5) - 0.5;
        gamma_ln1(t)
    } else if a < 10.0 {
        let n = (a - 1.25) as i32;
        let mut t = a;
        let mut w = 1.0;
        for _i in 1..=n {
            t = t - 1.0;
            w = t * w;
        }

        gamma_ln1(t - 1.0) + w.ln()
    } else {
        let t = pow2(1.0 / a);

        let w = (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a;

        (D + w) + (a - 0.5) * (a.ln() - 1.0)
    }
}

/// Evaluates the function ln(Γ(*a* + *b*)) (cdflib.f90:12489).
///
/// `gsumln` is used for 1 ≤ *a* ≤ 2 and 1 ≤ *b* ≤ 2.
#[inline]
pub fn gsumln(a: f64, b: f64) -> f64 {
    let x = a + b - 2.0;

    if x <= 0.25 {
        gamma_ln1(1.0 + x)
    } else if x <= 1.25 {
        gamma_ln1(x) + alnrel(x)
    } else {
        gamma_ln1(x - 1.0) + (x * (1.0 + x)).ln()
    }
}

/// Computes the Sterling remainder ln(Γ(*z*)) − Sterling(*z*)
/// (cdflib.f90:8657).
///
/// Sterling(*z*) is Sterling's approximation to ln(Γ(*z*)):
///
/// Sterling(*z*) = ln(√(2π)) + (*z* − 0.5) ln(*z*) − *z*.
///
/// If 6 < *z*, the routine uses 9 terms of a series in Bernoulli numbers,
/// with values calculated using Maple. Otherwise, the difference is
/// computed explicitly via [`gamma_log`].
///
/// # Panics
///
/// Panics if *z* ≤ 0 (CDFLIB prints a fatal error and stops).
///
/// # Example
///
/// ```
/// use cdflib::special::internal::dstrem;
///
/// // For large z, dstrem(z) ≈ 1/(12 z), the leading Bernoulli term.
/// let y = dstrem(100.0);
/// assert!((y - 1.0 / 1200.0).abs() < 1e-6);
/// ```
///
/// [`gamma_log`]: crate::special::gamma_log
#[inline]
pub fn dstrem(z: f64) -> f64 {
    const NCOEF: usize = 9;
    const COEF: [f64; NCOEF + 1] = [
        0.0,
        0.0833333333333333333333333333333,
        -0.00277777777777777777777777777778,
        0.000793650793650793650793650793651,
        -0.000595238095238095238095238095238,
        0.000841750841750841750841750841751,
        -0.00191752691752691752691752691753,
        0.00641025641025641025641025641026,
        -0.0295506535947712418300653594771,
        0.179644372368830573164938490016,
    ];
    const HLN2PI: f64 = 0.91893853320467274178;

    if z <= 0.0 {
        panic!("dstrem: argument z must be positive (got {z})");
    }

    if 6.0 < z {
        eval_pol(&COEF, 1.0 / pow2(z)) * z
    } else {
        let sterl = HLN2PI + (z - 0.5) * z.ln() - z;
        gamma_log(z) - sterl
    }
}

/// Evaluates exp(−*x*) · *x*ᵃ / Γ(*a*) (cdflib.f90:13806).
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn rcomp(a: f64, x: f64) -> f64 {
    // rt2pin = 1/sqrt(2 * pi).
    const RT2PIN: f64 = 0.398942280401433;

    if a < 20.0 {
        let t = a * x.ln() - x;

        if a < 1.0 {
            (a * t.exp()) * (1.0 + gam1(a))
        } else {
            t.exp() / gamma(a)
        }
    } else {
        let u = x / a;

        if u == 0.0 {
            0.0
        } else {
            let t = pow2(1.0 / a);
            let mut t1 = (((0.75 * t - 1.0) * t + 3.5) * t - 105.0) / (a * 1260.0);
            t1 = t1 - a * rlog(u);
            RT2PIN * a.sqrt() * t1.exp()
        }
    }
}

// =====================================================================
// gamma_inc: regularized incomplete gamma P(a,x), Q(a,x)
// =====================================================================

/// Accuracy request of [`gamma_inc_with_acc`] and [`try_gamma_inc_with_acc`].
///
/// This is the argument `ind` of CDFLIB's `gamma_inc`
/// (cdflib.f90:10537-10540), which selects `iop` = `ind` + 1. The
/// [`gamma_inc`] and [`try_gamma_inc`] entry points are equivalent to
/// passing [`GammaIncAcc::Max`].
///
/// [`gamma_inc`]: crate::special::gamma_inc
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
/// [`gamma_inc_with_acc`]: crate::special::gamma_inc_with_acc
/// [`try_gamma_inc_with_acc`]: crate::special::try_gamma_inc_with_acc
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GammaIncAcc {
    /// As much accuracy as possible (`ind` = 0).
    #[default]
    Max,
    /// To within 1 unit of the 6-th significant digit (`ind` = 1).
    Digits6,
    /// To within 1 unit of the 3rd significant digit (`ind` = 2).
    Digits3,
}

/// Errors of [`gamma_inc`].
///
/// CDFLIB's `gamma_inc` reports each of these conditions by setting `ans`
/// to 2 (cdflib.f90:10544-10548); [`try_gamma_inc`] and
/// [`try_gamma_inc_with_acc`] return them as errors.
///
/// [`gamma_inc`]: crate::special::gamma_inc
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
/// [`try_gamma_inc_with_acc`]: crate::special::try_gamma_inc_with_acc
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaIncError {
    /// *a* is negative.
    #[error("parameter a must be non-negative, got {0}")]
    ANegative(f64),
    /// *x* is negative.
    #[error("argument x must be non-negative, got {0}")]
    XNegative(f64),
    /// *a* and *x* are both 0.
    #[error("both a and x are zero")]
    BothZero,
    /// The answer is computationally indeterminate because *a* is extremely
    /// large and *x* is very close to *a*.
    #[error("indeterminate at a = {a}, x = {x} (deep asymptotic regime)")]
    Indeterminate { a: f64, x: f64 },
}

/// Evaluates the incomplete Γ ratio functions *P*(*a*, *x*) and
/// *Q*(*a*, *x*) with as much accuracy as possible.
///
/// *a* and *x* are the arguments of the incomplete Γ ratio. *a* and *x*
/// must be nonnegative, and they cannot both be zero. Returns
/// (*P*(*a*, *x*), *Q*(*a*, *x*)).
///
/// # Panics
///
/// Panics on a [`GammaIncError`]: if *a* or *x* is negative, if both are
/// 0, or when the answer is computationally indeterminate because *a* is
/// extremely large and *x* is very close to *a*. A NaN argument gives NaN
/// components, unless the other argument is negative. *x* = +∞ gives NaN
/// components too: as in the F90 for *a* ≥ 1, while for *a* < 1 the F90
/// never returns. The exception is *a* = 1/2, where the result is (1, 0);
/// *a* = +∞ with finite *x* gives (0, 1). Use [`try_gamma_inc`] for the
/// fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::gamma_inc;
///
/// let (p, q) = gamma_inc(2.5, 1.7);
/// assert!((p - 0.36143008).abs() < 1e-8);
/// assert!((q - 0.63856992).abs() < 1e-8);
/// ```
///
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
#[inline]
pub fn gamma_inc(a: f64, x: f64) -> (f64, f64) {
    gamma_inc_with_acc(a, x, GammaIncAcc::Max)
}

/// Evaluates the incomplete Γ ratio functions *P*(*a*, *x*) and
/// *Q*(*a*, *x*) with the accuracy requested by `accuracy`.
///
/// Lower accuracy uses shallower truncations of the Temme expansions and
/// different regime cutoffs. See [`gamma_inc`] for the arguments and the
/// result.
///
/// # Panics
///
/// Panics on a [`GammaIncError`], exactly where
/// [`try_gamma_inc_with_acc`] returns one. Use [`try_gamma_inc_with_acc`]
/// for the fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::{gamma_inc_with_acc, GammaIncAcc};
///
/// let (p, _) = gamma_inc_with_acc(2.5, 1.7, GammaIncAcc::Max);
/// let (p3, _) = gamma_inc_with_acc(2.5, 1.7, GammaIncAcc::Digits3);
/// assert!((p - p3).abs() < 1e-3);
/// ```
///
/// [`gamma_inc`]: crate::special::gamma_inc
/// [`try_gamma_inc_with_acc`]: crate::special::try_gamma_inc_with_acc
#[inline]
pub fn gamma_inc_with_acc(a: f64, x: f64, accuracy: GammaIncAcc) -> (f64, f64) {
    try_gamma_inc_with_acc(a, x, accuracy)
        .unwrap_or_else(|e| panic!("gamma_inc_with_acc({a}, {x}, {accuracy:?}): {e}"))
}

/// Fallible form of [`gamma_inc`]: returns [`GammaIncError`] where CDFLIB
/// sets `ans` to 2.
///
/// # Example
///
/// ```
/// use cdflib::special::{try_gamma_inc, GammaIncError};
///
/// let (p, q) = try_gamma_inc(2.5, 1.7).unwrap();
/// assert!((p - 0.36143008).abs() < 1e-8);
/// assert!((q - 0.63856992).abs() < 1e-8);
/// assert!(matches!(
///     try_gamma_inc(-1.0, 1.0),
///     Err(GammaIncError::ANegative(_)),
/// ));
/// ```
///
/// [`GammaIncError`]: crate::special::GammaIncError
#[inline]
pub fn try_gamma_inc(a: f64, x: f64) -> Result<(f64, f64), GammaIncError> {
    try_gamma_inc_with_acc(a, x, GammaIncAcc::Max)
}

// Data of gamma_inc (cdflib.f90:10564-10677). Rust arrays are 0-based, so
// F90 acc0(iop) is ACC0[iop - 1] and d0(1) is D0[0].
const ACC0: [f64; 3] = [5.0e-15, 5.0e-7, 5.0e-4];
const ALOG10: f64 = 2.30258509299405;
const RT2PIN: f64 = 0.398942280401433;
const RTPI: f64 = 1.77245385090552;
const BIG: [f64; 3] = [20.0, 14.0, 10.0];
const E00: [f64; 3] = [0.25e-3, 0.25e-1, 0.14];
const X00: [f64; 3] = [31.0, 17.0, 9.7];
const D0: [f64; 13] = [
    0.833333333333333e-1,
    -0.148148148148148e-1,
    0.115740740740741e-2,
    0.352733686067019e-3,
    -0.178755144032922e-3,
    0.391926317852244e-4,
    -0.218544851067999e-5,
    -0.185406221071516e-5,
    0.829671134095309e-6,
    -0.176659527368261e-6,
    0.670785354340150e-8,
    0.102618097842403e-7,
    -0.438203601845335e-8,
];
const D10: f64 = -0.185185185185185e-2;
const D1: [f64; 12] = [
    -0.347222222222222e-2,
    0.264550264550265e-2,
    -0.990226337448560e-3,
    0.205761316872428e-3,
    -0.401877572016461e-6,
    -0.180985503344900e-4,
    0.764916091608111e-5,
    -0.161209008945634e-5,
    0.464712780280743e-8,
    0.137863344691572e-6,
    -0.575254560351770e-7,
    0.119516285997781e-7,
];
const D20: f64 = 0.413359788359788e-2;
const D2: [f64; 10] = [
    -0.268132716049383e-2,
    0.771604938271605e-3,
    0.200938786008230e-5,
    -0.107366532263652e-3,
    0.529234488291201e-4,
    -0.127606351886187e-4,
    0.342357873409614e-7,
    0.137219573090629e-5,
    -0.629899213838006e-6,
    0.142806142060642e-6,
];
const D30: f64 = 0.649434156378601e-3;
const D3: [f64; 8] = [
    0.229472093621399e-3,
    -0.469189494395256e-3,
    0.267720632062839e-3,
    -0.756180167188398e-4,
    -0.239650511386730e-6,
    0.110826541153473e-4,
    -0.567495282699160e-5,
    0.142309007324359e-5,
];
const D40: f64 = -0.861888290916712e-3;
const D4: [f64; 6] = [
    0.784039221720067e-3,
    -0.299072480303190e-3,
    -0.146384525788434e-5,
    0.664149821546512e-4,
    -0.396836504717943e-4,
    0.113757269706784e-4,
];
const D50: f64 = -0.336798553366358e-3;
const D5: [f64; 4] = [
    -0.697281375836586e-4,
    0.277275324495939e-3,
    -0.199325705161888e-3,
    0.679778047793721e-4,
];
const D60: f64 = 0.531307936463992e-3;
const D6: [f64; 2] = [-0.592166437353694e-3, 0.270878209671804e-3];
const D70: f64 = 0.344367606892378e-3;

/// Fallible form of [`gamma_inc_with_acc`]; this is CDFLIB's `gamma_inc`
/// (cdflib.f90:10513-11305).
///
/// Evaluates the incomplete Γ ratio functions *P*(*a*, *x*) and
/// *Q*(*a*, *x*). *a* and *x* must be nonnegative, and they cannot both be
/// zero. On normal output returns (`ans`, `qans`) = (*P*(*a*, *x*),
/// *Q*(*a*, *x*)). Where CDFLIB sets `ans` to 2 (*a* or *x* is negative,
/// both are 0, or the answer is computationally indeterminate because *a*
/// is extremely large and *x* is very close to *a*) this function returns
/// the corresponding [`GammaIncError`].
///
/// `accuracy` is the accuracy request `ind`: as much accuracy as possible,
/// to within 1 unit of the 6-th significant digit, or to within 1 unit of
/// the 3rd significant digit.
///
/// [`gamma_inc_with_acc`]: crate::special::gamma_inc_with_acc
/// [`GammaIncError`]: crate::special::GammaIncError
#[allow(clippy::manual_range_contains)]
pub fn try_gamma_inc_with_acc(
    a: f64,
    x: f64,
    accuracy: GammaIncAcc,
) -> Result<(f64, f64), GammaIncError> {
    let e = f64::EPSILON;

    // F90 sets ans = 2 here (cdflib.f90:10681-10684); the error variant
    // names the negative argument.
    if a < 0.0 || x < 0.0 {
        if a < 0.0 {
            return Err(GammaIncError::ANegative(a));
        }
        return Err(GammaIncError::XNegative(x));
    }

    // F90 sets ans = 2 here (cdflib.f90:10686-10689).
    if a == 0.0 && x == 0.0 {
        return Err(GammaIncError::BothZero);
    }

    if a * x == 0.0 {
        if x <= a {
            return Ok((0.0, 1.0));
        } else {
            return Ok((1.0, 0.0));
        }
    }

    // iop = ind + 1, with GammaIncAcc standing for ind; F90 maps every ind
    // other than 0 and 1 to iop = 3 (cdflib.f90:10702-10703).
    let iop: usize = match accuracy {
        GammaIncAcc::Max => 1,
        GammaIncAcc::Digits6 => 2,
        GammaIncAcc::Digits3 => 3,
    };
    let acc = ACC0[iop - 1].max(e);
    let e0 = E00[iop - 1];
    let x0 = X00[iop - 1];

    // Select the appropriate algorithm.
    let r;
    'l40: {
        'l10: {
            if 1.0 <= a {
                break 'l10;
            }

            if a == 0.5 {
                return Ok(special_cases_390(x));
            }

            // Rust only: if a or x is NaN, or x is +inf, F90 loops forever
            // at label 160 or label 250; return NaN instead.
            if a.is_nan() || x.is_nan() || x == f64::INFINITY {
                return Ok((f64::NAN, f64::NAN));
            }

            if x < 1.1 {
                return Ok(taylor_series_for_p_over_x_pow_a(a, x, acc));
            }

            let t1 = a * x.ln() - x;
            let u = a * t1.exp();

            if u == 0.0 {
                return Ok((1.0, 0.0));
            }

            r = u * (1.0 + gam1(a));
            return Ok(continued_fraction_expansion(a, x, r, e, acc));
        }

        // Label 10.
        'l30: {
            if BIG[iop - 1] <= a {
                break 'l30;
            }

            'l20: {
                if x < a || x0 <= x {
                    break 'l20;
                }

                let twoa = a + a;
                let m = twoa as i32;

                if twoa == f64::from(m) {
                    let i = m / 2;
                    if a == f64::from(i) {
                        return Ok(finite_sums_for_q_210(x, i));
                    }
                    return Ok(finite_sums_for_q_220(x, i));
                }
            }

            // Label 20.
            let t1 = a * x.ln() - x;
            r = t1.exp() / gamma(a);
            break 'l40;
        }

        // Label 30.
        let l = x / a;

        if l == 0.0 {
            return Ok((0.0, 1.0));
        }

        let s = 0.5 + (0.5 - l);
        let z = rlog(l);
        if 700.0 / a <= z {
            return special_cases_410(a, x, s, e);
        }

        let y = a * z;
        let rta = a.sqrt();

        if s.abs() <= e0 / rta {
            return temme_expansion_for_l_eq_1(a, x, l, z, y, rta, e, iop);
        }

        if s.abs() <= 0.4 {
            return general_temme_expansion(a, x, l, s, z, y, rta, e, iop);
        }

        let t = pow2(1.0 / a);
        let mut t1 = (((0.75 * t - 1.0) * t + 3.5) * t - 105.0) / (a * 1260.0);
        t1 -= y;
        r = RT2PIN * rta * t1.exp();
    }

    // Label 40.
    if r == 0.0 {
        if x <= a {
            return Ok((0.0, 1.0));
        } else {
            return Ok((1.0, 0.0));
        }
    }

    // Rust only: when x is NaN, F90 reaches label 100 and loops forever,
    // unless big <= a, where it first computes z = rlog(NaN), which rlog
    // leaves unassigned (cdflib.f90:10772), so that its result is
    // undefined; the same happens when a = x = +inf. Return NaN instead.
    if x.is_nan() || (a == f64::INFINITY && x == f64::INFINITY) {
        return Ok((f64::NAN, f64::NAN));
    }

    if x <= a.max(ALOG10) {
        return Ok(taylor_series_for_p_over_r(a, x, r, acc));
    }

    if x < x0 {
        return Ok(continued_fraction_expansion(a, x, r, e, acc));
    }

    Ok(asymptotic_expansion(a, x, r, acc))
}

// Taylor series for P/R (label 50, cdflib.f90:10816-10861).
fn taylor_series_for_p_over_r(a: f64, x: f64, r: f64, acc: f64) -> (f64, f64) {
    let mut wk = [0.0_f64; 20];

    let mut apn = a + 1.0;
    let mut t = x / apn;
    wk[0] = t;

    let mut n: usize = 20;

    for i in 2..=20 {
        apn += 1.0;
        t *= x / apn;
        if t <= 1.0e-3 {
            n = i;
            break;
        }
        wk[i - 1] = t;
    }

    let mut sum1 = t;

    let tol = 0.5 * acc;

    loop {
        apn += 1.0;
        t *= x / apn;
        sum1 += t;

        if t <= tol {
            break;
        }
    }

    let n_max = n - 1;
    for _m in 1..=n_max {
        n -= 1;
        sum1 += wk[n - 1];
    }

    let ans = (r / a) * (1.0 + sum1);
    let qans = 0.5 + (0.5 - ans);
    (ans, qans)
}

// Asymptotic expansion (label 100, cdflib.f90:10862-10904).
fn asymptotic_expansion(a: f64, x: f64, r: f64, acc: f64) -> (f64, f64) {
    let mut wk = [0.0_f64; 20];

    let mut amn = a - 1.0;
    let mut t = amn / x;
    wk[0] = t;

    let mut n: usize = 20;

    for i in 2..=20 {
        amn -= 1.0;
        t *= amn / x;
        if t.abs() <= 1.0e-3 {
            n = i;
            break;
        }
        wk[i - 1] = t;
    }

    let mut sum1 = t;

    loop {
        if t.abs() <= acc {
            break;
        }

        amn -= 1.0;
        t *= amn / x;
        sum1 += t;
    }

    let n_max = n - 1;
    for _m in 1..=n_max {
        n -= 1;
        sum1 += wk[n - 1];
    }
    let qans = (r / x) * (1.0 + sum1);
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

// Taylor series for P(A,X)/X^A (label 160, cdflib.f90:10905-10972).
fn taylor_series_for_p_over_x_pow_a(a: f64, x: f64, acc: f64) -> (f64, f64) {
    let mut an: f64 = 3.0;
    let mut c = x;
    let mut sum1 = x / (a + 3.0);
    let tol = 3.0 * acc / (a + 1.0);

    loop {
        an += 1.0;
        c = -(c * (x / an));
        let t = c / (a + an);
        sum1 += t;

        if t.abs() <= tol {
            break;
        }
    }

    let j = a * x * ((sum1 / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));

    let z = a * x.ln();
    let h = gam1(a);
    let g = 1.0 + h;

    'l200: {
        'l190: {
            'l180: {
                if x < 0.25 {
                    break 'l180;
                }

                if a < x / 2.59 {
                    break 'l200;
                }

                break 'l190;
            }

            // Label 180.
            if -0.13394 < z {
                break 'l200;
            }
        }

        // Label 190.
        let w = z.exp();
        let ans = w * g * (0.5 + (0.5 - j));
        let qans = 0.5 + (0.5 - ans);
        return (ans, qans);
    }

    // Label 200.
    let l = rexp(z);
    let w = 0.5 + (0.5 + l);
    let qans = (w * j - l) * g - h;

    if qans < 0.0 {
        return (1.0, 0.0);
    }

    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

// Finite sums for Q when 1 <= A and 2*A is an integer: label 210, the
// entry for integer a (cdflib.f90:10973-10982).
fn finite_sums_for_q_210(x: f64, i: i32) -> (f64, f64) {
    let sum1 = (-x).exp();
    let t = sum1;
    let n = 1;
    let c = 0.0;
    finite_sums_for_q_230(x, i, sum1, t, n, c)
}

// Finite sums for Q when 1 <= A and 2*A is an integer: label 220, the
// entry for half-integer a (cdflib.f90:10984-10990).
fn finite_sums_for_q_220(x: f64, i: i32) -> (f64, f64) {
    let rtx = x.sqrt();
    let sum1 = error_fc(rtx);
    let t = (-x).exp() / (RTPI * rtx);
    let n = 0;
    let c = -0.5;
    finite_sums_for_q_230(x, i, sum1, t, n, c)
}

// Finite sums for Q when 1 <= A and 2*A is an integer: label 230, reached
// from labels 210 and 220 (cdflib.f90:10992-11003).
fn finite_sums_for_q_230(
    x: f64,
    i: i32,
    mut sum1: f64,
    mut t: f64,
    mut n: i32,
    mut c: f64,
) -> (f64, f64) {
    while n != i {
        n += 1;
        c += 1.0;
        t = (x * t) / c;
        sum1 += t;
    }

    let qans = sum1;
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

// Continued fraction expansion (label 250, cdflib.f90:11004-11035).
fn continued_fraction_expansion(a: f64, x: f64, r: f64, e: f64, acc: f64) -> (f64, f64) {
    let tol = (5.0 * e).max(acc);
    let mut a2nm1: f64 = 1.0;
    let mut a2n: f64 = 1.0;
    let mut b2nm1 = x;
    let mut b2n = x + (1.0 - a);
    let mut c: f64 = 1.0;
    let mut an0;

    loop {
        a2nm1 = x * a2n + c * a2nm1;
        b2nm1 = x * b2n + c * b2nm1;
        let am0 = a2nm1 / b2nm1;
        c += 1.0;
        let cma = c - a;
        a2n = a2nm1 + cma * a2n;
        b2n = b2nm1 + cma * b2n;
        an0 = a2n / b2n;

        if (an0 - am0).abs() < tol * an0 {
            break;
        }
    }

    let qans = r * an0;
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

// General Temme expansion (label 270, cdflib.f90:11036-11189), ending at
// label 310. F90 error_fc(1, .) is error_fc_scaled.
#[allow(clippy::too_many_arguments)]
fn general_temme_expansion(
    a: f64,
    x: f64,
    l: f64,
    s: f64,
    mut z: f64,
    y: f64,
    rta: f64,
    e: f64,
    iop: usize,
) -> Result<(f64, f64), GammaIncError> {
    // F90 sets ans = 2 here (cdflib.f90:11041-11044).
    if s.abs() <= 2.0 * e && 3.28e-3 < a * e * e {
        return Err(GammaIncError::Indeterminate { a, x });
    }

    let c = (-y).exp();
    let w = 0.5 * error_fc_scaled(y.sqrt());
    let u = 1.0 / a;
    z = (z + z).sqrt();

    if l < 1.0 {
        z = -z;
    }

    let t;
    if iop < 2 {
        if s.abs() <= 1.0e-3 {
            let c0 = ((((((D0[6] * z + D0[5]) * z + D0[4]) * z + D0[3]) * z + D0[2]) * z + D0[1])
                * z
                + D0[0])
                * z
                - 1.0 / 3.0;

            let c1 = (((((D1[5] * z + D1[4]) * z + D1[3]) * z + D1[2]) * z + D1[1]) * z + D1[0])
                * z
                + D10;

            let c2 = ((((D2[4] * z + D2[3]) * z + D2[2]) * z + D2[1]) * z + D2[0]) * z + D20;

            let c3 = (((D3[3] * z + D3[2]) * z + D3[1]) * z + D3[0]) * z + D30;

            let c4 = (D4[1] * z + D4[0]) * z + D40;
            let c5 = (D5[1] * z + D5[0]) * z + D50;
            let c6 = D6[0] * z + D60;

            t = ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0;
        } else {
            let c0 = ((((((((((((D0[12] * z + D0[11]) * z + D0[10]) * z + D0[9]) * z
                + D0[8])
                * z
                + D0[7])
                * z
                + D0[6])
                * z
                + D0[5])
                * z
                + D0[4])
                * z
                + D0[3])
                * z
                + D0[2])
                * z
                + D0[1])
                * z
                + D0[0])
                * z
                - 1.0 / 3.0;

            let c1 = (((((((((((D1[11] * z + D1[10]) * z + D1[9]) * z + D1[8]) * z + D1[7])
                * z
                + D1[6])
                * z
                + D1[5])
                * z
                + D1[4])
                * z
                + D1[3])
                * z
                + D1[2])
                * z
                + D1[1])
                * z
                + D1[0])
                * z
                + D10;

            let c2 = (((((((((D2[9] * z + D2[8]) * z + D2[7]) * z + D2[6]) * z + D2[5]) * z
                + D2[4])
                * z
                + D2[3])
                * z
                + D2[2])
                * z
                + D2[1])
                * z
                + D2[0])
                * z
                + D20;

            let c3 =
                (((((((D3[7] * z + D3[6]) * z + D3[5]) * z + D3[4]) * z + D3[3]) * z + D3[2]) * z
                    + D3[1])
                    * z
                    + D3[0])
                    * z
                    + D30;

            let c4 = (((((D4[5] * z + D4[4]) * z + D4[3]) * z + D4[2]) * z + D4[1]) * z + D4[0])
                * z
                + D40;

            let c5 = (((D5[3] * z + D5[2]) * z + D5[1]) * z + D5[0]) * z + D50;

            let c6 = (D6[1] * z + D6[0]) * z + D60;

            t = ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0;
        }
    } else if iop == 2 {
        let c0 = (((((D0[5] * z + D0[4]) * z + D0[3]) * z + D0[2]) * z + D0[1]) * z + D0[0]) * z
            - 1.0 / 3.0;

        let c1 = (((D1[3] * z + D1[2]) * z + D1[1]) * z + D1[0]) * z + D10;
        let c2 = D2[0] * z + D20;
        t = (c2 * u + c1) * u + c0;
    } else {
        // 2 < iop.
        t = ((D0[2] * z + D0[1]) * z + D0[0]) * z - 1.0 / 3.0;
    }

    Ok(label_310(l, c, w, t, rta))
}

// Label 310, the common end of labels 270 and 330 (cdflib.f90:11191-11201).
fn label_310(l: f64, c: f64, w: f64, t: f64, rta: f64) -> (f64, f64) {
    if 1.0 <= l {
        let qans = c * (w + RT2PIN * t / rta);
        let ans = 0.5 + (0.5 - qans);
        (ans, qans)
    } else {
        let ans = c * (w - RT2PIN * t / rta);
        let qans = 0.5 + (0.5 - ans);
        (ans, qans)
    }
}

// Temme expansion for L = 1 (label 330, cdflib.f90:11202-11273), ending at
// label 310.
#[allow(clippy::too_many_arguments, clippy::needless_late_init)]
fn temme_expansion_for_l_eq_1(
    a: f64,
    x: f64,
    l: f64,
    mut z: f64,
    y: f64,
    rta: f64,
    e: f64,
    iop: usize,
) -> Result<(f64, f64), GammaIncError> {
    // F90 sets ans = 2 here (cdflib.f90:11207-11210).
    if 3.28e-3 < a * e * e {
        return Err(GammaIncError::Indeterminate { a, x });
    }

    let c = 0.5 + (0.5 - y);
    let w = (0.5 - y.sqrt() * (0.5 + (0.5 - y / 3.0)) / RTPI) / c;
    let u = 1.0 / a;
    z = (z + z).sqrt();

    if l < 1.0 {
        z = -z;
    }

    let t;
    if iop < 2 {
        let c0 = ((((((D0[6] * z + D0[5]) * z + D0[4]) * z + D0[3]) * z + D0[2]) * z + D0[1]) * z
            + D0[0])
            * z
            - 1.0 / 3.0;

        let c1 =
            (((((D1[5] * z + D1[4]) * z + D1[3]) * z + D1[2]) * z + D1[1]) * z + D1[0]) * z + D10;

        let c2 = ((((D2[4] * z + D2[3]) * z + D2[2]) * z + D2[1]) * z + D2[0]) * z + D20;

        let c3 = (((D3[3] * z + D3[2]) * z + D3[1]) * z + D3[0]) * z + D30;

        let c4 = (D4[1] * z + D4[0]) * z + D40;
        let c5 = (D5[1] * z + D5[0]) * z + D50;
        let c6 = D6[0] * z + D60;

        t = ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0;
    } else if iop == 2 {
        let c0 = (D0[1] * z + D0[0]) * z - 1.0 / 3.0;
        let c1 = D1[0] * z + D10;
        t = (D20 * u + c1) * u + c0;
    } else {
        // 2 < iop.
        t = D0[0] * z - 1.0 / 3.0;
    }

    Ok(label_310(l, c, w, t, rta))
}

// Special cases: label 390, a = 0.5 (cdflib.f90:11274-11287).
fn special_cases_390(x: f64) -> (f64, f64) {
    if x < 0.25 {
        let ans = error_f(x.sqrt());
        let qans = 0.5 + (0.5 - ans);
        (ans, qans)
    } else {
        let qans = error_fc(x.sqrt());
        let ans = 0.5 + (0.5 - qans);
        (ans, qans)
    }
}

// Special cases: label 410 (cdflib.f90:11289-11304).
fn special_cases_410(a: f64, x: f64, s: f64, e: f64) -> Result<(f64, f64), GammaIncError> {
    // F90 sets ans = 2 here (cdflib.f90:11291-11294).
    if s.abs() <= 2.0 * e {
        return Err(GammaIncError::Indeterminate { a, x });
    }

    if x <= a {
        Ok((0.0, 1.0))
    } else {
        Ok((1.0, 0.0))
    }
}

// =====================================================================
// gamma_inc_inv: inverse of the incomplete gamma ratio function
// =====================================================================

/// Errors of [`gamma_inc_inv`].
///
/// Each variant corresponds to a negative value of the error flag `ierr`
/// of CDFLIB's `gamma_inc_inv`, except [`AtInfinity`], which corresponds
/// to the result *x* = `huge(x)` that the F90 routine gives with `ierr` =
/// 0 when *q* = 0. Variants for which the F90 routine gives a value for
/// *x* carry it as a field.
///
/// [`gamma_inc_inv`]: crate::special::gamma_inc_inv
/// [`AtInfinity`]: GammaIncInvError::AtInfinity
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaIncInvError {
    /// *a* ≤ 0 (`ierr` = −2).
    #[error("parameter a must be positive, got {0}")]
    ANotPositive(f64),
    /// No solution was obtained: the ratio *q*/*a* is too large (`ierr` =
    /// −3).
    #[error("no solution: q/a is too large")]
    NoSolution,
    /// *p* + *q* ≠ 1 (`ierr` = −4). Rust also returns it when *p* or *q* is
    /// NaN, which the F90 test lets through, or negative, which the F90
    /// requires but does not check.
    #[error("inconsistent inputs: p + q must equal 1")]
    InconsistentPq,
    /// 20 iterations were performed (`ierr` = −6). The F90 documentation
    /// says that this cannot occur if `x0` ≤ 0, but it does, for example
    /// when the solution is subnormal.
    #[error("iteration did not converge in 20 steps; last value: {partial}")]
    NotConverged {
        /// The most recent value obtained for *x*.
        partial: f64,
    },
    /// Iteration failed, and no value is given for *x*. This may occur
    /// when *x* is approximately 0 (`ierr` = −7).
    #[error("iteration failed: intermediate x went non-positive")]
    IterationFailed,
    /// A value for *x* has been obtained, but the routine is not certain
    /// of its accuracy. Iteration cannot be performed in this case. If
    /// `x0` is positive then this can occur when *a* is exceedingly close
    /// to *x* and *a* is extremely large (say *a* ≥ 10²⁰) (`ierr` = −8).
    /// The F90 documentation says that if `x0` ≤ 0 this can occur only
    /// when *p* or *q* is approximately 0, but it also occurs, for
    /// example, when the solution is subnormal.
    #[error("solution obtained but accuracy cannot be certified; value: {value}")]
    UncertainAccuracy {
        /// The value obtained for *x*; [`f64::MAX`] where the F90 routine
        /// sets *x* to `huge(x)`.
        value: f64,
    },
    /// *q* = 0: the solution is +∞, for which the F90 routine gives *x* =
    /// `huge(x)` with `ierr` = 0.
    #[error("inverse is +∞ (q = 0 means P(a, x) = 1 only as x → +∞)")]
    AtInfinity,
}

/// Computes the inverse incomplete Γ ratio function.
///
/// The routine is given positive *a*, and nonnegative *p* and *q* where
/// *p* + *q* = 1. The value *x* is computed with the property that
/// *P*(*a*, *x*) = *p* and *Q*(*a*, *x*) = *q*. Schröder iteration is
/// employed. The routine attempts to compute *x* to 10 significant digits
/// if this is possible for the particular computer arithmetic being used.
/// Algorithm by Alfred Morris.
///
/// `x0` is an optional initial approximation for the solution *x*. If the
/// user does not want to supply an initial approximation, then `x0` should
/// be set to 0, or a negative value.
///
/// Returns (*x*, `ierr`), where `ierr` is 0 if iteration was not used, and
/// otherwise the number of iterations performed.
///
/// # Panics
///
/// Panics on a [`GammaIncInvError`], exactly where [`try_gamma_inc_inv`]
/// returns one; a NaN or negative *p* or *q* is reported as
/// [`InconsistentPq`](GammaIncInvError::InconsistentPq). Use
/// [`try_gamma_inc_inv`] for the fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::gamma_inc_inv;
///
/// // For a = 2.0, p = 0.5, q = 0.5: the median of Γ(2, 1) is ≈ 1.6783.
/// let (x, _iters) = gamma_inc_inv(2.0, -1.0, 0.5, 0.5);
/// assert!((x - 1.6783).abs() < 1e-3);
/// ```
///
/// [`try_gamma_inc_inv`]: crate::special::try_gamma_inc_inv
#[inline]
pub fn gamma_inc_inv(a: f64, x0: f64, p: f64, q: f64) -> (f64, u32) {
    try_gamma_inc_inv(a, x0, p, q)
        .unwrap_or_else(|e| panic!("gamma_inc_inv(a={a}, x0={x0}, p={p}, q={q}): {e}"))
}

/// Fallible form of [`gamma_inc_inv`].
///
/// Returns (*x*, `ierr`) when the solution was obtained: `ierr` is 0 if
/// iteration was not used, and otherwise the number of iterations
/// performed. Each negative `ierr` of the F90 routine, and the result
/// *x* = `huge(x)` at *q* = 0, is returned as a [`GammaIncInvError`].
///
/// # Example
///
/// ```
/// use cdflib::special::{try_gamma_inc_inv, GammaIncInvError};
///
/// let (x, _iters) = try_gamma_inc_inv(2.0, -1.0, 0.5, 0.5).unwrap();
/// assert!((x - 1.6783).abs() < 1e-3);
/// assert_eq!(
///     try_gamma_inc_inv(2.0, -1.0, 1.0, 0.0),
///     Err(GammaIncInvError::AtInfinity),
/// );
/// ```
///
/// [`GammaIncInvError`]: crate::special::GammaIncInvError
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn try_gamma_inc_inv(a: f64, x0: f64, p: f64, q: f64) -> Result<(f64, u32), GammaIncInvError> {
    // Parameters, cdflib.f90:11371-11429. The tables indexed by iop keep
    // the F90 values iop = 1, 2 and are read at [iop - 1]. The table bmin
    // is read only at label 30, and tol only by the Schroder iterations,
    // so they live in label_30, schroder_p and schroder_q.
    const A0: f64 = 3.31125922108741;
    const A1: f64 = 11.6616720288968;
    const A2: f64 = 4.28342155967104;
    const A3: f64 = 0.213623493715853;
    const AMIN: [f64; 2] = [500.0, 100.0];
    const B1: f64 = 6.61053765625462;
    const B2: f64 = 6.40691597760039;
    const B3: f64 = 1.27364489782223;
    const B4: f64 = 0.036117081018842;
    const C: f64 = 0.577215664901533;
    const DMIN: [f64; 2] = [1.0e-6, 1.0e-4];
    const EMIN: [f64; 2] = [2.0e-3, 6.0e-3];
    const EPS0: [f64; 2] = [1.0e-10, 1.0e-8];
    const HALF: f64 = 0.5;
    const LN10: f64 = 2.302585;
    const TWO: f64 = 2.0;

    let e = f64::EPSILON;

    let mut x = 0.0;

    if a <= 0.0 {
        // ierr = -2.
        return Err(GammaIncInvError::ANotPositive(a));
    }

    let t = p + q - 1.0;

    if e < t.abs() {
        // ierr = -4.
        return Err(GammaIncInvError::InconsistentPq);
    }
    // Rust only: a NaN p or q makes t NaN, which passes the test above;
    // the F90 then returns a meaningless x, a negative ierr, or never
    // returns.
    if t.is_nan() {
        return Err(GammaIncInvError::InconsistentPq);
    }
    // Rust only: the F90 requires nonnegative p and q but does not check
    // it; a negative p or q gives a meaningless x, or the F90 never
    // returns.
    if p < 0.0 || q < 0.0 {
        return Err(GammaIncInvError::InconsistentPq);
    }

    let ierr: u32 = 0;

    if p == 0.0 {
        return Ok((x, ierr));
    }

    if q == 0.0 {
        // x = huge(x) with ierr = 0.
        return Err(GammaIncInvError::AtInfinity);
    }

    if a == 1.0 {
        if 0.9 <= q {
            x = -alnrel(-p);
        } else {
            x = -q.ln();
        }
        return Ok((x, ierr));
    }

    let e2 = TWO * e;
    let amax = 0.4e-10 / (e * e);

    let iop = if 1.0e-10 < e { 2 } else { 1 };

    let eps = EPS0[iop - 1];
    let mut xn = x0;

    'l160: {
        if 0.0 < x0 {
            break 'l160;
        }

        // Selection of the initial approximation XN of X when A < 1.
        'l80: {
            if 1.0 < a {
                break 'l80;
            }

            let g = gamma(a + 1.0);
            let qg = q * g;

            if qg == 0.0 {
                // x = huge(x), ierr = -8.
                return Err(GammaIncInvError::UncertainAccuracy { value: f64::MAX });
            }

            let b = qg / a;

            'l40: {
                if 0.6 * a < qg {
                    break 'l40;
                }

                if a < 0.30 && 0.35 <= b {
                    let t = (-(b + C)).exp();
                    let u = t * t.exp();
                    xn = t * u.exp();
                    break 'l160;
                }

                if 0.45 <= b {
                    break 'l40;
                }

                if b == 0.0 {
                    // x = huge(x), ierr = -8.
                    return Err(GammaIncInvError::UncertainAccuracy { value: f64::MAX });
                }

                let y = -b.ln();
                let s = HALF + (HALF - a);
                let z = y.ln();
                let t = y - s * z;

                if 0.15 <= b {
                    xn = y - s * t.ln() - (1.0 + s / (t + 1.0)).ln();
                    // go to 220
                    return schroder_q(a, xn, q, ierr, e2, amax, eps);
                }

                if 0.01 < b {
                    let u = ((t + TWO * (3.0 - a)) * t + (TWO - a) * (3.0 - a))
                        / ((t + (5.0 - a)) * t + TWO);
                    xn = y - s * t.ln() - u.ln();
                    // go to 220
                    return schroder_q(a, xn, q, ierr, e2, amax, eps);
                }

                // Fall through to label 30.
                return label_30(a, b, q, s, y, z, iop, ierr, e2, amax, eps);
            }

            // Label 40.
            if b * q <= 1.0e-8 {
                xn = (-(q / a + C)).exp();
            } else if 0.9 < p {
                xn = ((alnrel(-q) + gamma_ln1(a)) / a).exp();
            } else {
                xn = ((p * g).ln() / a).exp();
            }

            if xn == 0.0 {
                // ierr = -3.
                return Err(GammaIncInvError::NoSolution);
            }

            let t = HALF + (HALF - xn / (a + 1.0));
            xn = xn / t;
            break 'l160;
        }

        // Selection of the initial approximation XN of X when 1 < A.
        // Label 80.
        let mut w = if 0.5 < q { p.ln() } else { q.ln() };

        let t = (-(TWO * w)).sqrt();

        let mut s = t
            - (((A3 * t + A2) * t + A1) * t + A0) / ((((B4 * t + B3) * t + B2) * t + B1) * t + 1.0);

        if 0.5 < q {
            s = -s;
        }

        let rta = a.sqrt();
        let s2 = s * s;

        xn = a + s * rta + (s2 - 1.0) / 3.0 + s * (s2 - 7.0) / (36.0 * rta)
            - ((3.0 * s2 + 7.0) * s2 - 16.0) / (810.0 * a)
            + s * ((9.0 * s2 + 256.0) * s2 - 433.0) / (38880.0 * a * rta);

        xn = xn.max(0.0);

        if AMIN[iop - 1] <= a {
            x = xn;
            let d = HALF + (HALF - x / a);

            if d.abs() <= DMIN[iop - 1] {
                return Ok((x, ierr));
            }
        }

        'l130: {
            if p <= 0.5 {
                break 'l130;
            }

            if xn < 3.0 * a {
                // go to 220
                return schroder_q(a, xn, q, ierr, e2, amax, eps);
            }

            let y = -(w + gamma_log(a));
            let d = TWO.max(a * (a - 1.0));

            if LN10 * d <= y {
                let s = 1.0 - a;
                let z = y.ln();
                // go to 30. The F90 leaves b unset on this path: since
                // 1 < a, label 30 goes to 220 without reading it.
                return label_30(a, f64::NAN, q, s, y, z, iop, ierr, e2, amax, eps);
            }

            let t = a - 1.0;
            xn = y + t * xn.ln() - alnrel(-t / (xn + 1.0));
            xn = y + t * xn.ln() - alnrel(-t / (xn + 1.0));
            // go to 220
            return schroder_q(a, xn, q, ierr, e2, amax, eps);
        }

        // Label 130.
        let ap1 = a + 1.0;

        if 0.70 * ap1 < xn {
            // go to 170
            return schroder_p(a, xn, p, ierr, e2, amax, eps);
        }

        w = w + gamma_log(ap1);

        // The bare F90 literal 0.15 is binary64, as regenerate.sh compiles
        // with -fdefault-real-8.
        if xn <= 0.15 * ap1 {
            let ap2 = a + TWO;
            let ap3 = a + 3.0;
            x = ((w + x) / a).exp();
            x = ((w + x - (1.0 + (x / ap1) * (1.0 + x / ap2)).ln()) / a).exp();
            x = ((w + x - (1.0 + (x / ap1) * (1.0 + x / ap2)).ln()) / a).exp();
            x = ((w + x - (1.0 + (x / ap1) * (1.0 + (x / ap2) * (1.0 + x / ap3))).ln()) / a).exp();
            xn = x;

            if xn <= 1.0e-2 * ap1 {
                if xn <= EMIN[iop - 1] * ap1 {
                    return Ok((x, ierr));
                }
                // go to 170
                return schroder_p(a, xn, p, ierr, e2, amax, eps);
            }
        }

        let mut apn = ap1;
        let mut t = xn / apn;
        let mut sum1 = 1.0 + t;

        // Rust only: when xn is NaN or +inf, t never drops to 1.0e-4 and
        // the F90 loop below never ends; return NaN instead.
        if !xn.is_finite() {
            return Ok((f64::NAN, ierr));
        }

        loop {
            apn = apn + 1.0;
            t = t * (xn / apn);
            sum1 = sum1 + t;

            if t <= 1.0e-4 {
                break;
            }
        }

        let t = w - sum1.ln();
        xn = ((xn + t) / a).exp();
        xn = xn * (1.0 - (a * xn.ln() - xn - t) / (a - xn));
        // go to 170
        return schroder_p(a, xn, p, ierr, e2, amax, eps);
    }

    // Schroder iteration using P.
    // Label 160.
    if 0.5 < p {
        // go to 220
        return schroder_q(a, xn, q, ierr, e2, amax, eps);
    }

    schroder_p(a, xn, p, ierr, e2, amax, eps)
}

// Label 30 of gamma_inc_inv, cdflib.f90:11546-11577. It is entered by
// falling through from the selection for A < 1 (cdflib.f90:11544) and by
// the go to at cdflib.f90:11652 when 1 < A, with different s and z.
#[allow(clippy::too_many_arguments)]
fn label_30(
    a: f64,
    b: f64,
    q: f64,
    s: f64,
    y: f64,
    z: f64,
    iop: usize,
    ierr: u32,
    e2: f64,
    amax: f64,
    eps: f64,
) -> Result<(f64, u32), GammaIncInvError> {
    // bmin, cdflib.f90:11389-11390.
    const BMIN: [f64; 2] = [1.0e-28, 1.0e-13];
    const HALF: f64 = 0.5;
    const TWO: f64 = 2.0;

    let c1 = -(s * z);
    let c2 = -(s * (1.0 + c1));

    let c3 = s * ((HALF * c1 + (TWO - a)) * c1 + (2.5 - 1.5 * a));

    let c4 = -(s
        * (((c1 / 3.0 + (2.5 - 1.5 * a)) * c1 + ((a - 6.0) * a + 7.0)) * c1
            + ((11.0 * a - 46.0) * a + 47.0) / 6.0));

    let c5 = -(s
        * ((((-(c1 / 4.0) + (11.0 * a - 17.0) / 6.0) * c1 + ((-(3.0 * a) + 13.0) * a - 13.0))
            * c1
            + HALF * (((TWO * a - 25.0) * a + 72.0) * a - 61.0))
            * c1
            + (((25.0 * a - 195.0) * a + 477.0) * a - 379.0) / 12.0));

    let xn = ((((c5 / y + c4) / y + c3) / y + c2) / y + c1) + y;

    if 1.0 < a {
        // go to 220
        return schroder_q(a, xn, q, ierr, e2, amax, eps);
    }

    if BMIN[iop - 1] < b {
        // go to 220
        return schroder_q(a, xn, q, ierr, e2, amax, eps);
    }

    let x = xn;
    Ok((x, ierr))
}

// Schroder iteration using P of gamma_inc_inv, labels 170 to 210,
// cdflib.f90:11721-11813. Label 160 (cdflib.f90:11715-11719) is in
// try_gamma_inc_inv.
#[allow(clippy::assign_op_pattern)]
fn schroder_p(
    a: f64,
    mut xn: f64,
    p: f64,
    mut ierr: u32,
    e2: f64,
    amax: f64,
    eps: f64,
) -> Result<(f64, u32), GammaIncInvError> {
    const HALF: f64 = 0.5;
    // tol, cdflib.f90:11428.
    const TOL: f64 = 1.0e-5;

    if p <= 1.0e10 * f64::MIN_POSITIVE {
        // x = xn, ierr = -8.
        return Err(GammaIncInvError::UncertainAccuracy { value: xn });
    }

    let am1 = (a - HALF) - HALF;

    // Label 180.
    loop {
        if amax < a {
            let d = HALF + (HALF - xn / a);
            if d.abs() <= e2 {
                // x = xn, ierr = -8.
                return Err(GammaIncInvError::UncertainAccuracy { value: xn });
            }
        }

        if 20 <= ierr {
            // ierr = -6. The F90 x is the most recent value, equal to xn
            // since label 210.
            return Err(GammaIncInvError::NotConverged { partial: xn });
        }

        ierr = ierr + 1;
        // Rust only: the F90 gamma_inc signals failure with ans = 2 and
        // leaves qn unset; map the failure to ierr = -8 with x = xn.
        let Ok((pn, qn)) = try_gamma_inc(a, xn) else {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        };

        if pn == 0.0 || qn == 0.0 {
            // x = xn, ierr = -8.
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }

        let r = rcomp(a, xn);

        if r == 0.0 {
            // x = xn, ierr = -8.
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }

        let t = (pn - p) / r;
        let w = HALF * (am1 - xn);

        let x;
        let d;

        'l210: {
            'l200: {
                if t.abs() <= 0.1 && (w * t).abs() <= 0.1 {
                    break 'l200;
                }

                x = xn * (1.0 - t);

                if x <= 0.0 {
                    // ierr = -7.
                    return Err(GammaIncInvError::IterationFailed);
                }

                d = t.abs();
                break 'l210;
            }

            // Label 200.
            let h = t * (1.0 + w * t);
            x = xn * (1.0 - h);

            if x <= 0.0 {
                // ierr = -7.
                return Err(GammaIncInvError::IterationFailed);
            }

            if 1.0 <= w.abs() && w.abs() * t * t <= eps {
                return Ok((x, ierr));
            }

            d = h.abs();
        }

        // Label 210.
        xn = x;

        if d <= TOL {
            if d <= eps {
                return Ok((x, ierr));
            }

            if (p - pn).abs() <= TOL * p {
                return Ok((x, ierr));
            }
        }
    }
}

// Schroder iteration using Q of gamma_inc_inv, labels 220 to 260,
// cdflib.f90:11817-11909.
#[allow(clippy::assign_op_pattern)]
fn schroder_q(
    a: f64,
    mut xn: f64,
    q: f64,
    mut ierr: u32,
    e2: f64,
    amax: f64,
    eps: f64,
) -> Result<(f64, u32), GammaIncInvError> {
    const HALF: f64 = 0.5;
    // tol, cdflib.f90:11428.
    const TOL: f64 = 1.0e-5;

    if q <= 1.0e10 * f64::MIN_POSITIVE {
        // x = xn, ierr = -8.
        return Err(GammaIncInvError::UncertainAccuracy { value: xn });
    }

    let am1 = (a - HALF) - HALF;

    // Label 230.
    loop {
        if amax < a {
            let d = HALF + (HALF - xn / a);
            if d.abs() <= e2 {
                // x = xn, ierr = -8.
                return Err(GammaIncInvError::UncertainAccuracy { value: xn });
            }
        }

        if 20 <= ierr {
            // ierr = -6. The F90 x is the most recent value, equal to xn
            // since label 260.
            return Err(GammaIncInvError::NotConverged { partial: xn });
        }

        ierr = ierr + 1;
        // Rust only: the F90 gamma_inc signals failure with ans = 2 and
        // leaves qn unset; map the failure to ierr = -8 with x = xn.
        let Ok((pn, qn)) = try_gamma_inc(a, xn) else {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        };

        if pn == 0.0 || qn == 0.0 {
            // x = xn, ierr = -8.
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }

        let r = rcomp(a, xn);

        if r == 0.0 {
            // x = xn, ierr = -8.
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }

        let t = (q - qn) / r;
        let w = HALF * (am1 - xn);

        let x;
        let d;

        'l260: {
            'l250: {
                // The bare F90 literals 0.1 are binary64, as regenerate.sh
                // compiles with -fdefault-real-8.
                if t.abs() <= 0.1 && (w * t).abs() <= 0.1 {
                    break 'l250;
                }

                x = xn * (1.0 - t);

                if x <= 0.0 {
                    // ierr = -7.
                    return Err(GammaIncInvError::IterationFailed);
                }

                d = t.abs();
                break 'l260;
            }

            // Label 250.
            let h = t * (1.0 + w * t);
            x = xn * (1.0 - h);

            if x <= 0.0 {
                // ierr = -7.
                return Err(GammaIncInvError::IterationFailed);
            }

            if 1.0 <= w.abs() && w.abs() * t * t <= eps {
                return Ok((x, ierr));
            }

            d = h.abs();
        }

        // Label 260.
        xn = x;

        if TOL < d {
            continue;
        }

        if d <= eps {
            return Ok((x, ierr));
        }

        if (q - qn).abs() <= TOL * q {
            return Ok((x, ierr));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================================ dexpm1

    #[test]
    fn dexpm1_at_zero() {
        assert_eq!(dexpm1(0.0), 0.0);
    }

    #[test]
    fn dexpm1_small_argument_matches_rational() {
        // |x| ≤ 0.15: rational approximation. Should match exp(x) − 1
        // closely (better than naive subtraction near zero).
        for &x in &[-0.1_f64, -0.05, 0.0001, 0.05, 0.1, 0.15] {
            let got = dexpm1(x);
            let ref_val = x.exp() - 1.0;
            assert!(
                (got - ref_val).abs() < 1e-14,
                "x={x}: dexpm1={got}, ref={ref_val}"
            );
        }
    }

    #[test]
    fn dexpm1_large_argument_matches_exp_minus_one() {
        for &x in &[-3.0_f64, -1.0, 0.5, 1.0, 5.0] {
            let got = dexpm1(x);
            let ref_val = x.exp() - 1.0;
            let scale = ref_val.abs().max(1.0);
            assert!(
                (got - ref_val).abs() / scale < 1e-14,
                "x={x}: dexpm1={got}, ref={ref_val}"
            );
        }
    }

    // Exact 1e-15 agreement requires bit-identical libm exp; miri's
    // soft-float shim drifts by 1 ULP at x = 2. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn dexpm1_matches_rexp_in_overlap() {
        // Both routines exist; for the same x they should agree.
        for &x in &[-2.0_f64, -0.5, -0.01, 0.0, 0.01, 0.5, 2.0] {
            let a = dexpm1(x);
            let b = rexp(x);
            assert!((a - b).abs() < 1e-15, "x={x}: dexpm1={a}, rexp={b}");
        }
    }

    // ============================================================ dstrem

    #[test]
    fn dstrem_large_z_matches_bernoulli_lead() {
        // For large z, dstrem(z) ≈ 1/(12 z) − 1/(360 z³) + … .
        // At z = 100 the next term is about 3.3·10⁻⁶ of the leading one,
        // so the leading-term match is about 3.3·10⁻⁶.
        let r = dstrem(100.0);
        let lead = 1.0 / 1200.0;
        assert!((r - lead).abs() / lead < 1e-4, "r = {r}, leading = {lead}");
    }

    #[test]
    fn dstrem_small_z_matches_explicit_difference() {
        // For z ≤ 6, dstrem uses gamma_log(z) − Sterling(z) directly.
        // At z = 5: lnΓ(5) = ln 24 = 3.178053830347946...,
        // Sterling(5) ≈ ½ ln(2π) + 4.5·ln 5 − 5 = 3.161409...,
        // so dstrem(5) ≈ 0.016645...
        let r = dstrem(5.0);
        let sterl = 0.91893853320467274178 + 4.5 * 5.0_f64.ln() - 5.0;
        let expected = gamma_log(5.0) - sterl;
        assert!((r - expected).abs() < 1e-14, "r = {r}");
    }

    #[test]
    fn dstrem_continuous_across_z_eq_6() {
        // The two branches (z ≤ 6 and z > 6) should agree closely at z = 6.
        let just_below = dstrem(6.0);
        let just_above = dstrem(6.0 + 1.0e-9);
        assert!((just_below - just_above).abs() < 1e-7);
    }

    #[test]
    #[should_panic(expected = "argument z must be positive")]
    fn dstrem_panics_on_nonpositive() {
        let _ = dstrem(0.0);
    }

    // ============================================================ gamma_inc_inv

    /// Helper: round-trip residual at the returned x.
    fn round_trip_residual(a: f64, p: f64, q: f64) -> (f64, f64) {
        let (x, _iters) = gamma_inc_inv(a, -1.0, p, q);
        let (pn, qn) = gamma_inc(a, x);
        (pn - p, qn - q)
    }

    #[test]
    fn gamma_inc_inv_rejects_invalid_a() {
        assert_eq!(
            try_gamma_inc_inv(0.0, -1.0, 0.5, 0.5),
            Err(GammaIncInvError::ANotPositive(0.0))
        );
        assert_eq!(
            try_gamma_inc_inv(-1.0, -1.0, 0.5, 0.5),
            Err(GammaIncInvError::ANotPositive(-1.0))
        );
    }

    #[test]
    fn gamma_inc_inv_rejects_p_plus_q_not_one() {
        assert_eq!(
            try_gamma_inc_inv(2.0, -1.0, 0.5, 0.6),
            Err(GammaIncInvError::InconsistentPq)
        );
    }

    #[test]
    fn gamma_inc_inv_trivial_endpoints() {
        // p = 0 ⇒ x = 0 (P(a, 0) = 0 for all a > 0).
        assert_eq!(try_gamma_inc_inv(2.0, -1.0, 0.0, 1.0), Ok((0.0, 0)));
        // q = 0 means P(a, x) = 1, which only happens as x → +∞:
        // reported as the AtInfinity variant.
        assert_eq!(
            try_gamma_inc_inv(2.0, -1.0, 1.0, 0.0),
            Err(GammaIncInvError::AtInfinity)
        );
    }

    #[test]
    #[should_panic(expected = "inverse is +∞")]
    fn gamma_inc_inv_at_infinity_panics() {
        let _ = gamma_inc_inv(2.0, -1.0, 1.0, 0.0);
    }

    #[test]
    fn gamma_inc_inv_a_equals_one_closed_form() {
        // a = 1: incomplete gamma collapses to the exponential.
        // Q(1, x) = exp(-x) ⇒ x = -ln(q). The F90 also uses
        // -alnrel(-p) when q ≥ 0.9 (== p ≤ 0.1) for tail accuracy.
        for &(p, q) in &[(0.5, 0.5), (0.9, 0.1), (0.05, 0.95), (0.001, 0.999)] {
            let (x, iters) = gamma_inc_inv(1.0, -1.0, p, q);
            assert!(
                (-q.ln() - x).abs() / x.abs().max(1.0) < 1e-13,
                "p={p}: x={x}"
            );
            // The a = 1 path is closed-form (F90 ierr = 0): no iteration.
            assert_eq!(iters, 0, "p={p}: unexpected Schröder iteration");
        }
    }

    #[test]
    fn gamma_inc_inv_round_trip_small_a() {
        // a < 1 branch.
        for &a in &[0.05_f64, 0.2, 0.5, 0.95] {
            for &p in &[0.1_f64, 0.25, 0.5, 0.75, 0.9] {
                let q = 1.0 - p;
                let (dp, dq) = round_trip_residual(a, p, q);
                assert!(
                    dp.abs() < 1e-7 && dq.abs() < 1e-7,
                    "a={a}, p={p}: |p_n − p| = {}, |q_n − q| = {}",
                    dp.abs(),
                    dq.abs(),
                );
            }
        }
    }

    #[test]
    fn gamma_inc_inv_round_trip_large_a() {
        // a > 1 branch.
        for &a in &[1.5_f64, 5.0, 50.0, 500.0] {
            for &p in &[0.1_f64, 0.25, 0.5, 0.75, 0.9] {
                let q = 1.0 - p;
                let (dp, dq) = round_trip_residual(a, p, q);
                assert!(
                    dp.abs() < 1e-7 && dq.abs() < 1e-7,
                    "a={a}, p={p}: |p_n − p| = {}, |q_n − q| = {}",
                    dp.abs(),
                    dq.abs(),
                );
            }
        }
    }

    #[test]
    fn gamma_inc_inv_deep_tails() {
        // Verify tail behavior: p near 0 (small x) and p near 1 (large x).
        let (x, _) = gamma_inc_inv(3.0, -1.0, 1.0e-6, 1.0 - 1.0e-6);
        let (pn, _) = gamma_inc(3.0, x);
        assert!((pn - 1.0e-6).abs() / 1.0e-6 < 1e-4);

        let (x, _) = gamma_inc_inv(3.0, -1.0, 1.0 - 1.0e-6, 1.0e-6);
        let (_, qn) = gamma_inc(3.0, x);
        assert!((qn - 1.0e-6).abs() / 1.0e-6 < 1e-4);
    }

    #[test]
    fn gamma_inc_inv_with_caller_supplied_x0() {
        // x0 > 0 mode: caller supplies an initial approximation. The
        // routine should still converge to the same answer (within
        // tolerance) as the x0 ≤ 0 mode.
        let (x_auto, _) = gamma_inc_inv(5.0, -1.0, 0.7, 0.3);
        let (x_seeded, _) = gamma_inc_inv(5.0, x_auto * 1.1, 0.7, 0.3);
        assert!((x_seeded - x_auto).abs() / x_auto < 1e-7);
    }

    #[test]
    fn gamma_inc_inv_sweeps_all_regimes() {
        // Fine grid exercising the selection of the initial approximation:
        //   a <= 1 with various q levels (different b = qg/a ranges),
        //     reaching label 30 when b <= 0.01
        //   1 < a (label 80)
        //   amin(1) = 500 <= a (cdflib.f90:11627-11636; at these p the
        //     test |d| <= dmin(1) fails and the iteration follows)
        let grid_a = [
            0.01_f64, 0.05, 0.1, 0.25, 0.5, 0.95, // a <= 1
            1.5, 2.0, 5.0, 25.0, 99.0, 150.0, // 1 < a < amin(1)
            600.0, 2000.0, // amin(1) <= a
        ];
        let grid_p = [
            1.0e-9_f64,
            1.0e-5,
            1.0e-3,
            0.01,
            0.05,
            0.1,
            0.3,
            0.5,
            0.7,
            0.9,
            0.95,
            0.99,
            0.999,
            0.99999,
            1.0 - 1.0e-9,
        ];
        for &a in &grid_a {
            for &p in &grid_p {
                let q = 1.0 - p;
                let result = try_gamma_inc_inv(a, -1.0, p, q);
                // The F90 documents three "give-up" outcomes that are
                // part of its contract, not port regressions: no
                // solution (NoSolution), iterate went non-positive
                // (IterationFailed), and accuracy cannot be certified
                // (UncertainAccuracy).
                let (x, _iters) = match result {
                    Ok(pair) => pair,
                    Err(GammaIncInvError::NoSolution)
                    | Err(GammaIncInvError::IterationFailed)
                    | Err(GammaIncInvError::UncertainAccuracy { .. }) => continue,
                    Err(e) => panic!("a={a}, p={p}: unexpected error {e:?}"),
                };
                assert!(x.is_finite() && x > 0.0, "a={a}, p={p}: x = {x}");
                let (pn, qn) = gamma_inc(a, x);
                let dp = (pn - p).abs() / p.max(1e-300);
                let dq = (qn - q).abs() / q.max(1e-300);
                // Either tail should match to ~1e-5; CDFLIB's stated goal
                // is 10 significant digits when possible, but Schröder
                // iteration's tolerance constant tol = 1e-5 is the
                // floor.
                assert!(dp.min(dq) < 1e-4, "a={a}, p={p}, x={x}: dp={dp}, dq={dq}",);
            }
        }
    }

    #[test]
    fn gamma_inc_inv_caller_supplied_x0_p_le_half() {
        // x0 > 0 mode with p ≤ 0.5 routes through schroder_p; the auto
        // mode with p ≤ 0.5 also routes through it for small a.
        let (x, _) = gamma_inc_inv(3.0, 1.5, 0.3, 0.7);
        let (pn, _) = gamma_inc(3.0, x);
        assert!((pn - 0.3).abs() < 1e-7);
    }

    #[test]
    fn gamma_inc_inv_caller_supplied_x0_p_gt_half() {
        // x0 > 0 with p > 0.5 routes through schroder_q.
        let (x, _) = gamma_inc_inv(3.0, 5.0, 0.8, 0.2);
        let (_, qn) = gamma_inc(3.0, x);
        assert!((qn - 0.2).abs() < 1e-7);
    }

    // ============================================================ known values

    #[test]
    fn gamma_log_at_small_integer_arguments() {
        // ln Γ(1) = 0, ln Γ(2) = 0, ln Γ(3) = ln 2, ln Γ(4) = ln 6
        assert!((gamma_log(1.0)).abs() < 1e-14);
        assert!((gamma_log(2.0)).abs() < 1e-14);
        assert!((gamma_log(3.0) - 2.0_f64.ln()).abs() < 1e-14);
        assert!((gamma_log(4.0) - 6.0_f64.ln()).abs() < 1e-14);
        assert!((gamma_log(10.0) - 362880.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn gamma_at_small_integers() {
        assert!((gamma(1.0) - 1.0).abs() < 1e-14);
        assert!((gamma(2.0) - 1.0).abs() < 1e-14);
        assert!((gamma(3.0) - 2.0).abs() < 1e-14);
        assert!((gamma(4.0) - 6.0).abs() < 1e-14);
        assert!((gamma(5.0) - 24.0).abs() < 1e-13);
    }

    #[test]
    fn gamma_inc_basic_identities() {
        // P(a, 0) = 0, Q(a, 0) = 1 for a > 0.
        let (p, q) = gamma_inc(2.5, 0.0);
        assert_eq!(p, 0.0);
        assert_eq!(q, 1.0);

        // P(a, x) + Q(a, x) = 1.
        for &a in &[0.5, 1.0, 2.5, 7.0, 25.0] {
            for &x in &[0.1, 1.0, 5.0, 20.0] {
                let (p, q) = gamma_inc(a, x);
                assert!((p + q - 1.0).abs() < 1e-12, "a={a}, x={x}: p+q = {}", p + q);
            }
        }
    }

    #[test]
    fn psi_at_known_points() {
        // ψ(1) = -γ (Euler–Mascheroni)
        let γ = 0.5772156649015328606;
        assert!((psi(1.0) + γ).abs() < 1e-9, "psi(1) = {}", psi(1.0));
        // ψ(2) = 1 - γ
        assert!((psi(2.0) - (1.0 - γ)).abs() < 1e-9);
        // ψ(0.5) = -γ - 2 ln 2
        let expected = -γ - 2.0 * 2.0_f64.ln();
        assert!(
            (psi(0.5) - expected).abs() < 1e-9,
            "psi(0.5) = {}",
            psi(0.5)
        );
    }

    #[test]
    fn gamma_inc_at_a_half_uses_erf() {
        // P(1/2, x) = erf(sqrt(x)).
        for &x in &[0.1, 0.5, 1.0, 4.0, 9.0] {
            let (p, _q) = gamma_inc(0.5, x);
            let expected = error_f(x.sqrt());
            assert!(
                (p - expected).abs() < 1e-13,
                "x={x}: P = {p}, erf(√x) = {expected}"
            );
        }
    }

    // ===== Validation and error sentinels =====

    #[test]
    fn gamma_inc_rejects_negative_a() {
        assert!(matches!(
            try_gamma_inc(-1.0, 1.0),
            Err(GammaIncError::ANegative(_))
        ));
    }

    #[test]
    fn gamma_inc_rejects_negative_x() {
        assert!(matches!(
            try_gamma_inc(1.0, -1.0),
            Err(GammaIncError::XNegative(_))
        ));
    }

    #[test]
    fn gamma_inc_rejects_both_zero() {
        assert!(matches!(
            try_gamma_inc(0.0, 0.0),
            Err(GammaIncError::BothZero)
        ));
    }

    #[test]
    fn nan_arguments_reach_the_f90_error_checks() {
        // The F90 checks that come first still apply: gamma_inc sets
        // ans = 2 (cdflib.f90:10681-10684) and gamma_inc_inv sets ierr = -2
        // (cdflib.f90:11442-11445). The F90 test e < abs(t) is false for a
        // NaN t, so the F90 goes on to return x = 0 when p = 0
        // (cdflib.f90:11456-11458) and x = huge(x) when q = 0
        // (cdflib.f90:11460-11463); Rust reports InconsistentPq instead.
        let nan = f64::NAN;
        assert_eq!(
            try_gamma_inc(-1.0, nan),
            Err(GammaIncError::ANegative(-1.0))
        );
        assert_eq!(
            try_gamma_inc(nan, -1.0),
            Err(GammaIncError::XNegative(-1.0))
        );
        assert_eq!(
            try_gamma_inc_inv(0.0, -1.0, nan, 0.5),
            Err(GammaIncInvError::ANotPositive(0.0))
        );
        assert_eq!(
            try_gamma_inc_inv(2.0, -1.0, 0.0, nan),
            Err(GammaIncInvError::InconsistentPq)
        );
        assert_eq!(
            try_gamma_inc_inv(2.0, -1.0, nan, 0.0),
            Err(GammaIncInvError::InconsistentPq)
        );
        assert_eq!(
            try_gamma_inc_inv(2.0, -1.0, nan, 0.5),
            Err(GammaIncInvError::InconsistentPq)
        );
        // The panicking forms agree with the fallible ones.
        let (p, q) = gamma_inc(2.0, nan);
        assert!(p.is_nan() && q.is_nan());
    }

    #[test]
    #[should_panic(expected = "parameter a must be non-negative")]
    fn gamma_inc_negative_a_with_nan_x_panics() {
        let _ = gamma_inc(-1.0, f64::NAN);
    }

    #[test]
    #[should_panic(expected = "parameter a must be positive")]
    fn gamma_inc_inv_zero_a_with_nan_p_panics() {
        let _ = gamma_inc_inv(0.0, -1.0, f64::NAN, 0.5);
    }

    #[test]
    fn gamma_inc_a_zero_x_positive() {
        // a = 0, x > 0: P(0, x) = 1 (the limit).
        let (p, q) = gamma_inc(0.0, 1.0);
        assert_eq!(p, 1.0);
        assert_eq!(q, 0.0);
    }

    // ===== Regime switching points =====

    #[test]
    fn gamma_inc_regime_a_lt_1_x_lt_1() {
        // a < 1 and x < 1.1: Taylor series for P(A,X)/X^A.
        let (p, q) = gamma_inc(0.3, 0.2);
        assert!((p + q - 1.0).abs() < 1e-14);
        assert!(p > 0.0 && p < 1.0);
    }

    #[test]
    fn gamma_inc_regime_a_lt_1_x_ge_1() {
        // a < 1 and 1.1 <= x: continued fraction expansion.
        let (p, q) = gamma_inc(0.3, 5.0);
        assert!((p + q - 1.0).abs() < 1e-14);
        assert!(p > 0.99); // upper-tail saturation
    }

    #[test]
    fn gamma_inc_regime_a_eq_1() {
        // Boundary a == 1: exponential CDF.
        let (p, q) = gamma_inc(1.0, 2.0);
        let expected_p = 1.0 - (-2.0_f64).exp();
        assert!((p - expected_p).abs() < 1e-14);
        assert!((q - (-2.0_f64).exp()).abs() < 1e-14);
    }

    #[test]
    fn gamma_inc_regime_a_large_x_near_a() {
        // big(1) = 20 <= a and x = a: Temme expansion for L = 1.
        let (p, q) = gamma_inc(100.0, 100.0);
        assert!((p + q - 1.0).abs() < 1e-12);
        // Gamma(a, 1) at x=a: just above the median (which is at
        // ≈ a-1/3 for large a), so cdf is slightly above 0.5.
        assert!(p > 0.5 && p < 0.55, "p={p}");
    }

    #[test]
    fn gamma_inc_regime_a_very_large() {
        // Temme expansion for L = 1 at a very large a.
        let (p, q) = gamma_inc(1e6, 1e6);
        assert!((p + q - 1.0).abs() < 1e-12);
        assert!(p.is_finite() && q.is_finite());
        assert!(p > 0.499 && p < 0.501);
    }

    #[test]
    fn gamma_inc_half_integer_a_uses_finite_sum() {
        // a half-integer ≥ 1: finite-sum regime.
        for &a in &[1.5_f64, 2.5, 3.5, 10.5] {
            let (p, q) = gamma_inc(a, a);
            assert!((p + q - 1.0).abs() < 1e-12, "a={a}");
        }
    }

    #[test]
    fn gamma_inc_a_x_very_unbalanced() {
        // x >> a: deep right tail.
        let (p, q) = gamma_inc(2.0, 50.0);
        assert!(p > 0.9999);
        assert!(q > 0.0 && q < 1e-15);
        // x << a: deep left tail.
        let (p, q) = gamma_inc(50.0, 2.0);
        assert!(q > 0.9999);
        assert!(p > 0.0 && p < 1e-15);
    }

    #[test]
    fn gamma_at_half_integer() {
        // Γ(1/2) = √π.
        assert!((gamma(0.5) - std::f64::consts::PI.sqrt()).abs() < 1e-13);
        // Γ(3/2) = √π/2.
        assert!((gamma(1.5) - std::f64::consts::PI.sqrt() / 2.0).abs() < 1e-13);
    }

    #[test]
    fn try_gamma_overflow_and_underflow_are_err() {
        // gamma_user returns its sentinel for 1000 <= |a|; Rust reports
        // GammaDomainError::Overflow for a >= 1000 and Underflow for
        // a <= -1000.
        assert_eq!(try_gamma(1001.0), Err(GammaDomainError::Overflow(1001.0)));
        assert_eq!(try_gamma(1e10), Err(GammaDomainError::Overflow(1e10)));
        assert_eq!(
            try_gamma(-1000.5),
            Err(GammaDomainError::Underflow(-1000.5))
        );
        assert_eq!(
            try_gamma(-1000.0),
            Err(GammaDomainError::Underflow(-1000.0))
        );
    }

    #[test]
    #[should_panic(expected = "gamma(1001): Γ(1001) overflows f64")]
    fn gamma_overflow_panics() {
        let _ = gamma(1001.0);
    }

    #[test]
    fn gamma_at_negative_non_integer() {
        // Γ(-0.5) = -2√π via reflection.
        let expected = -2.0 * std::f64::consts::PI.sqrt();
        let got = gamma(-0.5);
        assert!(
            (got - expected).abs() < 1e-12,
            "got = {got}, expected = {expected}"
        );
        // Γ(-1.5) = 4√π/3.
        let expected = 4.0 / 3.0 * std::f64::consts::PI.sqrt();
        let got = gamma(-1.5);
        assert!(
            (got - expected).abs() < 1e-12,
            "got = {got}, expected = {expected}"
        );
    }

    #[test]
    fn gamma_at_negative_mid_range() {
        // Γ(-3.5) = 16√π/105 ≈ 0.2701. Check it against the reflection
        // formula Γ(z)Γ(1 - z) = π/sin(πz), which gives
        // Γ(-3.5) = π / (sin(-3.5π) Γ(4.5)).
        let g = gamma(-3.5);
        assert!(g.is_finite());
        let g_45 = gamma(4.5);
        let expected = std::f64::consts::PI / ((-3.5_f64 * std::f64::consts::PI).sin() * g_45);
        assert!((g - expected).abs() / expected.abs() < 1e-10);
    }

    #[test]
    fn gamma_at_negative_large_magnitude() {
        // Γ(-20.5): asymptotic-reflection branch for |a| ≥ 15.
        let g = gamma(-20.5);
        // Should be finite and tiny (~1e-19).
        assert!(g.is_finite() && g.abs() < 1e-15);
    }

    #[test]
    fn try_gamma_at_negative_integer_is_pole() {
        // Γ has a pole at every non-positive integer.
        assert_eq!(try_gamma(-3.0), Err(GammaDomainError::Pole(-3.0)));
        assert_eq!(try_gamma(-10.0), Err(GammaDomainError::Pole(-10.0)));
    }

    #[test]
    fn try_gamma_at_zero_is_pole() {
        // Γ has a pole at 0 (the mathematical limit is +∞).
        assert_eq!(try_gamma(0.0), Err(GammaDomainError::Pole(0.0)));
    }

    #[test]
    #[should_panic(expected = "gamma(0): Γ has a pole at 0")]
    fn gamma_at_zero_panics() {
        let _ = gamma(0.0);
    }

    #[test]
    fn gamma_at_large_positive() {
        // 15 <= a < 1000 uses the modified asymptotic sum.
        let ln_gamma_50 = gamma_log(50.0);
        let g_50 = gamma(50.0);
        // log(g_50) should match ln_gamma_50 to high precision.
        assert!((g_50.ln() - ln_gamma_50).abs() < 1e-9);
        // For very large a (but still < 1000), result should be huge but finite.
        let g_100 = gamma(100.0);
        assert!(g_100.is_finite() && g_100 > 1e150);
    }

    #[test]
    fn gamma_log_at_every_branch() {
        // gamma_log has branches for a <= 0.8, 0.8 < a <= 2.25,
        // 2.25 < a < 10 and 10 <= a; evaluate each and the boundaries.
        for &a in &[0.5_f64, 0.8, 1.0, 1.5, 2.0, 2.25, 5.0, 8.0, 10.0, 50.0] {
            let h = gamma_log(a);
            assert!(h.is_finite(), "a={a}");
        }
    }

    #[test]
    fn psi_at_a_lt_05_and_large() {
        // psi has different branches for x < 0.5, x in [0.5..3] and x > 3.
        assert!(psi(0.1).is_finite());
        assert!(psi(0.3).is_finite());
        assert!(psi(2.5).is_finite());
        assert!(psi(50.0).is_finite());
        assert!(psi(1000.0).is_finite());
    }
}
