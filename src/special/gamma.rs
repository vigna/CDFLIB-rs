//! Γ function family: ln Γ, Γ, regularized incomplete Γ ratios
//! *P*(*a*, *x*) and *Q*(*a*, *x*).
//!
//! This implementations dispatches across five computational regimes:
//!
//! - **Power series** for small *a*, *x*.
//!
//! - **Continued fraction** in the body when *x* > *a* is moderate.
//!
//! - **Tricomi–Temme asymptotic expansion** for large *a*.
//!
//! - **Finite sum** when *a* ≥ 1 and 2*a* is an integer.
//!
//! - **Special cases** for *a* = 1/2 (reduces to [`error_f`] / [`error_fc`])
//!   and for *a*·*x* = 0.
//!
//! [`gamma_inc`]: crate::special::gamma_inc
//! [`error_f`]: crate::special::error_f
//! [`error_fc`]: crate::special::error_fc

#![allow(clippy::approx_constant, clippy::excessive_precision)]

use super::erf::{error_f, error_fc, error_fc_scaled};
use super::eval_pol;

// =====================================================================
// Low-level helpers
// =====================================================================

/// Returns ln(1 + *a*). Rational approximation for |*a*| ≤ 0.375, otherwise fall
/// back to `(1 + a).ln()`.
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
        let num = ((P3 * t2 + P2) * t2 + P1) * t2 + 1.0;
        let den = ((Q3 * t2 + Q2) * t2 + Q1) * t2 + 1.0;
        2.0 * t * (num / den)
    } else {
        (1.0 + a).ln()
    }
}

/// Returns exp(*x*) − 1. Rational approximation for |*x*| ≤ 0.15 to preserve
/// precision near zero where naive `exp(x) - 1.0` cancels.
#[inline]
pub fn rexp(x: f64) -> f64 {
    const P1: f64 = 0.914041914819518e-9;
    const P2: f64 = 0.238082361044469e-1;
    const Q1: f64 = -0.499999999085958;
    const Q2: f64 = 0.107141568980644;
    const Q3: f64 = -0.119041179760821e-1;
    const Q4: f64 = 0.595130811860248e-3;
    if x.abs() <= 0.15 {
        let num = (P2 * x + P1) * x + 1.0;
        let den = (((Q4 * x + Q3) * x + Q2) * x + Q1) * x + 1.0;
        x * (num / den)
    } else {
        let w = x.exp();
        if x <= 0.0 {
            (w - 0.5) - 0.5
        } else {
            w * (0.5 + (0.5 - 1.0 / w))
        }
    }
}

/// Returns exp(*x*) − 1. Double-precision sibling of [`rexp`] (same algorithm:
/// rational approximation for |*x*| ≤ 0.15, fall back to `exp(x)` outside).
///
/// [`rexp`]: crate::special::internal::rexp
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

/// Returns *x* − 1 − ln(*x*). A precision-preserving form valid for *x* near 1.
#[inline]
pub fn rlog(x: f64) -> f64 {
    const A: f64 = 0.566749439387324e-1;
    const B: f64 = 0.456512608815524e-1;
    const P0: f64 = 0.333333333333333;
    const P1: f64 = -0.224696413112536;
    const P2: f64 = 0.620886815375787e-2;
    const Q1: f64 = -0.127408923933623e1;
    const Q2: f64 = 0.354508718369557;

    // Mirror F90 rlog (cdflib.f90:13660-13700) branch boundaries exactly:
    // strict less-than at 0.61, 0.82, 1.18, 1.57. Boundary points
    // 0.82, 1.18, 1.57 fall into the next branch (not the previous one).
    if x < 0.61 {
        let r = x - 0.5 - 0.5;
        return r - x.ln();
    }
    if x < 1.57 {
        let (u, w1) = if x < 0.82 {
            let u = (x - 0.7) / 0.7;
            (u, A - u * 0.3)
        } else if x < 1.18 {
            let u = x - 0.5 - 0.5;
            (u, 0.0)
        } else {
            let u = 0.75 * x - 1.0;
            (u, B + u / 3.0)
        };
        let r = u / (u + 2.0);
        let t = r * r;
        let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
        return 2.0 * t * (1.0 / (1.0 - r) - r * w) + w1;
    }
    // 1.57 <= x
    let r = x - 0.5 - 0.5;
    r - x.ln()
}

/// Returns *x* − ln(1 + *x*). Sibling of `rlog` for shifted input.
#[inline]
pub fn rlog1(x: f64) -> f64 {
    const A: f64 = 0.566749439387324e-1;
    const B: f64 = 0.456512608815524e-1;
    const P0: f64 = 0.333333333333333;
    const P1: f64 = -0.224696413112536;
    const P2: f64 = 0.620886815375787e-2;
    const Q1: f64 = -0.127408923933623e1;
    const Q2: f64 = 0.354508718369557;

    if !(-0.39..=0.57).contains(&x) {
        let w = x + 0.5 + 0.5;
        return x - w.ln();
    }
    let (h, w1) = if x < -0.18 {
        let h = (x + 0.3) / 0.7;
        (h, A - h * 0.3)
    } else if x > 0.18 {
        let h = 0.75 * x - 0.25;
        (h, B + h / 3.0)
    } else {
        (x, 0.0)
    };
    let r = h / (h + 2.0);
    let t = r * r;
    let w = ((P2 * t + P1) * t + P0) / ((Q2 * t + Q1) * t + 1.0);
    2.0 * t * (1.0 / (1.0 - r) - r * w) + w1
}

/// Returns 1/Γ(1 + *a*) − 1 for −0.5 ≤ *a* ≤ 1.5. Used to evaluate
/// Γ(1 + *a*) accurately when *a* is small (without computing the
/// near-1 value via subtraction).
#[inline]
pub fn gam1(a: f64) -> f64 {
    const S1: f64 = 0.273076135303957;
    const S2: f64 = 0.559398236957378e-1;
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
        1.0,
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

    let mut t = a;
    let d = a - 0.5;
    if d > 0.0 {
        t = d - 0.5;
    }
    if t == 0.0 {
        return 0.0;
    }
    if t > 0.0 {
        let top =
            ((((((P[6] * t + P[5]) * t + P[4]) * t + P[3]) * t + P[2]) * t + P[1]) * t) + P[0];
        let bot = (((Q[4] * t + Q[3]) * t + Q[2]) * t + Q[1]) * t + 1.0;
        let w = top / bot;
        if d <= 0.0 {
            a * w
        } else {
            (t / a) * ((w - 0.5) - 0.5)
        }
    } else {
        // t < 0
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
    }
}

// =====================================================================
// Γ(a): the Γ function itself
// =====================================================================

/// Errors of [`try_gamma`].
///
/// In CDFLIB these were signalled by returning the sentinel `0.0` (the
/// `gamma_user` body initializes the result to 0 and falls through on
/// every overflow/pole path). Here each failure mode is lifted into a
/// typed Rust error.
///
/// [`try_gamma`]: crate::special::try_gamma
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaDomainError {
    /// Argument is zero or a negative integer; Γ has a pole there.
    #[error("Γ has a pole at {0}")]
    Pole(f64),
    /// Result would overflow f64 (*a* ≥ 1000, or the asymptotic
    /// `exp(g)` would overflow).
    #[error("Γ({0}) overflows f64")]
    Overflow(f64),
    /// Argument is too negative (*a* ≤ −1000). F90 `gamma_user`
    /// (cdflib.f90:10162-10164) returns its sentinel 0 here; the true
    /// value of Γ underflows.
    #[error("Γ({0}) underflows f64")]
    Underflow(f64),
}

/// Returns Γ(*a*), the Γ function.
///
/// # Panics
///
/// Panics on a [`GammaDomainError`] (pole at a non-positive integer, or
/// overflow when |*a*| ≥ 1000). Use [`try_gamma`] for the fallible form.
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
    if a.is_nan() {
        return f64::NAN;
    }
    try_gamma(a).unwrap_or_else(|e| panic!("gamma({a}): {e}"))
}

/// Fallible form of [`gamma`]: returns [`GammaDomainError`] when the argument
/// lands on a pole or the result would overflow f64.
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
pub fn try_gamma(a: f64) -> Result<f64, GammaDomainError> {
    const D: f64 = 0.41893853320467274178;
    const PI: f64 = 3.1415926535898;
    const P_COEF: [f64; 7] = [
        0.539637273585445e-3,
        0.261939260042690e-2,
        0.204493667594920e-1,
        0.730981088720487e-1,
        0.279648642639792,
        0.553413866010467,
        1.0,
    ];
    const Q_COEF: [f64; 7] = [
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
    /// Largest exp arg for IEEE 754 binary64; corresponds to CDFLIB's
    /// `exparg(0) = 0.99999 * (1024 * lnb)`, with `lnb = 0.69314718055995`
    /// (the F90 truncation of ln 2). The `0.99999` safety factor is baked
    /// in here so that `0.99999 * POS_EXPARG` reproduces F90's call-site
    /// expression `0.99999 * exparg(0)` exactly.
    const POS_EXPARG: f64 = 709.775_615_066_259_888_4;

    let mut x = a;
    if a.abs() < 15.0 {
        // |a| < 15: shift to [0..1] and use rational approximation.
        let mut t = 1.0;
        let m = (a as i64) - 1;
        if m >= 0 {
            // a >= 1: peel off factors a-1, a-2, ..., a-m (loop empty when m == 0).
            for _ in 1..=m {
                x -= 1.0;
                t *= x;
            }
            x -= 1.0;
        } else {
            // a < 1: shift up.
            t = a;
            if a <= 0.0 {
                let mneg = -m - 1;
                for _ in 1..=mneg {
                    x += 1.0;
                    t *= x;
                }
                x = (x + 0.5) + 0.5;
                t *= x;
                if t == 0.0 {
                    // a is a non-positive integer (0, -1, -2, …): pole.
                    return Err(GammaDomainError::Pole(a));
                }
            }
            // For 1/t overflow check, CDFLIB compares fabs(t)*MAX <= 1.0001.
            if t.abs() < 1e-30 && t.abs() * f64::MAX <= 1.0001 {
                return Err(GammaDomainError::Overflow(a));
            }
            if t.abs() < 1e-30 {
                return Ok(1.0 / t);
            }
        }
        // Compute Γ(1 + x) for 0 <= x < 1 via the rational form.
        let mut top = P_COEF[0];
        let mut bot = Q_COEF[0];
        for i in 1..7 {
            top = P_COEF[i] + x * top;
            bot = Q_COEF[i] + x * bot;
        }
        let g = top / bot;
        return Ok(if a < 1.0 { g / t } else { g * t });
    }

    // |a| >= 15: asymptotic. F90 cdflib.f90:10162-10164 returns its
    // sentinel 0 for |a| >= 1000 regardless of sign; we split by sign
    // so the error variant carries the correct meaning.
    if a >= 1e3 {
        return Err(GammaDomainError::Overflow(a));
    }
    if a <= -1e3 {
        return Err(GammaDomainError::Underflow(a));
    }
    // For negative a, build the reflection coefficient s and switch x = -a.
    // F90 cdflib.f90:10166-10186.
    let mut s = 0.0;
    if a <= 0.0 {
        x = -a;
        let n = x as i64;
        let mut t = x - n as f64;
        if t > 0.9 {
            t = 1.0 - t;
        }
        s = (PI * t).sin() / PI;
        if n % 2 == 0 {
            s = -s;
        }
        if s == 0.0 {
            // a is a negative integer of magnitude ≥ 15: pole.
            return Err(GammaDomainError::Pole(a));
        }
    }
    // Modified asymptotic sum, F90 cdflib.f90:10190-10208.
    let t = 1.0 / (x * x);
    let g = ((((R1 * t + R2) * t + R3) * t + R4) * t + R5) / x;
    let lnx = x.ln();
    let g = D + g + (x - 0.5) * (lnx - 1.0);
    let w = g;
    let t = g - w;
    if w > 0.99999 * POS_EXPARG {
        return Err(GammaDomainError::Overflow(a));
    }
    let mut result = w.exp() * (1.0 + t);
    if a < 0.0 {
        result = (1.0 / (result * s)) / x;
    }
    Ok(result)
}

/// Returns ln Γ(1 + *a*) for −0.2 ≤ *a* ≤ 1.25.
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
        let num = ((((((P6 * a + P5) * a + P4) * a + P3) * a + P2) * a + P1) * a) + P0;
        let den = ((((((Q6 * a + Q5) * a + Q4) * a + Q3) * a + Q2) * a + Q1) * a) + 1.0;
        -(a * (num / den))
    } else {
        let x = a - 0.5 - 0.5;
        let num = (((((R5 * x + R4) * x + R3) * x + R2) * x + R1) * x) + R0;
        let den = (((((S5 * x + S4) * x + S3) * x + S2) * x + S1) * x) + 1.0;
        x * (num / den)
    }
}

/// Errors of [`try_psi`].
///
/// In CDFLIB these were signalled by returning the sentinel `0.0`; the
/// F90 doc header says “PSI is assigned the value 0 when the psi function
/// is undefined”. Here each failure mode is lifted into a typed Rust error.
///
/// [`try_psi`]: crate::special::try_psi
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum PsiError {
    /// Argument is zero or a non-positive integer; *ψ* has a pole there.
    #[error("ψ has a pole at {0}")]
    Pole(f64),
    /// Argument is too large in magnitude on the negative side for the
    /// cotangent reflection to remain accurate; *ψ* effectively
    /// overflows.
    #[error("ψ argument {0} is too large in magnitude (overflow)")]
    Overflow(f64),
}

/// Returns the digamma function *ψ*(*x*) = d/d*x* ln Γ(*x*).
///
/// Port of `psi` (FUNPACK / Cody–Strecok–Thacher, modified by Morris).
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
    if xx.is_nan() {
        return f64::NAN;
    }
    try_psi(xx).unwrap_or_else(|e| panic!("psi({xx}): {e}"))
}

/// Fallible form of [`psi`]: returns [`PsiError`] when the argument lands
/// on a pole or overflows the cotangent reduction.
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
pub fn try_psi(xx: f64) -> Result<f64, PsiError> {
    const DX0: f64 = 1.461632144968362341262659542325721325;
    const PIOV4: f64 = 0.785398163397448;
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

    // F90 cdflib.f90:13212-13213 builds xmax1 in two steps:
    //   xmax1 = real(ipmpar(3), kind=rk)
    //   xmax1 = min(xmax1, 1.0 / epsilon(xmax1))
    // With the IEEE 754 binary64 configuration of ipmpar bundled in the
    // F90 source, ipmpar(3) is the largest int32, 2_147_483_647 (about
    // 2.15e9). The second line caps that at 1/EPSILON (about 4.5e15),
    // but 2.15e9 < 4.5e15 so the cap never fires. The effective value
    // is therefore the constant below.
    let xmax1 = 2_147_483_647.0_f64;
    let xsmall = 1.0e-9;

    let mut x = xx;
    let mut aug = 0.0;
    if x == 0.0 {
        return Err(PsiError::Pole(xx));
    }
    if x < 0.5 {
        // Reflection: ψ(1-x) = ψ(x) + π cot(πx).
        if x.abs() <= xsmall {
            aug = -1.0 / x;
        } else {
            // Argument reduction for cotan.
            let mut w = -x;
            let mut sgn = PIOV4;
            if w <= 0.0 {
                w = -w;
                sgn = -sgn;
            }
            if w >= xmax1 {
                return Err(PsiError::Overflow(xx));
            }
            let mut nq = w as i64;
            w -= nq as f64;
            nq = (w * 4.0) as i64;
            w = 4.0 * (w - (nq as f64) * 0.25);
            let mut n = nq / 2;
            if n + n != nq {
                w = 1.0 - w;
            }
            let z = PIOV4 * w;
            let m = n / 2;
            if m + m != n {
                sgn = -sgn;
            }
            // -π cot(πx)
            n = (nq + 1) / 2;
            let m2 = n / 2;
            if m2 + m2 == n {
                if z == 0.0 {
                    return Err(PsiError::Pole(xx));
                }
                aug = 4.0 * sgn * (z.cos() / z.sin());
            } else {
                aug = 4.0 * sgn * (z.sin() / z.cos());
            }
        }
        x = 1.0 - x;
    }

    if x <= 3.0 {
        // 0.5 ≤ x ≤ 3: rational approximation around x = dx0.
        let mut den = x;
        let mut upper = P1[0] * x;
        for i in 1..=5 {
            den = (den + Q1[i - 1]) * x;
            upper = (upper + P1[i]) * x;
        }
        let den = (upper + P1[6]) / (den + Q1[5]);
        Ok(den * (x - DX0) + aug)
    } else if x < xmax1 {
        // 3 < x < xmax1: asymptotic.
        let w = (1.0 / x) / x;
        let mut den = w;
        let mut upper = P2[0] * w;
        for i in 1..=3 {
            den = (den + Q2[i - 1]) * w;
            upper = (upper + P2[i]) * w;
        }
        let aug = upper / (den + Q2[3]) - 0.5 / x + aug;
        Ok(aug + x.ln())
    } else {
        // xmax1 ≤ x
        Ok(aug + x.ln())
    }
}

/// Returns ln Γ(*a*) for *a* > 0.
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
pub fn gamma_log(a: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-1;
    const C1: f64 = -0.277777777760991e-2;
    const C2: f64 = 0.793650666825390e-3;
    const C3: f64 = -0.595202931351870e-3;
    const C4: f64 = 0.837308034031215e-3;
    const C5: f64 = -0.165322962780713e-2;
    const D: f64 = 0.418938533204673;

    if a <= 0.8 {
        return gamma_ln1(a) - a.ln();
    }
    if a <= 2.25 {
        let t = a - 0.5 - 0.5;
        return gamma_ln1(t);
    }
    if a < 10.0 {
        let n = (a - 1.25) as i64;
        let mut t = a;
        let mut w = 1.0;
        for _ in 1..=n {
            t -= 1.0;
            w *= t;
        }
        return gamma_ln1(t - 1.0) + w.ln();
    }
    // a >= 10: asymptotic. The F90 computes (1/a)**2, rounding the
    // reciprocal before squaring; 1/(a*a) differs in the last ulp.
    let t = (1.0 / a).powi(2);
    let w = (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a;
    D + w + (a - 0.5) * (a.ln() - 1.0)
}

/// Returns ln Γ(*a* + *b*) for 1 ≤ *a* ≤ 2 and 1 ≤ *b* ≤ 2.
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

/// Returns the Stirling remainder ln Γ(*z*) − Stirling(*z*) for *z* > 0, where
/// Stirling(*z*) = ln √(2π) + (*z* − ½) ln *z* − *z*.
///
/// For *z* > 6 uses a 9-term series in Bernoulli numbers; otherwise
/// computes the difference explicitly via [`gamma_log`].
///
/// # Panics
///
/// Panics if *z* ≤ 0 (the F90 routine prints and calls `stop`).
///
/// [`gamma_log`]: crate::special::gamma_log
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
#[inline]
pub fn dstrem(z: f64) -> f64 {
    const NCOEF: usize = 9;
    // F90 coef(0:9): a leading 0 followed by 9 Bernoulli coefficients.
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
    const HLN2PI: f64 = 0.91893853320467274178; // ½ ln(2π)

    if z <= 0.0 {
        panic!("dstrem: argument z must be positive (got {z})");
    }

    if z > 6.0 {
        eval_pol(&COEF, 1.0 / (z * z)) * z
    } else {
        let sterl = HLN2PI + (z - 0.5) * z.ln() - z;
        gamma_log(z) - sterl
    }
}

/// Returns exp(−*x*) · *xᵃ* / Γ(*a*). Used inside [`gamma_inc`] to multiply the
/// regularized integral.
///
/// [`gamma_inc`]: crate::special::gamma_inc
#[inline]
pub fn rcomp(a: f64, x: f64) -> f64 {
    const RT2PIN: f64 = 0.398942280401433;
    if a < 20.0 {
        let t = a * x.ln() - x;
        if a < 1.0 {
            a * t.exp() * (1.0 + gam1(a))
        } else {
            t.exp() / gamma(a)
        }
    } else {
        let u = x / a;
        if u == 0.0 {
            return 0.0;
        }
        let t = (1.0 / a).powi(2);
        let mut t1 = (((0.75 * t - 1.0) * t + 3.5) * t - 105.0) / (a * 1260.0);
        t1 -= a * rlog(u);
        RT2PIN * a.sqrt() * t1.exp()
    }
}

// =====================================================================
// gamma_inc: regularized incomplete gamma P(a,x), Q(a,x)
// =====================================================================

/// Accuracy regime for [`gamma_inc_with_acc`] and [`try_gamma_inc_with_acc`].
///
/// Mirrors CDFLIB's `ind` selector (cdflib.f90:10218). The free-standing
/// [`gamma_inc`] and [`try_gamma_inc`] entry points are equivalent to
/// passing [`GammaIncAcc::Max`].
///
/// [`gamma_inc`]: crate::special::gamma_inc
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
/// [`gamma_inc_with_acc`]: crate::special::gamma_inc_with_acc
/// [`try_gamma_inc_with_acc`]: crate::special::try_gamma_inc_with_acc
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GammaIncAcc {
    /// As much accuracy as possible (F90 ind = 0).
    #[default]
    Max,
    /// Within one unit of the 6th significant digit (F90 ind = 1).
    Digits6,
    /// Within one unit of the 3rd significant digit (F90 ind = 2).
    Digits3,
}

/// Errors of [`gamma_inc`].
///
/// CDFLIB's `gamma_inc` signals invalid input and indeterminate results
/// by returning a sentinel `(2.0, 0.0)` tuple; this enum lifts each
/// failure mode into a typed Rust error.
///
/// [`gamma_inc`]: crate::special::gamma_inc
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaIncError {
    /// The parameter *a* was negative.
    #[error("parameter a must be non-negative, got {0}")]
    ANegative(f64),
    /// The argument *x* was negative.
    #[error("argument x must be non-negative, got {0}")]
    XNegative(f64),
    /// Both *a* and *x* were zero, leaving the result undefined.
    #[error("both a and x are zero")]
    BothZero,
    /// The deep asymptotic regime where the Tricomi–Temme expansion cannot
    /// resolve *P* versus *Q* (triggered when *a* · ε² > 3.28·10⁻³ and
    /// *x* ≈ *a*).
    #[error("indeterminate at a = {a}, x = {x} (deep asymptotic regime)")]
    Indeterminate { a: f64, x: f64 },
}

/// Returns the regularized incomplete Γ function as (*P*(*a*, *x*),
/// *Q*(*a*, *x*)) with *P* + *Q* = 1 (computed independently, not via
/// subtraction).
///
/// The returned pair (*p*, *q*) is the (lower-tail, upper-tail)
/// probability, analogous to the (*w*, *w*₁) pair returned by
/// [`beta_inc`]; both are computed independently rather than one from
/// the other so that the small tail keeps its precision.
///
/// Five computational regimes are selected based on the location of
/// (*a*, *x*) in the parameter space.
///
/// # Panics
///
/// Panics on a [`GammaIncError`]. Use [`try_gamma_inc`] for the fallible
/// form.
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
/// [`beta_inc`]: crate::special::beta_inc
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
#[inline]
pub fn gamma_inc(a: f64, x: f64) -> (f64, f64) {
    gamma_inc_with_acc(a, x, GammaIncAcc::Max)
}

/// Accuracy-selectable form of [`gamma_inc`].
///
/// Lower accuracy uses shallower truncations of the Tricomi-Temme
/// expansion and looser regime cutoffs. Most callers want
/// [`GammaIncAcc::Max`]; the other levels exist for callers that
/// can trade precision for speed (see CDFLIB's `ind` parameter).
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
#[inline]
pub fn gamma_inc_with_acc(a: f64, x: f64, accuracy: GammaIncAcc) -> (f64, f64) {
    if a.is_nan() || x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    try_gamma_inc_with_acc(a, x, accuracy)
        .unwrap_or_else(|e| panic!("gamma_inc_with_acc({a}, {x}, {accuracy:?}): {e}"))
}

/// Fallible form of [`gamma_inc`]: returns [`GammaIncError`] on invalid
/// input or in the deep asymptotic regime where the Tricomi–Temme expansion
/// cannot resolve *P* and *Q*.
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

/// Accuracy-selectable form of [`try_gamma_inc`].
///
/// [`try_gamma_inc`]: crate::special::try_gamma_inc
pub fn try_gamma_inc_with_acc(
    a: f64,
    x: f64,
    accuracy: GammaIncAcc,
) -> Result<(f64, f64), GammaIncError> {
    if a < 0.0 {
        return Err(GammaIncError::ANegative(a));
    }
    if x < 0.0 {
        return Err(GammaIncError::XNegative(x));
    }
    if a == 0.0 && x == 0.0 {
        return Err(GammaIncError::BothZero);
    }
    let (p, q) = gamma_inc_core(a, x, accuracy);
    if p == 2.0 {
        return Err(GammaIncError::Indeterminate { a, x });
    }
    Ok((p, q))
}

/// Core dispatcher of [`gamma_inc`] without input validation.
///
/// An output of `(2.0, 0.0)` represents the indeterminate-result sentinel
/// that the public wrapper lifts into [`GammaIncError::Indeterminate`].
fn gamma_inc_core(a: f64, x: f64, accuracy: GammaIncAcc) -> (f64, f64) {
    use GammaIncAcc::{Digits3, Digits6, Max};
    // The Tricomi-Temme coefficient tables D0..D6 live inside the dedicated
    // helpers (temme_general, temme_for_l_eq_1); the dispatcher only needs
    // ALOG10 (the log(10) cutoff between series and continued fraction) and
    // the per-regime tolerances and cutoffs (cdflib.f90:10342).
    const ALOG10: f64 = 2.30258509299405;
    const RT2PIN: f64 = 0.398942280401433;

    let (acc_raw, e0, x0, big): (f64, f64, f64, f64) = match accuracy {
        Max => (5.0e-15, 0.25e-3, 31.0, 20.0),
        Digits6 => (5.0e-7, 0.25e-1, 17.0, 14.0),
        Digits3 => (5.0e-4, 0.14, 9.7, 10.0),
    };
    let acc = acc_raw.max(f64::EPSILON);
    let e = f64::EPSILON;

    if a * x == 0.0 {
        return if x <= a { (0.0, 1.0) } else { (1.0, 0.0) };
    }

    // ---------- Compute r and dispatch to the right tail formula. ----

    let r;
    if a < 1.0 {
        if a == 0.5 {
            // S390: reduces to erf/erfc.
            let rtx = x.sqrt();
            return if x < 0.25 {
                let ans = error_f(rtx);
                (ans, 0.5 + (0.5 - ans))
            } else {
                let qans = error_fc(rtx);
                (0.5 + (0.5 - qans), qans)
            };
        }
        if x < 1.1 {
            // S160: Taylor series for P(a, x) / x^a.
            return taylor_p_over_xa(a, x, acc);
        }
        let t1 = a * x.ln() - x;
        let u = a * t1.exp();
        if u == 0.0 {
            return (1.0, 0.0);
        }
        r = u * (1.0 + gam1(a));
        // Continued fraction (S250).
        return continued_fraction(a, x, r, acc, e);
    }

    // a >= 1
    if a < big {
        // Finite-sum branch for "a >= 1 and 2a is an integer", restricted
        // to a <= x < x0 where it's profitable.
        if a <= x && x < x0 {
            let twoa = a + a;
            let m = twoa as i64;
            if twoa == m as f64 {
                let ii = m / 2;
                let a_eq_i = a == ii as f64;
                return finite_q_half_integer(a, x, ii, a_eq_i);
            }
        }
        // S20: r = exp(t1) / Γ(a).
        let t1 = a * x.ln() - x;
        r = t1.exp() / gamma(a);
    } else {
        // a >= big
        let l = x / a;
        if l == 0.0 {
            return (0.0, 1.0); // S370
        }
        let s = 0.5 + (0.5 - l);
        let z = rlog(l);
        if z >= 700.0 / a {
            // S410
            if s.abs() <= 2.0 * e {
                return (2.0, 0.0);
            }
            return if x <= a { (0.0, 1.0) } else { (1.0, 0.0) };
        }
        let y = a * z;
        let rta = a.sqrt();
        if s.abs() <= e0 / rta {
            return temme_for_l_eq_1(a, l, z, y, e, accuracy);
        }
        if s.abs() <= 0.4 {
            return temme_general(a, l, z, y, rta, e, accuracy);
        }
        let t = (1.0 / a).powi(2);
        let mut t1 = (((0.75 * t - 1.0) * t + 3.5) * t - 105.0) / (a * 1260.0);
        t1 -= y;
        r = RT2PIN * rta * t1.exp();
    }

    // S40: post-r dispatch.
    if r == 0.0 {
        return if x <= a { (0.0, 1.0) } else { (1.0, 0.0) };
    }
    if x <= a.max(ALOG10) {
        // S50: Taylor series for P / r.
        return taylor_p_over_r(a, x, r, acc);
    }
    if x < x0 {
        return continued_fraction(a, x, r, acc, e);
    }
    // S100: asymptotic expansion (large x).
    asymptotic_expansion_q(a, x, r, acc)
}

// ---------- gamma_inc helper code paths ----------

fn taylor_p_over_xa(a: f64, x: f64, acc: f64) -> (f64, f64) {
    // S160 in CDFLIB.
    let mut an: f64 = 3.0;
    let mut c = x;
    let mut sum = x / (a + 3.0);
    let tol = 3.0 * acc / (a + 1.0);
    loop {
        an += 1.0;
        c = -(c * (x / an));
        let t = c / (a + an);
        sum += t;
        if t.abs() <= tol {
            break;
        }
    }
    let j = a * x * ((sum / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));
    let z = a * x.ln();
    let h = gam1(a);
    let g = 1.0 + h;
    let use_label_200 = if x < 0.25 { z > -0.13394 } else { a < x / 2.59 };
    if use_label_200 {
        let l = rexp(z);
        let w = 0.5 + (0.5 + l);
        let qans = (w * j - l) * g - h;
        if qans < 0.0 {
            return (1.0, 0.0);
        }
        let ans = 0.5 + (0.5 - qans);
        (ans, qans)
    } else {
        let w = z.exp();
        let ans = w * g * (0.5 + (0.5 - j));
        let qans = 0.5 + (0.5 - ans);
        (ans, qans)
    }
}

fn taylor_p_over_r(a: f64, x: f64, r: f64, acc: f64) -> (f64, f64) {
    // S50: Taylor series for P/r when r > 0.
    let mut wk = [0.0_f64; 20];
    let mut apn = a + 1.0;
    let mut t = x / apn;
    wk[0] = t;
    let mut n_filled = 20;
    for n in 2..=20 {
        apn += 1.0;
        t *= x / apn;
        if t <= 1e-3 {
            n_filled = n;
            break;
        }
        wk[n - 1] = t;
    }
    let mut sum = t;
    let tol = 0.5 * acc;
    loop {
        apn += 1.0;
        t *= x / apn;
        sum += t;
        if t <= tol {
            break;
        }
    }
    // wk[0..n_filled-1] backfill (reverse).
    let max = n_filled - 1;
    for m in 1..=max {
        let n = n_filled - m;
        sum += wk[n - 1];
    }
    let ans = r / a * (1.0 + sum);
    let qans = 0.5 + (0.5 - ans);
    (ans, qans)
}

fn continued_fraction(a: f64, x: f64, r: f64, acc: f64, e: f64) -> (f64, f64) {
    // S250.
    let tol = (5.0 * e).max(acc);
    let mut a2n: f64 = 1.0;
    let mut a2nm1: f64 = 1.0;
    let mut b2nm1 = x;
    let mut b2n = x + (1.0 - a);
    let mut c: f64 = 1.0;
    loop {
        a2nm1 = x * a2n + c * a2nm1;
        b2nm1 = x * b2n + c * b2nm1;
        let am0 = a2nm1 / b2nm1;
        c += 1.0;
        let cma = c - a;
        a2n = a2nm1 + cma * a2n;
        b2n = b2nm1 + cma * b2n;
        let an0 = a2n / b2n;
        if (an0 - am0).abs() < tol * an0 {
            let qans = r * an0;
            let ans = 0.5 + (0.5 - qans);
            return (ans, qans);
        }
    }
}

fn asymptotic_expansion_q(a: f64, x: f64, r: f64, acc: f64) -> (f64, f64) {
    // S100.
    let mut wk = [0.0_f64; 20];
    let mut amn = a - 1.0;
    let mut t = amn / x;
    wk[0] = t;
    let mut n_filled = 20;
    for n in 2..=20 {
        amn -= 1.0;
        t *= amn / x;
        if t.abs() <= 1e-3 {
            n_filled = n;
            break;
        }
        wk[n - 1] = t;
    }
    let mut sum = t;
    while t.abs() > acc {
        amn -= 1.0;
        t *= amn / x;
        sum += t;
    }
    let max = n_filled - 1;
    for m in 1..=max {
        let n = n_filled - m;
        sum += wk[n - 1];
    }
    let qans = r / x * (1.0 + sum);
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

fn finite_q_half_integer(_a: f64, x: f64, i: i64, a_eq_i: bool) -> (f64, f64) {
    // S210 / S220 in CDFLIB. "a >= 1 and 2a is an integer" path.
    let mut sum;
    let mut t;
    let mut n;
    let mut c;
    if a_eq_i {
        // S210: integer a.
        sum = (-x).exp();
        t = sum;
        n = 1_i64;
        c = 0.0_f64;
    } else {
        // S220: half-integer a.
        let rtx = x.sqrt();
        sum = error_fc(rtx);
        t = (-x).exp() / (1.77245385090552 * rtx);
        n = 0_i64;
        c = -0.5_f64;
    }
    while n != i {
        n += 1;
        c += 1.0;
        t = x * t / c;
        sum += t;
    }
    let qans = sum;
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

fn temme_general(
    a: f64,
    l: f64,
    z_in: f64,
    y: f64,
    rta: f64,
    e: f64,
    accuracy: GammaIncAcc,
) -> (f64, f64) {
    use GammaIncAcc::{Digits3, Digits6, Max};
    // S270: general Tricomi-Temme expansion (|s| <= 0.4, s = 1 - l). The
    // accuracy regime selects the truncation depth of the polynomial in z
    // (cdflib.f90:10758, L10873, L10888): Max uses the long expansion (with
    // a short sub-branch for |s| <= 1e-3), Digits6 a moderate expansion,
    // Digits3 the shortest.
    //
    // Indeterminate-sentinel check (cdflib.f90:10744, cdflib.f:840): when s ≈
    // 0 and a is so large that a*ε² > 3.28e-3, the Tricomi-Temme expansion
    // cannot resolve P/Q reliably. Returns the sentinel 2.0.
    let s = 0.5 + (0.5 - l);
    if s.abs() <= 2.0 * e && 3.28e-3 < a * e * e {
        return (2.0, 0.0);
    }
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
    const D4: [f64; 6] = [
        0.784039221720067e-3,
        -0.299072480303190e-3,
        -0.146384525788434e-5,
        0.664149821546512e-4,
        -0.396836504717943e-4,
        0.113757269706784e-4,
    ];
    const D5: [f64; 4] = [
        -0.697281375836586e-4,
        0.277275324495939e-3,
        -0.199325705161888e-3,
        0.679778047793721e-4,
    ];
    const D6: [f64; 2] = [-0.592166437353694e-3, 0.270878209671804e-3];
    const D10: f64 = -0.185185185185185e-2;
    const D20: f64 = 0.413359788359788e-2;
    const D30: f64 = 0.649434156378601e-3;
    const D40: f64 = -0.861888290916712e-3;
    const D50: f64 = -0.336798553366358e-3;
    const D60: f64 = 0.531307936463992e-3;
    const D70: f64 = 0.344367606892378e-3;
    const THIRD: f64 = 1.0 / 3.0;
    const RT2PIN: f64 = 0.398942280401433;

    let c = (-y).exp();
    let w = 0.5 * error_fc_scaled(y.sqrt());
    let u = 1.0 / a;
    let mut z = (z_in + z_in).sqrt();
    if l < 1.0 {
        z = -z;
    }
    let t = match accuracy {
        Max => {
            // S270 / S280 (max accuracy, F90 iop=1). F90 splits this into
            // two sub-branches at |s| <= 1e-3 (short) versus |s| > 1e-3
            // (long), cdflib.f90:10760.
            if s.abs() <= 1.0e-3 {
                let c0 = ((((((D0[6] * z + D0[5]) * z + D0[4]) * z + D0[3]) * z + D0[2]) * z
                    + D0[1])
                    * z
                    + D0[0])
                    * z
                    - THIRD;
                let c1 =
                    (((((D1[5] * z + D1[4]) * z + D1[3]) * z + D1[2]) * z + D1[1]) * z + D1[0]) * z
                        + D10;
                let c2 = ((((D2[4] * z + D2[3]) * z + D2[2]) * z + D2[1]) * z + D2[0]) * z + D20;
                let c3 = (((D3[3] * z + D3[2]) * z + D3[1]) * z + D3[0]) * z + D30;
                let c4 = (D4[1] * z + D4[0]) * z + D40;
                let c5 = (D5[1] * z + D5[0]) * z + D50;
                let c6 = D6[0] * z + D60;
                ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0
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
                    - THIRD;
                let c1 = (((((((((((D1[11] * z + D1[10]) * z + D1[9]) * z + D1[8]) * z
                    + D1[7])
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
                let c2 = (((((((((D2[9] * z + D2[8]) * z + D2[7]) * z + D2[6]) * z + D2[5])
                    * z
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
                let c3 = (((((((D3[7] * z + D3[6]) * z + D3[5]) * z + D3[4]) * z + D3[3]) * z
                    + D3[2])
                    * z
                    + D3[1])
                    * z
                    + D3[0])
                    * z
                    + D30;
                let c4 =
                    (((((D4[5] * z + D4[4]) * z + D4[3]) * z + D4[2]) * z + D4[1]) * z + D4[0]) * z
                        + D40;
                let c5 = (((D5[3] * z + D5[2]) * z + D5[1]) * z + D5[0]) * z + D50;
                let c6 = (D6[1] * z + D6[0]) * z + D60;
                ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0
            }
        }
        Digits6 => {
            // S290 (within 1 unit of the 6th significant digit, F90 iop=2).
            let c0 = (((((D0[5] * z + D0[4]) * z + D0[3]) * z + D0[2]) * z + D0[1]) * z + D0[0])
                * z
                - THIRD;
            let c1 = (((D1[3] * z + D1[2]) * z + D1[1]) * z + D1[0]) * z + D10;
            let c2 = D2[0] * z + D20;
            (c2 * u + c1) * u + c0
        }
        Digits3 => {
            // S300 (within 1 unit of the 3rd significant digit, F90 iop=3).
            ((D0[2] * z + D0[1]) * z + D0[0]) * z - THIRD
        }
    };

    // Tricomi-Temme normalization with the scaled erfc factor.
    if l < 1.0 {
        let ans = c * (w - RT2PIN * t / rta);
        let qans = 0.5 + (0.5 - ans);
        return (ans, qans);
    }
    let qans = c * (w + RT2PIN * t / rta);
    let ans = 0.5 + (0.5 - qans);
    (ans, qans)
}

fn temme_for_l_eq_1(
    a: f64,
    l: f64,
    z_in: f64,
    y: f64,
    e: f64,
    accuracy: GammaIncAcc,
) -> (f64, f64) {
    use GammaIncAcc::{Digits3, Digits6, Max};
    // S330 / S340 / S350 / S360 in CDFLIB (cdflib.f90:10908). The accuracy
    // regime selects the truncation depth: Max uses the full c0..c6 + d70
    // expansion, Digits6 a shallow expansion, Digits3 the shallowest.
    //
    // Indeterminate-sentinel check (cdflib.f90:10910, cdflib.f:915): when
    // a*ε² > 3.28e-3, the Tricomi-Temme L=1 expansion cannot resolve P/Q.
    if 3.28e-3 < a * e * e {
        return (2.0, 0.0);
    }
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
    const D4: [f64; 6] = [
        0.784039221720067e-3,
        -0.299072480303190e-3,
        -0.146384525788434e-5,
        0.664149821546512e-4,
        -0.396836504717943e-4,
        0.113757269706784e-4,
    ];
    const D5: [f64; 4] = [
        -0.697281375836586e-4,
        0.277275324495939e-3,
        -0.199325705161888e-3,
        0.679778047793721e-4,
    ];
    const D6: [f64; 2] = [-0.592166437353694e-3, 0.270878209671804e-3];
    const D10: f64 = -0.185185185185185e-2;
    const D20: f64 = 0.413359788359788e-2;
    const D30: f64 = 0.649434156378601e-3;
    const D40: f64 = -0.861888290916712e-3;
    const D50: f64 = -0.336798553366358e-3;
    const D60: f64 = 0.531307936463992e-3;
    const D70: f64 = 0.344367606892378e-3;
    const THIRD: f64 = 1.0 / 3.0;
    const RT2PIN: f64 = 0.398942280401433;
    const RTPI: f64 = 1.77245385090552;

    let c = 0.5 + (0.5 - y);
    let w = (0.5 - y.sqrt() * (0.5 + (0.5 - y / 3.0)) / RTPI) / c;
    let u = 1.0 / a;
    let mut z = (z_in + z_in).sqrt();
    if l < 1.0 {
        z = -z;
    }
    let t = match accuracy {
        Max => {
            // S330 / S340 (F90 iop=1).
            let c0 = ((((((D0[6] * z + D0[5]) * z + D0[4]) * z + D0[3]) * z + D0[2]) * z + D0[1])
                * z
                + D0[0])
                * z
                - THIRD;
            let c1 = (((((D1[5] * z + D1[4]) * z + D1[3]) * z + D1[2]) * z + D1[1]) * z + D1[0])
                * z
                + D10;
            let c2 = ((((D2[4] * z + D2[3]) * z + D2[2]) * z + D2[1]) * z + D2[0]) * z + D20;
            let c3 = (((D3[3] * z + D3[2]) * z + D3[1]) * z + D3[0]) * z + D30;
            let c4 = (D4[1] * z + D4[0]) * z + D40;
            let c5 = (D5[1] * z + D5[0]) * z + D50;
            let c6 = D6[0] * z + D60;
            ((((((D70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0
        }
        Digits6 => {
            // S350 (F90 iop=2).
            let c0 = (D0[1] * z + D0[0]) * z - THIRD;
            let c1 = D1[0] * z + D10;
            (D20 * u + c1) * u + c0
        }
        Digits3 => {
            // S360 (F90 iop=3).
            D0[0] * z - THIRD
        }
    };
    let rta = a.sqrt();
    if l < 1.0 {
        let ans = c * (w - RT2PIN * t / rta);
        (ans, 0.5 + (0.5 - ans))
    } else {
        let qans = c * (w + RT2PIN * t / rta);
        (0.5 + (0.5 - qans), qans)
    }
}

// =====================================================================
// gamma_inc_inv: inverse of the regularized incomplete gamma ratio
// =====================================================================

/// Errors of [`gamma_inc_inv`].
///
/// CDFLIB's `gamma_inc_inv` reports outcomes via a signed integer status;
/// this enum gives each one a named, descriptive form. Variants whose F90
/// counterpart still produces a usable approximation carry it as a field,
/// so callers that want to ignore the soft-failure can recover the value.
///
/// [`gamma_inc_inv`]: crate::special::gamma_inc_inv
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum GammaIncInvError {
    /// *a* ≤ 0 (CDFLIB `ierr = -2`).
    #[error("parameter a must be positive, got {0}")]
    ANotPositive(f64),
    /// No solution: *q*/*a* is too large (CDFLIB `ierr = -3`).
    #[error("no solution: q/a is too large")]
    NoSolution,
    /// *p* + *q* ≠ 1 within tolerance (CDFLIB `ierr = -4`).
    #[error("inconsistent inputs: p + q must equal 1")]
    InconsistentPq,
    /// 20 Schröder iterations did not converge; only fires when the
    /// caller supplied `x0 > 0` (CDFLIB `ierr = -6`).
    #[error("iteration did not converge in 20 steps; last value: {partial}")]
    NotConverged {
        /// Last running approximation when the iteration limit was hit.
        partial: f64,
    },
    /// Iteration failed: intermediate *x* went non-positive (CDFLIB
    /// `ierr = -7`).
    #[error("iteration failed: intermediate x went non-positive")]
    IterationFailed,
    /// A value was obtained but the routine cannot certify its accuracy
    /// (deep tails or *a* extremely close to *x*; CDFLIB `ierr = -8`).
    #[error("solution obtained but accuracy cannot be certified; value: {value}")]
    UncertainAccuracy {
        /// The approximate solution the routine produced; reliability is
        /// not guaranteed but the value is often usable for warm-starting.
        value: f64,
    },
    /// The mathematical answer is +∞ (this happens when `q == 0`, i.e.
    /// `p == 1` exactly: the incomplete-Γ ratio reaches 1 only in the
    /// limit *x* → +∞).
    #[error("inverse is +∞ (q = 0 means P(a, x) = 1 only as x → +∞)")]
    AtInfinity,
}

/// Returns the inverse of the regularized incomplete Γ function: *x* such
/// that *P*(*a*, *x*) = *p* and *Q*(*a*, *x*) = *q*.
///
/// Uses Schröder iteration; an asymptotic-series approximation provides
/// the initial guess (or the caller supplies `x0 > 0`). Algorithm by
/// Alfred Morris.
///
/// # Panics
///
/// Panics on a [`GammaIncInvError`]. Use [`try_gamma_inc_inv`] for the
/// fallible form. Returns `(value, iterations)`; see [`try_gamma_inc_inv`]
/// for the meaning of the second component.
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
    if a.is_nan() || x0.is_nan() || p.is_nan() || q.is_nan() {
        return (f64::NAN, 0);
    }
    try_gamma_inc_inv(a, x0, p, q)
        .unwrap_or_else(|e| panic!("gamma_inc_inv(a={a}, x0={x0}, p={p}, q={q}): {e}"))
}

/// Fallible form of [`gamma_inc_inv`]: returns `(value, iterations)` on
/// success or [`GammaIncInvError`] otherwise. `iterations` is the Schröder
/// iteration count, mirroring CDFLIB's positive `ierr`: 0 indicates the
/// closed-form path was taken (F90 ierr = 0), positive values mean K
/// Schröder iterations were performed (F90 ierr = K).
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
pub fn try_gamma_inc_inv(a: f64, x0: f64, p: f64, q: f64) -> Result<(f64, u32), GammaIncInvError> {
    // Constants from cdflib.f90:11070 onward.
    const A0: f64 = 3.31125922108741;
    const A1: f64 = 11.6616720288968;
    const A2: f64 = 4.28342155967104;
    const A3: f64 = 0.213623493715853;
    const B1: f64 = 6.61053765625462;
    const B2: f64 = 6.40691597760039;
    const B3: f64 = 1.27364489782223;
    const B4: f64 = 0.036117081018842;
    const C: f64 = 0.577215664901533; // Euler–Mascheroni γ
    const HALF: f64 = 0.5;
    const TWO: f64 = 2.0;
    const LN10: f64 = 2.302585;

    // iop-indexed tables. The F90 uses iop ∈ {1, 2} (one-based); we
    // use iop ∈ {0, 1} to fit Rust array indexing.
    const AMIN: [f64; 2] = [500.0, 100.0];
    const BMIN: [f64; 2] = [1.0e-28, 1.0e-13];
    const DMIN: [f64; 2] = [1.0e-6, 1.0e-4];
    const EMIN: [f64; 2] = [2.0e-3, 6.0e-3];
    const EPS0: [f64; 2] = [1.0e-10, 1.0e-8];

    let e = f64::EPSILON;

    // ---- validation ----------------------------------------------------
    if a <= 0.0 {
        return Err(GammaIncInvError::ANotPositive(a));
    }
    if (p + q - 1.0).abs() > e {
        return Err(GammaIncInvError::InconsistentPq);
    }
    if p == 0.0 {
        return Ok((0.0, 0));
    }
    if q == 0.0 {
        return Err(GammaIncInvError::AtInfinity);
    }
    if a == 1.0 {
        return Ok((if q >= 0.9 { -alnrel(-p) } else { -q.ln() }, 0));
    }

    // ---- setup ---------------------------------------------------------
    let e2 = TWO * e;
    let amax = 0.4e-10 / (e * e);
    let iop = if e > 1.0e-10 { 1 } else { 0 };
    let eps = EPS0[iop];

    // xn is the running approximation; use_q selects which Schröder
    // branch to use after the initial approximation is chosen.
    let xn;
    let use_q;

    if x0 > 0.0 {
        // F90 L11160 path: caller supplied an initial approximation.
        // F90 L11445: branch on p after go to 160.
        xn = x0;
        use_q = p > 0.5;
    } else if a < 1.0 {
        // ---- a < 1 branch (F90 L11185+) -----------------------------
        let g = gamma(a + 1.0);
        let qg = q * g;
        if qg == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: f64::MAX });
        }
        let b = qg / a;

        // F90 L11210: if 0.6*a < qg, jump to label 40 (small-b path).
        let go_to_40 = 0.6 * a < qg;

        if !go_to_40 && a < 0.30 && b >= 0.35 {
            // F90 L11214: closed-form for tiny a, moderate b.
            let t = (-(b + C)).exp();
            let u = t * t.exp();
            xn = t * u.exp();
            // L11217 go to 160: then L11445 → use P (since p must be ≤ 0.5
            // here: for a < 0.3, q is large, so p = 1 − q is small).
            use_q = p > 0.5;
        } else if !go_to_40 && b >= 0.45 {
            // Re-enter the small-b path.
            xn = initial_approx_small_b(a, p, q, g, b, C);
            if xn == 0.0 {
                return Err(GammaIncInvError::NoSolution);
            }
            use_q = p > 0.5;
        } else if !go_to_40 && b == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: f64::MAX });
        } else if !go_to_40 {
            // F90 L11243: y, s, z, t setup.
            let y = -(b.ln());
            let s = HALF + (HALF - a);
            let z = y.ln();
            let t = y - s * z;

            if b >= 0.15 {
                // F90 L11247
                xn = y - s * t.ln() - (1.0 + s / (t + 1.0)).ln();
                use_q = true; // go to 220
            } else if b > 0.01 {
                // F90 L11252
                let u = ((t + TWO * (3.0 - a)) * t + (TWO - a) * (3.0 - a))
                    / ((t + (5.0 - a)) * t + TWO);
                xn = y - s * t.ln() - u.ln();
                use_q = true; // go to 220
            } else {
                // F90 L11261: label 30, c1..c5 expansion.
                let xn_30 = label_30(a, s, z, y);
                // L11269: if a > 1: go to 220. (false here since a < 1.)
                // L11272: if bmin[iop] < b: go to 220.
                if BMIN[iop] < b {
                    xn = xn_30;
                    use_q = true;
                } else {
                    return Ok((xn_30, 0));
                }
            }
        } else {
            // F90 L11279: label 40.
            xn = initial_approx_small_b(a, p, q, g, b, C);
            if xn == 0.0 {
                return Err(GammaIncInvError::NoSolution);
            }
            // F90 L11296 → go to 160.
            use_q = p > 0.5;
        }
    } else {
        // ---- a > 1 branch (F90 L11298: label 80) --------------------
        let w = if q > 0.5 { p.ln() } else { q.ln() };
        let t = (-TWO * w).sqrt();
        let mut s = t
            - (((A3 * t + A2) * t + A1) * t + A0) / ((((B4 * t + B3) * t + B2) * t + B1) * t + 1.0);
        if q > 0.5 {
            s = -s;
        }
        let rta = a.sqrt();
        let s2 = s * s;
        let xn0 = a + s * rta + (s2 - 1.0) / 3.0 + s * (s2 - 7.0) / (36.0 * rta)
            - ((3.0 * s2 + 7.0) * s2 - 16.0) / (810.0 * a)
            + s * ((9.0 * s2 + 256.0) * s2 - 433.0) / (38880.0 * a * rta);
        let xn0 = xn0.max(0.0);

        // F90 L11322: if amin[iop] <= a: check the tolerance early-out.
        if AMIN[iop] <= a {
            let d = HALF + (HALF - xn0 / a);
            if d.abs() <= DMIN[iop] {
                return Ok((xn0, 0));
            }
        }

        if p <= 0.5 {
            // F90 L11340: label 130, a > 1, p ≤ 0.5. The F90 output
            // variable x starts at 0.0 (cdflib.f90:11141) and is only set
            // to xn when the amin[iop] ≤ a branch above (cdflib.f90:11330)
            // runs. label_130's exp-iteration sub-branch uses that x,
            // so we pass the same conditional initial value.
            let x_initial = if AMIN[iop] <= a { xn0 } else { 0.0 };
            let (refined_xn, early_return) = label_130(a, p, q, xn0, x_initial, EMIN[iop]);
            if early_return {
                return Ok((refined_xn, 0));
            }
            xn = refined_xn;
            use_q = false; // go to 170
        } else if xn0 < 3.0 * a {
            // F90 L11356: go to 220.
            xn = xn0;
            use_q = true;
        } else {
            // F90 L11358+
            let y = -(w + gamma_log(a));
            let d = TWO.max(a * (a - 1.0));
            if LN10 * d <= y {
                // F90 L11363: go to 30 (with s and z recomputed).
                let s = 1.0 - a;
                let z = y.ln();
                let xn_30 = label_30(a, s, z, y);
                // L11269: a > 1 here, so go to 220.
                xn = xn_30;
                use_q = true;
            } else {
                let t = a - 1.0;
                let xn1 = y + t * xn0.ln() - alnrel(-t / (xn0 + 1.0));
                let xn2 = y + t * xn1.ln() - alnrel(-t / (xn1 + 1.0));
                xn = xn2;
                use_q = true;
            }
        }
    }

    // ---- Schröder iteration (F90 L11414+ for P, L11506+ for Q) --------
    if use_q {
        schroder_q(a, xn, p, q, eps, amax, e2)
    } else {
        schroder_p(a, xn, p, q, eps, amax, e2)
    }
}

/// F90 L11261: label 30, small-*b* fallback initial approximation
/// using the c1..c5 expansion.
fn label_30(a: f64, s: f64, z: f64, y: f64) -> f64 {
    const HALF: f64 = 0.5;
    const TWO: f64 = 2.0;
    let c1 = -s * z;
    let c2 = -s * (1.0 + c1);
    let c3 = s * ((HALF * c1 + (TWO - a)) * c1 + (2.5 - 1.5 * a));
    let c4 = -s
        * (((c1 / 3.0 + (2.5 - 1.5 * a)) * c1 + ((a - 6.0) * a + 7.0)) * c1
            + ((11.0 * a - 46.0) * a + 47.0) / 6.0);
    let c5 = -s
        * ((((-c1 / 4.0 + (11.0 * a - 17.0) / 6.0) * c1 + ((-3.0 * a + 13.0) * a - 13.0)) * c1
            + HALF * (((TWO * a - 25.0) * a + 72.0) * a - 61.0))
            * c1
            + (((25.0 * a - 195.0) * a + 477.0) * a - 379.0) / 12.0);
    ((((c5 / y + c4) / y + c3) / y + c2) / y + c1) + y
}

/// F90 L11279: label 40, initial approximation for the small-*b* regime.
fn initial_approx_small_b(a: f64, p: f64, q: f64, g: f64, b: f64, gamma_eu: f64) -> f64 {
    const HALF: f64 = 0.5;
    let xn = if b * q <= 1.0e-8 {
        (-(q / a + gamma_eu)).exp()
    } else if p > 0.9 {
        ((alnrel(-q) + gamma_ln1(a)) / a).exp()
    } else {
        ((p * g).ln() / a).exp()
    };
    if xn == 0.0 {
        return 0.0;
    }
    let t = HALF + (HALF - xn / (a + 1.0));
    xn / t
}

/// F90 L11340: label 130, refinement for *a* > 1 and *p* ≤ 0.5.
///
/// Returns the refined `xn` and a flag set to `true` if an early-return
/// from the parent routine is warranted (sub-`emin` approximation).
///
/// `x_initial` matches the F90 output variable `x` at entry to label 130:
/// it is `xn0` if the caller's `amin[iop] <= a` branch ran (F90 cdflib.f90:11330
/// `x = xn`) and `0.0` otherwise (F90's `x = 0.0D+00` initialization at
/// cdflib.f90:11141 is never overwritten).
fn label_130(a: f64, p: f64, _q: f64, xn0: f64, x_initial: f64, emin_iop: f64) -> (f64, bool) {
    let ap1 = a + 1.0;
    if 0.70 * ap1 < xn0 {
        // F90 L11343: go to 170, no refinement needed.
        return (xn0, false);
    }
    let mut w = p.ln() + gamma_log(ap1);
    let mut xn = xn0;
    // The F90 literal 0.15 here carries no D exponent; regenerate.sh
    // compiles with -fdefault-real-8 so it is double precision.
    if xn <= 0.15 * ap1 {
        // F90 L11348: closed-form refinement via three corrective x = exp(...)
        // updates that match the F90 line by line. F90 starts from the
        // output variable x (= x_initial, see fn-doc above).
        let ap2 = a + 2.0;
        let ap3 = a + 3.0;
        let mut x = ((w + x_initial) / a).exp();
        x = ((w + x - (1.0 + (x / ap1) * (1.0 + x / ap2)).ln()) / a).exp();
        x = ((w + x - (1.0 + (x / ap1) * (1.0 + x / ap2)).ln()) / a).exp();
        x = ((w + x - (1.0 + (x / ap1) * (1.0 + (x / ap2) * (1.0 + x / ap3))).ln()) / a).exp();
        xn = x;
        if xn <= 1.0e-2 * ap1 {
            if xn <= emin_iop * ap1 {
                return (xn, true);
            }
            // F90 L11362: go to 170.
            return (xn, false);
        }
    }

    // F90 L11369+: series-sum refinement.
    let mut apn = ap1;
    let mut t = xn / apn;
    let mut sum1 = 1.0 + t;
    loop {
        apn += 1.0;
        t *= xn / apn;
        sum1 += t;
        if t <= 1.0e-4 {
            break;
        }
    }
    let t = w - sum1.ln();
    w = t; // for the next line of the F90 expression
    let mut xn_new = ((xn + w) / a).exp();
    xn_new *= 1.0 - (a * xn_new.ln() - xn_new - w) / (a - xn_new);
    (xn_new, false)
}

/// Schröder iteration in the P branch (F90 L11414: label 170/180).
fn schroder_p(
    a: f64,
    xn0: f64,
    p: f64,
    _q: f64,
    eps: f64,
    amax: f64,
    e2: f64,
) -> Result<(f64, u32), GammaIncInvError> {
    const HALF: f64 = 0.5;
    const TOL: f64 = 1.0e-5;
    if p <= 1.0e10 * f64::MIN_POSITIVE {
        return Err(GammaIncInvError::UncertainAccuracy { value: xn0 });
    }
    let am1 = (a - HALF) - HALF;
    let mut xn = xn0;
    let mut iter: u32 = 0;
    loop {
        // F90 L11432: amax check.
        if amax < a {
            let d = HALF + (HALF - xn / a);
            if d.abs() <= e2 {
                return Err(GammaIncInvError::UncertainAccuracy { value: xn });
            }
        }
        if iter >= 20 {
            return Err(GammaIncInvError::NotConverged { partial: xn });
        }
        iter += 1;
        let (pn, qn) =
            try_gamma_inc(a, xn).map_err(|_| GammaIncInvError::UncertainAccuracy { value: xn })?;
        if pn == 0.0 || qn == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }
        let r = rcomp(a, xn);
        if r == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }
        let t = (pn - p) / r;
        let w = HALF * (am1 - xn);
        let (x, d) = if t.abs() <= 0.1 && (w * t).abs() <= 0.1 {
            // F90 L11468: label 200, second-order Schröder update.
            let h = t * (1.0 + w * t);
            let x = xn * (1.0 - h);
            if x <= 0.0 {
                return Err(GammaIncInvError::IterationFailed);
            }
            if w.abs() >= 1.0 && w.abs() * t * t <= eps {
                return Ok((x, iter));
            }
            (x, h.abs())
        } else {
            // First-order update.
            let x = xn * (1.0 - t);
            if x <= 0.0 {
                return Err(GammaIncInvError::IterationFailed);
            }
            (x, t.abs())
        };
        xn = x;
        if d <= TOL {
            if d <= eps {
                return Ok((xn, iter));
            }
            if (p - pn).abs() <= TOL * p {
                return Ok((xn, iter));
            }
        }
    }
}

/// Schröder iteration in the Q branch (F90 L11506: label 220/230).
/// Structurally symmetric to [`schroder_p`].
fn schroder_q(
    a: f64,
    xn0: f64,
    _p: f64,
    q: f64,
    eps: f64,
    amax: f64,
    e2: f64,
) -> Result<(f64, u32), GammaIncInvError> {
    const HALF: f64 = 0.5;
    const TOL: f64 = 1.0e-5;
    if q <= 1.0e10 * f64::MIN_POSITIVE {
        return Err(GammaIncInvError::UncertainAccuracy { value: xn0 });
    }
    let am1 = (a - HALF) - HALF;
    let mut xn = xn0;
    let mut iter: u32 = 0;
    loop {
        if amax < a {
            let d = HALF + (HALF - xn / a);
            if d.abs() <= e2 {
                return Err(GammaIncInvError::UncertainAccuracy { value: xn });
            }
        }
        if iter >= 20 {
            return Err(GammaIncInvError::NotConverged { partial: xn });
        }
        iter += 1;
        let (pn, qn) =
            try_gamma_inc(a, xn).map_err(|_| GammaIncInvError::UncertainAccuracy { value: xn })?;
        if pn == 0.0 || qn == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }
        let r = rcomp(a, xn);
        if r == 0.0 {
            return Err(GammaIncInvError::UncertainAccuracy { value: xn });
        }
        let t = (q - qn) / r;
        let w = HALF * (am1 - xn);
        // The F90 q branch writes these two 0.1 thresholds without a D
        // exponent; regenerate.sh compiles with -fdefault-real-8 so they
        // are double precision like the p branch's 0.1D+00.
        let (x, d) = if t.abs() <= 0.1 && (w * t).abs() <= 0.1 {
            let h = t * (1.0 + w * t);
            let x = xn * (1.0 - h);
            if x <= 0.0 {
                return Err(GammaIncInvError::IterationFailed);
            }
            if w.abs() >= 1.0 && w.abs() * t * t <= eps {
                return Ok((x, iter));
            }
            (x, h.abs())
        } else {
            let x = xn * (1.0 - t);
            if x <= 0.0 {
                return Err(GammaIncInvError::IterationFailed);
            }
            (x, t.abs())
        };
        xn = x;
        if d > TOL {
            continue;
        }
        if d <= eps {
            return Ok((xn, iter));
        }
        if (q - qn).abs() <= TOL * q {
            return Ok((xn, iter));
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
        // At z = 100 the next term is ~3·10⁻⁸ relative, so the
        // leading-term match is ~4·10⁻⁵.
        let r = dstrem(100.0);
        let lead = 1.0 / 1200.0;
        assert!((r - lead).abs() / lead < 1e-4, "r = {r}, leading = {lead}");
    }

    #[test]
    fn dstrem_small_z_matches_explicit_difference() {
        // For z ≤ 6, dstrem uses gamma_log(z) − Stirling(z) directly.
        // At z = 5: lnΓ(5) = ln 24 = 3.178053830347946...,
        // Stirling(5) ≈ ½ ln(2π) + 4.5·ln 5 − 5 = 3.161549...,
        // so dstrem(5) ≈ 0.0165...
        let r = dstrem(5.0);
        let stirling = 0.91893853320467274178 + 4.5 * 5.0_f64.ln() - 5.0;
        let expected = gamma_log(5.0) - stirling;
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
        // Fine grid exercising the initial-approximation tree's regimes:
        //   a < 0.3 with various q levels (different b = qg/a ranges)
        //   0.3 < a < 1 (the label-30 c1..c5 fallback)
        //   1 < a < 100 (the rational a > 1 path)
        //   a > 500 (the early-return-via-dmin path)
        let grid_a = [
            0.01_f64, 0.05, 0.1, 0.25, 0.5, 0.95, // a < 1
            1.5, 2.0, 5.0, 25.0, 99.0, // moderate a > 1
            150.0, 600.0, 2000.0, // a above amin[iop]
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

    // ============================================================ existing

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
                if p == 2.0 {
                    continue; // error sentinel
                }
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

    // ===== Validation & error sentinels =====

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
    fn gamma_inc_a_zero_x_positive() {
        // a = 0, x > 0: P(0, x) = 1 (the limit).
        let (p, q) = gamma_inc(0.0, 1.0);
        assert_eq!(p, 1.0);
        assert_eq!(q, 0.0);
    }

    // ===== Regime switching points =====

    #[test]
    fn gamma_inc_regime_a_lt_1_x_lt_1() {
        // Power-series small-a small-x regime.
        let (p, q) = gamma_inc(0.3, 0.2);
        assert!((p + q - 1.0).abs() < 1e-14);
        assert!(p > 0.0 && p < 1.0);
    }

    #[test]
    fn gamma_inc_regime_a_lt_1_x_ge_1() {
        // Different branch: small a, larger x.
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
        // Tricomi–Temme asymptotic regime: a ≥ ~20.
        let (p, q) = gamma_inc(100.0, 100.0);
        assert!((p + q - 1.0).abs() < 1e-12);
        // Gamma(a, 1) at x=a: just above the median (which is at
        // ≈ a-1/3 for large a), so cdf is slightly above 0.5.
        assert!(p > 0.5 && p < 0.55, "p={p}");
    }

    #[test]
    fn gamma_inc_regime_a_very_large() {
        // Pure asymptotic regime, a beyond rmathlib's NaN cliff.
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
    fn try_gamma_overflow_is_err() {
        // Γ(a) overflows f64 for |a| ≥ 1000; reported as GammaDomainError::Overflow.
        assert_eq!(try_gamma(1001.0), Err(GammaDomainError::Overflow(1001.0)));
        assert_eq!(try_gamma(1e10), Err(GammaDomainError::Overflow(1e10)));
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
        // Γ(-3.5) = -8√π/15 via reflection identity.
        // Γ(n+0.5) = (2n)! √π / (4^n n!). For n=3: 6! / (4^3 3!) = 720/384 = 15/8.
        // So Γ(3.5) = (15/8)√π. And Γ(-3.5) = (-1)^4 π / (sin(3.5π) Γ(4.5))
        // ... easier: numerical check vs known-stable computation.
        let g = gamma(-3.5);
        assert!(g.is_finite());
        // The reflection formula: Γ(z)Γ(1-z) = π/sin(πz)
        // → Γ(-3.5) = π / (sin(-3.5π) Γ(4.5))
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
        // a in [15..1000] uses asymptotic.
        let ln_gamma_50 = gamma_log(50.0);
        let g_50 = gamma(50.0);
        // log(g_50) should match ln_gamma_50 to high precision.
        assert!((g_50.ln() - ln_gamma_50).abs() < 1e-9);
        // For very large a (but still < 1000), result should be huge but finite.
        let g_100 = gamma(100.0);
        assert!(g_100.is_finite() && g_100 > 1e150);
    }

    #[test]
    fn gamma_log_at_a_lt_8_branches() {
        // gamma_log has different branches for a in (0.8, 2.25) vs >= 8.
        // Verify continuity at the boundaries.
        for &a in &[0.5_f64, 0.8, 1.0, 1.5, 2.0, 2.25, 5.0, 8.0, 10.0, 50.0] {
            let h = gamma_log(a);
            assert!(h.is_finite(), "a={a}");
        }
    }

    #[test]
    fn psi_at_a_lt_05_and_large() {
        // psi has different branches for x < 0.5, in (0.5, 3], > 3.
        assert!(psi(0.1).is_finite());
        assert!(psi(0.3).is_finite());
        assert!(psi(2.5).is_finite());
        assert!(psi(50.0).is_finite());
        assert!(psi(1000.0).is_finite());
    }
}
