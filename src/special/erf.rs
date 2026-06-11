//! Error function and complementary error function.

#![allow(clippy::excessive_precision)]

/// Largest negative argument to `exp` for which the result is nonzero in
/// IEEE 754 binary64.
///
/// Matches F90's `exparg(1)` exactly: `0.99999 × (-1022) × 0.69314718055995`,
/// where -1022 = `ipmpar(9) - 1` for IEEE binary64 and `0.69314718055995` is
/// the literal `lnb` used by F90 for base b=2 (cdflib.f90:9544).
const NEG_EXPARG: f64 = -708.389_334_568_083_540_9;

// Coefficients for |x| ≤ 0.5.
const A: [f64; 5] = [
    7.71058495001320e-05,
    -1.33733772997339e-03,
    3.23076579225834e-02,
    4.79137145607681e-02,
    1.28379167095513e-01,
];
const B: [f64; 3] = [
    3.01048631703895e-03,
    5.38971687740286e-02,
    3.75795757275549e-01,
];

// Coefficients for 0.5 < |x| ≤ 4.
const P: [f64; 8] = [
    -1.36864857382717e-07,
    5.64195517478974e-01,
    7.21175825088309e+00,
    4.31622272220567e+01,
    1.52989285046940e+02,
    3.39320816734344e+02,
    4.51918953711873e+02,
    3.00459261020162e+02,
];
const Q: [f64; 8] = [
    1.00000000000000e+00,
    1.27827273196294e+01,
    7.70001529352295e+01,
    2.77585444743988e+02,
    6.38980264465631e+02,
    9.31354094850610e+02,
    7.90950925327898e+02,
    3.00459260956983e+02,
];

// Coefficients for |x| > 4.
const R: [f64; 5] = [
    2.10144126479064e+00,
    2.62370141675169e+01,
    2.13688200555087e+01,
    4.65807828718470e+00,
    2.82094791773523e-01,
];
const S: [f64; 4] = [
    9.41537750555460e+01,
    1.87114811799590e+02,
    9.90191814623914e+01,
    1.80124575948747e+01,
];

const C: f64 = 0.564189583547756;

/// Returns the error function erf(*x*) = (2/√π) ∫₀ˣ e⁻*ᵗ*² d*t*.
///
/// # Example
///
/// ```
/// use cdflib::special::error_f;
///
/// let y = error_f(0.8);
/// assert!((y - 0.74210096).abs() < 1e-8);
/// ```
#[inline]
pub fn error_f(x: f64) -> f64 {
    let ax = x.abs();

    let mut erf = if ax <= 0.5 {
        let t = x * x;
        let top = ((((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t) + A[4] + 1.0;
        let bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;
        ax * (top / bot)
    } else if ax <= 4.0 {
        let top = (((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5])
            * ax
            + P[6])
            * ax)
            + P[7];
        let bot = (((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5])
            * ax
            + Q[6])
            * ax)
            + Q[7];
        // erf = 1 - exp(-x²) * top/bot ; written as 0.5 + (0.5 - …) for
        // tail-precision-friendly assembly, matching CDFLIB.
        0.5 + (0.5 - (-(x * x)).exp() * top / bot)
    } else if ax < 5.8 {
        let x2 = x * x;
        let t = 1.0 / x2;
        let top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];
        let bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;
        let erf = (C - top / (x2 * bot)) / ax;
        0.5 + (0.5 - (-x2).exp() * erf)
    } else {
        // |x| >= 5.8: erf saturates to 1 before the sign fixup.
        1.0
    };
    if x < 0.0 {
        erf = -erf;
    }
    erf
}

/// Returns the complementary error function erfc(*x*) = 1 − erf(*x*).
///
/// Computed directly (not as `1 - error_f(x)`) so the small right-tail
/// values stay accurate to ~15 digits. CDFLIB exposes this and the scaled
/// variant [`error_fc_scaled`] via a single `int *ind` flag; we split them
/// into two Rust functions for clarity.
///
/// [`error_fc_scaled`]: crate::special::error_fc_scaled
///
/// # Example
///
/// ```
/// use cdflib::special::error_fc;
///
/// let y = error_fc(2.0);
/// assert!((y - 0.00467773).abs() < 1e-8);
/// ```
#[inline]
pub fn error_fc(x: f64) -> f64 {
    error_fc_inner(x, false)
}

/// Returns erfc(*x*) · exp(*x*²). Useful for very large |*x*| where erfc itself
/// underflows but its exponentially-scaled form does not.
#[inline]
pub fn error_fc_scaled(x: f64) -> f64 {
    error_fc_inner(x, true)
}

fn error_fc_inner(x: f64, scaled: bool) -> f64 {
    let ax = x.abs();

    // |x| ≤ 0.5
    if ax <= 0.5 {
        let t = x * x;
        let top = ((((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t) + A[4] + 1.0;
        let bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;
        let mut erfc = 0.5 + (0.5 - x * (top / bot));
        if scaled {
            erfc *= t.exp();
        }
        return erfc;
    }

    // 0.5 < |x| ≤ 4
    let mut erfc;
    if ax <= 4.0 {
        let top = (((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5])
            * ax
            + P[6])
            * ax)
            + P[7];
        let bot = (((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5])
            * ax
            + Q[6])
            * ax)
            + Q[7];
        erfc = top / bot;
    } else {
        // |x| > 4
        // Large-negative-x cutoff: erfc(x) → 2 as x → -∞.
        if x <= -5.6 {
            return if scaled { 2.0 * (x * x).exp() } else { 2.0 };
        }
        // For the unscaled form, also check the overflow boundary on the
        // positive side: when -x² ≤ NEG_EXPARG, exp(-x²) underflows.
        if !scaled && (x > 100.0 || x * x > -NEG_EXPARG) {
            return 0.0;
        }

        let t = (1.0 / x).powi(2);
        let top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];
        let bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;
        erfc = (C - t * top / bot) / ax;
    }

    // Final assembly.
    if scaled {
        if x < 0.0 {
            erfc = 2.0 * (x * x).exp() - erfc;
        }
    } else {
        erfc *= (-(x * x)).exp();
        if x < 0.0 {
            erfc = 2.0 - erfc;
        }
    }
    erfc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_zero() {
        assert_eq!(error_f(0.0), 0.0);
        assert_eq!(error_fc(0.0), 1.0);
    }

    #[test]
    fn erf_is_odd() {
        for &x in &[0.1, 0.7, 2.0, 4.5] {
            let a = error_f(x);
            let b = error_f(-x);
            assert!((a + b).abs() < 1e-15, "erf({x}) = {a}, erf(-{x}) = {b}");
        }
    }

    #[test]
    fn erf_saturates() {
        assert_eq!(error_f(10.0), 1.0);
        assert_eq!(error_f(-10.0), -1.0);
    }

    #[test]
    fn erfc_complement_relation() {
        for &x in &[-2.0, -0.5, 0.0, 0.3, 1.5, 3.7] {
            let s = error_f(x) + error_fc(x);
            assert!((s - 1.0).abs() < 1e-14, "erf({x}) + erfc({x}) = {s}");
        }
    }
}
