//! Error function and complementary error function (cdflib.f90:9294, :9450).

#![allow(clippy::excessive_precision)]

use super::{exparg, pow2};

// Rational approximation tables. The F90 declares the same digits in both
// error_f (cdflib.f90:9356-9386) and error_fc (cdflib.f90:9503-9534); they are
// shared here.
const A: [f64; 5] = [
    0.771058495001320e-4,
    -0.133733772997339e-2,
    0.323076579225834e-1,
    0.479137145607681e-1,
    0.128379167095513,
];
const B: [f64; 3] = [
    0.301048631703895e-2,
    0.538971687740286e-1,
    0.375795757275549,
];
const C: f64 = 0.564189583547756;
const P: [f64; 8] = [
    -1.36864857382717e-7,
    5.64195517478974e-1,
    7.21175825088309,
    4.31622272220567e1,
    1.52989285046940e2,
    3.39320816734344e2,
    4.51918953711873e2,
    3.00459261020162e2,
];
const Q: [f64; 8] = [
    1.00000000000000,
    1.27827273196294e1,
    7.70001529352295e1,
    2.77585444743988e2,
    6.38980264465631e2,
    9.31354094850610e2,
    7.90950925327898e2,
    3.00459260956983e2,
];
const R: [f64; 5] = [
    2.10144126479064,
    2.62370141675169e1,
    2.13688200555087e1,
    4.65807828718470,
    2.82094791773523e-1,
];
const S: [f64; 4] = [
    9.41537750555460e1,
    1.87114811799590e2,
    9.90191814623914e1,
    1.80124575948747e1,
];

/// Evaluates the error function.
///
/// The function is defined by
///
/// erf(*x*) = (2 / √π) ∫₀ˣ exp(−*t*²) d*t*.
///
/// Properties of the function include: erf(*x*) → −1 as *x* → −∞;
/// erf(0) = 0; erf(0.476936…) = 0.5; erf(*x*) → +1 as *x* → +∞; and
/// ½ (erf(*x* / √2) + 1) = Φ(*x*).
///
/// Since some compilers already supply a routine named `erf`, CDFLIB gives
/// this routine the distinct name `error_f` (cdflib.f90:9294).
///
/// A NaN argument gives 1, as in the F90.
///
/// Reference: Armido DiDinato, Alfred Morris, Algorithm 708: Significant
/// Digit Computation of the Incomplete Beta Function Ratios, ACM Transactions
/// on Mathematical Software, Volume 18, 1992, pages 360-373.
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
    let mut error_f;

    let ax = x.abs();

    if ax <= 0.5 {
        let t = x * x;

        let top = ((((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t + A[4]) + 1.0;

        let bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;
        error_f = ax * (top / bot);
    } else if ax <= 4.0 {
        let top = ((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5]) * ax
            + P[6])
            * ax
            + P[7];

        let bot = ((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5]) * ax
            + Q[6])
            * ax
            + Q[7];

        error_f = 0.5 + (0.5 - (-(x * x)).exp() * top / bot);
    } else if ax < 5.8 {
        let x2 = x * x;
        let t = 1.0 / x2;

        let top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];

        let bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;

        error_f = (C - top / (x2 * bot)) / ax;
        error_f = 0.5 + (0.5 - (-x2).exp() * error_f);
    } else {
        error_f = 1.0;
    }

    if x < 0.0 {
        error_f = -error_f;
    }

    error_f
}

/// Evaluates the complementary error function erfc(*x*) = 1 − erf(*x*).
///
/// This is CDFLIB's `error_fc(ind, x)` (cdflib.f90:9450) with *ind* = 0.
/// It is computed directly, not as 1 − [`error_f`]\(*x*\), so small values
/// keep full relative accuracy. A NaN argument gives NaN.
///
/// Reference: Armido DiDinato, Alfred Morris, Algorithm 708: Significant
/// Digit Computation of the Incomplete Beta Function Ratios, ACM Transactions
/// on Mathematical Software, Volume 18, 1992, pages 360-373.
///
/// # Example
///
/// ```
/// use cdflib::special::error_fc;
///
/// let y = error_fc(2.0);
/// assert!((y - 0.00467773).abs() < 1e-8);
/// ```
///
/// [`error_f`]: crate::special::error_f
#[inline]
pub fn error_fc(x: f64) -> f64 {
    error_fc_ind(0, x)
}

/// Evaluates the scaled complementary error function exp(*x*²) · erfc(*x*).
///
/// This is CDFLIB's `error_fc(ind, x)` (cdflib.f90:9450) with *ind* ≠ 0.
/// It stays finite for large positive *x*, where erfc(*x*) underflows. A
/// NaN argument gives NaN, and −∞ gives +∞.
///
/// Reference: Armido DiDinato, Alfred Morris, Algorithm 708: Significant
/// Digit Computation of the Incomplete Beta Function Ratios, ACM Transactions
/// on Mathematical Software, Volume 18, 1992, pages 360-373.
#[inline]
pub fn error_fc_scaled(x: f64) -> f64 {
    error_fc_ind(1, x)
}

// The body of error_fc (cdflib.f90:9450-9640). If ind is nonzero, the value
// returned has been multiplied by exp(x * x).
#[allow(clippy::assign_op_pattern)]
fn error_fc_ind(ind: i32, x: f64) -> f64 {
    let mut error_fc;

    // Case abs(x) <= 0.5.
    let ax = x.abs();

    if ax <= 0.5 {
        let t = x * x;

        let top = ((((A[0] * t + A[1]) * t + A[2]) * t + A[3]) * t + A[4]) + 1.0;

        let bot = ((B[0] * t + B[1]) * t + B[2]) * t + 1.0;

        error_fc = 0.5 + (0.5 - x * (top / bot));

        if ind != 0 {
            error_fc = t.exp() * error_fc;
        }

        return error_fc;
    }

    // Case 0.5 < abs(x) <= 4.
    if ax <= 4.0 {
        let top = ((((((P[0] * ax + P[1]) * ax + P[2]) * ax + P[3]) * ax + P[4]) * ax + P[5]) * ax
            + P[6])
            * ax
            + P[7];

        let bot = ((((((Q[0] * ax + Q[1]) * ax + Q[2]) * ax + Q[3]) * ax + Q[4]) * ax + Q[5]) * ax
            + Q[6])
            * ax
            + Q[7];

        error_fc = top / bot;
    } else {
        // Case 4 < abs(x).
        if x <= -5.6 {
            if ind == 0 {
                error_fc = 2.0;
            } else {
                error_fc = 2.0 * (x * x).exp();
            }

            return error_fc;
        }

        if ind == 0 {
            if 100.0 < x {
                error_fc = 0.0;
                return error_fc;
            }

            if -exparg(1) < x * x {
                error_fc = 0.0;
                return error_fc;
            }
        }

        let t = pow2(1.0 / x);

        let top = (((R[0] * t + R[1]) * t + R[2]) * t + R[3]) * t + R[4];

        let bot = (((S[0] * t + S[1]) * t + S[2]) * t + S[3]) * t + 1.0;

        error_fc = (C - t * top / bot) / ax;
    }

    // Final assembly.
    if ind != 0 {
        if x < 0.0 {
            error_fc = 2.0 * (x * x).exp() - error_fc;
        }
    } else {
        let w = x * x;
        let t = w;
        let e = w - t;
        error_fc = ((0.5 + (0.5 - e)) * (-t).exp()) * error_fc;

        if x < 0.0 {
            error_fc = 2.0 - error_fc;
        }
    }

    error_fc
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
