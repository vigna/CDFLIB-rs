//! Standard normal cumulative distribution function and its inverse
//! (cdflib.f90:7582, cdflib.f90:8041, cdflib.f90:14233 and cdflib.f90:8584).

#![allow(clippy::excessive_precision)]

use super::eval_pol;
use super::gamma::alnrel;

/// Computes the cumulative normal distribution.
///
/// Evaluates the normal distribution function
///
/// *P*(*x*) = (1 / √(2π)) ∫ exp(−*t*²/2) d*t*,
///
/// the integral running from −∞ to *x* = *arg*, the upper limit of
/// integration, and returns (*cum*, *ccum*), the normal CDF and the
/// complementary CDF. Both are returned because the smaller one is computed
/// directly, which preserves precision in either tail.
///
/// This transportable program uses rational functions that theoretically
/// approximate the normal distribution function to at least 18 significant
/// decimal digits. The accuracy achieved depends on the arithmetic system,
/// the compiler, the intrinsic functions, and proper selection of the machine
/// dependent constants.
///
/// This is CDFLIB's `cumnor` (cdflib.f90:7582).
///
/// As in the F90, the result is (NaN, NaN) when 16·|*x*| overflows
/// (|*x*| > `f64::MAX` / 16, including *x* = ±∞), and when *x* is NaN.
/// The distributions return the exact limits there.
///
/// References: William Cody, Rational Chebyshev approximations for the error
/// function, Mathematics of Computation, 1969, pages 631-637. William Cody,
/// Algorithm 715: SPECFUN - A Portable Fortran Package of Special Function
/// Routines and Test Drivers, ACM Transactions on Mathematical Software,
/// Volume 19, Number 1, 1993, pages 22-32.
///
/// # Example
///
/// ```
/// use cdflib::special::cumnor;
///
/// let (p, q) = cumnor(1.96);
/// assert!((p - 0.975).abs() < 1e-3);
/// assert!((q - 0.025).abs() < 1e-3);
/// ```
#[inline]
#[allow(unused_assignments, clippy::assign_op_pattern)]
pub fn cumnor(arg: f64) -> (f64, f64) {
    const A: [f64; 5] = [
        2.2352520354606839287,
        1.6102823106855587881e2,
        1.0676894854603709582e3,
        1.8154981253343561249e4,
        6.5682337918207449113e-2,
    ];
    const B: [f64; 4] = [
        4.7202581904688241870e1,
        9.7609855173777669322e2,
        1.0260932208618978205e4,
        4.5507789335026729956e4,
    ];
    const C: [f64; 9] = [
        3.9894151208813466764e-1,
        8.8831497943883759412,
        9.3506656132177855979e1,
        5.9727027639480026226e2,
        2.4945375852903726711e3,
        6.8481904505362823326e3,
        1.1602651437647350124e4,
        9.8427148383839780218e3,
        1.0765576773720192317e-8,
    ];
    const D: [f64; 8] = [
        2.2266688044328115691e1,
        2.3538790178262499861e2,
        1.5193775994075548050e3,
        6.4855582982667607550e3,
        1.8615571640885098091e4,
        3.4900952721145977266e4,
        3.8912003286093271411e4,
        1.9685429676859990727e4,
    ];
    const P: [f64; 6] = [
        2.1589853405795699e-1,
        1.274011611602473639e-1,
        2.2235277870649807e-2,
        1.421619193227893466e-3,
        2.9112874951168792e-5,
        2.307344176494017303e-2,
    ];
    const Q: [f64; 5] = [
        1.28426009614491121,
        4.68238212480865118e-1,
        6.59881378689285515e-2,
        3.78239633202758244e-3,
        7.29751555083966205e-5,
    ];
    const ROOT32: f64 = 5.656854248;
    const SIXTEN: f64 = 16.0;
    const SQRPI: f64 = 3.9894228040143267794e-1;
    const THRSH: f64 = 0.66291;

    let mut cum;
    let mut ccum;
    let mut xsq;

    // Machine dependent constants: eps is epsilon(1.0D+00) * 0.5D+00
    // (cdflib.f90:7713).
    let eps = f64::EPSILON * 0.5;

    let x = arg;
    let y = x.abs();

    if y <= THRSH {
        // Evaluate anorm for abs(x) <= 0.66291.
        if eps < y {
            xsq = x * x;
        } else {
            xsq = 0.0;
        }

        let mut xnum = A[4] * xsq;
        let mut xden = xsq;
        for i in 0..3 {
            xnum = (xnum + A[i]) * xsq;
            xden = (xden + B[i]) * xsq;
        }
        cum = x * (xnum + A[3]) / (xden + B[3]);
        let temp = cum;
        cum = 0.5 + temp;
        ccum = 0.5 - temp;
    } else if y <= ROOT32 {
        // Evaluate anorm for 0.66291 <= abs(x) <= sqrt(32).
        let mut xnum = C[8] * y;
        let mut xden = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        cum = (xnum + C[7]) / (xden + D[7]);
        xsq = (y * SIXTEN).trunc() / SIXTEN;
        let del = (y - xsq) * (y + xsq);
        cum = (-(xsq * xsq * 0.5)).exp() * (-(del * 0.5)).exp() * cum;
        ccum = 1.0 - cum;

        if 0.0 < x {
            std::mem::swap(&mut cum, &mut ccum);
        }
    } else {
        // Evaluate anorm for sqrt(32) < abs(x). cdflib.f90:7763 has the
        // dead store cum = 0.
        cum = 0.0;
        xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }

        cum = xsq * (xnum + P[4]) / (xden + Q[4]);
        cum = (SQRPI - cum) / y;
        xsq = (x * SIXTEN).trunc() / SIXTEN;
        let del = (x - xsq) * (x + xsq);
        cum = (-(xsq * xsq * 0.5)).exp() * (-(del * 0.5)).exp() * cum;
        ccum = 1.0 - cum;

        if 0.0 < x {
            std::mem::swap(&mut cum, &mut ccum);
        }
    }

    // The threshold is tiny(cum) (cdflib.f90:7786 and cdflib.f90:7790).
    if cum < f64::MIN_POSITIVE {
        cum = 0.0;
    }

    if ccum < f64::MIN_POSITIVE {
        ccum = 0.0;
    }

    (cum, ccum)
}

/// Computes the inverse of the normal distribution.
///
/// Returns *x* such that [`cumnor`]\(*x*\) = *p*, that is, so that
///
/// *p* = ∫ exp(−*u*²/2) / √(2π) d*u*,
///
/// the integral running from −∞ to *x*. The arguments *p* and *q* are the
/// probability and the complementary probability; the search runs on the
/// smaller of the two, which preserves precision when *p* is close to 1.
///
/// The rational function on page 95 of Kennedy and Gentle ([`stvaln`]) is
/// used as a starting value for the Newton method of finding roots. If the
/// Newton method does not converge in 100 iterations, the starting value is
/// returned.
///
/// This is CDFLIB's `dinvnr` (cdflib.f90:8041).
///
/// As in the F90, the result is NaN when *p* or *q* is 0, where the
/// answer is ∓∞, and when both are NaN; the distributions handle these
/// endpoints separately.
///
/// Reference: William Kennedy, James Gentle, Statistical Computing, Marcel
/// Dekker, NY, 1980.
///
/// # Example
///
/// ```
/// use cdflib::special::dinvnr;
///
/// let x = dinvnr(0.975, 0.025);
/// assert!((x - 1.95996).abs() < 1e-4);
/// ```
///
/// [`cumnor`]: crate::special::cumnor
/// [`stvaln`]: crate::special::internal::stvaln
#[inline]
pub fn dinvnr(p: f64, q: f64) -> f64 {
    const EPS: f64 = 1.0e-13;
    const MAXIT: i32 = 100;
    const R2PI: f64 = 0.3989422804014326;

    // As with gfortran's min, f64::min returns the other argument when one
    // of p and q is NaN.
    let pp = p.min(q);
    let strtx = stvaln(pp);
    let mut xcur = strtx;

    // Newton iterations.
    for _i in 1..=MAXIT {
        let (cum, _ccum) = cumnor(xcur);
        let dx = (cum - pp) / (R2PI * (-(0.5 * xcur * xcur)).exp());
        xcur -= dx;

        if (dx / xcur).abs() < EPS {
            return if p <= q { xcur } else { -xcur };
        }
    }

    if p <= q {
        strtx
    } else {
        -strtx
    }
}

/// Provides starting values for the inverse of the normal distribution.
///
/// Returns an *x* for which it is approximately true that
/// *p* = [`cumnor`]\(*x*\), that is,
///
/// *p* = ∫ exp(−*u*²/2) / √(2π) d*u*,
///
/// the integral running from −∞ to *x*; *p* is the probability whose normal
/// deviate is sought. This is the starting value of [`dinvnr`].
///
/// This is CDFLIB's `stvaln` (cdflib.f90:14233).
///
/// As in the F90, the result is NaN when *p* is 0 or 1, and when *p* is
/// NaN.
///
/// Reference: William Kennedy, James Gentle, Statistical Computing, Marcel
/// Dekker, NY, 1980, page 95.
///
/// [`cumnor`]: crate::special::cumnor
/// [`dinvnr`]: crate::special::dinvnr
#[allow(clippy::assign_op_pattern, clippy::needless_late_init)]
pub fn stvaln(p: f64) -> f64 {
    const XDEN: [f64; 5] = [
        0.993484626060e-1,
        0.588581570495,
        0.531103462366,
        0.103537752850,
        0.38560700634e-2,
    ];
    const XNUM: [f64; 5] = [
        -0.322232431088,
        -1.000000000000,
        -0.342242088547,
        -0.204231210245e-1,
        -0.453642210148e-4,
    ];

    let sgn;
    let z;

    if p <= 0.5 {
        sgn = -1.0;
        z = p;
    } else {
        sgn = 1.0;
        z = 1.0 - p;
    }

    let y = (-2.0 * z.ln()).sqrt();
    let mut stvaln = y + eval_pol(&XNUM, y) / eval_pol(&XDEN, y);
    stvaln = sgn * stvaln;

    stvaln
}

/// Evaluates the logarithm of the asymptotic normal CDF.
///
/// Computes the logarithm of the cumulative normal distribution from |*x*|
/// to infinity, that is, ln Pr[*X* > |*x*|] for a standard normal *X*, for
/// 5 ≤ |*x*|.
///
/// The relative error at *x* = 5 is about 0.5·10⁻⁵.
///
/// This is CDFLIB's `dlanor` (cdflib.f90:8584).
///
/// Reference: Milton Abramowitz, Irene Stegun, Handbook of Mathematical
/// Functions, 1966, Formula 26.2.12.
///
/// # Panics
///
/// Panics if |*x*| < 5. In this case the F90 routine prints a fatal-error
/// message and then continues with the asymptotic formula anyway.
///
/// # Example
///
/// ```
/// use cdflib::special::{cumnor, dlanor};
///
/// // At x = 8 the complementary probability is around 6.22e-16.
/// // dlanor returns its log; exp(dlanor) should match cumnor's ccum.
/// let log_q = dlanor(8.0);
/// let (_, q) = cumnor(8.0);
/// assert!((log_q.exp() / q - 1.0).abs() < 1e-5);
/// ```
#[inline]
pub fn dlanor(x: f64) -> f64 {
    // The coefficients (-1)^(k+1) (2k+1)!! of the asymptotic series in
    // 1 / x^2.
    const COEF: [f64; 12] = [
        -1.0,
        3.0,
        -15.0,
        105.0,
        -945.0,
        10395.0,
        -135135.0,
        2027025.0,
        -34459425.0,
        654729075.0,
        -13749310575.0,
        316234143225.0,
    ];
    const DLSQPI: f64 = 0.91893853320467274177;

    let xx = x.abs();

    if x.abs() < 5.0 {
        // Rust only: panic where the F90 prints a fatal-error message and
        // continues (cdflib.f90:8641-8645).
        panic!("dlanor: argument |x| must be ≥ 5 (got {x})");
    }

    let approx = -DLSQPI - 0.5 * x * x - x.abs().ln();

    let xx2 = xx * xx;
    let mut correc = eval_pol(&COEF, 1.0 / xx2) / xx2;
    correc = alnrel(correc);

    approx + correc
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================================ dlanor

    #[test]
    fn dlanor_matches_log_cumnor_tail() {
        // For each x ≥ 5, exp(dlanor(x)) should equal cumnor(x)'s ccum
        // to a few digits (the asymptotic formula's stated accuracy at
        // x = 5 is about 5·10⁻⁶, improving as |x| grows).
        for &x in &[5.0_f64, 6.0, 7.0, 8.0, 10.0, 15.0] {
            let log_q = dlanor(x);
            let (_, q) = cumnor(x);
            let rel = (log_q.exp() - q).abs() / q;
            assert!(
                rel < 1e-5,
                "x={x}: dlanor.exp={}, ccum={q}, rel={rel}",
                log_q.exp()
            );
        }
    }

    #[test]
    fn dlanor_symmetric_in_magnitude() {
        // dlanor depends on |x| only; sign should not change the result.
        for &x in &[5.5_f64, 8.0, 12.0] {
            assert_eq!(dlanor(x), dlanor(-x));
        }
    }

    #[test]
    fn dlanor_decreasing_in_x() {
        // Pr[X > x] decreases as x grows, so log Pr[...] also decreases.
        let a = dlanor(5.0);
        let b = dlanor(10.0);
        let c = dlanor(20.0);
        assert!(a > b && b > c);
    }

    #[test]
    #[should_panic(expected = "argument |x| must be ≥ 5")]
    fn dlanor_panics_below_threshold() {
        let _ = dlanor(4.99);
    }

    // ============================================================ cumnor / dinvnr

    #[test]
    fn cumnor_at_zero() {
        let (p, q) = cumnor(0.0);
        assert!((p - 0.5).abs() < 1e-15, "p = {p}");
        assert!((q - 0.5).abs() < 1e-15, "q = {q}");
    }

    #[test]
    fn cumnor_at_one_sigma() {
        let (p, q) = cumnor(1.0);
        // Reference: Φ(1) ≈ 0.8413447460685429
        assert!((p - 0.8413447460685429).abs() < 1e-14, "p = {p}");
        assert!((p + q - 1.0).abs() < 1e-15);
    }

    #[test]
    fn cumnor_symmetry() {
        for &x in &[0.1, 1.0, 2.5, 4.0, 8.0] {
            let (p_pos, q_pos) = cumnor(x);
            let (p_neg, q_neg) = cumnor(-x);
            assert!((p_pos - q_neg).abs() < 1e-15, "x = {x}");
            assert!((q_pos - p_neg).abs() < 1e-15, "x = {x}");
        }
    }

    #[test]
    fn cumnor_tail_accuracy() {
        // Φ(-10) ≈ 7.62e-24. A naive 1-Φ(10) would underflow to 0.
        let (_p, q) = cumnor(10.0);
        assert!(q > 0.0 && q < 1e-22, "q = {q}");
    }

    #[test]
    fn dinvnr_roundtrip() {
        for &x in &[-3.0, -1.0, -0.1, 0.5, 2.0, 4.0] {
            let (p, q) = cumnor(x);
            let back = dinvnr(p, q);
            assert!((back - x).abs() < 1e-12, "x = {x}, back = {back}");
        }
    }

    #[test]
    fn dinvnr_tail_accuracy() {
        // dinvnr should hit ~5.0 even when p is essentially 1.0 because we
        // route through the small tail q.
        let (_p, q) = cumnor(5.0);
        let back = dinvnr(1.0 - q, q);
        assert!((back - 5.0).abs() < 1e-9, "back = {back}");
    }
}
