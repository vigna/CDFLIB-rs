//! Β function family: Β(*a*, *b*), ln Β(*a*, *b*), and the regularized
//! incomplete Β function *Iₓ*(*a*, *b*) and its complement.

#![allow(clippy::approx_constant, clippy::excessive_precision)]

use super::erf::error_fc_scaled;
use super::gamma::{alnrel, gam1, gamma_ln1, gamma_log, gsumln, psi, rexp, rlog1};
use super::pow2;

/// Evaluates exp(*mu* + *x*) (cdflib.f90:9641).
///
/// The integer *mu* and the real *x* are both part of the argument.
#[inline]
#[allow(clippy::collapsible_if)]
pub fn esum(mu: i32, x: f64) -> f64 {
    if x <= 0.0 {
        if 0 <= mu {
            let w = mu as f64 + x;
            if w <= 0.0 {
                return w.exp();
            }
        }
    } else if 0.0 < x {
        if mu <= 0 {
            let w = mu as f64 + x;
            if 0.0 <= w {
                return w.exp();
            }
        }
    }

    let w = mu as f64;
    w.exp() * x.exp()
}

/// Computes ln(Γ(*b*) / Γ(*a* + *b*)) when 8 ≤ *b* (cdflib.f90:1).
///
/// In this algorithm, DEL(*x*) is the function defined by
/// ln Γ(*x*) = (*x* − 0.5) ln *x* − *x* + 0.5 ln(2π) + DEL(*x*).
#[inline]
#[allow(clippy::needless_late_init)]
pub fn algdiv(a: f64, b: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-1;
    const C1: f64 = -0.277777777760991e-2;
    const C2: f64 = 0.793650666825390e-3;
    const C3: f64 = -0.595202931351870e-3;
    const C4: f64 = 0.837308034031215e-3;
    const C5: f64 = -0.165322962780713e-2;

    let c;
    let d;
    let h;
    let x;

    if b < a {
        h = b / a;
        c = 1.0 / (1.0 + h);
        x = h / (1.0 + h);
        d = a + (b - 0.5);
    } else {
        h = a / b;
        c = h / (1.0 + h);
        x = 1.0 / (1.0 + h);
        d = b + (a - 0.5);
    }

    // Set SN = (1 - X^N)/(1 - X).
    let x2 = x * x;
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);

    // Set W = DEL(B) - DEL(A + B).
    let t = pow2(1.0 / b);
    let mut w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;

    w *= c / b;

    // Combine the results.
    let u = d * alnrel(a / b);
    let v = a * (b.ln() - 1.0);

    if v < u {
        (w - v) - u
    } else {
        (w - u) - v
    }
}

/// Evaluates the logarithm of the Β function, ln Β(*a0*, *b0*)
/// (cdflib.f90:1460).
///
/// *a0* and *b0* should be nonnegative.
///
/// # Example
///
/// ```
/// use cdflib::special::beta_log;
///
/// let y = beta_log(3.0, 4.0);
/// // Β(3, 4) = 1/60
/// assert!((y - (1.0/60.0_f64).ln()).abs() < 1e-14);
/// ```
#[inline]
pub fn beta_log(a0: f64, b0: f64) -> f64 {
    const E: f64 = 0.918938533204673;

    let mut a = a0.min(b0);
    let mut b = a0.max(b0);

    // 8 <= A.
    if 8.0 <= a {
        let w = bcorr(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -((a - 0.5) * c.ln());
        let v = b * alnrel(h);

        return if v < u {
            (((-(0.5 * b.ln()) + E) + w) - v) - u
        } else {
            (((-(0.5 * b.ln()) + E) + w) - u) - v
        };
    }

    // Procedure when A < 1.
    if a < 1.0 {
        return if b < 8.0 {
            gamma_log(a) + (gamma_log(b) - gamma_log(a + b))
        } else {
            gamma_log(a) + algdiv(a, b)
        };
    }

    // Procedure when 1 <= A < 8.
    let mut w;
    'l60: {
        'l40: {
            if 2.0 < a {
                break 'l40;
            }

            if b <= 2.0 {
                return gamma_log(a) + gamma_log(b) - gsumln(a, b);
            }

            w = 0.0;

            if b < 8.0 {
                break 'l60;
            }

            return gamma_log(a) + algdiv(a, b);
        }

        // Label 40: reduction of A when 1000 < B.
        if 1000.0 < b {
            let n = (a - 1.0) as i32;
            w = 1.0;
            for _ in 1..=n {
                a -= 1.0;
                w *= a / (1.0 + a / b);
            }

            return (w.ln() - n as f64 * b.ln()) + (gamma_log(a) + algdiv(a, b));
        }

        let n = (a - 1.0) as i32;
        w = 1.0;
        for _ in 1..=n {
            a -= 1.0;
            let h = a / b;
            w *= h / (1.0 + h);
        }
        w = w.ln();

        if 8.0 <= b {
            return w + gamma_log(a) + algdiv(a, b);
        }
    }

    // Label 60: reduction of B when B < 8.
    let n = (b - 1.0) as i32;
    let mut z = 1.0;
    for _ in 1..=n {
        b -= 1.0;
        z *= b / (a + b);
    }

    w + z.ln() + (gamma_log(a) + (gamma_log(b) - gsumln(a, b)))
}

/// Evaluates DEL(*a0*) + DEL(*b0*) − DEL(*a0* + *b0*) (cdflib.f90:282).
///
/// The function DEL(*a*) is a remainder term that is used in the expression
/// ln Γ(*a*) = (*a* − 0.5) ln *a* − *a* + 0.5 ln(2π) + DEL(*a*).
///
/// It is assumed that 8 ≤ *a0* and 8 ≤ *b0*.
#[inline]
pub fn bcorr(a0: f64, b0: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-1;
    const C1: f64 = -0.277777777760991e-2;
    const C2: f64 = 0.793650666825390e-3;
    const C3: f64 = -0.595202931351870e-3;
    const C4: f64 = 0.837308034031215e-3;
    const C5: f64 = -0.165322962780713e-2;

    let a = a0.min(b0);
    let b = a0.max(b0);

    let h = a / b;
    let c = h / (1.0 + h);
    let x = 1.0 / (1.0 + h);
    let x2 = x * x;

    // Set SN = (1 - X^N)/(1 - X).
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);

    // Set W = DEL(B) - DEL(A + B).
    let mut t = pow2(1.0 / b);

    let mut w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;

    w *= c / b;

    // Compute DEL(A) + W.
    t = pow2(1.0 / a);

    (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a + w
}

/// Evaluates the Β function, Β(*a*, *b*) (cdflib.f90:400).
///
/// # Example
///
/// ```
/// use cdflib::special::beta;
///
/// let y = beta(3.0, 4.0);
/// assert!((y - 1.0/60.0).abs() < 1e-14);
/// ```
#[inline]
pub fn beta(a: f64, b: f64) -> f64 {
    beta_log(a, b).exp()
}

/// Computes the Sterling remainder for the complete Β function
/// (cdflib.f90:7919).
///
/// ln Β(*a*, *b*) = ln Γ(*a*) + ln Γ(*b*) − ln Γ(*a* + *b*). Let *zz* be the
/// approximation obtained if each ln Γ is approximated by Sterling's formula,
/// Sterling(*z*) = ln √(2π) + (*z* − 0.5) ln *z* − *z*. The Sterling remainder
/// is ln Β(*a*, *b*) − *zz*.
///
/// # Example
///
/// ```
/// use cdflib::special::internal::dbetrm;
///
/// // Sterling remainder is small and decreasing in (a, b) for large args.
/// let r = dbetrm(50.0, 60.0);
/// assert!(r.abs() < 0.01);
/// ```
#[inline]
pub fn dbetrm(a: f64, b: f64) -> f64 {
    use super::gamma::dstrem;

    // Try to sum from smallest to largest.
    let mut dbetrm = -dstrem(a + b);
    dbetrm += dstrem(a.max(b));
    dbetrm += dstrem(a.min(b));

    dbetrm
}

/// Evaluates *Iₓ*(*a*, *b*) for very small *b* (cdflib.f90:10075).
///
/// This routine is appropriate for use when *b* < min(*eps*, *eps*·*a*) and
/// *x* ≤ 0.5.
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn fpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    use super::exparg;

    let mut fpser = 1.0;

    if 1.0e-3 * eps < a {
        fpser = 0.0;
        let t = a * x.ln();
        if t < exparg(1) {
            return fpser;
        }
        fpser = t.exp();
    }

    // 1/B(A,B) = B
    fpser = (b / a) * fpser;
    let tol = eps / a;
    let mut an = a + 1.0;
    let mut t = x;
    let mut s = t / an;

    loop {
        an += 1.0;
        t = x * t;
        let c = t / an;
        s += c;

        if c.abs() <= tol {
            break;
        }
        // Rust only: the F90 loop never exits once tol is NaN or a term is
        // NaN or infinite.
        if tol.is_nan() || !c.is_finite() {
            return f64::NAN;
        }
    }

    fpser * (1.0 + a * s)
}

/// Computes the incomplete Β ratio *I*₁₋ₓ(*b*, *a*) (cdflib.f90:188).
///
/// `apser` is used only for cases where *a* ≤ min(*eps*, *eps*·*b*),
/// *b*·*x* ≤ 1, and *x* ≤ 0.5.
#[inline]
pub fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    const G: f64 = 0.577215664901533;

    let bx = b * x;
    let mut t = x - bx;

    let c = if b * eps <= 0.02 {
        x.ln() + psi(b) + G + t
    } else {
        bx.ln() + G + t
    };

    let tol = 5.0 * eps * c.abs();
    let mut j = 1.0;
    let mut s = 0.0;

    loop {
        j += 1.0;
        t *= x - bx / j;
        let aj = t / j;
        s += aj;

        if aj.abs() <= tol {
            break;
        }
        // Rust only: the F90 loop never exits once tol is NaN or a term is
        // NaN or infinite.
        if tol.is_nan() || !aj.is_finite() {
            return f64::NAN;
        }
    }

    -(a * (c + s))
}

/// Uses a power series expansion to evaluate *Iₓ*(*a*, *b*)
/// (cdflib.f90:1625).
///
/// `beta_pser` is used when *b* ≤ 1 or *b*·*x* ≤ 0.7. *eps* is the tolerance.
#[inline]
#[allow(clippy::assign_op_pattern)]
pub fn beta_pser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let mut beta_pser = 0.0;

    if x == 0.0 {
        return beta_pser;
    }

    // Compute the factor X^A/(A*BETA(A,B)).
    let a0 = a.min(b);

    if 1.0 <= a0 {
        let z = a * x.ln() - beta_log(a, b);
        beta_pser = z.exp() / a;
    } else {
        let mut b0 = a.max(b);

        if b0 <= 1.0 {
            beta_pser = x.powf(a);
            if beta_pser == 0.0 {
                return beta_pser;
            }

            let apb = a + b;

            let z = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                let u = a + b - 1.0;
                (1.0 + gam1(u)) / apb
            };

            let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
            beta_pser = beta_pser * c * (b / apb);
        } else if b0 < 8.0 {
            let mut u = gamma_ln1(a0);
            let m = (b0 - 1.0) as i32;

            let mut c = 1.0;
            for _ in 1..=m {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }

            u = c.ln() + u;
            let z = a * x.ln() - u;
            b0 -= 1.0;
            let apb = a0 + b0;

            let t = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                u = a0 + b0 - 1.0;
                (1.0 + gam1(u)) / apb
            };

            beta_pser = z.exp() * (a0 / a) * (1.0 + gam1(b0)) / t;
        } else if 8.0 <= b0 {
            let u = gamma_ln1(a0) + algdiv(a0, b0);
            let z = a * x.ln() - u;
            beta_pser = (a0 / a) * z.exp();
        }
    }

    if beta_pser == 0.0 || a <= 0.1 * eps {
        return beta_pser;
    }

    // Compute the series.
    let mut sum1 = 0.0;
    let mut n = 0.0;
    let mut c = 1.0;
    let tol = eps / a;

    loop {
        n += 1.0;
        c = c * (0.5 + (0.5 - b / n)) * x;
        let w = c / (a + n);
        sum1 += w;

        if w.abs() <= tol {
            break;
        }
        // Rust only: the F90 loop never exits once tol is NaN or a term is
        // NaN or infinite.
        if tol.is_nan() || !w.is_finite() {
            return f64::NAN;
        }
    }

    beta_pser * (1.0 + a * sum1)
}

/// Evaluates *xᵃ* · *yᵇ* / Β(*a*, *b*) (cdflib.f90:1798).
///
/// *a* and *b* should be nonnegative; *x* and *y* define the numerator of the
/// fraction.
#[inline]
#[allow(clippy::assign_op_pattern, clippy::needless_late_init)]
pub fn beta_rcomp(a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST: f64 = 0.398942280401433;

    let mut beta_rcomp = 0.0;
    if x == 0.0 || y == 0.0 {
        return beta_rcomp;
    }

    let a0 = a.min(b);

    if a0 < 8.0 {
        let lnx;
        let lny;
        if x <= 0.375 {
            lnx = x.ln();
            lny = alnrel(-x);
        } else if y <= 0.375 {
            lnx = alnrel(-y);
            lny = y.ln();
        } else {
            lnx = x.ln();
            lny = y.ln();
        }

        let mut z = a * lnx + b * lny;

        if 1.0 <= a0 {
            z -= beta_log(a, b);
            beta_rcomp = z.exp();
            return beta_rcomp;
        }

        // Procedure for A < 1 or B < 1.
        let mut b0 = a.max(b);

        if b0 <= 1.0 {
            beta_rcomp = z.exp();
            if beta_rcomp == 0.0 {
                return beta_rcomp;
            }

            let apb = a + b;

            if apb <= 1.0 {
                z = 1.0 + gam1(apb);
            } else {
                let u = a + b - 1.0;
                z = (1.0 + gam1(u)) / apb;
            }

            let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
            beta_rcomp = beta_rcomp * (a0 * c) / (1.0 + a0 / b0);
        } else if b0 < 8.0 {
            let mut u = gamma_ln1(a0);
            let n = (b0 - 1.0) as i32;

            let mut c = 1.0;
            for _ in 1..=n {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u = c.ln() + u;

            z -= u;
            b0 -= 1.0;
            let apb = a0 + b0;

            let t = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                u = a0 + b0 - 1.0;
                (1.0 + gam1(u)) / apb
            };

            beta_rcomp = a0 * z.exp() * (1.0 + gam1(b0)) / t;
        } else if 8.0 <= b0 {
            let u = gamma_ln1(a0) + algdiv(a0, b0);
            beta_rcomp = a0 * (z - u).exp();
        }
    } else {
        let h;
        let x0;
        let y0;
        let lambda;
        if a <= b {
            h = a / b;
            x0 = h / (1.0 + h);
            y0 = 1.0 / (1.0 + h);
            lambda = a - (a + b) * x;
        } else {
            h = b / a;
            x0 = 1.0 / (1.0 + h);
            y0 = h / (1.0 + h);
            lambda = (a + b) * y - b;
        }

        let mut e = -(lambda / a);

        let u = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (x / x0).ln()
        };

        e = lambda / b;

        let v = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (y / y0).ln()
        };

        let z = (-(a * u + b * v)).exp();
        beta_rcomp = CONST * (b * x0).sqrt() * z * (-bcorr(a, b)).exp();
    }

    beta_rcomp
}

/// Evaluates exp(*mu*) · *xᵃ* · *yᵇ* / Β(*a*, *b*) (cdflib.f90:1992).
///
/// *a* and *b* should be nonnegative; *x* and *y* are the quantities whose
/// powers form part of the expression.
#[inline]
#[allow(clippy::assign_op_pattern, clippy::needless_late_init)]
pub fn beta_rcomp1(mu: i32, a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST: f64 = 0.398942280401433;

    let mut beta_rcomp1;

    let a0 = a.min(b);

    if 8.0 <= a0 {
        // Procedure for 8 <= A and 8 <= B.
        let h;
        let x0;
        let y0;
        let lambda;
        if a <= b {
            h = a / b;
            x0 = h / (1.0 + h);
            y0 = 1.0 / (1.0 + h);
            lambda = a - (a + b) * x;
        } else {
            h = b / a;
            x0 = 1.0 / (1.0 + h);
            y0 = h / (1.0 + h);
            lambda = (a + b) * y - b;
        }

        let mut e = -(lambda / a);

        let u = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (x / x0).ln()
        };

        e = lambda / b;

        let v = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (y / y0).ln()
        };

        let z = esum(mu, -(a * u + b * v));
        beta_rcomp1 = CONST * (b * x0).sqrt() * z * (-bcorr(a, b)).exp();
    } else {
        // Procedure for A < 8 or B < 8.
        let lnx;
        let lny;
        if x <= 0.375 {
            lnx = x.ln();
            lny = alnrel(-x);
        } else if y <= 0.375 {
            lnx = alnrel(-y);
            lny = y.ln();
        } else {
            lnx = x.ln();
            lny = y.ln();
        }

        let mut z = a * lnx + b * lny;

        if 1.0 <= a0 {
            z -= beta_log(a, b);
            beta_rcomp1 = esum(mu, z);
            return beta_rcomp1;
        }

        // Procedure for A < 1 or B < 1.
        let mut b0 = a.max(b);

        if 8.0 <= b0 {
            let u = gamma_ln1(a0) + algdiv(a0, b0);
            beta_rcomp1 = a0 * esum(mu, z - u);
            return beta_rcomp1;
        }

        if 1.0 < b0 {
            // Algorithm for 1 < B0 < 8.
            let mut u = gamma_ln1(a0);
            let n = (b0 - 1.0) as i32;

            let mut c = 1.0;
            for _ in 1..=n {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u = c.ln() + u;

            z -= u;
            b0 -= 1.0;
            let apb = a0 + b0;

            let t = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                u = a0 + b0 - 1.0;
                (1.0 + gam1(u)) / apb
            };

            beta_rcomp1 = a0 * esum(mu, z) * (1.0 + gam1(b0)) / t;
        } else {
            // Algorithm for B0 <= 1.
            beta_rcomp1 = esum(mu, z);
            if beta_rcomp1 == 0.0 {
                return beta_rcomp1;
            }

            let apb = a + b;

            if apb <= 1.0 {
                z = 1.0 + gam1(apb);
            } else {
                let u = a + b - 1.0;
                z = (1.0 + gam1(u)) / apb;
            }

            let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
            beta_rcomp1 = beta_rcomp1 * (a0 * c) / (1.0 + a0 / b0);
        }
    }

    beta_rcomp1
}

/// Evaluates *Iₓ*(*a*, *b*) − *Iₓ*(*a* + *n*, *b*) where *n* is a positive
/// integer (cdflib.f90:2196).
///
/// *a* and *b* should be nonnegative; *n* is the increment to the first
/// argument of *Iₓ*; *eps* is the tolerance.
#[inline]
#[allow(clippy::assign_op_pattern, clippy::collapsible_if)]
pub fn beta_up(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64) -> f64 {
    use super::exparg;

    // Obtain the scaling factor exp(-MU) and
    // exp(MU) * (X^A * Y^B / BETA(A,B)) / A.
    let apb = a + b;
    let ap1 = a + 1.0;
    let mut mu = 0;
    let mut d = 1.0;

    if n != 1 {
        if 1.0 <= a {
            if 1.1 * ap1 <= apb {
                mu = exparg(1).abs() as i32;
                let k = exparg(0) as i32;
                if k < mu {
                    mu = k;
                }
                let t = mu as f64;
                d = (-t).exp();
            }
        }
    }

    let beta_up = beta_rcomp1(mu, a, b, x, y) / a;

    if n == 1 || beta_up == 0.0 {
        return beta_up;
    }

    let mut w = d;

    // Let K be the index of the maximum term.
    let mut k = 0;

    if 1.0 < b {
        if y <= 0.0001 {
            k = n - 1;
        } else {
            let r = (b - 1.0) * x / y - a;

            if 1.0 <= r {
                k = n - 1;
                let t = (n - 1) as f64;
                if r < t {
                    k = r as i32;
                }
            }
        }

        // Add the increasing terms of the series.
        for i in 1..=k {
            let l = (i - 1) as f64;
            d = ((apb + l) / (ap1 + l)) * x * d;
            w += d;
        }
    }

    // Add the remaining terms of the series.
    for i in (k + 1)..=(n - 1) {
        let l = (i - 1) as f64;
        d = ((apb + l) / (ap1 + l)) * x * d;
        w += d;
        if d <= eps * w {
            return beta_up * w;
        }
    }

    beta_up * w
}

/// Evaluates the incomplete Γ ratio functions *P*(*a*, *x*) and
/// *Q*(*a*, *x*) (cdflib.f90:12208).
///
/// It is assumed that *a* ≤ 1. The argument *r* is the value
/// exp(−*x*) · *x*^*a* / Γ(*a*), and *eps* is the tolerance. Returns
/// (*p*, *q*), the values of *P*(*a*, *x*) and *Q*(*a*, *x*).
#[inline]
pub fn gamma_rat1(a: f64, x: f64, r: f64, eps: f64) -> (f64, f64) {
    use super::erf::{error_f, error_fc};

    let p;
    let mut q;

    if a * x == 0.0 {
        if x <= a {
            p = 0.0;
            q = 1.0;
        } else {
            p = 1.0;
            q = 0.0;
        }

        return (p, q);
    }

    if a == 0.5 {
        if x < 0.25 {
            p = error_f(x.sqrt());
            q = 0.5 + (0.5 - p);
        } else {
            q = error_fc(x.sqrt());
            p = 0.5 + (0.5 - q);
        }

        return (p, q);
    }

    // Rust only: with a or x NaN or infinite, the F90 Taylor series or
    // continued fraction below never exits, except for an infinite a with
    // x < 1.1, where the Taylor series exits at once and j = a * x * 0 is
    // NaN.
    if !a.is_finite() || !x.is_finite() {
        return (f64::NAN, f64::NAN);
    }

    // Taylor series for P(A,X)/X^A.
    if x < 1.1 {
        let mut an = 3.0;
        let mut c = x;
        let mut sum1 = x / (a + 3.0);
        let tol = 0.1 * eps / (a + 1.0);

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

        'l50: {
            'l40: {
                'l30: {
                    if x < 0.25 {
                        break 'l30;
                    }

                    if a < x / 2.59 {
                        break 'l50;
                    } else {
                        break 'l40;
                    }
                }

                // Label 30.
                if -0.13394 < z {
                    break 'l50;
                }
            }

            // Label 40.
            let w = z.exp();
            p = w * g * (0.5 + (0.5 - j));
            q = 0.5 + (0.5 - p);
            return (p, q);
        }

        // Label 50.
        let l = rexp(z);
        let w = 0.5 + (0.5 + l);
        q = (w * j - l) * g - h;

        if q < 0.0 {
            p = 1.0;
            q = 0.0;
        } else {
            p = 0.5 + (0.5 - q);
        }
    } else {
        // Continued fraction expansion.
        let mut a2nm1 = 1.0;
        let mut a2n = 1.0;
        let mut b2nm1 = x;
        let mut b2n = x + (1.0 - a);
        let mut c = 1.0;
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

            if (an0 - am0).abs() < eps * an0 {
                break;
            }
        }

        q = r * an0;
        p = 0.5 + (0.5 - q);
    }

    (p, q)
}

/// Failure modes of [`beta_grat`].
///
/// CDFLIB's `beta_grat` reports each of them as `ierr = 1` and returns
/// without changing *w*; callers recover with `result.unwrap_or(w)`, as
/// CDFLIB's `beta_inc` does by ignoring `ierr`.
///
/// [`beta_grat`]: crate::special::internal::beta_grat
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum BetaGratError {
    /// *b* · *z* evaluated to zero (CDFLIB `ierr = 1`).
    #[error("b·z evaluated to zero")]
    BzZero,
    /// The scale *r* · exp(−(algdiv(*b*, *a*) + *b* · ln *nu*))
    /// underflowed to zero (CDFLIB `ierr = 1`).
    #[error("u underflowed to zero")]
    UnderflowedScale,
    /// The partial sum of the expansion became nonpositive (CDFLIB
    /// `ierr = 1`).
    #[error("partial sum went non-positive")]
    NonPositiveSum,
}

/// Evaluates an asymptotic expansion for *Iₓ*(*a*, *b*) (cdflib.f90:763).
///
/// *a* and *b* are the parameters of the function and should be
/// nonnegative; it is assumed that 15 ≤ *a* and *b* ≤ 1, and that *b* is
/// less than *a*. *x* is the argument of the function and should satisfy
/// 0 ≤ *x* ≤ 1; *y* should equal 1 − *x*. *w* is a quantity to which the
/// result of the computation is added, and *eps* is a tolerance.
///
/// Returns *w* plus the expansion. CDFLIB's `ierr = 1` is returned as a
/// [`BetaGratError`]; in that case CDFLIB leaves *w* unchanged, so callers
/// recover with `result.unwrap_or(w)`.
///
/// [`BetaGratError`]: crate::special::internal::BetaGratError
#[inline]
pub fn beta_grat(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    mut w: f64,
    eps: f64,
) -> Result<f64, BetaGratError> {
    let mut c = [0.0_f64; 30];
    let mut d = [0.0_f64; 30];

    let bm1 = (b - 0.5) - 0.5;
    let nu = a + 0.5 * bm1;

    let lnx = if y <= 0.375 { alnrel(-y) } else { x.ln() };

    let z = -(nu * lnx);

    if b * z == 0.0 {
        // ierr = 1, mapped onto BetaGratError.
        return Err(BetaGratError::BzZero);
    }

    // Computation of the expansion. Set R = EXP(-Z)*Z^B/GAMMA(B).
    let mut r = b * (1.0 + gam1(b)) * (b * z.ln()).exp();
    r = r * (a * lnx).exp() * (0.5 * bm1 * lnx).exp();
    let mut u = algdiv(b, a) + b * nu.ln();
    u = r * (-u).exp();

    if u == 0.0 {
        // ierr = 1, mapped onto BetaGratError.
        return Err(BetaGratError::UnderflowedScale);
    }

    let (_p, q) = gamma_rat1(b, z, r, eps);

    let v = 0.25 * pow2(1.0 / nu);
    let t2 = 0.25 * lnx * lnx;
    let l = w / u;
    let mut j = q / r;
    let mut sum1 = j;
    let mut t = 1.0;
    let mut cn = 1.0;
    let mut n2 = 0.0;

    for n in 1..=30 {
        let bp2n = b + n2;
        j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
        n2 += 2.0;
        t *= t2;
        cn /= n2 * (n2 + 1.0);
        c[n - 1] = cn;
        let mut s = 0.0;

        let mut coef = b - n as f64;
        for i in 1..n {
            s += coef * c[i - 1] * d[n - i - 1];
            coef += b;
        }

        d[n - 1] = bm1 * cn + s / n as f64;
        let dj = d[n - 1] * j;
        sum1 += dj;

        if sum1 <= 0.0 {
            // ierr = 1, mapped onto BetaGratError.
            return Err(BetaGratError::NonPositiveSum);
        }

        if dj.abs() <= eps * (sum1 + l) {
            w += u * sum1;
            return Ok(w);
        }
    }

    w += u * sum1;
    Ok(w)
}

/// Computes an asymptotic expansion for *Iₓ*(*a*, *b*), for large *a* and
/// *b* (cdflib.f90:439).
///
/// *a* and *b* are the parameters of the function and should be
/// nonnegative; it is assumed that both *a* and *b* are greater than or
/// equal to 15. *lambda* is the value of (*a* + *b*) · *y* − *b*, assumed
/// nonnegative, and *eps* is the tolerance.
#[inline]
#[allow(clippy::assign_op_pattern, clippy::needless_late_init)]
pub fn beta_asym(a: f64, b: f64, lambda: f64, eps: f64) -> f64 {
    const NUM: usize = 20;
    const E0: f64 = 1.12837916709551;
    const E1: f64 = 0.353553390593274;

    let mut a0 = [0.0_f64; NUM + 1];
    let mut b0 = [0.0_f64; NUM + 1];
    let mut c = [0.0_f64; NUM + 1];
    let mut d = [0.0_f64; NUM + 1];

    let mut beta_asym = 0.0;

    let h;
    let r0;
    let r1;
    let w0;

    if a < b {
        h = a / b;
        r0 = 1.0 / (1.0 + h);
        r1 = (b - a) / b;
        w0 = 1.0 / (a * (1.0 + h)).sqrt();
    } else {
        h = b / a;
        r0 = 1.0 / (1.0 + h);
        r1 = (b - a) / a;
        w0 = 1.0 / (b * (1.0 + h)).sqrt();
    }

    let f = a * rlog1(-(lambda / a)) + b * rlog1(lambda / b);
    let t = (-f).exp();
    if t == 0.0 {
        return beta_asym;
    }

    let z0 = f.sqrt();
    let z = 0.5 * (z0 / E1);
    let z2 = f + f;

    a0[0] = (2.0 / 3.0) * r1;
    c[0] = -(0.5 * a0[0]);
    d[0] = -c[0];
    let mut j0 = (0.5 / E0) * error_fc_scaled(z0);
    let mut j1 = E1;
    let mut sum1 = j0 + d[0] * w0 * j1;

    let mut s = 1.0;
    let h2 = h * h;
    let mut hn = 1.0;
    let mut w = w0;
    let mut znm1 = z;
    let mut zn = z2;

    for n in (2..=NUM).step_by(2) {
        hn = h2 * hn;
        a0[n - 1] = 2.0 * r0 * (1.0 + h * hn) / (n as f64 + 2.0);
        let np1 = n + 1;
        s += hn;
        a0[np1 - 1] = 2.0 * r1 * s / (n as f64 + 3.0);

        for i in n..=np1 {
            let r = -(0.5 * (i as f64 + 1.0));
            b0[0] = r * a0[0];
            for m in 2..=i {
                let mut bsum = 0.0;
                let mm1 = m - 1;
                for j in 1..=mm1 {
                    let mmj = m - j;
                    bsum += (j as f64 * r - mmj as f64) * a0[j - 1] * b0[mmj - 1];
                }
                b0[m - 1] = r * a0[m - 1] + bsum / m as f64;
            }

            c[i - 1] = b0[i - 1] / (i as f64 + 1.0);

            let mut dsum = 0.0;
            for j in 1..i {
                dsum += d[i - j - 1] * c[j - 1];
            }
            d[i - 1] = -(dsum + c[i - 1]);
        }

        j0 = E1 * znm1 + (n as f64 - 1.0) * j0;
        j1 = E1 * zn + n as f64 * j1;
        znm1 = z2 * znm1;
        zn = z2 * zn;
        w = w0 * w;
        let t0 = d[n - 1] * w * j0;
        w = w0 * w;
        let t1 = d[np1 - 1] * w * j1;
        sum1 += t0 + t1;

        if (t0.abs() + t1.abs()) <= eps * sum1 {
            let u = (-bcorr(a, b)).exp();
            beta_asym = E0 * t * u * sum1;
            return beta_asym;
        }
    }

    let u = (-bcorr(a, b)).exp();
    beta_asym = E0 * t * u * sum1;
    beta_asym
}

/// Evaluates a continued fraction expansion for *Iₓ*(*a*, *b*)
/// (cdflib.f90:626).
///
/// *a* and *b* are the parameters of the function and should be
/// nonnegative; it is assumed that both *a* and *b* are greater than 1.
/// *x* is the argument of the function and should satisfy 0 ≤ *x* ≤ 1;
/// *y* should equal 1 − *x*. *lambda* is the value of
/// (*a* + *b*) · *y* − *b*, and *eps* is a tolerance.
#[inline]
pub fn beta_frac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64) -> f64 {
    let mut beta_frac = beta_rcomp(a, b, x, y);

    if beta_frac == 0.0 {
        return beta_frac;
    }

    let c = 1.0 + lambda;
    let c0 = b / a;
    let c1 = 1.0 + 1.0 / a;
    let yp1 = y + 1.0;

    let mut n = 0.0;
    let mut p = 1.0;
    let mut s = a + 1.0;
    let mut an = 0.0;
    let mut bn = 1.0;
    let mut anp1 = 1.0;
    let mut bnp1 = c / c1;
    let mut r = c1 / c;

    // Continued fraction calculation.
    loop {
        n += 1.0;
        let mut t = n / a;
        let w = n * (b - n) * x;
        let mut e = a / s;
        let alpha = (p * (p + c0) * e * e) * (w * x);
        e = (1.0 + t) / (c1 + t + t);
        let beta = n + w / s + e * (c + n * yp1);
        p = 1.0 + t;
        s += 2.0;

        // Update AN, BN, ANP1, and BNP1.
        t = alpha * an + beta * anp1;
        an = anp1;
        anp1 = t;
        t = alpha * bn + beta * bnp1;
        bn = bnp1;
        bnp1 = t;

        let r0 = r;
        r = anp1 / bnp1;

        if (r - r0).abs() <= eps * r {
            beta_frac *= r;
            break;
        }

        // Rust only: once r is NaN the F90 loop never exits.
        if r.is_nan() {
            return f64::NAN;
        }

        // Rescale AN, BN, ANP1, and BNP1.
        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0;
    }

    beta_frac
}

/// Errors of [`beta_inc`], one variant per nonzero `ierr` of CDFLIB's
/// `beta_inc`, in CDFLIB order.
///
/// On error CDFLIB returns *w* = *w*₁ = 0, except for `ierr = 6`, where
/// it returns *w* = 0 and *w*₁ = 1.
///
/// [`beta_inc`]: crate::special::beta_inc
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum BetaIncError {
    /// *a* or *b* is negative (CDFLIB `ierr = 1`).
    #[error("a or b is negative: a = {a}, b = {b}")]
    NegativeParameter { a: f64, b: f64 },
    /// *a* = *b* = 0 (CDFLIB `ierr = 2`).
    #[error("both a and b are zero")]
    BothZero,
    /// *x* < 0 or 1 < *x* (CDFLIB `ierr = 3`).
    #[error("x must be in [0..1], got {0}")]
    XOutOfRange(f64),
    /// *y* < 0 or 1 < *y* (CDFLIB `ierr = 4`).
    #[error("y must be in [0..1], got {0}")]
    YOutOfRange(f64),
    /// *x* + *y* ≠ 1, that is, 3*ε* < |*x* + *y* − 1| (CDFLIB `ierr = 5`).
    #[error("x + y must equal 1, got x = {x}, y = {y}")]
    InconsistentSum { x: f64, y: f64 },
    /// *x* = *a* = 0 (CDFLIB `ierr = 6`).
    #[error("degenerate: x = 0 and a = 0")]
    XZeroAndAZero,
    /// *y* = *b* = 0 (CDFLIB `ierr = 7`).
    #[error("degenerate: y = 0 and b = 0")]
    YZeroAndBZero,
}

/// Evaluates the incomplete Β function *Iₓ*(*a*, *b*) (cdflib.f90:928).
///
/// *a* and *b* are the parameters of the function and should be
/// nonnegative. *x* is the argument of the function and should satisfy
/// *x* ∈ [0 . . 1]; *y* should equal 1 − *x*. Returns (*w*, *w*₁), the values of
/// *Iₓ*(*a*, *b*) and 1 − *Iₓ*(*a*, *b*).
///
/// The caller supplies *y* rather than letting the routine compute 1 − *x*,
/// because in the deep tail that subtraction would lose digits. Each series
/// or expansion computes one of *w* and *w*₁ and derives the other as
/// 0.5 + (0.5 − ·), as with the (*p*, *q*) pair returned by [`gamma_inc`].
///
/// # Panics
///
/// Panics on a [`BetaIncError`] (a nonzero CDFLIB `ierr`). Use
/// [`try_beta_inc`] for the fallible form.
///
/// # Example
///
/// ```
/// use cdflib::special::beta_inc;
///
/// let (w, _w1) = beta_inc(2.0, 5.0, 0.3, 0.7);
/// assert!((w - 0.579825).abs() < 1e-6);
/// ```
///
/// [`gamma_inc`]: crate::special::gamma_inc
/// [`BetaIncError`]: crate::special::BetaIncError
/// [`try_beta_inc`]: crate::special::try_beta_inc
#[inline]
pub fn beta_inc(a: f64, b: f64, x: f64, y: f64) -> (f64, f64) {
    try_beta_inc(a, b, x, y).unwrap_or_else(|e| panic!("beta_inc({a}, {b}, {x}, {y}): {e}"))
}

/// Fallible form of [`beta_inc`]: returns a [`BetaIncError`] where CDFLIB's
/// `beta_inc` sets a nonzero `ierr`.
///
/// # Example
///
/// ```
/// use cdflib::special::{try_beta_inc, BetaIncError};
///
/// let (w, _w1) = try_beta_inc(2.0, 5.0, 0.3, 0.7).unwrap();
/// assert!((w - 0.579825).abs() < 1e-6);
/// assert!(matches!(
///     try_beta_inc(-1.0, 1.0, 0.5, 0.5),
///     Err(BetaIncError::NegativeParameter { .. }),
/// ));
/// ```
///
/// [`BetaIncError`]: crate::special::BetaIncError
#[inline]
// The F90 initialisation w = 0 is never read, since every error return is
// an Err; w1 = 0 is read by the direct go to 150 (cdflib.f90:1164).
#[allow(unused_assignments, clippy::manual_swap)]
pub fn try_beta_inc(a: f64, b: f64, x: f64, y: f64) -> Result<(f64, f64), BetaIncError> {
    let mut eps = f64::EPSILON;
    let mut w = 0.0;
    let mut w1 = 0.0;

    if a < 0.0 || b < 0.0 {
        // ierr = 1, mapped onto BetaIncError.
        return Err(BetaIncError::NegativeParameter { a, b });
    }

    if a == 0.0 && b == 0.0 {
        // ierr = 2, mapped onto BetaIncError.
        return Err(BetaIncError::BothZero);
    }

    // A NaN x or y passes these two tests, as in the F90.
    if x < 0.0 || 1.0 < x {
        // ierr = 3, mapped onto BetaIncError.
        return Err(BetaIncError::XOutOfRange(x));
    }

    if y < 0.0 || 1.0 < y {
        // ierr = 4, mapped onto BetaIncError.
        return Err(BetaIncError::YOutOfRange(y));
    }

    let z = ((x + y) - 0.5) - 0.5;

    if 3.0 * eps < z.abs() {
        // ierr = 5, mapped onto BetaIncError.
        return Err(BetaIncError::InconsistentSum { x, y });
    }

    if x == 0.0 {
        w = 0.0;
        w1 = 1.0;
        if a == 0.0 {
            // ierr = 6, mapped onto BetaIncError.
            return Err(BetaIncError::XZeroAndAZero);
        }
        return Ok((w, w1));
    }

    if y == 0.0 {
        if b == 0.0 {
            // ierr = 7, mapped onto BetaIncError.
            return Err(BetaIncError::YZeroAndBZero);
        }
        w = 1.0;
        w1 = 0.0;
        return Ok((w, w1));
    }

    if a == 0.0 {
        w = 1.0;
        w1 = 0.0;
        return Ok((w, w1));
    }

    if b == 0.0 {
        w = 0.0;
        w1 = 1.0;
        return Ok((w, w1));
    }

    eps = eps.max(1.0e-15);

    'l260: {
        if a.max(b) < 0.001 * eps {
            break 'l260;
        }

        // Rust only: past this point, given a NaN argument, the F90 either
        // never exits a loop, returns NaN, or returns a value computed from
        // the arguments that are not NaN (for example, from y alone at label
        // 120); return NaN.
        if a.is_nan() || b.is_nan() || x.is_nan() || y.is_nan() {
            return Ok((f64::NAN, f64::NAN));
        }

        let mut ind = 0;
        let mut a0 = a;
        let mut b0 = b;
        let mut x0 = x;
        let mut y0 = y;
        let n;

        // Each go to 90 ... 200 below calls the fn for that label, which
        // returns (w, w1); the go to 250 that ends every label is the break
        // out of 'l250 carrying that pair.
        (w, w1) = 'l250: {
            'l40: {
                if 1.0 < a0.min(b0) {
                    break 'l40;
                }

                // Procedure for A0 <= 1 or B0 <= 1.
                if 0.5 < x {
                    ind = 1;
                    a0 = b;
                    b0 = a;
                    x0 = y;
                    y0 = x;
                }

                if b0 < eps.min(eps * a0) {
                    break 'l250 label_90(a0, b0, x0, eps);
                }

                if a0 < eps.min(eps * b0) && b0 * x0 <= 1.0 {
                    break 'l250 label_100(a0, b0, x0, eps);
                }

                'l20: {
                    if 1.0 < a0.max(b0) {
                        break 'l20;
                    }

                    if 0.2_f64.min(b0) <= a0 {
                        break 'l250 label_110(a0, b0, x0, eps);
                    }

                    if x0.powf(a0) <= 0.9 {
                        break 'l250 label_110(a0, b0, x0, eps);
                    }

                    if 0.3 <= x0 {
                        break 'l250 label_120(a0, b0, y0, eps);
                    }

                    n = 20;
                    break 'l250 label_140(a0, b0, x0, y0, n, eps);
                }

                // Label 20.
                if b0 <= 1.0 {
                    break 'l250 label_110(a0, b0, x0, eps);
                }

                if 0.3 <= x0 {
                    break 'l250 label_120(a0, b0, y0, eps);
                }

                'l30: {
                    if 0.1 <= x0 {
                        break 'l30;
                    }

                    if (x0 * b0).powf(a0) <= 0.7 {
                        break 'l250 label_110(a0, b0, x0, eps);
                    }
                }

                // Label 30.
                if 15.0 < b0 {
                    // w1 is still the 0 of cdflib.f90:1007.
                    break 'l250 label_150(a0, b0, x0, y0, w1, eps);
                }

                n = 20;
                break 'l250 label_140(a0, b0, x0, y0, n, eps);
            }

            // Label 40: procedure for 1 < A0 and 1 < B0.
            let mut lambda = if a <= b {
                a - (a + b) * x
            } else {
                (a + b) * y - b
            };

            if lambda < 0.0 {
                ind = 1;
                a0 = b;
                b0 = a;
                x0 = y;
                y0 = x;
                lambda = lambda.abs();
            }

            if b0 < 40.0 && b0 * x0 <= 0.7 {
                break 'l250 label_110(a0, b0, x0, eps);
            }

            if b0 < 40.0 {
                break 'l250 label_160(a0, b0, x0, y0, eps);
            }

            'l80: {
                if b0 < a0 {
                    break 'l80;
                }

                if a0 <= 100.0 {
                    break 'l250 label_130(a0, b0, x0, y0, lambda, eps);
                }

                if 0.03 * a0 < lambda {
                    break 'l250 label_130(a0, b0, x0, y0, lambda, eps);
                }

                break 'l250 label_200(a0, b0, lambda, eps);
            }

            // Label 80.
            if b0 <= 100.0 {
                break 'l250 label_130(a0, b0, x0, y0, lambda, eps);
            }

            if 0.03 * b0 < lambda {
                break 'l250 label_130(a0, b0, x0, y0, lambda, eps);
            }

            break 'l250 label_200(a0, b0, lambda, eps);
        };

        // Evaluation of the appropriate algorithm.

        // Label 90 (cdflib.f90:1227-1231).
        fn label_90(a0: f64, b0: f64, x0: f64, eps: f64) -> (f64, f64) {
            let w = fpser(a0, b0, x0, eps);
            let w1 = 0.5 + (0.5 - w);
            (w, w1)
        }

        // Label 100 (cdflib.f90:1233-1237).
        fn label_100(a0: f64, b0: f64, x0: f64, eps: f64) -> (f64, f64) {
            let w1 = apser(a0, b0, x0, eps);
            let w = 0.5 + (0.5 - w1);
            (w, w1)
        }

        // Label 110 (cdflib.f90:1239-1243).
        fn label_110(a0: f64, b0: f64, x0: f64, eps: f64) -> (f64, f64) {
            let w = beta_pser(a0, b0, x0, eps);
            let w1 = 0.5 + (0.5 - w);
            (w, w1)
        }

        // Label 120 (cdflib.f90:1245-1249).
        fn label_120(a0: f64, b0: f64, y0: f64, eps: f64) -> (f64, f64) {
            let w1 = beta_pser(b0, a0, y0, eps);
            let w = 0.5 + (0.5 - w1);
            (w, w1)
        }

        // Label 130 (cdflib.f90:1251-1255).
        fn label_130(a0: f64, b0: f64, x0: f64, y0: f64, lambda: f64, eps: f64) -> (f64, f64) {
            let w = beta_frac(a0, b0, x0, y0, lambda, 15.0 * eps);
            let w1 = 0.5 + (0.5 - w);
            (w, w1)
        }

        // Label 140 (cdflib.f90:1257-1260), which falls through to label 150.
        fn label_140(a0: f64, mut b0: f64, x0: f64, y0: f64, n: i32, eps: f64) -> (f64, f64) {
            let w1 = beta_up(b0, a0, y0, x0, n, eps);
            b0 += n as f64;
            label_150(a0, b0, x0, y0, w1, eps)
        }

        // Label 150 (cdflib.f90:1262-1266). The F90 does not test ierr1: on
        // failure beta_grat leaves w1 unchanged.
        fn label_150(a0: f64, b0: f64, x0: f64, y0: f64, mut w1: f64, eps: f64) -> (f64, f64) {
            w1 = beta_grat(b0, a0, y0, x0, w1, 15.0 * eps).unwrap_or(w1);
            let w = 0.5 + (0.5 - w1);
            (w, w1)
        }

        // Label 160 (cdflib.f90:1268-1298). The F90 does not test ierr1: on
        // failure beta_grat leaves w unchanged.
        fn label_160(mut a0: f64, mut b0: f64, x0: f64, y0: f64, eps: f64) -> (f64, f64) {
            let mut n = b0 as i32;
            b0 -= n as f64;

            if b0 == 0.0 {
                n -= 1;
                b0 = 1.0;
            }

            let mut w = beta_up(b0, a0, y0, x0, n, eps);

            if x0 <= 0.7 {
                w += beta_pser(a0, b0, x0, eps);
                let w1 = 0.5 + (0.5 - w);
                return (w, w1);
            }

            if a0 <= 15.0 {
                n = 20;
                w += beta_up(a0, b0, x0, y0, n, eps);
                a0 += n as f64;
            }

            w = beta_grat(a0, b0, x0, y0, w, 15.0 * eps).unwrap_or(w);
            let w1 = 0.5 + (0.5 - w);
            (w, w1)
        }

        // Label 200 (cdflib.f90:1300-1304).
        fn label_200(a0: f64, b0: f64, lambda: f64, eps: f64) -> (f64, f64) {
            let w = beta_asym(a0, b0, lambda, 100.0 * eps);
            let w1 = 0.5 + (0.5 - w);
            (w, w1)
        }

        // Label 250: termination of the procedure.
        if ind != 0 {
            let t = w;
            w = w1;
            w1 = t;
        }

        return Ok((w, w1));
    }

    // Label 260: procedure for A and B < 0.001 * EPS.
    w = b / (a + b);
    w1 = a / (a + b);

    Ok((w, w1))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================================ dbetrm

    // The residual at (100, 100) sits right at the native-FPU 1e-12 edge;
    // miri's soft-float ln pushes it over. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn dbetrm_matches_beta_log_minus_sterling() {
        // For each (a, b), dbetrm should equal ln Β(a, b) − Sterling decomposition.
        const HLN2PI: f64 = 0.91893853320467274178;
        fn sterling(z: f64) -> f64 {
            HLN2PI + (z - 0.5) * z.ln() - z
        }
        for &(a, b) in &[(2.5_f64, 3.5), (10.0, 20.0), (50.0, 60.0), (100.0, 100.0)] {
            let r = dbetrm(a, b);
            let lnb = beta_log(a, b);
            let sterling_sum = sterling(a) + sterling(b) - sterling(a + b);
            let expected = lnb - sterling_sum;
            assert!(
                (r - expected).abs() < 1e-12,
                "a={a}, b={b}: dbetrm={r}, expected={expected}"
            );
        }
    }

    #[test]
    fn dbetrm_decreases_for_large_args() {
        // The Sterling remainder shrinks as a, b grow.
        let r10 = dbetrm(10.0, 10.0);
        let r100 = dbetrm(100.0, 100.0);
        assert!(r100.abs() < r10.abs());
        assert!(r100.abs() < 0.01);
    }

    // ============================================================ beta_log

    #[test]
    fn beta_log_at_integer_arguments() {
        // ln Β(1, 1) = 0
        assert!(beta_log(1.0, 1.0).abs() < 1e-14);
        // ln Β(2, 2) = ln(1/6) = -ln 6
        assert!((beta_log(2.0, 2.0) - (-6.0_f64.ln())).abs() < 1e-13);
        // ln Β(3, 4) = ln(Γ(3)Γ(4)/Γ(7)) = ln(2·6/720) = ln(1/60)
        assert!((beta_log(3.0, 4.0) - (-60.0_f64.ln())).abs() < 1e-13);
    }

    #[test]
    fn beta_inc_at_x_half_with_a_b_equal() {
        // The incomplete Β function at x = 0.5 with a = b is 0.5 by
        // symmetry.
        for &a in &[0.5, 1.0, 2.0, 5.0, 30.0] {
            let (w, w1) = beta_inc(a, a, 0.5, 0.5);
            assert!((w - 0.5).abs() < 1e-10, "a={a}: w={w}");
            assert!((w1 - 0.5).abs() < 1e-10);
        }
    }

    #[test]
    fn beta_inc_at_boundaries() {
        assert_eq!(try_beta_inc(2.0, 3.0, 0.0, 1.0), Ok((0.0, 1.0)));
        assert_eq!(try_beta_inc(2.0, 3.0, 1.0, 0.0), Ok((1.0, 0.0)));
    }

    #[test]
    fn beta_inc_p_plus_q_equals_one() {
        for &(a, b) in &[(1.0, 1.0), (2.0, 5.0), (10.0, 20.0), (0.5, 3.0)] {
            for x in [0.1, 0.3, 0.5, 0.7, 0.9] {
                let (w, w1) = beta_inc(a, b, x, 1.0 - x);
                assert!((w + w1 - 1.0).abs() < 1e-12, "a={a}, b={b}, x={x}");
            }
        }
    }

    // ===== Validation-error paths (each error variant) =====

    #[test]
    fn beta_inc_negative_parameter() {
        assert!(matches!(
            try_beta_inc(-1.0, 2.0, 0.5, 0.5),
            Err(BetaIncError::NegativeParameter { .. })
        ));
        assert!(matches!(
            try_beta_inc(2.0, -1.0, 0.5, 0.5),
            Err(BetaIncError::NegativeParameter { .. })
        ));
    }

    #[test]
    fn beta_inc_both_zero() {
        assert_eq!(
            try_beta_inc(0.0, 0.0, 0.5, 0.5),
            Err(BetaIncError::BothZero)
        );
    }

    #[test]
    fn beta_inc_x_out_of_range() {
        assert!(matches!(
            try_beta_inc(2.0, 3.0, -0.1, 1.1),
            Err(BetaIncError::XOutOfRange(-0.1))
        ));
        assert!(matches!(
            try_beta_inc(2.0, 3.0, 1.1, -0.1),
            Err(BetaIncError::XOutOfRange(_))
        ));
    }

    #[test]
    fn beta_inc_y_out_of_range() {
        // x in [0..1] but y not.
        assert!(matches!(
            try_beta_inc(2.0, 3.0, 0.5, -0.1),
            Err(BetaIncError::YOutOfRange(_))
        ));
        assert!(matches!(
            try_beta_inc(2.0, 3.0, 0.5, 1.1),
            Err(BetaIncError::YOutOfRange(_))
        ));
    }

    #[test]
    fn beta_inc_x_plus_y_not_one() {
        assert!(matches!(
            try_beta_inc(2.0, 3.0, 0.3, 0.5),
            Err(BetaIncError::InconsistentSum { .. })
        ));
    }

    #[test]
    fn beta_inc_x_zero_and_a_zero() {
        assert_eq!(
            try_beta_inc(0.0, 3.0, 0.0, 1.0),
            Err(BetaIncError::XZeroAndAZero)
        );
    }

    #[test]
    fn beta_inc_y_zero_and_b_zero() {
        assert_eq!(
            try_beta_inc(3.0, 0.0, 1.0, 0.0),
            Err(BetaIncError::YZeroAndBZero)
        );
    }

    #[test]
    fn beta_inc_a_zero_with_b_positive() {
        assert_eq!(try_beta_inc(0.0, 3.0, 0.5, 0.5), Ok((1.0, 0.0)));
    }

    #[test]
    fn beta_inc_b_zero_with_a_positive() {
        assert_eq!(try_beta_inc(3.0, 0.0, 0.5, 0.5), Ok((0.0, 1.0)));
    }

    #[test]
    fn beta_inc_both_tiny_a_b() {
        // max(a, b) < 0.001 * eps: label 260 returns (b/(a+b), a/(a+b)).
        let tiny = 1e-20;
        let (w, w1) = beta_inc(tiny, tiny, 0.5, 0.5);
        // Both ratios equal 0.5 by symmetry.
        assert!((w - 0.5).abs() < 1e-10);
        assert!((w1 - 0.5).abs() < 1e-10);
    }

    // ===== Regime switching points =====

    #[test]
    fn beta_inc_small_a_large_b_uses_grat_path() {
        // a0 ≤ 1 and 15 < b0: label 30 goes to label 150 (beta_grat).
        let (w, w1) = beta_inc(0.5, 30.0, 0.05, 0.95);
        // Sanity: numerically plausible.
        assert!(w > 0.0 && w < 1.0 && (w + w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn beta_inc_both_moderate_uses_frac_path() {
        // lambda > 0 at label 40, 40 ≤ b0 and a0 ≤ 100: label 130 (beta_frac,
        // through beta_rcomp's 8 ≤ a0 path).
        let (w, w1) = beta_inc(10.0, 60.0, 0.1, 0.9);
        assert!((w + w1 - 1.0).abs() < 1e-10);
    }

    // 1e-22 absolute tolerance on a tail probability is well below
    // miri's soft-float libm precision. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn beta_inc_extreme_skew_matches_high_precision_reference() {
        let (w, w1) = beta_inc(0.5, 100.0, 0.15, 0.85);
        assert!((w - 0.999_999_987_603_646_8).abs() < 1e-15);
        assert!((w1 - 1.239_635_319_310_601_4e-8).abs() < 1e-22);
    }

    #[test]
    fn beta_inc_a_large_b_moderate() {
        // b < a with 0 ≤ lambda: no swap at label 40, and b0 < 40 goes to
        // label 160 (beta_grat with 15 < a0).
        let (w, w1) = beta_inc(60.0, 10.0, 0.85, 0.15);
        assert!((w + w1 - 1.0).abs() < 1e-10);
        // Symmetric to the test above, since Ix(a, b) = 1 - Iy(b, a) with
        // y = 1 - x. Here lambda < 0, so label 40 swaps into the same
        // label 160 call.
        let (w2, _) = beta_inc(10.0, 60.0, 0.15, 0.85);
        assert!((w - (1.0 - w2)).abs() < 1e-10);
    }

    #[test]
    fn beta_inc_extreme_lambda_asym_path() {
        // 100 < a0, b0 with lambda ≤ 0.03 * min(a0, b0): label 200 (beta_asym).
        // Two parameter orderings to cover both branches of beta_asym.
        // lambda < 0 at label 40 swaps to b0 < a0, so label 80 goes to
        // label 200 with a0 > b0 (beta_asym else branch).
        let (w, w1) = beta_inc(150.0, 200.0, 150.0 / 350.0 + 0.001, 200.0 / 350.0 - 0.001);
        assert!((w + w1 - 1.0).abs() < 1e-8);
        // No swap, a0 ≤ b0, 100 < a0 and lambda ≤ 0.03 * a0.
        // mean = 150/550 ≈ 0.2727; x just below mean → lambda small positive.
        let (w, w1) = beta_inc(150.0, 400.0, 0.272, 0.728);
        assert!((w + w1 - 1.0).abs() < 1e-8);
    }

    // The bit-for-bit comparison needs two evaluations of the same libm
    // calls to agree, which miri's soft-float libm does not guarantee.
    // Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn beta_inc_a_below_eps_uses_apser() {
        // With eps = 1e-15, a0 < min(eps, eps * b0) and b0 * x0 <= 1
        // (cdflib.f90:1120): label 100 computes w1 with apser.
        let eps = f64::EPSILON.max(1e-15);
        let (a, b, x, y) = (1e-16, 5.0, 0.05, 0.95);
        let w1 = apser(a, b, x, eps);
        assert_eq!(beta_inc(a, b, x, y), (0.5 + (0.5 - w1), w1));
    }

    // The bit-for-bit comparison needs two evaluations of the same libm
    // calls to agree, which miri's soft-float libm does not guarantee.
    // Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn beta_inc_b_below_eps_uses_fpser() {
        // With eps = 1e-15, b0 < min(eps, eps * a0) (cdflib.f90:1116):
        // label 90 computes w with fpser.
        let eps = f64::EPSILON.max(1e-15);
        let (a, b, x, y) = (2.0, 3e-16, 0.5, 0.5);
        let w = fpser(a, b, x, eps);
        assert_eq!(beta_inc(a, b, x, y), (w, 0.5 + (0.5 - w)));
    }

    // ===== Direct helper-function tests =====

    #[test]
    fn esum_at_all_branches() {
        // x > 0, mu > 0: fallthrough to exp(mu)*exp(x).
        let r1 = esum(1, 2.0);
        assert!((r1 - (3.0_f64).exp()).abs() < 1e-12);
        // x > 0, mu < 0, mu+x < 0: fallthrough.
        let r2 = esum(-5, 1.0);
        assert!((r2 - (-4.0_f64).exp()).abs() < 1e-14);
        // x > 0, mu <= 0, mu+x >= 0: takes the early return.
        let r3 = esum(-1, 2.0);
        assert!((r3 - (1.0_f64).exp()).abs() < 1e-14);
        // x <= 0, mu >= 0, mu+x > 0: fallthrough.
        let r4 = esum(2, -1.0);
        assert!((r4 - (1.0_f64).exp()).abs() < 1e-14);
        // x <= 0, mu < 0: fallthrough.
        let r5 = esum(-1, -1.0);
        assert!((r5 - (-2.0_f64).exp()).abs() < 1e-14);
    }

    #[test]
    fn algdiv_in_b_le_a_branch() {
        // b < a: the first branch of cdflib.f90:71. The formula computes
        // ln(Γ(b)/Γ(a+b)) for b ≥ 8 (the precondition).
        // Compare against beta_log identity: ln Γ(b) - ln Γ(a+b).
        let a = 10.0;
        let b = 8.0;
        let direct = gamma_log(b) - gamma_log(a + b);
        let via_algdiv = algdiv(a, b);
        assert!((via_algdiv - direct).abs() < 1e-12);
    }

    #[test]
    fn beta_wrapper() {
        // Β(a,b) = exp(beta_log(a, b)). Verify at simple integer points.
        assert!((beta(1.0, 1.0) - 1.0).abs() < 1e-14);
        assert!((beta(2.0, 2.0) - 1.0 / 6.0).abs() < 1e-14);
        assert!((beta(3.0, 4.0) - 1.0 / 60.0).abs() < 1e-14);
    }

    #[test]
    fn beta_rcomp_at_extreme_b() {
        // a = 41 takes the 8 <= a0 path (cdflib.f90:1957-1986) with
        // b = 1e300, where x0 = h / (1 + h) is about 4e-299 and z
        // underflows to 0; the result must stay finite.
        let r = beta_rcomp(41.0, 1e300, 0.8, 0.2);
        assert!(r.is_finite(), "beta_rcomp returned non-finite: {r}");
    }

    #[test]
    fn beta_inc_a0_b0_at_most_one_corners() {
        // With max(a0, b0) ≤ 1 and a0 < min(0.2, b0), cdflib.f90:1132-1141
        // choose among three evaluations by x0 and x0^a0:
        //   x0^a0 ≤ 0.9: label 110 (beta_pser)
        //   0.9 < x0^a0 and 0.3 ≤ x0: label 120 (beta_pser(b0, a0, y0))
        //   0.9 < x0^a0 and x0 < 0.3: label 140 (beta_up + beta_grat)
        for &(a, b, x) in &[
            (0.1, 0.5, 0.3),  // x^a ≈ 0.887 ≤ 0.9: label 110
            (0.1, 0.5, 0.5),  // x^a ≈ 0.933 > 0.9, x ≥ 0.3: label 120
            (0.01, 0.5, 0.1), // x^a ≈ 0.977 > 0.9, x < 0.3: label 140
        ] {
            let (w, w1) = beta_inc(a, b, x, 1.0 - x);
            assert!((w + w1 - 1.0).abs() < 1e-10, "a={a}, b={b}, x={x}");
        }
    }

    #[test]
    fn beta_inc_label_160_x0_above_0_7() {
        // Label 160 with 0.7 < x0 and a0 ≤ 15: beta_up, then beta_up and
        // beta_grat with a0 + 20 (cdflib.f90:1288-1296).
        let (w, w1) = beta_inc(10.0, 3.0, 0.75, 0.25);
        assert!((w + w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fpser_apser_beta_pser_edge_cases() {
        let eps = f64::EPSILON.max(1e-15);

        // fpser: a*ln(x) < exparg(1) → return 0.
        // For a=1, x = 1e-308 (denormal) makes t ≈ -709, below exparg(1) ≈ -708.4.
        assert_eq!(fpser(1.0, 1e-20, 1e-308, eps), 0.0);

        // beta_pser at x == 0: explicit early return.
        assert_eq!(beta_pser(2.0, 3.0, 0.0, eps), 0.0);

        // apser exists and returns finite for small a, sensible b.
        let r = apser(1e-15, 5.0, 0.05, eps);
        assert!(r.is_finite() && r >= 0.0);
    }

    #[test]
    fn beta_grat_returns_ok_on_happy_path() {
        let eps = f64::EPSILON.max(1e-15);
        // Normal happy path through the 30-term expansion loop.
        let w = beta_grat(20.0, 0.5, 0.2, 0.8, 0.0, 15.0 * eps).unwrap();
        assert!(w.is_finite() && (0.0..=1.0).contains(&w));
    }

    #[test]
    fn beta_rcomp_a0_lt_1_unreached_via_beta_inc_but_safe() {
        // beta_inc calls beta_rcomp only through beta_frac, with a0 and b0
        // above 1, so the a0 < 1 paths of beta_rcomp are reached only by a
        // direct call. Verify those branches return a finite, non-negative
        // value at sensible inputs.
        // Path b0 ≥ 8 (large b, tiny a):
        let r = beta_rcomp(0.5, 30.0, 0.05, 0.95);
        assert!(r.is_finite() && r >= 0.0);
        // Path 1 < b0 < 8:
        let r = beta_rcomp(0.5, 5.0, 0.3, 0.7);
        assert!(r.is_finite() && r >= 0.0);
        // Path b0 ≤ 1:
        let r = beta_rcomp(0.5, 0.7, 0.5, 0.5);
        assert!(r.is_finite() && r >= 0.0);
    }

    #[test]
    fn beta_inc_nan_x_with_tiny_a_and_b_takes_label_260() {
        // max(a, b) < 0.001 eps (cdflib.f90:1092-1094): label 260 gives
        // w = b/(a+b) and w1 = a/(a+b) whatever x and y are; gfortran
        // returns (0.5, 0.5) here.
        assert_eq!(
            try_beta_inc(1e-20, 1e-20, f64::NAN, f64::NAN),
            Ok((0.5, 0.5))
        );
        assert_eq!(try_beta_inc(1e-20, 1e-20, 0.5, f64::NAN), Ok((0.5, 0.5)));
        let (w, w1) = try_beta_inc(200.0, 200.0, f64::NAN, f64::NAN).unwrap();
        assert!(w.is_nan() && w1.is_nan());
    }
}
