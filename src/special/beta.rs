//! Β function family: Β(*a*, *b*), ln Β(*a*, *b*), and the regularized
//! incomplete Β function *Iₓ*(*a*, *b*) and its complement.

#![allow(clippy::approx_constant, clippy::excessive_precision)]

use super::erf::error_fc_scaled;
use super::gamma::{alnrel, gam1, gamma_ln1, gamma_log, gsumln, psi, rexp, rlog1};

/// Largest negative argument to `exp` for which the result is nonzero in
/// IEEE 754 binary64; corresponds to CDFLIB's `exparg(1)`.
/// 0.99999 · (−1022) · 0.69314718055995 matching F90 cdflib.f90:9544, :9555.
const NEG_EXPARG: f64 = -708.389_334_568_083_540_9;
/// Largest positive argument to `exp`; F90's `exparg(0)`.
/// 0.99999 · 1024 · 0.69314718055995.
const POS_EXPARG: f64 = 709.775_615_066_259_888_4;

/// Returns exp(*μ* + *x*) where *μ* is a small integer scaling factor and *x* is a
/// double. Splits into pieces to avoid intermediate overflow.
#[inline]
pub fn esum(mu: i32, x: f64) -> f64 {
    if x <= 0.0 {
        if mu >= 0 {
            let w = mu as f64 + x;
            if w <= 0.0 {
                return w.exp();
            }
        }
    } else if mu <= 0 {
        let w = mu as f64 + x;
        if w >= 0.0 {
            return w.exp();
        }
    }
    (mu as f64).exp() * x.exp()
}

/// Returns ln(Γ(*b*) / Γ(*a* + *b*)) for *b* ≥ 8.
#[inline]
pub fn algdiv(a: f64, b: f64) -> f64 {
    const C0: f64 = 0.833333333333333e-1;
    const C1: f64 = -0.277777777760991e-2;
    const C2: f64 = 0.793650666825390e-3;
    const C3: f64 = -0.595202931351870e-3;
    const C4: f64 = 0.837308034031215e-3;
    const C5: f64 = -0.165322962780713e-2;

    let (c, x, d) = if b < a {
        let h = b / a;
        let c = 1.0 / (1.0 + h);
        let x = h / (1.0 + h);
        let d = a + (b - 0.5);
        (c, x, d)
    } else {
        let h = a / b;
        let c = h / (1.0 + h);
        let x = 1.0 / (1.0 + h);
        let d = b + (a - 0.5);
        (c, x, d)
    };

    let x2 = x * x;
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);

    let t = (1.0 / b).powi(2);
    let w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;
    let w = w * (c / b);

    let u = d * alnrel(a / b);
    let v = a * (b.ln() - 1.0);

    if v < u {
        w - v - u
    } else {
        w - u - v
    }
}

/// Returns ln Β(*a*, *b*) = ln Γ(*a*) + ln Γ(*b*) − ln Γ(*a* + *b*).
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
    let b = a0.max(b0);

    if a >= 8.0 {
        // Procedure for a ≥ 8.
        let w = bcorr(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -((a - 0.5) * c.ln());
        let v = b * alnrel(h);
        return if v < u {
            -(0.5 * b.ln()) + E + w - v - u
        } else {
            -(0.5 * b.ln()) + E + w - u - v
        };
    }

    if a < 1.0 {
        // a < 1
        if b < 8.0 {
            return gamma_log(a) + (gamma_log(b) - gamma_log(a + b));
        }
        return gamma_log(a) + algdiv(a, b);
    }

    // 1 ≤ a < 8
    if a <= 2.0 {
        if b <= 2.0 {
            return gamma_log(a) + gamma_log(b) - gsumln(a, b);
        }
        let w = 0.0;
        if b >= 8.0 {
            return gamma_log(a) + algdiv(a, b);
        }
        // fall through to S60 with w = 0
        return reduce_b(a, b, w);
    }

    // 2 < a < 8
    if b > 1000.0 {
        // S80: reduction of a
        let n = (a - 1.0) as i64;
        let mut w = 1.0;
        for _ in 0..n {
            a -= 1.0;
            w *= a / (1.0 + a / b);
        }
        return w.ln() - (n as f64) * b.ln() + (gamma_log(a) + algdiv(a, b));
    }

    // b ≤ 1000: reduce a
    let n = (a - 1.0) as i64;
    let mut w = 1.0;
    for _ in 0..n {
        a -= 1.0;
        let h = a / b;
        w *= h / (1.0 + h);
    }
    let w = w.ln();
    if b >= 8.0 {
        return w + gamma_log(a) + algdiv(a, b);
    }
    reduce_b(a, b, w)
}

fn reduce_b(a: f64, b0: f64, w: f64) -> f64 {
    let mut b = b0;
    let n = (b - 1.0) as i64;
    let mut z = 1.0;
    for _ in 0..n {
        b -= 1.0;
        z *= b / (a + b);
    }
    w + z.ln() + (gamma_log(a) + (gamma_log(b) - gsumln(a, b)))
}

/// Returns Δ(*a*) + Δ(*b*) − Δ(*a* + *b*), for *a* ≥ 8 and *b* ≥ 8.
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
    let s3 = 1.0 + (x + x2);
    let s5 = 1.0 + (x + x2 * s3);
    let s7 = 1.0 + (x + x2 * s5);
    let s9 = 1.0 + (x + x2 * s7);
    let s11 = 1.0 + (x + x2 * s9);
    let t = (1.0 / b).powi(2);
    let w = ((((C5 * s11 * t + C4 * s9) * t + C3 * s7) * t + C2 * s5) * t + C1 * s3) * t + C0;
    let w = w * (c / b);
    let t = (1.0 / a).powi(2);
    (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a + w
}

/// Returns Β(*a*, *b*) = Γ(*a*) Γ(*b*) / Γ(*a* + *b*).
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

/// Returns the Stirling remainder for the complete Β function:
/// ln Β(*a*, *b*) − [Stirling(*a*) + Stirling(*b*) − Stirling(*a* + *b*)],
/// where Stirling(*z*) = ln √(2π) + (*z* − ½) ln *z* − *z*.
///
/// Sums from smallest to largest argument for accuracy.
///
/// # Example
///
/// ```
/// use cdflib::special::internal::dbetrm;
///
/// // Stirling remainder is small and decreasing in (a, b) for large args.
/// let r = dbetrm(50.0, 60.0);
/// assert!(r.abs() < 0.01);
/// ```
#[inline]
pub fn dbetrm(a: f64, b: f64) -> f64 {
    use super::gamma::dstrem;
    let mut r = -dstrem(a + b);
    r += dstrem(a.max(b));
    r += dstrem(a.min(b));
    r
}

/// Returns *Iₓ*(*a*, *b*) when *b* < min(*ε*, *ε*·*a*) and *x* ≤ 0.5.
#[inline]
pub fn fpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    let mut result = 1.0;
    if a > 1e-3 * eps {
        // F90 cdflib.f90:9854-9863 has fpser = 0.0D+00 here before the
        // t = a * log(x) line, used as a dead-store before the t.exp()
        // overwrite or the return on underflow. Mirror it explicitly:
        result = 0.0;
        let t = a * x.ln();
        if t < NEG_EXPARG {
            return result;
        }
        result = t.exp();
    }
    result *= b / a;
    let tol = eps / a;
    let mut an = a + 1.0;
    let mut t = x;
    let mut s = t / an;
    loop {
        an += 1.0;
        t *= x;
        let c = t / an;
        s += c;
        if c.abs() <= tol {
            break;
        }
    }
    result * (1.0 + a * s)
}

/// Returns *I*₁ ₋ *ₓ*(*b*, *a*) when *a* is very small. Note the swapped
/// parameter convention: caller passes (*a*, *b*, *x*) where *a* is the
/// small parameter.
#[inline]
pub fn apser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    const G: f64 = 0.577215664901533;
    let bx = b * x;
    let mut t = x - bx;
    let c = if b * eps <= 2e-2 {
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
    }
    -(a * (c + s))
}

/// Returns *Iₓ*(*a*, *b*) by power series when *b* ≤ 1 or *b*·*x* ≤ 0.7.
#[inline]
pub fn beta_pser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }

    let a0 = a.min(b);
    let mut result;
    if a0 >= 1.0 {
        let z = a * x.ln() - beta_log(a, b);
        result = z.exp() / a;
    } else {
        let mut b0 = a.max(b);
        if b0 <= 1.0 {
            // a < 1, b ≤ 1
            result = x.powf(a);
            if result == 0.0 {
                return 0.0;
            }
            let apb = a + b;
            let z = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                let u = a + b - 1.0;
                (1.0 + gam1(u)) / apb
            };
            let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
            result = result * c * (b / apb);
        } else if b0 < 8.0 {
            // a < 1, 1 < b < 8
            let mut u = gamma_ln1(a0);
            let m = (b0 - 1.0) as i64;
            let mut c = 1.0;
            for _ in 1..=m {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u += c.ln();
            let z = a * x.ln() - u;
            b0 -= 1.0;
            let apb = a0 + b0;
            let t = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                let u = a0 + b0 - 1.0;
                (1.0 + gam1(u)) / apb
            };
            result = z.exp() * (a0 / a) * (1.0 + gam1(b0)) / t;
        } else {
            // a < 1, b ≥ 8
            let u = gamma_ln1(a0) + algdiv(a0, b0);
            let z = a * x.ln() - u;
            result = (a0 / a) * z.exp();
        }
    }

    if result == 0.0 || a <= 0.1 * eps {
        return result;
    }

    // Series.
    let mut sum = 0.0;
    let mut n = 0.0;
    let mut c = 1.0;
    let tol = eps / a;
    loop {
        n += 1.0;
        c = c * (0.5 + (0.5 - b / n)) * x;
        let w = c / (a + n);
        sum += w;
        if w.abs() <= tol {
            break;
        }
    }
    result * (1.0 + a * sum)
}

/// Returns *xᵃ* · *yᵇ* / Β(*a*, *b*).
#[inline]
pub fn beta_rcomp(a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST_VAL: f64 = 0.398942280401433; // 1/√(2π)
    if x == 0.0 || y == 0.0 {
        return 0.0;
    }
    let a0 = a.min(b);
    if a0 < 8.0 {
        let (lnx, lny) = if x <= 0.375 {
            (x.ln(), alnrel(-x))
        } else if y <= 0.375 {
            (alnrel(-y), y.ln())
        } else {
            (x.ln(), y.ln())
        };
        let z = a * lnx + b * lny;
        if a0 >= 1.0 {
            return (z - beta_log(a, b)).exp();
        }
        // Procedure for a < 1 or b < 1.
        let mut b0 = a.max(b);
        if b0 <= 1.0 {
            let result = z.exp();
            if result == 0.0 {
                return 0.0;
            }
            let apb = a + b;
            let z = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                let u = a + b - 1.0;
                (1.0 + gam1(u)) / apb
            };
            let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
            return result * (a0 * c) / (1.0 + a0 / b0);
        }
        if b0 < 8.0 {
            let mut u = gamma_ln1(a0);
            let n = (b0 - 1.0) as i64;
            let mut c = 1.0;
            for _ in 1..=n {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u += c.ln();
            let z = z - u;
            b0 -= 1.0;
            let apb = a0 + b0;
            let t = if apb <= 1.0 {
                1.0 + gam1(apb)
            } else {
                let u = a0 + b0 - 1.0;
                (1.0 + gam1(u)) / apb
            };
            return a0 * z.exp() * (1.0 + gam1(b0)) / t;
        }
        // 8 <= b0
        let u = gamma_ln1(a0) + algdiv(a0, b0);
        return a0 * (z - u).exp();
    }
    // a ≥ 8 and b ≥ 8.
    let (x0, y0, lambda) = if a <= b {
        let h = a / b;
        (h / (1.0 + h), 1.0 / (1.0 + h), a - (a + b) * x)
    } else {
        let h = b / a;
        (1.0 / (1.0 + h), h / (1.0 + h), (a + b) * y - b)
    };
    let e = -(lambda / a);
    let u = if e.abs() <= 0.6 {
        rlog1(e)
    } else {
        e - (x / x0).ln()
    };
    let e = lambda / b;
    // Use y0 directly, not 1.0 - x0. The two are mathematically equal but
    // 1.0 - x0 loses precision (down to exactly 0) when
    // h = min(a,b)/max(a,b) is below f64 epsilon, while y0 = h/(1+h)
    // preserves the small value. Matches F90's use of log(y/y0).
    let v = if e.abs() <= 0.6 {
        rlog1(e)
    } else {
        e - (y / y0).ln()
    };
    let z = (-(a * u + b * v)).exp();
    CONST_VAL * (b * x0).sqrt() * z * (-bcorr(a, b)).exp()
}

/// Returns exp(*μ*) · *xᵃ* · *yᵇ* / Β(*a*, *b*).
#[inline]
pub fn beta_rcomp1(mu: i32, a: f64, b: f64, x: f64, y: f64) -> f64 {
    const CONST_VAL: f64 = 0.398942280401433;
    let a0 = a.min(b);
    if a0 >= 8.0 {
        let (x0, y0, lambda) = if a <= b {
            let h = a / b;
            (h / (1.0 + h), 1.0 / (1.0 + h), a - (a + b) * x)
        } else {
            let h = b / a;
            (1.0 / (1.0 + h), h / (1.0 + h), (a + b) * y - b)
        };
        let e = -(lambda / a);
        let u = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (x / x0).ln()
        };
        let e = lambda / b;
        // Use y0 directly; see the comment in beta_rcomp for the rationale.
        let v = if e.abs() <= 0.6 {
            rlog1(e)
        } else {
            e - (y / y0).ln()
        };
        let t4 = -(a * u + b * v);
        let z = esum(mu, t4);
        return CONST_VAL * (b * x0).sqrt() * z * (-bcorr(a, b)).exp();
    }

    let (lnx, lny) = if x <= 0.375 {
        (x.ln(), alnrel(-x))
    } else if y <= 0.375 {
        (alnrel(-y), y.ln())
    } else {
        (x.ln(), y.ln())
    };
    let z = a * lnx + b * lny;
    if a0 >= 1.0 {
        return esum(mu, z - beta_log(a, b));
    }
    // Procedure for a < 1 or b < 1.
    let mut b0 = a.max(b);
    if b0 >= 8.0 {
        let u = gamma_ln1(a0) + algdiv(a0, b0);
        return a0 * esum(mu, z - u);
    }
    if b0 > 1.0 {
        // Algorithm for 1 < b0 < 8.
        let mut u = gamma_ln1(a0);
        let n = (b0 - 1.0) as i64;
        let mut c = 1.0;
        for _ in 1..=n {
            b0 -= 1.0;
            c *= b0 / (a0 + b0);
        }
        u += c.ln();
        let z = z - u;
        b0 -= 1.0;
        let apb = a0 + b0;
        let t = if apb <= 1.0 {
            1.0 + gam1(apb)
        } else {
            let u = a0 + b0 - 1.0;
            (1.0 + gam1(u)) / apb
        };
        return a0 * esum(mu, z) * (1.0 + gam1(b0)) / t;
    }
    // Algorithm for b0 ≤ 1.
    let result = esum(mu, z);
    if result == 0.0 {
        return 0.0;
    }
    let apb = a + b;
    let z = if apb <= 1.0 {
        1.0 + gam1(apb)
    } else {
        let u = a + b - 1.0;
        (1.0 + gam1(u)) / apb
    };
    let c = (1.0 + gam1(a)) * (1.0 + gam1(b)) / z;
    result * (a0 * c) / (1.0 + a0 / b0)
}

/// Returns *Iₓ*(*a*, *b*) − *Iₓ*(*a* + *n*, *b*) for positive integer *n*.
#[inline]
pub fn beta_up(a: f64, b: f64, x: f64, y: f64, n: i32, eps: f64) -> f64 {
    let apb = a + b;
    let ap1 = a + 1.0;
    let mut mu = 0;
    let mut d = 1.0;
    if n != 1 && a >= 1.0 && apb >= 1.1 * ap1 {
        // F90 (cdflib.f90:2267-2273): mu = abs(exparg(1)), k = exparg(0).
        // NEG_EXPARG = exparg(1) (negative bound), POS_EXPARG = exparg(0).
        mu = NEG_EXPARG.abs() as i32;
        let k = POS_EXPARG as i32;
        if k < mu {
            mu = k;
        }
        d = (-(mu as f64)).exp();
    }
    let bup = beta_rcomp1(mu, a, b, x, y) / a;
    if n == 1 || bup == 0.0 {
        return bup;
    }
    let nm1 = n - 1;
    let mut w = d;
    let mut k = 0_i32;
    if b > 1.0 {
        if y <= 1e-4 {
            k = nm1;
        } else {
            let r = (b - 1.0) * x / y - a;
            if r >= 1.0 {
                let t = nm1 as f64;
                k = nm1;
                if r < t {
                    k = r as i32;
                }
            }
        }
        // Add the increasing terms of the series.
        for i in 1..=k {
            let l = (i - 1) as f64;
            d *= ((apb + l) / (ap1 + l)) * x;
            w += d;
        }
    }
    // Add remaining terms.
    let kp1 = k + 1;
    for i in kp1..=nm1 {
        let l = (i - 1) as f64;
        d *= (apb + l) / (ap1 + l) * x;
        w += d;
        if d <= eps * w {
            break;
        }
    }
    bup * w
}

/// Returns the incomplete Γ ratios *P*(*a*, *x*), *Q*(*a*, *x*) specialized to
/// *a* ≤ 1. Used by [`beta_grat`].
///
/// [`beta_grat`]: crate::special::internal::beta_grat
#[inline]
pub fn gamma_rat1(a: f64, x: f64, r: f64, eps: f64) -> (f64, f64) {
    use super::erf::{error_f, error_fc};
    if a * x == 0.0 {
        return if x <= a { (0.0, 1.0) } else { (1.0, 0.0) };
    }
    if a == 0.5 {
        let rtx = x.sqrt();
        return if x < 0.25 {
            let p = error_f(rtx);
            (p, 0.5 + (0.5 - p))
        } else {
            let q = error_fc(rtx);
            (0.5 + (0.5 - q), q)
        };
    }

    if x < 1.1 {
        // Taylor series for P(a, x)/x^a.
        let mut an: f64 = 3.0;
        let mut c = x;
        let mut sum = x / (a + 3.0);
        let tol = 0.1 * eps / (a + 1.0);
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
        let use_label_50 = if x < 0.25 { z > -0.13394 } else { a < x / 2.59 };
        return if use_label_50 {
            let l = rexp(z);
            let w = 0.5 + (0.5 + l);
            let q = (w * j - l) * g - h;
            if q < 0.0 {
                (1.0, 0.0)
            } else {
                let p = 0.5 + (0.5 - q);
                (p, q)
            }
        } else {
            let w = z.exp();
            let p = w * g * (0.5 + (0.5 - j));
            let q = 0.5 + (0.5 - p);
            (p, q)
        };
    }

    // Continued fraction.
    let mut a2nm1: f64 = 1.0;
    let mut a2n: f64 = 1.0;
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
        if (an0 - am0).abs() < eps * an0 {
            let q = r * an0;
            return (0.5 + (0.5 - q), q);
        }
    }
}

/// Failure modes of [`beta_grat`].
///
/// All variants are *soft* failures: the routine cannot add a
/// correction to *w*, but the input *w* itself is a usable answer
/// (CDFLIB's documented fallback). Callers typically recover with
/// `result.unwrap_or(w_in)`.
///
/// [`beta_grat`]: crate::special::internal::beta_grat
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum BetaGratError {
    /// `b · z` evaluated to zero (b ≈ 0 or z = 0 with the other finite).
    #[error("b·z evaluated to zero")]
    BzZero,
    /// The exponential scale `u = r · exp(−u)` underflowed to zero.
    #[error("u underflowed to zero")]
    UnderflowedScale,
    /// The partial sum went non-positive during the 30-term expansion.
    #[error("partial sum went non-positive")]
    NonPositiveSum,
}

/// Returns *Iₓ*(*a*, *b*) by asymptotic expansion when 15 ≤ *a* and *b* ≤ 1.
///
/// Adds a correction to *w*; on success returns the updated value.
///
/// Each `Err` variant is a *soft* failure: the routine cannot add a
/// correction in that regime, but the input *w* itself is a usable
/// answer, the F90's documented fallback. Callers recover with
/// `result.unwrap_or(w_in)`.
#[inline]
pub fn beta_grat(
    a: f64,
    b: f64,
    x: f64,
    y: f64,
    w_in: f64,
    eps: f64,
) -> Result<f64, BetaGratError> {
    let bm1 = b - 0.5 - 0.5;
    let nu = a + 0.5 * bm1;
    let lnx = if y > 0.375 { x.ln() } else { alnrel(-y) };
    let z = -(nu * lnx);
    if b * z == 0.0 {
        return Err(BetaGratError::BzZero);
    }

    let mut r = b * (1.0 + gam1(b)) * (b * z.ln()).exp();
    r = r * (a * lnx).exp() * (0.5 * bm1 * lnx).exp();
    let mut u = algdiv(b, a) + b * nu.ln();
    u = r * (-u).exp();
    if u == 0.0 {
        return Err(BetaGratError::UnderflowedScale);
    }

    let (_, q) = gamma_rat1(b, z, r, eps);
    let v = 0.25 * (1.0 / nu).powi(2);
    let t2 = 0.25 * lnx * lnx;
    let l = w_in / u;
    let mut j = q / r;
    let mut sum = j;
    let mut t = 1.0;
    let mut cn = 1.0;
    let mut n2 = 0.0;
    let mut c_arr = [0.0_f64; 30];
    let mut d_arr = [0.0_f64; 30];

    for n in 1..=30 {
        let bp2n = b + n2;
        j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
        n2 += 2.0;
        t *= t2;
        cn /= n2 * (n2 + 1.0);
        c_arr[n - 1] = cn;
        let mut s = 0.0;
        if n != 1 {
            let mut coef = b - n as f64;
            for i in 1..n {
                s += coef * c_arr[i - 1] * d_arr[n - i - 1];
                coef += b;
            }
        }
        d_arr[n - 1] = bm1 * cn + s / (n as f64);
        let dj = d_arr[n - 1] * j;
        sum += dj;
        if sum <= 0.0 {
            return Err(BetaGratError::NonPositiveSum);
        }
        if dj.abs() <= eps * (sum + l) {
            return Ok(w_in + u * sum);
        }
    }
    Ok(w_in + u * sum)
}

/// Returns *Iₓ*(*a*, *b*) by asymptotic expansion when both *a* and *b* are ≥ 15.
#[inline]
pub fn beta_asym(a: f64, b: f64, lambda: f64, eps: f64) -> f64 {
    const E0: f64 = 1.12837916709551; // 2/√π
    const E1: f64 = 0.353553390593274; // 2^(-3/2)
    const NUM: usize = 20;

    let (h, r0, r1, w0) = if a < b {
        let h = a / b;
        let r0 = 1.0 / (1.0 + h);
        let r1 = (b - a) / b;
        let w0 = 1.0 / (a * (1.0 + h)).sqrt();
        (h, r0, r1, w0)
    } else {
        let h = b / a;
        let r0 = 1.0 / (1.0 + h);
        let r1 = (b - a) / a;
        let w0 = 1.0 / (b * (1.0 + h)).sqrt();
        (h, r0, r1, w0)
    };

    let f = a * rlog1(-(lambda / a)) + b * rlog1(lambda / b);
    let t = (-f).exp();
    if t == 0.0 {
        return 0.0;
    }
    let z0 = f.sqrt();
    let z = 0.5 * (z0 / E1);
    let z2 = f + f;

    let mut a0_arr = [0.0_f64; 21];
    let mut b0_arr = [0.0_f64; 21];
    let mut c_arr = [0.0_f64; 21];
    let mut d_arr = [0.0_f64; 21];

    a0_arr[0] = 2.0 / 3.0 * r1;
    c_arr[0] = -(0.5 * a0_arr[0]);
    d_arr[0] = -c_arr[0];
    let mut j0 = 0.5 / E0 * error_fc_scaled(z0);
    let mut j1 = E1;
    let mut sum = j0 + d_arr[0] * w0 * j1;
    let mut s = 1.0;
    let h2 = h * h;
    let mut hn = 1.0;
    let mut w = w0;
    let mut znm1 = z;
    let mut zn = z2;

    let mut n = 2;
    while n <= NUM {
        hn *= h2;
        a0_arr[n - 1] = 2.0 * r0 * (1.0 + h * hn) / (n as f64 + 2.0);
        let np1 = n + 1;
        s += hn;
        a0_arr[np1 - 1] = 2.0 * r1 * s / (n as f64 + 3.0);

        for i in n..=np1 {
            let r = -(0.5 * (i as f64 + 1.0));
            b0_arr[0] = r * a0_arr[0];
            for m in 2..=i {
                let mut bsum = 0.0;
                for j in 1..m {
                    let mmj = m - j;
                    bsum += (j as f64 * r - mmj as f64) * a0_arr[j - 1] * b0_arr[mmj - 1];
                }
                b0_arr[m - 1] = r * a0_arr[m - 1] + bsum / m as f64;
            }
            c_arr[i - 1] = b0_arr[i - 1] / (i as f64 + 1.0);
            let mut dsum = 0.0;
            for j in 1..i {
                let imj = i - j;
                dsum += d_arr[imj - 1] * c_arr[j - 1];
            }
            d_arr[i - 1] = -(dsum + c_arr[i - 1]);
        }
        j0 = E1 * znm1 + (n as f64 - 1.0) * j0;
        j1 = E1 * zn + (n as f64) * j1;
        znm1 *= z2;
        zn *= z2;
        w *= w0;
        let t0 = d_arr[n - 1] * w * j0;
        w *= w0;
        let t1 = d_arr[np1 - 1] * w * j1;
        sum += t0 + t1;
        if t0.abs() + t1.abs() <= eps * sum {
            break;
        }

        n += 2;
    }

    let u = (-bcorr(a, b)).exp();
    E0 * t * u * sum
}

/// Returns *Iₓ*(*a*, *b*) by continued fraction expansion when both
/// *a* and *b* are > 1.
#[inline]
pub fn beta_frac(a: f64, b: f64, x: f64, y: f64, lambda: f64, eps: f64) -> f64 {
    let bfrac_init = beta_rcomp(a, b, x, y);
    if bfrac_init == 0.0 {
        return 0.0;
    }
    let c = 1.0 + lambda;
    let c0 = b / a;
    let c1 = 1.0 + 1.0 / a;
    let yp1 = y + 1.0;
    let mut n = 0.0;
    let mut p = 1.0;
    let mut s = a + 1.0;
    let mut an = 0.0;
    let mut anp1 = 1.0;
    let mut bn = 1.0;
    let mut bnp1 = c / c1;
    let mut r = c1 / c;

    loop {
        n += 1.0;
        let t_local = n / a;
        let w = n * (b - n) * x;
        let e1 = a / s;
        let alpha = p * (p + c0) * e1 * e1 * (w * x);
        let e2 = (1.0 + t_local) / (c1 + t_local + t_local);
        let beta_v = n + w / s + e2 * (c + n * yp1);
        p = 1.0 + t_local;
        s += 2.0;

        let t_new = alpha * an + beta_v * anp1;
        an = anp1;
        anp1 = t_new;
        let t_new = alpha * bn + beta_v * bnp1;
        bn = bnp1;
        bnp1 = t_new;
        let r0 = r;
        r = anp1 / bnp1;

        if (r - r0).abs() <= eps * r {
            return bfrac_init * r;
        }

        // Rescale.
        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0;
    }
}

/// Errors of [`beta_inc`].
///
/// All variants correspond to invalid inputs: CDFLIB's `beta_inc` reports
/// them via a positive integer `ierr` [1 . . 7] and returns zeros for *w*,
/// *w*₁. This enum gives each one a named, descriptive form.
///
/// [`beta_inc`]: crate::special::beta_inc
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum BetaIncError {
    /// *a* or *b* is negative (CDFLIB `ierr = 1`).
    #[error("a or b is negative: a = {a}, b = {b}")]
    NegativeParameter { a: f64, b: f64 },
    /// Both *a* and *b* are zero (CDFLIB `ierr = 2`).
    #[error("both a and b are zero")]
    BothZero,
    /// *x* ∉ [0 . . 1] (CDFLIB `ierr = 3`).
    #[error("x must be in [0..1], got {0}")]
    XOutOfRange(f64),
    /// *y* ∉ [0 . . 1] (CDFLIB `ierr = 4`).
    #[error("y must be in [0..1], got {0}")]
    YOutOfRange(f64),
    /// *x* + *y* ≠ 1 within tolerance (CDFLIB `ierr = 5`).
    #[error("x + y must equal 1, got x = {x}, y = {y}")]
    InconsistentSum { x: f64, y: f64 },
    /// Degenerate: *x* = 0 and *a* = 0 (CDFLIB `ierr = 6`).
    #[error("degenerate: x = 0 and a = 0")]
    XZeroAndAZero,
    /// Degenerate: *y* = 0 and *b* = 0 (CDFLIB `ierr = 7`).
    #[error("degenerate: y = 0 and b = 0")]
    YZeroAndBZero,
}

/// Returns the regularized incomplete Β function *Iₓ*(*a*, *b*) and its
/// complement 1 − *Iₓ*(*a*, *b*).
///
/// The argument pair (*x*, *y*) is the (value, complement) of the
/// integration upper limit: the caller must supply *y* = 1 − *x* directly
/// rather than letting the routine subtract, because in the deep tail
/// the cancellation in `1.0 - x` would lose digits. The returned pair
/// (*w*, *w*₁) is the (lower-tail, upper-tail) probability with
/// *w* + *w*₁ = 1, analogous to the (*p*, *q*) pair returned by [`gamma_inc`];
/// both are computed independently rather than one from the other, for
/// the same precision reason.
///
/// # Panics
///
/// Panics on a [`BetaIncError`]. Use [`try_beta_inc`] for the fallible
/// form.
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

/// Fallible form of [`beta_inc`]: returns [`BetaIncError`] on invalid input.
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
#[allow(clippy::manual_range_contains)]
pub fn try_beta_inc(a: f64, b: f64, x: f64, y: f64) -> Result<(f64, f64), BetaIncError> {
    let eps = f64::EPSILON;
    if a < 0.0 || b < 0.0 {
        return Err(BetaIncError::NegativeParameter { a, b });
    }
    if a == 0.0 && b == 0.0 {
        return Err(BetaIncError::BothZero);
    }
    // Mirror CDFLIB's x < 0 || x > 1 form, not RangeInclusive::contains:
    // with NaN inputs, both comparisons return false, so NaN passes
    // through here and the x == 0 / y == 0 short-circuits below get a
    // chance to fire, matching CDFLIB's behavior for e.g. cumt with
    // extreme |t|.
    if x < 0.0 || x > 1.0 {
        return Err(BetaIncError::XOutOfRange(x));
    }
    if y < 0.0 || y > 1.0 {
        return Err(BetaIncError::YOutOfRange(y));
    }
    let z = x + y - 0.5 - 0.5;
    if z.abs() > 3.0 * eps {
        return Err(BetaIncError::InconsistentSum { x, y });
    }
    if x == 0.0 {
        if a == 0.0 {
            return Err(BetaIncError::XZeroAndAZero);
        }
        return Ok((0.0, 1.0));
    }
    if y == 0.0 {
        if b == 0.0 {
            return Err(BetaIncError::YZeroAndBZero);
        }
        return Ok((1.0, 0.0));
    }
    if a == 0.0 {
        return Ok((1.0, 0.0));
    }
    if b == 0.0 {
        return Ok((0.0, 1.0));
    }
    // A NaN that survives the short-circuits above would make the series
    // loops below diverge (their exit tolerances never compare true);
    // return the natural NaN instead.
    if a.is_nan() || b.is_nan() || x.is_nan() || y.is_nan() {
        return Ok((f64::NAN, f64::NAN));
    }

    let eps = eps.max(1e-15);
    if a.max(b) < 1e-3 * eps {
        return Ok((b / (a + b), a / (a + b)));
    }

    let mut ind = 0;
    let mut a0 = a;
    let mut b0 = b;
    let mut x0 = x;
    let mut y0 = y;

    let (w, w1) = if a0.min(b0) <= 1.0 {
        // Procedure for a0 ≤ 1 or b0 ≤ 1.
        if x > 0.5 {
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
        }
        // Note: variables now refer to (a0, b0, x0, y0).
        small_branch(a0, b0, x0, y0, eps)
    } else {
        // a0 > 1 and b0 > 1.
        let lambda = if a > b {
            (a + b) * y - b
        } else {
            a - (a + b) * x
        };
        let (lambda, ind_flip) = if lambda < 0.0 {
            (lambda.abs(), true)
        } else {
            (lambda, false)
        };
        if ind_flip {
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
        }
        large_branch(a0, b0, x0, y0, lambda, eps)
    };

    Ok(if ind == 0 { (w, w1) } else { (w1, w) })
}

fn small_branch(a0: f64, b0: f64, x0: f64, y0: f64, eps: f64) -> (f64, f64) {
    // S10 in CDFLIB.
    if b0 < eps.min(eps * a0) {
        let w = fpser(a0, b0, x0, eps);
        return (w, 0.5 + (0.5 - w));
    }
    if a0 < eps.min(eps * b0) && b0 * x0 <= 1.0 {
        let w1 = apser(a0, b0, x0, eps);
        return (0.5 + (0.5 - w1), w1);
    }
    if a0.max(b0) > 1.0 {
        // Falls into b0 > 1 path
        if b0 <= 1.0 {
            let w = beta_pser(a0, b0, x0, eps);
            return (w, 0.5 + (0.5 - w));
        }
        if x0 >= 0.3 {
            let w1 = beta_pser(b0, a0, y0, eps);
            return (0.5 + (0.5 - w1), w1);
        }
        if x0 < 0.1 && (x0 * b0).powf(a0) <= 0.7 {
            let w = beta_pser(a0, b0, x0, eps);
            return (w, 0.5 + (0.5 - w));
        }
        if b0 > 15.0 {
            // S150: beta_grat with b0 ≤ 1; but b0 > 15 here, so we have
            // to swap into beta_grat-on-(b0, a0, y0, x0) territory. The
            // CDFLIB code uses beta_up + beta_grat composition.
            let w1_grat = beta_grat(b0, a0, y0, x0, 0.0, 15.0 * eps).unwrap_or(0.0);
            let w = 0.5 + (0.5 - w1_grat);
            return (w, w1_grat);
        }
        let n = 20;
        let w1 = beta_up(b0, a0, y0, x0, n, eps);
        let b0_shifted = b0 + n as f64;
        let w1_total = beta_grat(b0_shifted, a0, y0, x0, w1, 15.0 * eps).unwrap_or(w1);
        return (0.5 + (0.5 - w1_total), w1_total);
    }
    // a0.max(b0) ≤ 1.
    if a0 >= 0.2_f64.min(b0) {
        let w = beta_pser(a0, b0, x0, eps);
        return (w, 0.5 + (0.5 - w));
    }
    if x0.powf(a0) <= 0.9 {
        let w = beta_pser(a0, b0, x0, eps);
        return (w, 0.5 + (0.5 - w));
    }
    if x0 >= 0.3 {
        let w1 = beta_pser(b0, a0, y0, eps);
        return (0.5 + (0.5 - w1), w1);
    }
    let n = 20;
    let w1 = beta_up(b0, a0, y0, x0, n, eps);
    let b0_shifted = b0 + n as f64;
    let w1_total = beta_grat(b0_shifted, a0, y0, x0, w1, 15.0 * eps).unwrap_or(w1);
    (0.5 + (0.5 - w1_total), w1_total)
}

fn large_branch(a0: f64, b0: f64, x0: f64, y0: f64, lambda: f64, eps: f64) -> (f64, f64) {
    // a0, b0 > 1.
    if b0 < 40.0 && b0 * x0 <= 0.7 {
        let w = beta_pser(a0, b0, x0, eps);
        return (w, 0.5 + (0.5 - w));
    }
    if b0 < 40.0 {
        // S160: reduce b0 to an integer + frac.
        let n = b0 as i32;
        let mut b0r = b0 - n as f64;
        let mut n_use = n;
        if b0r == 0.0 {
            n_use -= 1;
            b0r = 1.0;
        }
        let mut w = beta_up(b0r, a0, y0, x0, n_use, eps);
        if x0 <= 0.7 {
            w += beta_pser(a0, b0r, x0, eps);
            return (w, 0.5 + (0.5 - w));
        }
        let mut a0r = a0;
        if a0 <= 15.0 {
            let nn = 20;
            w += beta_up(a0r, b0r, x0, y0, nn, eps);
            a0r += nn as f64;
        }
        let w_total = beta_grat(a0r, b0r, x0, y0, w, 15.0 * eps).unwrap_or(w);
        return (w_total, 0.5 + (0.5 - w_total));
    }
    // b0 ≥ 40.
    if a0 <= b0 {
        if a0 <= 100.0 || lambda > 0.03 * a0 {
            let w = beta_frac(a0, b0, x0, y0, lambda, 15.0 * eps);
            return (w, 0.5 + (0.5 - w));
        }
        let w = beta_asym(a0, b0, lambda, 100.0 * eps);
        return (w, 0.5 + (0.5 - w));
    }
    if b0 <= 100.0 || lambda > 0.03 * b0 {
        let w = beta_frac(a0, b0, x0, y0, lambda, 15.0 * eps);
        return (w, 0.5 + (0.5 - w));
    }
    let w = beta_asym(a0, b0, lambda, 100.0 * eps);
    (w, 0.5 + (0.5 - w))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ============================================================ dbetrm

    // The residual at (100, 100) sits right at the native-FPU 1e-12 edge;
    // miri's soft-float ln pushes it over. Skipped under miri.
    #[cfg(not(miri))]
    #[test]
    fn dbetrm_matches_beta_log_minus_stirling() {
        // For each (a, b), dbetrm should equal ln Β(a, b) − Stirling decomposition.
        const HLN2PI: f64 = 0.91893853320467274178;
        fn stirling(z: f64) -> f64 {
            HLN2PI + (z - 0.5) * z.ln() - z
        }
        for &(a, b) in &[(2.5_f64, 3.5), (10.0, 20.0), (50.0, 60.0), (100.0, 100.0)] {
            let r = dbetrm(a, b);
            let lnb = beta_log(a, b);
            let stirling_sum = stirling(a) + stirling(b) - stirling(a + b);
            let expected = lnb - stirling_sum;
            assert!(
                (r - expected).abs() < 1e-12,
                "a={a}, b={b}: dbetrm={r}, expected={expected}"
            );
        }
    }

    #[test]
    fn dbetrm_decreases_for_large_args() {
        // The Stirling remainder shrinks as a, b grow.
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
        // I_{0.5}(a, a) = 0.5 by symmetry.
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
        // a.max(b) < 1e-3 * eps path → return (b/(a+b), a/(a+b)).
        let tiny = 1e-20;
        let (w, w1) = beta_inc(tiny, tiny, 0.5, 0.5);
        // Both ratios equal 0.5 by symmetry.
        assert!((w - 0.5).abs() < 1e-10);
        assert!((w1 - 0.5).abs() < 1e-10);
    }

    // ===== Regime switching points =====

    #[test]
    fn beta_inc_small_a_large_b_uses_grat_path() {
        // a0 ≤ 1, b0 > 15: small_branch's beta_grat composition.
        let (w, w1) = beta_inc(0.5, 30.0, 0.05, 0.95);
        // Sanity: numerically plausible.
        assert!(w > 0.0 && w < 1.0 && (w + w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn beta_inc_both_moderate_uses_frac_path() {
        // a, b both ≥ 8, b ≥ 40: large_branch → beta_frac → beta_rcomp's a0≥8 path.
        let (w, w1) = beta_inc(10.0, 60.0, 0.15, 0.85);
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
        // Swap territory: a > b, swap so b becomes the larger.
        let (w, w1) = beta_inc(60.0, 10.0, 0.85, 0.15);
        assert!((w + w1 - 1.0).abs() < 1e-10);
        // Symmetric to the previous: I_x(a,b) = 1 - I_{1-x}(b,a).
        let (w2, _) = beta_inc(10.0, 60.0, 0.15, 0.85);
        assert!((w - (1.0 - w2)).abs() < 1e-10);
    }

    #[test]
    fn beta_inc_extreme_lambda_asym_path() {
        // a, b ≥ 100 with lambda ≤ 0.03*a: triggers beta_asym.
        // Two parameter orderings to cover both branches of beta_asym.
        // a > b case (b0 > 100 after swap → uses else branch).
        let (w, w1) = beta_inc(150.0, 200.0, 150.0 / 350.0 + 0.001, 200.0 / 350.0 - 0.001);
        assert!((w + w1 - 1.0).abs() < 1e-8);
        // a < b case (a0 > 100, lambda < 0.03*a0).
        // mean = 150/550 ≈ 0.2727; x just below mean → lambda small positive.
        let (w, w1) = beta_inc(150.0, 400.0, 0.272, 0.728);
        assert!((w + w1 - 1.0).abs() < 1e-8);
    }

    #[test]
    fn beta_inc_a_below_eps_uses_apser() {
        // a < eps · max(a, b) AND b*x ≤ 1: triggers apser branch.
        let a = 1e-15;
        let b = 5.0;
        let (w, w1) = beta_inc(a, b, 0.05, 0.95);
        // For very small a, I_x(a, b) → 1.
        assert!(w > 0.99 && (w + w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn beta_inc_b_below_eps_uses_fpser() {
        // b < eps · max(a, b): triggers fpser branch.
        let a = 5.0;
        let b = 1e-15;
        let (w, w1) = beta_inc(a, b, 0.5, 0.5);
        // For very small b, I_x(a, b) → 0.
        assert!(w < 0.01 && (w + w1 - 1.0).abs() < 1e-10);
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
        // x = 0 path (covered via x <= 0).
        let r4 = esum(2, -1.0);
        assert!((r4 - (1.0_f64).exp()).abs() < 1e-14);
        // x <= 0, mu < 0: fallthrough.
        let r5 = esum(-1, -1.0);
        assert!((r5 - (-2.0_f64).exp()).abs() < 1e-14);
    }

    #[test]
    fn algdiv_in_b_le_a_branch() {
        // Force b <= a (the "swap" branch). The formula computes
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
        // At a = 41, b = 1e300, x = 0.8, y = 0.2 the intermediate 1 - x0
        // cancels to zero, so the result has to be computed via y0
        // directly to stay finite. Regression-guards that path.
        let r = beta_rcomp(41.0, 1e300, 0.8, 0.2);
        assert!(r.is_finite(), "beta_rcomp returned non-finite: {r}");
    }

    #[test]
    fn beta_inc_very_small_a_small_branch_corners() {
        // Lines 1022-1034 in small_branch require BOTH a0 ≤ 1 AND b0 ≤ 1
        // (so a0.max(b0) ≤ 1.0). Within that, three sub-branches based on
        // x and x^a:
        //   line 1022: x0.powf(a0) ≤ 0.9 → beta_pser
        //   line 1026: x0 ≥ 0.3 AND x0^a0 > 0.9 → beta_pser(b, a, y)
        //   line 1030: x0 < 0.3 AND x0^a0 > 0.9 → beta_up + beta_grat
        for &(a, b, x) in &[
            (0.1, 0.5, 0.3),  // x^a ≈ 0.887 ≤ 0.9 → line 1022
            (0.1, 0.5, 0.5),  // x^a ≈ 0.933 > 0.9, x ≥ 0.3 → line 1026
            (0.01, 0.5, 0.1), // x^a ≈ 0.977 > 0.9, x < 0.3 → line 1030
        ] {
            let (w, w1) = beta_inc(a, b, x, 1.0 - x);
            assert!((w + w1 - 1.0).abs() < 1e-10, "a={a}, b={b}, x={x}");
        }
    }

    #[test]
    fn beta_inc_large_branch_b_moderate_x_above_0_7() {
        // large_branch with b < 40 AND x > 0.7: triggers the a0 ≤ 15
        // beta_up + beta_grat path (lines 1052-1064).
        let (w, w1) = beta_inc(5.0, 10.0, 0.85, 0.15);
        assert!((w + w1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fpser_apser_beta_pser_edge_cases() {
        let eps = f64::EPSILON.max(1e-15);

        // fpser: a*ln(x) < NEG_EXPARG → return 0.
        // For a=1, x = 1e-308 (denormal) makes t ≈ -709, below NEG_EXPARG ≈ -708.4.
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
        // beta_inc routes a0 < 1 to small_branch (no beta_rcomp call).
        // beta_rcomp's a0 < 1 path is reachable only if a caller invokes
        // beta_rcomp directly. Verify those branches return a finite,
        // non-negative value at sensible inputs.
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
}
