//! Student's *t* special functions: [`dt1`] (cdflib.f90:8735), an
//! approximate inverse of the *t* CDF.
//!
//! [`dt1`]: dt1

use super::eval_pol;
use super::normal::dinvnr;

/// Computes an approximate inverse of the cumulative *t* distribution.
///
/// Returns the approximate value of *x* for which the *t* CDF with *df*
/// degrees of freedom has value *p*; *q* is 1 − *p*. This is CDFLIB's
/// `dt1` (cdflib.f90:8735), which `cdft` (cdflib.f90:6409-6412) uses as
/// the starting point of its search for *t*.
///
/// The result is NaN when *p* or *q* is 0. For small *df* it can have the
/// wrong sign (`dt1(0.6, 0.4, 0.1)` ≈ −36.3), be huge or ±∞, or be NaN
/// (`dt1(0.3, 0.7, 1e-300)`).
///
/// # Example
///
/// ```
/// use cdflib::special::dt1;
///
/// // At df = 10, p = 0.975 the exact t-quantile is ≈ 2.2281.
/// // dt1 is an approximation; for use as a starting value
/// // a few digits of accuracy suffice.
/// let t = dt1(0.975, 0.025, 10.0);
/// assert!((t - 2.228138851).abs() < 5e-3);
/// ```
#[inline]
#[allow(clippy::needless_late_init)]
pub fn dt1(p: f64, q: f64, df: f64) -> f64 {
    // The F90 coef(0:4,4) table (cdflib.f90:8776-8780): row i here is
    // column i there, the coefficients of a polynomial of degree ideg(i).
    const COEF: [[f64; 5]; 4] = [
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [3.0, 16.0, 5.0, 0.0, 0.0],
        [-15.0, 17.0, 19.0, 3.0, 0.0],
        [-945.0, -1920.0, 1482.0, 776.0, 79.0],
    ];
    const DENOM: [f64; 4] = [4.0, 96.0, 384.0, 92160.0];
    const IDEG: [usize; 4] = [1, 2, 3, 4];

    let x = dinvnr(p, q).abs();
    let xx = x * x;

    let mut sum1 = x;
    let mut denpow = 1.0;
    for i in 0..4 {
        let term = eval_pol(&COEF[i][..=IDEG[i]], xx) * x;
        denpow *= df;
        sum1 += term / (denpow * DENOM[i]);
    }

    let xp;
    if 0.5 <= p {
        xp = sum1;
    } else {
        xp = -sum1;
    }

    xp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dt1_at_median_is_essentially_zero() {
        // p = 0.5 ⇒ x = |dinvnr(0.5, 0.5)| ≈ 0 ⇒ all terms vanish.
        // dinvnr's Newton iteration leaves a sub-ULP residual (~8e-17)
        // rather than exact zero; the dt1.csv fixture records the
        // same residual, so the match against CDFLIB is exact.
        for &df in &[1.0_f64, 5.0, 50.0, 1.0e6] {
            let r = dt1(0.5, 0.5, df);
            assert!(r.abs() < 1e-15, "df = {df}: dt1(0.5, 0.5) = {r}");
        }
    }

    #[test]
    fn dt1_symmetric_around_median() {
        // dt1(p, 1-p, df) = -dt1(1-p, p, df) for any df > 0.
        for &df in &[1.0_f64, 5.0, 50.0] {
            for &p in &[0.1_f64, 0.3, 0.7, 0.9, 0.99] {
                let a = dt1(p, 1.0 - p, df);
                let b = dt1(1.0 - p, p, df);
                assert!(
                    (a + b).abs() < 1e-13,
                    "p={p}, df={df}: dt1(p,1-p)={a}, dt1(1-p,p)={b}, sum={}",
                    a + b
                );
            }
        }
    }

    #[test]
    fn dt1_approaches_normal_quantile_for_large_df() {
        // As df → ∞ the t distribution converges to the standard normal,
        // so dt1(p, q, df) → dinvnr(p, q). At df = 1e6 the leading
        // correction is O(1/df) ≈ 1e-6.
        for &(p, q) in &[(0.7_f64, 0.3), (0.9, 0.1), (0.975, 0.025), (0.999, 0.001)] {
            let z = dinvnr(p, q);
            let t = dt1(p, q, 1.0e6);
            assert!((t - z).abs() < 1e-4, "p={p}: t = {t}, z = {z}");
        }
    }

    #[test]
    fn dt1_known_value_at_df_10() {
        // Exact t-quantile at p = 0.975, df = 10: t ≈ 2.228138851... .
        // dt1 is accurate to ~3 digits at this df.
        let r = dt1(0.975, 0.025, 10.0);
        assert!((r - 2.228138851).abs() < 5e-3, "r = {r}");
    }

    #[test]
    fn dt1_known_value_at_df_5() {
        // Exact t-quantile at p = 0.975, df = 5: t ≈ 2.570581836... .
        let r = dt1(0.975, 0.025, 5.0);
        // Lower df means coarser approximation; ~1e-2 is plausible.
        assert!((r - 2.5705818356).abs() < 5e-2, "r = {r}");
    }

    #[test]
    fn dt1_sign_follows_p() {
        // p < 0.5 ⇒ result is negative; p > 0.5 ⇒ positive.
        assert!(dt1(0.3, 0.7, 5.0) < 0.0);
        assert!(dt1(0.7, 0.3, 5.0) > 0.0);
        assert!(dt1(0.001, 0.999, 5.0) < -1.0);
        assert!(dt1(0.999, 0.001, 5.0) > 1.0);
    }
}
