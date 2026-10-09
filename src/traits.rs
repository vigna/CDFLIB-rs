//! Traits for distribution capabilities.
//!
//! The traits are deliberately small and focused so generic code can require
//! exactly the capability it needs. Every distribution in this crate
//! implements [`ContinuousCdf`] or [`DiscreteCdf`] plus, where applicable,
//! [`Continuous`] / [`Discrete`], [`Mean`], [`Variance`], and [`Entropy`].
//!
//! [`ContinuousCdf`]: crate::traits::ContinuousCdf
//! [`DiscreteCdf`]: crate::traits::DiscreteCdf
//! [`Continuous`]: crate::traits::Continuous
//! [`Discrete`]: crate::traits::Discrete
//! [`Mean`]: crate::traits::Mean
//! [`Variance`]: crate::traits::Variance
//! [`Entropy`]: crate::traits::Entropy

/// Cumulative distribution function (CDF), complementary CDF, and inverse CDF for
/// a continuous distribution.
///
/// # Example
///
/// ```
/// use cdflib::Normal;
/// use cdflib::traits::ContinuousCdf;
///
/// let n = Normal::new(0.0, 1.0);
/// let p = n.cdf(0.0);       // 0.5
/// let x = n.inverse_cdf(p).unwrap(); // 0.0
/// ```
pub trait ContinuousCdf {
    /// Domain-specific error type returned by the inverse routines.
    type Error;

    /// Returns Pr\[*X* ≤ *x*\].
    fn cdf(&self, x: f64) -> f64;

    /// Returns Pr\[*X* > *x*\], the complementary CDF.
    ///
    /// Implementations compute this independently of [`cdf`] rather than as
    /// `1 − cdf(x)`, so the small tail keeps its precision deep into the
    /// tails where the subtraction would lose digits to cancellation. The
    /// exceptions are [`ChiSquaredNoncentral`] and
    /// [`FisherSnedecorNoncentral`], for which CDFLIB computes
    /// `1 − cdf(x)`.
    ///
    /// [`cdf`]: ContinuousCdf::cdf
    /// [`ChiSquaredNoncentral`]: crate::ChiSquaredNoncentral
    /// [`FisherSnedecorNoncentral`]: crate::FisherSnedecorNoncentral
    fn ccdf(&self, x: f64) -> f64;

    /// Returns the smallest *x* such that [cdf]\(*x*\) ≥ *p*, for *p* ∈ [0 . . 1].
    ///
    /// At *p* = 0 returns the infimum of support, at *p* = 1 the supremum
    /// (either may be infinite).
    ///
    /// [cdf]: ContinuousCdf::cdf
    fn inverse_cdf(&self, p: f64) -> Result<f64, Self::Error>;
}

/// Cumulative distribution function (CDF), complementary CDF, and inverse CDF for
/// a discrete distribution over the non-negative integers.
///
/// # Example
///
/// ```
/// use cdflib::Poisson;
/// use cdflib::traits::DiscreteCdf;
///
/// let p = Poisson::new(3.0);
/// let c = p.cdf(2);
/// let s = p.inverse_cdf(c).unwrap(); // 2
/// ```
pub trait DiscreteCdf {
    /// Domain-specific error type returned by the inverse routines.
    type Error;

    /// Returns Pr\[*X* ≤ *x*\].
    fn cdf(&self, x: u64) -> f64;

    /// Returns Pr\[*X* > *x*\] = 1 − [cdf]\(*x*\).
    ///
    /// Required method: implementors must compute the upper tail
    /// independently from the lower tail rather than as `1.0 - cdf(x)`,
    /// so the small tail keeps its precision deep into the tails (where
    /// the subtraction would lose digits to cancellation).
    ///
    /// [cdf]: DiscreteCdf::cdf
    fn ccdf(&self, x: u64) -> f64;

    /// Returns the smallest integer *x* such that [cdf]\(*x*\) ≥ *p*, or
    /// the largest admissible *x* if there is none. At *p* = 0 returns 0;
    /// at *p* = 1 returns the largest admissible *x*: the number of trials
    /// for the binomial distribution, and [`u64::MAX`] for the unbounded
    /// ones, also when the parameters reduce the support to {0}.
    ///
    /// [cdf]: DiscreteCdf::cdf
    fn inverse_cdf(&self, p: f64) -> Result<u64, Self::Error>;
}

/// Probability density function (and its log) for a continuous distribution.
///
/// Implemented only when the density admits a closed-form expression. The
/// terms of such an expression can cancel when the parameters are very
/// large: the χ² density, for example, has a relative error of about
/// 10⁻⁵ at *df* = 10¹⁰.
pub trait Continuous {
    /// Returns the density *f*(*x*) of the distribution at *x*.
    fn pdf(&self, x: f64) -> f64;
    /// Returns the logarithm of the density *f*(*x*). Computing in log-space
    /// avoids underflow in the tails.
    fn ln_pdf(&self, x: f64) -> f64;
}

/// Probability mass function (and its log) for a discrete distribution.
pub trait Discrete {
    /// Returns the mass Pr\[*X* = *x*\] at the support point *x*.
    fn pmf(&self, x: u64) -> f64;
    /// Returns the logarithm of the mass Pr\[*X* = *x*\]. Computing in log-space
    /// avoids underflow in the tails.
    fn ln_pmf(&self, x: u64) -> f64;
}

/// First moment, AKA the mean.
pub trait Mean {
    /// Returns the expected value E\[*X*\].
    ///
    /// Returns NaN when the mean is not defined for the distribution's
    /// parameters.
    fn mean(&self) -> f64;
}

/// Second central moment, AKA the variance.
pub trait Variance {
    /// Returns the variance Var(*X*) = E\[(*X* − E\[*X*\])²\].
    ///
    /// Returns NaN when the variance is not defined for the distribution's
    /// parameters.
    fn variance(&self) -> f64;
    /// Returns the standard deviation (the square root of the [variance]).
    ///
    /// [variance]: Variance::variance
    #[inline]
    fn std_dev(&self) -> f64 {
        self.variance().sqrt()
    }
}

/// Differential entropy (for continuous distributions) or Shannon
/// entropy (for discrete distributions), in nats.
///
/// Implemented only when the entropy admits a closed-form expression. The
/// terms of such an expression cancel when the parameters are very large,
/// and the result loses its accuracy beyond about 10⁹.
pub trait Entropy {
    /// Returns the entropy of the distribution in nats.
    fn entropy(&self) -> f64;
}
