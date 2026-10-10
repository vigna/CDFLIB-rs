# Change Log

## [0.4.4] - 2026-10-10

### Changed

- Every routine follows the Fortran 90 line by line: the `cum*` routines
  are crate-private functions with their Fortran names, the `cdf*` searches
  use the same `dstinv`/`dstzr` arguments and loops, argument checks run in
  the Fortran order, and the machine-constant routines `ipmpar` and
  `exparg` are ported. Because of the check order, a parameter of −∞ is
  now reported as not positive (or negative) rather than not finite, and
  some invalid combinations of arguments report a different error.

- The integer quantile of the discrete distributions brackets the answer by
  doubling instead of using a heuristic starting point.

- The Fortran squares `x**2` are computed as the exact product `x * x`, as
  gfortran does, instead of with `powi`, whose precision Rust does not
  specify.

- `SearchError::StartOutOfRange` displays its interval as `[a..b]`.

### Improved

- New Fortran 90 reference tables for every kernel, for the error exits of
  `gamma_user`, `psi`, `gamma_inc`, `gamma_inc_inv`, and `beta_inc`, for
  the call logs (results and status codes) of every `cdf*` routine,
  including the search failures and the status 10 that the public API can
  reach, and for the iteration traces of `dinvr` and `dzror`.
  `tests/regenerate/coverage.sh` checks that they execute every line of the
  Fortran source except those listed with a reason in
  `tests/regenerate/unreachable.txt`, and that they execute none of those.
  The tables also pin the documented results of the error-function, normal,
  `stvaln`, `dt1`, and Student's *t* routines at NaN, infinite, and huge
  arguments. Where the Fortran status is undefined, because `dzror` fails at
  its label 240 without assigning `qleft`, the Rust reports the answer as
  above the upper bound of the search, and tests pin this choice.
  `tests/regenerate/regenerate.sh` fails when a generator ends before writing
  its completion line, as after a Fortran `stop`, and the solver-trace test
  fails on a truncated trace.

- The documentation states the behavior of the special functions at NaN,
  infinite, and endpoint arguments, the precision limits of the noncentral
  distributions (whose series stop after a fixed number of terms, so that
  their values are badly wrong for a large noncentrality, and can exceed 1
  by a few ulps for a small one), the tolerances of the searches, which do
  not resolve answers below about 10⁻¹⁰, the searches that end without an
  error on a NaN function value, as in the Fortran (for example, for Β and
  *F* distributions with a parameter above about 2 · 10³⁰⁷, whose CDF is
  NaN), the searches for a probability below about 10⁻¹⁵⁸, which can end
  away from the answer or report the wrong bound, as in the Fortran, the
  cancellation in the densities, masses, and entropies for very large
  parameters, the rounding of integer arguments beyond 2⁵³, and how to build
  with Rust 1.71 to 1.76, for which the latest releases of `thiserror` are
  too recent. The examples in the API documentation check the values they
  compute.

- The continuous integration checks the minimum supported Rust version,
  the formatting and the documentation, and runs the tests in debug and
  release builds on Linux and, with bit-exact comparisons, on macOS on
  Apple silicon.

### Fixed

- The Fortran 90 reference is the current upstream `cdflib.f90`, verbatim,
  and the code cites its line numbers. It calls the library's own Γ function
  (`gamma_user`) in `gamma_inc`, `gamma_inc_inv`, and `rcomp`, as the original
  DCDFLIB does; the upstream Fortran 90 used to call the compiler's intrinsic
  instead (it has been fixed now, after our report). The reference tables are
  generated without fused multiply-adds, and the Rust port reproduces all of
  them bit for bit, in debug and release builds, on the platform that
  generates them (macOS on Apple silicon); elsewhere the tests allow the
  last-bit differences of the system libm.

- `psi` keeps the optimizer from merging the sine and cosine of its
  reflection formula into one `sincos` call, whose last bit can differ from
  the separate calls of the Fortran, so that release builds are bit-exact
  too.

- Inputs on which the Fortran loops forever now return NaN or an error:
  NaN arguments, and infinite arguments other than the tolerance, of
  `gamma_inc`, `beta_pser`, `beta_grat`, `gamma_rat1`, `beta_frac`, `fpser`,
  and `apser`, arguments of `gamma_rat1` beyond about 1.3 · 10¹⁵⁴, where its
  continued fraction overflows, and the rare finite arguments on which
  `beta_inc` never returns, give NaN; `gamma_inc_inv`
  gives NaN, or `NotConverged` when the Schröder iteration cannot proceed;
  `Gamma::search_shape` reports `AnswerAboveUpperBound` for *p* < 1 where
  *x* times the rate overflows. Several of these inputs never returned in
  0.4.3 either: for example, the `cdf` at +∞ of a χ² distribution with
  *df* ≤ 0.5, and these calls of `Gamma::search_shape`.

- Branches taken on NaN in `dinvr`, `dzror`, `gamma_inc_inv`, and the
  kernels now follow the Fortran; `gamma_inc_inv` keeps the Fortran `w` at
  label 130.

- `gamma_inc` and `gamma_inc_inv` check their other arguments in the
  Fortran order before handling a NaN argument, and `gamma_inc_inv`
  reports a NaN or negative *p* or *q*, for which the Fortran returns a
  meaningless *x*, an error, or never returns, as `InconsistentPq`.
  `beta_inc` with a NaN or infinite argument returns what the Fortran
  returns, whenever the Fortran returns.

- `try_gamma` reports `Underflow` instead of `Overflow` for negative
  arguments beyond about −171.6, where Γ underflows.

- `FisherSnedecor` rejects degrees of freedom whose half is 0, where the
  Fortran ignores the error of `beta_inc` inside `cumf`.

- `FisherSnedecorNoncentral` checks *dfn* ≥ 1 and *dfd* ≥ 1 after the
  `cdffnc` status checks, as the Fortran does.

- `Poisson` rejects a *λ* whose double overflows as `LambdaNotFinite`;
  `cumpoi` computes the χ² argument 2*λ*, and the Fortran returns NaN.

- The parameter searches return an error when the parameter they compute
  is not valid, where the Fortran returns it: rate 0 (at *p* = 0) or +∞
  (at *q* = 0, `RateNotFinite` instead of `GammaIncInv(AtInfinity)`) in
  `Gamma::search_rate`, a mean of ±∞ (at *p* = 0 or 1) in
  `Normal::search_mean`, a *σ* that is not positive (when no positive *σ*
  exists, when *x* = *μ*, and at *p* = 0 or 1) or not finite in
  `Normal::search_sd`, and a parameter of 0, at the lower end of the search
  interval, in `Gamma::search_shape`, `ChiSquared::search_df`,
  `Beta::search_a`, `Beta::search_b`, `ChiSquaredNoncentral::search_df`,
  `Binomial::search_trials`, `NegativeBinomial::search_r`, and
  `NegativeBinomial::search_pr`. `NegativeBinomial::search_pr` also
  rejects *r* = 0 as `RNotPositive`, as `NegativeBinomial::try_new` does.

- The noncentral χ² distribution at the smallest subnormal *x*, whose half
  is 0 and where the Fortran series gives NaN, sums the series at the next
  floating-point number, in `cdf`, `ccdf`, and the searches.

- `apser` uses the value 0 that the Fortran `psi` returns where it fails,
  which can happen only outside the domain of `apser`, instead of
  panicking.

- `Beta::search_a` with *x* = 0, `Beta::search_b` with *x* = 1, the `cdf`
  of the *t* and χ² distributions with the smallest subnormal *df*, and the
  `inverse_cdf` of a binomial distribution with `u64::MAX` trials no longer
  panic.

- `cdf` and `ccdf` are exactly 0 or 1 at ±∞ for the normal distribution,
  and at +∞ for the Γ, χ², noncentral χ² and noncentral *F*
  distributions, where the Fortran gives NaN or never returns (for the
  noncentral *F*, a sum truncated short of 1); for the normal and Γ
  distributions also where the standardized argument, or *x* times the
  rate, overflows. Where only *x* − *μ* overflows, the normal `cdf` and
  `ccdf` compute the standardized argument as *x*/*σ* − *μ*/*σ*, instead of
  returning NaN (with *μ* = −10³⁰⁸ and *σ* = `f64::MAX`, the `cdf` at
  10³⁰⁸ was NaN instead of 0.867).

- `inverse_ccdf` at *q* = 0 returns +∞ for the Poisson and negative
  binomial distributions, where the Fortran search stops at a finite value
  at which the upper tail is below its tolerance, and *n* for the binomial
  distribution, as the Fortran does, also where the Fortran search stops at
  its start value (*pr* = 0) or cannot start (*n* < 5); for the Poisson
  distribution with *λ* = 0 and the negative binomial distribution with
  *pr* = 1 the Fortran search stops at its start value too.

- The noncentral χ² and *F* distributions panic where the Fortran default
  integers would overflow (*λ* beyond about 4.3 · 10⁹), instead of
  overflowing, and the noncentral *F* distribution panics at once when its
  degrees of freedom are so large that its series sums to NaN, where the
  Fortran never returns.

- `Poisson::pmf` is 1 at 0 for *λ* = 0 and `NegativeBinomial::pmf` is 1 at
  0 for *pr* = 1, instead of NaN, and the densities of the Γ, χ², Β and *F*
  distributions at the ends of the support are their limits (for example,
  the rate for Γ with shape 1, +∞ for χ² with *df* < 2, and 0 at +∞)
  instead of 0 or NaN.

- The densities of the *t*, *F* and Β distributions use `ln_1p`, which
  keeps their precision for large *df*, *dfn*, *dfd*, or *b* (the
  *F*(1, 10¹⁷) density at 1 was 0.399, and the *F*(10¹⁶, 1) density at 1
  was 4.0 · 10⁷, instead of 0.242), and the *t* density no longer
  overflows for |*t*| beyond about 1.3 · 10¹⁵⁴, nor the *F* density where
  *dfn*/*dfd* or *dfn*·*x*/*dfd* overflows or underflows (the log-density
  of *F*(5, 2) at 10³⁰⁸ was −∞ instead of −1418.4). The masses of the
  binomial and negative binomial distributions compute ln(1 − *pr*) with
  `ln_1p` (the mass at 0 of a binomial distribution with *n* = 10¹² and
  *pr* = 10⁻¹³ had a relative error of 3 · 10⁻⁵).

- The means and variances of the Β, Γ, *F*, and noncentral *F*
  distributions no longer overflow or underflow for extreme parameters.

- The entropies of the χ², *t*, Β and *F* distributions are their limits
  for subnormal parameters instead of panicking or returning NaN (for the
  *F* distribution, unless both parameters are below about 10⁻³⁰⁸), and the
  entropy of the normal distribution no longer overflows or underflows for
  extreme *σ*.

- The error messages of `Poisson`, `FisherSnedecor`, and
  `NegativeBinomial` describe the condition they report, those of `Beta`
  and `NegativeBinomial` contain no backticks, those for a probability out
  of range say whether it is *p* or *q*, and those of the noncentral
  distributions and `Poisson` say “positive”, “nonnegative”, or “at least 1”
  instead of using comparison symbols. The error messages of the special
  functions write “nonnegative” and “nonpositive” instead of “non-negative”
  and “non-positive”.

## [0.4.3] - 2026-06-11

### Improved

- Lower tolerance margins, closer code, more Fortran 90 tests.

## [0.4.2] - 2026-05-24

### Improved

- More tests.

- Cleaned up code.

## [0.4.1] - 2026-05-21

### Changed

- A few further exact alignment to Fortran 90 code.

## [0.4.0] - 2026-05-21

### New

- `try_new` constructors return errors on wrong parameters,
  `new` constructors panic.

- `const` getters.

### Changed

- Error diagnostics has been made uniform.

- The parameter searches `solve_*` were renamed `search_*`, and
  `NegativeBinomial::solve_trials` became `NegativeBinomial::search_r`.

- `cdflib::special::GammaError` -> `cdflib::special::GammaDomainError`
  to avoid clash with `cdflib::GammaError`.

- `distribution` module was renamed `dist`.

- Lowered MSRV to 1.71 and switched to Rust 2021 edition.

- Shortened names of error variants.

- `inverse_sf` is no longer a trait method.

- "Survival function" -> "Complementary cumulative distribution function".

### Fixed

- Fixed erratic behavior at endpoints 0 and 1.

- Restored accuracy selection for `gamma_inc` via `gamma_inc_with_acc`.

- Solvers now take _P_ and _Q_ instead of computing _Q_ = 1 - _P_.

## [0.3.1] - 2026-05-20

### Fixed

- Examples about `statrs` were not correct.

## [0.3.0] - 2026-05-20

### Changed

- Major API surface change: fallible functions get `try_` prefix, and
  infallible variants panic. There are no longer silent errors returned as
  special values.

## [0.2.0] - 2026-05-20

### Changed

- `gamma_x` (the Γ function) has been renamed `gamma`.

## [0.1.1] - 2026-05-20

### Fixed

- Fixed repo name and a few mathematical typos.

## [0.1.0] - 2026-05-19

### New

- Initial release.
