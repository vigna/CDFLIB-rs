# Change Log

## 0.4.4 - 2026-10-09

### Fixed

- The Fortran 90 reference is the current upstream `cdflib.f90`, verbatim,
  and the code cites its line numbers. It calls the library's own Γ function
  (`gamma_user`) in `gamma_inc`, `gamma_inc_inv`, and `rcomp`, as the original
  DCDFLIB does; the upstream Fortran 90 used to call the compiler's intrinsic
  instead (it has been fixed now, after our report). The reference tables are
  generated without fused multiply-adds, and the Rust port reproduces all of
  them bit for bit.

- Inputs that make the Fortran loop forever (NaN or infinite arguments of
  `gamma_inc`, `gamma_inc_inv`, `beta_pser`, `beta_grat`, `gamma_rat1`,
  `beta_frac`, `fpser`, `apser`) now return NaN.

- Branches taken on NaN in `dinvr`, `dzror`, `gamma_inc_inv`, and the
  kernels now follow the Fortran; `gamma_inc_inv` keeps the Fortran `w` at
  label 130.

- `gamma_inc` and `gamma_inc_inv` with a NaN argument follow the Fortran
  error checks instead of returning NaN, and `gamma_inc_inv` reports a NaN
  *p* or *q*, for which the Fortran returns a meaningless *x*, as
  `InconsistentPq`. `beta_inc` with a NaN *x* or *y* and tiny *a* and *b*
  returns the Fortran values of label 260.

- `FisherSnedecor` rejects degrees of freedom whose half is 0, where the
  Fortran ignores the error of `beta_inc` inside `cumf`.

- `FisherSnedecorNoncentral` checks *dfn* ≥ 1 and *dfd* ≥ 1 after the
  `cdffnc` status checks, as the Fortran does.

- The closed-form parameter searches return an error when the parameter
  they compute is not valid, where the Fortran returns it: rate 0 (at
  *p* = 0) or +∞ (at *q* = 0) in `Gamma::search_rate`, a mean of ±∞ (at
  *p* = 0 or 1) in `Normal::search_mean`, and a *σ* that is not positive
  (when no positive *σ* exists, when *x* = *μ*, and at *p* = 0 or 1) in
  `Normal::search_sd`.

- `Poisson::pmf` is 1 at 0 for *λ* = 0 instead of NaN, and the densities
  of the Γ, χ², Β and *F* distributions at the ends of the support are
  their limits (for example, the rate for Γ with shape 1) instead of 0.

- The noncentral χ² and *F* distributions panic where the Fortran default
  integers would overflow (*λ* beyond about 4.3 · 10⁹), instead of
  overflowing.

### Changed

- Every routine follows the Fortran 90 line by line: the `cum*` routines
  are crate-private functions with their Fortran names, the `cdf*` searches
  use the same `dstinv`/`dstzr` arguments and loops, argument checks run in
  the Fortran order, and the machine-constant routines `ipmpar` and
  `exparg` are ported.

- The integer quantile of the discrete distributions brackets the answer by
  doubling instead of using a heuristic starting point.

- The Fortran squares `x**2` are computed as the exact product `x * x`, as
  gfortran does, instead of with `powi`, whose precision Rust does not
  specify.

### Improved

- New Fortran 90 reference tables for every kernel, for the error exits of
  `gamma_user`, `psi`, `gamma_inc`, `gamma_inc_inv`, and `beta_inc`, for
  the call logs (results and status codes) of every `cdf*` routine, and for
  the iteration traces of `dinvr` and `dzror`. `tests/regenerate/coverage.sh`
  checks that they execute every line of the Fortran source except those
  listed with a reason in `tests/regenerate/unreachable.txt`.

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

- `NegativeBinomial::solve_trials` -> `NegativeBinomial::solve_r`.

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
