# Change Log

## [0.4.4] - 2026-10-10

### Changed

- Every routine follows the Fortran 90 line by line, including the order of
  the argument checks, so some invalid inputs (e.g., a parameter of −∞)
  report a different error.

- The integer quantile of the discrete distributions brackets the answer by
  doubling.

- Error messages are more descriptive and print floating-point values in
  short form.

### Improved

- Fortran 90 reference tables for every kernel, error exit, `cdf*` call and
  status, and for the iteration traces of `dinvr` and `dzror`;
  `tests/regenerate/coverage.sh` checks that they execute every reachable
  line of `cdflib.f90`.

- The documentation states the behavior at NaN, infinite and endpoint
  arguments, the tolerances of the searches, and the precision limits of
  the noncentral distributions and of the closed-form densities and masses.

- CI also runs Miri, the tests with the MSRV, semver checks, the Fortran
  coverage check, and the regeneration of the reference tables.

### Fixed

- The reference is the current upstream `cdflib.f90`, which calls
  `gamma_user` instead of the intrinsic Γ (fixed upstream after our report).
  The port reproduces every reference table bit for bit, in debug and
  release builds, on macOS on Apple silicon.

- Inputs on which the Fortran loops forever return NaN or an error, and the
  noncentral distributions panic where the Fortran integers would overflow.

- Constructors and searches reject parameters on which the Fortran returns
  meaningless values, and several inputs no longer panic.

- `cdf`, `ccdf` and the inverses are exact at the ends of the support, and
  a NaN argument gives NaN.

- Densities, masses, moments and entropies no longer overflow, underflow or
  lose precision for extreme parameters.

- `try_gamma` reports `Underflow` instead of `Overflow` below about −171.6.

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
