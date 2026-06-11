# Change Log

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
