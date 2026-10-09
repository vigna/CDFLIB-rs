//! Root finders for the parameter searches of the `cdf*` routines.
//!
//! CDFLIB's `dinvr` (with its entry `dstinv`) and `dzror` (with its entry
//! `dstzr`) are reverse-communication state machines: the caller loops
//! while `status == 1`, evaluating the function at `x` and calling the
//! routine again with the value `fx`. The state machines live in
//! [`dinvr`] and [`dzror`]; this module wraps them in [`Dinvr`] and
//! [`Dzror`], whose methods mirror the F90 calls one to one. A `cdf*`
//! search written in F90 as
//!
//! ```text
//! call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )
//! status = 0
//! a = 5.0D+00
//! fx = 0.0D+00
//! call dinvr ( status, a, fx, qleft, qhi )
//! do while ( status == 1 )
//!   call cumbet ( x, y, a, b, cum, ccum )
//!   fx = cum - p
//!   call dinvr ( status, a, fx, qleft, qhi )
//! end do
//! ```
//!
//! reads in Rust as
//!
//! ```text
//! let mut d = dstinv(0.0, INF, 0.5, 0.5, 5.0, ATOL, TOL).dinvr(5.0)?;
//! while d.status() == 1 {
//!     let a = d.x();
//!     let (cum, ccum) = cumbet(x, y, a, b);
//!     let fx = cum - p;
//!     d.dinvr(fx);
//! }
//! ```
//!
//! [`Dinvr`]: self::Dinvr
//! [`Dzror`]: self::Dzror

mod dinvr;
mod dzror;

use crate::error::SearchError;
use dinvr::{InvrAction, InvrConfig, InvrState};
use dzror::{ZrorAction, ZrorConfig, ZrorState};

/// The arguments of F90 `dstinv` (cdflib.f90:8494): the search
/// interval [*small* . . *big*], the step parameters of the bracketing
/// phase, and the tolerances of the final `dzror` phase.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Dstinv {
    small: f64,
    big: f64,
    absstp: f64,
    relstp: f64,
    stpmul: f64,
    abstol: f64,
    reltol: f64,
}

/// F90 `dstinv`: sets the parameters of the next `dinvr` search.
#[inline]
pub(crate) fn dstinv(
    small: f64,
    big: f64,
    absstp: f64,
    relstp: f64,
    stpmul: f64,
    abstol: f64,
    reltol: f64,
) -> Dstinv {
    Dstinv {
        small,
        big,
        absstp,
        relstp,
        stpmul,
        abstol,
        reltol,
    }
}

impl Dstinv {
    /// F90 `status = 0; fx = 0; call dinvr ( status, x, fx, qleft, qhi )`:
    /// starts a `dinvr` search at *x*.
    ///
    /// Returns [`SearchError::StartOutOfRange`] where F90 `dinvr` stops
    /// with the fatal error "The values SMALL, X, BIG are not monotone"
    /// (cdflib.f90:8258-8263).
    #[inline]
    pub(crate) fn dinvr(self, x: f64) -> Result<Dinvr, SearchError> {
        // Rust only: the test of cdflib.f90:8258, made here so that the F90
        // fatal stop becomes an error; InvrState::step repeats it as a
        // panic, which this check makes unreachable.
        if !(self.small <= x && x <= self.big) {
            return Err(SearchError::StartOutOfRange {
                start: x,
                small: self.small,
                big: self.big,
            });
        }
        let cfg = InvrConfig {
            small: self.small,
            big: self.big,
            absstp: self.absstp,
            relstp: self.relstp,
            stpmul: self.stpmul,
            abstol: self.abstol,
            reltol: self.reltol,
        };
        let mut d = Dinvr {
            state: InvrState::new(cfg, x),
            status: 0,
            x,
            qleft: false,
            qhi: false,
        };
        d.dinvr(0.0);
        Ok(d)
    }
}

/// A `dinvr` search in progress, holding the F90 in/out arguments
/// `status`, `x`, `qleft` and `qhi`.
#[derive(Debug)]
pub(crate) struct Dinvr {
    state: InvrState,
    status: i32,
    x: f64,
    qleft: bool,
    qhi: bool,
}

impl Dinvr {
    /// F90 `call dinvr ( status, x, fx, qleft, qhi )` with *fx* the value
    /// of the function at [`x`](Self::x).
    #[inline]
    pub(crate) fn dinvr(&mut self, fx: f64) {
        match self.state.step(fx) {
            InvrAction::NeedEval(x) => {
                self.status = 1;
                self.x = x;
            }
            InvrAction::Converged(x) => {
                self.status = 0;
                self.x = x;
            }
            InvrAction::Failed { qleft, qhi, x } => {
                self.status = -1;
                self.qleft = qleft;
                self.qhi = qhi;
                self.x = x;
            }
        }
    }

    /// The F90 `status`: 1 if the function must be evaluated at
    /// [`x`](Self::x), 0 if [`x`](Self::x) is the answer, −1 if the search
    /// failed (see [`qleft`](Self::qleft) and [`qhi`](Self::qhi)).
    #[inline]
    pub(crate) fn status(&self) -> i32 {
        self.status
    }

    /// The F90 `x`: the next evaluation point, the answer, or on failure
    /// the value `dinvr` leaves in `x`.
    #[inline]
    pub(crate) fn x(&self) -> f64 {
        self.x
    }

    /// The F90 `qleft`: on failure, true if the stepping search terminated
    /// unsuccessfully at *small*, false if it terminated unsuccessfully at
    /// *big*.
    #[inline]
    pub(crate) fn qleft(&self) -> bool {
        self.qleft
    }

    /// The F90 `qhi`: on failure, true if *Y* < F(*X*) at the termination
    /// of the search, false if F(*X*) < *Y*.
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn qhi(&self) -> bool {
        self.qhi
    }
}

/// The arguments of F90 `dstzr` (cdflib.f90:9130): the interval
/// [*xlo* . . *xhi*] and the tolerances of a `dzror` search.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Dstzr {
    xlo: f64,
    xhi: f64,
    abstol: f64,
    reltol: f64,
}

/// F90 `dstzr`: sets the parameters of the next `dzror` search.
#[inline]
pub(crate) fn dstzr(xlo: f64, xhi: f64, abstol: f64, reltol: f64) -> Dstzr {
    Dstzr {
        xlo,
        xhi,
        abstol,
        reltol,
    }
}

impl Dstzr {
    /// F90 `status = 0; fx = 0; call dzror ( status, x, fx, xlo, xhi,
    /// qleft, qhi )`: starts a `dzror` search.
    #[inline]
    pub(crate) fn dzror(self) -> Dzror {
        let cfg = ZrorConfig {
            xlo: self.xlo,
            xhi: self.xhi,
            abstol: self.abstol,
            reltol: self.reltol,
        };
        let mut z = Dzror {
            state: ZrorState::new(cfg),
            status: 0,
            x: 0.0,
            qleft: false,
            qhi: false,
        };
        z.dzror(0.0);
        z
    }
}

/// A `dzror` search in progress, holding the F90 in/out arguments
/// `status`, `x`, `qleft` and `qhi`.
#[derive(Debug)]
pub(crate) struct Dzror {
    state: ZrorState,
    status: i32,
    x: f64,
    qleft: bool,
    qhi: bool,
}

impl Dzror {
    /// F90 `call dzror ( status, x, fx, xlo, xhi, qleft, qhi )` with *fx*
    /// the value of the function at [`x`](Self::x).
    #[inline]
    pub(crate) fn dzror(&mut self, fx: f64) {
        match self.state.step(fx) {
            ZrorAction::NeedEval(x) => {
                self.status = 1;
                self.x = x;
            }
            ZrorAction::Converged { x, .. } => {
                self.status = 0;
                self.x = x;
            }
            ZrorAction::Failed { qleft, qhi, .. } => {
                self.status = -1;
                self.qleft = qleft;
                self.qhi = qhi;
            }
        }
    }

    /// The F90 `status`: 1 if the function must be evaluated at
    /// [`x`](Self::x), 0 if [`x`](Self::x) is the answer, −1 if the search
    /// failed (see [`qleft`](Self::qleft) and [`qhi`](Self::qhi)).
    #[inline]
    pub(crate) fn status(&self) -> i32 {
        self.status
    }

    /// The F90 `x`: the next evaluation point, or the answer.
    #[inline]
    pub(crate) fn x(&self) -> f64 {
        self.x
    }

    /// The F90 `qleft`: on failure, true if the search terminated
    /// unsuccessfully at *xlo*, false if it terminated unsuccessfully at
    /// *xhi*.
    #[inline]
    pub(crate) fn qleft(&self) -> bool {
        self.qleft
    }

    /// The F90 `qhi`: on failure, true if *Y* < F(*X*) at the termination
    /// of the search, false if F(*X*) < *Y*.
    #[inline]
    #[allow(dead_code)]
    pub(crate) fn qhi(&self) -> bool {
        self.qhi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Drives a dinvr search on [small..big] from x, returning the answer or
    // the qleft flag of a failed search.
    fn solve(small: f64, big: f64, x: f64, f: impl Fn(f64) -> f64) -> Result<f64, bool> {
        let mut d = dstinv(small, big, 0.5, 0.5, 5.0, 1.0e-10, 1.0e-8)
            .dinvr(x)
            .unwrap();
        while d.status() == 1 {
            let fx = f(d.x());
            d.dinvr(fx);
        }
        if d.status() == -1 {
            return Err(d.qleft());
        }
        Ok(d.x())
    }

    #[test]
    fn solves_increasing_function() {
        // f(x) = x³ - 8; root at x = 2.
        let r = solve(0.0, 100.0, 1.0, |x| x.powi(3) - 8.0).unwrap();
        assert!((r - 2.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn solves_decreasing_function() {
        // f(x) = 1/x - 0.25; root at x = 4.
        let r = solve(0.01, 1000.0, 10.0, |x| 1.0 / x - 0.25).unwrap();
        assert!((r - 4.0).abs() < 1e-10, "r = {r}");
    }

    #[test]
    fn solves_root_at_moderate_value() {
        // f(x) = ln(x) - 1, root at x = e.
        let r = solve(1e-10, 1000.0, 1.0, |x| x.ln() - 1.0).unwrap();
        let e = std::f64::consts::E;
        assert!((r - e).abs() / e < 1e-8, "r = {r}, e = {e}");
    }

    #[test]
    fn reports_qleft_and_qhi_failures() {
        // Increasing, already positive at small: fails at small.
        assert_eq!(solve(1.0, 10.0, 5.0, |x| x + 1.0), Err(true));
        // Increasing, still negative at big: fails at big.
        assert_eq!(solve(1.0, 10.0, 5.0, |x| x - 100.0), Err(false));
        // Decreasing, already negative at small: fails at small.
        assert_eq!(solve(1.0, 10.0, 5.0, |x| -x - 1.0), Err(true));
        // Decreasing, still positive at big: fails at big.
        assert_eq!(solve(1.0, 10.0, 5.0, |x| 100.0 - x), Err(false));
    }

    #[test]
    fn converges_immediately_when_start_is_the_root() {
        let r = solve(0.0, 10.0, 3.0, |x| x - 3.0).unwrap();
        assert!((r - 3.0).abs() < 1e-15);
    }

    #[test]
    fn nan_objective_fails_at_small() {
        // As in F90, a NaN function value takes neither branch of the
        // monotonicity test and the stepping search fails at small.
        assert_eq!(solve(0.0, 1.0, 0.5, |_| f64::NAN), Err(true));
    }

    #[test]
    fn start_outside_range_is_the_f90_fatal_error() {
        let err = dstinv(0.0, 4.0, 0.5, 0.5, 5.0, 1.0e-10, 1.0e-8)
            .dinvr(5.0)
            .unwrap_err();
        assert_eq!(
            err,
            SearchError::StartOutOfRange {
                start: 5.0,
                small: 0.0,
                big: 4.0
            }
        );
    }

    #[test]
    fn dzror_finds_bracketed_zero_and_reports_failures() {
        let run = |f: &dyn Fn(f64) -> f64| {
            let mut z = dstzr(0.0, 1.0, 1.0e-10, 1.0e-8).dzror();
            while z.status() == 1 {
                let fx = f(z.x());
                z.dzror(fx);
            }
            (z.status(), z.x(), z.qleft())
        };
        let (status, x, _) = run(&|x| x * x - 0.25);
        assert_eq!(status, 0);
        assert!((x - 0.5).abs() < 1e-8, "x = {x}");
        // Same sign at both ends.
        let (status, _, _) = run(&|x| x + 1.0);
        assert_eq!(status, -1);
    }

    // Replays tests/data/solver_traces.txt, written by
    // tests/regenerate/gen_solver_traces.f90: each case feeds the recorded
    // function values to Dinvr or Dzror and checks every requested x and
    // the outcome bit for bit against cdflib.f90.
    #[test]
    #[cfg(not(miri))]
    fn solver_traces_match_f90() {
        let text = std::fs::read_to_string("tests/data/solver_traces.txt").unwrap();
        let hex = |t: &str| f64::from_bits(u64::from_str_radix(t, 16).unwrap());
        let flag = |t: &str| t == "T";
        let mut lines = text.lines().peekable();
        let mut cases = 0;
        while let Some(head) = lines.next() {
            let h: Vec<&str> = head.split_whitespace().collect();
            // Status, x, qleft and qhi after each call, for either solver.
            let (mut status, mut x, mut qleft, mut qhi);
            let mut dinvr = None;
            let mut dzror = None;
            match h[0] {
                "dinvr" => {
                    let d = dstinv(hex(h[1]), hex(h[2]), 0.5, 0.5, 5.0, hex(h[4]), hex(h[5]))
                        .dinvr(hex(h[3]))
                        .unwrap();
                    (status, x, qleft, qhi) = (d.status(), d.x(), d.qleft(), d.qhi());
                    dinvr = Some(d);
                }
                "dzror" => {
                    let z = dstzr(hex(h[1]), hex(h[2]), hex(h[3]), hex(h[4])).dzror();
                    (status, x, qleft, qhi) = (z.status(), z.x(), z.qleft(), z.qhi());
                    dzror = Some(z);
                }
                _ => panic!("unexpected line {head}"),
            }
            cases += 1;
            while let Some(line) = lines.next_if(|l| !l.starts_with('d')) {
                let t: Vec<&str> = line.split_whitespace().collect();
                match t[0] {
                    "need" => {
                        assert_eq!(status, 1, "{head}: {line}");
                        assert_eq!(x.to_bits(), hex(t[1]).to_bits(), "{head}: {line}");
                        let fx = hex(t[2]);
                        if let Some(d) = dinvr.as_mut() {
                            d.dinvr(fx);
                            (status, x, qleft, qhi) = (d.status(), d.x(), d.qleft(), d.qhi());
                        } else if let Some(z) = dzror.as_mut() {
                            z.dzror(fx);
                            (status, x, qleft, qhi) = (z.status(), z.x(), z.qleft(), z.qhi());
                        }
                    }
                    "conv" => {
                        assert_eq!(status, 0, "{head}: {line}");
                        assert_eq!(x.to_bits(), hex(t[1]).to_bits(), "{head}: {line}");
                    }
                    "fail" => {
                        assert_eq!(status, -1, "{head}: {line}");
                        assert_eq!(qleft, flag(t[1]), "{head}: {line}");
                        // F90 leaves qleft and qhi unset when dzror fails at
                        // label 240; the Rust state machine reports false.
                        if dinvr.is_some() {
                            assert_eq!(qhi, flag(t[2]), "{head}: {line}");
                        }
                        assert_eq!(x.to_bits(), hex(t[3]).to_bits(), "{head}: {line}");
                    }
                    "cap" => {}
                    _ => panic!("unexpected line {line}"),
                }
            }
        }
        assert!(cases >= 150, "{cases}");
    }
}
