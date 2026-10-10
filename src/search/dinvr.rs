//! The F90 subroutine `dinvr` and its entry `dstinv` (cdflib.f90:8135-8583).
//!
//! `dinvr` bounds the zero of the function and invokes `dzror`. It is
//! given a monotone function F and seeks an argument value *X* such that
//! F(*X*) = *Y*. It compares F(*X*) with *Y* for the input value of *X*,
//! then uses `qincr` to determine whether to step left or right to bound
//! the desired *X*; the initial step size is max(`absstp`, `relstp` ·
//! |*X*|). It iteratively steps right or left until it bounds *X*,
//! multiplying the step size by `stpmul` at each step which does not bound
//! *X*, and never steps beyond `small` or `big`. If *X* is successfully
//! bounded, Algorithm R of Bus and Dekker ([`ZrorState`]) is employed to
//! find the zero of F(*X*) − *Y*.
//!
//! The F90 keeps its state in `save` variables and resumes at the label
//! stored in `i99999`; here the saved variables are the fields of
//! [`InvrState`], with the F90 names, and the resume label is a private
//! enum.
//!
//! [`ZrorState`]: super::dzror::ZrorState

use super::dzror::{ZrorAction, ZrorConfig, ZrorState};

/// The arguments of the F90 entry `dstinv` (cdflib.f90:8494).
///
/// The fields carry the names of the `save` variables that `dstinv`
/// sets (cdflib.f90:8574-8580), not those of its dummy arguments
/// (`zsmall`, `zbig`, …). `small` and `big` are the left and right
/// endpoints of the interval to be searched for a solution; the initial
/// step size in the search is max(`absstp`, `relstp` · |*X*|); when a
/// step does not bound the zero, the step size is multiplied by `stpmul`
/// and another step taken; `abstol` and `reltol` are two numbers that
/// determine the accuracy of the solution.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InvrConfig {
    pub small: f64,
    pub big: f64,
    pub absstp: f64,
    pub relstp: f64,
    pub stpmul: f64,
    pub abstol: f64,
    pub reltol: f64,
}

/// The F90 resume label `i99999` (cdflib.f90:8235-8254).
#[derive(Debug, Clone, Copy)]
enum Resume {
    /// No evaluation pending: the F90 call with `status` = 0.
    Entry,
    /// Resume at label 10 with F(`small`).
    Label10,
    /// Resume at label 20 with F(`big`).
    Label20,
    /// Resume at label 90 with F(`xsave`).
    Label90,
    /// Resume at label 130 with F(`xub`).
    Label130,
    /// Resume at label 200 with F(`xlb`).
    Label200,
    /// Resume at label 270 with the value requested by `dzror`.
    Label270,
}

/// The F90 `status` on return from `dinvr`.
#[derive(Debug, Clone, Copy)]
pub(crate) enum InvrAction {
    /// `status` = 1: the function must be evaluated at *x*, and the value
    /// passed as `fx` to the next [`InvrState::step`] call.
    NeedEval(f64),
    /// `status` = 0: *x* is an approximate root of F(*X*).
    Converged(f64),
    /// `status` = -1: the routine cannot bound the function. `qleft` is
    /// true if the stepping search terminated unsuccessfully at `small`,
    /// and false if the search terminated unsuccessfully at `big`; `qhi`
    /// is true if *Y* < F(*X*) at the termination of the search and false
    /// if F(*X*) < *Y*; `x` is the F90 `x` on return.
    Failed { qleft: bool, qhi: bool, x: f64 },
}

/// The `save` variables of the F90 `dinvr`, with the F90 names, together
/// with the dummy argument `x`, which the F90 caller passes back unchanged
/// on every call.
#[derive(Debug)]
pub(crate) struct InvrState {
    // Set by the entry dstinv.
    small: f64,
    big: f64,
    absstp: f64,
    relstp: f64,
    stpmul: f64,
    abstol: f64,
    reltol: f64,
    i99999: Resume,
    // The reverse-communication variable: the caller's x on the first
    // call, then the point handed out for evaluation.
    x: f64,
    xsave: f64,
    fsmall: f64,
    fbig: f64,
    qincr: bool,
    step: f64,
    xlb: f64,
    xub: f64,
    // The save variables of dzror, set by the dstzr call at label 240.
    zror: ZrorState,
}

impl InvrState {
    /// The F90 entry `dstinv` (cdflib.f90:8574-8582), with the starting
    /// point `x` the caller passes to its first `dinvr` call. That call,
    /// with `status` = 0, is the first [`step`].
    ///
    /// [`step`]: Self::step
    #[inline]
    pub(crate) fn new(cfg: InvrConfig, x: f64) -> Self {
        Self {
            small: cfg.small,
            big: cfg.big,
            absstp: cfg.absstp,
            relstp: cfg.relstp,
            stpmul: cfg.stpmul,
            abstol: cfg.abstol,
            reltol: cfg.reltol,
            i99999: Resume::Entry,
            x,
            xsave: 0.0,
            fsmall: 0.0,
            fbig: 0.0,
            qincr: false,
            step: 0.0,
            xlb: 0.0,
            xub: 0.0,
            // Rust only: a placeholder, since the field needs a value before
            // the dstzr call at label 240 sets it; no dzror step runs before
            // label 240. Its xxlo = xxhi = 0 are the F90 initial values
            // (cdflib.f90:8919-8920); abstol and reltol have none in the F90.
            zror: ZrorState::new(ZrorConfig {
                xlo: 0.0,
                xhi: 0.0,
                abstol: 0.0,
                reltol: 0.0,
            }),
        }
    }

    /// One call of the F90 `dinvr` (cdflib.f90:8135).
    ///
    /// The first call is the F90 call with `status` = 0, and `fx` is
    /// ignored. If the routine needs the function to be evaluated, it
    /// returns [`InvrAction::NeedEval`] (`status` = 1); the value of the
    /// function must be passed as `fx` to the next call. If the routine
    /// finishes without error, it returns [`InvrAction::Converged`]
    /// (`status` = 0) with an approximate root of F(*X*). If it cannot
    /// bound the function, it returns [`InvrAction::Failed`]
    /// (`status` = -1) with `qleft` and `qhi`.
    ///
    /// # Panics
    ///
    /// If the starting point is not in [*small* . . *big*], where the F90
    /// stops with a fatal error.
    #[allow(clippy::assign_op_pattern)]
    #[inline]
    pub(crate) fn step(&mut self, fx: f64) -> InvrAction {
        // cdflib.f90:8235-8254: if 0 < status, resume at label i99999.
        // Rust only: the F90 stop 1 for an illegal i99999 cannot occur, as
        // Resume has no other values.
        match self.i99999 {
            Resume::Entry => {
                // cdflib.f90:8256 also stores this test in qcond, which
                // label 130 or label 200 overwrites before label 110 or
                // label 180 reads it.
                if !(self.small <= self.x && self.x <= self.big) {
                    // F90 stop 1 (cdflib.f90:8258-8263). Dstinv::dinvr in
                    // src/search/mod.rs rejects this case with
                    // SearchError::StartOutOfRange before driving the state.
                    panic!("dinvr: the values small, x, big are not monotone");
                }

                self.xsave = self.x;
                // See that SMALL and BIG bound the zero and set QINCR.
                self.x = self.small;
                // GET-function-VALUE
                self.i99999 = Resume::Label10;
                InvrAction::NeedEval(self.x)
            }
            Resume::Label10 => {
                // Label 10.
                self.fsmall = fx;
                self.x = self.big;
                // GET-function-VALUE
                self.i99999 = Resume::Label20;
                InvrAction::NeedEval(self.x)
            }
            Resume::Label20 => {
                // Label 20.
                self.fbig = fx;

                self.qincr = self.fsmall < self.fbig;

                if self.fsmall <= self.fbig {
                    if 0.0 < self.fsmall {
                        // status = -1
                        return InvrAction::Failed {
                            qleft: true,
                            qhi: true,
                            x: self.x,
                        };
                    }

                    if self.fbig < 0.0 {
                        // status = -1
                        return InvrAction::Failed {
                            qleft: false,
                            qhi: false,
                            x: self.x,
                        };
                    }
                } else if self.fbig < self.fsmall {
                    if self.fsmall < 0.0 {
                        // status = -1
                        return InvrAction::Failed {
                            qleft: true,
                            qhi: false,
                            x: self.x,
                        };
                    }

                    if 0.0 < self.fbig {
                        // status = -1
                        return InvrAction::Failed {
                            qleft: false,
                            qhi: true,
                            x: self.x,
                        };
                    }
                }

                self.x = self.xsave;
                self.step = self.absstp.max(self.relstp * self.x.abs());
                // YY = F(X) - Y
                // GET-function-VALUE
                self.i99999 = Resume::Label90;
                InvrAction::NeedEval(self.x)
            }
            Resume::Label90 => {
                // Label 90.
                let yy = fx;

                if yy == 0.0 {
                    // status = 0
                    return InvrAction::Converged(self.x);
                }

                let qup = (self.qincr && (yy < 0.0)) || (!self.qincr && (0.0 < yy));
                // Handle case in which we must step higher.
                if !qup {
                    return self.label_170();
                }

                self.xlb = self.xsave;
                self.xub = (self.xlb + self.step).min(self.big);
                self.label_120()
            }
            Resume::Label130 => {
                // Label 130.
                let yy = fx;
                let qbdd = (self.qincr && (0.0 <= yy)) || (!self.qincr && (yy <= 0.0));
                let qlim = self.big <= self.xub;
                let qcond = qbdd || qlim;

                if !qcond {
                    self.step = self.stpmul * self.step;
                    self.xlb = self.xub;
                    self.xub = (self.xlb + self.step).min(self.big);
                }

                // Label 110.
                if qcond {
                    return self.label_150(qbdd, qlim);
                }

                self.label_120()
            }
            Resume::Label200 => {
                // Label 200.
                let yy = fx;
                let qbdd = (self.qincr && (yy <= 0.0)) || (!self.qincr && (0.0 <= yy));
                let qlim = self.xlb <= self.small;
                let qcond = qbdd || qlim;

                if !qcond {
                    self.step = self.stpmul * self.step;
                    self.xub = self.xlb;
                    self.xlb = (self.xub - self.step).max(self.small);
                }

                // Label 180.
                if qcond {
                    return self.label_220(qbdd, qlim);
                }

                self.label_190()
            }
            Resume::Label270 => {
                // Label 270: go to 250, where status = 1 falls through to
                // label 260.
                self.label_260(fx)
            }
        }
    }

    // Label 120 (cdflib.f90:8368-8378): YY = F(XUB) - Y.
    #[inline]
    fn label_120(&mut self) -> InvrAction {
        self.x = self.xub;
        // GET-function-VALUE
        self.i99999 = Resume::Label130;
        InvrAction::NeedEval(self.x)
    }

    // Label 150 (cdflib.f90:8396-8406).
    #[inline]
    fn label_150(&mut self, qbdd: bool, qlim: bool) -> InvrAction {
        if qlim && !qbdd {
            // status = -1
            self.x = self.big;
            return InvrAction::Failed {
                qleft: false,
                qhi: !self.qincr,
                x: self.x,
            };
        }

        self.label_240()
    }

    // Label 170 (cdflib.f90:8408-8414): handle the case in which we must
    // step lower.
    #[inline]
    fn label_170(&mut self) -> InvrAction {
        self.xub = self.xsave;
        self.xlb = (self.xub - self.step).max(self.small);
        self.label_190()
    }

    // Label 190 (cdflib.f90:8422-8432): YY = F(XLB) - Y.
    #[inline]
    fn label_190(&mut self) -> InvrAction {
        self.x = self.xlb;
        // GET-function-VALUE
        self.i99999 = Resume::Label200;
        InvrAction::NeedEval(self.x)
    }

    // Label 220 (cdflib.f90:8450-8458), falling through to label 240.
    #[inline]
    fn label_220(&mut self, qbdd: bool, qlim: bool) -> InvrAction {
        if qlim && !qbdd {
            // status = -1
            self.x = self.small;
            return InvrAction::Failed {
                qleft: true,
                qhi: self.qincr,
                x: self.x,
            };
        }

        self.label_240()
    }

    // Label 240 (cdflib.f90:8460-8467).
    #[inline]
    fn label_240(&mut self) -> InvrAction {
        // call dstzr ( xlb, xub, abstol, reltol )
        self.zror = ZrorState::new(ZrorConfig {
            xlo: self.xlb,
            xhi: self.xub,
            abstol: self.abstol,
            reltol: self.reltol,
        });
        // If we reach here, XLB and XUB bound the zero of F.
        // status = 0, so dzror starts a new search and ignores fx.
        self.label_260(0.0)
    }

    // Labels 250 and 260 (cdflib.f90:8469-8489).
    #[inline]
    fn label_260(&mut self, fx: f64) -> InvrAction {
        // call dzror ( status, x, fx, xlo, xhi, qdum1, qdum2 )
        match self.zror.step(fx) {
            ZrorAction::NeedEval(x) => {
                self.x = x;
                // GET-function-VALUE
                self.i99999 = Resume::Label270;
                InvrAction::NeedEval(self.x)
            }
            // status /= 1: go to 250. A dzror failure (status = -1, with
            // its qleft and qhi in qdum1 and qdum2) ends here too.
            ZrorAction::Converged { xlo, .. } | ZrorAction::Failed { xlo, .. } => {
                // Label 250.
                self.x = xlo;
                // status = 0
                InvrAction::Converged(self.x)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> InvrConfig {
        InvrConfig {
            small: 0.0,
            big: 1.0,
            absstp: 0.5,
            relstp: 0.5,
            stpmul: 5.0,
            abstol: 1.0e-50,
            reltol: 1.0e-8,
        }
    }

    #[test]
    fn rejects_increasing_range_when_small_is_already_positive() {
        let mut state = InvrState::new(cfg(), 0.5);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(1.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(
            state.step(2.0),
            InvrAction::Failed {
                qleft: true,
                qhi: true,
                x: 1.0
            }
        ));
    }

    #[test]
    fn rejects_increasing_range_when_big_is_still_negative() {
        let mut state = InvrState::new(cfg(), 0.5);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(-2.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(
            state.step(-1.0),
            InvrAction::Failed {
                qleft: false,
                qhi: false,
                x: 1.0
            }
        ));
    }

    #[test]
    fn rejects_decreasing_range_when_small_is_already_negative() {
        let mut state = InvrState::new(cfg(), 0.5);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(-1.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(
            state.step(-2.0),
            InvrAction::Failed {
                qleft: true,
                qhi: false,
                x: 1.0
            }
        ));
    }

    #[test]
    fn rejects_decreasing_range_when_big_is_still_positive() {
        let mut state = InvrState::new(cfg(), 0.5);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(2.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(
            state.step(1.0),
            InvrAction::Failed {
                qleft: false,
                qhi: true,
                x: 1.0
            }
        ));
    }

    #[test]
    fn reports_upper_bound_failure_when_search_runs_out_of_room() {
        let mut state = InvrState::new(cfg(), 0.9);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(-1.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(state.step(1.0), InvrAction::NeedEval(0.9)));
        assert!(matches!(state.step(-0.1), InvrAction::NeedEval(1.0)));
        assert!(matches!(
            state.step(-0.05),
            InvrAction::Failed {
                qleft: false,
                qhi: false,
                x: 1.0
            }
        ));
    }

    #[test]
    fn reports_lower_bound_failure_when_search_runs_out_of_room() {
        let mut state = InvrState::new(cfg(), 0.1);
        assert!(matches!(state.step(0.0), InvrAction::NeedEval(0.0)));
        assert!(matches!(state.step(-1.0), InvrAction::NeedEval(1.0)));
        assert!(matches!(state.step(1.0), InvrAction::NeedEval(0.1)));
        assert!(matches!(state.step(0.1), InvrAction::NeedEval(0.0)));
        assert!(matches!(
            state.step(0.05),
            InvrAction::Failed {
                qleft: true,
                qhi: true,
                x: 0.0
            }
        ));
    }
}
