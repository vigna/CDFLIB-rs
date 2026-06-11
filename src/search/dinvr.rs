//! CDFLIB's `dinvr` / `E0000` range-and-invert search.
//!
//! Like [`super::dzror::ZrorState`] this is a state machine driven by
//! reverse communication. The expansion phase doubles the step until
//! the search interval contains a sign change, then hands off to
//! [`super::dzror::ZrorState`] for refinement.

use super::dzror::{ZrorAction, ZrorConfig, ZrorState};

/// Configuration mirroring CDFLIB's `dstinv`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct InvrConfig {
    pub small: f64,
    pub big: f64,
    pub abs_step: f64,
    pub rel_step: f64,
    pub stp_mul: f64,
    pub abs_tol: f64,
    pub rel_tol: f64,
}

#[derive(Debug, Clone, Copy)]
enum Stage {
    /// About to evaluate F(small).
    Start,
    AwaitFsmall,
    AwaitFbig,
    AwaitInitial,
    AwaitUpper,
    AwaitLower,
    InZror,
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum InvrAction {
    NeedEval(f64),
    Converged(f64),
    /// CDFLIB-style failure flags. `qleft` is true if the stepping
    /// search terminated unsuccessfully at the lower bound `small`.
    /// `qhi` and `x` mirror CDFLIB's interface; only `qleft` is consulted
    /// by the current search_monotone driver, but the others are retained
    /// for callers that drive `InvrState` directly.
    #[allow(dead_code)]
    Failed {
        qleft: bool,
        qhi: bool,
        x: f64,
    },
}

#[derive(Debug)]
pub(crate) struct InvrState {
    cfg: InvrConfig,
    stage: Stage,
    xsave: f64,
    fsmall: f64,
    fbig: f64,
    qincr: bool,
    step: f64,
    xlb: f64,
    xub: f64,
    zror: Option<ZrorState>,
}

impl InvrState {
    #[inline]
    pub(crate) fn new(cfg: InvrConfig, x_initial: f64) -> Self {
        Self {
            cfg,
            stage: Stage::Start,
            xsave: x_initial,
            fsmall: 0.0,
            fbig: 0.0,
            qincr: false,
            step: 0.0,
            xlb: 0.0,
            xub: 0.0,
            zror: None,
        }
    }

    /// Returns the next action of the root-finder after driving one
    /// iteration. On the first call, `fx` is ignored (no evaluation has
    /// happened yet). On subsequent calls, `fx` must be the value of *f*
    /// at the *x* from the previous `NeedEval`.
    #[inline]
    pub(crate) fn step(&mut self, fx: f64) -> InvrAction {
        match self.stage {
            Stage::Start => {
                // x = small; request F(small).
                self.stage = Stage::AwaitFsmall;
                InvrAction::NeedEval(self.cfg.small)
            }
            Stage::AwaitFsmall => {
                self.fsmall = fx;
                self.stage = Stage::AwaitFbig;
                InvrAction::NeedEval(self.cfg.big)
            }
            Stage::AwaitFbig => {
                self.fbig = fx;
                self.qincr = self.fbig > self.fsmall;
                // Check that the zero lies in [small, big]. F90 dinvr
                // (cdflib.f90:8054-8083) splits on fsmall <= fbig
                // (inclusive) versus fbig < fsmall (strict). The tie
                // case fsmall == fbig takes the <= branch. qincr above
                // is strict, used only for search-direction logic below.
                // F90 does not assign x in these four returns, so x
                // keeps the last requested evaluation point, big.
                if self.fsmall <= self.fbig {
                    if self.fsmall > 0.0 {
                        return InvrAction::Failed {
                            qleft: true,
                            qhi: true,
                            x: self.cfg.big,
                        };
                    }
                    if self.fbig < 0.0 {
                        return InvrAction::Failed {
                            qleft: false,
                            qhi: false,
                            x: self.cfg.big,
                        };
                    }
                } else {
                    if self.fsmall < 0.0 {
                        return InvrAction::Failed {
                            qleft: true,
                            qhi: false,
                            x: self.cfg.big,
                        };
                    }
                    if self.fbig > 0.0 {
                        return InvrAction::Failed {
                            qleft: false,
                            qhi: true,
                            x: self.cfg.big,
                        };
                    }
                }
                // S70/S80: start search from xsave.
                let x = self.xsave;
                self.step = self.cfg.abs_step.max(self.cfg.rel_step * x.abs());
                self.stage = Stage::AwaitInitial;
                InvrAction::NeedEval(x)
            }
            Stage::AwaitInitial => {
                let yy = fx;
                if yy == 0.0 {
                    return InvrAction::Converged(self.xsave);
                }
                let qup = (self.qincr && yy < 0.0) || (!self.qincr && yy > 0.0);
                if qup {
                    // Step higher.
                    self.xlb = self.xsave;
                    self.xub = (self.xlb + self.step).min(self.cfg.big);
                    self.stage = Stage::AwaitUpper;
                    InvrAction::NeedEval(self.xub)
                } else {
                    // Step lower.
                    self.xub = self.xsave;
                    self.xlb = (self.xub - self.step).max(self.cfg.small);
                    self.stage = Stage::AwaitLower;
                    InvrAction::NeedEval(self.xlb)
                }
            }
            Stage::AwaitUpper => {
                let yy = fx;
                let qbdd = (self.qincr && yy >= 0.0) || (!self.qincr && yy <= 0.0);
                let qlim = self.xub >= self.cfg.big;
                if qbdd {
                    // Sign change found: start dzror on (xlb, xub).
                    return self.start_zror();
                }
                if qlim {
                    return InvrAction::Failed {
                        qleft: false,
                        qhi: !self.qincr,
                        x: self.cfg.big,
                    };
                }
                // Double step and continue upward.
                self.step *= self.cfg.stp_mul;
                self.xlb = self.xub;
                self.xub = (self.xlb + self.step).min(self.cfg.big);
                InvrAction::NeedEval(self.xub)
            }
            Stage::AwaitLower => {
                let yy = fx;
                let qbdd = (self.qincr && yy <= 0.0) || (!self.qincr && yy >= 0.0);
                let qlim = self.xlb <= self.cfg.small;
                if qbdd {
                    return self.start_zror();
                }
                if qlim {
                    return InvrAction::Failed {
                        qleft: true,
                        qhi: self.qincr,
                        x: self.cfg.small,
                    };
                }
                self.step *= self.cfg.stp_mul;
                self.xub = self.xlb;
                self.xlb = (self.xub - self.step).max(self.cfg.small);
                InvrAction::NeedEval(self.xlb)
            }
            Stage::InZror => {
                let z = self.zror.as_mut().expect("zror");
                match z.step(fx) {
                    ZrorAction::NeedEval(x) => InvrAction::NeedEval(x),
                    ZrorAction::Converged { xlo, .. } => InvrAction::Converged(xlo),
                    // F90 dinvr (cdflib.f90:8233-8237 label 250) treats a
                    // dzror failure (status = -1) identically to convergence:
                    // it executes "x = xlo; status = 0; return". Mirror that
                    // silent-success behavior here.
                    ZrorAction::Failed { xlo, .. } => InvrAction::Converged(xlo),
                }
            }
        }
    }

    #[inline]
    fn start_zror(&mut self) -> InvrAction {
        let mut z = ZrorState::new(ZrorConfig {
            xlo: self.xlb,
            xhi: self.xub,
            abstol: self.cfg.abs_tol,
            reltol: self.cfg.rel_tol,
        });
        // Drive the first step (which always requests an eval).
        let first = z.step(0.0);
        self.zror = Some(z);
        self.stage = Stage::InZror;
        match first {
            ZrorAction::NeedEval(x) => InvrAction::NeedEval(x),
            ZrorAction::Converged { xlo, .. } => InvrAction::Converged(xlo),
            // Same F90 silent-success-at-xlo behavior as the InZror handler
            // above (cdflib.f90:8233-8237).
            ZrorAction::Failed { xlo, .. } => InvrAction::Converged(xlo),
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
            abs_step: 0.5,
            rel_step: 0.5,
            stp_mul: 5.0,
            abs_tol: 1.0e-50,
            rel_tol: 1.0e-8,
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
