//! CDFLIB's `dzror` / `E0001` zero-finder.
//!
//! Algorithm R of Bus & Dekker (ACM TOMS 1975): a hybrid of linear
//! and inverse quadratic interpolation. CDFLIB uses a
//! reverse-communication idiom driven by a `static`-local switch on a
//! “where to resume” integer `i99999`; this module uses an explicit
//! `Stage` enum carried inside [`ZrorState`] instead.
//!
//! The exact iteration trace matches CDFLIB. Variable names follow the
//! CDFLIB source (`a`, `b`, `c`, `d`, `fa`, `fb`, `fc`, `fd`, `fda`,
//! `fdb`, `m`, `mb`, `p`, `q`, `w`, `tol`, `ext`, `first`) for
//! line-by-line cross-referencing.

/// Configuration mirroring CDFLIB's `dstzr`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ZrorConfig {
    pub xlo: f64,
    pub xhi: f64,
    pub abstol: f64,
    pub reltol: f64,
}

#[derive(Debug, Clone, Copy)]
enum Stage {
    /// Initial entry: about to request F(xlo).
    Start,
    /// Awaiting F(xlo); next is to request F(xhi).
    AwaitFb,
    /// Awaiting F(xhi); next is to validate the sign change.
    AwaitFa,
    /// Awaiting F(b) after an interpolation step.
    AwaitFbStep,
}

/// One step of the [`ZrorState`] machine.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ZrorAction {
    /// Caller must evaluate *f*(*x*) and pass the result to the next `step` call.
    NeedEval(f64),
    /// Successful convergence: a root lies in [`xlo`, `xhi`]. `x` is the
    /// F90 reverse-communication variable: the point handed out for the
    /// last evaluation. The swap at F90 label 80 updates `xlo` but not
    /// `x`, so the two can differ by up to the convergence tolerance.
    /// The direct-`dzror` dispatchers (`cdfbet`, `cdfbin`, `cdfnbn`)
    /// report `x`, while the `dinvr`-embedded path overwrites it with
    /// `xlo` (cdflib.f90:8233).
    Converged {
        x: f64,
        xlo: f64,
        #[allow(dead_code)]
        xhi: f64,
    },
    /// *f*(`xlo`) and *f*(`xhi`) do not straddle zero. `qleft` / `qhi`
    /// follow the CDFLIB sign-flag convention. `xlo` is the last lower
    /// bound dzror maintained — F90 cdflib.f90:8233 returns it as the
    /// approximate root on failure.
    #[allow(dead_code)]
    Failed { xlo: f64, qleft: bool, qhi: bool },
}

#[derive(Debug)]
pub(crate) struct ZrorState {
    cfg: ZrorConfig,
    stage: Stage,
    // working state, names from the CDFLIB source
    xlo: f64,
    xhi: f64,
    // the reverse-communication variable: last point handed out for
    // evaluation
    x: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    fa: f64,
    fb: f64,
    fc: f64,
    fd: f64,
    w: f64,
    mb: f64,
    ext: i32,
    first: bool,
}

impl ZrorState {
    #[inline]
    pub(crate) fn new(cfg: ZrorConfig) -> Self {
        Self {
            cfg,
            stage: Stage::Start,
            xlo: 0.0,
            xhi: 0.0,
            x: 0.0,
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            fa: 0.0,
            fb: 0.0,
            fc: 0.0,
            fd: 0.0,
            w: 0.0,
            mb: 0.0,
            ext: 0,
            first: true,
        }
    }

    /// Returns the next action of the root-finder after driving one
    /// iteration. On the first call, `fx` is ignored (no evaluation has
    /// happened yet). On subsequent calls, `fx` must be the value of *f*
    /// at the *x* from the previous `NeedEval`.
    #[inline]
    pub(crate) fn step(&mut self, fx: f64) -> ZrorAction {
        loop {
            match self.stage {
                Stage::Start => {
                    self.xlo = self.cfg.xlo;
                    self.xhi = self.cfg.xhi;
                    self.b = self.xlo;
                    self.x = self.xlo;
                    self.stage = Stage::AwaitFb;
                    return ZrorAction::NeedEval(self.b);
                }
                Stage::AwaitFb => {
                    self.fb = fx;
                    self.xlo = self.xhi;
                    self.a = self.xlo;
                    self.x = self.xlo;
                    self.stage = Stage::AwaitFa;
                    return ZrorAction::NeedEval(self.a);
                }
                Stage::AwaitFa => {
                    // Validate sign change.
                    if self.fb < 0.0 {
                        if fx < 0.0 {
                            return ZrorAction::Failed {
                                xlo: self.xlo,
                                qleft: fx < self.fb,
                                qhi: false,
                            };
                        }
                    } else if self.fb > 0.0 && fx > 0.0 {
                        return ZrorAction::Failed {
                            xlo: self.xlo,
                            qleft: fx > self.fb,
                            qhi: true,
                        };
                    }
                    self.fa = fx;
                    self.first = true;
                    // Enter S70: c = a; fc = fa; ext = 0.
                    self.restart_c_from_a();
                    // Fall through to the swap+interpolation step.
                    if let Some(action) = self.refine_iteration() {
                        return action;
                    }
                }
                Stage::AwaitFbStep => {
                    self.fb = fx;
                    if self.fc * self.fb >= 0.0 {
                        // Sign-change lost: restart with c <- a.
                        self.restart_c_from_a();
                    } else if self.w == self.mb {
                        self.ext = 0;
                    } else {
                        self.ext += 1;
                    }
                    if let Some(action) = self.refine_iteration() {
                        return action;
                    }
                }
            }
        }
    }

    #[inline]
    fn restart_c_from_a(&mut self) {
        self.c = self.a;
        self.fc = self.fa;
        self.ext = 0;
    }

    /// Drive one round of the swap → check-convergence → step loop.
    /// Returns `Some(action)` if we need to leave (with eval request or
    /// terminal result); `None` if we've internally looped to the next
    /// iteration without external evaluation.
    #[inline]
    fn refine_iteration(&mut self) -> Option<ZrorAction> {
        // S80: swap so |fb| is the smaller residual.
        if self.fc.abs() < self.fb.abs() {
            if self.c == self.a {
                self.d = self.a;
                self.fd = self.fa;
            }
            self.a = self.b;
            self.fa = self.fb;
            self.xlo = self.c;
            self.b = self.xlo;
            self.fb = self.fc;
            self.c = self.a;
            self.fc = self.fa;
        }
        // S100: check convergence.
        let tol = 0.5 * self.cfg.abstol.max(self.cfg.reltol * self.xlo.abs());
        let m = 0.5 * (self.c + self.b);
        let mb = m - self.b;
        self.mb = mb;
        if mb.abs() <= tol {
            // Convergence section S240.
            self.xhi = self.c;
            let qrzero = (self.fc >= 0.0 && self.fb <= 0.0) || (self.fc < 0.0 && self.fb >= 0.0);
            if qrzero {
                return Some(ZrorAction::Converged {
                    x: self.x,
                    xlo: self.xlo,
                    xhi: self.xhi,
                });
            }
            return Some(ZrorAction::Failed {
                xlo: self.xlo,
                qleft: false,
                qhi: false,
            });
        }

        // S110 / step-size selection.
        let w;
        if self.ext > 3 {
            w = mb;
        } else {
            let tol_signed = tol.copysign(mb);
            let mut p = (self.b - self.a) * self.fb;
            let q;
            if self.first {
                q = self.fa - self.fb;
                self.first = false;
            } else {
                let fdb = if self.d == self.b {
                    1.0
                } else {
                    (self.fd - self.fb) / (self.d - self.b)
                };
                let fda = if self.d == self.a {
                    1.0
                } else {
                    (self.fd - self.fa) / (self.d - self.a)
                };
                p *= fda;
                q = fdb * self.fa - fda * self.fb;
            }
            let (mut p, q) = if p < 0.0 { (-p, -q) } else { (p, q) };
            if self.ext == 3 {
                p *= 2.0;
            }
            if p == 0.0 || p <= q * tol_signed {
                w = tol_signed;
            } else if p < mb * q {
                w = p / q;
            } else {
                w = mb;
            }
        }
        self.w = w;

        // S170: update history and step b.
        self.d = self.a;
        self.fd = self.fa;
        self.a = self.b;
        self.fa = self.fb;
        self.b += w;
        self.xlo = self.b;
        self.x = self.xlo;
        self.stage = Stage::AwaitFbStep;
        Some(ZrorAction::NeedEval(self.b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> ZrorConfig {
        ZrorConfig {
            xlo: 0.0,
            xhi: 1.0,
            abstol: 1.0e-50,
            reltol: 1.0e-8,
        }
    }

    #[test]
    fn swap_branch_preserves_history_when_c_equals_a() {
        let mut z = ZrorState {
            cfg: cfg(),
            stage: Stage::AwaitFbStep,
            xlo: 0.0,
            xhi: 0.0,
            x: 0.0,
            a: 1.0,
            b: 2.0,
            c: 1.0,
            d: 99.0,
            fa: 3.0,
            fb: 2.0,
            fc: 1.0,
            fd: 77.0,
            w: 0.0,
            mb: 0.0,
            ext: 0,
            first: false,
        };

        let action = z.refine_iteration();
        match action {
            Some(ZrorAction::NeedEval(x)) => assert!((x - 4.0 / 3.0).abs() < 1e-15),
            other => panic!("unexpected action: {other:?}"),
        }
    }

    #[test]
    fn guarded_divided_difference_when_d_equals_a() {
        let mut z = ZrorState {
            cfg: cfg(),
            stage: Stage::AwaitFbStep,
            xlo: 2.0,
            xhi: 4.0,
            x: 0.0,
            a: 1.0,
            b: 2.0,
            c: 4.0,
            d: 1.0,
            fa: 4.0,
            fb: -1.0,
            fc: 5.0,
            fd: 5.0,
            w: 0.0,
            mb: 0.0,
            ext: 0,
            first: false,
        };

        let action = z.refine_iteration();
        match action {
            Some(ZrorAction::NeedEval(x)) => assert!((x - (2.0 + 1.0 / 23.0)).abs() < 1e-15),
            other => panic!("unexpected action: {other:?}"),
        }
    }

    #[test]
    fn guarded_divided_difference_when_d_equals_b() {
        let mut z = ZrorState {
            cfg: cfg(),
            stage: Stage::AwaitFbStep,
            xlo: 2.0,
            xhi: 4.0,
            x: 0.0,
            a: 1.0,
            b: 2.0,
            c: 4.0,
            d: 2.0,
            fa: 4.0,
            fb: -1.0,
            fc: 5.0,
            fd: 5.0,
            w: 0.0,
            mb: 0.0,
            ext: 0,
            first: false,
        };

        let action = z.refine_iteration();
        match action {
            Some(ZrorAction::NeedEval(x)) => assert_eq!(x, 3.0),
            other => panic!("unexpected action: {other:?}"),
        }
    }

    #[test]
    fn fails_when_both_initial_values_are_negative() {
        let mut z = ZrorState::new(cfg());
        assert!(matches!(z.step(0.0), ZrorAction::NeedEval(0.0)));
        assert!(matches!(z.step(-2.0), ZrorAction::NeedEval(1.0)));
        assert!(matches!(
            z.step(-1.0),
            ZrorAction::Failed {
                qleft: false,
                qhi: false,
                ..
            }
        ));
    }

    #[test]
    fn fails_when_both_initial_values_are_positive() {
        let mut z = ZrorState::new(cfg());
        assert!(matches!(z.step(0.0), ZrorAction::NeedEval(0.0)));
        assert!(matches!(z.step(1.0), ZrorAction::NeedEval(1.0)));
        assert!(matches!(
            z.step(2.0),
            ZrorAction::Failed {
                qleft: true,
                qhi: true,
                ..
            }
        ));
    }

    #[test]
    fn reports_failed_convergence_if_interval_no_longer_straddles_zero() {
        let mut z = ZrorState {
            cfg: cfg(),
            stage: Stage::AwaitFbStep,
            xlo: 1.0,
            xhi: 1.0,
            x: 0.0,
            a: 1.0,
            b: 1.0,
            c: 1.0,
            d: 0.0,
            fa: 1.0,
            fb: 1.0,
            fc: 1.0,
            fd: 0.0,
            w: 0.0,
            mb: 0.0,
            ext: 0,
            first: false,
        };

        assert!(matches!(
            z.refine_iteration(),
            Some(ZrorAction::Failed {
                qleft: false,
                qhi: false,
                ..
            })
        ));
    }
}

#[cfg(test)]
mod protocol_tests {
    use super::*;

    // F90 dzror returns with the reverse-communication variable x still
    // holding the last point it handed out for evaluation; the swap at
    // label 80 updates xlo but not x. The direct-dzror dispatchers
    // (cdfbet, cdfbin, cdfnbn) report x, so Converged must carry the
    // last evaluation point separately from xlo. A step function makes
    // the |fc| < |fb| swap fire on the converging iteration, so the two
    // genuinely differ.
    #[test]
    fn converged_x_is_last_eval_point_not_xlo_after_final_swap() {
        let f = |x: f64| if x < 0.7 { -1.0 } else { 2.0 };
        let mut z = ZrorState::new(ZrorConfig {
            xlo: 0.0,
            xhi: 1.0,
            abstol: 1.0e-10,
            reltol: 1.0e-8,
        });
        let mut fx = 0.0;
        let mut last = f64::NAN;
        loop {
            match z.step(fx) {
                ZrorAction::NeedEval(x) => {
                    last = x;
                    fx = f(x);
                }
                ZrorAction::Converged { x, xlo, xhi } => {
                    assert_eq!(x, last, "x must be the last evaluation point");
                    assert_ne!(x, xlo, "the final swap must have moved xlo off x");
                    assert!((xlo - 0.7).abs() < 1e-7 && (xhi - 0.7).abs() < 1e-7);
                    return;
                }
                ZrorAction::Failed { .. } => panic!("no sign change reported"),
            }
        }
    }
}
