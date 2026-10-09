//! The F90 subroutine `dzror` and its entry `dstzr` (cdflib.f90:8819-9209).
//!
//! `dzror` seeks a zero of a function, using reverse communication, by
//! Algorithm R of Bus and Dekker: a hybrid of linear interpolation,
//! inverse quadratic interpolation and bisection. The F90 keeps its state
//! in `save` variables and resumes at the label stored in `i99999`; here
//! the saved variables are the fields of [`ZrorState`], with the F90
//! names, and the resume label is a private enum.
//!
//! Reference: J. C. P. Bus and T. J. Dekker, Two Efficient Algorithms with
//! Guaranteed Convergence for Finding a Zero of a Function, *ACM
//! Transactions on Mathematical Software*, 1(4):330-345, 1975.

/// The arguments of the F90 entry `dstzr` (cdflib.f90:9130).
///
/// `xlo` and `xhi` (`zxlo`, `zxhi`) are the left and right endpoints of
/// the interval to be searched for a solution; `abstol` and `reltol`
/// (`zabstl`, `zreltl`) are two numbers that determine the accuracy of
/// the solution.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ZrorConfig {
    pub xlo: f64,
    pub xhi: f64,
    pub abstol: f64,
    pub reltol: f64,
}

/// The F90 resume label `i99999` (cdflib.f90:9194-9207).
#[derive(Debug, Clone, Copy)]
enum Resume {
    /// No evaluation pending: the F90 call with `status` = 0.
    Entry,
    /// Resume at label 10 with F(`xlo`).
    Label10,
    /// Resume at label 20 with F(`xhi`).
    Label20,
    /// Resume at label 200 with F(`b`).
    Label200,
}

/// The F90 `status` on return from `dzror`.
#[derive(Debug, Clone, Copy)]
pub(crate) enum ZrorAction {
    /// `status` = 1: the function must be evaluated at *x*, and the value
    /// passed as `fx` to the next [`ZrorState::step`] call.
    NeedEval(f64),
    /// `status` = 0: `xlo` and `xhi` bound the answer. `x` is the F90
    /// reverse-communication variable, the last point handed out for
    /// evaluation; label 80 can move `xlo` without changing `x`. The F90
    /// direct callers of `dzror` (`cdfbet`, `cdfbin`, `cdfnbn`) report
    /// `x`, while `dinvr` reports `xlo` (cdflib.f90:8472).
    Converged {
        x: f64,
        xlo: f64,
        // xhi mirrors the F90 output; no caller reads it.
        #[allow(dead_code)]
        xhi: f64,
    },
    /// `status` = -1: F(`xlo`) and F(`xhi`) have the same sign. `qleft`
    /// is true if the search terminated unsuccessfully at `xlo`, false if
    /// it terminated unsuccessfully at `xhi`; `qhi` is true if
    /// *Y* < F(*X*) at the termination of the search and false if
    /// F(*X*) < *Y*. `xlo` is the F90 `xlo` on return, which `dinvr`
    /// reports as its approximate root (cdflib.f90:8472).
    Failed { xlo: f64, qleft: bool, qhi: bool },
}

/// The `save` variables of the F90 `dzror`, with the F90 names, together
/// with the dummy arguments `x`, `xlo` and `xhi`, which the F90 caller
/// passes back unchanged on every call.
#[derive(Debug)]
pub(crate) struct ZrorState {
    // Set by the entry dstzr.
    xxlo: f64,
    xxhi: f64,
    abstol: f64,
    reltol: f64,
    i99999: Resume,
    // The dummy arguments: x is the point handed out for evaluation, and
    // xlo and xhi bound the answer.
    x: f64,
    xlo: f64,
    xhi: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    fa: f64,
    fb: f64,
    fc: f64,
    fd: f64,
    mb: f64,
    w: f64,
    ext: i32,
    first: bool,
}

/// The F90 intrinsic `sign ( a, b )`: the magnitude of *a* with the sign
/// of *b*.
#[inline]
fn sign(a: f64, b: f64) -> f64 {
    a.copysign(b)
}

impl ZrorState {
    /// The F90 entry `dstzr` (cdflib.f90:9183-9187). The first `dzror`
    /// call, with `status` = 0, is the first [`step`](Self::step).
    ///
    /// Given a function F, find `xlo` such that F(`xlo`) = 0. Input
    /// condition: F is a function of a single argument and `xlo` and `xhi`
    /// are such that F(`xlo`) · F(`xhi`) ≤ 0. If the input condition is
    /// met, `dzror` returns `status` = 0 and the output values of `xlo` and
    /// `xhi` satisfy F(`xlo`) · F(`xhi`) ≤ 0, |F(`xlo`)| ≤ |F(`xhi`)| and
    /// |`xlo` − `xhi`| ≤ TOL(*X*), where TOL(*X*) = max(`abstol`,
    /// `reltol` · |*X*|).
    #[inline]
    pub(crate) fn new(cfg: ZrorConfig) -> Self {
        Self {
            xxlo: cfg.xlo,
            xxhi: cfg.xhi,
            abstol: cfg.abstol,
            reltol: cfg.reltol,
            i99999: Resume::Entry,
            x: 0.0,
            xlo: 0.0,
            xhi: 0.0,
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            fa: 0.0,
            fb: 0.0,
            fc: 0.0,
            fd: 0.0,
            mb: 0.0,
            w: 0.0,
            ext: 0,
            first: false,
        }
    }

    /// The F90 statement function `ftol` (cdflib.f90:8929).
    #[inline]
    fn ftol(&self, zx: f64) -> f64 {
        0.5 * self.abstol.max(self.reltol * zx.abs())
    }

    /// One call of the F90 `dzror` (cdflib.f90:8819).
    ///
    /// The first call is the F90 call with `status` = 0, and `fx` is
    /// ignored. When `dzror` needs the function evaluated, it returns
    /// [`ZrorAction::NeedEval`] (`status` = 1); the value of the function
    /// must be passed as `fx` to the next call. When `dzror` has finished
    /// without error, it returns [`ZrorAction::Converged`] (`status` = 0);
    /// if it finds an error (which implies that F(`xlo`) − *Y* and
    /// F(`xhi`) − *Y* have the same sign), it returns
    /// [`ZrorAction::Failed`] (`status` = -1).
    #[allow(clippy::collapsible_if, clippy::assign_op_pattern)]
    #[inline]
    pub(crate) fn step(&mut self, fx: f64) -> ZrorAction {
        // cdflib.f90:8931-8933: if 0 < status, go to label 280, which
        // resumes at label i99999 (cdflib.f90:9194-9207). Rust only: the
        // F90 stop 1 for an illegal i99999 cannot occur, as Resume has no
        // other values.
        match self.i99999 {
            Resume::Entry => {
                self.xlo = self.xxlo;
                self.xhi = self.xxhi;
                self.b = self.xlo;
                self.x = self.xlo;
                // GET-function-VALUE
                self.i99999 = Resume::Label10;
                self.label_270()
            }
            Resume::Label10 => {
                // Label 10.
                self.fb = fx;
                self.xlo = self.xhi;
                self.a = self.xlo;
                self.x = self.xlo;
                // GET-function-VALUE
                self.i99999 = Resume::Label20;
                self.label_270()
            }
            Resume::Label20 => {
                // Label 20. Check that F(ZXLO) < 0 < F(ZXHI) or
                // F(ZXLO) > 0 > F(ZXHI).
                if self.fb < 0.0 {
                    if fx < 0.0 {
                        // status = -1
                        return ZrorAction::Failed {
                            xlo: self.xlo,
                            qleft: fx < self.fb,
                            qhi: false,
                        };
                    }
                }

                if 0.0 < self.fb {
                    if 0.0 < fx {
                        // status = -1
                        return ZrorAction::Failed {
                            xlo: self.xlo,
                            qleft: self.fb < fx,
                            qhi: true,
                        };
                    }
                }

                self.fa = fx;
                self.first = true;
                self.label_70()
            }
            Resume::Label200 => {
                // Label 200.
                self.fb = fx;

                if 0.0 <= self.fc * self.fb {
                    self.label_70()
                } else {
                    if self.w == self.mb {
                        self.ext = 0;
                    } else {
                        self.ext = self.ext + 1;
                    }

                    self.label_80()
                }
            }
        }
    }

    // Label 70 (cdflib.f90:8982-8986), falling through to label 80.
    #[inline]
    fn label_70(&mut self) -> ZrorAction {
        self.c = self.a;
        self.fc = self.fa;
        self.ext = 0;
        self.label_80()
    }

    // Label 80 (cdflib.f90:8988-9094): the interpolation step, through
    // labels 150, 180 and 190, up to the request for F(b).
    #[allow(
        clippy::assign_op_pattern,
        clippy::needless_late_init,
        clippy::neg_cmp_op_on_partial_ord
    )]
    #[inline]
    fn label_80(&mut self) -> ZrorAction {
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

        let mut tol = self.ftol(self.xlo);
        let m = (self.c + self.b) * 0.5;
        self.mb = m - self.b;

        if !(tol < self.mb.abs()) {
            return self.label_240();
        }

        'l190: {
            'l180: {
                if 3 < self.ext {
                    self.w = self.mb;
                    break 'l190;
                }

                tol = sign(tol, self.mb);
                let mut p = (self.b - self.a) * self.fb;
                let mut q;
                // The F90 guards both divided differences against a zero
                // denominator (cdflib.f90:9025-9026).
                if self.first {
                    q = self.fa - self.fb;
                    self.first = false;
                } else {
                    let fdb;
                    let fda;

                    if self.d == self.b {
                        fdb = 1.0;
                    } else {
                        fdb = (self.fd - self.fb) / (self.d - self.b);
                    }

                    if self.d == self.a {
                        fda = 1.0;
                    } else {
                        fda = (self.fd - self.fa) / (self.d - self.a);
                    }

                    p = fda * p;
                    q = fdb * self.fa - fda * self.fb;
                }

                if p < 0.0 {
                    p = -p;
                    q = -q;
                }

                if self.ext == 3 {
                    p = p * 2.0;
                }

                'l150: {
                    if !((p * 1.0) == 0.0 || p <= (q * tol)) {
                        break 'l150;
                    }

                    self.w = tol;
                    break 'l180;
                }

                // Label 150.
                if p < self.mb * q {
                    self.w = p / q;
                } else {
                    self.w = self.mb;
                }
            }
        }

        // Labels 180 and 190.
        self.d = self.a;
        self.fd = self.fa;
        self.a = self.b;
        self.fa = self.fb;
        self.b = self.b + self.w;
        self.xlo = self.b;
        self.x = self.xlo;
        // GET-function-VALUE
        self.i99999 = Resume::Label200;
        self.label_270()
    }

    // Label 240 (cdflib.f90:9116-9128).
    #[inline]
    fn label_240(&mut self) -> ZrorAction {
        self.xhi = self.c;
        let qrzero = (0.0 <= self.fc && self.fb <= 0.0) || (self.fc < 0.0 && self.fb >= 0.0);

        if qrzero {
            // status = 0
            ZrorAction::Converged {
                x: self.x,
                xlo: self.xlo,
                xhi: self.xhi,
            }
        } else {
            // status = -1. Rust only: the F90 leaves qleft and qhi unset
            // here, and the Rust reports both as false.
            ZrorAction::Failed {
                xlo: self.xlo,
                qleft: false,
                qhi: false,
            }
        }
    }

    // Label 270 (cdflib.f90:9191-9192): TO GET-function-VALUE, that is,
    // status = 1 and return.
    #[inline]
    fn label_270(&self) -> ZrorAction {
        ZrorAction::NeedEval(self.x)
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
            i99999: Resume::Label200,
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
            ..ZrorState::new(cfg())
        };

        let action = z.label_80();
        match action {
            ZrorAction::NeedEval(x) => assert!((x - 4.0 / 3.0).abs() < 1e-15),
            other => panic!("unexpected action: {other:?}"),
        }
    }

    #[test]
    fn guarded_divided_difference_when_d_equals_a() {
        let mut z = ZrorState {
            i99999: Resume::Label200,
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
            ..ZrorState::new(cfg())
        };

        let action = z.label_80();
        match action {
            ZrorAction::NeedEval(x) => assert!((x - (2.0 + 1.0 / 23.0)).abs() < 1e-15),
            other => panic!("unexpected action: {other:?}"),
        }
    }

    #[test]
    fn guarded_divided_difference_when_d_equals_b() {
        let mut z = ZrorState {
            i99999: Resume::Label200,
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
            ..ZrorState::new(cfg())
        };

        let action = z.label_80();
        match action {
            ZrorAction::NeedEval(x) => assert_eq!(x, 3.0),
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
            i99999: Resume::Label200,
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
            ..ZrorState::new(cfg())
        };

        assert!(matches!(
            z.label_80(),
            ZrorAction::Failed {
                qleft: false,
                qhi: false,
                ..
            }
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
