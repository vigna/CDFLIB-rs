#![cfg(not(miri))]

//! Replays the call logs of every cdf* dispatcher written by
//! `tests/regenerate/gen_dispatcher_calls.f90`.
//!
//! Each row records `which`, the arguments before the call, the F90
//! `status` and `bound`, and the arguments after the call. For every row
//! the Rust API can express, the test calls the method corresponding to
//! `which` and checks either the computed value, bit for bit, or that the
//! error maps back to the same `status` (and, for statuses 1 and 2, to the
//! same `bound`).
//!
//! Rows the API cannot express are skipped: a pair (*p*, *q*) that is not
//! exactly complementary passed to a method that derives one from the other,
//! a negative count passed to a `u64` argument, a number of successes above
//! the number of trials passed to the binomial `cdf` (which returns 1
//! instead of an error), the Rust-only rejections of a zero *r* or *pr* in
//! `NegativeBinomial`, and the NaN that `cdfnor` returns at *p* = 0, where
//! the Rust inverse returns −∞. Where the F90 returns a parameter of 0 with
//! status 0, at the lower end of a search interval, the Rust returns the
//! error of the constructor for that parameter instead.

mod common;

use cdflib::traits::{ContinuousCdf, DiscreteCdf};
use cdflib::{
    Beta, BetaError, Binomial, BinomialError, ChiSquared, ChiSquaredError, ChiSquaredNoncentral,
    ChiSquaredNoncentralError, FisherSnedecor, FisherSnedecorError, FisherSnedecorNoncentral,
    FisherSnedecorNoncentralError, Gamma, GammaError, NegativeBinomial, NegativeBinomialError,
    Normal, NormalError, Poisson, PoissonError, SearchError, StudentsT, StudentsTError,
};
use common::{assert_exact, read_csv};
use std::fmt::Debug;

/// The F90 `status` (and `bound`, for statuses 1 and 2) an error stands for.
trait F90Status {
    fn f90_status(&self) -> (i32, f64);
}

fn search_status(e: &SearchError) -> (i32, f64) {
    match *e {
        SearchError::AnswerBelowLowerBound { bound } => (1, bound),
        SearchError::AnswerAboveUpperBound { bound } => (2, bound),
        // F90 stops: never present in a call log.
        SearchError::StartOutOfRange { .. } => (-99, f64::NAN),
    }
}

impl F90Status for BetaError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            BetaError::PNotInRange(_) => (-2, f64::NAN),
            BetaError::QNotInRange(_) => (-3, f64::NAN),
            BetaError::XOutOfRange(_) => (-4, f64::NAN),
            BetaError::ANotPositive(_) => (-6, f64::NAN),
            BetaError::BNotPositive(_) => (-7, f64::NAN),
            BetaError::PQSumNotOne { .. } => (3, f64::NAN),
            BetaError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for BinomialError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            BinomialError::PNotInRange(_) => (-2, f64::NAN),
            BinomialError::QNotInRange(_) => (-3, f64::NAN),
            BinomialError::SuccessesExceedTrials { .. } => (-4, f64::NAN),
            BinomialError::TrialsZero => (-5, f64::NAN),
            BinomialError::PrOutOfRange(_) => (-6, f64::NAN),
            BinomialError::PQSumNotOne { .. } => (3, f64::NAN),
            BinomialError::Search(s) => search_status(s),
        }
    }
}

impl F90Status for ChiSquaredError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            ChiSquaredError::PNotInRange(_) => (-2, f64::NAN),
            ChiSquaredError::QNotInRange(_) => (-3, f64::NAN),
            ChiSquaredError::XNotPositive(_) => (-4, f64::NAN),
            ChiSquaredError::DfNotPositive(_) => (-5, f64::NAN),
            ChiSquaredError::PQSumNotOne { .. } => (3, f64::NAN),
            ChiSquaredError::Search(s) => search_status(s),
            ChiSquaredError::GammaInc(_) => (10, f64::NAN),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for ChiSquaredNoncentralError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            ChiSquaredNoncentralError::PNotInRange(_) => (-2, f64::NAN),
            ChiSquaredNoncentralError::XNotPositive(_) => (-4, f64::NAN),
            ChiSquaredNoncentralError::DfNotPositive(_) => (-5, f64::NAN),
            ChiSquaredNoncentralError::NcpNegative(_) => (-6, f64::NAN),
            ChiSquaredNoncentralError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for FisherSnedecorError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            FisherSnedecorError::PNotInRange(_) => (-2, f64::NAN),
            FisherSnedecorError::QNotInRange(_) => (-3, f64::NAN),
            FisherSnedecorError::FNotPositive(_) => (-4, f64::NAN),
            FisherSnedecorError::DfnNotPositive(_) => (-5, f64::NAN),
            FisherSnedecorError::DfdNotPositive(_) => (-6, f64::NAN),
            FisherSnedecorError::PQSumNotOne { .. } => (3, f64::NAN),
            FisherSnedecorError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for FisherSnedecorNoncentralError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            FisherSnedecorNoncentralError::PNotInRange(_) => (-2, f64::NAN),
            FisherSnedecorNoncentralError::FNotPositive(_) => (-4, f64::NAN),
            FisherSnedecorNoncentralError::DfnNotPositive(_) => (-5, f64::NAN),
            FisherSnedecorNoncentralError::DfdNotPositive(_) => (-6, f64::NAN),
            FisherSnedecorNoncentralError::NcpNegative(_) => (-7, f64::NAN),
            FisherSnedecorNoncentralError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for GammaError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            GammaError::PNotInRange(_) => (-2, f64::NAN),
            GammaError::QNotInRange(_) => (-3, f64::NAN),
            GammaError::XNotPositive(_) => (-4, f64::NAN),
            GammaError::ShapeNotPositive(_) => (-5, f64::NAN),
            GammaError::RateNotPositive(_) => (-6, f64::NAN),
            GammaError::PQSumNotOne { .. } => (3, f64::NAN),
            GammaError::Search(s) => search_status(s),
            GammaError::GammaIncInv(_) | GammaError::GammaInc(_) => (10, f64::NAN),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for NegativeBinomialError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            NegativeBinomialError::PNotInRange(_) => (-2, f64::NAN),
            NegativeBinomialError::QNotInRange(_) => (-3, f64::NAN),
            NegativeBinomialError::PrOutOfRange(_) => (-6, f64::NAN),
            NegativeBinomialError::PQSumNotOne { .. } => (3, f64::NAN),
            NegativeBinomialError::Search(s) => search_status(s),
            NegativeBinomialError::RNotPositive => (-99, f64::NAN),
        }
    }
}

impl F90Status for NormalError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            NormalError::PNotInRange(_) => (-2, f64::NAN),
            NormalError::QNotInRange(_) => (-3, f64::NAN),
            NormalError::SdNotPositive(_) => (-6, f64::NAN),
            NormalError::PQSumNotOne { .. } => (3, f64::NAN),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for PoissonError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            PoissonError::PNotInRange(_) => (-2, f64::NAN),
            PoissonError::QNotInRange(_) => (-3, f64::NAN),
            PoissonError::LambdaNegative(_) => (-5, f64::NAN),
            PoissonError::PQSumNotOne { .. } => (3, f64::NAN),
            PoissonError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

impl F90Status for StudentsTError {
    fn f90_status(&self) -> (i32, f64) {
        match self {
            StudentsTError::PNotInRange(_) => (-2, f64::NAN),
            StudentsTError::QNotInRange(_) => (-3, f64::NAN),
            StudentsTError::DfNotPositive(_) => (-5, f64::NAN),
            StudentsTError::PQSumNotOne { .. } => (3, f64::NAN),
            StudentsTError::Search(s) => search_status(s),
            _ => (-99, f64::NAN),
        }
    }
}

/// One row of a call log: which, the arguments before the call, status,
/// bound, and the arguments after the call.
struct Call {
    which: i32,
    before: Vec<f64>,
    status: i32,
    bound: f64,
    after: Vec<f64>,
    row: Vec<f64>,
}

fn calls(name: &str, nargs: usize) -> Vec<Call> {
    read_csv(&format!("tests/data/{name}_calls.csv"))
        .into_iter()
        .map(|row| {
            assert_eq!(row.len(), 2 * nargs + 3, "{name}: {row:?}");
            Call {
                which: row[0] as i32,
                before: row[1..=nargs].to_vec(),
                status: row[nargs + 1] as i32,
                bound: row[nargs + 2],
                after: row[nargs + 3..].to_vec(),
                row,
            }
        })
        .collect()
}

/// Checks a Rust result against the F90 status and the expected value.
#[track_caller]
fn check<E: F90Status + Debug>(got: Result<f64, E>, c: &Call, expected: f64) {
    match got {
        Ok(v) => {
            assert_eq!(0, c.status, "{:?}: got Ok({v:e})", c.row);
            assert_exact(v, expected, &c.row);
        }
        Err(e) => {
            let (status, bound) = e.f90_status();
            assert_eq!(status, c.status, "{:?}: got {e:?}", c.row);
            if status == 1 || status == 2 {
                assert_exact(bound, c.bound, &c.row);
            }
        }
    }
}

/// Checks a search for a parameter that must be positive: where the F90
/// returns 0 with status 0, at the lower end of the search interval, the
/// Rust reports the error that `zero` accepts instead.
#[track_caller]
fn check_positive<E: F90Status + Debug>(
    got: Result<f64, E>,
    c: &Call,
    expected: f64,
    zero: impl Fn(&E) -> bool,
) {
    if c.status == 0 && expected == 0.0 {
        match got {
            Err(e) if zero(&e) => {}
            r => panic!("{:?}: expected the error for 0, got {r:?}", c.row),
        }
    } else {
        check(got, c, expected);
    }
}

/// Checks a constructor error against the F90 status.
#[track_caller]
fn check_new<T: Debug, E: F90Status + Debug>(got: Result<T, E>, c: &Call) {
    match got {
        Ok(d) => panic!("{:?}: constructor accepted {d:?}", c.row),
        Err(e) => assert_eq!(e.f90_status().0, c.status, "{:?}: got {e:?}", c.row),
    }
}

/// Whether *q* is exactly 1 − *p*, as derived by a method taking *p* alone.
fn q_derived(p: f64, q: f64) -> bool {
    (1.0 - p).to_bits() == q.to_bits()
}

/// Whether *p* is exactly 1 − *q*, as derived by a method taking *q* alone.
fn p_derived(p: f64, q: f64) -> bool {
    (1.0 - q).to_bits() == p.to_bits()
}

/// Calls the inverse taking *p* or the one taking *q*, whichever receives
/// exactly the F90 pair; returns None when neither can express the row, or
/// when *p* or *q* is 0 or 1 and the F90 result is NaN, where the Rust
/// inverses return the exact endpoint of the support instead.
fn by_p_or_q<T>(c: &Call, by_p: impl FnOnce(f64) -> T, by_q: impl FnOnce(f64) -> T) -> Option<T> {
    let (p, q) = (c.before[0], c.before[1]);
    if c.status == 0 && (p == 0.0 || p == 1.0 || q == 0.0 || q == 1.0) && c.after[2].is_nan() {
        return None;
    }
    match c.status {
        -2 => Some(by_p(p)),
        -3 => Some(by_q(q)),
        _ if q_derived(p, q) => Some(by_p(p)),
        _ if p_derived(p, q) => Some(by_q(q)),
        _ => None,
    }
}

/// Whether *v* can be passed as a `u64` count.
fn count(v: f64) -> Option<u64> {
    (v >= 0.0 && v.fract() == 0.0).then_some(v as u64)
}

#[test]
fn cdfbet_calls() {
    let mut n = 0;
    for c in calls("cdfbet", 6) {
        let [p, q, x, y, a, b] = c.before[..] else {
            unreachable!()
        };
        // The Rust API derives y from x where x is an input.
        if c.which != 2 && (1.0 - x).to_bits() != y.to_bits() {
            continue;
        }
        let d = Beta::try_new(a, b);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(x), c.after[0], &c.row);
                    assert_exact(d.ccdf(x), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<Beta, _>(e), &c),
            },
            2 => match d {
                Ok(d) => {
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                }
                Err(e) => check_new(Err::<Beta, _>(e), &c),
            },
            3 => check_positive(Beta::search_a(p, q, x, b), &c, c.after[4], |e| {
                *e == BetaError::ANotPositive(0.0)
            }),
            4 => check_positive(Beta::search_b(p, q, x, a), &c, c.after[5], |e| {
                *e == BetaError::BNotPositive(0.0)
            }),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 37);
}

#[test]
fn cdfbin_calls() {
    let mut n = 0;
    for c in calls("cdfbin", 6) {
        let [p, q, s, xn, pr, ompr] = c.before[..] else {
            unreachable!()
        };
        if (1.0 - pr).to_bits() != ompr.to_bits() {
            continue;
        }
        match c.which {
            1 => {
                let (Some(s), Some(xn)) = (count(s), count(xn)) else {
                    continue;
                };
                match Binomial::try_new(xn, pr) {
                    Ok(d) if c.status == 0 => {
                        assert_exact(d.cdf(s), c.after[0], &c.row);
                        assert_exact(d.ccdf(s), c.after[1], &c.row);
                    }
                    // CDFLIB rejects s > xn; cdf returns 1 there.
                    Ok(_) => continue,
                    Err(e) => check_new(Err::<Binomial, _>(e), &c),
                }
            }
            2 => {
                let Some(xn) = count(xn) else { continue };
                match Binomial::try_new(xn, pr) {
                    Ok(d) => {
                        if c.status == -2 {
                            assert_eq!(d.inverse_cdf(p).unwrap_err().f90_status().0, -2);
                        } else if c.status == -3 || p_derived(p, q) {
                            check(d.inverse_ccdf(q), &c, c.after[2]);
                        } else {
                            continue;
                        }
                    }
                    Err(e) => check_new(Err::<Binomial, _>(e), &c),
                }
            }
            3 => {
                let Some(s) = count(s) else { continue };
                check_positive(Binomial::search_trials(p, q, pr, s), &c, c.after[3], |e| {
                    *e == BinomialError::TrialsZero
                });
            }
            4 => {
                let (Some(s), Some(xn)) = (count(s), count(xn)) else {
                    continue;
                };
                check(Binomial::search_pr(p, q, xn, s), &c, c.after[4]);
            }
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 38);
}

#[test]
fn cdfchi_calls() {
    let mut n = 0;
    for c in calls("cdfchi", 4) {
        let [p, q, x, df] = c.before[..] else {
            unreachable!()
        };
        let d = ChiSquared::try_new(df);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(x), c.after[0], &c.row);
                    assert_exact(d.ccdf(x), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<ChiSquared, _>(e), &c),
            },
            2 => match d {
                Ok(d) => {
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                }
                Err(e) => check_new(Err::<ChiSquared, _>(e), &c),
            },
            3 => check_positive(ChiSquared::search_df(p, q, x), &c, c.after[3], |e| {
                *e == ChiSquaredError::DfNotPositive(0.0)
            }),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 28);
}

#[test]
fn cdfchn_calls() {
    let mut n = 0;
    for c in calls("cdfchn", 5) {
        let [p, _q, x, df, pnonc] = c.before[..] else {
            unreachable!()
        };
        let d = ChiSquaredNoncentral::try_new(df, pnonc);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(x), c.after[0], &c.row);
                    assert_exact(d.ccdf(x), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<ChiSquaredNoncentral, _>(e), &c),
            },
            2 => match d {
                Ok(d) => check(d.inverse_cdf(p), &c, c.after[2]),
                Err(e) => check_new(Err::<ChiSquaredNoncentral, _>(e), &c),
            },
            3 => check_positive(
                ChiSquaredNoncentral::search_df(p, x, pnonc),
                &c,
                c.after[3],
                |e| *e == ChiSquaredNoncentralError::DfNotPositive(0.0),
            ),
            4 => check(ChiSquaredNoncentral::search_ncp(p, x, df), &c, c.after[4]),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 28);
}

#[test]
fn cdff_calls() {
    let mut n = 0;
    for c in calls("cdff", 5) {
        let [p, q, f, dfn, dfd] = c.before[..] else {
            unreachable!()
        };
        let d = FisherSnedecor::try_new(dfn, dfd);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(f), c.after[0], &c.row);
                    assert_exact(d.ccdf(f), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<FisherSnedecor, _>(e), &c),
            },
            2 => match d {
                Ok(d) => {
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                }
                Err(e) => check_new(Err::<FisherSnedecor, _>(e), &c),
            },
            3 => check(FisherSnedecor::search_dfn(p, q, f, dfd), &c, c.after[3]),
            4 => check(FisherSnedecor::search_dfd(p, q, f, dfn), &c, c.after[4]),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 32);
}

#[test]
fn cdffnc_calls() {
    let mut n = 0;
    for c in calls("cdffnc", 6) {
        let [p, _q, f, dfn, dfd, pnonc] = c.before[..] else {
            unreachable!()
        };
        let d = FisherSnedecorNoncentral::try_new(dfn, dfd, pnonc);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(f), c.after[0], &c.row);
                    assert_exact(d.ccdf(f), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<FisherSnedecorNoncentral, _>(e), &c),
            },
            2 => match d {
                Ok(d) => check(d.inverse_cdf(p), &c, c.after[2]),
                Err(e) => check_new(Err::<FisherSnedecorNoncentral, _>(e), &c),
            },
            3 => check(
                FisherSnedecorNoncentral::search_dfn(p, f, dfd, pnonc),
                &c,
                c.after[3],
            ),
            4 => check(
                FisherSnedecorNoncentral::search_dfd(p, f, dfn, pnonc),
                &c,
                c.after[4],
            ),
            5 => check(
                FisherSnedecorNoncentral::search_ncp(p, f, dfn, dfd),
                &c,
                c.after[5],
            ),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 32);
}

#[test]
fn cdfgam_calls() {
    let mut n = 0;
    for c in calls("cdfgam", 5) {
        let [p, q, x, shape, scale] = c.before[..] else {
            unreachable!()
        };
        let d = Gamma::try_new(shape, scale);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(x), c.after[0], &c.row);
                    assert_exact(d.ccdf(x), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<Gamma, _>(e), &c),
            },
            2 => match d {
                Ok(d) => {
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                }
                Err(e) => check_new(Err::<Gamma, _>(e), &c),
            },
            3 => check_positive(Gamma::search_shape(p, q, x, scale), &c, c.after[3], |e| {
                *e == GammaError::ShapeNotPositive(0.0)
            }),
            4 => check(Gamma::search_rate(p, q, x, shape), &c, c.after[4]),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 31);
}

#[test]
fn cdfnbn_calls() {
    let mut n = 0;
    for c in calls("cdfnbn", 6) {
        // CDFLIB's F is the number of failures, S the number of successes.
        let [p, q, f, s, pr, ompr] = c.before[..] else {
            unreachable!()
        };
        if (1.0 - pr).to_bits() != ompr.to_bits() {
            continue;
        }
        match c.which {
            1 => {
                let (Some(f), Some(s)) = (count(f), count(s)) else {
                    continue;
                };
                match NegativeBinomial::try_new(s, pr) {
                    Ok(d) => {
                        assert_eq!(0, c.status, "{:?}", c.row);
                        assert_exact(d.cdf(f), c.after[0], &c.row);
                        assert_exact(d.ccdf(f), c.after[1], &c.row);
                    }
                    Err(NegativeBinomialError::RNotPositive) => continue,
                    Err(e) => check_new(Err::<NegativeBinomial, _>(e), &c),
                }
            }
            2 => {
                let Some(s) = count(s) else { continue };
                match NegativeBinomial::try_new(s, pr) {
                    Ok(d) => {
                        if c.status == -2 {
                            assert_eq!(d.inverse_cdf(p).unwrap_err().f90_status().0, -2);
                        } else if c.status == -3 || p_derived(p, q) {
                            check(d.inverse_ccdf(q), &c, c.after[2]);
                        } else {
                            continue;
                        }
                    }
                    Err(NegativeBinomialError::RNotPositive) => continue,
                    Err(e) => check_new(Err::<NegativeBinomial, _>(e), &c),
                }
            }
            3 => {
                let Some(f) = count(f) else { continue };
                // Rust also rejects pr = 0.
                if pr == 0.0 {
                    continue;
                }
                check_positive(
                    NegativeBinomial::search_r(p, q, pr, f),
                    &c,
                    c.after[3],
                    |e| *e == NegativeBinomialError::RNotPositive,
                );
            }
            4 => {
                let (Some(f), Some(s)) = (count(f), count(s)) else {
                    continue;
                };
                // Rust also rejects s = 0 successes, which cdfnbn accepts.
                if s == 0 {
                    assert_eq!(
                        NegativeBinomial::search_pr(p, q, s, f),
                        Err(NegativeBinomialError::RNotPositive),
                        "{:?}",
                        c.row
                    );
                    continue;
                }
                check_positive(
                    NegativeBinomial::search_pr(p, q, s, f),
                    &c,
                    c.after[4],
                    |e| *e == NegativeBinomialError::PrOutOfRange(0.0),
                );
            }
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 30);
}

#[test]
fn cdfnor_calls() {
    let mut n = 0;
    for c in calls("cdfnor", 5) {
        let [p, q, x, mean, sd] = c.before[..] else {
            unreachable!()
        };
        let d = Normal::try_new(mean, sd);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(x), c.after[0], &c.row);
                    assert_exact(d.ccdf(x), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<Normal, _>(e), &c),
            },
            2 => {
                // cdfnor checks P, Q and P + Q before SD.
                if matches!(c.status, -2 | -3 | 3) || d.is_ok() {
                    let d = d.unwrap_or(Normal::standard());
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                } else {
                    check_new(d, &c);
                }
            }
            3 => check(Normal::search_mean(p, q, x, sd), &c, c.after[3]),
            4 => match Normal::search_sd(p, q, x, mean) {
                // Rust only: an sd that is not positive, which F90 returns
                // with status 0, is reported as SdNotPositive.
                Err(NormalError::SdNotPositive(sd)) if c.status == 0 => {
                    assert!(c.after[4] <= 0.0, "{:?}: got sd = {sd:e}", c.row);
                    assert_exact(sd, c.after[4], &c.row);
                }
                r => check(r, &c, c.after[4]),
            },
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 20);
}

#[test]
fn cdfpoi_calls() {
    let mut n = 0;
    for c in calls("cdfpoi", 4) {
        let [p, q, s, xlam] = c.before[..] else {
            unreachable!()
        };
        match c.which {
            1 => {
                let Some(s) = count(s) else { continue };
                match Poisson::try_new(xlam) {
                    Ok(d) => {
                        assert_eq!(0, c.status, "{:?}", c.row);
                        assert_exact(d.cdf(s), c.after[0], &c.row);
                        assert_exact(d.ccdf(s), c.after[1], &c.row);
                    }
                    Err(e) => check_new(Err::<Poisson, _>(e), &c),
                }
            }
            2 => match Poisson::try_new(xlam) {
                Ok(d) => {
                    if c.status == -2 {
                        assert_eq!(d.inverse_cdf(p).unwrap_err().f90_status().0, -2);
                    } else if c.status == -3 || p_derived(p, q) {
                        check(d.inverse_ccdf(q), &c, c.after[2]);
                    } else {
                        continue;
                    }
                }
                Err(e) => check_new(Err::<Poisson, _>(e), &c),
            },
            3 => {
                let Some(s) = count(s) else { continue };
                check(Poisson::search_lambda(p, q, s), &c, c.after[3]);
            }
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 21);
}

#[test]
fn cdft_calls() {
    let mut n = 0;
    for c in calls("cdft", 4) {
        let [p, q, t, df] = c.before[..] else {
            unreachable!()
        };
        let d = StudentsT::try_new(df);
        match c.which {
            1 => match d {
                Ok(d) => {
                    assert_eq!(0, c.status, "{:?}", c.row);
                    assert_exact(d.cdf(t), c.after[0], &c.row);
                    assert_exact(d.ccdf(t), c.after[1], &c.row);
                }
                Err(e) => check_new(Err::<StudentsT, _>(e), &c),
            },
            2 => {
                // cdft checks P and Q before DF.
                if matches!(c.status, -2 | -3) || d.is_ok() {
                    let d = d.unwrap_or(StudentsT::new(1.0));
                    let Some(r) = by_p_or_q(&c, |p| d.inverse_cdf(p), |q| d.inverse_ccdf(q)) else {
                        continue;
                    };
                    check(r, &c, c.after[2]);
                } else {
                    check_new(d, &c);
                }
            }
            3 => check(StudentsT::search_df(p, q, t), &c, c.after[3]),
            _ => unreachable!(),
        }
        n += 1;
    }
    assert_eq!(n, 23);
}

/// cdfbin with which = 2 starts its search at *s* = 5 in [0 . . *xn*]; for
/// xn < 5 F90 dinvr stops the program. Rust reports it as an error.
#[test]
fn cdfbin_start_outside_search_range_is_an_error() {
    let err = Binomial::new(3, 0.5).inverse_ccdf(0.5).unwrap_err();
    assert!(matches!(
        err,
        BinomialError::Search(SearchError::StartOutOfRange { .. })
    ));
}
