#![cfg(not(miri))]

//! Per-dispatcher tests pinning the literal `bound` value that each F90
//! `cdf*` dispatcher writes when its inner search fails at `status = 1`
//! (qleft, "answer below lower range") or `status = 2` (qhi, "answer above
//! upper range").
//!
//! The F90 call logs in tests/dispatcher_calls.rs check the same bounds
//! against cdflib.f90; these tests name the cases explicitly. Each
//! assertion cites the F90 source line where the `bound = …` literal is
//! written. Three dispatchers (cdft which=3, cdffnc which=3, cdffnc
//! which=4) write `bound = 0.0D+00` even though the lower end of their
//! search interval is 1.0.
//!
//! The last three tests pin the Rust status where the F90 status is
//! undefined: a `dzror` failure at label 240 leaves `qleft` unassigned,
//! so no call log can record it (see tests/regenerate/unreachable.txt).

use cdflib::{
    Beta, BetaError, Binomial, BinomialError, FisherSnedecorNoncentral,
    FisherSnedecorNoncentralError, NegativeBinomial, NegativeBinomialError, SearchError, StudentsT,
    StudentsTError,
};

// ---------------------------------------------------------------------------
// Drift sites: F90 writes bound = 0.0D+00 for qleft even though small > 0.
// ---------------------------------------------------------------------------

#[test]
fn students_t_search_df_qleft_bound_is_f90_zero_not_small() {
    // cdflib.f90:6475 writes bound = 0.0D+00 for cdft which=3 qleft,
    // despite small = 1.0 (cdflib.f90:6450 sets dstinv(1.0, maxdf, ...)).
    //
    // Trigger qleft with t = -2.0, p = q = 0.5: the t-CDF at t = -2 is
    // decreasing in df, with cum(-2, df=1) ≈ 0.148 and cum(-2, df→∞) → 0.023.
    // The target p = 0.5 lies above both endpoint values, so f = cum - p
    // is negative across the whole range, so qleft fires.
    let err = StudentsT::search_df(0.5, 0.5, -2.0).unwrap_err();
    assert!(
        matches!(
            err,
            StudentsTError::Search(SearchError::AnswerBelowLowerBound { bound }) if bound == 0.0
        ),
        "expected AnswerBelowLowerBound {{ bound: 0.0 }} per cdflib.f90:6475, got {err:?}"
    );
}

#[test]
fn students_t_search_df_qhi_bound_is_f90_maxdf() {
    // cdflib.f90:6482 writes bound = maxdf = 1.0D+10 for cdft which=3 qhi.
    //
    // Trigger qhi with t = -2.0, p = 0.001, q = 0.999: the search-residual
    // pivot uses cum - p (since p ≤ q is true). cum(-2, df) decreases
    // from ≈ 0.148 to ≈ 0.023; f = cum - 0.001 is positive everywhere,
    // and since f is decreasing the answer would lie above big, so qhi fires.
    let err = StudentsT::search_df(0.001, 0.999, -2.0).unwrap_err();
    assert!(
        matches!(
            err,
            StudentsTError::Search(SearchError::AnswerAboveUpperBound { bound }) if bound == 1.0e10
        ),
        "expected AnswerAboveUpperBound {{ bound: 1.0e10 }} per cdflib.f90:6482, got {err:?}"
    );
}

#[test]
fn fisher_snedecor_noncentral_search_dfn_qleft_bound_is_f90_zero_not_small() {
    // cdflib.f90:4758 writes bound = 0.0D+00 for cdffnc which=3 qleft,
    // despite small = 1.0 (cdflib.f90:4738).
    //
    // Trigger qleft with f = 2.0, dfd = 10, ncp = 0, p = 0.5: the central-F
    // CDF at f=2, (dfn=1, dfd=10) is ≈ 0.83. Increasing dfn pushes cum
    // even higher (the F distribution concentrates near 1, and 2 > 1
    // sits in the upper tail). So cum ≥ 0.83 everywhere in [1..1e30],
    // making f = cum - p positive, so qleft fires.
    let err = FisherSnedecorNoncentral::search_dfn(0.5, 2.0, 10.0, 0.0).unwrap_err();
    assert!(
        matches!(
            err,
            FisherSnedecorNoncentralError::Search(SearchError::AnswerBelowLowerBound { bound })
                if bound == 0.0
        ),
        "expected AnswerBelowLowerBound {{ bound: 0.0 }} per cdflib.f90:4758, got {err:?}"
    );
}

#[test]
fn fisher_snedecor_noncentral_search_dfd_qleft_bound_is_f90_zero_not_small() {
    // cdflib.f90:4796 writes bound = 0.0D+00 for cdffnc which=4 qleft,
    // despite small = 1.0 (cdflib.f90:4777).
    //
    // Trigger qleft with f = 0.5, dfn = 10, ncp = 0, p = 0.99: cum is
    // decreasing in dfd, from about 0.19 at dfd = 1 to about 0.11, the
    // χ²(10)/10 limit, as dfd grows, so cum - p is negative at both ends
    // of the search and decreasing, and dinvr reports qleft.
    let err = FisherSnedecorNoncentral::search_dfd(0.99, 0.5, 10.0, 0.0).unwrap_err();
    assert!(
        matches!(
            err,
            FisherSnedecorNoncentralError::Search(SearchError::AnswerBelowLowerBound { bound })
                if bound == 0.0
        ),
        "expected AnswerBelowLowerBound {{ bound: 0.0 }} per cdflib.f90:4796, got {err:?}"
    );
}

// ---------------------------------------------------------------------------
// dzror failures at label 240, where the F90 reads an uninitialized qleft.
// ---------------------------------------------------------------------------

// dzror assigns qleft and qhi only when it fails at label 20
// (cdflib.f90:8959). When it fails at label 240 (cdflib.f90:9125) the
// cdf* routine tests a qleft that nothing has assigned, so the F90 status
// is undefined; the Rust reports both flags as false, which selects the
// status 2 branch and its bound of 1.

#[test]
fn beta_inverse_ccdf_dzror_label_240_failure_is_status_2() {
    // cdfbet which=2, cdflib.f90:2778-2794: the search for y with
    // q = 1e-300 ends at label 240.
    let err = Beta::new(0.5, 1e5).inverse_ccdf(1e-300).unwrap_err();
    assert_eq!(
        err,
        BetaError::Search(SearchError::AnswerAboveUpperBound { bound: 1.0 })
    );
}

#[test]
fn binomial_search_pr_dzror_label_240_failure_is_status_2() {
    // cdfbin which=4, cdflib.f90:3318-3334: cdf(0) = (1 - pr)^100 = 1e-300
    // holds for pr near 0.999, but the search on [0..1] ends at label 240.
    let err = Binomial::search_pr(1e-300, 1.0, 100, 0).unwrap_err();
    assert_eq!(
        err,
        BinomialError::Search(SearchError::AnswerAboveUpperBound { bound: 1.0 })
    );
}

#[test]
fn negative_binomial_search_pr_dzror_label_240_failure_is_status_2() {
    // cdfnbn which=4, cdflib.f90:5616-5634: the search for ompr with
    // q = 1e-300 ends at label 240.
    let err = NegativeBinomial::search_pr(1.0, 1e-300, 5, 100).unwrap_err();
    assert_eq!(
        err,
        NegativeBinomialError::Search(SearchError::AnswerAboveUpperBound { bound: 1.0 })
    );
}
