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

use cdflib::{
    FisherSnedecorNoncentral, FisherSnedecorNoncentralError, SearchError, StudentsT, StudentsTError,
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
    // Trigger qleft with f = 0.5, dfn = 10, ncp = 0, p = 0.5: the central-F
    // CDF at f = 0.5, dfn = 10, dfd → ∞ approaches χ²(10)/10 ≤ 0.5 ≈ 0.0083,
    // and at small dfd it's larger. cum(dfd=1) > 0.5, so f = cum - p is
    // positive at small. Increasing dfd from 1 reduces cum (concentrates
    // around 1.0 ≥ 0.5), so f stays positive somewhere along the way then
    // can cross. Actually: cum here is decreasing in dfd. We need both
    // endpoints same-signed. With p = 0.99: cum(f=0.5, dfn=10, dfd=1)
    // is well below 0.99; cum(f=0.5, dfn=10, dfd→∞) approaches ≈ 0.008.
    // So f = cum - 0.99 is always negative. f is decreasing → qleft.
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
