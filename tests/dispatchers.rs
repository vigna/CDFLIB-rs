#![cfg(not(miri))]

//! Reference-table tests for every cdf* dispatcher in cdflib.f90.
//!
//! Each `tests/data/cdf<dist>_<solved-var>.csv` file is the output of
//! calling the corresponding Fortran dispatcher in `gen_dispatchers.f90`. The
//! columns are the inputs followed by the C-computed answer. These
//! tests assert that the corresponding Rust `search_*` /
//! `inverse_cdf` method returns the same value.
//!
//! Tolerance: both sides drive `dinvr` with CDFLIB's `tol = 1e-8`
//! (see `src/search/mod.rs`). Converged values therefore agree at the
//! 1e-8 floor of the algorithm, amplified by ~1/pdf at tail quantiles,
//! not to machine precision. See `REL` below.

mod common;

use cdflib::traits::ContinuousCdf;
use cdflib::{
    Beta, Binomial, ChiSquared, ChiSquaredNoncentral, FisherSnedecor, FisherSnedecorNoncentral,
    Gamma, NegativeBinomial, Normal, Poisson, StudentsT,
};
use common::{assert_close_eps, read_csv};

// Search converges at 1e-8 relative; tail quantiles inflate that by
// ~1/pdf. Measured max on this grid: 9.8e-7 for the most amplifying
// cases. 3e-6 leaves ~3x margin.
const REL: f64 = 3e-6;
// Some cdf* dispatchers compute parameters that are exactly 0 (e.g.,
// cdfnor's search_mean recovering mean=0). Rust converges to a denormal
// near-zero. 1e-14 absorbs that without masking real divergences.
const ABS: f64 = 1e-14;

// ============================================================== Beta

#[test]
fn cdfbet_x_matches_beta_inverse_cdf() {
    for row in read_csv("tests/data/cdfbet_x.csv") {
        let [p, _q, a, b, x_ref] = row[..] else {
            panic!("width")
        };
        let got = Beta::new(a, b).inverse_cdf(p).unwrap();
        assert_close_eps(got, x_ref, REL, ABS);
    }
}

#[test]
fn cdfbet_a_matches_beta_solve_a() {
    for row in read_csv("tests/data/cdfbet_a.csv") {
        let [p, q, x, b, a_ref] = row[..] else {
            panic!("width")
        };
        let got = Beta::search_a(p, q, x, b).unwrap();
        assert_close_eps(got, a_ref, REL, ABS);
    }
}

#[test]
fn cdfbet_b_matches_beta_solve_b() {
    for row in read_csv("tests/data/cdfbet_b.csv") {
        let [p, q, x, a, b_ref] = row[..] else {
            panic!("width")
        };
        let got = Beta::search_b(p, q, x, a).unwrap();
        assert_close_eps(got, b_ref, REL, ABS);
    }
}

// ============================================================ Binomial

#[test]
fn cdfbin_s_matches_binomial_inverse_ccdf() {
    for row in read_csv("tests/data/cdfbin_s.csv") {
        let [_p, q, xn, pr, s_ref] = row[..] else {
            panic!("width")
        };
        let got = Binomial::new(xn as u64, pr).inverse_ccdf(q).unwrap();
        assert_close_eps(got, s_ref, REL, ABS);
    }
}

#[test]
fn cdfbin_xn_matches_binomial_solve_trials() {
    for row in read_csv("tests/data/cdfbin_xn.csv") {
        let [p, q, s, pr, xn_ref] = row[..] else {
            panic!("width")
        };
        let got = Binomial::search_trials(p, q, pr, s as u64).unwrap();
        assert_close_eps(got, xn_ref, REL, ABS);
    }
}

#[test]
fn cdfbin_pr_matches_binomial_solve_pr() {
    for row in read_csv("tests/data/cdfbin_pr.csv") {
        let [p, q, s, xn, pr_ref] = row[..] else {
            panic!("width")
        };
        let got = Binomial::search_pr(p, q, xn as u64, s as u64).unwrap();
        assert_close_eps(got, pr_ref, REL, ABS);
    }
}

// ============================================================ ChiSquared

#[test]
fn cdfchi_x_matches_chi_squared_inverse_cdf() {
    for row in read_csv("tests/data/cdfchi_x.csv") {
        let [p, _q, df, x_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquared::new(df).inverse_cdf(p).unwrap();
        assert_close_eps(got, x_ref, REL, ABS);
    }
}

#[test]
fn cdfchi_df_matches_chi_squared_solve_df() {
    for row in read_csv("tests/data/cdfchi_df.csv") {
        let [p, q, x, df_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquared::search_df(p, q, x).unwrap();
        assert_close_eps(got, df_ref, REL, ABS);
    }
}

// =================================================== ChiSquaredNoncentral

#[test]
fn cdfchn_x_matches_chi_squared_noncentral_inverse_cdf() {
    for row in read_csv("tests/data/cdfchn_x.csv") {
        let [p, _q, df, pnonc, x_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquaredNoncentral::new(df, pnonc).inverse_cdf(p).unwrap();
        assert_close_eps(got, x_ref, REL, ABS);
    }
}

#[test]
fn cdfchn_df_matches_chi_squared_noncentral_solve_df() {
    for row in read_csv("tests/data/cdfchn_df.csv") {
        let [p, _q, x, pnonc, df_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquaredNoncentral::search_df(p, x, pnonc).unwrap();
        assert_close_eps(got, df_ref, REL, ABS);
    }
}

#[test]
fn cdfchn_pnonc_matches_chi_squared_noncentral_solve_ncp() {
    for row in read_csv("tests/data/cdfchn_pnonc.csv") {
        let [p, _q, x, df, ncp_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquaredNoncentral::search_ncp(p, x, df).unwrap();
        assert_close_eps(got, ncp_ref, REL, ABS);
    }
}

// ========================================================== FisherSnedecor

#[test]
fn cdff_f_matches_fisher_snedecor_inverse_cdf() {
    for row in read_csv("tests/data/cdff_f.csv") {
        let [p, _q, dfn, dfd, f_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecor::new(dfn, dfd).inverse_cdf(p).unwrap();
        assert_close_eps(got, f_ref, REL, ABS);
    }
}

#[test]
fn cdff_dfn_matches_fisher_snedecor_solve_dfn() {
    for row in read_csv("tests/data/cdff_dfn.csv") {
        let [p, q, f, dfd, dfn_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecor::search_dfn(p, q, f, dfd).unwrap();
        assert_close_eps(got, dfn_ref, REL, ABS);
    }
}

#[test]
fn cdff_dfd_matches_fisher_snedecor_solve_dfd() {
    for row in read_csv("tests/data/cdff_dfd.csv") {
        let [p, q, f, dfn, dfd_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecor::search_dfd(p, q, f, dfn).unwrap();
        assert_close_eps(got, dfd_ref, REL, ABS);
    }
}

// =================================================== FisherSnedecorNoncentral

#[test]
fn cdffnc_f_matches_fisher_snedecor_noncentral_inverse_cdf() {
    for row in read_csv("tests/data/cdffnc_f.csv") {
        let [p, _q, dfn, dfd, phonc, f_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::new(dfn, dfd, phonc)
            .inverse_cdf(p)
            .unwrap();
        assert_close_eps(got, f_ref, REL, ABS);
    }
}

#[test]
fn cdffnc_dfn_matches_fisher_snedecor_noncentral_solve_dfn() {
    for row in read_csv("tests/data/cdffnc_dfn.csv") {
        let [p, _q, f, dfd, phonc, dfn_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_dfn(p, f, dfd, phonc).unwrap();
        assert_close_eps(got, dfn_ref, REL, ABS);
    }
}

#[test]
fn cdffnc_dfd_matches_fisher_snedecor_noncentral_solve_dfd() {
    for row in read_csv("tests/data/cdffnc_dfd.csv") {
        let [p, _q, f, dfn, phonc, dfd_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_dfd(p, f, dfn, phonc).unwrap();
        assert_close_eps(got, dfd_ref, REL, ABS);
    }
}

#[test]
fn cdffnc_phonc_matches_fisher_snedecor_noncentral_solve_ncp() {
    for row in read_csv("tests/data/cdffnc_phonc.csv") {
        let [p, _q, f, dfn, dfd, ncp_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_ncp(p, f, dfn, dfd).unwrap();
        assert_close_eps(got, ncp_ref, REL, ABS);
    }
}

// ============================================================== Gamma

// CDFLIB's cdfgam names its second parameter scale, but the code
// computes P(shape, x * scale), so it's mathematically the rate.
// This crate calls the parameter rate (see src/distribution/gamma.rs);
// we pass the CSV's scale column directly to Gamma::new as rate.

#[test]
fn cdfgam_x_matches_gamma_inverse_cdf() {
    for row in read_csv("tests/data/cdfgam_x.csv") {
        let [p, _q, shape, rate, x_ref] = row[..] else {
            panic!("width")
        };
        let got = Gamma::new(shape, rate).inverse_cdf(p).unwrap();
        assert_close_eps(got, x_ref, REL, ABS);
    }
}

#[test]
fn cdfgam_shape_matches_gamma_solve_shape() {
    for row in read_csv("tests/data/cdfgam_shape.csv") {
        let [p, q, x, rate, shape_ref] = row[..] else {
            panic!("width")
        };
        let got = Gamma::search_shape(p, q, x, rate).unwrap();
        assert_close_eps(got, shape_ref, REL, ABS);
    }
}

#[test]
fn cdfgam_scale_matches_gamma_solve_rate() {
    for row in read_csv("tests/data/cdfgam_scale.csv") {
        let [p, q, x, shape, rate_ref] = row[..] else {
            panic!("width")
        };
        let got = Gamma::search_rate(p, q, x, shape).unwrap();
        assert_close_eps(got, rate_ref, REL, ABS);
    }
}

// ========================================================= NegativeBinomial

#[test]
fn cdfnbn_s_matches_negative_binomial_inverse_ccdf() {
    for row in read_csv("tests/data/cdfnbn_s.csv") {
        let [_p, q, r, pr, s_ref] = row[..] else {
            panic!("width")
        };
        let got = NegativeBinomial::new(r as u64, pr).inverse_ccdf(q).unwrap();
        assert_close_eps(got, s_ref, REL, ABS);
    }
}

#[test]
fn cdfnbn_xn_matches_negative_binomial_solve_r() {
    for row in read_csv("tests/data/cdfnbn_xn.csv") {
        let [p, q, s, pr, xn_ref] = row[..] else {
            panic!("width")
        };
        let got = NegativeBinomial::search_r(p, q, pr, s as u64).unwrap();
        assert_close_eps(got, xn_ref, REL, ABS);
    }
}

#[test]
fn cdfnbn_pr_matches_negative_binomial_solve_pr() {
    for row in read_csv("tests/data/cdfnbn_pr.csv") {
        let [p, q, s, xn, pr_ref] = row[..] else {
            panic!("width")
        };
        let got = NegativeBinomial::search_pr(p, q, xn as u64, s as u64).unwrap();
        assert_close_eps(got, pr_ref, REL, ABS);
    }
}

// ============================================================== Normal

#[test]
fn cdfnor_x_matches_normal_inverse_cdf() {
    for row in read_csv("tests/data/cdfnor_x.csv") {
        let [p, _q, mean, sd, x_ref] = row[..] else {
            panic!("width")
        };
        let got = Normal::new(mean, sd).inverse_cdf(p).unwrap();
        assert_close_eps(got, x_ref, REL, ABS);
    }
}

#[test]
fn cdfnor_mean_matches_normal_solve_mean() {
    for row in read_csv("tests/data/cdfnor_mean.csv") {
        let [p, q, x, sd, mean_ref] = row[..] else {
            panic!("width")
        };
        let got = Normal::search_mean(p, q, x, sd).unwrap();
        assert_close_eps(got, mean_ref, REL, ABS);
    }
}

#[test]
fn cdfnor_sd_matches_normal_solve_sd() {
    for row in read_csv("tests/data/cdfnor_sd.csv") {
        let [p, q, x, mean, sd_ref] = row[..] else {
            panic!("width")
        };
        let got = Normal::search_sd(p, q, x, mean).unwrap();
        assert_close_eps(got, sd_ref, REL, ABS);
    }
}

// ============================================================== Poisson

#[test]
fn cdfpoi_s_matches_poisson_inverse_ccdf() {
    for row in read_csv("tests/data/cdfpoi_s.csv") {
        let [_p, q, lambda, s_ref] = row[..] else {
            panic!("width")
        };
        let got = Poisson::new(lambda).inverse_ccdf(q).unwrap();
        assert_close_eps(got, s_ref, REL, ABS);
    }
}

#[test]
fn cdfpoi_xlam_matches_poisson_solve_lambda() {
    for row in read_csv("tests/data/cdfpoi_xlam.csv") {
        let [p, q, s, xlam_ref] = row[..] else {
            panic!("width")
        };
        let got = Poisson::search_lambda(p, q, s as u64).unwrap();
        assert_close_eps(got, xlam_ref, REL, ABS);
    }
}

// ============================================================== StudentsT

#[test]
fn cdft_t_matches_students_t_inverse_cdf() {
    for row in read_csv("tests/data/cdft_t.csv") {
        let [p, _q, df, t_ref] = row[..] else {
            panic!("width")
        };
        let got = StudentsT::new(df).inverse_cdf(p).unwrap();
        assert_close_eps(got, t_ref, REL, ABS);
    }
}

#[test]
fn cdft_df_matches_students_t_solve_df() {
    for row in read_csv("tests/data/cdft_df.csv") {
        let [p, q, t, df_ref] = row[..] else {
            panic!("width")
        };
        let got = StudentsT::search_df(p, q, t).unwrap();
        assert_close_eps(got, df_ref, REL, ABS);
    }
}
