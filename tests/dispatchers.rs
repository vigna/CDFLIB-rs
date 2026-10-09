#![cfg(not(miri))]

//! Reference-table tests for every cdf* dispatcher in cdflib.f90.
//!
//! Each `tests/data/cdf<dist>_<solved-var>.csv` file is the output of
//! calling the corresponding Fortran dispatcher in `gen_dispatchers.f90`. The
//! columns are the inputs followed by the F90-computed answer. These
//! tests assert that the corresponding Rust `search_*` /
//! `inverse_cdf` / `inverse_ccdf` method returns the same value, bit for
//! bit: the solver and the cumulative distribution functions are exact
//! ports, so the whole search follows the F90 iteration.

mod common;

use cdflib::traits::ContinuousCdf;
use cdflib::{
    Beta, Binomial, ChiSquared, ChiSquaredNoncentral, FisherSnedecor, FisherSnedecorNoncentral,
    Gamma, NegativeBinomial, Normal, NormalError, Poisson, StudentsT,
};
use common::{assert_exact, read_csv};

// ============================================================== Beta

#[test]
fn cdfbet_x_matches_beta_inverse_cdf() {
    for row in read_csv("tests/data/cdfbet_x.csv") {
        let [p, q, a, b, x_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = Beta::new(a, b);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, x_ref, &row);
    }
}

#[test]
fn cdfbet_a_matches_beta_search_a() {
    for row in read_csv("tests/data/cdfbet_a.csv") {
        let [p, q, x, b, a_ref] = row[..] else {
            panic!("width")
        };
        let got = Beta::search_a(p, q, x, b).unwrap();
        assert_exact(got, a_ref, &row);
    }
}

#[test]
fn cdfbet_b_matches_beta_search_b() {
    for row in read_csv("tests/data/cdfbet_b.csv") {
        let [p, q, x, a, b_ref] = row[..] else {
            panic!("width")
        };
        let got = Beta::search_b(p, q, x, a).unwrap();
        assert_exact(got, b_ref, &row);
    }
}

// ============================================================ Binomial

#[test]
fn cdfbin_s_matches_binomial_inverse_ccdf() {
    for row in read_csv("tests/data/cdfbin_s.csv") {
        let [p, q, xn, pr, s_ref] = row[..] else {
            panic!("width")
        };
        // The Rust API derives p = 1 - q; rows where F90 received a
        // different p are covered with exact complements by
        // tests/dispatcher_calls.rs.
        if (1.0 - q).to_bits() != p.to_bits() && p <= q {
            continue;
        }
        let got = Binomial::new(xn as u64, pr).inverse_ccdf(q).unwrap();
        assert_exact(got, s_ref, &row);
    }
}

#[test]
fn cdfbin_xn_matches_binomial_search_trials() {
    for row in read_csv("tests/data/cdfbin_xn.csv") {
        let [p, q, s, pr, xn_ref] = row[..] else {
            panic!("width")
        };
        let got = Binomial::search_trials(p, q, pr, s as u64).unwrap();
        assert_exact(got, xn_ref, &row);
    }
}

#[test]
fn cdfbin_pr_matches_binomial_search_pr() {
    for row in read_csv("tests/data/cdfbin_pr.csv") {
        let [p, q, s, xn, pr_ref] = row[..] else {
            panic!("width")
        };
        let got = Binomial::search_pr(p, q, xn as u64, s as u64).unwrap();
        assert_exact(got, pr_ref, &row);
    }
}

// ============================================================ ChiSquared

#[test]
fn cdfchi_x_matches_chi_squared_inverse_cdf() {
    for row in read_csv("tests/data/cdfchi_x.csv") {
        let [p, q, df, x_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = ChiSquared::new(df);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, x_ref, &row);
    }
}

#[test]
fn cdfchi_df_matches_chi_squared_search_df() {
    for row in read_csv("tests/data/cdfchi_df.csv") {
        let [p, q, x, df_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquared::search_df(p, q, x).unwrap();
        assert_exact(got, df_ref, &row);
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
        assert_exact(got, x_ref, &row);
    }
}

#[test]
fn cdfchn_df_matches_chi_squared_noncentral_search_df() {
    for row in read_csv("tests/data/cdfchn_df.csv") {
        let [p, _q, x, pnonc, df_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquaredNoncentral::search_df(p, x, pnonc).unwrap();
        assert_exact(got, df_ref, &row);
    }
}

#[test]
fn cdfchn_pnonc_matches_chi_squared_noncentral_search_ncp() {
    for row in read_csv("tests/data/cdfchn_pnonc.csv") {
        let [p, _q, x, df, ncp_ref] = row[..] else {
            panic!("width")
        };
        let got = ChiSquaredNoncentral::search_ncp(p, x, df).unwrap();
        assert_exact(got, ncp_ref, &row);
    }
}

// ========================================================== FisherSnedecor

#[test]
fn cdff_f_matches_fisher_snedecor_inverse_cdf() {
    for row in read_csv("tests/data/cdff_f.csv") {
        let [p, q, dfn, dfd, f_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = FisherSnedecor::new(dfn, dfd);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, f_ref, &row);
    }
}

#[test]
fn cdff_dfn_matches_fisher_snedecor_search_dfn() {
    for row in read_csv("tests/data/cdff_dfn.csv") {
        let [p, q, f, dfd, dfn_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecor::search_dfn(p, q, f, dfd).unwrap();
        assert_exact(got, dfn_ref, &row);
    }
}

#[test]
fn cdff_dfd_matches_fisher_snedecor_search_dfd() {
    for row in read_csv("tests/data/cdff_dfd.csv") {
        let [p, q, f, dfn, dfd_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecor::search_dfd(p, q, f, dfn).unwrap();
        assert_exact(got, dfd_ref, &row);
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
        assert_exact(got, f_ref, &row);
    }
}

#[test]
fn cdffnc_dfn_matches_fisher_snedecor_noncentral_search_dfn() {
    for row in read_csv("tests/data/cdffnc_dfn.csv") {
        let [p, _q, f, dfd, phonc, dfn_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_dfn(p, f, dfd, phonc).unwrap();
        assert_exact(got, dfn_ref, &row);
    }
}

#[test]
fn cdffnc_dfd_matches_fisher_snedecor_noncentral_search_dfd() {
    for row in read_csv("tests/data/cdffnc_dfd.csv") {
        let [p, _q, f, dfn, phonc, dfd_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_dfd(p, f, dfn, phonc).unwrap();
        assert_exact(got, dfd_ref, &row);
    }
}

#[test]
fn cdffnc_pnonc_matches_fisher_snedecor_noncentral_search_ncp() {
    for row in read_csv("tests/data/cdffnc_phonc.csv") {
        let [p, _q, f, dfn, dfd, ncp_ref] = row[..] else {
            panic!("width")
        };
        let got = FisherSnedecorNoncentral::search_ncp(p, f, dfn, dfd).unwrap();
        assert_exact(got, ncp_ref, &row);
    }
}

// ============================================================== Gamma

// CDFLIB's cdfgam names its second parameter scale, but the code
// computes P(shape, x * scale), so it's mathematically the rate.
// This crate calls the parameter rate (see src/dist/gamma.rs);
// we pass the CSV's scale column directly to Gamma::new as rate.

#[test]
fn cdfgam_x_matches_gamma_inverse_cdf() {
    for row in read_csv("tests/data/cdfgam_x.csv") {
        let [p, q, shape, rate, x_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = Gamma::new(shape, rate);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, x_ref, &row);
    }
}

#[test]
fn cdfgam_shape_matches_gamma_search_shape() {
    for row in read_csv("tests/data/cdfgam_shape.csv") {
        let [p, q, x, rate, shape_ref] = row[..] else {
            panic!("width")
        };
        let got = Gamma::search_shape(p, q, x, rate).unwrap();
        assert_exact(got, shape_ref, &row);
    }
}

#[test]
fn cdfgam_scale_matches_gamma_search_rate() {
    for row in read_csv("tests/data/cdfgam_scale.csv") {
        let [p, q, x, shape, rate_ref] = row[..] else {
            panic!("width")
        };
        let got = Gamma::search_rate(p, q, x, shape).unwrap();
        assert_exact(got, rate_ref, &row);
    }
}

// ========================================================= NegativeBinomial

#[test]
fn cdfnbn_s_matches_negative_binomial_inverse_ccdf() {
    for row in read_csv("tests/data/cdfnbn_s.csv") {
        let [p, q, r, pr, s_ref] = row[..] else {
            panic!("width")
        };
        // The Rust API derives p = 1 - q; rows where F90 received a
        // different p are covered with exact complements by
        // tests/dispatcher_calls.rs.
        if (1.0 - q).to_bits() != p.to_bits() && p <= q {
            continue;
        }
        let got = NegativeBinomial::new(r as u64, pr).inverse_ccdf(q).unwrap();
        assert_exact(got, s_ref, &row);
    }
}

#[test]
fn cdfnbn_xn_matches_negative_binomial_search_r() {
    for row in read_csv("tests/data/cdfnbn_xn.csv") {
        let [p, q, s, pr, xn_ref] = row[..] else {
            panic!("width")
        };
        let got = NegativeBinomial::search_r(p, q, pr, s as u64).unwrap();
        assert_exact(got, xn_ref, &row);
    }
}

#[test]
fn cdfnbn_pr_matches_negative_binomial_search_pr() {
    for row in read_csv("tests/data/cdfnbn_pr.csv") {
        let [p, q, s, xn, pr_ref] = row[..] else {
            panic!("width")
        };
        let got = NegativeBinomial::search_pr(p, q, xn as u64, s as u64).unwrap();
        assert_exact(got, pr_ref, &row);
    }
}

// ============================================================== Normal

#[test]
fn cdfnor_x_matches_normal_inverse_cdf() {
    for row in read_csv("tests/data/cdfnor_x.csv") {
        let [p, q, mean, sd, x_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = Normal::new(mean, sd);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, x_ref, &row);
    }
}

#[test]
fn cdfnor_mean_matches_normal_search_mean() {
    for row in read_csv("tests/data/cdfnor_mean.csv") {
        let [p, q, x, sd, mean_ref] = row[..] else {
            panic!("width")
        };
        let got = Normal::search_mean(p, q, x, sd).unwrap();
        assert_exact(got, mean_ref, &row);
    }
}

#[test]
fn cdfnor_sd_matches_normal_search_sd() {
    for row in read_csv("tests/data/cdfnor_sd.csv") {
        let [p, q, x, mean, sd_ref] = row[..] else {
            panic!("width")
        };
        match Normal::search_sd(p, q, x, mean) {
            Ok(got) => assert_exact(got, sd_ref, &row),
            // Rust only: an sd that is not positive, which F90 returns, is
            // reported as SdNotPositive.
            Err(NormalError::SdNotPositive(got)) => {
                assert!(sd_ref <= 0.0, "{row:?}");
                assert_exact(got, sd_ref, &row);
            }
            Err(e) => panic!("{row:?}: {e}"),
        }
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
        assert_exact(got, s_ref, &row);
    }
}

#[test]
fn cdfpoi_xlam_matches_poisson_search_lambda() {
    for row in read_csv("tests/data/cdfpoi_xlam.csv") {
        let [p, q, s, xlam_ref] = row[..] else {
            panic!("width")
        };
        let got = Poisson::search_lambda(p, q, s as u64).unwrap();
        assert_exact(got, xlam_ref, &row);
    }
}

// ============================================================== StudentsT

#[test]
fn cdft_t_matches_students_t_inverse_cdf() {
    for row in read_csv("tests/data/cdft_t.csv") {
        let [p, q, df, t_ref] = row[..] else {
            panic!("width")
        };
        // F90 uses p in the search when p <= q and q otherwise.
        let d = StudentsT::new(df);
        let got = if p <= q {
            d.inverse_cdf(p).unwrap()
        } else {
            d.inverse_ccdf(q).unwrap()
        };
        assert_exact(got, t_ref, &row);
    }
}

#[test]
fn cdft_df_matches_students_t_search_df() {
    for row in read_csv("tests/data/cdft_df.csv") {
        let [p, q, t, df_ref] = row[..] else {
            panic!("width")
        };
        let got = StudentsT::search_df(p, q, t).unwrap();
        assert_exact(got, df_ref, &row);
    }
}
