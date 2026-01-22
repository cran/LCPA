#include <RcppArmadillo.h>
#include <unordered_map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

std::unordered_map<std::string, arma::umat> path_cache;

// [[Rcpp::export]]
arma::umat make_latent_paths_cpp(int L, int times) {
  std::string key = "L" + std::to_string(L) + "_T" + std::to_string(times);

  auto it = path_cache.find(key);
  if (it != path_cache.end()) {
    return it->second;
  }

  size_t n_paths = 1;
  for (int i = 0; i < times; i++) {
    n_paths *= L;
  }

  arma::umat paths(n_paths, times, arma::fill::zeros);
  arma::uvec counter(times, arma::fill::zeros);

  for (size_t i = 0; i < n_paths; i++) {
    paths.row(i) = counter.t();

    int j = 0;
    while (j < times && counter(j) == static_cast<unsigned int>(L - 1)) {
      counter(j) = 0;
      j++;
    }
    if (j < times) {
      counter(j)++;
    }
  }

  path_cache[key] = paths;
  return paths;
}

// [[Rcpp::export]]
List lta_vector_to_parameters_cpp(NumericVector params,
                                  List covariates_list,
                                  int L,
                                  int ref_class) {
  int times = covariates_list.size();
  NumericMatrix cov1 = Rcpp::as<NumericMatrix>(covariates_list[0]);
  int p1 = cov1.ncol();
  int num_beta = p1 * (L - 1);

  arma::mat beta(p1, L, arma::fill::zeros);
  std::vector<int> non_ref_classes;
  for (int i = 1; i <= L; i++) {
    if (i != ref_class) non_ref_classes.push_back(i - 1);
  }

  if (num_beta > 0 && non_ref_classes.size() > 0) {
    for (size_t j = 0; j < non_ref_classes.size(); j++) {
      int col = non_ref_classes[j];
      for (int i = 0; i < p1; i++) {
        if (i + static_cast<int>(j) * p1 < params.size()) {
          beta(i, col) = params[i + j * p1];
        }
      }
    }
  }

  List gamma(times - 1);
  int idx = num_beta;

  for (int t_idx = 1; t_idx < times; t_idx++) {
    NumericMatrix cov_t = Rcpp::as<NumericMatrix>(covariates_list[t_idx]);
    int pt = cov_t.ncol();
    int size = L * (L - 1) * pt;

    List g_t(L);
    NumericVector slice(size > 0 ? size : 0);

    if (size > 0 && idx + size <= params.size()) {
      for (int i = 0; i < size; i++) {
        slice[i] = params[idx + i];
      }
    }
    idx += size;

    int pos = 0;
    for (int l = 0; l < L; l++) {
      List g_t_l(L);
      for (int ll = 0; ll < L; ll++) {
        if (ll + 1 == ref_class) {
          g_t_l[ll] = NumericVector(pt, 0.0);
        } else {
          if (pt > 0 && pos + pt <= slice.size()) {
            NumericVector coef(pt);
            for (int k = 0; k < pt; k++) {
              coef[k] = slice[pos + k];
            }
            g_t_l[ll] = coef;
            pos += pt;
          } else {
            g_t_l[ll] = NumericVector(pt, 0.0);
          }
        }
      }
      g_t[l] = g_t_l;
    }

    if (t_idx - 1 < gamma.size()) {
      gamma[t_idx - 1] = g_t;
    }
  }

  return List::create(_["beta"] = beta, _["gamma"] = gamma);
}

// [[Rcpp::export]]
double get_log_lik_lta_optim_cpp(NumericVector init_par,
                                 List CEP_list,
                                 List Zs_list,
                                 List covariates_list,
                                 bool covariates_timeCross,
                                 int ref_class,
                                 arma::umat latent_paths) {

  NumericMatrix CEP0 = Rcpp::as<NumericMatrix>(CEP_list[0]);
  int L = CEP0.ncol();
  int N = Rcpp::as<IntegerVector>(Zs_list[0]).size();
  int times = covariates_list.size();

  if ((int)latent_paths.n_cols != times) {
    stop("latent_paths dimensions do not match number of time points");
  }

  List params = lta_vector_to_parameters_cpp(init_par, covariates_list, L, ref_class);
  arma::mat beta = Rcpp::as<arma::mat>(params["beta"]);
  List gamma = params["gamma"];

  if (covariates_timeCross && times > 2) {
    for (int t = 2; t < times; t++) {
      if (t-1 < gamma.size()) {
        gamma[t-1] = gamma[t-2];
      }
    }
  }

  int n_paths = latent_paths.n_rows;
  double log_lik = 0.0;
  const double log_min = std::log(1e-200);

  std::vector<arma::mat> CEP_vec(times);
  std::vector<arma::uvec> Zs_vec(times);
  std::vector<arma::mat> cov_vec(times);

  for (int t = 0; t < times; t++) {
    NumericMatrix CEP_t = Rcpp::as<NumericMatrix>(CEP_list[t]);
    CEP_vec[t] = arma::mat(CEP_t.begin(), CEP_t.nrow(), CEP_t.ncol(), false);

    IntegerVector Zs_t = Rcpp::as<IntegerVector>(Zs_list[t]);
    arma::uvec tmp(Zs_t.size());
    for (int i = 0; i < Zs_t.size(); i++) {
      tmp(i) = static_cast<unsigned int>(Zs_t[i] - 1);
    }
    Zs_vec[t] = tmp;

    NumericMatrix cov_t = Rcpp::as<NumericMatrix>(covariates_list[t]);
    cov_vec[t] = arma::mat(cov_t.begin(), cov_t.nrow(), cov_t.ncol(), false);
  }

  arma::mat zeta(L, L);
  std::vector<arma::mat> P_t(times);

  for (int n = 0; n < N; n++) {
    arma::rowvec Xt1 = cov_vec[0].row(n);
    arma::rowvec eta1 = Xt1 * beta;

    double max_eta1 = eta1.max();
    eta1 = eta1 - max_eta1;
    arma::rowvec exp_eta1 = arma::exp(eta1);
    double sum_exp_eta1 = arma::sum(exp_eta1);

    if (sum_exp_eta1 <= 0 || !std::isfinite(sum_exp_eta1)) {
      sum_exp_eta1 = std::numeric_limits<double>::min();
    }

    arma::rowvec P1 = exp_eta1 / sum_exp_eta1;
    P_t[0] = P1;

    for (int t_idx = 1; t_idx < times; t_idx++) {
      arma::rowvec Xt = cov_vec[t_idx].row(n);
      int current_time_idx = t_idx - 1;

      zeta.fill(0.0);
      for (int l = 0; l < L; l++) {
        if (current_time_idx >= 0 && current_time_idx < static_cast<int>(gamma.size())) {
          List g_t = Rcpp::as<List>(gamma[current_time_idx]);
          if (l >= 0 && l < static_cast<int>(g_t.size())) {
            List g_t_l = Rcpp::as<List>(g_t[l]);
            for (int ll = 0; ll < L; ll++) {
              if (ll >= 0 && ll < static_cast<int>(g_t_l.size())) {
                NumericVector coef = Rcpp::as<NumericVector>(g_t_l[ll]);
                if (coef.size() > 0) {
                  arma::rowvec coef_vec(coef.begin(), coef.size(), false);
                  zeta(l, ll) = arma::dot(Xt, coef_vec);
                }
              }
            }
          }
        }
      }

      arma::mat P_t_mat(L, L);
      for (int l = 0; l < L; l++) {
        arma::rowvec row = zeta.row(l);
        double max_val = row.max();
        row = row - max_val;
        arma::rowvec exp_row = arma::exp(row);
        double sum_exp_row = arma::sum(exp_row);

        if (sum_exp_row <= 0 || !std::isfinite(sum_exp_row)) {
          sum_exp_row = std::numeric_limits<double>::min();
        }

        P_t_mat.row(l) = exp_row / sum_exp_row;
      }
      P_t[t_idx] = P_t_mat;
    }

    double total_prob = 0.0;

    for (int k = 0; k < n_paths; k++) {
      arma::urowvec path = latent_paths.row(k);
      double prob = 1.0;
      double weight = 1.0;

      int class1 = path(0);
      if (class1 < 0 || class1 >= L) continue;

      prob *= P_t[0](class1);

      if (static_cast<size_t>(Zs_vec[0](n)) < CEP_vec[0].n_rows &&
          static_cast<size_t>(class1) < CEP_vec[0].n_cols) {
        weight *= CEP_vec[0](Zs_vec[0](n), class1);
      } else {
        weight = 0.0;
      }

      for (int t_idx = 1; t_idx < times; t_idx++) {
        int prev_class = path(t_idx - 1);
        int curr_class = path(t_idx);

        if (prev_class < 0 || prev_class >= L || curr_class < 0 || curr_class >= L) {
          prob = 0.0;
          weight = 0.0;
          break;
        }

        if (t_idx < static_cast<int>(P_t.size()) &&
            static_cast<size_t>(prev_class) < P_t[t_idx].n_rows &&
            static_cast<size_t>(curr_class) < P_t[t_idx].n_cols) {
          prob *= P_t[t_idx](prev_class, curr_class);
        } else {
          prob = 0.0;
        }

        if (static_cast<size_t>(Zs_vec[t_idx](n)) < CEP_vec[t_idx].n_rows &&
            static_cast<size_t>(curr_class) < CEP_vec[t_idx].n_cols) {
          weight *= CEP_vec[t_idx](Zs_vec[t_idx](n), curr_class);
        } else {
          weight = 0.0;
        }

        if (prob <= 0 || weight <= 0) break;
      }

      total_prob += prob * weight;
    }

    if (total_prob < 1e-200 || !std::isfinite(total_prob)) {
      log_lik += log_min;
    } else {
      log_lik += std::log(total_prob);
    }
  }

  if (!std::isfinite(log_lik)) {
    warning("Non-finite log-likelihood detected: %f", log_lik);
    return 1e10;
  }

  return -log_lik;
}
