#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace {

std::unordered_map<std::string, arma::umat> path_cache;

arma::rowvec softmax(const arma::rowvec& linear_predictor) {
  double maximum = linear_predictor.max();
  arma::rowvec probability = arma::exp(linear_predictor - maximum);
  double total = arma::sum(probability);
  if (!std::isfinite(total) || total <= 0.0) {
    probability.fill(1.0 / linear_predictor.n_elem);
    return probability;
  }
  return probability / total;
}

}

// [[Rcpp::export]]
arma::umat make_latent_paths_cpp(int L, int times) {
  std::string key = "L" + std::to_string(L) + "_T" + std::to_string(times);
  auto cached = path_cache.find(key);
  if (cached != path_cache.end()) return cached->second;

  std::size_t n_paths = 1;
  for (int t = 0; t < times; ++t) n_paths *= L;

  arma::umat paths(n_paths, times, arma::fill::zeros);
  arma::uvec counter(times, arma::fill::zeros);
  for (std::size_t path = 0; path < n_paths; ++path) {
    paths.row(path) = counter.t();
    int t = 0;
    while (t < times && counter[t] == static_cast<unsigned int>(L - 1)) {
      counter[t] = 0;
      ++t;
    }
    if (t < times) ++counter[t];
  }

  path_cache[key] = paths;
  return paths;
}

// [[Rcpp::export]]
List lta_vector_to_parameters_cpp(const NumericVector& params,
                                  const List& covariates_list,
                                  int L,
                                  int ref_class) {
  int times = covariates_list.size();
  NumericMatrix covariates_t1 = as<NumericMatrix>(covariates_list[0]);
  int p1 = covariates_t1.ncol();
  std::vector<int> free_classes;
  for (int l = 0; l < L; ++l) {
    if (l != ref_class - 1) free_classes.push_back(l);
  }

  arma::mat beta(p1, L, arma::fill::zeros);
  int position = 0;
  for (int current_class : free_classes) {
    for (int j = 0; j < p1; ++j) beta(j, current_class) = params[position++];
  }

  List gamma(std::max(times - 1, 0));
  for (int t = 1; t < times; ++t) {
    NumericMatrix covariates_t = as<NumericMatrix>(covariates_list[t]);
    int p = covariates_t.ncol();
    List gamma_t(L);
    for (int from_class = 0; from_class < L; ++from_class) {
      List gamma_from(L);
      for (int to_class = 0; to_class < L; ++to_class) {
        NumericVector coefficient(p);
        if (to_class != ref_class - 1) {
          for (int j = 0; j < p; ++j) coefficient[j] = params[position++];
        }
        gamma_from[to_class] = coefficient;
      }
      gamma_t[from_class] = gamma_from;
    }
    gamma[t - 1] = gamma_t;
  }

  if (position != params.size()) stop("Parameter vector has an incompatible length");
  return List::create(_["beta"] = beta, _["gamma"] = gamma);
}

// [[Rcpp::export]]
List get_log_lik_lta_optim_cpp(const NumericVector& params,
                               const List& CEP_list,
                               const List& Zs_list,
                               const List& covariates_list,
                               bool covariates_time_cross,
                               int ref_class,
                               bool compute_gradient) {
  int times = covariates_list.size();
  if (times < 1 || CEP_list.size() != times || Zs_list.size() != times) {
    stop("CEP, Zs, and covariates must contain the same time points");
  }

  NumericMatrix CEP_t1 = as<NumericMatrix>(CEP_list[0]);
  IntegerVector Z_t1 = as<IntegerVector>(Zs_list[0]);
  int L = CEP_t1.nrow();
  int N = Z_t1.size();
  if (CEP_t1.ncol() != L || ref_class < 1 || ref_class > L) {
    stop("Invalid LTA class dimensions");
  }
  static_cast<void>(covariates_time_cross);

  std::vector<int> free_classes;
  for (int l = 0; l < L; ++l) {
    if (l != ref_class - 1) free_classes.push_back(l);
  }

  std::vector<arma::mat> CEP(times);
  std::vector<arma::uvec> observed(times);
  std::vector<arma::mat> covariates(times);
  for (int t = 0; t < times; ++t) {
    NumericMatrix CEP_current = as<NumericMatrix>(CEP_list[t]);
    IntegerVector Z_current = as<IntegerVector>(Zs_list[t]);
    NumericMatrix covariates_current = as<NumericMatrix>(covariates_list[t]);
    if (CEP_current.nrow() != L || CEP_current.ncol() != L ||
        Z_current.size() != N || covariates_current.nrow() != N) {
      stop("Incompatible LTA dimensions across time points");
    }

    CEP[t] = arma::mat(CEP_current.begin(), L, L, false);
    covariates[t] = arma::mat(
      covariates_current.begin(), N, covariates_current.ncol(), false
    );
    observed[t].set_size(N);
    for (int n = 0; n < N; ++n) {
      int current = Z_current[n] - 1;
      if (current < 0 || current >= L) stop("Zs contains an invalid class");
      observed[t][n] = current;
    }
  }

  int position = 0;
  int p1 = covariates[0].n_cols;
  arma::mat beta(p1, L, arma::fill::zeros);
  for (int current_class : free_classes) {
    for (int j = 0; j < p1; ++j) beta(j, current_class) = params[position++];
  }

  std::vector<arma::cube> gamma(std::max(times - 1, 0));
  std::vector<int> gamma_position(std::max(times - 1, 0));
  for (int t = 1; t < times; ++t) {
    int p = covariates[t].n_cols;
    gamma[t - 1].zeros(p, L, L);
    gamma_position[t - 1] = position;
    for (int from_class = 0; from_class < L; ++from_class) {
      for (int to_class : free_classes) {
        for (int j = 0; j < p; ++j) {
          gamma[t - 1](j, to_class, from_class) = params[position++];
        }
      }
    }
  }
  if (position != params.size()) stop("Parameter vector has an incompatible length");

  double log_likelihood = 0.0;
  arma::vec score(params.size(), arma::fill::zeros);
  std::vector<arma::vec> emission(times);
  std::vector<arma::vec> forward(times);
  std::vector<arma::vec> backward(times);
  std::vector<arma::mat> transition(std::max(times - 1, 0));
  std::vector<double> scale(times);

  for (int n = 0; n < N; ++n) {
    arma::rowvec initial_probability = softmax(covariates[0].row(n) * beta);
    for (int t = 0; t < times; ++t) emission[t] = CEP[t].col(observed[t][n]);

    for (int t = 1; t < times; ++t) {
      transition[t - 1].set_size(L, L);
      for (int from_class = 0; from_class < L; ++from_class) {
        arma::rowvec linear_predictor(L, arma::fill::zeros);
        for (int to_class : free_classes) {
          linear_predictor[to_class] = arma::dot(
            covariates[t].row(n), gamma[t - 1].slice(from_class).col(to_class).t()
          );
        }
        transition[t - 1].row(from_class) = softmax(linear_predictor);
      }
    }

    forward[0] = initial_probability.t() % emission[0];
    scale[0] = arma::sum(forward[0]);
    bool valid = std::isfinite(scale[0]) && scale[0] > 0.0;
    if (valid) forward[0] /= scale[0];

    for (int t = 1; t < times && valid; ++t) {
      forward[t] = emission[t] % (transition[t - 1].t() * forward[t - 1]);
      scale[t] = arma::sum(forward[t]);
      valid = std::isfinite(scale[t]) && scale[t] > 0.0;
      if (valid) forward[t] /= scale[t];
    }

    if (!valid) {
      arma::vec invalid_gradient(params.size(), arma::fill::zeros);
      return List::create(
        _["objective"] = R_PosInf,
        _["gradient"] = compute_gradient ? wrap(invalid_gradient) : R_NilValue
      );
    }
    for (int t = 0; t < times; ++t) log_likelihood += std::log(scale[t]);
    if (!compute_gradient) continue;

    backward[times - 1].ones(L);
    for (int t = times - 2; t >= 0; --t) {
      backward[t] = transition[t] * (emission[t + 1] % backward[t + 1]);
      backward[t] /= scale[t + 1];
    }

    arma::vec initial_posterior = forward[0] % backward[0];
    initial_posterior /= arma::sum(initial_posterior);
    for (std::size_t j = 0; j < free_classes.size(); ++j) {
      int current_class = free_classes[j];
      double residual = initial_posterior[current_class] - initial_probability[current_class];
      for (int v = 0; v < p1; ++v) {
        score[static_cast<int>(j) * p1 + v] += covariates[0](n, v) * residual;
      }
    }

    for (int t = 1; t < times; ++t) {
      arma::vec destination = emission[t] % backward[t];
      arma::mat transition_posterior =
        (forward[t - 1] * destination.t()) % transition[t - 1];
      transition_posterior /= scale[t];
      transition_posterior /= arma::accu(transition_posterior);
      arma::vec origin_posterior = arma::sum(transition_posterior, 1);
      int p = covariates[t].n_cols;

      for (int from_class = 0; from_class < L; ++from_class) {
        for (std::size_t j = 0; j < free_classes.size(); ++j) {
          int to_class = free_classes[j];
          double residual = transition_posterior(from_class, to_class) -
            origin_posterior[from_class] * transition[t - 1](from_class, to_class);
          int block = from_class * (L - 1) + static_cast<int>(j);
          for (int v = 0; v < p; ++v) {
            score[gamma_position[t - 1] + block * p + v] +=
              covariates[t](n, v) * residual;
          }
        }
      }
    }
  }

  if (!std::isfinite(log_likelihood)) stop("Non-finite log-likelihood detected");
  return List::create(
    _["objective"] = -log_likelihood,
    _["gradient"] = compute_gradient ? wrap(-score) : R_NilValue
  );
}
