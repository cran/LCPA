#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {

int draw_logit(const std::vector<double>& eta) {
  double maximum = -std::numeric_limits<double>::infinity();
  for (double value : eta) maximum = std::max(maximum, value);

  double total = 0.0;
  std::vector<double> probability(eta.size());
  for (std::size_t k = 0; k < eta.size(); ++k) {
    probability[k] = std::exp(eta[k] - maximum);
    total += probability[k];
  }

  double draw = R::runif(0.0, total);
  double cumulative = 0.0;
  for (std::size_t k = 0; k < probability.size(); ++k) {
    cumulative += probability[k];
    if (draw <= cumulative) return static_cast<int>(k) + 1;
  }
  return static_cast<int>(probability.size());
}

}

// [[Rcpp::export]]
IntegerVector sample_multinomial_logit_cpp(const NumericMatrix& covariates,
                                            const NumericMatrix& coefficients) {
  int N = covariates.nrow();
  int p = covariates.ncol();
  int L = coefficients.ncol();
  if (coefficients.nrow() != p) stop("Incompatible multinomial-logit dimensions");

  IntegerVector classes(N);
  std::vector<double> eta(L);
  for (int n = 0; n < N; ++n) {
    for (int l = 0; l < L; ++l) {
      eta[l] = 0.0;
      for (int j = 0; j < p; ++j) eta[l] += covariates(n, j) * coefficients(j, l);
    }
    classes[n] = draw_logit(eta);
  }
  return classes;
}

// [[Rcpp::export]]
NumericVector mean_multinomial_logit_probability_cpp(
    const NumericMatrix& covariates,
    const NumericMatrix& coefficients) {
  int N = covariates.nrow();
  int p = covariates.ncol();
  int L = coefficients.ncol();
  if (N < 1 || coefficients.nrow() != p) {
    stop("Incompatible multinomial-logit dimensions");
  }

  NumericVector mean_probability(L);
  std::vector<double> eta(L);
  std::vector<double> probability(L);
  for (int n = 0; n < N; ++n) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (int l = 0; l < L; ++l) {
      eta[l] = 0.0;
      for (int j = 0; j < p; ++j) {
        eta[l] += covariates(n, j) * coefficients(j, l);
      }
      maximum = std::max(maximum, eta[l]);
    }

    double total = 0.0;
    for (int l = 0; l < L; ++l) {
      probability[l] = std::exp(eta[l] - maximum);
      total += probability[l];
    }
    if (!std::isfinite(total) || total <= 0.0) {
      stop("Invalid multinomial-logit probabilities");
    }
    for (int l = 0; l < L; ++l) {
      mean_probability[l] += probability[l] / total;
    }
  }
  return mean_probability / N;
}

// [[Rcpp::export]]
IntegerVector sample_transition_logit_cpp(const IntegerVector& previous,
                                           const NumericMatrix& covariates,
                                           const List& gamma) {
  int N = covariates.nrow();
  int p = covariates.ncol();
  int L = gamma.size();
  if (previous.size() != N) stop("Incompatible transition dimensions");

  std::vector<double> coefficient(static_cast<std::size_t>(L) * L * p);
  for (int from = 0; from < L; ++from) {
    List gamma_from = gamma[from];
    if (gamma_from.size() != L) stop("gamma must contain L destination classes");
    for (int to = 0; to < L; ++to) {
      NumericVector current = gamma_from[to];
      if (current.size() != p) stop("gamma coefficient length does not match covariates");
      for (int j = 0; j < p; ++j) {
        coefficient[(from * L + to) * p + j] = current[j];
      }
    }
  }

  IntegerVector classes(N);
  std::vector<double> eta(L);
  for (int n = 0; n < N; ++n) {
    int from = previous[n] - 1;
    if (from < 0 || from >= L) stop("previous contains an invalid class");
    for (int to = 0; to < L; ++to) {
      eta[to] = 0.0;
      for (int j = 0; j < p; ++j) {
        eta[to] += covariates(n, j) * coefficient[(from * L + to) * p + j];
      }
    }
    classes[n] = draw_logit(eta);
  }
  return classes;
}

// [[Rcpp::export]]
IntegerVector sample_markov_cpp(const IntegerVector& previous,
                                const NumericMatrix& rate) {
  int N = previous.size();
  int L = rate.nrow();
  if (rate.ncol() != L) stop("rate must be square");

  IntegerVector classes(N);
  for (int n = 0; n < N; ++n) {
    int from = previous[n] - 1;
    if (from < 0 || from >= L) stop("previous contains an invalid class");
    double draw = R::runif(0.0, 1.0);
    double cumulative = 0.0;
    classes[n] = L;
    for (int to = 0; to < L; ++to) {
      cumulative += rate(from, to);
      if (draw <= cumulative) {
        classes[n] = to + 1;
        break;
      }
    }
  }
  return classes;
}
