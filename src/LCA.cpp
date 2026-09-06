#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List lca_expectation_cpp(const IntegerMatrix& response,
                         const NumericVector& par,
                         const NumericVector& P_Z) {
  IntegerVector dimensions = par.attr("dim");
  if (dimensions.size() != 3) stop("par must be a three-dimensional array");

  int N = response.nrow();
  int I = response.ncol();
  int L = dimensions[0];
  int categories = dimensions[2];
  if (dimensions[1] != I || P_Z.size() != L) stop("Incompatible LCA dimensions");

  NumericMatrix posterior(N, L);
  NumericVector log_likelihood_observation(N);
  bool valid = true;

  for (int n = 0; n < N; ++n) {
    std::vector<double> log_joint(L);
    double maximum = -std::numeric_limits<double>::infinity();

    for (int l = 0; l < L; ++l) {
      double value = P_Z[l] > 0.0
        ? std::log(P_Z[l])
        : -std::numeric_limits<double>::infinity();
      for (int i = 0; i < I; ++i) {
        int category = response(n, i);
        if (category < 0 || category >= categories) stop("response category is outside par");
        double probability = par[l + i * L + category * L * I];
        if (!std::isfinite(probability) || probability < 0.0) {
          stop("par contains an invalid probability for an observed category");
        }
        if (probability == 0.0) {
          value = -std::numeric_limits<double>::infinity();
          break;
        }
        value += std::log(probability);
      }
      log_joint[l] = value;
      maximum = std::max(maximum, value);
    }

    if (!std::isfinite(maximum)) {
      valid = false;
      for (int l = 0; l < L; ++l) posterior(n, l) = 1.0 / L;
      log_likelihood_observation[n] = std::log(std::numeric_limits<double>::min());
      continue;
    }

    double denominator = 0.0;
    for (int l = 0; l < L; ++l) {
      posterior(n, l) = std::exp(log_joint[l] - maximum);
      denominator += posterior(n, l);
    }
    for (int l = 0; l < L; ++l) posterior(n, l) /= denominator;
    log_likelihood_observation[n] = maximum + std::log(denominator);
  }

  double log_likelihood = Rcpp::sum(log_likelihood_observation);
  return List::create(
    _["posterior"] = posterior,
    _["Log.Lik"] = log_likelihood,
    _["log.likelihood.observation"] = log_likelihood_observation,
    _["valid"] = valid
  );
}

// [[Rcpp::export]]
List lca_maximization_cpp(const IntegerMatrix& response,
                          const NumericMatrix& posterior,
                          const IntegerVector& poly_value,
                          double smoothing) {
  int N = response.nrow();
  int I = response.ncol();
  int L = posterior.ncol();
  int poly_max = max(poly_value);
  if (posterior.nrow() != N || poly_value.size() != I) stop("Incompatible LCA dimensions");

  NumericVector par(L * I * poly_max, NA_REAL);
  par.attr("dim") = Dimension(L, I, poly_max);
  NumericVector P_Z(L);
  if (!std::isfinite(smoothing) || smoothing <= 0.0) stop("smoothing must be positive");

  for (int l = 0; l < L; ++l) {
    double class_mass = 0.0;
    for (int n = 0; n < N; ++n) class_mass += posterior(n, l);
    P_Z[l] = class_mass / N;

    for (int i = 0; i < I; ++i) {
      for (int category = 0; category < poly_value[i]; ++category) {
        double count = 0.0;
        for (int n = 0; n < N; ++n) {
          if (response(n, i) == category) count += posterior(n, l);
        }
        par[l + i * L + category * L * I] =
          (count + smoothing) / (class_mass + poly_value[i] * smoothing);
      }
    }
  }

  return List::create(_["par"] = par, _["P.Z"] = P_Z);
}

// [[Rcpp::export]]
IntegerMatrix sample_lca_response_cpp(const IntegerVector& classes,
                                      const NumericVector& par,
                                      const IntegerVector& poly_value) {
  IntegerVector dimensions = par.attr("dim");
  if (dimensions.size() != 3) stop("par must be a three-dimensional array");
  int N = classes.size();
  int L = dimensions[0];
  int I = dimensions[1];
  if (poly_value.size() != I) stop("poly_value must contain one value per indicator");

  IntegerMatrix response(N, I);
  for (int n = 0; n < N; ++n) {
    int l = classes[n] - 1;
    if (l < 0 || l >= L) stop("classes contains an invalid class");
    for (int i = 0; i < I; ++i) {
      double draw = R::runif(0.0, 1.0);
      double cumulative = 0.0;
      response(n, i) = poly_value[i] - 1;
      for (int category = 0; category < poly_value[i]; ++category) {
        cumulative += par[l + i * L + category * L * I];
        if (draw <= cumulative) {
          response(n, i) = category;
          break;
        }
      }
    }
  }
  return response;
}
