#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

namespace {

constexpr double lpa_covariance_jitter = 1e-4;

bool covariance_cholesky(arma::mat& lower, arma::mat covariance,
                         bool /*repair*/) {
  covariance = 0.5 * (covariance + covariance.t());
  covariance.diag() += lpa_covariance_jitter;
  if (arma::chol(lower, covariance, "lower")) return true;
  return false;
}

bool multivariate_log_density(arma::vec& result,
                              const arma::mat& response,
                              const arma::rowvec& mean,
                              const arma::mat& covariance,
                              bool repair) {
  arma::mat lower;
  if (!covariance_cholesky(lower, covariance, repair)) return false;
  arma::mat deviation = response.each_row() - mean;
  arma::mat standardized;
  bool solved = arma::solve(
    standardized, arma::trimatl(lower), deviation.t(),
    arma::solve_opts::no_approx
  );
  if (!solved) {
    arma::mat inverse;
    if (!arma::pinv(inverse, lower)) return false;
    standardized = inverse * deviation.t();
  }
  arma::rowvec quadratic = arma::sum(arma::square(standardized), 0);
  double constant = -0.5 * (
    response.n_cols * std::log(2.0 * arma::datum::pi) +
    2.0 * arma::sum(arma::log(lower.diag()))
  );
  result = constant - 0.5 * quadratic.t();
  return result.is_finite();
}

}

// [[Rcpp::export]]
List lpa_expectation_cpp(const arma::mat& response,
                         const arma::mat& means,
                         const arma::cube& covs,
                         const arma::vec& P_Z,
                         bool repair) {
  arma::uword N = response.n_rows;
  arma::uword I = response.n_cols;
  arma::uword L = means.n_rows;
  if (means.n_cols != I || covs.n_rows != I || covs.n_cols != I ||
      covs.n_slices != L || P_Z.n_elem != L) {
    stop("Incompatible LPA dimensions");
  }

  constexpr double model_eps = 1e-6;
  arma::mat log_joint(N, L);
  for (arma::uword l = 0; l < L; ++l) {
    if (!std::isfinite(P_Z[l]) || P_Z[l] <= 0.0) {
      return List::create(_["valid"] = false);
    }
    arma::vec log_density;
    if (!multivariate_log_density(
          log_density, response, means.row(l), covs.slice(l), repair
        )) {
      return List::create(_["valid"] = false);
    }
    log_joint.col(l) = log_density + std::log(P_Z[l] + model_eps);
  }

  arma::mat posterior(N, L);
  double log_likelihood = 0.0;
  for (arma::uword n = 0; n < N; ++n) {
    double maximum = log_joint.row(n).max();
    arma::rowvec relative = arma::exp(log_joint.row(n) - maximum);
    double denominator = arma::sum(relative);
    if (!std::isfinite(denominator) || denominator <= 0.0) {
      return List::create(_["valid"] = false);
    }
    posterior.row(n) = relative / denominator;
    log_likelihood += maximum + std::log(denominator);
  }

  return List::create(
    _["posterior"] = posterior,
    _["Log.Lik"] = log_likelihood,
    _["valid"] = true
  );
}

// [[Rcpp::export]]
NumericVector mvn_log_density_cpp(const arma::mat& response,
                                  const arma::rowvec& mean,
                                  const arma::mat& covariance,
                                  bool repair) {
  if (mean.n_elem != response.n_cols ||
      covariance.n_rows != response.n_cols ||
      covariance.n_cols != response.n_cols) {
    stop("Incompatible multivariate-normal dimensions");
  }
  arma::vec result;
  if (!multivariate_log_density(result, response, mean, covariance, repair)) {
    stop("Unable to obtain a positive-definite covariance matrix");
  }
  return wrap(result);
}
