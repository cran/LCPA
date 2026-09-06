#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// BCH-weighted multinomial logit kernel.
//
// The parameter vector is ordered by non-reference destination class, with all
// design coefficients for one class stored contiguously. The score returned for
// observation i and class c is
//
//   x_i * {w_ic - w_i+ * p_ic},
//
// and the analytic information block for free classes c and d is
//
//   sum_i w_i+ * p_ic * {I(c=d) - p_id} * x_i x_i'.
//
// BCH weights may be negative. They are deliberately retained because clipping
// changes the BCH estimating equation.
//
// [[Rcpp::export]]
List bch_multinomial_cpp(const arma::vec& params,
                         const arma::mat& design,
                         const arma::mat& weight,
                         int ref_class,
                         bool compute_gradient,
                         bool compute_information) {
  const arma::uword N = design.n_rows;
  const arma::uword p = design.n_cols;
  const arma::uword L = weight.n_cols;
  if (weight.n_rows != N || L < 2 || ref_class < 1 ||
      ref_class > static_cast<int>(L)) {
    stop("Invalid BCH multinomial dimensions or reference class");
  }

  std::vector<arma::uword> free_classes;
  free_classes.reserve(L - 1);
  for (arma::uword l = 0; l < L; ++l) {
    if (l != static_cast<arma::uword>(ref_class - 1)) {
      free_classes.push_back(l);
    }
  }
  const arma::uword q = p * (L - 1);
  if (params.n_elem != q) {
    stop("The BCH multinomial parameter vector has incompatible length");
  }

  arma::mat coefficient(p, L, arma::fill::zeros);
  for (arma::uword j = 0; j < free_classes.size(); ++j) {
    coefficient.col(free_classes[j]) = params.subvec(j * p, (j + 1) * p - 1);
  }

  arma::mat eta = design * coefficient;
  arma::mat probability(N, L);
  double objective = 0.0;
  for (arma::uword i = 0; i < N; ++i) {
    const double maximum = eta.row(i).max();
    arma::rowvec value = arma::exp(eta.row(i) - maximum);
    const double total = arma::accu(value);
    if (!std::isfinite(total) || total <= 0.0) {
      stop("Non-finite BCH multinomial probability");
    }
    probability.row(i) = value / total;
    for (arma::uword l = 0; l < L; ++l) {
      objective -= weight(i, l) *
        std::log(std::max(probability(i, l), 1e-300));
    }
  }
  if (!compute_gradient && !compute_information) {
    return List::create(_["objective"] = objective);
  }

  const arma::vec row_mass = arma::sum(weight, 1);
  arma::mat score(N, q, arma::fill::zeros);
  for (arma::uword j = 0; j < free_classes.size(); ++j) {
    const arma::uword l = free_classes[j];
    score.cols(j * p, (j + 1) * p - 1) = design.each_col() %
      (weight.col(l) - row_mass % probability.col(l));
  }
  arma::vec gradient = -arma::sum(score, 0).t();
  NumericVector gradient_vector(gradient.begin(), gradient.end());
  if (!compute_information) {
    return List::create(
      _["objective"] = objective,
      _["gradient"] = gradient_vector,
      _["score"] = score,
      _["probability"] = probability
    );
  }

  arma::mat information(q, q, arma::fill::zeros);
  for (arma::uword j = 0; j < free_classes.size(); ++j) {
    const arma::uword class_j = free_classes[j];
    for (arma::uword k = 0; k < free_classes.size(); ++k) {
      const arma::uword class_k = free_classes[k];
      arma::vec multiplier = row_mass % probability.col(class_j) %
        ((j == k ? 1.0 : 0.0) - probability.col(class_k));
      information.submat(j * p, k * p, (j + 1) * p - 1,
                         (k + 1) * p - 1) =
        design.t() * (design.each_col() % multiplier);
    }
  }

  return List::create(
    _["objective"] = objective,
    _["gradient"] = gradient_vector,
    _["score"] = score,
    _["probability"] = probability,
    _["information"] = information
  );
}
