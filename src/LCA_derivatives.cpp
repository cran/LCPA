#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List lca_score_information_cpp(IntegerMatrix response, NumericVector par,
                               NumericVector P_Z, IntegerVector poly_value,
                               bool compute_information) {
  const int N = response.nrow();
  const int I = response.ncol();
  const int L = P_Z.size();
  int q = L - 1;
  for (int i = 0; i < I; ++i) {
    q += L * (poly_value[i] - 1);
  }

  IntegerMatrix block_start(L, I);
  int position = L - 1;
  for (int i = 0; i < I; ++i) {
    for (int l = 0; l < L; ++l) {
      block_start(l, i) = position;
      position += poly_value[i] - 1;
    }
  }

  arma::mat posterior(N, L, arma::fill::zeros);
  arma::mat score_observation(N, q, arma::fill::zeros);
  arma::vec class_mass(L, arma::fill::zeros);
  arma::mat missing_information;
  if (compute_information) {
    missing_information.zeros(q, q);
  }

  for (int n = 0; n < N; ++n) {
    arma::vec log_joint(L);
    for (int l = 0; l < L; ++l) {
      double value = std::log(P_Z[l]);
      for (int i = 0; i < I; ++i) {
        const int category = response(n, i);
        const int index = l + i * L + category * L * I;
        value += std::log(par[index]);
      }
      log_joint[l] = value;
    }
    const double maximum = log_joint.max();
    arma::vec tau = arma::exp(log_joint - maximum);
    tau /= arma::accu(tau);
    posterior.row(n) = tau.t();
    class_mass += tau;

    if (compute_information) {
      arma::mat complete_score(q, L, arma::fill::zeros);
      for (int l = 0; l < L; ++l) {
        for (int j = 0; j < L - 1; ++j) {
          complete_score(j, l) = (j == l ? 1.0 : 0.0) - P_Z[j];
        }
        for (int i = 0; i < I; ++i) {
          const int category = response(n, i);
          const int start = block_start(l, i);
          for (int k = 0; k < poly_value[i] - 1; ++k) {
            const int par_index = l + i * L + k * L * I;
            complete_score(start + k, l) =
              (category == k ? 1.0 : 0.0) - par[par_index];
          }
        }
      }
      arma::vec score_n = complete_score * tau;
      score_observation.row(n) = score_n.t();
      arma::mat weighted_score = complete_score;
      weighted_score.each_row() %= arma::sqrt(tau).t();
      missing_information += weighted_score * weighted_score.t() -
        score_n * score_n.t();
    } else {
      for (int j = 0; j < L - 1; ++j) {
        score_observation(n, j) = tau[j] - P_Z[j];
      }
      for (int i = 0; i < I; ++i) {
        const int category = response(n, i);
        for (int l = 0; l < L; ++l) {
          const int start = block_start(l, i);
          for (int k = 0; k < poly_value[i] - 1; ++k) {
            const int par_index = l + i * L + k * L * I;
            score_observation(n, start + k) = tau[l] *
              ((category == k ? 1.0 : 0.0) - par[par_index]);
          }
        }
      }
    }
  }

  List result = List::create(
    _["posterior"] = posterior,
    _["score.observation"] = score_observation,
    _["score"] = arma::sum(score_observation, 0).t()
  );

  if (compute_information) {
    arma::mat expected_information(q, q, arma::fill::zeros);
    if (L > 1) {
      arma::vec probability(L - 1);
      for (int j = 0; j < L - 1; ++j) {
        probability[j] = P_Z[j];
      }
      expected_information.submat(0, 0, L - 2, L - 2) =
        N * (arma::diagmat(probability) - probability * probability.t());
    }
    for (int i = 0; i < I; ++i) {
      const int categories = poly_value[i] - 1;
      if (categories == 0) {
        continue;
      }
      for (int l = 0; l < L; ++l) {
        arma::vec probability(categories);
        for (int k = 0; k < categories; ++k) {
          const int par_index = l + i * L + k * L * I;
          probability[k] = par[par_index];
        }
        const int start = block_start(l, i);
        expected_information.submat(
          start, start, start + categories - 1, start + categories - 1
        ) = class_mass[l] *
          (arma::diagmat(probability) - probability * probability.t());
      }
    }
    arma::mat information = expected_information - missing_information;
    result["information"] = (information + information.t()) / 2.0;
  }

  return result;
}
