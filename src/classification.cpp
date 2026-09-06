#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List classification_error_counts_cpp(const NumericMatrix& posterior) {
  int N = posterior.nrow();
  int L = posterior.ncol();
  NumericMatrix count(L, L);
  NumericVector class_mass(L);

  for (int n = 0; n < N; ++n) {
    int predicted = 0;
    double maximum = posterior(n, 0);
    for (int l = 1; l < L; ++l) {
      if (posterior(n, l) > maximum) {
        maximum = posterior(n, l);
        predicted = l;
      }
    }
    for (int true_class = 0; true_class < L; ++true_class) {
      double value = posterior(n, true_class);
      class_mass[true_class] += value;
      count(true_class, predicted) += value;
    }
  }

  return List::create(_["count"] = count, _["class.mass"] = class_mass);
}
