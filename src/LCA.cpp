#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List em_e_step(IntegerMatrix Y, NumericVector par_pre, NumericVector P_Z_pre, int N, int I, int L, int poly_max) {
  NumericMatrix L_Xi_Z(L, N);
  NumericMatrix P_Z_Xn(N, L);
  NumericVector L_Xi_vec(N);

  // Initialize L_Xi_Z to 1
  for (int idx = 0; idx < L*N; idx++) {
    L_Xi_Z[idx] = 1.0;
  }

  // E-step calculations
  for (int p = 0; p < N; p++) {
    for (int i = 0; i < I; i++) {
      int y_val = Y(p, i);
      for (int l = 0; l < L; l++) {
        int idx_par = l + i*L + y_val*L*I;
        L_Xi_Z(l, p) *= par_pre[idx_par];
      }
    }

    double L_Xi = 0.0;
    for (int l = 0; l < L; l++) {
      L_Xi += L_Xi_Z(l, p) * P_Z_pre[l];
    }

    if (L_Xi < 1e-300) L_Xi = 1e-300;
    L_Xi_vec[p] = L_Xi;

    for (int l = 0; l < L; l++) {
      P_Z_Xn(p, l) = (L_Xi_Z(l, p) * P_Z_pre[l]) / L_Xi;
    }
  }

  return List::create(
    _["L_Xi_Z"] = L_Xi_Z,
    _["P_Z_Xn"] = P_Z_Xn,
    _["L_Xi_vec"] = L_Xi_vec
  );
}

// [[Rcpp::export]]
List em_m_step(NumericMatrix P_Z_Xn, NumericVector Y_hot, IntegerVector poly_value, int N, int I, int L, int poly_max) {
  NumericVector Rlj_vec(L * I * poly_max, 0.0);
  NumericVector Ilj_vec(L * I * poly_max, 0.0);
  NumericVector S_l(L, 0.0);

  // Precompute S_l = sum_p P_Z_Xn(p, l)
  for (int l = 0; l < L; l++) {
    for (int p = 0; p < N; p++) {
      S_l[l] += P_Z_Xn(p, l);
    }
  }

  // M-step calculations
  for (int l = 0; l < L; l++) {
    for (int i = 0; i < I; i++) {
      int max_po = poly_value[i];
      for (int po = 0; po < max_po; po++) {
        double sum_val = 0.0;
        for (int p = 0; p < N; p++) {
          int idx_yhot = p + i*N + po*N*I;
          sum_val += P_Z_Xn(p, l) * Y_hot[idx_yhot];
        }

        int idx_rlj = l + i*L + po*L*I;
        Rlj_vec[idx_rlj] = sum_val;
        Ilj_vec[idx_rlj] = S_l[l];
      }
    }
  }

  // Set dimensions
  Dimension dim(L, I, poly_max);
  Rlj_vec.attr("dim") = dim;
  Ilj_vec.attr("dim") = dim;

  return List::create(
    _["Rlj"] = Rlj_vec,
    _["Ilj"] = Ilj_vec
  );
}
