// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List adjust_response_cpp(NumericMatrix response) {
  int N = response.nrow();
  int I = response.ncol();

  // 存储每列的唯一非 NA 值
  std::vector<std::vector<double>> poly_temp(I);
  std::vector<int> poly_value(I, 0);

  // 第一步：提取每列非 NA 唯一值并排序
  for (int i = 0; i < I; ++i) {
    std::set<double> unique_vals;
    for (int n = 0; n < N; ++n) {
      double val = response(n, i);
      if (!NumericVector::is_na(val)) {
        unique_vals.insert(val);
      }
    }
    poly_temp[i] = std::vector<double>(unique_vals.begin(), unique_vals.end());
    poly_value[i] = poly_temp[i].size();
  }

  // 计算最大类别数
  int poly_max = *std::max_element(poly_value.begin(), poly_value.end());

  // 构造 poly.orig: I x poly_max 矩阵
  NumericMatrix poly_orig(I, poly_max);
  poly_orig.fill(NA_REAL);
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < poly_value[i]; ++j) {
      poly_orig(i, j) = poly_temp[i][j];
    }
  }

  // 构造 response.adjusted
  NumericMatrix response_adjusted(N, I);
  response_adjusted.fill(NA_REAL);

  // 使用映射加速：每列构建值 → 编码 的 map
  for (int i = 0; i < I; ++i) {
    const auto& levels = poly_temp[i];
    std::map<double, int> value_to_code;
    for (size_t po = 0; po < levels.size(); ++po) {
      value_to_code[levels[po]] = static_cast<int>(po);  // 0-based
    }

    for (int n = 0; n < N; ++n) {
      double val = response(n, i);
      if (!NumericVector::is_na(val)) {
        auto it = value_to_code.find(val);
        if (it != value_to_code.end()) {
          response_adjusted(n, i) = it->second;
        }
      }
    }
  }

  // 恢复原始 NA 位置
  for (int n = 0; n < N; ++n) {
    for (int i = 0; i < I; ++i) {
      if (NumericVector::is_na(response(n, i))) {
        response_adjusted(n, i) = NA_REAL;
      }
    }
  }

  return List::create(
    Named("poly.orig") = poly_orig,
    Named("poly.value") = wrap(poly_value),
    Named("poly.max") = poly_max,
    Named("response") = response_adjusted
  );
}
