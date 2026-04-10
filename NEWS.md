
# LCPA 1.0.2

-   Added -  Added a new function `adjust.model` for aligning two `LCA` or `LPA` objects.
-   Added -  Added visualization plots (S3 method `plot`) for `LCA` and `LPA`.
-   Fixed -  Resolved the case sensitivity issue of observed indicator names when calling Mplus
-   Fixed -  Improved Python environment selection in the `install_python_dependencies` function.
-   Fixed -  Fixed the description issues in the comment documentation.
-   Fixed -  Enhanced the stability of the `LRT.test.Bootstrap` function.
-   Fixed -  Removed unnecessary Mplus statements from `LCA` and `LPA`.
-   Fixed -  Improved the efficiency of Python code calling NVIDIA GPUs (Cuda + cuDNN + Pytorch) for deep network training.

# LCPA 1.0.1

-   Fixed - Resolved the configuration issue with `control.NNE`.
-   Fiexd - Resolved return parameter `res$params$par` naming issue in LCA.
-   Added - For `LCA`, `LPA`, `LCPA`, and `LTA`, users can control whether to enable the self-attention mechanism (i.e., transformer encoder) through the `use.attention` arguments.

# LCPA 1.0.0

-   Initial release
