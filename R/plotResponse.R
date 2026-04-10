#' @title Visualize Response Distributions with Density Plots
#'
#' @description
#' Creates a publication-quality density plot showing the distribution of responses
#' across multiple indicators/items/variables. Automatically handles variable ordering, color
#' scaling, and legend layout based on the number of variables.
#'
#' @param response A matrix or data frame containing response data where:
#' \itemize{
#'   \item Rows represent respondents, samples or observations
#'   \item Columns represent variables, indicators or questions (must have numeric suffixes, e.g., "indicator1", "Q2")
#'   \item Cell values contain numeric responses
#' }
#' Non-numeric columns (except row identifiers) will cause errors.
#'
#' @return A \code{ggplot} object containing:
#' \itemize{
#'   \item Density curves for each variable colored by indicator
#'   \item Adaptive color schemes based on number of variables
#'   \item Optimized legend layout for large numbers of indicators
#'   \item Publication-ready theme with grid lines and clean styling
#' }
#' The plot can be further customized using standard ggplot2 syntax.
#'
#' @section Theming Details:
#' The plot uses a minimal theme with:
#' \itemize{
#'   \item Light grey grid lines for readability
#'   \item Black axis lines and ticks (0.8pt thickness)
#'   \item White background with no panel border
#'   \item Optimized font sizes (13pt axis titles, 11pt tick labels)
#'   \item Legend positioned on right with adaptive sizing
#' }
#'
#' @examples
#' # Simulate response data for 5 indicators
#' set.seed(42)
#' resp_data <- data.frame(
#'   indicator1 = rnorm(100, mean = 3, sd = 1),
#'   indicator2 = rnorm(100, mean = 2, sd = 0.8),
#'   indicator3 = rnorm(100, mean = 4, sd = 1.2),
#'   indicator4 = rnorm(100, mean = 3.5, sd = 0.9),
#'   indicator5 = rnorm(100, mean = 2.5, sd = 1.1)
#' )
#'
#' library(LCPA)
#' # Generate and display plot
#' p <- plotResponse(resp_data)
#' print(p)
#'
#' # For data with many indicators (18 indicators example)
#' many_indicators <- as.data.frame(replicate(18, rnorm(50, mean = runif(1, 1, 5), sd = 1)))
#' names(many_indicators) <- paste0("Q", 1:18)
#' p_large <- plotResponse(many_indicators)
#' print(p_large)
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr n
#' @importFrom grDevices hcl
#' @export
plotResponse <- function(response) {
  response <- as.data.frame(response)
  response$row_id <- 1:nrow(response)
  response_long <- response %>%
    pivot_longer(cols = -all_of("row_id"), names_to = "variable", values_to = "value") %>%
    mutate(
      var_num = as.numeric(gsub("\\D+", "", .data[["variable"]])),
      variable = factor(.data[["variable"]],
                        levels = unique(.data[["variable"]][order(.data[["var_num"]])]),
                        ordered = TRUE)
    )
  response_long <- response_long %>%
    mutate(
      var_num = as.numeric(gsub("\\D+", "", .data[["variable"]])),
      variable = factor(.data[["variable"]],
                        levels = unique(.data[["variable"]][order(.data[["var_num"]])]),
                        ordered = TRUE)
    )

  n_vars <- ncol(response) - 1

  # Use .data pronoun to avoid NOTEs
  p <- ggplot(response_long, aes(x = .data[["value"]], fill = .data[["variable"]])) +
    geom_density(alpha = ifelse(n_vars > 15, 0.2, 0.4), color = "black", linewidth = 0.5) +
    labs(x = "Value", y = "Density") +
    theme_minimal() +
    theme(
      axis.ticks = element_line(color = "#333333", linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      axis.line = element_line(color = "#333333", linewidth = 0.8),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 1.0),
      axis.text = element_text(size = 11, color = "#333333", margin = margin(t = 2, r = 3, b = 2, l = 3)),
      axis.title = element_text(size = 13, color = "#333333"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right",
      panel.grid.minor = element_line(linewidth = 0.3, color = "grey90"),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey90")
    )

  if (n_vars <= 9) {
    p <- p + scale_fill_brewer(palette = "Set1")
  } else if (n_vars <= 12) {
    p <- p + scale_fill_brewer(palette = "Set3")
  } else {
    message("Number of variables exceeds 12. Using default color wheel.")
    p <- p + scale_fill_hue(h = c(15, 375), l = 80, c = 100)  # Fixed invalid parameters
  }

  if (n_vars > 15) {
    p <- p +
      guides(fill = guide_legend(
        ncol = ceiling(n_vars / 15),
        title.theme = element_text(size = 13),
        label.theme = element_text(size = 11)
      )) +
      theme(
        legend.direction = "vertical",
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.title = element_blank()
      )
  } else {
    p <- p + guides(fill = guide_legend(
      title = NULL,
      title.theme = element_text(size = 13),
      label.theme = element_text(size = 11)
    ))
  }

  return(p)
}
