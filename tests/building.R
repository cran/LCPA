#' script_path <- rstudioapi::getActiveDocumentContext()$path
#' print(script_path)
#' working_directory <- dirname(dirname(script_path))
#' setwd(working_directory)
#'
#' pkgdown::build_site()
#'
#' devtools::document()
#'
#' devtools::build_manual()
#'
#' devtools::build()
#' devtools::install()
#' devtools::check(manual = TRUE)
#'
#' a <- devtools::spell_check(vignettes = TRUE, use_wordlist = TRUE)
#' a
#' a[1]
#'
#' devtools::check_rhub()
#'
#' devtools::check_win_devel()
#'
#' devtools::release()
#' devtools::submit_cran()
#'
#'
#' NEWS.md
#' DESCRIPTION
#' cran-comments.md
#'
#'
#' pack <- available.packages()
#'
#' which(rownames(pack) == "Qval")
#'
#'
#' library(cranlogs)
#' downloads <- cran_downloads(packages = "LCAP", from = "2024-08-01", to = "2026-4-09")
#' print(downloads)
#' sum(downloads$count)
#'
#'
