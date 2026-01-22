
.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

StartWelcomeMessage <- function(){
  paste("LCPA R Package ",
        "(version ", utils::packageDescription("LCPA")$Version,
        "; ", utils::packageDescription("LCPA")$Date, ")\n",
        sep="")
}

printPackageInfo <- function() {
  packageinfo <- utils::packageDescription("LCPA")
  cat(paste("LCPA version ", packageinfo$Version, " (", packageinfo$Date, ")", sep = ""), "\n")
}
