.onAttach <- function(libname, pkgname) {
  tv <- utils::packageVersion("sdStaf")
  m <- paste0("This is version ", tv, " of the \"sdStaf\" package, for testing only\n")
  packageStartupMessage(m)
}