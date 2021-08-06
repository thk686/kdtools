#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname){
  if (!has_cxx17())
    packageStartupMessage("The kdtools package was compiled without c++17 and will have reduced functionality\n")
}
