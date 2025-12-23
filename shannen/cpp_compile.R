# shannen/helpers.R
library(pkgbuild)
library(devtools)
library(Rcpp)

reload_cpp <- function(path=".", compile_attrs=FALSE) {
  path <- normalizePath(path)
  
  if (compile_attrs) Rcpp::compileAttributes(path)
  
  pkgbuild::compile_dll(path, quiet = FALSE)
  devtools::load_all(path, recompile = FALSE, quiet = FALSE)
  
  invisible(TRUE)
}
