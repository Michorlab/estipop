###################################################################################----
#'  Utility functions for time-dependent rates
#'  taken with minor adaptions from:
#'
#'  McDonald, T. O., & Michor, F. (2017).
#'  SIApopr: a computational method to simulate evolutionary branching trees for
#'  analysis of tumor clonal evolution. Bioinformatics, 33(14), 2221-2223.
#'
#'  Original code author: T. O. McDonald (mcdonald(at)jimmy.harvard.edu)
#'

##------------------------------------------------------------------------
#' reload
#'
#' Recompiles and reinstalls the package for use after adjusting the model
#' code. WARNING: For advanced users. Adjusting the C++ code can lead to
#' errors in compilation.
#'
#' @param path path to siapopr package source
#'
#' @export
reload <- function( path ){

  detach("package:estipop", unload = TRUE)
  library.dynam.unload("estipop", system.file(package = "estipop"))

  path <- paste( "--vanilla  CMD INSTALL ", path )

  system2( 'R', path  )
  require("estipop")
}


#' .pop_off
#'
#' pops off last instance of pattern in a string
#'
#' @param string - character string where matches are sought.
#' @param pattern - character string to search for.
#'
#' @return .pop_off - string with everything after final pattern removed.
#' @export
.pop_off <- function(string, pattern = ">", ...){
  string <- unlist(strsplit(string, pattern, ...))
  string <- string[1:length(string)-1]
  paste(string, collapse = pattern)
}

##------------------------------------------------------------------------
#' .pop
#'
#' pops off last element of a string after final pattern in a vector
#'
#' @param string - character string where matches are sought.
#' @param pattern - character string to search for.
#'
#' @return .pop - string of everything after final pattern.
#' @export
.pop <- function(string, pattern = ">", ...){
  string <- unlist(strsplit(string, pattern, ...))
  string[length(string)]
}


#'
#' ##------------------------------------------------------------------------
#' #' .replace
#' #'
#' #' replaces an element of a vector with a particular pattern with the
#' #' replacement given
#' #'
#' #' @param vec - character vector
#' #' @param pattern - character element to search for
#' #' @param replacement - character to replace pattern with
#' #'
#' #' @return .replace - returns a vector with the replaced elements
#' #' @export
#' .replace <- function(vec, pattern, replacement) {
#'   vec[vec == pattern] <- replacement
#'   vec
#' }


make_function <- function(args, body, env = parent.frame()) {
args <- as.pairlist(args)

eval(call("function", args, body), env)
}
