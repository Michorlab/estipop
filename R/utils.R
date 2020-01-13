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
  string <- string[1:lengthgth(string)-1]
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
  string[lengthgth(string)]
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


#' Generate CPP code for computing a custom rate
#' @export
generateCpp <- function(exprn, params)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.atomic(first) && is.numeric(first) && !is.array(first))
  {
    return(first)
  }
  if (is.name(first) && deparse(first) == "t"){
    return("t")
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^", "<", ">", "<=", ">="))
  {
    # allowed operations
    if (length(exprn) > 2)
    {
      return(paste("(", generateCpp(exprn[[2]], params), 
                   ")", deparse(first), "(", generateCpp(exprn[[3]], params), ")", sep = " "))
    } else
    {
      return(paste(deparse(first), "(", generateCpp(exprn[[2]], params), ")", sep = ""))
    }
  }
  if(deparse(first) %in% c("exp", "log", "sin", "cos"))
  {
    if (length(exprn) > 2)
    {
      stop("log and exp take only one parameter")  
    } else
    {
      return(paste(deparse(first), "(", generateCpp(exprn[[2]], params), ")", sep = ""))
    }
    
  }
  if (deparse(first) == "(")
  {
    return(generateCpp(exprn[[2]], params))
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^params\\[[0-9]+\\]$")))
    {
      # function parameters
      s <- stringr::str_extract(deparse(exprn), "^params\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 8, nchar(s) - 1))
      if (num > length(params) || num <= 0)
      {
        print(num)
        print(params)
        stop(paste("Parameter", s, "goes beyond the number of parameters specified.", 
                   sep = " "))
      }
      return(sprintf("%f", params[num]))
    }
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  } 
  else
  {
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  }
}