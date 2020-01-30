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



##------------------------------------------------------------------------
#' isConst
#' 
#' determines whether an expression depends on time 
#' 
#' @export
isConst <- function(exprn)
{
  if (is.atomic(exprn) || is.name(exprn))
  {
    first <- exprn
  } else
  {
    first <- exprn[[1]]
  }
  
  if (is.numeric(first) && !is.array(first))
  {
    return(T)
  }
  if (is.name(first) && deparse(first) == "t"){
    return(F)
  }
  if (deparse(first) %in% c("+", "-", "/", "*", "^", "<", ">", "<=", ">="))
  {
    # allowed operations
    if (length(exprn) > 2)
    {
      r1 = isConst(exprn[[2]])
      r2 = isConst(exprn[[3]])
      return(r1 && r2)
    } else
    {
      return(isConst(exprn[[2]]))
    }
  }
  if(deparse(first) %in% c("exp", "log", "sin", "cos", "("))
  {    
    if (length(exprn) != 2)
    {
      stop("unary function applied to wrong number of parameters")  
    } 
    return(isConst(exprn[[2]]))
  }
  if (deparse(first) == "[")
  {
    if (!is.na(stringr::str_extract(deparse(exprn), "^params\\[[0-9]+\\]$")))
    {
      s <- stringr::str_extract(deparse(exprn), "^params\\[[0-9]+\\]$")
      num <- strtoi(substr(s, 8, nchar(s) - 1))
      if (num <= 0)
      {
        stop(paste("parameter indices must be positive!", 
                   sep = " "))
      }
      return(T)
    }
    stop(paste("Invalid expression:", deparse(exprn), sep = " "))
  } 
  stop(paste("Invalid expression:", deparse(exprn), sep = " "))
}

##------------------------------------------------------------------------
#' generateCpp
#' 
#' generates custom C++ code for rate computation
#' 
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
  
  if (is.numeric(first) && !is.array(first))
  {
    return(first)
  }
  if (is.name(first) && deparse(first) %in% c("t")){
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
    if (length(exprn) != 2)
    {
      stop("unary function applied to wrong number of parameters")  
    } 
    else
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

##------------------------------------------------------------------------
#' formatSimData
#' 
#' coerces data from simulation fuction into correct format for estimation function
#' 
#' @export
formatSimData = function(sim_data, ntypes){
  for (i in 1:ntypes)
  {
    cellname <- sprintf("t%d_cells", i)
    names(sim_data)[2 +  i] <- cellname
    prevname <- paste(cellname, "prev", sep = "_")
    sim_data <- dplyr::group_by(sim_data, rep)
    sim_data <- dplyr::mutate(sim_data, prev = dplyr::lag(!!dplyr::sym(cellname)))
    names(sim_data)[2 + ntypes + i] <- prevname
  }
  sim_data <- dplyr::mutate(sim_data, prev_time = dplyr::lag(time))
  sim_data <- dplyr::filter(sim_data,  !is.na(t1_cells_prev))
  sim_data
}