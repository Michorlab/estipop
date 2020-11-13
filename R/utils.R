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
#' walk_exprn
#'  
#' helper for walking ast in R
#' 
walk_ast <- function(ast, base_fn, combine_fn, is_base_case) {
  if (is.atomic(ast) || is.name(ast) || is_base_case(ast)) {
    return(base_fn(ast))
  } else if (is.call(ast)) {
    return(combine_fn(ast[[1]], lapply(ast[-1], walk_ast, base_fn, combine_fn, is_base_case))) #recursively compute on arguments and return results
  } else {
    # User supplied incorrect input
    stop("invalid expression!")
  }
}

##------------------------------------------------------------------------
#' check_valid
#'  
#' helper for determining whether an expression is valid
#' 
#' @export
check_valid <- function(ast) {
  base_fn <- function(x){
    if ((is.numeric(x) && !is.array(x)) || (is.name(x) && deparse(x) == "t"))
    {
      return(T)
    }
    else if (is.name(x)){
      stop(sprintf("invalid name %s", deparse(x)))
    }
    if (is.call(x) && deparse(x[[1]]) == "[" && deparse(x[[2]]) == "params" && is.numeric(x[[3]]) && x[[2]] > 0)
    {
      return(T)
    }
    stop(sprintf("invalid expression %s", deparse(x)))
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec){
    if (deparse(fname) %in% c("/", "*", "^", "<", ">", "<=", ">="))
    {
      if(length(rec) == 2){
        return(rec[[1]] && rec[[2]])
      }
      stop(sprintf("binary operator %s takes 2 argments!", deparse(fname)))
    }
    if(deparse(fname) %in% c("exp", "log", "log10", "sin", "cos", "tan", "asin", "acos", "atan", "tanh", "cosh", "sinh", "atanh", "acosh", "asinh", "("))
    {    
      if(length(rec) == 1){
        return(rec[[1]])
      }
      stop(sprintf("unary operator %s takes 1 argment!", deparse(fname)))
    }
    if(deparse(fname) %in% c("+","-")){
      if(length(rec) == 1){
        return(rec[[1]])
      }
      if(length(rec) == 2){
        return(rec[[1]] && rec[[2]])
      }
      stop("operator %s takes 1 or 2 arguments!", deparse(fname))
    }
    stop(sprintf("invalid operator %s!", deparse(fname)))
  }
  
  is_base_case <- function (ast){
    return(is.call(ast) && ast[[1]] == "[")
  }
  walk_ast(ast, base_fn, combine_fn, is_base_case)
}




##------------------------------------------------------------------------
#' is_const
#'  
#' helper for determining whether an expression depends on time
#' 
#' @export
is_const <- function(ast) {
  check_valid(ast)
  base_fn <- function(x){
    if(deparse(x) == "t"){
      return(F)
    }
    return(T)
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec){
    res = T
    for(i in 1:length(rec)){res <- res && rec[[i]]}
    return(res)
  }
  
  is_base_case <- function(ast){
    return(is.call(ast) && ast[[1]] == "[")
  }
  walk_ast(ast, base_fn, combine_fn, is_base_case)
}


##------------------------------------------------------------------------
#' generate_cpp
#'  
#' generate_cpp C++ code from R expression
#' 
#' @export
generate_cpp <- function(ast, params) {
  check_valid(ast)
  base_fn <- function(x){
    if (is.call(x) && deparse(x[[1]]) == "[")
    {
      idx = as.numeric(x[[3]])
      if(idx > length(params) || idx <= 0){
        stop(sprintf("Parameter params[%d] goes beyond the number of parameters provided!", idx))
      }
      return(sprintf("%f", params[idx]))
    }
    return(deparse(x))
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec){
    if (deparse(fname) %in% c("+", "-", "/", "*", "^", "<", ">", "<=", ">=") && length(rec) == 2)
    {
      return(paste(rec[[1]], deparse(fname), rec[[2]],  sep = " "))    
    }
    if(deparse(fname) %in% c("+", "-", "exp", "log", "sin", "cos", "(") && length(rec) == 1)
    {    
      if(deparse(fname) != "("){
        return(paste(deparse(fname), "(", rec[[1]], ")", sep = ""))
      }
      return(paste("(", rec[[1]], ")", sep=""))
    }
  }
  
  is_base_case <- function(ast){
    return(is.call(ast) && ast[[1]] == "[")
  }
  walk_ast(ast, base_fn, combine_fn, is_base_case)
}



##------------------------------------------------------------------------
#' formatSimData
#' 
#' coerces data from simulation fuction into correct format for estimation function
#' 
#' @export
format_sim_data = function(sim_data, ntypes){
  for (i in 1:ntypes)
  {
    cellname <- sprintf("type%d", i)
    prevname <- paste(cellname, "prev", sep = "_")
    sim_data <- dplyr::group_by(sim_data, rep)
    sim_data <- dplyr::mutate(sim_data, prev = dplyr::lag(!!dplyr::sym(cellname)))
    names(sim_data)[2 + ntypes + i] <- prevname
  }
  sim_data <- dplyr::mutate(sim_data, prev_time = dplyr::lag(time))
  sim_data <- dplyr::mutate(sim_data, dtime = time - prev_time)
  if(all(is.na(sim_data$type1_prev))){
    warning("There are no initial data points in this dataset! Make sure to record at least 2 timepoints per simulation")
  }
  sim_data <- dplyr::filter(sim_data,  !is.na(type1_prev))
  return(sim_data)
}