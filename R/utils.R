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
#' @export
walk_ast <- function(ast, base_fn, combine_fn) {
  if (is.atomic(ast) || is.name(ast)) {
    return(base_fn(ast))
  } else if (is.call(ast)) {
    return(combine_fn(ast[[1]], lapply(ast[-1], walk_ast, base_fn, combine_fn), ast[-1]))
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
    return(F)
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec, args){
    if (deparse(fname) %in% c("+", "-", "/", "*", "^", "<", ">", "<=", ">=") && length(rec) == 2)
    {
      return(rec[[1]] && rec[[2]])
    }
    if(deparse(fname) %in% c("+", "-", "exp", "log", "sin", "cos", "(") && length(rec) == 1)
    {    
      return(rec[[1]])
    }
    if (deparse(fname) == "[" && deparse(args[[1]]) == "params" && is.numeric(args[[2]]) > 0 && args[[2]] > 0)
    {
        return(T)
    }
    return(F)
  }
  walk_ast(ast, base_fn, combine_fn)
}




##------------------------------------------------------------------------
#' is_const
#'  
#' helper for determining whether an expression depends on time
#' 
#' @export
is_const <- function(ast) {
  base_fn <- function(x){
    if(deparse(x) == "t"){
      return(F)
    }
    return(T)
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec, args){
    res = T
    for(i in 1:length(args)){res <- rec && args[[i]]}
    return(res)
  }
  walk_ast(ast, base_fn, combine_fn)
}


##------------------------------------------------------------------------
#' generate_cpp
#'  
#' generate_cpp C++ code from R expression
#' 
#' @export
generate_cpp <- function(ast, params) {
  base_fn <- function(x){
    return(deparse(x))
  }
  
  # fname = name of function being applies, rec = resluts of recusively computing function on arguments, args = values of arguments
  combine_fn <- function(fname, rec, args){
    if (deparse(fname) %in% c("+", "-", "/", "*", "^", "<", ">", "<=", ">=") && length(args) == 2)
    {
      return(paste("(", rec[[1]], ")", deparse(fname), "(",rec[[2]], ")", sep = " "))    
    }
    if(deparse(fname) %in% c("+", "-", "exp", "log", "sin", "cos", "(") && length(rec) == 1)
    {    
      return(paste(deparse(fname), "(", rec[[1]], ")", sep = ""))
    }
    if (deparse(fname) == "[")
    {
      idx = as.numeric(rec[[2]])
      if(idx > length(params) || idx <= 0){
        stop(sprintf("Parameter params[%d] goes beyond the number of parameters provided!", idx))
      }
      return(sprintf("%f", params[idx]))
    }
    stop("tried to generate from invalid expression")
  }
  walk_ast(ast, base_fn, combine_fn)
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
    cellname <- sprintf("t%d_cells", i)
    names(sim_data)[2 +  i] <- cellname
    prevname <- paste(cellname, "prev", sep = "_")
    sim_data <- dplyr::group_by(sim_data, rep)
    sim_data <- dplyr::mutate(sim_data, prev = dplyr::lag(!!dplyr::sym(cellname)))
    names(sim_data)[2 + ntypes + i] <- prevname
  }
  sim_data <- dplyr::mutate(sim_data, prev_time = dplyr::lag(time))
  sim_data <- dplyr::filter(sim_data,  !is.na(t1_cells_prev))
  return(sim_data)
}