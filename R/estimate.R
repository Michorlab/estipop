#' estimate_td
#'
#' Estimates the rate parameters for a general multitype branching process with time-dependent rates
#'
#' @param model the \code{process_model} object representing the process generating the data
#' @param params the vector of parameters for which we are computing the likelihood
#' @param init_pop a \code{nobs x ntype} matrix with initial population for each observation
#' @param start_times the \code{nobs} length vector of times at which the initial populations were observed
#' @param end_times the \code{nobs} length vector of times at which the final populations were observed
#' @param final_pop the \code{nobs x mtype} matrix of final populations observed
#' @param initial vector of initial parameters estimates for MLE optimization
#' @param known boolean vector of known parameter rates, if NULL, all rates will be estimated
#' @param lower vector of lower bounds on rate parameters for optimization
#' @param upper vector of upper bounds on rate parameters for optimization
#' @param trace level of output for optimizer - see optim function, control - trace for "L-BFGS-B" method
#'
#' @export
estimate_td = function(model, init_pop, start_times, end_times, final_pop, initial, known = NULL, lower = NULL, upper = NULL, trace = 1){
  if(class(model) != "estipop_process_model"){
    stop("model must be a process_model object!")
  }
  if (is.vector(final_pop) | is.vector(init_pop))
  {
    warning("final_pop or init_pop is a vector, not a matrix. Attempting conversion to one-column matrices")
    final_pop <- matrix(final_pop, ncol = 1)
    init_pop <- matrix(init_pop, ncol = 1)
  }
  nobs = nrow(init_pop)
  if(length(start_times) != nobs || length(end_times) != nobs || nrow(final_pop) != nobs){
    stop("init_pop, start_times, end_times, and final_pop must all have the same number of rows!")
  }
  if(ncol(init_pop) != model$ntypes || ncol(final_pop) != model$ntypes){
    stop("init_pop, final_pop, and model must all have the same number of types!")
  }
  if(!is.numeric(init_pop) || !is.numeric(start_times) || !is.numeric(end_times) || !is.numeric(final_pop) || !is.numeric(initial)){
    stop("all time and population inputs must be numeric!")
  }


  

  # MLE
  loglik <- function(params){ -1*bp_loglik(model, params, init_pop, start_times, end_times, final_pop)}
  control <-  list(trace = trace, factr=10, pgtol=1e-20, fnscale = 1e7)  
  if(is.null(lower)){
    lower <- 1e-10*1:length(initial)
  }
  
  if(is.null(upper)){
    upper <- 1.75 + 1e-10*1:length(initial)
  }
  
  if(is.null(known)){
    mle <- optim(initial,
                  loglik, method = "L-BFGS-B",
                  lower = lower, upper = upper,
                  control = control)
  } else {
    mle <- bossMaps:::optifix(initial,
                               known,
                               loglik, method = "L-BFGS-B",
                               lower = lower, upper = upper,
                               control = control)
  }
  return(mle)
}

#' estimate
#'
#' Estimates the rate parameters for a general multitype branching process with constant rates
#'
#' @param model the \code{process_model} object representing the process generating the data
#' @param params the vector of parameters for which we are computing the likelihood
#' @param init_pop a \code{nobs x ntype} matrix with initial population for each observation
#' @param time the \code{nobs} length vector containing the time between the initial and final population observations
#' @param final_pop the \code{nobs x mtype} matrix of final populations observed
#' @param initial vector of initial parameters estimates for MLE optimization
#' @param known boolean vector of known parameter rates, if NULL, all rates will be estimated
#' @param lower vector of lower bounds on rate parameters for optimization
#' @param upper vector of upper bounds on rate parameters for optimization
#' @param trace level of output for optimizer - see optim function, control - trace for "L-BFGS-B" method
#'
#' @export
estimate = function(model, init_pop, times, final_pop, initial, known = NULL, lower = NULL, upper = NULL, trace = 1){
  for(trans in model$transition_list){
    if(!is_const(trans$rate$exp)){
      stop("for time-dependent models, use estimate_td")
    }
  }
  return(estimate_td(model, init_pop, max(times) - times, rep(max(times), length(times)), final_pop, initial, known, lower, upper, trace))
}
