#' estimate_td
#'
#' Estimates the rate parameters for a general multitype branching process with time-dependent rates
#'
#' @param model the \code{process_model} object representing the process generating the data
#' @param init_pop a \code{nobs x ntype} matrix with initial population for each observation
#' @param final_pop the \code{nobs x mtype} matrix of final populations observed
#' @param start_times the \code{nobs} length vector of times at which the initial populations were observed
#' @param end_times the \code{nobs} length vector of times at which the final populations were observed
#' @param initial_params vector of initial parameters estimates for MLE optimization
#'
#' @export
estimate_td = function(model, init_pop, final_pop, start_times, end_times,initial_params, control = list(), ...){
  if(class(model) != "estipop_process_model"){
    stop("model must be a process_model object!")
  }
  if( !is.numeric(init_pop) || !is.numeric(start_times) ||  !is.numeric(end_times) || !is.numeric(final_pop)){
    stop("all time and population inputs must be numeric!")
  }
  if (is.vector(final_pop) | is.vector(init_pop))
  {
    warning("final_pop or init_pop is a vector, not a matrix. Attempting conversion to matrices")
    final_pop <- matrix(final_pop, ncol = model$ntypes)
    init_pop <- matrix(init_pop, ncol = model$ntypes)
  }
  nobs = nrow(init_pop)
  if(length(start_times) != nobs || length(end_times) != nobs || nrow(final_pop) != nobs){
    stop("init_pop, start_times, end_times, and final_pop must all have the same number of rows!")
  }
  if(ncol(init_pop) != model$ntypes || ncol(final_pop) != model$ntypes){
    stop("init_pop, final_pop, and model must all have the same number of types!")
  }
  if(!is.numeric(init_pop) || !is.numeric(start_times) || !is.numeric(end_times) || !is.numeric(final_pop) || !is.numeric(initial_params)){
    stop("all time and population inputs must be numeric!")
  }
  if(any(init_pop < 0) || any(final_pop < 0) || any(start_times < 0) || any(end_times < 0) || any(end_times - start_times < 0)){
    stop("population and time variables must be nonnegative!")
  }

  

  # MLE
  loglik <- function(params){ -1*bp_loglik(model, params, init_pop, start_times, end_times, final_pop)}
  
  # Allowing the user to specify control variables and adding our own in
  default_control <-  list(trace = 1, factr=10, pgtol=1e-20, fnscale = 1e7)
  if(length(control) == 0) {
    control <- default_control
  } else {
    control <- c(control, default_control[setdiff(names(default_control), names(control))])
  }

  mle <- optim(initial_params,
                loglik,
                control = control, ...)
  return(mle)
}

#' estimate
#'
#' Estimates the rate parameters for a general multitype branching process with constant rates
#'
#' @param model the \code{process_model} object representing the process generating the data
#' @param params the vector of parameters for which we are computing the likelihood
#' @param init_pop a \code{nobs x ntype} matrix with initial population for each observation
#' @param final_pop the \code{nobs x mtype} matrix of final populations observed
#' @param time the \code{nobs} length vector containing the time between the initial and final population observations
#' @param initial_params vector of initial parameters estimates for MLE optimization
#'
#' @export
estimate = function(model, init_pop,  final_pop, times, initial_params, control = list(), ...){
  if(class(model) != "estipop_process_model"){
    stop("model must be a process_model object!")
  }
  for(trans in model$transition_list){
    if(!is_const(trans$rate$exp)){
      stop("for time-dependent models, use estimate_td")
    }
  }
  return(estimate_td(model, init_pop, final_pop, max(times) - times, rep(max(times), length(times)), initial_params, control = control, ...))
}
