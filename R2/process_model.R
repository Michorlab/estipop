#' new_rate
#' constructor for object of class \code{rate}
#' 
#' @param exprn an expression dictating the functional dependence of the rate on parameters and on time 
#' 
#' @return A new branching process rate object
new_rate <- function(exprn){
  r <- list("exp" = exprn)
  class(r) <- "estipop_rate"
  return(r)
}


#' validate_rate
#' verifies the correctness of an object of class \code{rate}
#' 
#' @return The \code{rate} object if it is valid, throws an error otherwise
validate_rate <- function(rate_obj){
  if(class(rate_obj) != "estipop_rate"){
    stop("not a valid estipop rate object!")
  }
  if(check_valid(rate_obj$exp)){
    return(rate_obj)
  }
  stop("invalid expression!")
}

#' rate
#' constructs and validates an object of class \code{rate}
#' 
#' @param exprn an expression definintion a functional dependence of rate on parameters and on time
#'
#' @return The \code{rate} object if it is valid, throws an error otherwise
#' @export
rate <- function(exprn){
  return(validate_rate(new_rate(substitute(exprn))))
}

#' new_transition
#' constructor for object of class \code{transition}
#' 
#' @param rate a rate for the transition 
#' @param parent the parent type for the transition
#' @param offspring the update to the system after the transition
#' 
#' @return A new branching process rate object
new_transition <- function(rate, parent, offspring){
  trans <- list("rate"= rate, "parent" = parent, "offspring" = offspring)
  class(trans) <- "estipop_transition"
  return(trans)
}

#' validate_transition
#' verifies the correctness of an object of class \code{transition}
#' 
#' @return The \code{transition} object if it is valid, throws an error otherwise
validate_transition <- function(trans_obj){
  if(class(trans_obj) != "estipop_transition"){
    stop("not a valid estipop transition object!")
  }
  if(class(trans_obj$rate) != "estipop_rate"){
    stop("rate is not a valid estipop rate object!")
  }
  if(!is.numeric(parent) || !is.numeric(offspring)){
    stop("parent and offpspring must be numeric!")
  }
  return(trans_obj)
}

#' transition 
#' constructs and validates an object of class \code{transition}
#' 
#' @param rate a rate for the transition 
#' @param parent the parent type for the transition
#' @param offspring the update to the system after the transition
#' 
#' @return The \code{transition} object if it is valid, throws an error otherwise
#' @export
transition <- function(rate, parent, offspring){
  return(validate_transition(new_transition(rate, parent, offspring)))
}

#' new_process_model
#' constructor for object of class \code{process_model}
#' 
#' @param transition_list a list of estipop transition objects
#' 
#' @return A new branching process model object
new_process_model <- function(transition_list){
  pm <- list(transition_list = transition_list, ntypes = length(transition_list[[1]]$offspring))
  class(pm) <- "estipop_process_model"
  return(pm)
}

#' validate_process_model
#' verifies the correctness of an object of class \code{process_model}
#' 
#' @return The \code{process_model} object if it is valid, throws an error otherwise
validate_process_model <- function(proc_model){
  if(!is.list(proc_model$transition_list) || length(proc_model$transition_list) == 0){
    stop("transition_list must be a list with a positive number of elements!")
  }
  lapply(proc_model$transition_list, function(b){if(class(b) != "estipop_transition"){stop("invalid transition object!")}})
  lapply(proc_model$transition_list, function(b){if(length(b$offpspring) != proc_model$ntypes){stop("transitions have inconsistent number of types!")}})
  if(class(proc_model) != "estipop_process_model"){
    stop("invalid process model!")
  }
  return(proc_model)
}


#' process_model
#' constructs and validates an object of class \code{process_model}
#' 
#' @param transitions a list of estipop transition objects
#'
#' @return The \code{process_model} object if it is valid, throws an error otherwise
#' @export
process_model <- function(transitions){
  return(validate_process_model(new_process_model(transitions)))
}
