#' Transition
#'
#' Designates a Transition, that is, a transition vector that may contain or may not contain random elements
#'
#' @param population index of population for which transition applies
#' @param prob probability that this transition occurs
#' @param fixed
#' @param random
#'
#' @export
#' @examples
#' \dontrun{
#' Transition(population = 0, prob = 0.5, fixed = c(1, 1, 0), random = c(FALSE, FALSE, TRUE))
#' }
Transition = function(population, prob, fixed, random = NULL){
  rlist = list()
  # If there is no random component
  if(length(random) == 0){
    rlist = list(population, FALSE, prob, fixed)
    names(rlist) = c("pop", "is_random", "prob", "fixed")
  }

  # If there is a random component
  else {

    # check to make sure the fixed and random are of the same size
    if(length(random) != length(fixed))
      stop("fixed and random vector must be of same length")
    rlist = list(population, TRUE, prob, fixed, which(random)-1)
    names(rlist) = c("population", "is_random", "prob", "fixed", "rand_indices")
  }
  return(rlist)
}

#' TransitionList
#'
#' Makes a single TransitionList object out of multiple Transition objects
#'
#' @param ... Transition objects
#'
#' @export
#' @examples
#' \dontrun{
#' TransitionList(Transition(prob = 0.5, fixed = c(1, 1, 0), random = c(FALSE, FALSE, TRUE)),
#'                Transition(prob = 0.5, fixed = c(0, 0, 0), random = c(TRUE, FALSE, FALSE)))
#' }
TransitionList = function(...){
  ts = list(...)
  return(ts)
}


