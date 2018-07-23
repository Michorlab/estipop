#' FixedTransition
#'
#' Designates a fixed transition, that is, a transition vector that does not contain random elements, such as a split event
#'
#' @param population index of population for which transition applies
#' @param rate rate at which this transition occurs
#' @param fixed
#'
#' @export
#' @examples
#' \dontrun{
#' FixedTransition(population = 0, rate = 0.5, fixed = c(1, 1, 0))
#' }
FixedTransition = function(population, rate, fixed){
  rlist = list(population, FALSE, rate, fixed)
  names(rlist) = c("pop", "is_random", "prob", "fixed")
  return(rlist)
}

#' RandomTransition
#'
#' Designates a fixed transition, that is, a transition vector that generates the number of offspring according to the distribution oDist with
#' parameters oParams, and then divides these according to a multinomial distribution with probability vector oVec
#'
#' @param population index of population for which transition applies
#' @param rate rate at which this transition occurs
#' @param oVec
#' @param oDist
#' @param oParams
#'
#' @export
#' @examples
#' \dontrun{
#' RandomTransition(population = 0, rate = 0.5, fixed = c(1, 1, 0))
#' }
RandomTransition = function(population, rate, oVec, oDist, oParams){
  rlist = list(population, FALSE, rate, oVec, oDist, oParams)
  names(rlist) = c("pop", "is_random", "prob", "oVec", "oDist", "oParams")
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


