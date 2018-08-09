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
  rlist = list(population, TRUE, rate, oVec, oDist, oParams)
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

#' StopCriterion
#'
#' Designates a condition upon which the simulation should stop
#'
#' If the sum of the specified states meets the inequality specified by the particular inequality and value, the simulation stops
#'
#' @param indices indices of states to be included in the sum
#' @param inequality ">", ">=
#' @param value value to compare the sum of states to
#'
#' @export
#' @examples
#' \dontrun{
#' StopCriterion(indices = c(0), inequality = ">=", value = 1000)
#' }
StopCriterion = function(indices, inequality, value){
  rlist = list(indices, inequality, value)
  names(rlist) = c("indices", "inequality", "value")
  return(rlist)
}

#' StopList
#'
#' Makes a single StopList object out of multiple StopCriterion objects
#'
#' @param ... StopCriterion objects
#'
#' @export
#' @examples
#' \dontrun{
#' StopList(StopCriterion(indices = c(0), inequality = ">=", value = 1000),
#'          StopCriterion(indices = c(0, 1), inequality = ">=", value = 10000))
#' }
StopList = function(...){
  ts = list(...)
  return(ts)
}

#' branch
#'
#' Makes a single StopList object out of multiple StopCriterion objects
#'
#' @param time number of time units to simulate
#' @param initial intial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param stopList StopList object specifying stopping conditions for the system
#'
#' @export
#' @examples
#' \dontrun{
#' branch(time = 100,
#'        initial = c(1,0),
#'        transistionList = TransitionList(Transition(prob = 0.5, fixed = c(1, 1, 0), random = c(FALSE, FALSE, TRUE)),
#'                                         Transition(prob = 0.5, fixed = c(0, 0, 0), random = c(TRUE, FALSE, FALSE))),
#'        stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 1000),
#'                   StopCriterion(indices = c(0, 1), inequality = ">=", value = 10000)))
#' }
branch = function(time, intial, transitionList, stopList){
  f = getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  gmbp3(time, f, initial, transitionList, stopList)
  res = read.csv(f, header = F)
  names(res)[1] = "time"
  return(res)
}






