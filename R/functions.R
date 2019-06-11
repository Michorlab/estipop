#### Model Specification Functions ####

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
FixedTransition = function(population, rate = 1.0, fixed){
  rlist = list(population, FALSE, rate, fixed)
  names(rlist) = c("pop", "is_random", "rate", "fixed")
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
  names(rlist) = c("pop", "is_random", "rate", "oVec", "oDist", "oParams")
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

#' Rate
#'
#' Designates a rate
#'
#' Specifies a rate.  User can specify using built-in types and accompanying parameters, or use a custom C++ function of time
#'
#' @param type type of rate function, see README
#' @param params list of parameters
#'
#' @export
#' @examples
#' \dontrun{
#' Rate(type = 0, params = c(0.01))
#' Rate(type = 1, params = c(0.1, -0.01))
#' }
Rate = function(type, params){
  rlist = list(type, params)
  names(rlist) = c("type", "params")
  return(rlist)
}


#### Simulation Functions ####

#' branch
#'
#' Makes a single StopList object out of multiple StopCriterion objects
#'
#' @param time number of time units to simulate
#' @param initial intial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param stopList StopList object specifying stopping conditions for the system
#' @param silent if true, verbose output will be shown.  Default: false
#' @param keep if true, the temporary comma-separated file generated while simulating while be kept.  if false, it will be deleted.  Default: false
#' @param seed seed for the random number generator.  If NULL, will use computer clock to set a random seed
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
branch = function(time, initial, transitionList, stopList, silent = FALSE, keep = FALSE, seed = NULL){
  f = R.utils::getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  if(is.null(seed)){
    gmbp3(time, f, initial, transitionList, stopList, silent)
  } else {
    gmbp3(time, f, initial, transitionList, stopList, silent, seed)
  }
  res = read.csv(f, header = F)
  names(res)[1] = "time"

  if(!keep)
    file.remove(f)

  return(res)
}

#' branchTD
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
#' branchTD(time = 100,
#'        initial = c(1,0),
#'        transistionList = TransitionList(Transition(prob = 0.5, fixed = c(1, 1, 0), random = c(FALSE, FALSE, TRUE)),
#'                                         Transition(prob = 0.5, fixed = c(0, 0, 0), random = c(TRUE, FALSE, FALSE))),
#'        stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 1000),
#'                   StopCriterion(indices = c(0, 1), inequality = ">=", value = 10000)))
#' }
branchTD = function(time, intial, transitionList, stopList, silent = FALSE, keep = FALSE){
  f = R.utils::getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  timeDepBranch(time, f, initial, transitionList, stopList, silent)
  res = read.csv(f, header = F)
  names(res)[1] = "time"

  if(!keep)
    file.remove(f)
  return(res)
}

#### Estimation Functions ####


#' estimateBP
#'
#' Estimates the rate parameters for a general multitype branching process
#'
#' @param time numeric or vector of time units for data observations
#' @param N intial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param data n x k data matrix
#' @param initial vector of initial estimates for MLE optimization
#' @param ... additional parameters to pass to optimizer
#'
#' @export
#' @examples
#' \dontrun{
#' estimateBP(time = 100,
#'        initial = c(1,0),
#'        transistionList = TransitionList(Transition(prob = 0.5, fixed = c(1, 1, 0), random = c(FALSE, FALSE, TRUE)),
#'                                         Transition(prob = 0.5, fixed = c(0, 0, 0), random = c(TRUE, FALSE, FALSE))),
#'        stopList = StopList(StopCriterion(indices = c(0), inequality = ">=", value = 1000),
#'                   StopCriterion(indices = c(0, 1), inequality = ">=", value = 10000)))
#' }
estimateBP = function(time, N, transitionList, data, initial){

  # Set up quantities for estimation, related to pgf
  parent = c()
  rate = c()
  offspring = matrix(nrow = length(transitionList), ncol = length(transitionList[[1]]$fixed))
  for(i in 1:length(transitionList)){
    parent = c(parent, transitionList[[i]]$pop)
    rate = c(rate, transitionList[[i]]$rate)
    offspring[i,] = as.matrix(transitionList[[i]]$fixed)
  }

  parent = parent + 1

  t = time

  # MLE
  loglik_ex2 <- function(rates) -1 * loglik_full2_time(data, t, N, parent, rates, offspring)

  # if (!is.na(...)){
  #   control = list(...)
  # } else {
  #   control = NULL
  # }
  control = NULL
  rMLE <- optim(initial,#b+runif(n = 1, min = -.3, max = .3),d+runif(n = 1, min = -.3, max = .3)),
                loglik_ex2, method = "L-BFGS-B",
                lower = 1e-10*1:length(initial), upper = 4 + 1e-10*1:length(initial),
                control = control)
  return(rMLE)
}





