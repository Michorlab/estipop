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
#' Designates a random transition, that is, a transition vector that generates the number of offspring according to the distribution oDist with
#' parameters oParams, and then divides these according to a multinomial distribution with probability vector oVec
#'
#' @param population index of population for which transition applies
#' @param rate rate at which this transition occurs
#' @param oVec offspring distribution probabilities
#' @param oDist offspring distribution function
#' @param oParams offspring distriution function parameters
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
#' Specifies a rate.  User writes an R expression which can involve time (t) or params[1..n] 
#'
#' @param expr the expression representing the function of time
#'
#' @export
#' @examples
#' \dontrun{
#' Rate(params[0]*t + params[1])
#' }
#' 
Rate = function(expr){
  return(substitute(expr))
}


#### Simulation Functions ####

#' branch
#'
#' Simulates a continuous-time markov branching process using the specified parameters.  Used for rates that are constant throughout time.
#'
#' @param time number of time units to simulate
#' @param initial intial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param stopList StopList object specifying stopping conditions for the system
#' @param silent if true, verbose output will be shown.  Default: false
#' @param keep if true, the temporary comma-separated file generated while simulating while be kept.  if false, it will be deleted.  Default: false
#' @param seed seed for the random number generator.  If NULL, will use computer clock to set a random seed
#' @param approx if true, simulation will proceed from draws from the asymptotic distribution
#' @param dtime if dtime and time are not null, observations will be made from time 0 to time in units of dtime
#' @param observations numeric vector specifying custom observation times
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
branch = function(time, initial, transitionList, stopList, silent = FALSE, keep = FALSE, seed = NULL, approx = FALSE, dtime = NULL, observations = NULL){
  if(approx){
    return(estipop:::sim_approx_full(1, time, initial, transitionList, observations, dtime))
  }

  time_obs = c()

  if(!is.null(observations)){
    time_obs = c(time_obs, observations)
  } else if (!is.null(dtime) & !is.null(time)){
    time_obs = seq(0, time, dtime)
  } else {
    time_obs = seq(0, time)
  }

  time_obs = as.matrix(time_obs)

  f = R.utils::getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  if(is.null(seed)){
    gmbp3(time_obs, f, initial, transitionList, stopList, silent)
  } else {
    gmbp3(time_obs, f, initial, transitionList, stopList, silent, seed)
  }
  res = read.csv(f, header = F)
  names(res)[1] = "time"

  if(!keep)
    file.remove(f)

  return(res)
}

#' branchTD
#' TODO: Parse the R expression into CPP for fast sim
#' Simulates a continuous-time markov branching process using the specified parameters.  Used for time-dependent rate simulation.
#'
#' @param time number of time units to simulate
#' @param initial intial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param stopList StopList object specifying stopping conditions for the system
#' @param silent if true, verbose output will be shown.  Default: false
#' @param keep if true, the temporary comma-separated file generated while simulating while be kept.  if false, it will be deleted.  Default: false
#' @param seed seed for the random number generator.  If NULL, will use computer clock to set a random seed
#' @param dtime if dtime and time are not null, observations will be made from time 0 to time in units of dtime
#' @param observations numeric vector specifying custom observation times
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
branchTD = function(time, initial, transitionList, transitionParams, stopList, reps, silent = FALSE, keep = FALSE, seed = NULL, dtime = NULL, observations = NULL){

  time_obs = c()

  if(!is.null(observations)){
    time_obs = c(time_obs, observations)
  } else if (!is.null(dtime) & !is.null(time)){
    time_obs = seq(0, time, dtime)
  } else {
    time_obs = seq(0, time)
  }
  
  time_obs = as.matrix(time_obs)
  for(i in 1:length(transitionList)){
    if(isConst(transitionList[[i]]$rate)){
      transitionList[[i]]$rate = eval(transitionList[[i]]$rate, list(params = transitionParams))
    }
    else{
      fname =  paste("custom_rate_", digest::digest(deparse(transitionList[[i]]$rate),"md5"), sep="")
      create_timedep_template(transitionList[[i]]$rate, transitionParams, paste(fname, ".cpp", sep=""))
      compile_timedep( paste(fname, ".cpp", sep=""))
      transitionList[[i]]$rate = list(type=3,params = c(paste(fname, ".so", sep=""), "rate"))
    }
  }

  f = R.utils::getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  if(is.null(seed)){
    timeDepBranch(time_obs, reps, f, initial, transitionList, stopList, silent)
  } else {
    timeDepBranch(time_obs, reps, f, initial, transitionList, stopList, silent, seed)
  }
  res = read.csv(f, header = F)
  names(res)[1:2] = c("rep","time")

  if(!keep)
    file.remove(f)
    for(trans in transitionList){
      if("type" %in% names(trans$rate)){
        fname = .pop(trans$rate$params[1],".so")
        file.remove(trans$rate$params[1])
        file.remove(paste(fname, ".cpp", sep=""))
        file.remove(paste(fname, ".h", sep=""))
        file.remove(paste(fname, ".o", sep=""))
        file.remove(paste(fname, ".cpp.backup", sep=""))
      }
    }

  return(res)
}

#### Estimation Functions ####


#' estimateBP
#'
#' Estimates the rate parameters for a general multitype branching process with time-dependent rates
#'
#' @param time numeric or vector of time units for data observations
#' @param N intial state vector
#' @param initTime time of observation for initial state vector
#' @param transitionList TransitionList object specifying transitions in system
#' @param data n x k data matrix
#' @param initial vector of initial estimates for MLE optimization
#' @param known boolean vector of known parameter rates, if NULL, all rates will be estimated
#' @param lower vector of lower bounds on rate parameters for optimization
#' @param upper vector of upper bounds on rate parameters for optimization
#' @param trace level of output for optimizer - see optim function, control - trace for "L-BFGS-B" method
#'
#' @export
estimateBP = function(time, N, initTime, transitionList, data, initial, known = NULL, lower = NULL, upper = NULL, trace = 1){
  
  if(length(transitionList) < 1){
    stop("No model specified by transitionList.")
  }
  
  if (is.vector(data))
  {
    warning("data is a vector, not a matrix. Converting to a one-column matrix")
    data <- matrix(data, ncol = 1)
    N <- matrix(N, ncol = 1)
  }
  
  for(i in 1:length(transitionList)){
    if(length(transitionList[[i]]$fixed) != ncol(data)){
      stop("Model specified by transitionList has different number of types (columns) than the data matrix.")
    }
  }
  
  # Set up quantities for estimation, related to pgf
  parent = c()
  rate_list = c() #list of all of the rates. Can include both numerics and rate objects
  offspring = matrix(nrow = length(transitionList), ncol = length(transitionList[[1]]$fixed))
  
  for(i in 1:length(transitionList)){
    parent = c(parent, transitionList[[i]]$pop)
    rate_list = c(rate_list, transitionList[[i]]$rate)
    offspring[i,] = as.matrix(transitionList[[i]]$fixed)
  }
  
  rate_func = function(t, params){
    sapply(rate_list, function(r)eval(r, list(t = t, params = params)))
  }

  parent = parent + 1
  
  t = time
  
  # MLE
  loglik <- function(params){ loglik_time_dependent(data, t, initTime, N, parent, rate_func, params, offspring)}
  control =  list(trace = trace, fnscale = 1e7)
  
  if(is.null(lower)){
    lower = 1e-10*1:length(initial)
  }
  
  if(is.null(upper)){
    upper = 1.75 + 1e-10*1:length(initial)
  }
  
  if(is.null(known)){
    rMLE <- optim(initial,
                  loglik, method = "L-BFGS-B",
                  lower = lower, upper = upper,
                  control = control)
  } else {
    rMLE <- bossMaps:::optifix(initial,
                               known,
                               loglik_ex2, method = "L-BFGS-B",
                               lower = lower, upper = upper,
                               control = control)
  }
  return(rMLE)
}




