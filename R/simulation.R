#' branch
#' Simulates a continuous-time time-inhomogenous markov branching process using the specified parameters. Uses C++ code for faster simulation.
#'
#' @param model the \code{process_model} object representing the process being simulates
#' @param params the vector of parameters for which we are simulating the model
#' @param time_obs the vector of times at which to record the process state
#' @param stops \code{stop_list} object specifying stopping conditions for the system
#' @param reps the number of replicates to simulate
#' @param silent if true, verbose output will be shown.  Default: false
#' @param keep if true, the temporary comma-separated file generated while simulating while be kept.  if false, it will be deleted.  Default: false
#' @param seed seed for the random number generator.  If NULL, will use computer clock to set a random seed
#'
#' @export
branch <- function(model, params, init_pop, time_obs, stops, reps, silent = FALSE, keep = FALSE, seed = NULL){
  if(class(model) != "estipop_process_model"){
    stop("model must be a process_model object!")
  }
  if(ncol(init_pop) != model$ntypes){
    stop("init_pop and model must have same number of types ")
  }
  if(!is.numeric(params) || !is.numeric(init_pop) || !is.numeric(time_obs) || !is.numeric(reps)){
    stop("all time, population, and parameter inputs must be numeric!")
  }
  
  #modify the transition list with metadata about whether the rate is constant
  for(i in 1:length(model$transition_list)){
    if(is_const(model$transition_list[[i]]$rate$exp)){
      #evaluate constant rates
      model$transition_list[[i]]$rate <- eval(transitionList[[i]]$rate, list(params = transitionParams))
      model$transition_list[[i]]$type <- 1
    }
    else{
      fname <-  paste("custom_rate_", digest::digest(deparse(transitionList[[i]]$rate),"md5"), sep="")
      create_timedep_template(transitionList[[i]]$rate, transitionParams, paste(fname, ".cpp", sep=""))
      compile_timedep( paste(fname, ".cpp", sep=""))
      model$transition_list[[i]]$rate <- c(paste(fname, ".so", sep=""), "rate")
      model$transition_list[[i]]$type <- 2
    }
  }
  
  f <- R.utils::getAbsolutePath(tempfile(pattern = paste("system_", format(Sys.time(), "%d-%m-%Y-%H%M%S"), "_", sep = ""), fileext = ".csv", tmpdir = getwd()))
  if(is.null(seed)){
    timeDepBranch(time_obs, reps, f, initial, model$transition_list, stops, silent)
  } else {
    timeDepBranch(time_obs, reps, f, initial, model$transition_list, stops, silent, seed)
  }
  res <- read.csv(f, header = F)
  names(res)[1:2] <- c("rep","time")
  
  if(!keep){
    file.remove(f)
  }
  
  #remove the temporary dynamic libraries
  for(trans in model$transition_list){
    if(trans$rate$type == 2){
      fname = .pop(trans$rate,".so")
      file.remove(trans$rate)
      file.remove(paste(fname, ".cpp", sep=""))
      file.remove(paste(fname, ".h", sep=""))
      file.remove(paste(fname, ".o", sep=""))
      file.remove(paste(fname, ".cpp.backup", sep=""))
    }
  }
  
  return(res)
}

#' branch_approx
#' Approximately simulated a time-inhomogenous Markov branching process by sampling from the asymptotic multivariate normal distribution.
#'
#' @param model the \code{process_model} object representing the process being simulates
#' @param params the vector of parameters for which we are simulating the model
#' @param time_obs the vector of times at which to record the process state
#' @param stops \code{stop_list} object specifying stopping conditions for the system
#' @param reps the number of replicates to simulate
#'
#' @export
branch_approx = function(model, params, init_pop, time_obs, stops, reps){
  if(class(model) != "estipop_process_model"){
    stop("model must be a process_model object!")
  }
  if(ncol(init_pop) != model$ntypes){
    stop("init_pop and model must have same number of types ")
  }
  if(!is.numeric(params) || !is.numeric(init_pop) || !is.numeric(time_obs) || !is.numeric(reps)){
    stop("all time, population, and parameter inputs must be numeric!")
  }
  
  ntype = model$ntypes  
  time_obs <- sort(unique(c(0, time_obs)))
  res <- matrix(rep(0, reps*ntype*length(time_obs)), ncol = ntype)
  mom <- moments(model, params, time_obs[1:(length(time_obs - 1))], time_obs[-1])
  
  for(r in 1:reps){
    curr_pop <- init_pop
    res[(r-1)*length(time_obs)+ 1,] <- c(time = 0, curr_pop)
    
    for(i in 2:length(time_obs)){
      desol <- as.matrix(mom[(mom$tf == time_obs[i]), -c(1,2)])
      
      m_mat <- matrix(desol[1:ntype**2],nrow= ntype)
      dt_mat <- matrix(desol[(ntype**2 + 1):(ntype**3 + ntype**2)],nrow = ntype) #second moment array
      mu_vec <- t(m_mat)%*%curr_pop #final mean population vector
      sigma_mat <- matrix(curr_pop%*%dt_mat,c(ntype,ntype*ntype)) #final population covariance matrix
      
      for(i in 1:ntype){
        for(j in 1:ntype){
          sigma_mat[i,j] <-  sigma_mat[i,j] - curr_pop%*%(m_mat[,i]*m_mat[,j])
        }
      }
      
      curr_pop <- round(MASS:::mvrnorm(n = 1, mu_vec, sigma_mat))
      res[(r-1)*length(time_obs) + i,] <- c(time = time_obs[i], rep = r, curr_pop)
    }
  }
}