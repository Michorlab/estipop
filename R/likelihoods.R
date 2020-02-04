#' moments
#' 
#' compute the moment matrices for a time-inhomogenous branching process model with a certain set of parameters. Helps for computing likelihoods.
#' 
#' @param model the \code{process_model} object representing the process whose moments will be computed
#' @param params the vector of parameters to plug into the process model during moment computation
#' @param start_times the various start times of the moments to be computed
#' @param end_times the various end times of the moments to be computed
#' 
#' @return a dataframe with the moments vector for each unique start and end combination
moments = function(model, params, start_times, end_times){
  ntype = model$ntypes
  
  rate_func <- function(t, params){
    lapply(model$transition_list, function(trans)eval(trans$rate$exp, list(t = t, params = params)))
  }
  
  # state is a vector containing all first and second moments evolving over time 
  # due to the nature of the mathematical derivation provided in the paper, we have to integrate the equation
  # *backward* in time, from the end time to the start time.
  moment_de <- function(curr_t, state, params){
    s = params[1] - curr_t
    rate_vec = rate_func(s, rate_params)
    
    #Expand rate vector to be a matrix where rate is in a colum corresponding to the parent
    rate_mat <- matrix(rep(0,nrow(offspring)*ntype), nrow = nrow(offspring), ncol = ntype)
    rate_mat[cbind(1:nrow(offspring), parent)] <- rate_vec
    lamb <- colSums(rate_mat)
    b_mat <- t(t(offspring)%*%rate_mat)/lamb
    
    a_mat <-  lamb*(b_mat - diag(ntype))
    
    c_mat <- array(rep(0,ntype**3), c(ntype, ntype, ntype)); #matrix of second derivatives of offspring PGF
    
    for(i in 1:ntype){
      i_mat <- t(t(offspring*offspring[,i])%*%rate_mat)/lamb
      i_mat[,i] <- i_mat[,i] - b_mat[,i]
      c_mat[,i,] <- t(i_mat)
      c_mat[i,,] <- t(i_mat)
    }
    
    mt_mat = matrix(state[1:ntype**2], c(ntype, ntype)) 
    dt_mat = matrix(state[(ntype**2+1):(ntype**3+ntype**2)], c(ntype, ntype*ntype)) 
    
    beta_mat <- matrix(rep(0,ntype**3), nrow = ntype, ncol = ntype*ntype); #beta_mat array for second moment ODE
    for(i in 1:ntype){
      ai <- lamb[i]
      beta_mat[i,] <- c(ai*t(mt_mat)%*%c_mat[,,i]%*%mt_mat) #vectorized computation of beta_mat
    }
    mt_prime <- c(a_mat%*%mt_mat)
    dt_prime <- c(a_mat%*%dt_mat + beta_mat)
    return(list(c(mt_prime, dt_prime)))
  }
  
  init_dt <- array(rep(0,ntype**3),c(ntype,ntype,ntype)) #inital values of second moments
  init_mt = array(rep(0, ntype**2), c(ntype, ntype))
  for(i in 1:ntype){
    init_dt[i,i,i] <- 1;
    init_mt[i,i] <- 1;
  }
  
  
  init_state <- c(c(init_mt),c(init_dt))
  ends = unique(end_times)
  out <- c()
  for(i in 1:length(ends)){
    starts = unique(start_times[end_times == ends[i]])
    sol <- deSolve::ode(y = init_state, c(0,sort(ends[i] - starts)), func = moment_de, parms = ends[i])
    out <- rbind(out, cbind(tf = ends[i], sol)[-1,])
  }
  return(data.frame(out))
}


#' bp_loglik
#' 
#' compute the log-likelihood of observing a set of data generated under a time-homogenous branching process model with some setting of parameters
#' 
#' @param model the \code{process_model} object representing the process generating the data
#' @param params the vector of parameters for which we are computing the likelihood
#' @param init_pop a \code{nobs x ntype} matrix with initial population for each observation
#' @param start_times the \code{nobs} length vector of times at which the initial populations were observed
#' @param end_times the \code{nobs} length vector of times at which the final populations were observed
#' @param final_pop the \code{nobs x mtype} matrix of final populations observed
#' 
#' @return a dataframe with the moments vector for each timepoint
bp_loglik = function(model, params, init_pop, start_times, end_times, final_pop){
  ntype = ncol(init_pop)
  mom = moments(model, params, start_times, end_times) #compute moments
  
  #Compute likelihood
  ll = 0
  for(obs in 1:nrow(final_pop)){
    pop <- final_pop[obs,]
    endtime <- end_times[obs]
    starttime <- start_times[obs]
    init <- init_pop[obs,]
    
    
    #get the state vector at the correct time
    desol <- as.matrix(mom[(mom$tf == endtime & out$time == endtime - starttime), -c(1,2)]) # = don't include time data in state vector 
    
    m_mat <- matrix(desol[1:ntype**2],nrow= ntype)
    dt_mat <- matrix(desol[(ntype**2 + 1):(ntype**3 + ntype**2)],nrow = ntype) #second moment array
    mu_vec <- t(m_mat)%*%init #final mean population vector
    sigma_mat <- matrix(init%*%dt_mat,c(ntype,ntype*ntype)) #final population covariance matrix
    
    for(i in 1:ntype){
      for(j in 1:ntype){
        sigma_mat[i,j] <-  sigma_mat[i,j] - init%*%(m_mat[,i]*m_mat[,j])
      }
    }
    
    #compute multivatiate normal likelihood
    ll <- ll - ntype/2*log(2*pi) - 1/2*log(det(sigma_mat))
    ll <- ll - 1/2*(mu_vec - pop)%*%solve(sigma_mat)%*%(mu_vec - pop)
  }  
  return(ll)
}
