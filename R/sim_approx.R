#' sim_approx
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t N-length vector of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
sim_approx <- function(numSamples, t, N, transitionList)
{

  parent = c()
  rate = c()
  offspring = matrix(nrow = length(transitionList), ncol = length(transitionList[[1]]$fixed))
  for(j in 1:length(transitionList)){
    parent = c(parent, transitionList[[j]]$pop)
    rate = c(rate, transitionList[[j]]$rate)
    offspring[j,] = as.matrix(transitionList[[j]]$fixed)
  }

  parent = parent + 1
  ntypes <- ncol(offspring)

  A <- matrix(ncol = ntypes, nrow = ntypes)
  for(i in 1:ntypes){
    ei <- rep(0, ntypes)
    ei[i] <- 1
    #b[i,] <- colSums(offspring[parent==i,] * prob[parent==i]) - ei
    A[i,] <- colSums(offspring[parent==i,,drop=F] * rate[parent==i]) - ei * sum(rate[parent == i])
  }
  #  A = diag(a) %*% b


  C <- array(rep(NA, ntypes^3), c(ntypes, ntypes, ntypes))
  for(i in 1:dim(C)[1])
  {
    tmp <- expand.grid(1:ntypes, 1:ntypes)
    for(j in 1:nrow(tmp))
    {
      #C[i,tmp[j,1], tmp[j,2]] <- sum(offspring[parent==i,tmp[j, 1]] * (offspring[parent==i,tmp[j, 2]] - I(tmp[j,1]==tmp[j,2])) * prob[parent==i]) * a[i]
      C[i,tmp[j,1], tmp[j,2]] <- sum(offspring[parent==i,tmp[j, 1],drop=F] * (offspring[parent==i,tmp[j, 2],drop=F] - I(tmp[j,1]==tmp[j,2])) * rate[parent==i])
    }
  }

  # Moments of the process --------------
  VL <- eigen(A)
  V <- VL$vectors
  L <- VL$values
  Vi <- solve(V)

  # Mean growth function
  m <- function(t, i, j)
  {
    val <- 0
    for(k in 1:ntypes)
    {
      val <- val + V[i,k] * exp(L[k]*t) * Vi[k,j]
    }
    val
  }

  beta <- function(t, i, j, k)
  {
    val <- 0;
    for(l in 1:ntypes)
    {
      for(n in 1:ntypes)
      {
        val <- val + C[i,l,n] * m(t, l, k) * m(t, n, j)
      }
    }
    val
  }

  d <- function(t, i, j, k)
  {
    val <- 0
    if(j == k) val <- val + m(t, i, j)

    for(n in 1:ntypes)
    {
      val <- val + integrate(function(s) m(t - s, i, n) * beta(s, n, j, k), 0, t)$value
    }
    val
  }

  nonSingular <- function(m) class(try(solve(m),silent=T))=="matrix"

  var_t <- function(t, j, k){
    val <- 0;
    for(n in 1:ntypes)
    {
      val <- val + N[n] * (d(t, n, j, k) - m(t, n, j) * m(t, n, k))
    }
    val
  }

  Mt <- matrix(ncol = ntypes, nrow = ntypes)
  Sigmat <- matrix(0, ncol = ntypes, nrow = ntypes)
  for(n1 in 1:ntypes)
  {
    for(n2 in 1:ntypes)
    {
      Mt[n1, n2] <- m(t, n1, n2)
      Sigmat[n1, n2] <- var_t(t, n1, n2)
    }
  }

  #print(N %*% Mt)
  #print(Sigmat)

  if(numSamples == 1){
    ret = t(as.matrix(round(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat)), ncol = ntypes))
    return(cbind(rep(t, numSamples), ret))
  } else{
    ret = t(as.matrix(round(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat)), ncol = ntypes))
    return(cbind(rep(t, numSamples), ret))
  }
}


#' sim_approx_timeseries
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t time to simulate until
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
#' @param observations a vector of observation times
#' @param dt make observations from time 0 until time t in units of time dt
#' 
#' @export
sim_approx_timeseries <- function(numSamples, t, N, transitionList, observations = NULL, dt = NULL){

  if(is.null(observations) & is.null(dt)){

    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))

    time = 1
    for(i in 1:t){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = 1, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + 1
    }

    ret_matrix <- na.omit(ret)
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  } else if (!is.null(observations)){
    d_obs = c(observations[1], diff(observations))

    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))
    
    # Remove 0 from our list of observations
    d_obs = d_obs[d_obs != 0]

    time = d_obs[1]
    for(i in 1:length(d_obs)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = time, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + d_obs[i+1]
    }
    ret_matrix <- na.omit(ret)
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  } else if(!is.null(dt)){
    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))

    # Change time = 0 to time = dt
    time = dt
    for(i in seq(0, t, dt)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = dt, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + dt
    }
    ret_matrix <- as.matrix(na.omit(ret))
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  }
}

#' sim_approx_timedep
#'
#' method to simulate from the asymptotic distriution of a time-inhomogenous branching process
#'
#' @param numSamples number of samples
#' @param t n-length vector of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
sim_approx_timedep <- function(numSamples, tf, t0, init_pop, parent, offspring, rate_func)
{
  
  ntype = ncol(offspring)
  
  #state is a vector containing all first and second moments evolving over time 
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
  
  sol <- deSolve::ode(y = init_state, c(0,tf - t0), func = moment_de, parms = tf)
  m_mat <- array(out[nrow(out),2:(ntype**2 + 1)],c(ntype,ntype))
  dt_mat <- array(out[nrow(out),(ntype**2 + 2):(ntype**3 + ntype**2 + 1)],c(ntype,ntype*ntype)) #second moment array
  
  mu_vec <- t(m_mat)%*%z0_vec #final mean population matrix
  sigma_mat <- matrix(z0_vec%*%dt_mat,c(ntype,ntype*ntype)) #final population covariance matrix
  
  for(i in 1:ntype){
    for(j in 1:ntype){
      sigma_mat[i,j] <-  sigma_mat[i,j] - init_pop%*%(m_mat[,i]*m_mat[,j])
    }
  }
  
  ret = t(as.matrix(round(MASS:::mvrnorm(n = numSamples, mu_vec, sigma_mat)), ncol = ntypes))
  return(cbind(rep(t, numSamples), ret))
}

#' sim_approx_timeseries_timedep
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t time to simulate until
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
#' @param observations a vector of observation times
#' @param dt make observations from time 0 until time t in units of time dt
#' 
#' @export
sim_approx_timeseries_timedep <- function(numSamples, t, N, transitionList, observations = NULL, dt = NULL){
  
  if(is.null(observations) & is.null(dt)){
    
    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)
    
    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))
    
    time = 1
    for(i in 1:t){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = 1, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + 1
    }
    
    ret_matrix <- na.omit(ret)
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  } else if (!is.null(observations)){
    d_obs = c(observations[1], diff(observations))
    
    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)
    
    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))
    
    # Remove 0 from our list of observations
    d_obs = d_obs[d_obs != 0]
    
    time = d_obs[1]
    for(i in 1:length(d_obs)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = time, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + d_obs[i+1]
    }
    ret_matrix <- na.omit(ret)
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  } else if(!is.null(dt)){
    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)
    
    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))
    
    # Change time = 0 to time = dt
    time = dt
    for(i in seq(0, t, dt)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx(numSamples = 1, t = dt, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + dt
    }
    ret_matrix <- as.matrix(na.omit(ret))
    if(any(ret_matrix[,-1] < 30)) warning("Population goes below 30. Normal approximation conditions may not hold.")
    return(ret_matrix)
  }
}
