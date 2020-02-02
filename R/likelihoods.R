
#' bploglikelihood
#'
#' log-likelihood function for N replicates various time observations from a general k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
#' 
#' @export
bploglikelihood <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }
  #print(t)
  ntypes <- ncol(offspring)

  full_dat = as.data.frame(cbind(dat, t))
  full_dat = full_dat[order(full_dat[,ntypes+1]),]
  dat = as.matrix(full_dat[,1:ntypes])
  t = as.matrix(full_dat[,ntypes+1])

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
  m <- function(t_, i, j)
  {
    val <- 0
    for(k in 1:ntypes)
    {
      val <- val + V[i,k] * exp(L[k]*t_) * Vi[k,j]
    }
    val
  }

  beta <- function(t_, i, j, k)
  {
    val <- 0;
    for(l in 1:ntypes)
    {
      for(n in 1:ntypes)
      {
        val <- val + C[i,l,n] * m(t_, l, k) * m(t_, n, j)
      }
    }
    val
  }

  d <- function(t_, i, j, k)
  {
    val <- 0
    if(j == k) val <- val + m(t_, i, j)

    for(n in 1:ntypes)
    {
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_, subdivisions = 1000)$value
    }
    val
  }

  nonSingular <- function(m) class(try(solve(m),silent=T))=="matrix"

  var_t <- function(t_, j, k){
    val <- 0;
    for(n in 1:ntypes)
    {
      val <- val + N[n] * (d(t_, n, j, k) - m(t_, n, j) * m(t_, n, k))
    }
    val
  }

  ll = 0
  t_curr = -1
  Mt <- matrix(ncol = ntypes, nrow = ntypes)
  Sigmat <- matrix(0, ncol = ntypes, nrow = ntypes)
  for(ind in 1:nrow(dat)){
    if(ind == 1 || t[ind] != t[ind-1]){
      for(n1 in 1:ntypes){
        for(n2 in 1:ntypes){
          Mt[n1, n2] <- m(t[ind], n1, n2)
          Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
          t_curr = t[ind]
        }
      }
    }

    Sigma_inv <- solve(Sigmat)
    ll = ll - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - N%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N%*%Mt)
  }
  ll
}


#' bploglikelihood_TD
#'
#' log-likelihood function for N replicates from various time observations and possibly different initial counts from a general time-inhomogenous k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param tf times when observations in dat were made
#' @param t0 times when initial states were observed
#' @param N Nxk matrix of initial ancestor counts for each type for each observation
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate_func a function of time returning a d-length vector specifying the rate at which each offspring transition is occurring
#' @param rate_params a vector of parameters for the rate functions. These are the parameters than will be estimated 
#' @param offspring a dxk matrix specifying the offspring transitions
bploglikelihood_TD <- function(dat, tf, t0, N, parent, rate_func, rate_params, offspring){
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
  
  
  init_state <- c(c(init_mt),c(init_dt))
  ends = unique(tf)
  out <- c()
  for(i in 1:length(ends)){
    starts = unique(t0[tf == ends[i]])
    sol <- deSolve::ode(y = init_state, c(0,sort(ends[i] - starts)), func = moment_de, parms = ends[i])
    out <- rbind(out, cbind(tf = ends[i], sol)[-1,])
  }
  out <- data.frame(out)

  #Compute likelihood
  ll = 0
  for(obs in 1:nrow(dat)){
    pop <- dat[obs,]
    endtime <- tf[obs]
    starttime <- t0[obs]
    init_pop <- N[obs,]
    

    #get the state vector at the correct time
    desol <- as.matrix(out[(out$tf == endtime & out$time == endtime - starttime), -c(1,2)]) # = don't include time data in state vector 
    
    m_mat <- matrix(desol[1:ntype**2],nrow= ntype)
    dt_mat <- matrix(desol[(ntype**2 + 1):(ntype**3 + ntype**2)],nrow = ntype) #second moment array
    mu_vec <- t(m_mat)%*%init_pop #final mean population vector
    sigma_mat <- matrix(init_pop%*%dt_mat,c(ntype,ntype*ntype)) #final population covariance matrix

    for(i in 1:ntype){
      for(j in 1:ntype){
        sigma_mat[i,j] <-  sigma_mat[i,j] - init_pop%*%(m_mat[,i]*m_mat[,j])
      }
    }
    
    #compute multivatiate normal likelihood
    ll <- ll - ntype/2*log(2*pi) - 1/2*log(det(sigma_mat))
    ll <- ll - 1/2*(mu_vec - pop)%*%solve(sigma_mat)%*%(mu_vec - pop)
  }  
  ll
}
  

