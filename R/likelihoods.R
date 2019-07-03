#### Estimation Log-likelihood functions ####

#' loglik_est
#'
#' log-likelihood function for N replicates from the same time observation from a general k-type branching process with d transitions
#'
#' @param dat nxK matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_est <- function(dat, t, N, parent, rate, offspring)
{
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

  if(!nonSingular(Sigmat)){
    write.csv(Sigmat, "singular_sigmat.csv")
    write.csv(dat, "singular_dat.csv")
    write.csv(rate, "singular_rates.csv")
    print("SINGULAR!")
    print(rate)
    print(A)
    return(-205561842*2)
  }

  Sigma_inv <- solve(Sigmat)
  ### ADDED LINE ###
  dat = as.matrix(dat)
  ll = - nrow(dat) / 2 * log(det(Sigmat))
  for(ind in 1:nrow(dat))
  {
    ll <- ll - 1/2 * (dat[ind,] - N%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N%*%Mt)
  }
  ll
}

#' loglik_est_time
#'
#' log-likelihood function for N replicates various time observations from a general k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_est_time <- function(dat, t, N, parent, rate, offspring)
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

#' loglik_est_time_piece
#'
#' log-likelihood function for N replicates from various time observations and possibly different initial counts from a general k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N Nxk matrix of initial ancestor counts for each type for each observation
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_est_time_piece <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }


  #print(t)
  ntypes <- ncol(offspring)

  full_dat = as.data.frame(cbind(dat, N, t))
  full_dat = full_dat[order(full_dat[,(2*ntypes)+1]),]
  dat = as.matrix(full_dat[,1:ntypes])
  N = as.matrix(full_dat[,(ntypes+1):(2*ntypes)])
  t = as.matrix(full_dat[,(2*ntypes)+1])

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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_, stop.on.error = F)$value
      #val = val + quadl(function(s) m(t_ - s, i, n) * beta(s, n, j, k),0,t_)
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
  for(ind in 1:nrow(dat))
  {
    if(t_curr == -1 || t[ind] != t[ind-1]){
      for(n1 in 1:ntypes)
      {
        for(n2 in 1:ntypes)
        {
          Mt[n1, n2] <- m(t[ind], n1, n2)
          Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
        }
      }
    }

    Sigma_inv <- solve(Sigmat)
    ll = ll - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - N[ind,]%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N[ind,]%*%Mt)
  }
  ll
}

#' loglik_est_time_z
#'
#' special log-likelihood function for the observing the convolution of two populations
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_est_time_z <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }
  #print(t)
  ntypes <- ncol(offspring)

  full_dat = as.data.frame(cbind(dat, t))
  full_dat = full_dat[order(full_dat[,ntypes]),]
  dat = as.matrix(full_dat[,1])
  t = as.matrix(full_dat[,ntypes])

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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
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

  for(ind in 1:nrow(dat))
  {
    if(t_curr == -1 || t[ind] != t[ind-1]){
      for(n1 in 1:ntypes)
      {
        for(n2 in 1:ntypes)
        {
          Mt[n1, n2] <- m(t[ind], n1, n2)
          Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
        }
      }
    }

    Sigmat = as.matrix(Sigmat[1,1] + Sigmat[2, 2] + 2 * Sigmat[1,2])
    Sigma_inv <- solve(Sigmat)
    ll = ll - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - sum(N%*%Mt)) %*% Sigma_inv %*% t(dat[ind,] - sum(N%*%Mt))
  }
  ll
}

#' loglik_est_z
#'
#' special log-likelihood function for the observing the convolution of two populations
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_est_z <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }
  #print(t)
  ntypes <- ncol(offspring)
  #a <- rep(NA, ntypes)
  # for(i in 1:ntypes)
  # {
  #   a[i] <- sum(rate[parent == i])
  # }

  #prob <- unlist(tapply(rate, parent, function(x) x / sum(x), simplify = F))

  #b <- matrix(rep(NA, ntypes^2), nrow = ntypes)
  #b2 <- matrix(rep(NA, ntypes^2), nrow = ntypes)
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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
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
  for(ind in 1:nrow(dat))
  {
    Mt <- matrix(ncol = ntypes, nrow = ntypes)
    Sigmat <- matrix(0, ncol = ntypes, nrow = ntypes)
    for(n1 in 1:ntypes)
    {
      for(n2 in 1:ntypes)
      {
        Mt[n1, n2] <- m(t[ind], n1, n2)
        Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
      }
    }

    #Mt = as.matrix(sum(Mt))
    Sigmat = as.matrix(Sigmat[1,1] + Sigmat[2, 2] + 2 * Sigmat[1,2])
    #print(Sigmat)
    #print(N)

    # if(!nonSingular(Sigmat)){
    #   print("SINGULAR!")
    #   print(rate)
    #   print(A)
    #   return(-205561842*2)
    # }

    Sigma_inv <- solve(Sigmat)
    ll = ll - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - sum(N%*%Mt)) %*% Sigma_inv %*% t(dat[ind,] - sum(N%*%Mt))
  }
  ll
}

#### Full Log-likelihood functions ####

#' loglik_full
#'
#' log-likelihood function for N replicates from the same time observation from a general k-type branching process with d transitions
#'
#' @param dat nxK matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_full <- function(dat, t, N, parent, rate, offspring)
{
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

  Sigma_inv <- solve(Sigmat)
  ### ADDED LINE ###
  dat = as.matrix(dat)
  ll = -nrow(dat) * ntypes / 2 * log(2*pi) - nrow(dat) / 2 * log(det(Sigmat))
  for(ind in 1:nrow(dat))
  {
    ll <- ll - 1/2 * (dat[ind,] - N%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N%*%Mt)
  }
  ll
}

#' loglik_full_time
#'
#' log-likelihood function for N replicates various time observations from a general k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_full_time <- function(dat, t, N, parent, rate, offspring)
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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
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
    #print(paste(ind, ": ", sep = ""))
    #print(Mt)
    #print(dat[ind,])
    #print(N)
    ll = ll - ntypes / 2 * log(2*pi) - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - N%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N%*%Mt)
    #print((dat[ind,] - N%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N%*%Mt))
    #print((dat[ind,] - N%*%Mt))
    #print(Sigma_inv)
    #print(t(dat[ind,] - N%*%Mt))
  }
  ll
}

#' loglik_full_time_piece
#'
#' log-likelihood function for N replicates from various time observations and possibly different initial counts from a general k-type branching process with d transitions
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N Nxk matrix of initial ancestor counts for each type for each observation
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_full_time_piece <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }


  #print(t)
  ntypes <- ncol(offspring)

  full_dat = as.data.frame(cbind(dat, N, t))
  full_dat = full_dat[order(full_dat[,(2*ntypes)+1]),]
  dat = as.matrix(full_dat[,1:ntypes])
  N = as.matrix(full_dat[,(ntypes+1):(2*ntypes)])
  t = as.matrix(full_dat[,(2*ntypes)+1])

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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
      #val = val + quadl(function(s) m(t_ - s, i, n) * beta(s, n, j, k),0,t_)
    }
    val
  }

  nonSingular <- function(m) class(try(solve(m),silent=T))=="matrix"

  var_t <- function(t_, j, k,ind){
    val <- 0;
    for(n in 1:ntypes)
    {
      val <- val + N[ind, n] * (d(t_, n, j, k) - m(t_, n, j) * m(t_, n, k))
    }
    val
  }

  ll = 0
  t_curr = -1
  Mt <- matrix(ncol = ntypes, nrow = ntypes)
  Sigmat <- matrix(0, ncol = ntypes, nrow = ntypes)
  for(ind in 1:nrow(dat))
  {
    if(ind == 1 || t[ind] != t[ind-1] || N[ind,] != N[ind-1,]){
      for(n1 in 1:ntypes)
      {
        for(n2 in 1:ntypes)
        {
          Mt[n1, n2] <- m(t[ind], n1, n2)
          Sigmat[n1, n2] <- var_t(t[ind], n1, n2, ind)
        }
      }
    }
    #print(paste(ind, ": ", sep = ""))
    #print(Mt)
    #print(dat[ind,])
    #print(N[ind,])
    Sigma_inv <- solve(Sigmat)
    ll = ll - ntypes / 2 * log(2*pi) - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - N[ind,]%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N[ind,]%*%Mt)
    #print((dat[ind,] - N[ind,]%*%Mt) %*% Sigma_inv %*% t(dat[ind,] - N[ind,]%*%Mt))
    #print((dat[ind,] - N[ind,]%*%Mt))
    #print(Sigma_inv)
    #print(t(dat[ind,] - N[ind,]%*%Mt))
  }
  ll
}

#' loglik_full_time_z
#'
#' special log-likelihood function for the observing the convolution of two populations
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_full_time_z <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }
  #print(t)
  ntypes <- ncol(offspring)

  full_dat = as.data.frame(cbind(dat, t))
  full_dat = full_dat[order(full_dat[,ntypes]),]
  dat = as.matrix(full_dat[,1])
  t = as.matrix(full_dat[,ntypes])

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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
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

  for(ind in 1:nrow(dat))
  {
    if(t_curr == -1 || t[ind] != t[ind-1]){
      for(n1 in 1:ntypes)
      {
        for(n2 in 1:ntypes)
        {
          Mt[n1, n2] <- m(t[ind], n1, n2)
          Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
        }
      }
    }

    Sigmat = as.matrix(Sigmat[1,1] + Sigmat[2, 2] + 2 * Sigmat[1,2])
    Sigma_inv <- solve(Sigmat)
    ll = ll - ntypes / 2 * log(2*pi) - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - sum(N%*%Mt)) %*% Sigma_inv %*% t(dat[ind,] - sum(N%*%Mt))
  }
  ll
}

#' loglik_full_z
#'
#' special log-likelihood function for the observing the convolution of two populations
#'
#' @param dat Nxk matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
loglik_full_z <- function(dat, t, N, parent, rate, offspring)
{
  if(length(t) == 1){
    t = rep(t, nrow(dat))
  }
  #print(t)
  ntypes <- ncol(offspring)
  #a <- rep(NA, ntypes)
  # for(i in 1:ntypes)
  # {
  #   a[i] <- sum(rate[parent == i])
  # }

  #prob <- unlist(tapply(rate, parent, function(x) x / sum(x), simplify = F))

  #b <- matrix(rep(NA, ntypes^2), nrow = ntypes)
  #b2 <- matrix(rep(NA, ntypes^2), nrow = ntypes)
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
      val <- val + integrate(function(s) m(t_ - s, i, n) * beta(s, n, j, k), 0, t_)$value
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
  for(ind in 1:nrow(dat))
  {
    Mt <- matrix(ncol = ntypes, nrow = ntypes)
    Sigmat <- matrix(0, ncol = ntypes, nrow = ntypes)
    for(n1 in 1:ntypes)
    {
      for(n2 in 1:ntypes)
      {
        Mt[n1, n2] <- m(t[ind], n1, n2)
        Sigmat[n1, n2] <- var_t(t[ind], n1, n2)
      }
    }

    #Mt = as.matrix(sum(Mt))
    Sigmat = as.matrix(Sigmat[1,1] + Sigmat[2, 2] + 2 * Sigmat[1,2])
    #print(Sigmat)
    #print(N)

    # if(!nonSingular(Sigmat)){
    #   print("SINGULAR!")
    #   print(rate)
    #   print(A)
    #   return(-205561842*2)
    # }

    Sigma_inv <- solve(Sigmat)
    ll = ll - ntypes / 2 * log(2*pi) - 1 / 2 * log(det(Sigmat))
    ll <- ll - 1/2 * (dat[ind,] - sum(N%*%Mt)) %*% Sigma_inv %*% t(dat[ind,] - sum(N%*%Mt))
  }
  ll
}

#### User-facing functions ####

#' bploglikelihood
#'
#' log-likelihood function for N replicates from a general k-type branching process with d transitions
#'
#' @param dat nxK matrix of observed data where each column is a count from a type and each row is an observation
#' @param t time from system initialization until observations in dat were made
#' @param N initial ancestor counts for each type
#' @param transitionList A TransitionList object to specify the model form
#'
#' @export
bploglikelihood = function(dat, t, N, transitionList){

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

  ndata = nrow(dat)
  ntypes = ncol(offspring)

  if(length(t) == 1 & is.null(dim(N))){
    #print("Sit1")
    return(estipop:::loglik_full(dat, t, N, parent, rate, offspring))
  } else if (length(t) != 1 & is.null(dim(N))){
    #print("Sit2")
    return(estipop:::loglik_full_time(dat, t, N, parent, rate, offspring))
  } else{
    #print("Sit3")
    return(estipop:::loglik_full_time_piece(dat, t, N, parent, rate, offspring))
  }

}

