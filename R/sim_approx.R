#' sim_approx
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param parent a d-length vector specifying which population type is associated with a specific offspring transition
#' @param rate a d-length vector specifying the rate at which each offspring transition is occurring
#' @param offspring a dxk matrix specifying the offspring transitions
sim_approx <- function(numSamples, t, N, parent, rate, offspring)
{

  # print("Parameters supplied:")
  # print(paste("t:", t))
  # print(paste("N:", N))
  # print(paste("parent:", parent))
  # print(paste("rate:", rate))
  # print(paste("offspring:", offspring))
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
    ret = t(as.matrix(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat), ncol = ntypes))
    return(cbind(rep(t, numSamples), ret))
  } else{

    return(cbind(rep(t, numSamples), as.matrix(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat))))
  }
}

#' sim_approx2
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t N-length vecotr of time from system initialization until each observation in dat were made
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
sim_approx2 <- function(numSamples, t, N, transitionList)
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
    ret = t(as.matrix(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat), ncol = ntypes))
    return(cbind(rep(t, numSamples), ret))
  } else{

    return(cbind(rep(t, numSamples), as.matrix(MASS:::mvrnorm(n = numSamples, N %*% Mt, Sigmat))))
  }
}


#' sim_approx_full
#'
#' method to simulate from the asymptotic distriution of a general multitype branching process
#'
#' @param numSamples number of samples
#' @param t time to simulate until
#' @param N k-length vector of the initial ancestor counts for each type
#' @param transitionList transitionList object to specify the model
#' @param observations a vector of observation times
#' @param dt make observations from time 0 until time t in units of time dt
sim_approx_full <- function(numSamples, t, N, transitionList, observations = NULL, dt = NULL){

  if(is.null(observations) & is.null(dt)){

    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))

    time = 1
    for(i in 1:t){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx2(numSamples = 1, t = 1, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + 1
    }

    return(na.omit(ret))
  } else if (!is.null(observations)){
    d_obs = c(observations[1], diff(observations))

    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))

    time = d_obs[1]
    for(i in 1:length(d_obs)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx2(numSamples = 1, t = time, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + d_obs[i+1]
    }

    return(na.omit(ret))
  } else if(!is.null(dt)){

    prev_values = matrix(rep(N, numSamples), ncol = length(N), byrow = T)

    ret = matrix(ncol = length(transitionList[[1]]$fixed )+1)
    ret = rbind(ret, c(0, N))

    time = 0
    for(i in seq(0, t, dt)){
      for(j in 1:nrow(prev_values)){
        samp = sim_approx2(numSamples = 1, t = dt, N = prev_values[j,], transitionList)
        prev_values[j,] = samp[,2:ncol(samp)]
        ret = rbind(ret, c(time, samp[,2:ncol(samp)]))
      }
      time = time + dt
    }

    return(as.matrix(na.omit(ret)))
  }
}
