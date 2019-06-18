library(MASS)

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
#'
#' @export
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
#'
#' @export
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

# library(gmbp)
#
# transitionList = TransitionList(FixedTransition(population = 0, rate = .1, fixed = c(2, 0)),
#                                 FixedTransition(population = 0, rate = .2, fixed = c(1, 1)),
#                                 FixedTransition(population = 0, rate = .3, fixed = c(0, 0)),
#                                 FixedTransition(population = 1, rate = .4, fixed = c(0, 2)),
#                                 FixedTransition(population = 1, rate = .5, fixed = c(0, 0)))
#
# parent = c()
# rate = c()
# offspring = matrix(nrow = length(transitionList), ncol = length(transitionList[[1]]$fixed))
# for(j in 1:length(transitionList)){
#   parent = c(parent, transitionList[[j]]$pop)
#   rate = c(rate, transitionList[[j]]$rate)
#   offspring[j,] = as.matrix(transitionList[[j]]$fixed)
# }
#
# parent = parent + 1
#
# time = 1
#
# N = c(500, 500)
#
# samp = sim_approx(n = 10000, t = time, N = N, parent = parent, rate = rate, offspring = offspring)
#
# library(gmbp)
# library(doParallel)
# library(foreach)
#
# # Set up parallelization
# cores=detectCores()
# cl <- makeCluster(cores[1]-1)
# registerDoParallel(cl)
#
# # Set up time to simulate - 1 unit
# time = 1
#
# # Set our initial population sizes
# initial = c(500,500)
#
#
# # No StopList items
# stopList = StopList()
#
# # Simulate 1000 samples
# ntrials = 10000
#
# full_res = foreach(j_=1:ntrials, .combine = "rbind", .packages = "gmbp") %dopar%{
#   res = branch(time, initial, transitionList, stopList, keep = FALSE, silent = TRUE)
#   as.matrix(res)
# }
# full_res = na.omit(full_res)
# data = as.matrix(full_res[full_res[,1] == 1,2:3])
#
# print(colMeans(data))
# print(colMeans(samp))
# print(var(data[,1]))
# print(var(samp[,1]))
