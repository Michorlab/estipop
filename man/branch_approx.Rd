% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{branch_approx}
\alias{branch_approx}
\title{branch_approx
Approximately simulated a time-inhomogenous Markov branching process by sampling from the asymptotic multivariate normal distribution.}
\usage{
branch_approx(model, params, init_pop, time_obs, reps)
}
\arguments{
\item{model}{the \code{process_model} object representing the process being simulates}

\item{params}{the vector of parameters for which we are simulating the model}

\item{time_obs}{the vector of times at which to record the process state}

\item{reps}{the number of replicates to simulate}
}
\description{
branch_approx
Approximately simulated a time-inhomogenous Markov branching process by sampling from the asymptotic multivariate normal distribution.
}
