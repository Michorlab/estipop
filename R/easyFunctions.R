#' simBirth
#'
#' Simulates a 1-type linear birth process
#'
#' @param birth constant birth rate
#' @param init initial population size
#' @param time number of time units to simulate
#' @param approx boolean if TRUE, approximate simulation is used, if FALSE, exact simulation is performed
#'
#' @export
#' @examples
#' \dontrun{
#' simBirth(birth = 1, init = 100, time = 10)
#' }
simBirth = function(birth = 1.0, init = 100, time = 1, approx = FALSE){
  transitionList = TransitionList(FixedTransition(population = 0, rate = birth, fixed = c(2)))

  stopList = StopList()

  return(estipop:::branch(time, init, transitionList, stopList, keep = FALSE, silent = TRUE, approx))
}

#' simBirthDeath
#'
#' Simulates a 1-type linear birth-death process
#'
#' @param birth constant birth rate
#' @param death constant death rate
#' @param init initial population size
#' @param time number of time units to simulate
#' @param approx boolean if TRUE, approximate simulation is used, if FALSE, exact simulation is performed
#'
#' @export
#' @examples
#' \dontrun{
#' simBirthDeath(birth = 1, death = 0.7, init = 100, time = 10)
#' }
simBirthDeath = function(birth = 1.0, death = 0.0, init = 100, time = 1, approx = FALSE){
  transitionList = TransitionList(FixedTransition(population = 0, rate = birth, fixed = c(2)),
                                  FixedTransition(population = 0, rate = death, fixed = c(0)))

  stopList = StopList()

  return(estipop:::branch(time, init, transitionList, stopList, keep = FALSE, silent = TRUE, approx))
}

#' simpBirthDeath
#'
#' Simulates a 1-type linear birth-death process
#'
#' @param b constant birth rate
#' @param d constant death rate
#' @param init initial population size
#' @param t simulate at which time
#' @param N how many samples
#'
#' @export
#' @examples
#' \dontrun{
#' simpBirthDeath(birth = 1, death = 0.7, init = 100, time = 10)
#' }
simpBirthDeath = function(b = 1.0, d = 0.0, init = 100, t = 1, N = 1){
  l = b - d
  alpha = (d * exp(l * t) - d) / (b * exp(l*t) - d)
  beta = (b * exp(l*t) - b) / (b * exp(l*t) - d)

  Y = rbinom(n = N, size = init, prob = 1 - alpha)
  add_x = function(y){
    return(y + rnbinom(n = 1, size = y, prob = 1 - beta))
  }
  X = sapply(Y, add_x)
  return(X)
}

#' simBirthMutation
#'
#' Simulates a 2-type birth-mutation process
#'
#' @param birth1 constant birth rate of population 1
#' @param mutation constant mutation rate
#' @param birth2 constant birth rate of population 2
#' @param init vector of initial population sizes
#' @param time number of time units to simulate
#' @param approx boolean if TRUE, approximate simulation is used, if FALSE, exact simulation is performed
#'
#' @export
#' @examples
#' \dontrun{
#' simBirthMutation(birth1 = 1, death1 = 0.7, mutation = 0.1, birth2 = 1.0, death2 = 0.7, init = c(100,0), time = 10)
#' }
simBirthMutation = function(birth1 = 1.0, mutation = 0.0, birth2 = 1.0, init = c(100,0), time = 1, approx = FALSE){
  transitionList = TransitionList(FixedTransition(population = 0, rate = birth1, fixed = c(2,0)),
                                  FixedTransition(population = 0, rate = mutation, fixed = c(1,1)),
                                  FixedTransition(population = 1, rate = birth2, fixed = c(0,2)))

  stopList = StopList()

  return(estipop:::branch(time, init, transitionList, stopList, keep = FALSE, silent = TRUE, approx))
}

#' simBirthDeathMutation
#'
#' Simulates a 2-type birth-death-mutation process
#'
#' @param birth1 constant birth rate of population 1
#' @param death1 constant death rate of population 1
#' @param mutation constant mutation rate
#' @param birth2 constant birth rate of population 2
#' @param death2 constant death rate of population 2
#' @param init vector of initial population sizes
#' @param time number of time units to simulate
#' @param approx boolean if TRUE, approximate simulation is used, if FALSE, exact simulation is performed
#'
#' @export
#' @examples
#' \dontrun{
#' simBirthDeathMutation(birth1 = 1, death1 = 0.7, mutation = 0.1, birth2 = 1.0, death2 = 0.7, init = c(100,0), time = 10)
#' }
simBirthDeathMutation = function(birth1 = 1.0, death1 = 0.0, mutation = 0.0, birth2 = 1.0, death2 = 0.7, init = c(100,0), time = 1, approx = FALSE){
  transitionList = TransitionList(FixedTransition(population = 0, rate = birth1, fixed = c(2,0)),
                                  FixedTransition(population = 0, rate = death1, fixed = c(0,0)),
                                  FixedTransition(population = 0, rate = mutation, fixed = c(1,1)),
                                  FixedTransition(population = 1, rate = birth2, fixed = c(0,2)),
                                  FixedTransition(population = 1, rate = death2, fixed = c(0,0)))

  stopList = StopList()

  return(estipop:::branch(time, init, transitionList, stopList, keep = FALSE, silent = TRUE, approx))
}

#' simTwoTypeResistance
#'
#' Simulates 2-type Resistance process
#'
#' @param birth0 constant birth rate of population 0
#' @param death0 constant death rate of population 0
#' @param mutation1 constant mutation rate to type-1 resistance
#' @param mutation2 constant mutation rate to type-2 resistance
#' @param mutation12 constant mutation rate from type-1 resistance to include type-2 resistance
#' @param mutation21 constant mutation rate from type-2 restistance to include type-1 resistance
#' @param birth1 constant birth rate of population 1
#' @param death1 constant death rate of population 1
#' @param birth2 constant birth rate of population 2
#' @param death2 constant death rate of population 2
#' @param birth12 constant birth rate of population 1/2
#' @param death12 constant death rate of population 1/2
#' @param init vector of initial population sizes
#' @param time number of time units to simulate
#' @param approx boolean if TRUE, approximate simulation is used, if FALSE, exact simulation is performed
#'
#' @export
#' @examples
#' \dontrun{
#' simTwoTypeResistance(birth1 = 1, death1 = 0.7, mutation = 0.1, birth2 = 1.0, death2 = 0.7, init = c(100,0), time = 10)
#' }
simTwoTypeResistance = function(mutation1 = 0.0, mutation2 = 0.0, mutation12 = 0.0, mutation21 = 0.0,
                                birth0 = 1.0, death0 = 0.7,
                                birth1 = 1.0, death1 = 0.7,
                                birth2 = 1.0, death2 = 0.7,
                                birth12 = 1.0, death12 = 0.7,
                                init = c(100,0,0,0), time = 1, approx = FALSE){
  transitionList = TransitionList(FixedTransition(population = 0, rate = mutation1, fixed = c(1, 1, 0, 0)),
                                  FixedTransition(population = 0, rate = mutation2, fixed = c(1, 0, 1, 0)),
                                  FixedTransition(population = 1, rate = mutation12, fixed = c(0, 1, 0, 1)),
                                  FixedTransition(population = 2, rate = mutation21, fixed = c(0, 0, 1, 1)),
                                  FixedTransition(population = 0, rate = birth0, fixed = c(2, 0, 0, 0)),
                                  FixedTransition(population = 0, rate = death0, fixed = c(0, 0, 0, 0)),
                                  FixedTransition(population = 1, rate = birth1, fixed = c(0, 2, 0, 0)),
                                  FixedTransition(population = 1, rate = death1, fixed = c(0, 0, 0, 0)),
                                  FixedTransition(population = 2, rate = birth2, fixed = c(0, 0, 2, 0)),
                                  FixedTransition(population = 2, rate = death2, fixed = c(0, 0, 0, 0)),
                                  FixedTransition(population = 3, rate = birth12, fixed = c(0, 0, 0, 2)),
                                  FixedTransition(population = 3, rate = death12, fixed = c(0, 0, 0, 0)))

  stopList = StopList()

  return(estipop:::branch(time, init, transitionList, stopList, keep = FALSE, silent = TRUE, approx))
}

