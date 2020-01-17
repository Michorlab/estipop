test_sim <- function(){

  # Specify how many units of time to simulate
  time = 20

  # Specify that the simulation will initiate with a single type with size 100
  initial = c(500)

  compile_timedep("custom_rates/test_custom.cpp")
  compile_timedep("custom_rates/test_custom2.cpp")
  # Specify two fixed transitions, birth and death
  transitionList = TransitionList(FixedTransition(population = 0, rate = Rate(3, c("custom_rates/test_custom.so","rate")), fixed = c(2)),
                                FixedTransition(population = 0, rate = Rate(3, c("custom_rates/test_custom2.so","rate")), fixed = c(0)))
                                
  branchTD(time, initial, transitionList, list(), dtime = .1)
}