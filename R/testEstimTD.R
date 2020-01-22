
test_sim2 = function(){  
  # Specify how many units of time to simulate
  time = 20
  
  # Specify that the simulation will initiate with a single type with size 500
  initial = c(500)
  
  # Specify two fixed transitions, birth and death
  transitionList = TransitionList(FixedTransition(population = 0, rate = Rate(params[1]*exp(-params[3]*t)), fixed = c(2)),
                                  FixedTransition(population = 0,  rate = Rate(params[2]*exp(-params[3]*t)), fixed = c(0)))
  params = c(.3,.1,.6)  
  plot(branchTD(time, initial, transitionList, params, list(), dtime = .1))
}