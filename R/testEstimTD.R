
test_sim2 = function(){  
  # Specify how many units of time to simulate
  time = 20
  
  # Specify that the simulation will initiate with a single type with size 500
  initial = c(500)
  
  # Specify two fixed transitions, birth and death
  transitionList = TransitionList(FixedTransition(population = 0, rate = Rate(params[1] + params[1]*sin(3.14159*params[2]*t)), fixed = c(2)),
                                  FixedTransition(population = 0,  rate = Rate(params[3]), fixed = c(0)))
  params = c(.3,2,.2)  
  res = branchTD(time, initial, transitionList, params, list(), dtime = .1, 5)
  plot(res[res$rep == 1,]$V3)
}