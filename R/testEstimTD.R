
test_sim2 = function(){  
  # Specify how many units of time to simulate
  time = 20
  
  # Specify that the simulation will initiate with a single type with size 500
  initial = c(5000)
  
  # Specify two fixed transitions, birth and death
  transitionList = TransitionList(FixedTransition(population = 0, rate = Rate(params[1] - params[2]*exp(-params[3]*t)), fixed = c(2)),
                                  FixedTransition(population = 0,  rate = Rate(params[4]), fixed = c(0)))
  params = c(.3,.25,.1, .2)  
  res = branchTD(time, initial, transitionList, params, list(), dtime = .1, 1)
  plot(res[res$rep == 1,]$V3)
}

test_estim = function(){  
  # Specify how many units of time to simulate
  time = 20
  
  # Specify that the simulation will initiate with a single type with size 500
  initial = c(5000)
  
  params = c(.3,.25,.1, .2)  
  # Specify two fixed transitions, birth and death
  transitionList = TransitionList(FixedTransition(population = 0, rate = Rate(params[1] - params[2]*exp(-params[3]*t)), fixed = c(2)),
                                  FixedTransition(population = 0,  rate = Rate(params[4]), fixed = c(0)))
  simdat = formatSimData(branchTD(time, initial, transitionList, params, list(), dtime = .1, 100),1)
  est = estimateBP(time = simdat$time, N = simdat$t1_cells_prev, initTime = simdat$prev_time, data = simdat$t1_cells, transitionList = transitionList, lower = c(0,0,0,0), upper = c(.5,.5,.5,.5), initial = runif(4,0,.5))
}