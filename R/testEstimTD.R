
  #plot(res[res$rep == 1,]$V3)

test_estim = function(){  
  # Specify how many units of time to simulate
  time = 20
  
  # Specify that the simulation will initiate with a single type with size 500
  initial = c(5000)
  
  params = c(.3,.25,.1, .2)  
  # Specify two fixed transitions, birth and death
  model = process_model(transition(rate = rate(params[1] - params[2]*exp(-params[3]*t)), parent = 1, offspring = 2),
                        transition(rate = rate(params[4]), parent = 1, offspring = 0))
  res = branch(model, params, initial, seq(0,20,.1), 1)
  simdat = formatSimData(branchTD(time, initial, transitionList, params, list(), dtime = .1, 20),1)
  compute_mu_sigma(model, params, 0, 20, initial)
  est = estimateBP_timedep(time = simdat$time, N = simdat$t1_cells_prev, initTime = simdat$prev_time, data = simdat$t1_cells, transitionList = transitionList, lower = c(0,0,0,0), upper = c(.5,.5,.5,.5), initial = runif(4,0,.5))
}