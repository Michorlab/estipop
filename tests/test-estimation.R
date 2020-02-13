#Estimation tests. Not part of the unit test suite because they take a while, but run these to make sure estimtaion is accurate

#test that estimation results are accurate -- one type inhomogenous case
time = seq(0,20,.5)
initial = c(500)

params = c(.3,.25,.1, .2)  
model = process_model(transition(rate = rate(params[1] - params[2]*exp(-params[3]*t)), parent = 1, offspring = 2),
                      transition(rate = rate(params[4]), parent = 1, offspring = 0))
res = branch(model, params, initial, time, 10)
simdata = format_sim_data(res, model$ntypes)
estimate_td(model, init_pop = simdata$type1_prev,  start_times = simdata$prev_time, final_pop =simdata$type1, end_times = simdata$time, initial = runif(4,0,.5), lower = rep(0,4), upper = rep(.5,4))
