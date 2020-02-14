#Estimation tests. Not part of the unit test suite because they take a while, but run these to make sure estimtaion is accurate

#test that estimation results are accurate -- one type inhomogenous case
time = seq(0,20,.5)
initial = c(500)

params = c(.3,.25,.1, .2)  
model = process_model(transition(rate = rate(params[1] - params[2]*exp(-params[3]*t)), parent = 1, offspring = 2),
                      transition(rate = rate(params[4]), parent = 1, offspring = 0))
res = branch(model, params, initial, time, 20)
simdata = format_sim_data(res, model$ntypes)
estimate_td(model, init_pop = simdata$type1_prev,  start_times = simdata$prev_time, final_pop =simdata$type1, end_times = simdata$time, initial = runif(4,0,.5), lower = rep(0,4), upper = rep(.5,4))

#test that estimation results are correct - two type homogenous

time = seq(0,5,1) 
initial = c(100, 0)

params = c(.4,.1,.7,.1, .3)

model = process_model(transition(rate = rate(params[1]), parent = 1, offspring = c(2,0)),
                      transition(rate = rate(params[2]), parent = 1, offspring = c(0,0)),
                      transition(rate = rate(params[3]), parent = 2, offspring = c(0,2)),
                      transition(rate = rate(params[4]), parent = 2, offspring = c(0,0)),
                      transition(rate = rate(params[5]), parent = 1, offspring = c(1,1)))



res = branch(model, params, initial, time, 20)
simdata = format_sim_data(res, model$ntypes)
estimate(model, init_pop = cbind(simdata$type1_prev, simdata$type2_prev),  final_pop = cbind(simdata$type1, simdata$type2), times = simdata$dtime, initial = runif(5,0,.5), lower = rep(0,5), upper = rep(1,5))
