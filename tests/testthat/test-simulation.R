#Simulation tests. Not part of the unit test suite because they take a while, but run these to make sure simulation is accurate

#test that simulation results match moments -- one type cae
time = 20
initial = c(500)

params = c(.3,.25,.1, .2)  
model = process_model(transition(rate = rate(params[1] - params[2]*exp(-params[3]*t)), parent = 1, offspring = 2),
                      transition(rate = rate(params[4]), parent = 1, offspring = 0))
res = branch(model, params, initial, time, 1000)
mom = compute_mu_sigma(model, params, 0, 20, initial)



  #test that simulation results match moments -- two type case
  
  #test that simulation results match moments -- three type case
  
  #test that simulation rejects invalid inputs
