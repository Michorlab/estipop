#Simulation tests. Not part of the unit test suite because they take a while, but run these to make sure simulation is accurate

#test that simulation results match moments -- one type inhomogenous case
time = 20
initial = c(500)

params = c(.3,.25,.1, .2)  
model = process_model(transition(rate = rate(params[1] - params[2]*exp(-params[3]*t)), parent = 1, offspring = 2),
                      transition(rate = rate(params[4]), parent = 1, offspring = 0))
res = branch(model, params, initial, time, 1000)
mom = compute_mu_sigma(model, params, 0, 20, initial)
print(paste("Sample Mean: ", mean(res$V3)))
print(paste("True Mean: ", mom$mu))
print(paste("Sample Variance: ", var(res$V3)))
print(paste("True Variance: ", mom$Sigma))

#test that simulation results match moments -- two type inhomogenous case
  
time = 20
initial = c(500, 300)

model = process_model(transition(rate = rate(.6*exp(-.1*t)), parent = 1, offspring = c(2,0)),
                      transition(rate = rate(.4*exp(-.1*t)), parent = 1, offspring = c(0,0)),
                      transition(rate = rate(.3*exp(-.1*t)), parent = 2, offspring = c(0,1)),
                      transition(rate = rate(.2*exp(-.1*t)), parent = 2, offspring = c(0,0)),
                      transition(rate = rate(.1*exp(-.1*t)), parent = 1, offspring = c(1,1)))
                      

res = branch(model, params, initial, time, 1000)
mom = compute_mu_sigma(model, params, 0, 20, initial)
print(paste("Sample Mean: ", c(mean(res$V3), mean(res$V4))))
print(paste("True Mean: ", mom$mu))
print(paste("Sample Variance: ", c(var(res$V3), var(res$V4))))
print(paste("Sample Covariance:", cov(res$V3, res$V4)))
print(paste("True Variance: ", diag(mom$Sigma)))
print(paste("True Covariance: ", mom$Sigma[1,2]))

#test that simulation results match moments -- two type homogenous case
  
#test that simulation rejects invalid inputs