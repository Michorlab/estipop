context("Test that we are computing process moments correctly")

test_that("moments work in cases where we have analytical solutions", {
  
  #simple birth-death model
  process = process_model(transition(rate(params[1]), 1, 2),
                          transition(rate(params[2]), 1, 0))
  init_pop = 50
  b = 3
  d = 2
  t0 = 5
  tf = 10
  mom = compute_mu_sigma(process, c(b,d), t0, tf, init_pop)
  mu_real = exp((b-d)*(tf - t0))*init_pop
  sigma_real = (exp((tf - t0)*(b-d))*(b*(2*exp((tf - t0)*(b-d))-1)-d)/(b-d) - exp(2*(tf - t0)*(b-d)))*init_pop; #analytical variance solution
  
  #verify that we have < .001% disagreement
  expect_lt(abs(mom$mu - mu_real)/mu_real, .00001)
  expect_lt(abs(mom$Sigma - sigma_real)/sigma_real, .00001)
  
  #inhomogenous birth-death model
  process = process_model(transition(rate(params[1]*t), 1, 2),
                          transition(rate(params[2]*t), 1, 0))
  init_pop = 50
  b = 3
  d = 2
  t0 = 0
  tf = 3
  mom = compute_mu_sigma(process, c(b,d), t0, tf, init_pop)
  mu_real = exp(.5*(tf^2 - t0^2))*init_pop #integral formula for the mean

  #verify that we have < .001% disagreement
  expect_lt(abs(mom$mu - mu_real)/mu_real, .00001)
})
