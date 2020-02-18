context("verify that estimation functions throw errors on incorrect inputs")

test_that("estimation rejects incorrect inputs", {
  expect_error(estimate("a","b","c","d","e","f"), "model must be a process_model object!")
  model = process_model(transition(rate=rate(.5),parent=1,offspring=3), transition(rate = rate(.3), parent = 1, offspring = 0))
  expect_error(estimate(model,"b","c","d","e",c(0)),"all time and population inputs must be numeric!")
  expect_error(estimate(model, matrix(c(1,2),nrow=1), matrix(c(1,2),nrow=1) , c(4), 10, c(0)),"init_pop, final_pop, and model must all have the same number of types!")
  expect_error(estimate(model,1, 5, c(1,2,3,4), 10,c(0)), "init_pop, start_times, end_times, and final_pop must all have the same number of rows!")
  expect_error(estimate(model, c(3,3,1), c(1,4,5), c(0,-5,10),c(0)), "population and time variables must be nonnegative!")
  expect_error(estimate(model, c(3,3,-1), c(1,4,5), c(0,5,10),c(0)), "population and time variables must be nonnegative!")
  expect_error(estimate(model, c(3,3,1), c(1,4,-5), c(0,5,10),c(0)), "population and time variables must be nonnegative!")
})
