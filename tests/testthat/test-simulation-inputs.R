context("verify that simulation functions throw errors on incorrect inputs")

test_that("exact simulation rejects incorrect inputs", {
  expect_error(branch("a","b","c","d","e"), "model must be a process_model object!")
  model = process_model(transition(rate=rate(.5),parent=1,offspring=3), transition(rate = rate(.3), parent = 1, offspring = 0))
  expect_error(branch(model,"b","c","d","e"),"all time, population, and parameter inputs must be numeric!")
  expect_error(branch(model,NULL, c(1,2),c(1,2,3,4),10),"init_pop and model must have same number of types ")
  expect_error(branch(model,NULL, -1,c(1,2,3,4),10),"population must be nonnegative and reps must be positive!")
  expect_error(branch(model,NULL, 1,c(1,2,3,-1),10), "all observation times must be nonnegative.")
  expect_error(branch(model,NULL, 1,c(1,2,3,5),-10), "population must be nonnegative and reps must be positive!")
  expect_error(branch(model,"c", 1,c(1,2,3,5),-10), "all time, population, and parameter inputs must be numeric!")
  
  model = process_model(transition(rate = rate(.3), parent = 1, offspring = c(2,0)),
                        transition(rate = rate(.2), parent = 1, offspring = c(0,0)),
                        transition(rate = rate(.25), parent = 2, offspring = c(0,1)),
                        transition(rate = rate(.15), parent = 2, offspring = c(0,0)),
                        transition(rate = rate(.05), parent = 1, offspring = c(1,1)))
  
  expect_error(branch(model,NULL, 1,c(1,2,3,4),10),"init_pop and model must have same number of types ")
})

test_that("approximate simulation rejects incorrect inputs", {
  expect_error(branch_approx("a","b","c","d","e"), "model must be a process_model object!")
  model = process_model(transition(rate=rate(.5),parent=1,offspring=3), transition(rate = rate(.3), parent = 1, offspring = 0))
  expect_error(branch_approx(model,"b","c","d","e"),"all time, population, and parameter inputs must be numeric!")
  expect_error(branch_approx(model,NULL, c(1,2),c(1,2,3,4),10),"init_pop and model must have same number of types ")
  expect_error(branch_approx(model,NULL, -1,c(1,2,3,4),10),"population must be nonnegative and reps must be positive!")
  expect_error(branch_approx(model,NULL, 1,c(1,2,3,-1),10), "all observation times must be nonnegative.")
  expect_error(branch_approx(model,NULL, 1,c(1,2,3,5),-10), "population must be nonnegative and reps must be positive!")
  expect_error(branch_approx(model,"c", 1,c(1,2,3,5),-10), "all time, population, and parameter inputs must be numeric!")
  
  model = process_model(transition(rate = rate(.3), parent = 1, offspring = c(2,0)),
                        transition(rate = rate(.2), parent = 1, offspring = c(0,0)),
                        transition(rate = rate(.25), parent = 2, offspring = c(0,1)),
                        transition(rate = rate(.15), parent = 2, offspring = c(0,0)),
                        transition(rate = rate(.05), parent = 1, offspring = c(1,1)))
  
  expect_error(branch_approx(model,NULL, 1,c(1,2,3,4),10),"init_pop and model must have same number of types ")
})
