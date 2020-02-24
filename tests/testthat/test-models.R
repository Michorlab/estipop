context("test the construction of objects in the pacakge, including rate, transition, process_model, stop_criteron, and stoplist")

test_that("incorrect rate objects are not instantiated", {
  expect_error(rate(params[1] + t[1]), "invalid expression t\\[1\\]")
  expect_error(rate(params[-1] + 2), "invalid expression params\\[-1\\]")
  expect_error(rate(c(1)), "invalid operator c!")
  expect_error(rate(exp(params[1],params[3])), "unary operator exp takes 1 argment!")
})

test_that("correct rate objects are instantiated", {
  r = rate(params[1]*t + log(t^5.0/4.0))
  expect_equal(class(r), "estipop_rate")
  expect_equal(r$exp, quote(params[1]*t + log(t^5.0/4.0)))
  expect_silent(rate(8.005*params[3] + 0.0001*log(t)))
  expect_silent(rate(sin(4*t)))
  expect_silent(rate(params[1] + params[2]*t + params[3]*t^2))
})

test_that("incorrect transition objects are not instantiated", {
  expect_error(transition("a","b","c"),"rate is not a valid estipop rate object!")
  expect_error(transition(rate(t*params[1]),"b","c"), "parent and offpspring must be numeric!")
  expect_error(transition(rate(t*params[1]),8,"c"), "parent and offpspring must be numeric!")
  expect_error(transition(rate(t*params[1]),0,c(1,2)), "parent must be positive and offspring must be nonnegative!")
  expect_error(transition(rate(t*params[1]),2,c(-1,2)), "parent must be positive and offspring must be nonnegative!")
  expect_error(transition(rate(t*params[1]),2,c(-1,2)), "parent must be positive and offspring must be nonnegative!")
  expect_error(transition(rate(t*params[1]),5,c(1,2)), "parent index is larger than length of offspring vector!")
  expect_error(transition(rate(t*params[1]),c(1,2),c(1,2)), "parent should not be a vector!")
})

test_that("correct transition objects are instantiated", {
  trans = transition(rate(sin(4*t)), 3, c(0,0,1,2,3))
  expect_equal(class(trans), "estipop_transition")
  expect_equal(trans$parent, 3)
  expect_equal(trans$offspring, c(0,0,1,2,3))
  expect_equal(trans$rate$exp, quote(sin(4*t)))
  expect_silent(transition(rate(8.005*params[3] + 0.0001*log(t)), 5, c(1,2,3,4,5)))
  expect_silent(transition(rate(sin(4*t)), 1, 1))
  expect_silent(transition(rate(params[1] + params[2]*t + params[3]*t^2), 2, c(0,0,0)))
})

test_that("incorrect process_model objects are not instantiated", {
  expect_error(process_model(1,2,3,4),"invalid transition object!")
  expect_error(process_model(transition(rate(t*params[1]), 1, c(1,2,3)),
                             transition(rate(t*params[2]), 3, c(1,2,3)),
                             "hello!"),"invalid transition object!")
  expect_error(process_model(),"transition_list must be a list with a positive number of elements!")
  expect_error(process_model(transition(rate(t*params[1]), 4, c(1,2,3))),"parent index is larger than length of offspring vector!")
  expect_error(process_model(transition(rate(t*params[1]), 1, c(1,2,3)),
                             transition(rate(t*params[2]), 3, c(1,2,3,4)),
                             transition(rate(t*params[1]), 2, c(1,2))),"transitions have inconsistent number of types!")
})

test_that("correct process_model objects are instantiated", {
  pm <- process_model(transition(rate(8.005*params[3] + 0.0001*log(t)), 2, c(1,2,3)),
                    transition(rate(sin(4*t)), 1, c(0,0,1)),
                    transition(rate(params[1] + params[2]*t + params[3]*t^2), 2, c(0,0,0)))
  expect_equal(pm$ntypes, 3)
  expect_equal(length(pm$transition_list), 3)
})

