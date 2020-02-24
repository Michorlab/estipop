context("test parsing of R expressions into C++ code")

test_that("generate_cpp flags expressions with incorrect numbers of parameters",{
  expect_error(generate_cpp(quote(params[10]), 2), "Parameter params\\[10\\] goes beyond the number of parameters provided.")
  expect_error(generate_cpp(quote(params[7]), c(1,2,3,4,5)), "Parameter params\\[7\\] goes beyond the number of parameters provided.")
  expect_error(generate_cpp(quote(params[0]), 1), "Parameter params\\[0\\] goes beyond the number of parameters provided.")
})

test_that("Valid expressions are correctly translated to C++ by generate_cpp",{
  expect_equal(generate_cpp(quote(params[1]), 5), "5.000000")
  expect_equal(generate_cpp(quote(t*5), c(1,2,3)), "t * 5")
  expect_equal(generate_cpp(quote(exp(params[3]) + log(t)), c(1,2,17)), "exp(17.000000) + log(t)")
  expect_equal(generate_cpp(quote(log(t^4)*(params[5]*(sin(t*params[4])*params[1]))),c(0,.5,1.0,1.75,3.6)), "log(t ^ 4) * (3.600000 * (sin(t * 1.750000) * 0.000000))")
  expect_equal(generate_cpp(quote(params[1]*exp(t*exp(1)*params[2])), c(5,6)), "5.000000 * exp(t * exp(1) * 6.000000)")
})

test_that("Inavlid expressions are detected by the validator",{
  expect_error(check_valid(quote(params[1] + t[1])), "invalid expression t\\[1\\]")
  expect_error(check_valid(quote(t1)), "invalid name t1")
  expect_error(check_valid(quote(params[-1] + 2)), "invalid expression params\\[-1\\]")
  expect_error(check_valid(quote(t*5 + params[-6])), "invalid expression params\\[-6\\]")
  expect_error(check_valid(quote(params[abcd])), "invalid expression params\\[abcd\\]")
  expect_error(check_valid(quote(params[])), "invalid expression params\\[\\]")
  expect_error(check_valid(quote(params[t])), "invalid expression params\\[t\\]")
  expect_error(check_valid(quote(params[1]*t[5])), "invalid expression t\\[5\\]")
  expect_error(check_valid(quote(params[1+1])), "invalid expression params\\[1 \\+ 1\\]")
  expect_error(check_valid(quote(params[3] + params12)), "invalid name params12")
  expect_error(check_valid(quote(params[[3]])), "invalid operator \\[\\[!")
  expect_error(check_valid(quote("hi mom")), "invalid expression \"hi mom\"")
  expect_error(check_valid(quote(params[2] & t)), "invalid operator &!")
  expect_error(check_valid(quote(c(1))), "invalid operator c!")
  expect_error(check_valid(quote(exp(params[1],params[3]))), "unary operator exp takes 1 argment!")
  expect_error(check_valid(quote(log(params[1],params[3], params[2]))), "unary operator log takes 1 argment!")
  expect_error(check_valid(quote(sin())), "unary operator sin takes 1 argment!")
})


test_that("is_const correctly recognizes whether expressions depend on time",{
  expect_equal(is_const(quote(params[1])), T)
  expect_equal(is_const(quote(params[1]*t)), F)
  expect_equal(is_const(quote(t)), F)
  expect_equal(is_const(quote(sin(4*t))), F)
  expect_equal(is_const(quote(4*params[8])), T)
  expect_equal(is_const(quote(8.005*params[3] + 0.0001*log(t))), F)
  expect_equal(is_const(quote(2*params[1] + 4*params[2]^2 + 8*params[3]^3)), T)
  expect_equal(is_const(quote(params[1] + params[2]*t + params[3]*t^2)), F)
})

