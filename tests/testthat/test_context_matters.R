x1 <- tibble(x = 1)

test_that("context_matters valid input", {
  expect_error(context_matters(muts_df = x1),
               "Column names of input data are not correct 
              in context_matters")
})