x1 <- tibble(x = 1,
             age = 2,
             indvar = FALSE)
x2 <- tibble(sample_id = 1,
             x = 2,
             indvar = FALSE)
x3 <- tibble(sample_id = 1,
             age = 2,
             x = FALSE)

test_that("get_signature valid input", {
  expect_error(get_signature(x1, factor = "f"),
               "Input data frame missing sample_id column")
  expect_error(get_signature(x2, factor = "f"),
               "Input data frame missing AGE column")
  expect_error(get_signature(x3, factor = "f"),
               "Input data frame missing IndVar column")
  
})