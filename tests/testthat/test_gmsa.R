context("generate_min_sima_algebra")

x1 <- list(A = c("A[C>A]T"))
x2 <- list(c("A[C>A]T"), c("A[C>A]T"))
x3 <- list(A = c("A[C>A]T"), B = c("A[C>A]T"))

test_that("gmsa valid input", {
  expect_error(generate_min_sigma_algebra(x1),
               "input_ls has length not equal to 2")
  expect_error(generate_min_sigma_algebra(x2), 
               "input_ls is not named")
})

input_ls = list(Exposed = c("[C>A]C", "T[T>G]T"), 
                Unexposed = c("C>A", "T[T>G]T"))
ref <- list(F1 = c("[C>A]A", "[C>A]G", "[C>A]T"),
            F2 = c("[C>A]C"),
            F3 = c("T[T>G]T"),
            Remaining = c("C>G", "C>T", "T>A", "T>C", "[T>G]A",
                          "[T>G]C", "[T>G]G", "A[T>G]T", "C[T>G]T",
                          "G[T>G]T"))
test_that("gmsa works", {
  out <- generate_min_sigma_algebra(input_ls, condense = T)
  expect_equal(out, ref)
})
