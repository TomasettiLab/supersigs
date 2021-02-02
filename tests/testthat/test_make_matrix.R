test_that("make_matrix valid output", {
  out <- make_matrix(example_dt)
  expect_equal(nrow(out), length(unique(example_dt$sample_id)))
  expect_equal(ncol(out), 98)
})