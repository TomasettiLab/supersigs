y <- c(rep(0, 5), rep(1, 5))
x1 <- 1:10
x2 <- 10:1
x3 <- c(1, 1, 5, 7, 9, 2, 4, 6, 8, 10)
x4 <- c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10)
x5 <- c(NA, 2:10)
x6 <- rep(NA, 10)
    
test_that("my_auc works", {
    expect_equal(my_auc(y, x1), 1)
    expect_equal(my_auc(y, x2), 0)
    expect_equal(my_auc(y, x3), 0.64)
    expect_equal(my_auc(y, x4), 0.6)
})

test_that("my_auc handles missing values", {
    expect_equal(my_auc(y, x5), 1)
    expect_equal(my_auc(y, x6), NA)
})
