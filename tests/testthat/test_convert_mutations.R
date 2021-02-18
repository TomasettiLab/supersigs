x <- "[C>A]T"
x_div <- c("A[C>A]T", "C[C>A]T", "G[C>A]T", "T[C>A]T")
y <- "A[C>A]T"

test_that("condense_mutations works", {
    expect_equal(condense_mutations(x_div), x)
    expect_equal(y, y)
})

test_that("convert_to_level3 works", {
    expect_equal(convert_to_level3(x), x_div)
    expect_equal(convert_to_level3(y), y)
    expect_equal(convert_to_level3(c(y, x)),
                             c(y, x_div))
})
