context("sgd batches")

test_that("sgd batches is linear with E", {
	expect_equal(sgdBatches(1e6, 10), sgdBatches(1e6, 1e4))
})

test_that("sgd batches has a small bump", {
	expect_lt(sgdBatches(5000), sgdBatches(10001))
	expect_lt(sgdBatches(9999) / 1.5, sgdBatches(10001))
})
