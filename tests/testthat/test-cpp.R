context("C++")
test_that("Catch unit tests pass", {
		skip_on_cran()
		skip_on_travis()
		skip_on_appveyor()
    expect_cpp_tests_pass("largeVis")
})
