test_that("Check error catching", {
  expect_error(PCA(c(12,11,2,3,4,5,76,45,4,4,)))
  expect_error(PCA(matrix(c("12","11","2","3","4","5","76","45","4"), nrow=3, ncol=3)))
})
