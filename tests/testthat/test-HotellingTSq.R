test_that("Check error catching", {
  expect_error(HotellingTSq(c(12,11,2,3,4,5,76,45,4,4,), c(0,0,0)))
  expect_error(HotellingTSq(matrix(c("-1","2","3","4","5","6","7","8","9"), nrow=3, ncol=3, byrow=TRUE), c(0,0,0)))
})
