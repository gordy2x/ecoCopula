test_that("output of stackedsdm", {
  X <- as.data.frame(spider$x)
  abund <- spider$abund 
  
  pa <-(abund>0)*1
  
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  expect_s3_class(spider_pa, "stackedsdm")
  expect_visible(spider_pa)
  expect_named(spider_pa)
  expect_equal(length(spider_pa), 9)
})
