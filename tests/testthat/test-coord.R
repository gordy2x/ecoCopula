test_that("Output", {
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  set.seed(5)
  spid_lv <- cord(spider_pa) 
  
  expect_s3_class(spid_lv, "cord")
  expect_visible(spid_lv) 
  expect_named(spid_lv) 
})


test_that("print and summary", {
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  expect_visible(spider_pa)
  expect_visible(summary(spider_pa))
})
