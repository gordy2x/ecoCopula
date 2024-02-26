test_that("plotting successfully", {
  skip_on_ci()
  X <- as.data.frame(spider$x)
  abund <- spider$abund 
  
  pa <-(abund>0)*1
  
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  residual_path <- save_png(plot(spider_pa))
  
  expect_snapshot_file(residual_path, "spider_residuals.png")
  
  set.seed(5)
  spid_lv <- cord(spider_pa) 
  
  path <- save_png(plot(spid_lv,biplot = TRUE))
  
  expect_snapshot_file(path, "spider_biplot.png")
})
