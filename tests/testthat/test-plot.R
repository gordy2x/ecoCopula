save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  set.seed(5)
  png(path, width = width, height = height)
  on.exit(dev.off())
  code
  
  path
}


test_that("plotting successfully", {
  X <- as.data.frame(spider$x)
  abund <- spider$abund 
  
  pa <-(abund>0)*1
  
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  set.seed(5)
  spid_lv <- cord(spider_pa) 
  
  path <- save_png(plot(spid_lv,biplot = TRUE))
  
  expect_snapshot_file(path, "spider_biplot.png")
})
