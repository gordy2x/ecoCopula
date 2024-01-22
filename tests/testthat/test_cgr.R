test_that("testing cgr so it error", {
  expect_error(cgr(4))
})

test_that("obj file type", {
  spider.mod=manyglm(mvabund(spider$abund)~.,data=as.data.frame(spider$x))
  expect_s3_class(cgr(spider.mod), "cgr")

  })
