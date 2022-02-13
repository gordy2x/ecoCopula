context("my first test file")

library(ecoCopula)
library(mvabund)
data(spider)
spider.mod=manyglm(mvabund(spider$abund)~.,data=spider$x)

test_that("error", {
  expect_error(cgr(4))
})

test_that("obj file type", {
  expect_is(cgr(spider.mod), "cgr")
})
