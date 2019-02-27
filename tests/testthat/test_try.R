context("my first test file")

library(saint)
data(spider)
spider.mod=manyglm(mvabund(spider$abund)~spider$x)

test_that("error", {
  expect_error(saint(4))
})

test_that("obj file type", {
  expect_is(saint(spider.mod), "saint")
})