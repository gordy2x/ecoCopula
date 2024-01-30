test_that("output of stackedsdm + residuals", {
  set.seed(5)
  spider_pa <- stackedsdm(pa,~1, data = X, family="binomial",ncores = 2)
  
  expect_s3_class(spider_pa, "stackedsdm")
  expect_visible(spider_pa)
  expect_named(spider_pa)
  expect_equal(length(spider_pa), 9)
  
  expect_visible(residuals(spider_pa))
  expect_equal(residuals(spider_pa) |> nrow(), 28)
  expect_equal(spider_pa$family |> unique(), "binomial")
  
  expect_visible(summary(spider_pa))
  expect_visible(predict(spider_pa, type = "response"))
})

test_that("Family argument = neg.binomial",{
  set.seed(5)
  #marginal model
  spider_nb <- stackedsdm(abund,~., data = X, family="negative.binomial", ncores = 2) #eqiv. manyglm()
  
  expect_s3_class(spider_nb, "stackedsdm")
  expect_visible(spider_nb)
  expect_named(spider_nb)
  expect_equal(spider_nb$family |> unique(), "negative.binomial")
}
)


test_that("Family argument = ordinal",{
  set.seed(5)
  #marginal model
  bryce_marg <- stackedsdm(bryceord, formula_X = ~ 1, data = brycesite, family="ordinal",ncores = 2)
  
  expect_s3_class(bryce_marg, "stackedsdm")
  expect_visible(bryce_marg)
  expect_named(bryce_marg)
  expect_equal(bryce_marg$family |> unique(), "ordinal")
}
)

