source("fixtures.R")

context("limited support for stackedsdm")

test_that("error if multiple distributions", {

  expect_error(
    simulate(ssdm_multi.cord)
  )

})

context("class of output")

test_that("matrix output", {

  returned = simulate(nb0.cord)
  expect_is(returned, "matrix")

})

context("manylm")

test_that("continuous and matrix size", {

  returned = simulate(fit1.cord, nsim=2)
  expect_false(any(returned%%1 == 0))
  expect_equal(nrow(returned), 56)
  expect_equal(colnames(returned), colnames(abund))

})

context("presence absence")

test_that("ones and zeros", {

  returned = simulate(bin0.cord)
  expect_true(all(returned == 0 | returned == 1))

})

test_that("relative change in proportion of ones", {

  returned = simulate(bin0.cord)
  expect_identical(
    round(abs(mean(returned)-mean(pa))/mean(pa)), 0
  )

})

context("manyglm - various family and specs")

test_that("discrete and matrix size", {

  for (fit in list(
    nb0.cord, nb1.cord, bin0.cord,
    poi1.cord, nb2.cord, nb_mix.cord,
    nb_fac2.cord, nb_fac4.cord, nb_mth.cord
  )) {

    returned = simulate(fit)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 28)
    expect_equal(colnames(returned), colnames(abund))

  }

})

context("stackedsdm - various specs")

test_that("discrete and matrix size", {

  for (fit in list(
    nb0.ssdm.cord, nb1.ssdm.cord, nb_fac4.ssdm.cord
  )) {

    returned = simulate(fit, nsim=2)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 56)
    expect_equal(colnames(returned), colnames(abund))
  }

})

context("newdata - new factor")

test_that("discrete and matrix size", {

  for (fit in list(
    nb_fac2.cord, nb_mth.cord
  )) {

    returned = simulate(fit, nsim=3, newdata=Xnew)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 84)
    expect_equal(colnames(returned), colnames(abund))
  }

})

context("newdata - subset")

test_that("discrete and matrix size", {

  for (fit in list(
    poi1.cord, nb2.cord, nb_mix.cord, nb_mth.cord,
    nb0.ssdm.cord, nb1.ssdm.cord, nb_fac4.ssdm.cord
  )) {

    returned = simulate(fit, newdata=Xnew_sub)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 10)
    expect_equal(colnames(returned), colnames(abund))
  }

})

context("newdata - diff size")

test_that("discrete and matrix size", {

  for (fit in list(
    nb_fac4.cord, nb_fac4.ssdm.cord
  )) {

    returned = simulate(fit, nsim=2, newdata=Xvec)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 30)
    expect_equal(colnames(returned), colnames(abund))
  }

})
