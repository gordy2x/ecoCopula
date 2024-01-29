test_that("Throws error if multiple distributions for stackedsdm", {

  expect_error(
    stats::simulate(ssdm_multi.cord)
  )

})

test_that("class of output", {

  returned = stats::simulate(nb0.cord)
  expect_equal(class(returned), c("matrix", "array"))

})

test_that("manylm continuous and matrix size", {

  returned = stats::simulate(fit1.cord, nsim=2)
  expect_false(any(returned%%1 == 0))
  expect_equal(nrow(returned), 56)
  expect_equal(colnames(returned), colnames(abund))

})


test_that("presence absence/ones and zeros", {

  returned = stats::simulate(bin0.cord)
  expect_true(all(returned == 0 | returned == 1))

})

test_that("relative change in proportion of ones", {

  returned = stats::simulate(bin0.cord)
  expect_identical(
    round(abs(mean(returned)-mean(pa))/mean(pa)), 0
  )

})


test_that("manyglm - various family and specs, discrete and matrix size", {

  for (fit in list(
    nb0.cord, nb1.cord, bin0.cord,
    poi1.cord, nb2.cord, nb_mix.cord,
    nb_fac2.cord, nb_fac4.cord, nb_mth.cord
  )) {

    returned = stats::simulate(fit)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 28)
    expect_equal(colnames(returned), colnames(abund))

  }

})

test_that("stackedsdm - various specs discrete and matrix size", {

  for (fit in list(
    nb0.ssdm.cord, nb1.ssdm.cord, nb_fac4.ssdm.cord
  )) {

    returned = stats::simulate(fit, nsim=2)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 56)
    expect_equal(colnames(returned), colnames(abund))
  }

})


test_that("newdata - new factor discrete and matrix size", {

  for (fit in list(
    nb_fac2.cord, nb_mth.cord
  )) {

    returned = stats::simulate(fit, nsim=3, newdata=Xnew)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 84)
    expect_equal(colnames(returned), colnames(abund))
  }

})

test_that("newdata - subset, discrete and matrix size", {

  for (fit in list(
    poi1.cord, nb2.cord, nb_mix.cord, nb_mth.cord,
    nb0.ssdm.cord, nb1.ssdm.cord, nb_fac4.ssdm.cord
  )) {

    returned = stats::simulate(fit, newdata=Xnew_sub)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 10)
    expect_equal(colnames(returned), colnames(abund))
  }

})


test_that("newdata - diff size, discrete and matrix size", {

  for (fit in list(
    nb_fac4.cord, nb_fac4.ssdm.cord
  )) {

    returned = stats::simulate(fit, nsim=2, newdata=Xvec)
    expect_true(all(returned%%1 == 0))
    expect_equal(nrow(returned), 30)
    expect_equal(colnames(returned), colnames(abund))
  }

})
