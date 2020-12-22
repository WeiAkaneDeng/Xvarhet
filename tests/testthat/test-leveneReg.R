context("leveneReg")

testthat::test_that("errors", {
  testthat::expect_error(
    leveneReg(GENO=NULL, SEX=NA, Y=NA, COV=NULL)
  )

  testthat::expect_error(
    leveneReg(GENO=rbinom(100, 2, 0.3), SEX = rbinom(100, 1, 0.3), PLINK = 0, Y=rnorm(100), COV=rbinom(200, 2, 0.3))
  )
})
