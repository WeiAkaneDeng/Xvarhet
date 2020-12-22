context("leveneTests")

testthat::test_that("errors", {
  testthat::expect_error(
   leveneTests(GENO=NULL, SEX=NA, Y=NA, COV=NULL)
  )

  testthat::expect_error(
    leveneTests(GENO=NULL, Y=NA, COV=NULL)
  )
  })
