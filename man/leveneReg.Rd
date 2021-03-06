% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LeveneReg.R
\name{leveneReg}
\alias{leveneReg}
\title{Levene's regression tests for variance homogeneity by SNP genotype}
\usage{
leveneReg(GENO, SEX, PLINK = FALSE, Y, COV = NULL, test_type = "REC")
}
\arguments{
\item{GENO}{the genotype of a SNP, must be a vector of 0, 1, 2's indicating the number of
reference alleles. The length of \code{GENO} should match that of \code{SEX}, \code{Y},
and any covariates \code{COV}.}

\item{SEX}{the genetic sex of individuals in the sample population, must be a vector of 1 and 2 or 0 and 1, depending on whether the PLINK sex code is used. Note that the default sex code is 1 for males and 2 for females in PLINK.}

\item{PLINK}{a logical indicating whether the SEX is coded following PLINK, i.e.
female coded as 2 and male colded as 1; if not, SEX is assumed to be coded 0 for females
and 1 for males.}

\item{Y}{a vector of quantitative traits, such as human height.}

\item{COV}{a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.}

\item{test_type}{a character of either "ALL", printing all strategies, or "M1",
"M2", or "M3" for the 8 stratigies given the mean stage models, or "REC" printing
only the recommended tests.}
}
\value{
a vector of Levene's test regression p-values according to the models
specified.
}
\description{
The function takes as input the genotype of a SNP (\code{GENO}), the SEX (\code{SEX}), and
 a quantitative trait (\code{Y}) in a sample population, and possibly additional
 covariates. It should be noted that these variables must be of the same length.
 The function then returns the variance heterogeneity \emph{p}-values using the
 generalized Levene's test. The residual function could alternatively be replaced with the quantile regression quantreg::rq following the "fn" algorithm, for more details see \code{?quantreg::rq}.
}
\note{
We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
the variance values tend to be proportional to mean values when stratified by \code{GENO}.
}
\examples{
N <- 5000
GENO <- rbinom(N, 2, 0.3)
sex <- rbinom(N, 1, 0.5)
Y <- rnorm(N)
cov <- matrix(rnorm(N*10), ncol=10)
leveneReg(GENO=GENO, SEX=sex, Y=Y, COV=cov, test_type="ALL")
leveneReg(GENO=GENO, SEX=sex, Y=Y, COV=cov, test_type="REC")

}
\author{
Wei Q. Deng \email{deng@utstat.toronto.edu}
}
