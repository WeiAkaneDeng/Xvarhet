#' Levene's test for variance homogeneity by SNP genotypes
#'
#' The function takes as input the genotypes of a SNP (\code{GENO}), the sex (\code{SEX}), and a quantitative
#'    trait (\code{Y}) in a sample population, and possibly additional covariates. It should be noted
#'    that these variables must be of the same length. The function then
#'    returns the variance heterogeneity \emph{p}-values for the model \eqn{Y\sim G},
#'    \eqn{Y \sim G\times S}, the sex-specific results based on model \eqn{Y\sim G}, as well as that
#'    using Fisher's method to combine sex-specific results.
#'
#'
#' @param GENO the genotype of a SNP, must be a vector of 0, 1, 2's indicating the number of reference alleles. The length of \code{GENO} should match that of \code{SEX}, \code{Y}, and any covariates.
#' @param SEX the genetic sex of individuals in the sample population, must be a vector of 1 and 2 or 0 and 1, depending on whether the PLINK sex code is used. Note that the default sex code is 1 for male and 2 for female in PLINK.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COV a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param PLINK a logical indicating whether the SEX is coded following PLINK or not.
#' @param centre a character indicating whether the absolute deviation should be calculated with respect to ``median'' or ``mean''.
#'
#' @return a vector of Levene's test p-values according to levels specified by \code{GENO}, the
#' interaction between \code{GENO} and \code{SEX}, sex-specific Levene's test stratified by
#' \code{GENO}, and the Fisher's method to combine the sex-specific Levene's test \emph{p}-values.
#'
#' @importFrom stats pchisq
#' @importFrom stats lm
#' @importFrom stats resid
#' @importFrom stats complete.cases
#' @importFrom stats anova
#'
#' @examples
#' N <- 5000
#' geno <- rbinom(N, 2, 0.3)
#' sex <- rbinom(N, 1, 0.5)
#' y <- rnorm(N)
#' cov <- matrix(rnorm(N*10), ncol=10)
#' leveneTests(GENO=geno, SEX=sex, Y=y, COV=cov)
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}
#'
#' @export leveneTests
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{G}.
#'

leveneTests <- function(GENO, SEX, PLINK = FALSE, Y, centre = "median", COV = NULL) {

  if (missing(GENO))
    stop("The genotype input is missing.")

  if (missing(SEX))
    stop("The sex input is missing.")

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (length(GENO)!=length(SEX)|length(GENO)!=length(Y)|length(Y)!=length(SEX))
    stop("Make sure the inputs have the same length.")


  if (!is.null(COV)){
    Y <- stats::resid(stats::lm(Y ~ COV, na.action = na.exclude))
  }

  if (length(table(GENO))==1) {

    warning("Monomorphic SNP detected");
    PVAL <- c(NA, NA, NA, NA, NA)
    names(PVAL) <- c("Lev3", "Lev5", "LevF", "LevM", "Fisher")

  } else {
  N <- length(GENO)
  GENO <- factor(GENO)
  group <- factor(interaction(GENO, SEX))
  meds <- tapply(Y, group, centre, na.rm = TRUE)
  meds3 <- tapply(Y, GENO, centre, na.rm = TRUE)
  resp <- abs(Y - meds[group])
  resp3 <- abs(Y - meds3[GENO])
  Lp_3G <- anova(lm(resp3 ~ GENO))[, c(1, 4, 5)][1, 3]
  Lp_5G <- anova(lm(resp ~ group))[, c(1, 4, 5)][1, 3]

  if (PLINK) {

  if (length(table(SEX))==1) {

    if (sum(SEX==1) == N){
      warning("Only Males detected")
      Lp_3G_F <- NA
      Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[, c(1, 4, 5)][1, 3]
      }else{
      warning("Only Females detected")
      Lp_2G_M <- NA
      Lp_3G_F <- anova(lm(resp[SEX == 2] ~ factor(GENO[SEX == 2])))[, c(1, 4, 5)][1, 3]
      PVAL <- c(NA, NA, NA, NA, NA)
      }

  } else {

  Lp_3G_F <- anova(lm(resp[SEX == 2] ~ factor(GENO[SEX == 2])))[, c(1, 4, 5)][1, 3]
  Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[, c(1, 4, 5)][1, 3]
  df_sum <- (anova(lm(resp[SEX == 2] ~ factor(GENO[SEX == 2])))[1, 1] +
               anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[1, 1])

  df_full <- ifelse(sum(table(interaction(GENO, SEX))>0) < df_sum, df_sum + 2,
                    sum(table(interaction(GENO, SEX))>0))
  Lp_Fisher5 <- pchisq(-2 * (log(Lp_3G_F) + log(Lp_2G_M)),
                       df = ifelse(df_sum < 3, df_sum + 1, df_full-1), lower.tail = FALSE)

  }
  }
    else {

      if (length(table(SEX))==1) {

        if (sum(SEX==1) == N){
          warning("Only Males detected")
          Lp_3G_F <- NA
          Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[, c(1, 4, 5)][1, 3]
        } else if (sum(SEX==1) == 0) {
          warning("Only Females detected")
          Lp_2G_M <- NA
          Lp_3G_F <- anova(lm(resp[SEX == 0] ~ factor(GENO[SEX == 0])))[, c(1, 4, 5)][1, 3]
          PVAL <- c(NA, NA, NA, NA, NA)
        }

      } else {

        Lp_3G_F <- anova(lm(resp[SEX == 0] ~ factor(GENO[SEX == 0])))[, c(1, 4, 5)][1, 3]
        Lp_2G_M <- anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[, c(1, 4, 5)][1, 3]
        df_sum <- (anova(lm(resp[SEX == 0] ~ factor(GENO[SEX == 0])))[1, 1] +
                     anova(lm(resp[SEX == 1] ~ factor(GENO[SEX == 1])))[1, 1])

        df_full <- ifelse(sum(table(interaction(GENO, SEX))>0) < df_sum, df_sum + 2,
                          sum(table(interaction(GENO, SEX))>0))
        Lp_Fisher5 <- pchisq(-2 * (log(Lp_3G_F) + log(Lp_2G_M)),
                             df = ifelse(df_sum < 3, df_sum + 1, df_full-1), lower.tail = FALSE)

      }
    }


PVAL <- c(Lp_3G, Lp_5G, Lp_3G_F, Lp_2G_M, Lp_Fisher5)
names(PVAL) <- c("Lev3", "Lev5", "LevF", "LevM", "Fisher")
return(PVAL)
}
}
