#' Levene's regression tests for variance homogeneity by SNP genotype
#'
#'  The function takes as input the genotype of a SNP (\code{GENO}), the SEX (\code{SEX}), and
#'  a quantitative trait (\code{Y}) in a sample population, and possibly additional
#'  covariates. It should be noted that these variables must be of the same length.
#'  The function then returns the variance heterogeneity \emph{p}-values using the
#'  generalized Levene's test. The residual from median is computed using the "fn"
#'  algorithm, for more details see \code{?rq}.
#'
#'
#' @param GENO the genotype of a SNP, must be a vector of 0, 1, 2's indicating the number of
#' reference alleles. The length of \code{GENO} should match that of \code{SEX}, \code{Y},
#' and any covariates \code{COV}.
#' @param SEX the genetic sex of individuals in the sample population, must be a vector of 1 and 2 or 0 and 1, depending on whether the PLINK sex code is used. Note that the default sex code is 1 for males and 2 for females in PLINK.
#' @param Y a vector of quantitative traits, such as human height.
#' @param COV a vector or matrix of covariates that are used to reduce bias due to confounding, such as age.
#' @param PLINK a logical indicating whether the SEX is coded following PLINK, i.e.
#' female coded as 2 and male colded as 1; if not, SEX is assumed to be coded 0 for females
#' and 1 for males.
#' @param test_type a character of either "ALL", printing all strategies, or "M1",
#' "M2", or "M3" for the 8 stratigies given the mean stage models, or "M3V3.2" printing
#' only the recommended test.
#'
#' @importFrom quantreg rq
#' @importFrom methods is
#' @importFrom stats resid
#' @importFrom stats lm
#' @importFrom stats complete.cases
#' @importFrom stats na.exclude
#' @importFrom stats anova
#'
#' @return a vector of Levene's test regression p-values according to the models
#' specified.
#' @export leveneReg
#' @note We recommend to quantile-normally transform \code{Y} to avoid ‘scale-effect’ where
#' the variance values tend to be proportional to mean values when stratified by \code{GENO}.
#'
#' @examples
#' N <- 5000
#' GENO <- rbinom(N, 2, 0.3)
#' sex <- rbinom(N, 1, 0.5)
#' y <- rnorm(N)
#' cov <- matrix(rnorm(N*10), ncol=10)
#' leveneReg(GENO=GENO, SEX=sex, Y=y, COV=cov, test_type="ALL")
#' leveneReg(GENO=GENO, SEX=sex, Y=y, COV=cov, test_type="M3V3.2")
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}
#'
#'


leveneReg <- function(GENO, SEX, PLINK = FALSE, Y, COV = NULL, test_type = "ALL"){

  if (missing(GENO))
    stop("The GENOtype input is missing.")

  if (missing(SEX))
    stop("The sex input is missing.")

  if (missing(Y))
    stop("The quantitative trait input is missing.")

  if (length(GENO)!=length(SEX)|length(GENO)!=length(Y)|length(Y)!=length(SEX))
    stop("Make sure the inputs have the same length.")


  if (!is.null(COV)){
    Y <- resid(lm(Y ~ COV, na.action = na.exclude))
  }

  if (PLINK){
    SEX[SEX==2] = 0
  }


  use_dat <- complete.cases(GENO, SEX, Y)
  GENO = GENO[use_dat]
  SEX = SEX[use_dat]
  Y = Y[use_dat]

  N = length(Y)

  if (sum(SEX==1) == N){
    stop("Only Males detected")

    } else if (sum(SEX==2) == N) {

    stop("Only Females detected")

  }


  if (length(table(GENO))==1) {

    stop("Monomorphic SNP detected");


  } else if (length(table(factor(interaction(GENO, SEX)))) < 4) {

    stop("Monomorphic SNP detected in either females, males");

  }


  if (test_type == "ALL"){

    res_GENO  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="fn"))));

    res_GENO_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="fn"))));
    res_GENO_sexInt <- try(as.numeric(resid(quantreg::rq(Y~GENO*SEX,na.action=na.exclude, method="fn"))));

     if (methods::is(res_GENO, "try-error")){

      Model1_G <- rep(NA, 3)
      Model1_2df <- NA
      Model1_G_Dummy <- rep(NA, 4)

    } else{

      Model1_G <- c(summary(lm(abs(res_GENO)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO)~GENO*SEX))$coef[2,4]);
      Model1_2df <- c(anova(lm(abs(res_GENO)~SEX), lm(abs(res_GENO)~GENO*SEX))$Pr[2])

      Model1_G_Dummy <- c(anova(lm(abs(res_GENO)~factor(GENO)), lm(abs(res_GENO)~1))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX+factor(GENO)), lm(abs(res_GENO)~SEX))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX*factor(GENO)), lm(abs(res_GENO)~SEX+SEX:I(GENO==1)))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX*factor(GENO)), lm(abs(res_GENO)~SEX))$Pr[2])
    }

    if (methods::is(res_GENO_Sex, "try-error")){

      Model2_G <- rep(NA, 3)
      Model2_2df <- NA
      Model2_G_Dummy <- rep(NA, 4)

    }else{


      Model2_G <- c(summary(lm(abs(res_GENO_Sex)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO_Sex)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO_Sex)~GENO*SEX))$coef[2,4]);
      Model2_2df <- c(anova(lm(abs(res_GENO_Sex)~SEX), lm(abs(res_GENO_Sex)~GENO*SEX))$Pr[2])

      Model2_G_Dummy <- c(anova(lm(abs(res_GENO_Sex)~factor(GENO)), lm(abs(res_GENO_Sex)~1))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX+factor(GENO)), lm(abs(res_GENO_Sex)~SEX))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX*factor(GENO)), lm(abs(res_GENO_Sex)~SEX+SEX:I(GENO==1)))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX*factor(GENO)), lm(abs(res_GENO_Sex)~SEX))$Pr[2])

    }

    if (methods::is(res_GENO_sexInt, "try-error")){

      Model3_G <- rep(NA, 3)
      Model3_2df <- NA
      Model3_G_Dummy <- rep(NA, 4)

    } else{


      Model3_G <- c(summary(lm(abs(res_GENO_sexInt)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO_sexInt)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO_sexInt)~GENO*SEX))$coef[2,4]);

      Model3_2df <- c(anova(lm(abs(res_GENO_sexInt)~SEX), lm(abs(res_GENO_sexInt)~GENO*SEX))$Pr[2])

      Model3_G_Dummy <-  c(anova(lm(abs(res_GENO_sexInt)~factor(GENO)), lm(abs(res_GENO_sexInt)~1))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX+factor(GENO)), lm(abs(res_GENO_sexInt)~SEX))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX*factor(GENO)), lm(abs(res_GENO_sexInt)~SEX+SEX:I(GENO==1)))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX*factor(GENO)), lm(abs(res_GENO_sexInt)~SEX))$Pr[2])

    }



    PVAL <- c(Model1_G, Model1_2df, Model1_G_Dummy, Model2_G, Model2_2df, Model2_G_Dummy, Model3_G, Model3_2df, Model3_G_Dummy)

    names(PVAL) <- c(paste("M1", c("V1", "V2", "V3"), sep=""), "M1V3.2", paste("M1", c("VNA1", "VNA2", "VNA3"), sep=""), "M1VNA3.2", paste("M2", c("V1", "V2", "V3"), sep=""), "M2V3.2", paste("M2", c("VNA1", "VNA2", "VNA3"), sep=""), "M2VNA3.2", paste("M3", c("V1", "V2", "V3"), sep=""), "M3V3.2", paste("M3", c("VNA1", "VNA2", "VNA3"), sep=""), "M3VNA3.2")


  } else if (test_type == "M1") {

    res_GENO  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="fn"))));

    if (methods::is(res_GENO, "try-error")){

      Model1_G <- rep(NA, 3)
      Model1_2df <- NA
      Model1_G_Dummy <- rep(NA, 4)

    } else{

      Model1_G <- c(summary(lm(abs(res_GENO)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO)~GENO*SEX))$coef[2,4]);
      Model1_2df <- c(anova(lm(abs(res_GENO)~SEX), lm(abs(res_GENO)~GENO*SEX))$Pr[2])

      Model1_G_Dummy <- c(anova(lm(abs(res_GENO)~factor(GENO)), lm(abs(res_GENO)~1))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX+factor(GENO)), lm(abs(res_GENO)~SEX))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX*factor(GENO)), lm(abs(res_GENO)~SEX+SEX:I(GENO==1)))$Pr[2],
                          anova(lm(abs(res_GENO)~SEX*factor(GENO)), lm(abs(res_GENO)~SEX))$Pr[2])
    }

    PVAL <- c(Model1_G, Model1_2df, Model1_G_Dummy)
    names(PVAL) <- c(paste("M1", c("V1", "V2", "V3"), sep=""), "M1V3.2", paste("M1", c("VNA1", "VNA2", "VNA3"), sep=""), "M1VNA3.2")

  } else if (test_type == "M2") {


    res_GENO_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="fn"))));

    if (methods::is(res_GENO_Sex, "try-error")){

      Model2_G <- rep(NA, 3)
      Model2_2df <- NA
      Model2_G_Dummy <- rep(NA, 4)

    }else{


      Model2_G <- c(summary(lm(abs(res_GENO_Sex)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO_Sex)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO_Sex)~GENO*SEX))$coef[2,4]);
      Model2_2df <- c(anova(lm(abs(res_GENO_Sex)~SEX), lm(abs(res_GENO_Sex)~GENO*SEX))$Pr[2])

      Model2_G_Dummy <- c(anova(lm(abs(res_GENO_Sex)~factor(GENO)), lm(abs(res_GENO_Sex)~1))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX+factor(GENO)), lm(abs(res_GENO_Sex)~SEX))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX*factor(GENO)), lm(abs(res_GENO_Sex)~SEX+SEX:I(GENO==1)))$Pr[2],
                          anova(lm(abs(res_GENO_Sex)~SEX*factor(GENO)), lm(abs(res_GENO_Sex)~SEX))$Pr[2])

    }


    PVAL <- c(Model2_G, Model2_2df, Model2_G_Dummy)
    names(PVAL) <- c(paste("M2", c("V1", "V2", "V3"), sep=""), "M2V3.2", paste("M2", c("VNA1", "VNA2", "VNA3"), sep=""), "M2VNA3.2")


  } else if (test_type == "M3") {


   res_GENO_sexInt <- try(as.numeric(resid(quantreg::rq(y~GENO*SEX,na.action=na.exclude, method="fn"))));

    if (methods::is(res_GENO_sexInt, "try-error")){

      Model3_G <- rep(NA, 3)
      Model3_2df <- NA
      Model3_G_Dummy <- rep(NA, 4)

    } else{


      Model3_G <- c(summary(lm(abs(res_GENO_sexInt)~GENO))$coef[2,4],
                    summary(lm(abs(res_GENO_sexInt)~GENO+SEX))$coef[2,4],
                    summary(lm(abs(res_GENO_sexInt)~GENO*SEX))$coef[2,4]);

      Model3_2df <- c(anova(lm(abs(res_GENO_sexInt)~SEX), lm(abs(res_GENO_sexInt)~GENO*SEX))$Pr[2])

      Model3_G_Dummy <-  c(anova(lm(abs(res_GENO_sexInt)~factor(GENO)), lm(abs(res_GENO_sexInt)~1))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX+factor(GENO)), lm(abs(res_GENO_sexInt)~SEX))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX*factor(GENO)), lm(abs(res_GENO_sexInt)~SEX+SEX:I(GENO==1)))$Pr[2],
                           anova(lm(abs(res_GENO_sexInt)~SEX*factor(GENO)), lm(abs(res_GENO_sexInt)~SEX))$Pr[2])

    }


    PVAL <- c(Model3_G, Model3_2df, Model3_G_Dummy)
    names(PVAL) <- c(paste("M3", c("V1", "V2", "V3"), sep=""), "M3V3.2", paste("M3", c("VNA1", "VNA2", "VNA3"), sep=""), "M3VNA3.2")


  } else if (test_type == "M3V3.2") {

    res_GENO_sexInt <- try(as.numeric(resid(quantreg::rq(y~GENO*SEX,na.action=na.exclude, method="fn"))));

    if (methods::is(res_GENO_sexInt, "try-error")){

      Model3_2df <- NA

    } else{

      Model3_2df <- c(anova(lm(abs(res_GENO_sexInt)~SEX), lm(abs(res_GENO_sexInt)~GENO*SEX))$Pr[2])
    }
    PVAL <- Model3_2df
    names(PVAL) <- c("M3V3.2")

  }

  return(PVAL)
}
