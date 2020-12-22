#' Levene's regression tests for variance homogeneity by SNP genotype
#'
#'  The function takes as input the genotype of a SNP (\code{GENO}), the SEX (\code{SEX}), and
#'  a quantitative trait (\code{Y}) in a sample population, and possibly additional
#'  covariates. It should be noted that these variables must be of the same length.
#'  The function then returns the variance heterogeneity \emph{p}-values using the
#'  generalized Levene's test. The residual function could alternatively be replaced with the quantile regression quantreg::rq following the "fn" algorithm, for more details see \code{?quantreg::rq}.
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
#' "M2", or "M3" for the 8 stratigies given the mean stage models, or "REC" printing
#' only the recommended tests.
#'
#' @import stats
#' @import quantreg
#' @import methods
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
#' Y <- rnorm(N)
#' cov <- matrix(rnorm(N*10), ncol=10)
#' leveneReg(GENO=GENO, SEX=sex, Y=Y, COV=cov, test_type="ALL")
#' leveneReg(GENO=GENO, SEX=sex, Y=Y, COV=cov, test_type="REC")
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}
#'
#'


leveneReg <- function(GENO, SEX, PLINK = FALSE, Y, COV = NULL, test_type = "REC"){

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


  if (class(PLINK)!="logical"){
  warning("PLINK only takes arguments TRUE or FALSE, analysis will be completed ignoring this input value (assuming SEX=0 for females and SEX=1 for males)")
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

    #res_GENO  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="fn"))));
	 res_geno  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="pfn"))));

    #res_GENO_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="fn"))));
    res_geno_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="pfn"))));

    #res_GENO_sexInt <- try(as.numeric(resid(quantreg::rq(Y~GENO*SEX,na.action=na.exclude, method="fn"))));
	res_geno_sexInt <- try(as.numeric(resid(quantreg::rq(Y~GENO*SEX,na.action=na.exclude, method="pfn"))));


	VAR <- tapply(Y, SEX, sd)
	w <- ifelse(SEX == as.numeric(names(VAR)[1]), VAR[1], VAR[2])

	r1_geno <- res_geno/w
	r1_geno_Sex <- res_geno_Sex/w
	r1_geno_SexInt <- res_geno_sexInt/w


     if (methods::is(res_geno, "try-error")){

      wModel1_G <- rep(NA, 3)
      wModel1_2df <- NA
      wModel1_G_Dummy <- rep(NA, 4)

    } else{

     wModel1_G <- c(summary(lm(abs(r1_geno)~GENO))$coef[2,4],
				summary(lm(abs(r1_geno)~GENO+SEX))$coef[2,4],
				summary(lm(abs(r1_geno)~GENO*SEX))$coef[2,4]);

	wModel1_2df <- c(anova(lm(abs(r1_geno)~SEX), lm(abs(r1_geno)~GENO*SEX))$Pr[2])

	wModel1_G_Dummy <- c(anova(lm(abs(r1_geno)~factor(GENO)), lm(abs(r1_geno)~1))$Pr[2],
	anova(lm(abs(r1_geno)~SEX+factor(GENO)), lm(abs(r1_geno)~SEX))$Pr[2],
	anova(lm(abs(r1_geno)~SEX*factor(GENO)), lm(abs(r1_geno)~SEX+SEX:I(GENO==1)))$Pr[2],
	anova(lm(abs(r1_geno)~SEX*factor(GENO)), lm(abs(r1_geno)~SEX))$Pr[2])


    }

     if (methods::is(res_geno_Sex, "try-error")){

      wModel2_G <- rep(NA, 3)
      wModel2_2df <- NA
      wModel2_G_Dummy <- rep(NA, 4)

    }else{


	wModel2_G <- c(summary(lm(abs(r1_geno_Sex)~GENO))$coef[2,4],
				summary(lm(abs(r1_geno_Sex)~GENO+SEX))$coef[2,4],
				summary(lm(abs(r1_geno_Sex)~GENO*SEX))$coef[2,4]);

	wModel2_2df <- c(anova(lm(abs(r1_geno_Sex)~SEX), lm(abs(r1_geno_Sex)~GENO*SEX))$Pr[2])

	wModel2_G_Dummy <- c(anova(lm(abs(r1_geno_Sex)~factor(GENO)), lm(abs(r1_geno_Sex)~1))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX+factor(GENO)), lm(abs(r1_geno_Sex)~SEX))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX*factor(GENO)), lm(abs(r1_geno_Sex)~SEX+SEX:I(GENO==1)))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX*factor(GENO)), lm(abs(r1_geno_Sex)~SEX))$Pr[2])


    }

     if (methods::is(res_geno_sexInt, "try-error")){

      wModel3_G <- rep(NA, 3)
      wModel3_2df <- NA
      wModel3_G_Dummy <- rep(NA, 4)

    } else{

	wModel3_G <- c(summary(lm(abs(r1_geno_SexInt)~GENO))$coef[2,4],
			summary(lm(abs(r1_geno_SexInt)~GENO+SEX))$coef[2,4],
			summary(lm(abs(r1_geno_SexInt)~GENO*SEX))$coef[2,4]);

	wModel3_2df <- c(anova(lm(abs(r1_geno_SexInt)~SEX), lm(abs(r1_geno_SexInt)~GENO*SEX))$Pr[2])

	wModel3_G_Dummy <-  c(anova(lm(abs(r1_geno_SexInt)~factor(GENO)), lm(abs(r1_geno_SexInt)~1))$Pr[2],
anova(lm(abs(r1_geno_SexInt)~SEX+factor(GENO)), lm(abs(r1_geno_SexInt)~SEX))$Pr[2],
anova(lm(abs(r1_geno_SexInt)~SEX*factor(GENO)), lm(abs(r1_geno_SexInt)~SEX+SEX:I(GENO==1)))$Pr[2],
anova(lm(abs(r1_geno_SexInt)~SEX+factor(GENO)+SEX:I(GENO==1)), lm(abs(r1_geno_SexInt)~SEX))$Pr[2])

    }

    PVAL <- c(wModel1_G, wModel1_2df, wModel1_G_Dummy, wModel2_G, wModel2_2df, wModel2_G_Dummy, wModel3_G, wModel3_2df, wModel3_G_Dummy)

    names(PVAL) <- c(paste("wM1", c("V1", "V2", "V3"), sep=""), "wM1V3.2", paste("wM1", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM1VNA3.3", paste("wM2", c("V1", "V2", "V3"), sep=""), "wM2V3.2", paste("W_M2", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM2VNA3.3", paste("wM3", c("V1", "V2", "V3"), sep=""), "wM3V3.2", paste("wM3", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM3VNA3.3")


  } else if (test_type == "M1") {

    #res_GENO  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="fn"))));
	res_geno  <- try(as.numeric(resid(quantreg::rq(Y~GENO,na.action=na.exclude, method="fn"))));

	VAR <- tapply(Y, SEX, sd)
	w <- ifelse(SEX == as.numeric(names(VAR)[1]), VAR[1], VAR[2])

	r1_geno <- res_geno/w

    if (methods::is(res_geno, "try-error")){

      wModel1_G <- rep(NA, 3)
      wModel1_2df <- NA
      wModel1_G_Dummy <- rep(NA, 4)

    } else{

     wModel1_G <- c(summary(lm(abs(r1_geno)~GENO))$coef[2,4],
				summary(lm(abs(r1_geno)~GENO+SEX))$coef[2,4],
				summary(lm(abs(r1_geno)~GENO*SEX))$coef[2,4]);

	wModel1_2df <- c(anova(lm(abs(r1_geno)~SEX), lm(abs(r1_geno)~GENO*SEX))$Pr[2])

	wModel1_G_Dummy <- c(anova(lm(abs(r1_geno)~factor(GENO)), lm(abs(r1_geno)~1))$Pr[2],
	anova(lm(abs(r1_geno)~SEX+factor(GENO)), lm(abs(r1_geno)~SEX))$Pr[2],
	anova(lm(abs(r1_geno)~SEX*factor(GENO)), lm(abs(r1_geno)~SEX+SEX:I(GENO==1)))$Pr[2],
	anova(lm(abs(r1_geno)~SEX*factor(GENO)), lm(abs(r1_geno)~SEX))$Pr[2])

    PVAL <- c(wModel1_G, wModel1_2df, wModel1_G_Dummy)
    names(PVAL) <- c(paste("wM1", c("V1", "V2", "V3"), sep=""), "wM1V3.2", paste("wM1", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM1VNA3.3")

    }

   } else if (test_type == "M2") {

	#res_GENO_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="fn"))));
    res_geno_Sex <- try(as.numeric(resid(quantreg::rq(Y~GENO+SEX,na.action=na.exclude, method="fn"))));

	VAR <- tapply(Y, SEX, sd)
	w <- ifelse(SEX == as.numeric(names(VAR)[1]), VAR[1], VAR[2])

	r1_geno_Sex <- res_geno_Sex/w

    if (methods::is(res_geno_Sex, "try-error")){

      wModel2_G <- rep(NA, 3)
      wModel2_2df <- NA
      wModel2_G_Dummy <- rep(NA, 4)

    } else {


	wModel2_G <- c(summary(lm(abs(r1_geno_Sex)~GENO))$coef[2,4],
				summary(lm(abs(r1_geno_Sex)~GENO+SEX))$coef[2,4],
				summary(lm(abs(r1_geno_Sex)~GENO*SEX))$coef[2,4]);

	wModel2_2df <- c(anova(lm(abs(r1_geno_Sex)~SEX), lm(abs(r1_geno_Sex)~GENO*SEX))$Pr[2])

	wModel2_G_Dummy <- c(anova(lm(abs(r1_geno_Sex)~factor(GENO)), lm(abs(r1_geno_Sex)~1))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX+factor(GENO)), lm(abs(r1_geno_Sex)~SEX))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX*factor(GENO)), lm(abs(r1_geno_Sex)~SEX+SEX:I(GENO==1)))$Pr[2],
	anova(lm(abs(r1_geno_Sex)~SEX*factor(GENO)), lm(abs(r1_geno_Sex)~SEX))$Pr[2])


    PVAL <- c(wModel2_G, wModel2_2df, wModel2_G_Dummy)
    names(PVAL) <- c(paste("wM2", c("V1", "V2", "V3"), sep=""), "wM2V3.2", paste("W_M2", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM2VNA3.3")

	}

  } else if (test_type == "M3") {


	res_geno_sexInt <- try(as.numeric(resid(quantreg::rq(Y~GENO*SEX,na.action=na.exclude, method="fn"))));


    if (methods::is(res_geno_sexInt, "try-error")){

      wModel3_G <- rep(NA, 3)
      wModel3_2df <- NA
      wModel3_G_Dummy <- rep(NA, 4)

    } else {

	wModel3_G <- c(summary(lm(abs(res_geno_sexInt)~GENO))$coef[2,4],
			summary(lm(abs(res_geno_sexInt)~GENO+SEX))$coef[2,4],
			summary(lm(abs(res_geno_sexInt)~GENO*SEX))$coef[2,4]);

	wModel3_2df <- c(anova(lm(abs(res_geno_sexInt)~SEX), lm(abs(res_geno_sexInt)~GENO*SEX))$Pr[2])

	wModel3_G_Dummy <-  c(anova(lm(abs(res_geno_sexInt)~factor(GENO)), lm(abs(r1_geno_SexInt)~1))$Pr[2],
anova(lm(abs(res_geno_sexInt)~SEX+factor(GENO)), lm(abs(res_geno_sexInt)~SEX))$Pr[2],
anova(lm(abs(res_geno_sexInt)~SEX*factor(GENO)), lm(abs(res_geno_sexInt)~SEX+SEX:I(GENO==1)))$Pr[2],
anova(lm(abs(res_geno_sexInt)~SEX+factor(GENO)+SEX:I(GENO==1)), lm(abs(res_geno_sexInt)~SEX))$Pr[2])

    }

    PVAL <- c(wModel3_G, wModel3_2df, wModel3_G_Dummy)
    names(PVAL) <- c(paste("wM3", c("V1", "V2", "V3"), sep=""), "wM3V3.2", paste("wM3", c("VNA1", "VNA2", "VNA3.2"), sep=""), "wM3VNA3.3")


  } else if (test_type == "REC") {

	res_geno_sexInt <- try(as.numeric(resid(quantreg::rq(Y~GENO*SEX,na.action=na.exclude, method="fn"))));

    if (methods::is(res_geno_sexInt, "try-error")){

      wModel3_2df <- NA
      wModel3_G_3df <- NA

      } else{

	wModel3_2df <- c(anova(lm(abs(res_geno_sexInt)~SEX), lm(abs(res_geno_sexInt)~GENO*SEX))$Pr[2])
    wModel3_G_3df <- anova(lm(abs(res_geno_sexInt)~SEX+factor(GENO)+SEX:I(GENO==1)), lm(abs(res_geno_sexInt)~SEX))$Pr[2]
    }
    PVAL <- c(wModel3_2df, wModel3_G_3df)
    names(PVAL) <- c("wM3V3.2", "wM3VNA3.3")

  }



  return(PVAL)
}
