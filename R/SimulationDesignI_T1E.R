#' Type I Error (T1E) Rates under Simulation Design I
#'
#' This function returns the type I error rates of testing strategies
#' under the simulation conditions specified by a linear regression model.
#' The quantitative trait $Y$ is simulated from the following model without any genetic effect:
#' \eqn{Y =\beta_{E} E + \beta_S S + \beta_{SE} SE + \epsilon}, where \eqn{E} and \eqn{\epsilon} are assumed to follow a standard normal distribution.
#'
#'
#' @param Nf an integer for the number of females in the sample population.
#' @param Nm an integer for the number of males in the sample population.
#' @param MAF_F a numerical between 0 and 0.5 for minor allele frequency of the simulated SNP in females.
#' @param MAF_M a numerical between 0 and 0.5 for minor allele frequency of the simulated SNP in males.
#' @param betaS a numerical value for the marginal sex effect on the simulated trait.
#' @param betaE a numerical value for the marginal environmental covariate effect on the simulated trait.
#' @param betaSE a numerical value for the interaction effect between sex and an environmental covariate.
#' @param nbSim an integer for the number of simulations.
#' @param alpha a numerical between 0 and 1 for the significance threshold, the default value is 0.05.
#' @param nbDigits an integer indicating the number of printed digits for empirical T1E rates.
#' @param test_type a character of either "ALL", printing all strategies, or "M1",
#' "M2", or "M3" for the 8 stratigies given the mean stage models, or "M3V3.2" printing
#' only the recommended test.
#' @param MisspecX a logic indicating whether the T1E should be calculated under a
#' mis-specified X-inactivation status, i.e. model generated under X-inactivation
#' but tested under absence of X-inactivation.
#'
#' @return a vector of type I error rates of the specified testing strategies.
#'
#' @importFrom stats rnorm
#' @importFrom stats resid
#' @importFrom stats lm
#' @importFrom stats complete.cases
#' @importFrom stats na.exclude
#' @importFrom stats anova
#'
#' @note We recommend to run the simulation in the background as it might take
#' long for more than 500 simulations on a 1.6 GHz Core.
#' @export T1EsimDesign1
#'
#' @examples
#' T1EsimDesign1(Nf=5000, Nm=5000, MAF_F=0.2, MAF_M=0.2, betaE=0.1,
#' betaS=0.2, betaSE=0, nbSim=5, alpha=0.05, nbDigits = 3, MisspecX =FALSE,
#' test_type = "ALL")
#'
#' T1EsimDesign1(Nf=5000, Nm=5000, MAF_F=0.2, MAF_M=0.2, betaE=0.1,
#' betaS=0.2, betaSE=0, nbSim=50, alpha=0.05, nbDigits = 3, MisspecX =FALSE,
#' test_type = "M3V3.2")
#'
#' T1EsimDesign1(Nf=5000, Nm=5000, MAF_F=0.2, MAF_M=0.2, betaE=0.1,
#' betaS=0.2, betaSE=0, nbSim=100, alpha=0.05, nbDigits = 3, Misspe
#' =TRUE,test_type = "ALL") # longer run time with 100 simulations
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}
#'
#'

T1EsimDesign1 <- function(Nf, Nm, MAF_F, MAF_M, betaE, betaS, betaSE,
                          nbSim=100, alpha=0.05, nbDigits = 3,
                          MisspecX = FALSE, test_type="ALL"){

  if (missing(Nf)|missing(Nm))
    stop("The sample size input is missing.")

  if (missing(MAF_F)| missing(MAF_M))
    stop("The minor allele frequency input is missing.")

  if (missing(betaE)| missing(betaSE)|missing(betaS))
    stop("The quantitative trait input is missing.")

  if (is.integer(Nf) | is.integer(Nm))
    stop("Sample size must be an integer.")

  if ((MAF_F>0.5 | MAF_F < 0.05) | (MAF_M>0.5 | MAF_M < 0.05))
    stop("Minor allele frequency must be a numerical value betweem 0.05 and 0.5.")

  generate_pval <- function(y, geno, sex , cov=NULL, PLINK=FALSE, tests = "ALL"){
    g11 <- leveneTests(GENO=geno, SEX=sex, Y=y, COV=cov, PLINK=PLINK)
    g16 <- leveneReg(GENO=geno, SEX=sex, PLINK = PLINK, Y=y, COV = cov, test_type = tests)
    return(c(g11, g16))
  }


  SEX =c(rep(0, Nf), rep(1, Nm))

  female_geno <- c(rep(0, round(Nf*(1-MAF_F)^2)), rep(1, round(Nf*(1-MAF_F)*2* MAF_F)), rep(2, Nf-round(Nf*(1-MAF_F)*2* MAF_F)-round(Nf*(1-MAF_F)^2)))

  male_geno <- c(rep(0, round(Nm*(1-MAF_M))), rep(1, Nm-round(Nm*(1-MAF_M))))

  GENOs <- c(female_geno, male_geno);

  run_sim <- function(test_type){
    e2 <- rnorm(Nf + Nm);
    yy <- rnorm(Nf + Nm) + betaS*SEX + betaE*e2 + betaSE*e2*SEX
    pval1 <- generate_pval(y=yy, geno=GENOs, sex=SEX , cov=NULL, PLINK=FALSE, tests = test_type)
    return(pval1)
    }

    All_Pval_Sets <- replicate(nbSim, run_sim(test_type=test_type))
    typeIerror <- apply(All_Pval_Sets, 1, function(dd) sum(dd < alpha)/length(dd))

    out_table <- data.frame("FemaleN" = Nf, "MaleN" = Nm,
                                       "MAF_Female" = MAF_F, "MAF_Male" = MAF_M,
                                       "betaEnvironment" = betaE, "betaSEX" = betaS,
                                       "beta_SE" = betaSE, as.data.frame(t(typeIerror)))

 return(out_table)
}
