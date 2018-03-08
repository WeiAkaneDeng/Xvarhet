#' Generating Conceptual Sex-Stratified Distributions as Shown in Xvarhet Manuscript
#'
#' This function generates the conceptual distribution of a quantitative trait
#' stratified by sex and genotypes under the null hypothesis and alternative
#' hypothesis.
#'
#' @param RAF a numeric between 0 and 1, the reference allele frequency of the
#' SNP, assumed to be equal in females and males.
#' @param prp_m a numeric between 0 and 1, giving the porportion of males
#' in the sample.
#' @param length.OUT an integer for the number of quantitative trait values
#' at which to evaluate; the default is set to 500.
#' @param Yrange a numeric for the maximum value of the quantitative trait to
#' be evaluated assuming the distribution is symmetric about zero;
#' the default is set to 3.
#' @param y1_vect a vector for the sex-stratified mean and standard deviation
#'  in females and males under no sexual dimorphism; the default is c(0,1), with zero
#'  mean and standard deviation one in both sexes.
#' @param y2_vect a vector for the sex-stratified mean and standard deviation
#' in females and males under sexual dimorphism in mean only; the default is
#' c(-0.5,0.5,1,1), with quantitative trait in female following a normal
#' distribution \eqn{N(\mu=-0.5, \sigma=1)}.
#' @param y3_vect a vector for the sex-stratified mean and standard deviation
#' in females and males under sexual dimorphism in variance only; the default
#' is c(0,0,0.8,1.2), with quantitative trait in female following a normal
#' distribution \eqn{N(\mu=0, \sigma=0.8)}.
#' @param y4_vect a vector for the sex-stratified mean and standard deviation
#'  in females and males under sexual dimorphism in both mean and variance;
#'  the default is c(-0.5,0.5,0.8,1.2), with quantitative trait in female following
#'  a normal distribution \eqn{N(\mu=-0.5, \sigma=0.8)}.
#' @param geno_sd a numeric between 0 and the minimum of the sex-stratified
#' standard deviation, giving the contribution to standard deviation per
#' reference allele under no sexual dimorphism assuming an additive genetic
#' variance effect. The default is 0.2.
#' @param cols a logic indicating whether the outputted graphical object should be
#' in colour or gray style.
#'
#' @importFrom ggplot2 theme
#' @import ggplot2
#' @importFrom stats dnorm
#'
#' @return a ggplot2 graph object
#'
#' @export ConceptFigure
#' @examples
#' #ConceptFigure(RAF=0.2, prp_m=0.1); # not run
#' #ConceptFigure(RAF=0.2, prp_m=0.1, cols = FALSE); # not run
#'
#' @author Wei Q. Deng \email{deng@utstat.toronto.edu}
#'


ConceptFigure <- function(RAF, prp_m, length.OUT=500, Yrange = 3,
                          y1_vect = c(0, 1), y2_vect = c(-0.5, 0.5, 1, 1),
                          y3_vect = c(0,0,0.8,1.2), y4_vect = c(-0.5,0.5,0.8,1.2),
                          geno_sd = 0.2, cols = T){

	Names <- c("Female G = 0", "Female G = 1", "Female G = 2",
	"Male G = 0","Male G = 2")

	if (RAF > 1 | RAF < 0)
	  stop("Please check the reference allele frequency is between 0 and 1")


	if (prp_m > 1 | prp_m < 0)
	  stop("Please check the proportion of males is between 0 and 1")

	RAF_PROP <- c((1-prp_m)*(1-RAF)^2, (1-prp_m)*(1-RAF)* RAF*2, (1-prp_m)*(RAF)^2, (1-RAF)* prp_m, RAF* prp_m)


	xx <- seq(-Yrange, Yrange, length.out=length.OUT)

	y_1_mix <- cbind(dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]))%*% RAF_PROP

	y_1 <- c(dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]))


	y_2_mix <- cbind(dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]))%*% RAF_PROP

	y_2 <- c(dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]))


	y_3_mix <- cbind(dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]),dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]), dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]))%*% RAF_PROP

	y_3 <- c(dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]),dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]), dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]))

	y_4_mix <- cbind(dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]),dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]), dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]))%*% RAF_PROP

	y_4 <- c(dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]),dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]), dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]))



	y_5_mix <- cbind(dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd*2), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd))%*% RAF_PROP

	y_5 <- c(dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd*2), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]), dnorm(xx, mean=y1_vect[1], sd=y1_vect[2]+geno_sd))


	y_6_mix <- cbind(dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]+geno_sd), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]+2*geno_sd), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]+geno_sd))%*% RAF_PROP

	y_6 <- c(dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]+geno_sd), dnorm(xx, mean=y2_vect[1], sd=y2_vect[3]+2*geno_sd), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]), dnorm(xx, mean=y2_vect[2], sd=y2_vect[4]+geno_sd))


	y_7_mix <- cbind(dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]),dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]+geno_sd), dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]+2*geno_sd), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]+geno_sd))%*% RAF_PROP

	y_7 <- c(dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]),dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]+geno_sd), dnorm(xx, mean=y3_vect[1], sd=y3_vect[3]+2*geno_sd), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]), dnorm(xx, mean=y3_vect[2], sd=y3_vect[4]+geno_sd))


	y_8_mix <- cbind(dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]),dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]+geno_sd), dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]+2*geno_sd), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]+geno_sd))%*% RAF_PROP

	y_8 <- c(dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]),dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]+geno_sd), dnorm(xx, mean=y4_vect[1], sd=y4_vect[3]+2*geno_sd), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]), dnorm(xx, mean=y4_vect[2], sd=y4_vect[4]+geno_sd))


cases <- rep(rep(paste("case", 1:4, sep=""), each = length.OUT*5), 2);
sex <- 	rep(rep(c("F", "F", "F", "M", "M"), each = length.OUT), 8)
geno <- rep(rep(c(0,1,2,0,1), each=length.OUT), 8)
hypo <-	rep(c("Null","Alt"), each = length.OUT*5*4)


case_lab <- c(expression(paste("Sex effect: ", mu[f]==mu[m]~ "," ~sigma[f]==sigma[m], collapse="")), expression(paste("Sex effect: ", mu[f]!=mu[m]~ "," ~sigma[f]==sigma[m], collapse="")), expression(paste("Sex effect: ", mu[f]==mu[m]~ "," ~sigma[f]!=sigma[m], collapse="")), expression(paste("Sex effect: ", paste(mu[f]!=mu[m]~ "," ~sigma[f]!=sigma[m], collapse=""))))


Hypo_lab <- c(expression(atop(H[o]~": no genotype variance effect.")),
              expression(atop(H[1]~": genotype variance effect.")))


	graph_data_frame <- data.frame("X" = rep(xx, 40), "Y" = c(y_1, y_2, y_3, y_4, y_5, y_6, y_7, y_8), "Cases" = factor(cases), "Hypothese" = factor(hypo), "SEX" = factor(sex), "GENO" = factor(geno), "data" = rep("singled", 40*length.OUT))

	graph_data_frameCombined <- data.frame("X" = rep(xx, 8), "Y" = c(y_1_mix, y_2_mix, y_3_mix, y_4_mix, y_5_mix, y_6_mix, y_7_mix, y_8_mix), "Cases" = factor(rep(rep(paste("case", 1:4, sep=""), each = length.OUT), 2)), "Hypothese" = factor(rep(c("Null", "Alt"), each = length.OUT*4)), "SEX" = factor(rep("Combined", 8*length.OUT)), "GENO" = rep(NA, length.OUT*8),  "data" = rep("combined", length.OUT*8))


	graph_data_frame$hypo_label = factor(graph_data_frame$Hypothese, labels= Hypo_lab, levels=c("Null","Alt"))
	graph_data_frame$case_label = factor(graph_data_frame$Cases, labels=case_lab, levels = c("case1", "case2", "case3", "case4"))

	graph_data_frameCombined$hypo_label = factor(graph_data_frameCombined$Hypothese, labels= Hypo_lab, levels=c("Null","Alt"))
	graph_data_frameCombined$case_label = factor(graph_data_frameCombined$Cases, labels=case_lab, levels = c("case1", "case2", "case3", "case4"))

	graph_data = rbind(graph_data_frameCombined, graph_data_frame)

p_cases <- ggplot2::ggplot(graph_data,  ggplot2::aes(x = X, y= Y), group=interaction(SEX, GENO)) + ggplot2::facet_wrap(hypo_label~ case_label, scales="free_x", nrow=2, labeller = ggplot2::label_parsed) + ggplot2::geom_line(data= graph_data[graph_data $data!="combined", ], ggplot2::aes(linetype=GENO, color=SEX)) + ggplot2::geom_line(data= graph_data[graph_data$data=="combined", ], size=1) + ggplot2::geom_segment(ggplot2::aes(x = -3, y = 0, xend = 3, yend = 0))+ ggplot2::scale_linetype_manual(values=c("solid", "twodash", "dotted")) + ggplot2::scale_fill_discrete(name="SEX", breaks= c("F","M"),labels =c("Female", "Male")) + ggplot2::theme(axis.title=ggplot2::element_text(face="bold.italic", size=14), plot.title = ggplot2::element_text(face="bold", hjust=0.5), legend.text = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size=13), strip.text.x = ggplot2::element_text(size = 14), strip.text.y = ggplot2::element_text(size = 14)) + ggplot2::scale_shape_identity() + ggplot2::labs(x="Phenotype", y="Density")+ ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10)) + ggplot2::scale_colour_grey()



cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p_cases_col <- ggplot2::ggplot(graph_data, ggplot2::aes(x = X, y= Y), group=interaction(SEX, GENO)) + ggplot2::facet_wrap(hypo_label~ case_label, scales="free_x", nrow=2, labeller = ggplot2::label_parsed) + ggplot2::geom_line(data= graph_data[graph_data $data!="combined", ], ggplot2::aes(linetype=GENO, color=SEX)) + ggplot2::geom_line(data= graph_data[graph_data$data=="combined", ], size=1) + ggplot2::geom_segment(ggplot2::aes(x = -3, y = 0, xend = 3, yend = 0))+ ggplot2::scale_linetype_manual(values=c("solid", "twodash", "dotted")) + ggplot2::scale_fill_discrete(name="SEX", breaks= c("F","M"),labels =c("Female", "Male")) + ggplot2::theme(axis.title=ggplot2::element_text(face="bold.italic", size=14), plot.title = ggplot2::element_text(face="bold", hjust=0.5), legend.text = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size=13), strip.text.x = ggplot2::element_text(size = 14), strip.text.y = ggplot2::element_text(size = 14)) + ggplot2::scale_shape_identity() + ggplot2::labs(x="Phenotype", y="Density")+ ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10), panel.spacing = ggplot2::unit(2, "lines")) + ggplot2::scale_colour_manual(values= cbbPalette)


if (cols == TRUE){

  return(p_cases_col)

  } else {

  return(p_cases)
}

}




