# Data and code for **SEX-SPECIFIC ADDITIVE GENETIC VARIANCES AND CORRELATIONS FOR FITNESS IN A SONG SPARROW (MELOSPIZA MELODIA) POPULATION SUBJECT TO NATURAL IMMIGRATION AND INBREEDING**
# [biorXiv](https://www.biorxiv.org/content/early/2018/02/26/272138)
# Wolak, Arcese, Keller, Nietlisbach, & Reid 

# These data come from the long-term song sparrow field study on Mandarte Island, BC, Canada.

# The data provided here are sufficient to replicate the analyses presented in the above paper, and are therefore a restricted subset of the full Mandarte dataset.

# If you are interested in running additional analyses that require further data then please get in touch with at least one (preferably all) of the following project leaders:
##  - Prof Peter Arcese (University of British Columbia): peter.arcese<at>ubc.ca
##  - Prof Lukas Keller (University of Zurich): lukas.keller<at>ieu.uzh.ch
##  - Prof Jane Reid (University of Aberdeen): jane.reid<at>abdn.ac.uk

# We are always happy to develop collaborations with researchers who have good ideas for new analyses.

# We would also appreciate it if you could let us know if you are intending to make use of the dataset below in order to facilitate coordination of different ongoing research efforts and allow us to keep track of all outputs from the long-term field study.








################################################################################


#		BIVARIATE FITNESS


################################################################################
rm(list = ls())
setwd("~/Dropbox/AberdeenPostdoc/Fitness/fitnessMS/dryadGithub")#setwd("<<Insert correct path here>>")

library(MCMCglmm)
library(wolakR) #<-- install from GitHub with devtools: `devtools::install_github("matthewwolak/wolakR")`
library(QGglmm)
########################################
#####################
# Read in data
#####################
fitListAinv <- read.table("./Wolak_et_al_fitness_listAinv.txt", header = TRUE)
  fitAinv <- as(sparseMatrix(i = fitListAinv[, "row"], j = fitListAinv[, "column"],
	x = fitListAinv[, "Ainv"],
	dims = rep(max(fitListAinv[, 1]), 2),
	dimnames = list(as.character(seq(1, max(fitListAinv[, 1]), 1)), NULL),
	symmetric = TRUE, index1 = TRUE), "dgCMatrix")

fitDat <- read.table("./Wolak_et_al_fitness_Data.txt", header = TRUE)
  fitDat$id <- as.factor(fitDat$id)
  fitDat$cohort <- as.factor(fitDat$cohort)

#####################
## Prior
#####################
bivPEpriorNukp1 <- list(R = list(V = diag(2), nu = 2),
	G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
		G2 = list(V = diag(2)*0.02, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

#####################
## Alternative Priors
## "Prior 2" in Supporting Information
bivPEpriorNuk <- list(R = list(V = diag(2), nu = 2),
	G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
		G2 = list(V = diag(2), nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

## "Prior 3" in Supporting Information Appendix S6
bivPEpriorIW <- list(R = list(V = diag(2), nu = 1.002),
	G = list(G1 = list(V = diag(2), nu = 1.002, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
		G2 = list(V = diag(2), nu = 1.002, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))
# End Alternative Priors
########################

nsamp <- 5000
BURN <- 5000; THIN <- 2500
(NITT <- BURN + THIN*nsamp)

fitMod <- MCMCglmm(fitness ~ sex-1 + sex:(f_coeff + immGG),
	random = ~ us(sex):id + us(sex):cohort,
	rcov = ~ idh(sex):units,
	ginverse = list(id = fitAinv),
	data = fitDat,
	prior = bivPEpriorNukp1,
	family = "poisson",
	nitt = NITT, thin = THIN, burnin = BURN)

# To run the above model in parallel MCMC chains, see code in the `hpc` repository
## [hpc - TorqueMoab](https://github.com/matthewwolak/hpc/tree/master/TorqueMoab/R-MCMCglmm)
## Particularly the files: `parMCMC-R.sh` and `parMCMC_postProcess.R`


###########################
# replace lower triangle with correlations
## First, save covariance version for use below
fitModCov <- fitMod
fitMod$VCV[, 2] <- posterior.cor(fitMod$VCV[, 1:4])[, 2]
fitMod$VCV[, 6] <- posterior.cor(fitMod$VCV[, 5:8])[, 2]



###############################

# LATENT SCALE HERITABILITIES

###############################
# de Villemereuil et al. 2016
fitH2latent <- with(fitMod, as.mcmc(matrix(c(VCV[, 1] / rowSums(VCV[, c(1,5,9)]),
	VCV[, 4] / rowSums(VCV[, c(4,8,10)])),
		ncol = 2, dimnames = list(NULL, c("femaleH2latent", "maleH2latent")))))


###############################

# OBSERVED SCALE COMPUTATIONS

###############################
# Assume f_coeff=0 and immigrant genetic group=0 (i.e., the 'founder' population)
fitMVparams_wrap <- function(i, itwrite = FALSE){
 if(itwrite) cat(i, " of ", nrow(fitModCov$VCV), "\n")
  Gi <- matrix(fitModCov$VCV[i, 1:4], 2, 2)
  Pi <- Gi + matrix(fitModCov$VCV[i, 5:8], 2, 2) + diag(fitModCov$VCV[i, 9:10])
  ovcvi <- QGmvparams(mu = fitMod$Sol[i, 1:2], vcv.G = Gi, vcv.P = Pi,
	models = rep("Poisson.log", 2),
	rel.acc = 0.01, width = 10, verbose = FALSE, mask = NULL)
 c(unlist(ovcvi), CORAobs = cov2cor(ovcvi[[3]])[1, 2], H = diag(ovcvi[[3]]) / diag(ovcvi[[2]]), CVA = 100*sqrt(diag(ovcvi[[3]])) / ovcvi[[1]], IA = diag(ovcvi[[3]]) / (ovcvi[[1]]^2))
}



#XXX XXX takes ~2.3 HOURS for 5,000 MCMC samples
system.time(fitObs <- as.mcmc(t(sapply(seq(nrow(fitMod$VCV)),
	FUN = fitMVparams_wrap, itwrite = TRUE))))



# Full posterior distribution of the relevant model parameters
fitPost <- as.mcmc(cbind(fitMod$Sol[, 1:6], fitMod$VCV,
	fitH2latent, fitObs))
attr(fitPost, "mcpar") <- attr(fitMod$VCV, "mcpar")
  dimnames(fitPost)[[2L]][19:35] <- c("mean.obsFemale", "mean.obsMale",
    "vcv.P.obsFemale.Female", "vcv.P.obsMale.Female", "vcv.P.obsFemale.Male", "vcv.P.obsMale.Male",
    "vcv.G.obsFemale.Female", "vcv.G.obsMale.Female", "vcv.G.obsFemale.Male", "vcv.G.obsMale.Male",
    "CORA.obsFemale.male", "h2.obsFemale", "h2.obsMale",
    "CVA.obsFemale", "CVA.obsMale", "IA.obsFemale", "IA.obsMale")




##########################

# FIGURE 2

##########################

fitnessPostHisBreaks <- 50
fitnessdf1 <- 1
par(mfcol = c(1, 3), mar = c(8, 8, 2, 1), mgp = c(4,1.5,0), cex.lab = 2, cex.axis = 1.9)
  fitppfva <- postPlot(fitPost[, "sexFemale:sexFemale.id"],
	xlab = expression("V"[" A"]~" female fitness"),
	histbreaks = fitnessPostHisBreaks,
	xlim = c(0, 7), ylim = c(0, 0.53),
	ylab = "\nDensity",
	at2 = seq(0, 0.5, by = 0.1))
     # Plot prior densities below y-axis maximum
     poh <- fitppfva$postDensity$histogram
     # evaluate densities at some fill-in of histogram
     prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10), 2*poh$mids[1]/10)
     prd <- df(prdx / 1000, df1 = fitnessdf1, df2 = 3, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(2*poh$mids[1]/10 * prd)# prd area
     sprd <- prd / prda
   lines(sprd[which(sprd <= 0.5)] ~ prdx[which(sprd <= 0.5)], lwd = 3, col = "blue")
   text(x = 6.9, y = 0.52, "A", cex = 3)

  fitppmva <- postPlot(fitPost[, "sexMale:sexMale.id"],
	xlab = expression("V"[" A"]~" male fitness"),
	histbreaks = fitnessPostHisBreaks,
	xlim = c(0, 7), ylim = c(0, 0.53),
	ylab = "",
	at2 = seq(0, 0.5, by = 0.1))
     # Plot prior densities below y-axis maximum
     poh <- fitppmva$postDensity$histogram
     # evaluate densities at some fill-in of histogram
     prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10), 2*poh$mids[1]/10)
     prd <- df(prdx / 1000, df1 = fitnessdf1, df2 = 3, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(2*poh$mids[1]/10 * prd)# prd area
     sprd <- prd / prda
   lines(sprd[which(sprd <= 0.5)] ~ prdx[which(sprd <= 0.5)], lwd = 3, col = "blue")
    text(x = 6.9, y = 0.52, "B", cex = 3)

  postPlot(fitPost[, "sexMale:sexFemale.id"], bw = 0.15,
	xlab = expression("Cross-sex r"[" A"]~" fitness"),
	histbreaks = fitnessPostHisBreaks,
	xlim = c(-1, 1.11), ylim = c(0, 1.58),
	ylab = "")
   # Prior is effectively flat across entire range, so to have area=1, y=0.5
   lines(rep(0.5, 100) ~ seq(-1, 1, length.out = 100), lwd = 3, col = "blue")
    text(x = 1.11, y = 1.55, "C", cex = 3)




#	XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX
################################################################################















################################################################################


#		TRIVARIATE FITNESS COMPONENTS


################################################################################
rm(list = ls())
setwd("~/Dropbox/AberdeenPostdoc/Fitness/fitnessMS/dryadGithub")#setwd("<<Insert correct path here>>")

library(MCMCglmm)
library(wolakR) #<-- install from GitHub with devtools: `devtools::install_github("matthewwolak/wolakR")`
library(QGglmm)
########################################
#####################
# Read in data
#####################
fitcompListAinv <- read.table("./Wolak_et_al_fitnessComp_listAinv.txt", header = TRUE)
  fitcompAinv <- as(sparseMatrix(i = fitcompListAinv[, "row"],
	j = fitcompListAinv[, "column"],
	x = fitcompListAinv[, "Ainv"],
	dims = rep(max(fitcompListAinv[, 1]), 2),
	dimnames = list(as.character(seq(1, max(fitcompListAinv[, 1]), 1)), NULL),
	symmetric = TRUE, index1 = TRUE), "dgCMatrix")

fitcompDat <- read.table("./Wolak_et_al_fitnessComp_Data.txt", header = TRUE)
  fitcompDat$id <- as.factor(fitcompDat$id)
  fitcompDat$cohort <- as.factor(fitcompDat$cohort)
  fitcompDat$socdam <- as.factor(fitcompDat$socdam)
  fitcompDat$socsire <- as.factor(fitcompDat$socsire)
  fitcompDat$brood <- as.factor(fitcompDat$brood)
  fitcompDat$pair <- as.factor(fitcompDat$pair)

#####################
## Prior
#####################
PEpriorNukp1 <- list(R = list(V = diag(3), nu = 2, fix = 3),
  G = list(G1 = list(V = diag(3)*0.02, nu = 4, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G2 = list(V = diag(3)*0.02, nu = 4, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G3 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G4 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G5 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G6 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10))))
########################
# Alternative Priors:
triPEpriorNuk <- list(R = list(V = diag(3), nu = 2, fix = 3),
  G = list(G1 = list(V = diag(3), nu = 3, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G2 = list(V = diag(3), nu = 3, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G3 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G4 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G5 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G6 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10))))

triPEpriorIW <- list(R = list(V = diag(3), nu = 2, fix = 3),
  G = list(G1 = list(V = diag(3), nu = 2.002, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G2 = list(V = diag(3), nu =  2.002, alpha.mu = rep(0, 3), alpha.V = diag(3)*c(1000, 1000, 10)),
	G3 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G4 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G5 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G6 = list(V = diag(1), nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10))))
# End Alternative Priors
########################



nsamp <- 5000
BURN <- 25000; THIN <- 8000
(NITT <- BURN + THIN*nsamp)

fitcompMod <- MCMCglmm(cbind(femLRS, maleLRS, jSurv) ~ at.level(trait, 1:2)-1 + at.level(trait, 3):sex-1 + at.level(trait, 1:2):(f_coeff + immGG) + at.level(trait, 3):(f_coeff + DFEs + immGG),
	random = ~ us(trait):id + us(trait):cohort + idh(at.level(trait, 3)):socdam + idh(at.level(trait, 3)):socsire + idh(at.level(trait, 3)):brood + idh(at.level(trait, 3)):pair,
	rcov = ~ idh(trait):units,
	ginverse = list(id = fitcompAinv),
	prior = PEpriorNukp1,
	family = c("poisson", "poisson", "categorical"),
	data = fitcompDat,
	slice = TRUE,
	nitt = NITT, thin = THIN, burnin = BURN)


# To run the above model in parallel MCMC chains, see code in the `hpc` repository
## [hpc - TorqueMoab](https://github.com/matthewwolak/hpc/tree/master/TorqueMoab/R-MCMCglmm)
## Particularly the files: `parMCMC-R.sh` and `parMCMC_postProcess.R`




###########################
# replace lower triangle with correlations
## First, save covariance version for use below
fitcompModCov <- fitcompMod
fitcompMod$VCV[, c(2,3,6)] <- posterior.cor(fitcompMod$VCV[, 1:9])[, c(2,3,6)]
fitcompMod$VCV[, c(11,12,15)] <- posterior.cor(fitcompMod$VCV[, 10:18])[, c(2,3,6)]




###########################
# Difference between female and male VA for LRS
## Male - female
VaLRSDiff <- with(triMod, VCV[,5]-VCV[,1])
postTable(VaLRSDiff)

# Difference between female and male inbreeding depression
## Male - female
fLRSDiff <- with(triMod, Sol[,6]-Sol[,5])
postTable(fLRSDiff)




###############################

# LATENT SCALE HERITABILITIES

###############################
# de Villemereuil et al. 2016
fitcompH2latent <- with(fitcompModCov, as.mcmc(matrix(c(VCV[, 1] / rowSums(VCV[, c(1,10,23)]),
	VCV[, 5] / rowSums(VCV[, c(5,14,24)]),
	VCV[, 9] / rowSums(VCV[, c(9,18:22,25)])), ncol = 3,
		dimnames = list(NULL, c("femaleLRSh2latent", "maleLRSh2latent", "jSurvH2latent")))))



###############################

# OBSERVED SCALE COMPUTATIONS

###############################
# Assume f_coeff=0 and immigrant genetic group=0 (i.e., the 'founder' population)
# juvenile survival: DFEs (date first egg laid in a nest) is centered. To evaluate this at the mean=0 allows us to drop it.
# juvenile survival: female vs. male intercept - predict for each individual then take mean of population
## Create the `X` design matrix just for sex and juvenile survival
Xsex <- model.matrix(~ sex - 1, data = subset(fitcompDat, subset = !is.na(jSurv)))
## Calculate proportion female and male
cmXsex <- matrix(colMeans(Xsex), ncol = 2)
## Calculate intercepts by weighting posterior samples of intercept by proportion female and male
jSurvMeanInt <- sapply(seq(nrow(fitcompMod$VCV)), FUN = function(i){c(cmXsex %*% matrix(fitcompMod$Sol[i, c("at.level(trait, 3):sexFemale", "at.level(trait, 3):sexMale")], ncol = 1))})



fitcompMVparams_wrap <- function(i, itwrite = FALSE){
 if(itwrite && i == 1) st <<- Sys.time()
 if(itwrite) cat(i, " of ", nrow(fitcompModCov$VCV), "\t", format(Sys.time()-st, format = "%H:%M:%S"), "\n")
  Gi <- matrix(fitcompModCov$VCV[i, 1:9], 3, 3)
  Pi <- Gi + matrix(fitcompModCov$VCV[i, 10:18], 3, 3) + diag(c(0,0,sum(fitcompModCov$VCV[i, 19:22]))) + diag(fitcompModCov$VCV[i, 23:25])
  ovcvi <- QGmvparams(mu = c(fitcompMod$Sol[i, 1:2], jSurvMeanInt[i]),
	vcv.G = Gi, vcv.P = Pi,
	models = c(rep("Poisson.log", 2), "binom1.logit"),
	rel.acc = 0.01, width = 10, verbose = FALSE, mask = NULL)

 c(unlist(ovcvi), CORA = cov2cor(ovcvi[[3]])[cbind(c(2, 3, 3), c(1, 1, 2))], H = diag(ovcvi[[3]]) / diag(ovcvi[[2]]), CVA = 100*sqrt(diag(ovcvi[[3]])[1:2]) / ovcvi[[1]][1:2], IA = diag(ovcvi[[3]])[1:2] / (ovcvi[[1]][1:2]^2))
}



#XXX XXX takes ~11.2 hours for 5,000 iterations
system.time(fitcompObs <- as.mcmc(t(sapply(seq(nrow(fitcompMod$VCV)),
	FUN = fitcompMVparams_wrap, itwrite = TRUE))))



# Full posterior distribution of the relevant model parameters
fitcompPost <- as.mcmc(cbind(fitcompMod$Sol[, 1:11],
		fitcompMod$VCV, fitcompH2latent, fitcompObs))
attr(fitcompPost, "mcpar") <- attr(fitcompMod$VCV, "mcpar")
  dimnames(fitcompPost)[[2L]][40:70] <- c("mean.obsFemLRS", "mean.obsMaleLRS", "mean.obsJSurv",
    "vcv.P.obsFemLRS.FemLRS", "vcv.P.obsMaleLRS.FemLRS", "vcv.P.obsJSurv.FemLRS",
    "vcv.P.obsFemLRS.MaleLRS", "vcv.P.obsMaleLRS.MaleLRS", "vcv.P.obsJSurv.MaleLRS",
    "vcv.P.obsFemLRS.JSurv", "vcv.P.obsMaleLRS.JSurv", "vcv.P.obsJSurv.JSurv",
    "vcv.G.obsFemLRS.FemLRS", "vcv.G.obsMaleLRS.FemLRS", "vcv.G.obsJSurv.FemLRS",
    "vcv.G.obsFemLRS.MaleLRS", "vcv.G.obsMaleLRS.MaleLRS", "vcv.G.obsJSurv.MaleLRS",
    "vcv.G.obsFemLRS.JSurv", "vcv.G.obsMaleLRS.JSurv", "vcv.G.obsJSurv.JSurv",
    "CORA.obsFemLRS.MaleLRS", "CORA.obsFemLRS.JSurv", "CORA.obsMaleLRS.JSurv",
    "h2.obs.femaleLRS", "h2.obs.maleLRS", "h2.obs.jSurv",
    "CVA.obs.femaleLRS", "CVA.obs.maleLRS", "IA.obs.femaleLRS", "IA.obs.maleLRS")


##########################

# FIGURE 3

##########################
jSurvRSPostHistbreaks <- 50

par(mfrow = c(2, 3), mar = c(9, 7, 2, 1), mgp = c(4.75,1.5,0), cex.lab = 2, cex.axis = 2)
  # first row
  trippjva <- postPlot(fitcompPost[, 20],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "V"[" A"]~" Juv. Survival")),
	xlim = c(0, 1.5), ylim = c(0, 3.2),
	ylab = "Density",
	at2 = seq(0, 3.2, 0.4))
     # Plot prior densities below y-axis maximum
     poh <- trippjva$postDensity$histogram
     # evaluate densities at some fill-in of histogram
     prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10), 2*poh$mids[1]/10)
     prd <- df(prdx / 10, df1 = 1, df2 = 4, ncp = (0^2)/10)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(2*poh$mids[1]/10 * prd)# prd area
     sprd <- prd / prda
   lines(sprd[which(sprd <= 3.3)] ~ prdx[which(sprd <= 3.3)], lwd = 3, col = "blue")
    text(x = 1.5, y = 3.2, "A", cex = 3)


  trippfva <- postPlot(fitcompPost[, 12],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "V"[" A"]~" female LRS")),
	xlim = c(0, 0.8), ylim = c(0, 25),
	ylab = "")
     # Plot prior densities below y-axis maximum
     poh <- trippfva$postDensity$histogram
     # evaluate densities at some fill-in of histogram
     prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10), 2*poh$mids[1]/10)
     prd <- df(prdx / 1000, df1 = 1, df2 = 4, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(2*poh$mids[1]/10 * prd)# prd area
     sprd <- prd / prda
   lines(sprd[which(sprd <= 25)] ~ prdx[which(sprd <= 25)], lwd = 3, col = "blue")
    text(x = 0.8, y = 25, "B", cex = 3)

  trippmva <- postPlot(fitcompPost[, 16],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "V"[" A"]~" male LRS")),
	xlim = c(0, 4), ylim = c(0, 1),
	ylab = "")
     # Plot prior densities below y-axis maximum
     poh <- trippmva$postDensity$histogram
     # evaluate densities at some fill-in of histogram
     prdx <- seq(poh$mids[1]/10, max(poh$breaks)-(poh$mids[1]/10), 2*poh$mids[1]/10)
     prd <- df(prdx / 1000, df1 = 1, df2 = 4, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(2*poh$mids[1]/10 * prd)# prd area
     sprd <- prd / prda
   lines(sprd[which(sprd <= 1)] ~ prdx[which(sprd <= 1)], lwd = 3, col = "blue")
    text(x = 4, y = 1, "C", cex = 3)

############
# second row
  postPlot(fitcompPost[, 14],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "r"[" A"]~" Juv. Survival,"~" female LRS")),
	xlim = c(-1, 1), ylim = c(0, 1.5),
	ylab = "Density",
	at2 = seq(0, 1.5, 0.5))
   # Prior is effectively flat across entire range, so to have area=1, y=0.5
   lines(rep(0.5, 100) ~ seq(-1, 1, length.out = 100), lwd = 3, col = "blue")
    text(x = 1.0, y = 1.5, "D", cex = 3)

  postPlot(fitcompPost[, 17],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "r"[" A"]~" Juv. Survival,"~" male LRS")),
	xlim = c(-1, 1), ylim = c(0, 1.5),
	ylab = "Density",
	at2 = seq(0, 1.5, 0.5))
   # Prior is effectively flat across entire range, so to have area=1, y=0.5
   lines(rep(0.5, 100) ~ seq(-1, 1, length.out = 100), lwd = 3, col = "blue")
    text(x = 1.0, y = 1.5, "E", cex = 3)

  postPlot(fitcompPost[, 13],
	histbreaks = jSurvRSPostHistbreaks,
	xlab = expression(atop("", "r"[" A"]~" female LRS,"~" male LRS")),
	xlim = c(-1, 1), ylim = c(0, 1.5),
	ylab = "",
	at2 = seq(0, 1.5, 0.5))
   # Prior is effectively flat across entire range, so to have area=1, y=0.5
   lines(rep(0.5, 100) ~ seq(-1, 1, length.out = 100), lwd = 3, col = "blue")
    text(x = 0.95, y = 1.5, "F", cex = 3)



#	XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX
################################################################################
















################################################################################


#		ANNUAL REPRODUCTIVE SUCCESS


################################################################################
rm(list = ls())
setwd("~/Dropbox/AberdeenPostdoc/Fitness/fitnessMS/dryadGithub")#setwd("<<Insert correct path here>>")

library(MCMCglmm)
library(wolakR) #<-- install from GitHub with devtools: `devtools::install_github("matthewwolak/wolakR")`
library(QGglmm)
########################################
#####################
# Read in data
#####################
arsListAinv <- read.table("./Wolak_et_al_ARS_listAinv.txt", header = TRUE)
  arsAinv <- as(sparseMatrix(i = arsListAinv[, "row"],
	j = arsListAinv[, "column"],
	x = arsListAinv[, "Ainv"],
	dims = rep(max(arsListAinv[, 1]), 2),
	dimnames = list(as.character(seq(1, max(arsListAinv[, 1]), 1)), NULL),
	symmetric = TRUE, index1 = TRUE), "dgCMatrix")

arsDat <- read.table("./Wolak_et_al_ARS_Data.txt", header = TRUE)
  arsDat$id <- arsDat$idr <- as.factor(arsDat$id)
  arsDat$year <- as.factor(arsDat$year)

#####################
## Prior
#####################
PEpriorNukp1 <- list(R = list(V = diag(2), nu = 2),
  G = list(G1 = list(V = diag(2)*0.02, nu = 3, alpha.mu = rep(0, 2), alpha.V = diag(2)*c(1000)),
	G2 = list(V = diag(2), nu = 2, alpha.mu = rep(0, 2), alpha.V = diag(2)*c(1000)),
	G3 = list(V = diag(2)*0.02, nu = 3, alpha.mu = rep(0, 2), alpha.V = diag(2)*c(1000))))




nsamp <- 5000
BURN <- 15000; THIN <- 5000
(NITT <- BURN + nsamp*THIN)

arsMod <- MCMCglmm(ARS ~ sex-1 + sex:(ageCat + f_coeff + immGG),
	random = ~ us(sex):id + idh(sex):idr + us(sex):year,
	rcov = ~ idh(sex):units,
	ginverse = list(id = arsAinv),
	prior = PEpriorNukp1,
	family = c("poisson"),
	data = arsDat,
	saveX = TRUE,  #<-- need design matrix for step below
	nitt = NITT, thin = THIN, burnin = BURN) 


# To run the above model in parallel MCMC chains, see code in the `hpc` repository
## [hpc - TorqueMoab](https://github.com/matthewwolak/hpc/tree/master/TorqueMoab/R-MCMCglmm)
## Particularly the files: `parMCMC-R.sh` and `parMCMC_postProcess.R`




###########################
# replace lower triangle with correlations
## First, save covariance version for use below
arsModCov <- arsMod 
arsMod$VCV[, 2] <- posterior.cor(arsMod$VCV[, 1:4])[, 2]
arsMod$VCV[, 8] <- posterior.cor(arsMod$VCV[, 7:10])[, 2]




###########################
# Difference between female and male inbreeding depression
## Male - female
fARSDiff <- with(arsMod, Sol[,"sexMale:f"]-Sol[,"sexFemale:f"])
postTable(fARSDiff)



###############################

# LATENT SCALE HERITABILITIES

###############################
# de Villemereuil et al. 2016
arsH2latent <- with(arsMod, as.mcmc(matrix(c(VCV[, 1] / rowSums(VCV[, c(1,5,7,11)]),
	VCV[, 4] / rowSums(VCV[, c(4,6,10,12)])),
		ncol = 2, dimnames = list(NULL, c("femaleARSh2latent", "maleARSh2latent")))))



###############################

# OBSERVED SCALE COMPUTATIONS

###############################
# Posterior prediction of the means
## Create an alternative dataset with covariates set to values below
noCov_arsDat <- arsDat
# Assume f=0 and immigrant genetic group=0 (i.e., the base population)
noCov_arsDat[, c("f_coeff", "immGG")] <- 0

# female vs. male intercept - predict for each individual then take mean of population
Xsex <- arsMod$X[, 1:2]
arsSexPredFun <- function(i){
  ipred <- predict.MCMCglmm(arsModCov, newdata = noCov_arsDat,
	type = "terms", interval = "none", it = i)

 (crossprod(Xsex, ipred) / matrix(colSums(Xsex), ncol = 1))@x
}


#XXX takes ~20 minutes for 5,000 iterations
system.time(arsSexInt <- t(sapply(seq(nrow(arsMod$VCV)), FUN = arsSexPredFun)))



arsMVparams_wrap <- function(i, itwrite = FALSE){
 if(itwrite) cat(i, " of ", nrow(arsModCov$VCV), "\n")
  Gi <- matrix(arsModCov$VCV[i, 1:4], 2, 2)
  Pi <- Gi + diag(arsModCov$VCV[i, 5:6]) + matrix(arsModCov$VCV[i, 7:10], 2, 2) + diag(arsModCov$VCV[i, 11:12])
  ovcvi <- QGmvparams(mu = arsSexInt[i,], vcv.G = Gi, vcv.P = Pi,
	models = rep("Poisson.log", 2),
	rel.acc = 0.01, width = 10, verbose = FALSE, mask = NULL)
 c(unlist(ovcvi), CORA = cov2cor(ovcvi[[3]])[1, 2], H = diag(ovcvi[[3]]) / diag(ovcvi[[2]]), CVA = 100*sqrt(diag(ovcvi[[3]])) / ovcvi[[1]], IA = diag(ovcvi[[3]]) / (ovcvi[[1]]^2))
}

#XXX XXX takes ~1.5 HOURS for 5,000 iterations
system.time(arsObs <- as.mcmc(t(sapply(seq(nrow(arsMod$VCV)),
	FUN = arsMVparams_wrap, itwrite = TRUE))))


# Full posterior distribution of the relevant model parameters
arsPost <- as.mcmc(cbind(arsMod$Sol[, 1:12], arsMod$VCV,
	arsH2latent, arsObs))
attr(arsPost, "mcpar") <- attr(arsMod$VCV, "mcpar")
  dimnames(arsPost)[[2L]][27:43] <- c("mean.obsFemaleARS", "mean.obsMaleARS",
    "vcv.P.obsFemaleARS.FemaleARS", "vcv.P.obsMaleARS.FemaleARS", "vcv.P.obsFemaleARS.MaleARS", "vcv.P.obsMaleARS.MaleARS",
    "vcv.G.obsFemaleARS.FemaleARS", "vcv.G.obsMaleARS.FemaleARS", "vcv.G.obsFemaleARS.MaleARS", "vcv.G.obsMaleARS.MaleARS",
    "CORA.obsFemaleARS.MaleARS", "h2.obsFemaleARS", "h2.obsMaleARS",
    "CVA.obsFemaleARS", "CVA.obsMaleARS", "IA.obsFemaleARS", "IA.obsMaleARS")





##########################

# FIGURE 4

##########################
arsVAx <- 0.6
fitnessPostHisBreaks <- 50

par(mfcol = c(1, 3), mar = c(8, 8, 2, 1), mgp = c(4,1.5,0), cex.lab = 2, cex.axis = 1.9)
oldparmar <- par()$mar
  postPlot(arsPost[, 13], plotHist = FALSE,
	xlab = expression("V"[" A"]~" female ARS"),
	ylab = "Density",
	xlim = c(0, arsVAx), ylim = c(0, 40),
	denslwd = 4,
	hpdlty = 0)
  text(x = arsVAx, y = 40, labels = "A", cex = 3)
    u <- par("usr")
    v <- c(grconvertX(u[1:2], "user", "ndc"),
  	   grconvertY(u[3:4], "user", "ndc"))
    v <- c(v[1]+((v[2]-v[1])/3), v[2]-((v[2]-v[1])/10), v[3]+((v[4]-v[3])/3), v[4])
    par(fig = v, new = TRUE, mar = c(0,0,0,0), cex.lab = 0.7, cex.axis = 1.4)
    arsppfva <- postPlot(arsPost[, 13],
	histbreaks = 25,
	denslwd = 5, hpdlwd = 5, meanlwd = 6,
	xlab = "",
	xlim = c(0, arsVAx*0.1), ylim = c(0, 34),
	ylab = "")
     # Plot prior densities
     poh <- arsppfva$postDensity$histogram
     prd <- df(poh$mids / 1000, df1 = 1, df2 = 3, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(diff(poh$breaks) * prd)# prd area
     sprd <- prd / prda
   lines(sprd ~ poh$mids, lwd = 5, col = "blue")

######
  par(fig = c(0.33, 0.66, 0, 1), mar = oldparmar, cex.lab = 2, cex.axis = 2, new = TRUE)
  arsppmva <- postPlot(arsPost[, 16],
	histbreaks = 25,
	xlab = expression("V"[" A"]~" male ARS"),
	xlim = c(0, arsVAx), ylim = c(0, 8),
	ylab = "")
     # Plot prior densities
     poh <- arsppmva$postDensity$histogram
     prd <- df(poh$mids / 1000, df1 = 1, df2 = 3, ncp = (0^2)/1000)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(diff(poh$breaks) * prd)# prd area
     sprd <- prd / prda
   lines(sprd ~ poh$mids, lwd = 5, col = "blue")
  text(x = arsVAx, y = 8, labels = "B", cex = 3)

######
  par(fig = c(0.66, 1, 0, 1), mar = oldparmar, cex.lab = 2, cex.axis = 2, new = TRUE)
  postPlot(arsPost[, 14],
	histbreaks = 50,
	xlab = expression("Cross-sex r"[" A"]~" ARS"),
	xlim = c(-1, 1), ylim = c(0, 1.07),
	ylab = "")
   # Prior is effectively flat across entire range, so to have area=1, y=0.5
   lines(rep(0.5, 100) ~ seq(-1, 1, length.out = 100), lwd = 3, col = "blue")
  text(x = 1, y = 1.06, labels = "C", cex = 3)



#	XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX					XXX XXX XXX XXX XXX
################################################################################
















################################################################################


#		ANNUAL ADULT SURVIVAL


################################################################################
rm(list = ls())
setwd("~/Dropbox/AberdeenPostdoc/Fitness/fitnessMS/dryadGithub")#setwd("<<Insert correct path here>>")

library(MCMCglmm)
library(wolakR) #<-- install from GitHub with devtools: `devtools::install_github("matthewwolak/wolakR")`
library(QGglmm)
########################################
#####################
# Read in data
#####################
adSurvListAinv <- read.table("./Wolak_et_al_AdultSurvival_listAinv.txt", header = TRUE)
  adSurvAinv <- as(sparseMatrix(i = adSurvListAinv[, "row"],
	j = adSurvListAinv[, "column"],
	x = adSurvListAinv[, "Ainv"],
	dims = rep(max(adSurvListAinv[, 1]), 2),
	dimnames = list(as.character(seq(1, max(adSurvListAinv[, 1]), 1)), NULL),
	symmetric = TRUE, index1 = TRUE), "dgCMatrix")

adSurvDat <- read.table("./Wolak_et_al_AdultSurvival_Data.txt", header = TRUE)
  adSurvDat$id <- adSurvDat$idr <- as.factor(adSurvDat$id)
  adSurvDat$year <- as.factor(adSurvDat$year)

#####################
## Prior
#####################
PEpriorNuk <- list(R = list(V = 1, fix = 1),
  G = list(G1 = list(V = 1, nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G2 = list(V = 1, nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10)),
	G3 = list(V = 1, nu = 1, alpha.mu = rep(0, 1), alpha.V = diag(1)*c(10))))





nsamp <- 5000
BURN <- 5000; THIN <- 2000
(NITT <- BURN + nsamp*THIN)


adSurvMod <- MCMCglmm(survive ~ sex-1 + ageCat + f_coeff + immGG,
	random = ~ id + idr + year,
	ginverse = list(id = adSurvAinv),
	data = adSurvDat,
	prior = PEpriorNuk,
	family = "categorical",
	slice = TRUE,
	saveX = TRUE,
	nitt = NITT, thin = THIN, burnin = BURN)



# To run the above model in parallel MCMC chains, see code in the `hpc` repository
## [hpc - TorqueMoab](https://github.com/matthewwolak/hpc/tree/master/TorqueMoab/R-MCMCglmm)
## Particularly the files: `parMCMC-R.sh` and `parMCMC_postProcess.R`






###############################

# LATENT SCALE HERITABILITIES

###############################
# de Villemereuil et al. 2016
adSurvH2latent <- with(adSurvMod, as.mcmc(matrix(c(VCV[, 1] / rowSums(VCV)),
		ncol = 1, dimnames = list(NULL, c("adSurvH2latent")))))




###############################

# OBSERVED SCALE COMPUTATIONS

###############################
# Posterior prediction of the means
## Create an alternative dataset with covariates set to values below
noCov_adSurvDat <- adSurvDat
# Assume f=0 and immigrant genetic group=0 (i.e., the base population)
noCov_adSurvDat[, c("f_coeff", "immGG")] <- 0

adSurvQGparams_wrap <- function(i){
 cat(i, "\n")
  yhat <- predict.MCMCglmm(adSurvMod, newdata = noCov_adSurvDat,
	type = "terms", interval = "none", it = i)
  QGparams(predict = yhat, var.a = adSurvMod$VCV[i, 1],
	var.p = sum(adSurvMod$VCV[i, ]),
	model = "binom1.logit", verbose = FALSE)
}


#XXX takes ~76.9 minutes for 5,000 iterations
system.time(adSurvObs <- do.call(rbind, lapply(seq(nrow(adSurvMod$VCV)), FUN = adSurvQGparams_wrap)))



# Full posterior distribution of the relevant model parameters
adSurvPost <- as.mcmc(cbind(adSurvMod$Sol[, 1:7], adSurvMod$VCV,
	adSurvH2latent, adSurvObs))
attr(adSurvPost, "mcpar") <- attr(adSurvMod$VCV, "mcpar")
  dimnames(adSurvPost)[[2L]][13:16] <- c("mean.obsAdSurv", 
	"var.obsAdSurv", "var.a.obsAdSurv", "h2.obsAdSurv")



##########################

# FIGURE 5C

##########################

par(mar = c(7.5, 7, 2, 2), mgp = c(4.5,1.5,0), cex.lab = 2, cex.axis = 2)
  adsurvppva <- postPlot(adSurvPost[, "id"], bw = 0.008,
	histbreaks = 50,
	xlab = expression(paste("V"[" A"], " Annual adult survival")),
	xlim = c(0, 0.3), ylim = c(0, 42),
	ylab = "Density")
     # Plot prior densities
     poh <- adsurvppva$postDensity$histogram
     prd <- df(poh$mids / 10, df1 = 1, df2 = 1, ncp = (0^2)/10)
     # scale so total of: density * width of histogram bar, summed over all bars = 1
       prda <- sum(diff(poh$breaks) * prd)# prd area
     sprd <- prd / prda
   lines(sprd ~ poh$mids, lwd = 5, col = "blue")
  text(x = 0.3, y = 42.5, labels = "C", cex = 3)





