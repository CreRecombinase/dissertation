rm(list = ls())
source("TADA.R")

### to guarantee the same answer for this example use the following line, do not use this in your actual analysis
set.seed(100)

### read mutation data
tada.file <- "TADA_demo_counts_de-novo_only.txt"
tada.data <- read.table(tada.file, header = T)

### specify the number of families included in the analysis
n.family <- 4500
n <- data.frame(dn = n.family, ca = NA, cn = NA)
sample.counts <- list(cls1 = n, cls2 = n)

# create the mutational data used by TADA-Denovo
cls1.counts <- data.frame(dn = tada.data$dn.cls1, ca = NA, cn = NA)
rownames(cls1.counts) <- tada.data$gene.id
cls2.counts <- data.frame(dn = tada.data$dn.cls2, ca = NA, cn = NA)
rownames(cls2.counts) <- tada.data$gene.id
tada.counts <- list(cls1 = cls1.counts, cls2 = cls2.counts)

### set up mutation rates
mu <- data.frame(cls1 = tada.data$mut.cls1, cls2 = tada.data$mut.cls2)

### set up denovo only TRUE/FALSE, here we do not want to restrict ourselves to de novo only analyses
denovo.only <- data.frame(cls1 = TRUE, cls2 = TRUE)

# set up parameters
cls1 <- data.frame(gamma.mean.dn = 20.0, beta.dn = 1, gamma.mean.CC = NA, beta.CC = NA, rho1 = NA, nu1 = NA, rho0 = NA, nu0 = NA)
cls2 <- data.frame(gamma.mean.dn = 4.7, beta.dn = 1, gamma.mean.CC = NA, beta.CC = NA, rho1 = NA, nu1 = NA, rho0 = NA, nu0 = NA)
hyperpar <- list(cls1 = cls1, cls2 = cls2)

# running TADA-Denovo
re.TADA <- do.call(cbind.data.frame, TADA(tada.counts = tada.counts, sample.counts = sample.counts, mu = mu, hyperpar = hyperpar, denovo.only = denovo.only))

# Bayesian FDR control
re.TADA$qval <- Bayesian.FDR(re.TADA$BF.total, pi0 = 0.95)

# run permutation to get the null distributions to use for calculating p-values for TADA
re.TADA.null <- do.call(cbind.data.frame, TADAnull(tada.counts = tada.counts, sample.counts = sample.counts, mu = mu, hyperpar = hyperpar, denovo.only = denovo.only, nrep = 100))
re.TADA$pval <- bayesFactor.pvalue(re.TADA$BF.total, re.TADA.null$BFnull.total)

# display top 10 genes based on BF.total
re.TADA[order(-re.TADA$BF.total)[1:10], ]
