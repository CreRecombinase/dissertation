rm(list = ls())
source("TADA.R")

### to guarantee the same answer for this example
#use the following line, do not use this in your actual analysis
set.seed(100)

### read mutation data
tada.file <- "TADA_demo_counts_de-novo_and_inherited.txt"
tada.data <- read.table(tada.file, header = T)

### specify the number of families and
#the number of cases and control samples included in the analysis
n.family <- 4500
n.case <- 1000
n.ctrl <- 3000
n <- data.frame(dn = n.family, ca = n.case + n.family, cn = n.ctrl + n.family)
sample.counts <- list(cls1 = n, cls2 = n)

# create the mutational data used by TADA
cls1.counts <- data.frame(dn = tada.data$dn.cls1,
                          ca = tada.data$trans.cls1 + tada.data$case.cls1,
                          cn = tada.data$ntrans.cls1 + tada.data$ctrl.cls1)
rownames(cls1.counts) <- tada.data$gene.id
cls2.counts <- data.frame(dn = tada.data$dn.cls2,
                          ca = tada.data$trans.cls2 + tada.data$case.cls2,
                          cn = tada.data$ntrans.cls2 + tada.data$ctrl.cls2)
rownames(cls2.counts) <- tada.data$gene.id
tada.counts <- list(cls1 = cls1.counts, cls2 = cls2.counts)

### set up mutation rates
mu <- data.frame(cls1 = tada.data$mut.cls1, cls2 = tada.data$mut.cls2)

### set up denovo only TRUE/FALSE,
#here we do not want to restrict ourselves to de novo only analyses
denovo.only <- data.frame(cls1 = FALSE, cls2 = FALSE)

# set up parameters
cls1 <- data.frame(gamma.mean.dn = 20.0,
                   beta.dn = 1,
                   gamma.mean.CC = 2.3,
                   beta.CC = 4.0,
                   rho1 = 0.1,
                   nu1 = 100,
                   rho0 = 0.1,
                   nu0 = 100)
cls2 <- data.frame(gamma.mean.dn = 4.7, beta.dn = 1, gamma.mean.CC = 1.0, beta.CC = 1000, rho1 = 0.15, nu1 = 100, rho0 = 0.15, nu0 = 100)
hyperpar <- list(cls1 = cls1, cls2 = cls2)

# running TADA
re.TADA <- do.call(cbind.data.frame, TADA(tada.counts = tada.counts, sample.counts = sample.counts, mu = mu, hyperpar = hyperpar, denovo.only = denovo.only))

# Bayesian FDR control
re.TADA$qval <- Bayesian.FDR(re.TADA$BF.total, pi0 = 0.95)

# run permutation to get the null distributions to use for calculating p-values for TADA
re.TADA.null <- do.call(cbind.data.frame, TADAnull(tada.counts = tada.counts, sample.counts = sample.counts, mu = mu, hyperpar = hyperpar, denovo.only = denovo.only, nrep = 100))
re.TADA$pval <- bayesFactor.pvalue(re.TADA$BF.total, re.TADA.null$BFnull.total)

# top 10 genes based on BF.total
re.TADA[order(-re.TADA$BF.total)[1:10], ]
