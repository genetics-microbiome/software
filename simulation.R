
require(dirmult) # for simulating from DM distribution. 
require(cluster) # for Partitioning Around medoid (PAM) clustering 
require(GUniFrac) # contains throat data (real data) to be used to obtain realistic sim params. 
require(picante) # to prune phylogenetic tree. 
require(ape) # to simulate phylogenetic tree for the OTUs. 
require(MASS) # simulate from multivariate Gaussian
library(gplots) # Needed for plotting (heatmap)
library(RColorBrewer) # Colors in plotting, say heatmaps

## Step 1: Simulating microbiome meta data.
# ----------------------------------------------------------------------------
# ----------------------------------
# Step 1(a): simulate OTU abundance
# ----------------------------------
n <- 500 # number of samples to simulate
p <- 200 # number of OTUs to simulate
conc.param = 24.21 # concentrating parameter.

# Real data for simulation parameters
data(throat.otu.tab) # from GUniFrac package
data(throat.tree)    # from GUniFrac package

# Only p [in 200, 400, 800] most abundant taxa will be used for parameter estimation.
otu.pruned <- throat.otu.tab[,order(colSums(throat.otu.tab), decreasing = TRUE)[1:p]]
DM.params <- dirmult(data=otu.pruned,epsilon=10^(-4), trace=TRUE) #obtain 'pi' and 'theta'.

# simulate OTU counts
otu.counts <- simPop(J=n, # number of samples to simulate
                     n=p,  # number of OTUs to simulate (this n is not sample size!!)
                     pi=DM.params$pi,
                     theta=DM.params$theta)$data
#dim(otu.counts)

# OTU counts normalzed by read counts, sampled from negative binom with mean 1000 and dispersion 25:
otu.abun <- apply(otu.counts,2,function(x){return(x/rnbinom(n=1,mu=5000,size=25) )})
colnames(otu.abun)=paste0("OTU",1:p)


# ------------------------------------------------------------------------------
# prune OTU tree to the considered-OTUs. 
otu.tree <- picante::prune.sample(throat.otu.tab[,colnames(otu.pruned)], throat.tree) #


# ------------------------------------------------------------------------------
# Use PAM algorithm, to partition OTUs into clusters based on patristic distances.
# Partition into k=50
otu.patr.dist.mat <- cophenetic(otu.tree) # compute patristic distances.
otu.abun.clust    <- pam(scale(otu.patr.dist.mat), k=50, metric = "manhattan", 
                         cluster.only = FALSE) # k=clusters

clus.abun=order(table(otu.abun.clust$clustering), decreasing = TRUE) #order from most abundant cluster
otu.in.cluster=otu.abun[, otu.abun.clust$clustering==clus.abun[1]] # cluster most abundant cluster
#dim(otu.in.cluster) #  OTUs in the most abundant cluster.
num.assoc.otu = 10 # set number of associated ('causal') OTUs
assoc.otus = otu.in.cluster[,order(colMeans(otu.in.cluster), decreasing = TRUE)][,1:num.assoc.otu]
#dim(assoc.otus)
print(assoc.otus) # get names of the causal OTUs.


# ------------------------------------------------------------------------------
#Simulate random effects, Psi. Use polya urn scheme.
polya_urn_model = function(base_color_distribution, num_balls, alpha) {
  balls = c()
  
  for (i in 1:num_balls) {
    if (runif(1) < alpha / (alpha + length(balls))) {
      # Add a new ball color.
      new_color = base_color_distribution()
      balls = c(balls, new_color)
    } else {
      # Pick out a ball from the urn, and add back a
      # ball of the same color.
      ball = balls[sample(1:length(balls), 1)]
      balls = c(balls, ball)
    }
  }
  
  balls
} 

Psi <- polya_urn_model(function() rnorm(1), num_balls = n, alpha = conc.param)

# ------------------------------------------------------------------------------
## Simulate phenotype.
assoc.beta = c(0.1,0.1,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0) # effects for associated OTUs.


names(assoc.beta) = colnames(assoc.otus) # locations of known true effects.

Beta = rep(0,p) 
names(Beta) <- paste0("OTU",1:p)
Beta[names(assoc.beta)] <- assoc.beta # set effects for causal OTUs, effects for non-causal ones set to 0.
#replace(true.beta, names(assoc.beta), assoc.beta) # alternative way.

z = otu.abun %*% Beta  + Psi
pr = exp(z)/(1+exp(z)) # 
repeat{
  y=rbinom(n, size = 1, prob = pr) #
  if(length(which(y==1))==0.5*n) break # simulate equal number of cases and controls.
}



out2<-phyloDPM(y,otu.abun,mTree=otu.tree,mu=0,alpha.par=NULL,Nsim=10000,post.plots=TRUE)















