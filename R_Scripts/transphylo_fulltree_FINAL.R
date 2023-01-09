### Library ---------------------------------------------------------------------

library(castor)
library(rlist)
library(TransPhylo)
#devtools::install_github('xavierdidelot/Transphylo')
library(ape)
library(coda)
library(stringr)
library(lubridate)
library(plyr)
library(tibble)
library(reshape2)
library(igraph)
library(rgdal) 
library(sp)
library(lubridate)
library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(ggridges)
library(wesanderson)





### Load files ------------------------------------------------------------------

#MCC tree
phyMCC<-read.nexus("")


### Transphylo MCC tree ---------------------------------------------------------

plot(phyMCC, show.tip.label=F)
axisPhylo(backward = F)


### Dated tree ------------------------------------------------------------------
# Gives the date of the last sample otherwise all values of dates are relative to each other 

ptreeMCC <- ptreeFromPhylo(phyMCC,dateLastSample=2017.351) 
plot(ptreeMCC, show.tip.label=F)


### Infer the transmission tree -------------------------------------------------

w.shape=2           # alpha of gamma distribution (generation time distribution)  (from primary to secondary infection)
w.scale=0.025          # beta of gamma distribution 
ws.shape=2          # alpha of gamma distribution (sampling time distribution)  (from primary infection to sampling)
ws.scale=0.025          # beta of gamma distribution   
dateT=Inf
iter=10000000           # MCMC iterations
thin=10000            # MCMC thinning interval between two sampled iterations

res<-inferTTree(ptreeMCC, w.shape=w.shape, w.scale=w.scale, ws.shape=ws.shape, ws.scale=ws.scale, dateT=dateT, 
                updatePi=T, updateOff.p = F,
                mcmcIterations=iter, thinning = thin)
#set.seed
save(res, file = "")


### Checking convergence --------------------------------------------------------

plot(res)


### Obtain the effective sample size (ESS) of paramaters ------------------------

mcmc=convertToCoda(res)
effectiveSize(mcmc)



### Plot the sampling proportion distribution -----------------------------------

hist(sapply(res,function(x) x$pi), xlab='Sampling proportion', main=NULL)


### Find the most representative (aka medoid) colored tree ----------------------
#Returns the transmission tree from the posterior that is least different from 
#all others according to a well-defined distance metric

med=medTTree(res, burnin=0.5)

ctree=medTTree(res)
plotCTree(ctree,showLabels = F,cols = c('black','black'))

#This function extract all the transmission events from a colored tree
extractTransEvents=function(ctree) {
  nam=ctree$nam
  ctree=ctree$ctree
  w=which(ctree[,2]>0&ctree[,3]==0)
  out=c()
  nsam=length(which(ctree[,2]==0&ctree[,3]==0))
  for (i in 1:length(w)) {
    cur=w[i]
    while (ctree[cur,2]>0) cur=ctree[cur,2]
    time=ctree[cur,1]-ctree[w[i],1]
    out=rbind(out,c(nam[cur],time))
  }
  return(out)
}
ete=extractTransEvents(ctree)


### Plot the corresponding transmission tree ------------------------------------

ttree=extractTTree(med)

plot(ttree,type='summarised', w.shape, w.scale, showLabels=F)


## Plot of sampled and unsampled cases over time -------------------------------

a=getIncidentCases(res, numBins=50, show.plot = T)

source("Script/Transphylo_R/my_getIncidentCases.R")

my_getIncidentCases(res, numBins=50, burnin=0.5, show.plot = T)


### Distribution of realised sampling times -------------------------------------

a=getSamplingTimeDist(res, maxi=0.2, numBins= 20, show.plot = T)


### Distribution of infection time for the individuals labelled i ---------------

a=getInfectionTimeDist(res, k=i, burnin=0.5, show.plot = T)

getInfectionTimeDist <- function(record,burnin=0.5,k,numBins=10,show.plot=F) {
  record=record[max(1,round(length(record)*burnin)):length(record)]
  for (i in 1:length(record)) record[[i]]=extractTTree(record[[i]]$ctree)
  times=matrix(NA,length(k),length(record))
  for (i in 1:length(k)) for (j in 1:length(record)) {
    ii=which(record[[j]]$nam==k[i])
    times[i,j]=record[[j]]$ttree[ii,1]
  }
  if (length(k)==1) times=as.vector(times)
  return(times)
}


### Distribution of realised generation times -----------------------------------

a=getGenerationTimeDist(res, maxi=0.2, numBins= 20, show.plot = F) 


### Matrix of probability of direct transmission for all pairs------------------ 

mat=computeMatWIW(res, burnin=0.5)

lattice::levelplot(mat,xlab='',ylab='',
                   scales=list(x=list(rot=90, cex=0.5),
                               y=list(cex=0.5)))
lattice::levelplot(mat+t(mat),xlab='',ylab='',
                   scales=list(x=list(rot=90, cex=0.5),
                               y=list(cex=0.5)))

