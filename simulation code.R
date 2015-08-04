
############################################################################
# Code to reproduce simulations from Waples & Ward 
# author: Eric Ward, eric.ward@noaa.gov  July 1 2015
# 
# Input for this code is EASYPOP files. For all scenarios, we used 250 f/m 
# from 2 popualtions. Migration = 0.0005. 
# 
############################################################################
library(adegenet)
library(genetics)

FST = "0.05" # We explored 3 levels of FST, set to 0.025, 0.05, 0.1

setwd(paste("/users/eric.ward/Wahlund-Waples-Ward/simulated data/Fst_",FST,sep=""))

# For specified value of FST, we ran 2 microsatellie scenarios (10, 20 loci) and 1
# SNP scenario. This loop for iii iterates over all 3.
for(iii in 1:3) {

if(iii==1) {
	if(FST=="0.025") run = "allele10_55"
	if(FST=="0.05") run = "allele10_125"
	if(FST=="0.1") run = "allele10_375"
}
if(iii==2) {
	if(FST=="0.025") run = "allele20_55"	
	if(FST=="0.05") run = "allele20_125"	
	if(FST=="0.1") run = "allele20_375"
}
if(iii==3) {
	if(FST=="0.025") run = "snp_50"
	if(FST=="0.05") run = "snp_140"	
	if(FST=="0.1") run = "snp_300"
}

SNP = FALSE
if(substr(run,1,3) == "snp") SNP = TRUE

library(hierfstat)
detach(package:hierfstat)

datfile = paste(run,".dat",sep="")
A <- read.fstat(datfile)

N = length(A@ind.names) # number of animals from both pops
g2h = genind2hierfstat(A,pop=c(rep(1,N/2), rep(2,N/2)))
gen =genind2genotype(A) # create genotype object

library(hierfstat)
genotypes = g2h[,-1]
LOCI = dim(genotypes)[2]
STATES = c(10, 20, 2)[iii]

library(pegas)
gt = genind2loci(A)

# Score the genotypes so they can be used to calculate r^2
# Calculate correlation over each allelic state. Average across states
g1 = matrix(0,dim(gen)[1],dim(gen)[2])
g2 = matrix(0,dim(gen)[1],dim(gen)[2])
gMat = array(0,dim=c(dim(gen)[1],dim(gen)[2],STATES))
for(i in 1:dim(gen)[1]) {
  for(j in 1:dim(gen)[2]) {
    g1[i,j] = as.numeric(strsplit(as.character(gen[i,j]),"/")[[1]][1])
    g2[i,j] = as.numeric(strsplit(as.character(gen[i,j]),"/")[[1]][2]) 
    gMat[i,j,g1[i,j]] = gMat[i,j,g1[i,j]] + 1
    gMat[i,j,g2[i,j]] = gMat[i,j,g2[i,j]] + 1
  }
}

# We want to iterate over (1) sample size of individuals, (2) loci, (3) mixing fractions, (4) number of loci
scenarios = expand.grid('samples' = c(10, 20,50,100,200,400), 'mix' = c(0.1,0.3,0.5),'nloci' = c(10,20))
if(SNP == TRUE) scenarios = expand.grid('samples' = c(10, 20,50,100,200,400), 'mix' = c(0.1,0.3,0.5),'nloci' = c(100,200))

# This loops over all the scenarios created with expand.grid()
for(ii in 1:dim(scenarios)[1]) {

mixPop1 = scenarios$mix[ii]# mixing fraction
nSampled = scenarios$samples[ii] # total animals in the mixture (both populations)
sampledLoci = scenarios$nloci[ii] # number of loci to sample 

weights = c(mixPop1, 1-mixPop1) 
sampledAnimals = ceiling(nSampled * weights) # numbers of animals to sample from each population

# Create a large number ~ 1000 replicates for this scenario
SIMS = 1000
output = list()
for(i in 1:SIMS) {

   # Generate constant samples of same size from pop 1 and 2
   # draw N/2 animals from each subpopulation
   animals = sample(seq(1,N/2), size = nSampled, replace=F)
   animals2 = sample(seq(N/2+1,N), size = nSampled, replace=F)    
   
   # Generate random sample of loci from the populations -- same loci, without replacement
   loci = sample(seq(1,LOCI), size = sampledLoci, replace=F)

   # calculate Fst following Nei & Chesser 1983
   # Fst is calculated based on mixture of 2 pops. 
   subft = g2h[c(animals,animals2),c(1,loci+1)] # pop ID in first column
   
   output[[i]] = matrix(NA, sampledLoci, 7) # matrix to store results
   rownames(output[[i]]) = names(subft)[-1]
   colnames(output[[i]]) = c("Fst","Fis","Zh1","Zh2","delta","Zh3","CorFst1Fst2_R2")
   
   # Calculate FST for each locus
   Fstvec = basic.stats(subft)$perloc$Fst
   output[[i]][,1] = Fstvec

   # Calculate pairwise products of Fst across loci
   Fst1_Fst2 = Fstvec%o%Fstvec # 

   # Generate 2nd sample for Fis / R2 calculation. This sample is the same as that 
   # specified in expand.grid() above. 
   animals = sample(seq(1,N/2), size = sampledAnimals[1], replace=F)
   animals2 = sample(seq(N/2+1,N), size = sampledAnimals[2], replace=F)  
   # Subset genotype data for these animals
   subft = g2h[c(animals,animals2),c(1,loci+1)]
   
   # Calculate gene diversity, using Hs by population
   div1 = basic.stats(subft)$Hs[,1]
   div2 = basic.stats(subft)$Hs[,2]
   output[[i]][,5] = (div1 - div2)/(div1+div2) # gene diveristy
   
   # Change popID to the same for all animals, because Fis is calculated 
   # including all animals as same pop
   subft$pop=1
   rownames(subft) = paste(seq(1,nSampled)) # character vector of animal ID

   # Calculate Fis for each locus
   output[[i]][,2] = basic.stats(subft)$perloc$Fis
   
   # Calculate pieces of theoretical expectations (L. Zhivotovsky 2015)
   output[[i]][,3] = log(1/output[[i]][,2] - 1)
   output[[i]][,4] = log(1/output[[i]][,1] - 1)
   output[[i]][,6] = log((1 - output[[i]][,5]*(1 - 2*scenarios$mix[ii]))/(4*scenarios$mix[ii]*(1-scenarios$mix[ii])))
   
   # Calculate correlation over each allelic state. Average across states
   # there are 10*9/2 = 45 pairs of loci and 45 Fst1*Fst2 values.  consider only the comparisons of locus1 and locus2
   # there are 10x10 = 100 comparisons of alleles at locus 1 vs locus 2
   # recode to 0,1,2 and calculate r^2 for each of the 100 comparisons
   # mean of the 100 r^2 values is the mean r^2 for that pair of loci
   # repeat for the other 44 comparisons of loci and compute the correlation of r^2 and Fst1*Fst2
   subgMat = gMat[c(animals,animals2),loci,]
   corMat = matrix(0, length(loci), length(loci))
   for(iii in 1:length(loci)) {
     for(jjj in 1:length(loci)) {
       # pairwise comparaison of all ellelic correlatinos 
       corMat[iii,jjj] = mean(cor(cbind(subgMat[,iii,],subgMat[,jjj,]))[1:STATES,((STATES+1):(2*STATES))]^2, na.rm=T)
     }    
   }
   
   # Calculate correlation between pairwise Fst products and R2
   output[[i]][,7]=cor(corMat[lower.tri(corMat)], Fst1_Fst2[lower.tri(Fst1_Fst2)], use = "pairwise.complete.obs")

} # end i sims

# save this output list to a file
save(output,file = paste("output/",run,"_",paste("mix",mixPop1,sep=""),paste("_N",nSampled,sep=""),paste("_n",sampledLoci,sep=""),".Rdata",sep=""))
}

} # End iii

