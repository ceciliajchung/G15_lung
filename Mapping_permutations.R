library("DOQTL")
library("QTLRel")
source("scanone_mi.R")
set.seed(1234)

# Meta data
# : 'g15.fam' with cage information as random variable
# : sex, cage information as fixed variables
# 'pheno.data'
# : data frame of phenotypes
# : column 1 is sex
# Sample names or IDs must be matching the rownames in all phenotypic and genotypic data.


## Perform Mapping ####
filename.lod <- vector();
filename.coef <- vector();

for(i in 2:ncol(pheno.data)){ #column 1 is sex 
  G15.qtl <- NULL;
  G15.qtl <- scanone_mi(pheno.data.miseq, pheno.col=i, probs=array.G.mice.g15, K=kin.G.mice.g15, snps=snps.G.mice, addcovar=covar);
  a <- NULL;
  a <- rbind(G15.qtl$lod$A,G15.qtl$lod$X);
  b <- NULL;
  b <- rbind(G15.qtl$coef$A,G15.qtl$coef$X);
  filename.lod[i] <- paste(colnames(pheno.data)[i],"lod_G15.csv",sep="_");
  filename.coef[i] <- paste(colnames(pheno.data)[i],"coef_G15.csv",sep="_");
  write.csv(a,file=filename.lod[i],quote=F,row.names=T);
  write.csv(b,file=filename.coef[i],quote=F,row.names=T);
  filename.lod <- c(filename.lod);
  filename.coef <- c(filename.coef);      
}


## Perform permutations ####
permute=10000;
a.mat <- data.frame(pheno.data);
rownames(a.mat) <- rownames(pheno.data)
names(a.mat) <- colnames(pheno.data)
perm.mat <- vector();
thresholds <- data.frame();

for (j in 1:ncol(a.mat)){
  for(i in 1:permute){
    a.perm <- NULL;
    a.temp <- apply(a.mat,2,sample); # shuffle without replacement
    a.temp <- as.data.frame(a.temp);
    a.temp <- apply(a.temp,2,as.numeric); # sex can be either numeric or factor  
    rownames(a.temp) <- rownames(a.mat); # force phenotypes to have the same row.names
    a.perm <- scanone_mi(a.temp,pheno.col=j,probs=array.G.mice.g15,K=kin.G.mice.g15,snps=snps.G.mice,addcovar=covar);
    print(max(c(a.perm$lod$A$lod,a.perm$lod$X$lod)));
    perm.mat[i] <- as.numeric(max(c(a.perm$lod$A$lod,a.perm$lod$X$lod)));
    perm.mat <- c(perm.mat);
    cat(paste(i,"..... permutation done\n",sep=""));
  }
  thr[j] <- round(perm.mat, digits=4);
  thr.95[j] <- round(quantile(perm.mat, probs=0.95),digits = 4); # Significant association
  thr.90[j] <- round(quantile(perm.mat, probs=0.90),digits = 4); # Suggestive association
  data <- cbind(colnames(a.mat)[j],thr.95[j],thr.90[j]);
  thresholds <-  rbind(thresholds,data);
}

colnames(thresholds) <- c("trait","thr.95","thr.90");
rownames(thresholds) <- colnames(a.mat)[-1];
write.csv(thresholds,"thresholds.csv");


# END
