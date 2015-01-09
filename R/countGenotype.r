# Count occurance of each genotype
countGenotypeCl<-function(g, samples=NA, max=2, num.threads=8) {
  # g         A genotype matrix
  # samples   A subset of samples/columns if specified
  # max       Maximum value of genotype
  
  library(snow);
  
  # sample subset
  samples<-samples[samples %in% colnames(g)];
  if (length(samples) == 0) samples<-colnames(g); # use all samples if not specified
  g<-g[, samples, drop=FALSE];
  
  # split matrix for parrallel
  sz<-ceiling(nrow(g)/num.threads);
  ind<-seq(1, nrow(g), sz);
  ind<-c(ind, nrow(g)+1);
  ind<-cbind(ind[-length(ind)], ind[-1]-1);
  gs<-apply(ind, 1, function(i) g[i[1]:i[2], , drop=FALSE]);
  
  cl<-makeCluster(num.threads, type='SOCK');
  ct<-clusterApplyLB(cl, gs, countGenotype, max=max);  
  try(stopCluster(cl));
  
  do.call('rbind', ct);
}


countGenotype<-function(g, samples=NA, max=2) {
  # g         A genotype matrix
  # samples   A subset of samples/columns if specified
  # max       Maximum value of genotype
  
  # sample subset
  samples<-samples[samples %in% colnames(g)];
  if (length(samples) == 0) samples<-colnames(g); # use all samples if not specified
  
  g<-g[, samples, drop=FALSE];
  
  if (is.na(max) | max<0) max<-max(g, na.rm=TRUE);
  
  ct<-sapply(0:max, function(i) apply(g, 1, function(g) {
    g<-g[!is.na(g) & g>=0 & g<=max];
    length(g[g==i]);
  }));
  
  colnames(ct)<-paste('count_', 0:max, sep='');
  
  ct;
}