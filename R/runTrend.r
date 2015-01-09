# Run Cochran-Armitage trend test on each variant of a genotype matrix in parallel
runTrendCl<-function(g, fct, values=0:max(g, na.rm=TRUE), cov=NA, num.threads=8) {
  # g             A matrix of genotype calls, each row is a variant and each column is a sample
  # fct           A factor variable or integer vector corresponding to sample grouping, much match the columns of genotype matrix
  # values        Accepted values of genotype calls
  # cov           A matrix whose rows matching the samples, each column is a covariant of test
  # num.threads   Number of threads
  library(snow);
  
  num.threads<-max(num.threads, 1);
  num.threads<-min(num.threads, nrow(g));
  
  if (is.null(cov)) cov<-NA;
  
  # split matrix
  split<-rep(1:num.threads, each=ceiling(nrow(g)/num.threads))[1:nrow(g)];
  gs<-lapply(1:num.threads, function(i) g[split==i, , drop=FALSE]);
    
  cl<-makeCluster(num.threads, type='SOCK');
  ps<-clusterApplyLB(cl, gs, runTrend, fct=fct, values=values, cov=cov);  
  try(stopCluster(cl));
  
  do.call('c', ps);
}


# Run Cochran-Armitage trend test on each variant of a genotype matrix
runTrend<-function(g, fct, values=0:max(g, na.rm=TRUE), cov=NA) {
  # g             A matrix of genotype calls, each row is a variant and each column is a sample
  # fct           A factor variable or integer vector corresponding to sample grouping, much match the columns of genotype matrix
  # values        Accepted values of genotype calls
  # cov           A matrix whose rows matching the samples, each column is a covariant of test
  
  library(coin);
  
  fct<-as.integer(fct);
  
  # limited to accepted values
  g[! (g %in% values)]<-NA;
  
  # create fomula
  fm<-'fct ~ d';
  if (class(cov) == 'matrix' | class(cov) == 'data.frame') if (nrow(cov) == length(fct)) {
    fm<-paste(fm, paste(paste('+ ', 'cov[, ', 1:ncol(cov), ']', sep=''), collapse=' '));
  }
  fm<-as.formula(fm);
  
  # run test
  p<-apply(g, 1, function(d) pvalue(independence_test(fm, teststat = "quad"))[[1]]);
  p[is.na(p)]<-1;
  names(p)<-rownames(g); 
  
  p;
}