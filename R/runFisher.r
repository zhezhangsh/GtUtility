# run Fisher test in parallel on each variant (allelic, recessive, dominant tests), based on counts of each genotype type
runFisherCl<-function(ct, num.threads=8) {
  # ct    A matrix of 6 columns in orders of AA0, AB0, BB0, AA1, AB1, and BB1, each is the number of samples having the genotype, 0 and 1 represent control and case groups
  
  library(snow);
  
  num.threads<-max(num.threads, 1);
  num.threads<-min(num.threads, nrow(ct));
  
  split<-rep(1:num.threads, each=ceiling(nrow(ct)/num.threads))[1:nrow(ct)];
  cts<-lapply(1:num.threads, function(i) ct[split==i, , drop=FALSE]);
  
  cl<-makeCluster(num.threads, type='SOCK');
  stat<-clusterApplyLB(cl, cts, runFisher);  
  try(stopCluster(cl));
  
  do.call('rbind', stat);
}


# run Fisher test on each variant (allelic, recessive, dominant tests), based on counts of each genotype type
runFisher<-function(ct) {
  # ct    A matrix of 6 columns in orders of AA0, AB0, BB0, AA1, AB1, and BB1, each is the number of samples having the genotype, 0 and 1 represent control and case groups
  
  runF<-function(cnt) {
    # prepare matrixes
    mtr<-list(
      allelic=cbind(control=as.vector(c(cnt[1]*2+cnt[2], cnt[3]*2+cnt[2])), case=as.vector(as.vector(c(cnt[4]*2+cnt[5], cnt[6]*2+cnt[5])))), 
      recessive=cbind(control=as.vector(c(cnt[1]+cnt[2], cnt[3])), case=as.vector(c(cnt[4]+cnt[5], cnt[6]))),
      dominant=cbind(control=as.vector(c(cnt[1], cnt[2]+cnt[3])), case=as.vector(c(cnt[4], cnt[5]+cnt[6])))
    )
    
    # Run fisher test
    stat<-as.vector(sapply(mtr, function(mtr) {
      fsh<-fisher.test(mtr); 
      c(as.vector(mtr), fsh$estimate[[1]], fsh$p.value[[1]]); 
    }));
    names(stat)<-paste(rep(c('Allelic', 'Recessive', 'Dominant'), each=6), rep(c('A0', 'B0', 'A1', 'B1', 'OR', 'p')), sep='_');
    names(cnt)<-c('AA0', 'AB0', 'BB0', 'AA1', 'AB1', 'BB1');
    
    c(Total=sum(cnt), cnt, stat);
  };
 
  t(apply(ct, 1, runF));
}