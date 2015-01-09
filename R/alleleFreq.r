# Summarize allele frequency of variants
alleleFreq<-function(gt) {
  # gt  Genotype call matrix, only 3 values will be counted: 0 = A/A; 1 = A/B; and 2 = B/B
  ct<-sapply(0:2, function(i) apply(gt, 1, function(gt) length(gt[!is.na(gt) & gt == i])));
  
  N<-rowSums(ct); 
  
  frq<-apply(ct, 2, function(x) 100*x/N);
  colnames(frq)<-c('AA%', 'AB%', 'BB%');
  
  frq<-cbind(Total=N, 'MAF(%)'=(ct[,2]+2*ct[,3])/2/N*100, frq);
  
  frq<-cbind(frq, Count_A=ct[,1]*2+ct[,2], Count_B=ct[,2]+ct[,3]*2);
  
  frq[N==0,]<-0;
  
  frq;
}