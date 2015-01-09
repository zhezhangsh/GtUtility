# recode genotype calls to callapse alternative variants when there are 2 or more of them
recodeCalls<-function(gt, coded) {
  # gt      The original genotype call matrix
  # coded   Coded value of each allele, rows must match rows of genotype call matrix, and columns match columns of alleles matrix except the first one (the reference allele)
  
  combs<-lapply(1:(1+ncol(coded)), function(i) cbind(0:(i-1), rep(i-1, length(i)))); 
  combs<-do.call('rbind', combs);
  c1<-combs[,1];
  c2<-combs[,2];
  
  A<-B<-gt;
  A[!is.na(A)&A>=0&A<nrow(combs)]<-c1[(A+1)[!is.na(A)&A>=0&A<nrow(combs)]];
  B[!is.na(B)&B>=0&B<nrow(combs)]<-c2[(B+1)[!is.na(B)&B>=0&B<nrow(combs)]];
  
  for (i in 1:nrow(gt)) {
    c<-coded[i,];
    A[i, !is.na(A[i,]) & A[i, ]>0]<-c[A[i, !is.na(A[i,]) & A[i, ]>0]];
    B[i, !is.na(B[i,]) & B[i, ]>0]<-c[B[i, !is.na(B[i,]) & B[i, ]>0]];    
  }
  
  a<-pmin(A, B);
  b<-pmax(A, B);
  
  combs<-lapply(1:(1+max(b, na.rm=TRUE)), function(i) cbind(0:(i-1), rep(i-1, length(i)))); 
  combs<-do.call('rbind', combs);
  
  c<-b;
  for (i in 1:nrow(combs)) {
    c[!is.na(a)&!is.na(b)&a==combs[i,1]&b==combs[i,2]]<-i-1;
  }; 
  
  c;
}