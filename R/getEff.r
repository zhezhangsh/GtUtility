# Retrieve EFF field from a VCF object
getEFF<-function(vcf, codes=c('LOW', 'MODERATE', 'HIGH', 'MODIFIER')) {
  
  info<-info(vcf);
  eff<-info$EFF;
  eff.all<-unlist(eff, use.names=FALSE);
  cd<-sapply(strsplit(eff.all, '[(|]'), function(x) x[2]);
  
  id<-rownames(info);
  n<-sapply(eff, length);
  ids<-rep(id, n);
  
  out<-matrix(0, nr=length(id), nc=length(codes), dimnames=list(id, codes));
  
  for (i in 1:length(codes)) out[unique(ids[cd==codes[i]]), i]<-1;
  
  out;
}