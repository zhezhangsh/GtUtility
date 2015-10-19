# Give a "CollapsedVCF" object, retrieve it's AD field to get sequencing depth of both alleles
obtainAlleleDepth<-function(vcf) {
  ad<-geno(vcf)$AD;
  
  cnm<-colnames(ad); 
  ct<-lapply(cnm, function(nm) {
    cat(nm, '\n');
    d<-ad[, nm];
    n<-elementLengths(d);
    d<-unlist(d, use.names=FALSE);
    sm<-cumsum(n);
    ind<-c(1, sm+1)[1:length(n)];
    cbind(REF=d[ind], ALT=d[ind+1]);
  });
  
  c<-lapply(1:2, function(i) sapply(ct, function(ct) ct[, i]));
  colnames(c[[1]])<-colnames(c[[2]])<-cnm;
  rownames(c[[1]])<-rownames(c[[2]])<-rownames(ad);
  names(c)<-c('REF', 'ALT');
  c;
}