# Merge regions of genotype matrixes, all regions must have the same samples in the same order
mergeRegions<-function(fn=NA, gt=NA) {
  # fn  if filenames specified, load matrix from files first before merging
  # gt  genotype matrixes with annotation. when filenames not given, merge the matrix directly from input
  
  if (!identical(fn, NA)) gt<-lapply(fn, function(fn) eval(parse(text=load(fn))));
  
  # merging genotype matrix
  mtrx<-lapply(gt, function(gt) gt$genotype);
  mtrx<-do.call('rbind', mtrx);
  
  # merging annotation
  anno<-lapply(gt, function(gt) gt$annotation);
  anno<-lapply(anno, function(anno) {
    alt<-anno$ALT;
    if (class(alt)=='DNAStringSetList') anno else {
      n<-elementLengths(alt);
      a<-unlist(alt);
      vr<-gsub('[^ACGTN]', '', a);
      a[a!=vr]<-'N';
      a<-DNAStringSet(a);
      anno$ALT<-split(a, rep(1:length(n), n));
      anno;
    }
    
  })
  if (length(anno)>1) anno<-do.call('c', anno) else anno<-anno[[1]];
  
  
  # merging alleles
  alll<-lapply(gt, function(gt) gt$allele);
  nc<-max(sapply(alll, ncol)); # maximum number of columns
  alll<-lapply(alll, function(alll) if (ncol(alll)==nc) alll else {cbind(alll, matrix(rep('', (nc-ncol(alll))*nrow(alll)), nr=nrow(alll)))});
  alll<-do.call('rbind', alll);
  
  list(genotype=mtrx, annotation=anno, allele=alll);
}