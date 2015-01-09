# select variant subsets suiting IBD analysis as sparsed variants reduced IBD accuracy
selectIBD<-function(v, ws=10^6, step=10^6, min.variants=250) {
  # v     GRanges, location of all variants
  # ws    Window size of the region for variant density
  # step  Step size to move the windows along chromosome
  # min.variants  Minimal number of variants within each window. Only window includes the same number or more variants will be used.
  
  chr<-as.vector(seqnames(v));
  pos<-end(v);
  
  # chromosome length
  len<-sapply(split(pos, chr), max);
  
  # Windows of all chromosomes
  gr<-lapply(names(len), function(c) {
    l<-as.numeric(len[[c]]);
    stt<-seq(1, l, step);
    GRanges(c, IRanges(stt, width=ws));
  });
  
  gr<-do.call('c', gr);
  
  ct<-countOverlaps(gr, v);
  
  gr<-gr[ct>=min.variants];
  
  countOverlaps(v, gr)>0;
}