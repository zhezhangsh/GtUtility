# Add annotation to a set of variants
mapAnnotation<-function(gr1, gr0, ref.alt=1:2) {
  # gr1     The location of a set of variants
  # gr0     The existing annotation of variants
  # ref.alt Column index of REF and ALT alleles
  library(GenomicRanges);
  
  gr0<-gr0[countOverlaps(gr0, gr1, type='equal', ignore.strand=TRUE)>0];
  
  ol<-as.matrix(findOverlaps(gr1, gr0, type='equal', ignore.strand=TRUE)); 
  
  if (nrow(ol)==0) NA else {
    ol<-ol[!duplicated(ol[, 1]), , drop=FALSE]; 
    
    g1<-gr1[ol[,1]];
    g0<-gr0[ol[,2]];
    
    meta1<-elementMetadata(g1);
    meta0<-elementMetadata(g0);
    
    ref.alt<-ref.alt[1:2];
    ref.alt<-ref.alt[!is.na(ref.alt) & ref.alt>0];
    if (length(ref.alt)==2 & max(ref.alt)<=ncol(meta1) & max(ref.alt)<=ncol(meta0)) {
      ref1<-meta1[, ref.alt[1]];
      alt1<-meta1[, ref.alt[2]];
      ref0<-meta0[, ref.alt[1]];
      alt0<-meta0[, ref.alt[2]]; 
      
      lll<-cbind(ref1, alt1, ref0, alt0); 
      for (i in 1:ncol(lll)) {
        lll[, i]<-toupper(lll[, i]);
        lll[, i]<-gsub('[^ACGT]', '', lll[, i]);
      }
      
      flg<-lll[, 1]==lll[, 3] & lll[,2]==lll[,4]; 
      g1<-g1[flg];
      g0<-g0[flg];
    }
    
    anno<-as.data.frame(elementMetadata(g0)); 
    rownames(anno)<-names(g1); 
    
    ex<-names(gr1)[!(names(gr1) %in% rownames(anno))];
    ex<-matrix(NA, nr=length(ex), nc=ncol(anno), dimnames=list(ex, colnames(anno))); 
    anno<-rbind(anno, ex); 
    
    anno[names(gr1), ];
  }
}