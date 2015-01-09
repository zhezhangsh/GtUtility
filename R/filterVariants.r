# Filter variants by call rate, ID, etc.
filterVariants<-function(gt, call.rate=1/ncol(gt[[1]]), removeMono=TRUE, haveID=FALSE, targeted=NA, targeted.xpd=0, qcScore=c(-Inf, Inf), passGATK=FALSE) {
  # gt            genotype matrixes with annotation, a list with 3 elements: the matrix, annotation as a GRanges object, and an allele matrix
  # call.rate     minimal call rate of a variant, default is at least one valid call
  # removeMono    whether to remove monophonic variants
  # haveID        remove variants without dbSNP IDs (rs) if TRUE
  # targeted      remove variants not within targeted region if a GRanges object
  # targeted.xpd  number of bp to expand from ends of targeted regions
  # passGATK      remove variants not pass GATK filtering if TRUE
  # qcScore       remove variants with QUAL score out of given range
  
  flg<-rep(TRUE, length(gt$annotation));
  
  g<-gt[[1]];
  
  # count occurance of each allele
  mx<-max(g, na.rm=TRUE);
  ct<-matrix(nr=nrow(g), nc=mx+1);
  ct[is.na(ct)]<-0;
  for (i in 0:mx) {
    ind<-which(!is.na(g) & g==i);
    ri<-ind %% nrow(g);
    n<-table(ri); 
    rows<-as.numeric(names(n));
    rows[rows==0]<-nrow(g);
    ct[rows, i+1]<-ct[rows, i+1]+as.integer(n);
  }
  
  # remove variants with call rates lower than cutoff
  flg[rowSums(ct)<call.rate]<-FALSE;
  
  # remove variants having only one genotype
  if (removeMono) {
    ct0<-ct;
    ct0[ct0>1]<-1;
    flg[rowSums(ct0)<2]<-FALSE;
  }
  
  # remove variants without dbSNP ID
  if (haveID) {
    ids<-rownames(g);
    ind<-grep('^rs', ids, ignore.case=TRUE);
    if (length(ind) != nrow(g)) flg[-ind]<-FALSE;
  }
  
  # select variants within targeted regions, allow for expanding at the ends
  if (class(targeted) == 'GRanges') {
    targeted.xpd<-max(0, targeted.xpd);
    c<-countOverlaps(gt[[2]], targeted, maxgap=targeted.xpd);
    flg[c==0]<-FALSE;
  }
  
  # select variants passed GATK filtering
  if (passGATK) {
    flt<-as.vector(gt[[2]]$FILTER);
    if (length(flt) == nrow(g)) flg[flt !='PASS']<-FALSE;
  }
  
  # select variants with QUAL score within given range
  if (min(qcScore) > -Inf | max(qcScore) < Inf) {
    mn<-min(qcScore);
    mx<-max(qcScore);
    sc<-as.vector(gt[[2]]$QUAL);
    if (length(sc) == nrow(g)) flg[is.na(sc) | sc<mn | sc>mx]<-FALSE;
  }
  
  gt[[1]]<-gt[[1]][flg, , drop=FALSE];
  gt[[2]]<-gt[[2]][flg];
  gt[[3]]<-gt[[3]][flg, , drop=FALSE];
  
  gt;
}