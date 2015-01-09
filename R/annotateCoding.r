# annotate variants with information about coding changes
# when there are multiple alternatives, perform annotation separately on each of them
annotateCoding<-function(rows, chromosome.mapping, species='human',  bads=c('frameshift', 'nonsense', 'nonsynonymous')) {
  # rows      A GRanges object created from the rowData field of a VCF file, must include a metadata column called 'ALT' with sequences of alternatives alleles
  # chromosome.mapping    Chromosome name mapping
  # species   Currently accept only human
  # bads      Types of 'bad' alternative alleles, usually those causing protein coding changes
  
  library(VariantAnnotation);
 
  if (species=='human') {
    library(BSgenome.Hsapiens.UCSC.hg19);
    ref<-Hsapiens
    
    library(TxDb.Hsapiens.UCSC.hg19.knownGene);
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;

    # make chromosome names consistent
    chr<-seqlevels(rows);
    mp<-unlist(sapply(chromosome.mapping, function(c) c[c %in% txdb$.chrom][1]))[chr];
    mp[is.na(mp)]<-chr[is.na(mp)];
    seqlevels(rows)<-as.vector(mp);
    seqlengths(rows)<-NA;
  }
  
  alt<-rows$ALT;
  len<-elementLengths(alt); # number of alternative alleles
  
  gr<-rep(rows, len);
  var<-unlist(alt);
  if (class(var) != "DNAStringSet") {
    vr<-gsub('[^ACGTN]', '', var); # remove illegal base character
    var[var!=vr]<-''; 
    var<-DNAStringSet(var);
  }
  # index of each alternative allele
  gr$ind0<-unlist(rep(1:length(rows), len)); # index of the variant
  gr$ind1<-unlist(lapply(len, function(len) 1:len)); # index within the alternative alleles of the same variant
  
  ################ predicting coding changes ################
  prd<-predictCoding(gr, txdb, seqSource=ref, varAllele=var);
  ###########################################################
  
  # full annotation
  meta0<-elementMetadata(gr);
  meta1<-elementMetadata(prd);
  meta1<-as.data.frame(meta1[, !(colnames(meta1) %in% colnames(meta0))]);
  meta1<-cbind(varID=names(rows)[prd$ind0], varIndex=prd$ind1, meta1);
  
  # re-code genotype based on mutation type; 0 if NA, 1 if not a 'bad' mutation, 2 if a 'bad' mutation such as missense
  coded<-matrix(nr=length(len), nc=max(len));
  coded[is.na(coded)]<-0;
  for (i in 1:max(len)) coded[len>=i, i]<-1;
  chg<-as.vector(prd$CONSEQUENCE);
  ind<-sapply(c('ind0', 'ind1'), function(x) elementMetadata(prd)[, x][chg %in% bads]);
  for (i in 1:max(ind[,2])) coded[unique(ind[ind[,2]==i,1, drop=FALSE]), i]<-2;
  rownames(coded)<-1:nrow(coded);

  list(annotated=meta1, coded=coded);
}