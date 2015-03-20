# Use the outputs of the <parseVcf> function as inputs to annotate variants with genes and intragenic regions
# map variants to genes and intragenic regions, as well as mutation type
map2gene<-function(rows, chromosome.mapping, species='human', targeted=NA, cadd.cnm=NA) {
  # rows      A GRanges object created from the rowData field of a VCF file, must include a metadata column called 'ALT' with sequences of alternatives alleles
  # chromosome.mapping    Chromosome name mapping
  # species   Currently accept only human
  # targeted  GRanges object of targeted regions on Exome sequencing  
  # cadd.cnm  ElementMetadata column name of cadd scores
  
  # annotated with genes and intragenic region
  anno1<-annotateLocation(rows, chromosome.mapping, species);
  
  # annotated with mutation type
  anno2<-annotateCoding(rows, chromosome.mapping, species);
  
  if (species=='human') {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene);
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
    
    # make chromosome names consistent
    chr<-seqlevels(rows);
    mp<-unlist(sapply(chromosome.mapping, function(c) c[c %in% txdb$.chrom][1]))[chr];
    mp[is.na(mp)]<-chr[is.na(mp)];
    seqlevels(rows)<-as.vector(mp);
    seqlengths(rows)<-NA;
  }
  
  # whether variants are within given targeted regions
  if (class(targeted)=='GRanges') ct<-countOverlaps(rows, targeted) else ct<-NA;
  ct[ct>1]<-1;
  
  # CADD scores of variants when available
  if (class(cadd.cnm)=='character' & cadd.cnm %in% colnames(elementMetadata(rows))) cadd<-elementMetadata(rows)[, cadd.cnm] else cadd<-NA;
  
  # column names of different types of regions
  nm1<-c('intergenic', 'promoter', 'fiveUTR', 'coding', 'spliceSite', 'intron', 'threeUTR');
  rg<-matrix(nr=length(rows), nc=length(nm1), data=0, dimnames=list(1:length(rows), nm1));
  id<-lapply(nm1, function(x) as.vector(anno1[[1]][[1]])[anno1[[1]]['LOCATION']==x]);
  for (i in 1:length(nm1)) rg[names(rows) %in% id[[i]], i]<-1;
  
  # column names of different types of mutations
  nm2<-c("synonymous", "nonsynonymous", "nonsense", "frameshift");
  mt<-matrix(nr=length(rows), nc=length(nm2), data=0, dimnames=list(1:length(rows), nm2));
  id<-lapply(nm2, function(x) as.vector(anno2[[1]][[1]])[anno2[[1]]['CONSEQUENCE']==x]);
  for (i in 1:length(nm2)) mt[names(rows) %in% id[[i]], i]<-1;
  
  anno<-data.frame(varID=names(rows), rg, mt);
  
  anno$gene=anno1$variant2gene[names(rows)];
  
  if (!identical(ct, NA)) anno<-cbind(anno, targeted=ct); 
  if (!identical(cadd, NA)) anno<-cbind(anno, cadd=cadd);
  
  # return values
  list(
    summary=anno,   
    location=anno1$annotated,
    mutation=anno2$annotated,
    gene2variant=anno1$gene2variant,
    variant2gene=anno1$variant2gene,
    alternative.code=anno2$coded
  )
}