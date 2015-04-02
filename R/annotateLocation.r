# annotate variants with information related to known genes
annotateLocation<-function(rows, chromosome.mapping, species='human') {
  # rows      A GRanges object created from the rowData field of a VCF file, must include a metadata column called 'ALT' with sequences of alternatives alleles
  # chromosome.mapping    Chromosome name mapping
  # species   Currently accept only human
  
  library(VariantAnnotation);
  library(AnnotationDbi);
  
  if (species=='human') {
    library(BSgenome.Hsapiens.UCSC.hg19);
    ref<-Hsapiens
    
    # known transcripts
    library(TxDb.Hsapiens.UCSC.hg19.knownGene);
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
    
    # known gene IDs and names
    library(org.Hs.eg.db); 
    sym<-org.Hs.egSYMBOL; 
    keys<-mappedkeys(sym); 
    gn<-unlist(as.list(sym[keys])); 
    
    # make chromosome names consistent
    chr<-seqlevels(rows);
    mp<-unlist(sapply(chromosome.mapping, function(c) c[c %in% txdb$.chrom][1]))[chr];
    mp[is.na(mp)]<-chr[is.na(mp)];
    seqlevels(rows)<-as.vector(mp);
    seqlengths(rows)<-NA;
  }
  
  loc<-locateVariants(rows, txdb, AllVariants());
  
  meta<-as.data.frame(elementMetadata(loc));
  meta<-cbind(varID=names(rows)[meta[,2]], GENE=gn[as.character(meta$GENEID)], meta[,-2]);
  
  gn2v<-split(as.vector(meta[[1]]), as.vector(meta[[2]]));
  gn2v<-lapply(gn2v, unique);
  gn2v<-lapply(gn2v, sort);
  
  v2gn<-split(as.vector(meta[[2]]), as.vector(meta[[1]]));
  v2gn<-lapply(v2gn, function(x) sort(unique(x[!is.na(x)])));
  
  list(annotated=meta, gene2variant=gn2v, variant2gene=v2gn, genes=gn);
}
