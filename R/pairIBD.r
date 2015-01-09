# IBD analysis of samples
runIBD<-function(g, gr, fn='pairIBD') {
  # g     Genotype matrix
  # gr    GRanges object of the variant locations, must match rows of the genotype matrix
  # fn    Name of gds file
  library(SNPRelate);
  
  g[g<0 | g>2]<-NA; # only use values 0, 1, 2 (AA, AB, BB)
  
  # remove null values
  #flg<-!is.na(g[,1]) & !is.na(g[,2]);
  #flg[g[,1]==g[,2]]<-FALSE;
  #g<-g[flg,];
  #gr<-gr[flg];
  
  ch<-as.vector(seqnames(gr));
  flg<-ch!='X' & ch!='Y';
  g<-g[flg,];
  gr<-gr[flg];
  
  ch<-as.integer(as.vector(seqnames(gr)));
  lc<-start(gr);
  
  # remove low density regions
  #flg<-selectIBD(gr, step=500000, min.variants=100);
  #g<-g[flg, ];
  #gr<-gr[flg];

  fn<-paste(fn, '.gds', sep='');
  
  # create gds file
  snpgdsCreateGeno(fn, genmat=g, sample.id=colnames(g), snp.id=1:length(gr), snp.chromosome=ch, snp.position=lc);
  
  # connet to gds file
  fl<-openfn.gds(fn); 
 
  snpset<-snpgdsLDpruning(fl, ld.threshold=0.1);
  snps<-as.vector(unlist(snpset));
  
  # run IBD using MLE method
  ibd<-snpgdsIBDMoM(fl, snp.id=snps, remove.monosnp=FALSE, kinship=TRUE, autosome.only=TRUE);
  
  # get kinship
  coeff<-snpgdsIBDSelection(ibd)
  
  closefn.gds(fl);
  
  list(IBD=ibd, Coefficients=coeff, Genotype=g);
}