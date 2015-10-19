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

# Parse the VCF file directly to retrieve snpEFF annotation
parseEFF<-function(fn.vcf, fn.out, n.total=-1, n.batch=1000000, n.cluster=8) {
  # fn.vcf    Input VCF file
  # fn.out    Name of output file
  # n.total   Total number of lines to read in; all lines if -1
  # n.batch   Number of lines to read in each batch
  # n.cluster Number of clusters for parallel computing
  
  if (n.total<=0) n.total<-Inf;
 
  cur.line<-1;

  write.table(matrix(c('LOW', 'MODERATE', 'HIGH', 'MODIFIER'), nr=1), fn.out, row=FALSE, col=FALSE, qu=FALSE)
  
  ########################################################################################
  # Parse a number of lines
  parser<-function(ls) {
    ########################################################################################
    # ID Parser
    parserId<-function(ls) {
      if (length(ls) == 0) c() else {      
        fld<-t(sapply(strsplit(ls, '\t'), function(x) x[1:5]));
        id<-fld[, 3];
        id[is.na(id)]<-'';
        id[id=='.']<-paste(fld[id=='.', 1], ':', fld[id=='.', 2], '_', fld[id=='.', 4], '/', sapply(strsplit(fld[id=='.', 5], ','), function(x) x[1]), sep='');
        id;
      }
    }
    ########################################################################################
    
    codes=c('LOW', 'MODERATE', 'HIGH', 'MODIFIER');
    #ls<-ls[grep(';EFF=', ls)];
    if (length(ls) == 0) matrix(nr=0, nc=4) else {
      snp.id<-parserId(ls); print(head(snp.id));
      ls<-sapply(strsplit(ls, ';EFF='), function(l) l[2]);
      ls<-sapply(strsplit(ls, '[;\t]'), function(l) l[1]);
      ls<-strsplit(ls, ',');
      eff<-unlist(ls, use.names=FALSE);
      cd<-sapply(strsplit(eff, '[(|]'), function(x) x[2]);
      out<-matrix(0, nr=length(ls), nc=length(codes), dimnames=list(1:length(ls), codes));
      
      ids<-rep(1:length(ls), sapply(ls, length));
      for (i in 1:length(codes)) out[unique(ids[cd==codes[i]]), i]<-1;
      rownames(out)<-snp.id;
      out;
    }
  }
  
  ########################################################################################
  # Parse a number of lines in parallel
  parseLines<-function(ls, n.cluster) {
    n.cluster<-max(1, n.cluster);
    ls<-split(ls,  rep(1:n.cluster, each=ceiling(length(ls)/n.cluster))[1:length(ls)]);
    
    library(snow);
    cl<-snow::makeCluster(n.cluster, type='SOCK');
    out<-clusterApplyLB(cl, ls, parser);
    stopCluster(cl);
    do.call('rbind', out);
  }
  
  ########################################################################################  
  con<-file(fn.vcf, open='r');
  while (length(ls<-readLines(con, n=n.batch, warn = FALSE)) > 0 & cur.line<n.total) {
    cat(cur.line, '\n');
    eff<-parseLines(ls, n.cluster);
    write.table(eff, fn.out, append=TRUE, quote=FALSE, sep='\t', row.names=TRUE, col.names=FALSE);
    cur.line <- cur.line + n.batch;
  } 
  
  close(con);
}