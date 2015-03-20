# Perform trio-based analysis, such as de novo variants, autosomal recessive analysis, and TDT
analyzeTrio<-function(geno, trio, analyses=list(
    de.novo = TRUE, # whether to summarize de novo mutations
    ar = TRUE # whether to summarize autosomal recessive variants
  )) {
  # gene        Coded genotype matrix (-1, 0, 1, 2, ...)
  # trio        A data.frame of three ID columns: child, parent1, parent2; corresponding to column names in the <geno> matrix
  # analyses    Whether to run a specific analysis, such as summarize de novo mutations
  
  # all results
  results<-list();
  
  cn<-colnames(geno);
  find.sample<-sapply(1:3, function(i) as.vector(trio[, i]) %in% cn & !is.na(as.vector(trio[, i])));
  trio<-trio[rowSums(find.sample)==3, , drop=FALSE];
  if (nrow(trio) == 0) {
    warnings("No complete trio was found in the genotyping data");
    results;
  } else {
    g<-lapply(1:3, function(i) geno[, as.vector(trio[, i]), drop=FALSE]);
    
    ##################################################################################################
    # summarize de novo mutations
    if (exists(analyses$de.novo)) if (analyses$de.novo[1] == TRUE) {
      ttl<-!is.na(g[[1]]) & !is.na(g[[2]]) & !is.na(g[[3]]) & g[[2]]==0 & g[[3]]==0;
      yes<-ttl & g[[1]]==1;
      no<-ttl & g[[1]]==0;
      results$de.novo<-list(Yes=yes, No=no, Summary=list(By.Variant=cbind(Yes=rowSums(yes), No=rowSums(no)), By.Sample=cbind(Yes=colSums(yes), No=colSums(no))));
    }
    
    results;
  }
}