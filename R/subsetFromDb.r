# Query a table in the genotype database to get a subset of variants and/or samples
subsetFromDb<-function(db.name, tbl.name, vs=c(), ss=c(), as.data.frame=TRUE) {
  # db.name     Full path to the SQLite db
  # tbl.name    Table name
  # vs          IDs of a subset of variants to retrieve. Retrieve all variants if not specified
  # ss          IDs of a subset of samples to retrieve. Retrieve all samples and annotation columns if not specified
  
  library(dplyr);
  
  db<-src_sqlite(db.name);
  tbl<-tbl(db, 'dep20_filtered');
  
  # get variant subset
  if (length(vs) > 0) {
    tbl<-dplyr::filter(tbl, id %in% vs);
  }
  
  # get sample subset
  if (length(ss) > 0) {
    ind<-which(colnames(tbl) %in% ss);
    tbl<-dplyr::select(tbl, ind);
  }
  
  tbl<-collect(tbl);
  
  if(as.data.frame) tbl<-as.data.frame(tbl);
  
  tbl;
}