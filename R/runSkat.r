# run SKAT test in parallel
runSkatCl<-function(ids=list(), gt, samples0, samples1, cov=NA, resampling=0, method='optimal.adj', num.threads=8, trimming=TRUE, max.call=10^6, weights=NULL, verbose=FALSE) {
  # ids         A list of id sets, run SKAT test on each set of IDs. Test on the complete genotype matrix if empty
  # gt          Genotype call matrix
  # samples0    Sample IDs of control group
  # samples1    Sample IDs of case group
  # cov         A matrix whose rows matching the samples, each column is a covariant of test
  # resampling  Number of re-sampling to evaluate observed p values; no resampling if equals 0
  # trimming    Trimming variants to fit max.call if TRUE
  # max.call    The maximum number of calls that is allowed (n.variant * n.sample). Trim the variants if the number is over the limit. Default is one million
  # method      Method of SKAT test

  library(snow);
  
  num.threads<-max(num.threads, 1);
  
  # remove variants no in the matrix
  id<-unlist(ids, use.names=FALSE);
  nm<-rep(names(ids), sapply(ids, length));
  ids<-split(id[id %in% rownames(gt)], nm[id %in% rownames(gt)]);
  ids<-lapply(ids, unique);
  ids<-ids[sapply(ids, length)>0]; 
  
  # trim variants if the total number of variants is over the limits
  #if (max.call>0 & max.call<Inf) ids<-runSkatTrim(ids, gt, ceiling(max.call/ncol(gt)));
  
  # split ID sets 
  n<-sapply(ids, length);
  ids<-ids[order(n)][length(ids):1];
  n<-sort(n)[length(n):1];
  sm<-cumsum(n);
  sz<-ceiling(sm[length(n)]/num.threads);
  sz<-min(sz, 10^8/ncol(g)); 
  sz<-round(max(sz, max(n))); 
  ind<-sapply(1:ceiling(sum(n)/sz), function(i) min(which(sm >= min(sm[length(n)], i*sz))));
  ind<-c(0, unique(ind));
  id.sets<-lapply(2:length(ind), function(i) ids[(ind[i-1]+1):ind[i]]);

  if (resampling<0) resampling<-0 else resampling<-round(resampling); 
  
  cl<-makeCluster(num.threads, type='SOCK');
  stat<-clusterApplyLB(cl, id.sets, runSkat, gt=gt, samples0=samples0, samples1=samples1, cov=cov, resampling=resampling, trimming=trimming, max.call=max.call, method=method, weights=weights, verbose=verbose);  
  try(stopCluster(cl));

  if (resampling==0) do.call('rbind', stat) else {
    list(stat=do.call('rbind', lapply(stat, function(stat) stat[[1]])), resampling=do.call('cbind', lapply(stat, function(stat) stat[[2]])));
  }
}


# run SKAT test
runSkat<-function(ids=list(), gt, samples0, samples1, cov=NA, method='optimal.adj', resampling=0, trimming=TRUE, max.call=10^6, weights=NULL, verbose=FALSE) {
  # ids         A list of id sets, run SKAT test on each set of IDs. Test on the complete genotype matrix if empty
  # gt          Genotype call matrix
  # samples0    Sample IDs of control group
  # samples1    Sample IDs of case group
  # cov         A matrix whose rows matching the samples, each column is a covariant of test
  # resampling  Number of re-sampling to evaluate observed p values; no resampling if equals 0
  # trimming    Trimming variants to fit max.call if TRUE
  # max.call    The maximum number of calls that is allowed (n.variant * n.sample). Trim the variants if the number is over the limit. Default is one million
  # method      Method of SKAT test
  
  library(SKAT);
  
  
#########################################
  # trim variant set if its total size (n.variant * n.sample) is larger than given number
  runSkatTrim<-function(ids, gt, mx, min.callrate=0.85, method=c('maf', 'callrate', 'first', 'random')) {
    # ids         A list of id sets, run SKAT test on each set of IDs. Test on the complete genotype matrix if empty
    # gt          Genotype call matrix
    # max         The maximum number of variants that is allowed. Trim the variants if the number is over the limit
    
    if (!is.numeric(mx) | mx<2) mx<-2 else mx<-round(mx);
    
    if (!is.numeric(min.callrate)) min.callrate<-0;
    
    n<-sapply(ids, length); 
    
    if (max(n) <=mx) ids else {
      ids[n>mx]<-lapply(ids[n>mx], function(id) { 
        id<-unique(id);
        id<-id[id %in% rownames(gt)];
        
        if (length(id) <= mx) id else {
          g<-gt[id, , drop=FALSE];
          ct<-sapply(0:2, function(i) apply(g, 1, function(g) length(g[!is.na(g) & g>=0 & g<=2])));
          cr<-rowSums(ct)/ncol(gt);
          mf<-(ct[,2]+ct[,3]*2)/2/rowSums(ct);
          
          id<-id[cr>min.callrate & mf>0 & cr>0 & apply(ct, 1, max) < rowSums(ct)]; 
          if (length(id) <= mx) id else {
            if (method[1] == 'first') {
              id[1:mx];
            } else if (method[1] == 'random') {
              sample(id, mx);
            } else if (method[1] == 'maf') {
              id[order(mf[id])][1:mx]; # select ones with lowest MAF
            } else {
              id[order(cr[id])][length(id):1][1:mx];
            }
          }
        }
      })
      ids;
    } 
  }
#########################################
  
  if (length(ids)==0) ids<-list(all=rownames(gt));
  
  if (verbose) cat("removing variants no in the matrix ...\n");
  # remove variants no in the matrix
  id<-unlist(ids, use.names=FALSE);
  nm<-rep(names(ids), sapply(ids, length));
  ids<-split(id[id %in% rownames(gt)], nm[id %in% rownames(gt)]);
  ids<-lapply(ids, unique);
  ids<-ids[sapply(ids, length)>0]; 
  
  # trim variants if the total number of variants is over the limits
  
  if (max.call>0 & max.call<Inf & trimming) {
    if (verbose) cat("trimming variants if the total number of variants is over the limits ...\n");
    ids<-runSkatTrim(ids, gt, ceiling(max.call/ncol(gt)));
  }
  
  samples0<-samples0[samples0 %in% colnames(gt)];
  samples1<-samples1[samples1 %in% colnames(gt)];
  
  if (resampling<0) resampling<-0 else resampling<-round(resampling); 
  
  # create null model
  if (verbose) cat('creating null model ...\n');
  f<-rep(0:1, c(length(samples0), length(samples1)));
  if (class(cov) == 'matrix' | class(cov) == 'data.frame') {
    if (nrow(cov) == length(f)) {
    cov<-as.matrix(cov);
    obj<-SKAT_Null_Model(f~cov, out_type='D', n.Resampling=resampling); 
    } else obj<-SKAT_Null_Model(f~1, out_type='D', n.Resampling=resampling); 
  } else obj<-SKAT_Null_Model(f~1, out_type='D', n.Resampling=resampling); 
  
  gt<-gt[, c(samples0, samples1)];

  gt[gt<0 | gt>2]<-NA;

  skat<-lapply(ids, function(ids) {
    if (verbose) cat('Testing', length(ids), 'variants.\n');
    
    if (length(weights)>0) if(weights[1]==1) weights<-rep(1, length(ids)) else weights<-NULL;

    tryCatch({
      SKAT(t(gt[ids, , drop=FALSE]), obj, method=method, weights=weights);
      }, error=function(e) {list(p.value=1, param=list(n.marker=0, n.marker.test=0));}
    );
  });
  names(skat)<-names(ids);

  p<-sapply(skat, function(x) x$p.value[[1]]);
  N<-sapply(skat, function(x) x$param$n.marker);
  n<-sapply(skat, function(x) x$param$n.marker.test);
  
  stat<-cbind(Skat_P=p, Variant_Total=N, Variant_Tested=n);
  rownames(stat)<-names(ids);

  if (resampling > 0) {
    res<-sapply(skat, function(skat) skat$p.value.resampling);
    list(stat=stat, resampling=res);
  } else stat;
  
}
