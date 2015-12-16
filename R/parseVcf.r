# parse information stored in an R VCF object to get more accessible values for downstream analysis, especially the genotype calls (replace 0/0, 0/1, 1/1, ... with 0, 1, 2, ...) 
parseVcf<-function(vcf, variants=NA, samples=NA, chromosome=NA, range=NA, snp.only=FALSE, info.fields=names(info(vcf)), filter=list(), keep.no.geno.qual=FALSE,
                   score=c(20, Inf), depth=c(10, Inf), depth.alleles=c(3, Inf), hetero.alleles.ratio=c(0.1, 10), other.alleles.pct=c(0, 0.05), missing.default=NA,
                   phased=FALSE, outputs=1) {
  # vcf:   A VCF or CollapsedVCF object
  # chromosome, range, variants, samples, snp.only:   Parameters used to filter variants; failed variants will be completely removed.
      # chromosome:   limit the variants to one or more chromosomes; no filitering if NA
      # range:        an GRanges object defining the genomic range within which the selected variants will be in; no filtering if NA
      # variants:     the IDs of one or more variants; no filtering if NA      # variants:     the IDs of one or more variants; no filtering if NA
      # samples:      the IDs of one or more samples; no filtering if NA  
      # snp.only:     select SNPs only (both REF and ALT are single base); no filtering if NA
      # info.field:   fields in the INFO table to be included in the outputs
  # filter:   Extra filtering criteria based on the elementMetadata fields, for example:
      # 'QUAL' = c(100, Inf) means the QUAL field in metadata must have values 100 or higher; and
      # 'FILTER' = c('PASS', 'pass') means the field called 'FILTER' must equal to 'PASS' or 'pass'
  # score, depth, depth.alleles, depth.alleles.ratio:     Cutoff values used to decide whether a single genotype call will be selected
      # score:                A selected genotype call should have GQ score within the given range, such as c(25, 99)
      # depth:                A selected genetype call should have total sequencing depth of alleles contributing to the genotype within the given range, such as (10, 200)
      # depth.alleles:        A selected heterozygous call should have sequencing depth of both alleles within the given range, such as (3, 200)
      # hetero.alleles.ratio: A selected heterozygous call should have the ratio of sequencing depths of two alleles within the given range, such as (1/10, 10)
      # other.alleles.pct:    A selected heterozygous call should have total reads from alleles not contributing to the genotype within the given range, such as (0, 0.1)
      # missing.default:      The default value to replace calls failed filtering 
  # phased:   If TRUE the genotype call is phased; 0/1 will be coded as 1 and 1/0 will be coded as -1; <missing.default> will be set to NA
  # outputs:  What to output; Onlfy the genotype matrix if 0; Variant annotation and alleles if 1.
  # keep.no.geno.qual: whether to keep the genotype call if it has no information associated to the call, such as depth, quality score, ...
  ###############################################################################################################
  #Filter variants by location, ID, etc.
  # selected variants and samples by ID
  variants<-variants[variants %in% rownames(vcf)];
  samples<-samples[samples %in% colnames(vcf)];
  if (length(variants)>0) vcf<-vcf[variants,];
  if (length(samples)>0) vcf<-vcf[, samples];

  rows<-rowRanges(vcf);
  
  # merge specified fields in INFO into row annotation
  info<-info(vcf);
  info.fields<-info.fields[info.fields %in% colnames(info)];
  if (length(info.fields)>0) {
    for (i in 1:length(info.fields)) elementMetadata(rows)[[info.fields[i]]]<-info[[info.fields[i]]];
  }

  # select variants on specified chromosome(s)
  chromosome<-chromosome[chromosome %in% seqlevels(rows)]; # chromosome names that have been included in the VCF
  if (length(chromosome)>0) rows<-rows[seqnames(rows) %in% chromosome];
  
  # select variants within given range
  if (class(range)=='GRanges')  rows<-rrows[countOverlaps(rows, range)>0];

  # select variants based on metadata fields
  names(filter)<-toupper(names(filter));
  names(elementMetadata(rows))<-toupper(names(elementMetadata(rows)));
  filter<-filter[names(filter) %in% names(elementMetadata(rows))];
  if (length(filter)>0) {
    for (i in 1:length(filter)) {
      nm<-names(filter)[i];
      vl<-filter[[i]];
      mt<-elementMetadata(rows)[[nm]];
      
      if (is.numeric(mt)) {
        if (length(vl)==1) rows<-rows[mt>=vl] else if (length(vl)>1) rows<-rows[mt>=min(vl[1:2]) & mt<=max(vl[1:2])];
      } else {
        rows<-rows[mt %in% vl];
      }
    }    
  }

  # select variants with just one base replacement
  if (snp.only) {
    rows<-rows[nchar(elementMetadata(rows)$REF)==1];
    alt<-elementMetadata(rows)$ALT; # alternative alleles
    mx<-sapply(nchar(alt), max); # max length of alternative alleles    
    rows<-rows[mx==1];
  }
  
  # filter out unwanted variants
  vcf<-vcf[rownames(vcf) %in% names(rows), ];
  
  if (nrow(vcf)==0) { # no more variants left
    gt<-matrix(nr=0, nc=ncol(vcf));
    colnames(gt)<-colnames(vcf);
    gt; # return an empty matrix
  } else { # some varaints left
    ref<-elementMetadata(rows)$REF; # reference allele
    alt<-elementMetadata(rows)$ALT; # alternative alleles
    nl<-elementLengths(alt); # number of altermative alleles
    alt<-split(as.character(unlist(alt)), rep(1:length(nl), nl)); # convert to characters
    
    # convert alleles (reference and alternatives) into a matrix, number of columns equal to the largest number of alleles of all variants
    lll<-matrix(nr=length(rows), nc=max(nl));
    for (i in 1:ncol(lll)) lll[nl>=i, i]<-sapply(alt[nl>=i], function(x) x[i]);
    lll<-cbind(as.character(ref), lll);
    rownames(lll)<-names(rows);
    lll[is.na(lll)]<-'';
    
    # codes to convert genotype calls to integer
    combs<-lapply(1:ncol(lll), function(i) cbind(0:(i-1), rep(i-1, length(i)))); 
    combs<-do.call('rbind', combs);
    code<-(1:nrow(combs))-1;
    names(code)<-paste(combs[,1], combs[,2], sep='/');
    rownames(combs)<-names(code);
    code0<-code; 
    names(code0)<-paste(combs[,2], combs[,1], sep='/');
    if (phased) code0<--1*code0;
    code<-c(code, code0);
    code1<-code;
    names(code1)<-sub('/', '|', names(code1));
    code<-c(code, code1);
    code<-code[!duplicated(names(code))];
    combs<-do.call('rbind', strsplit(names(code), '[/|]')); 
    combs<-cbind(as.numeric(combs[,1]), as.numeric(combs[,2]));
    rownames(combs)<-names(code);
    cmb<-combs[combs[,1]==0 & grepl('/', rownames(combs)), ];
    rownames(cmb)<-0:(ncol(lll)-1);
    combs<-rbind(combs, cmb);
    c1<-code[paste('0', rownames(cmb), sep='/')];
    names(c1)<-rownames(cmb);
    code<-c(code, c1);
    
    # Retrieve fields
    geno<-geno(vcf); 
    gq<-as.vector(geno$GQ); # genotype call quality
    gt<-as.vector(geno$GT); # genotype calls
    ad<-as.list(geno$AD); # allele depth
    
    ##################################################################
    # filtering based on genotyping call quality and depth
    if (length(gt) == length(gq) & length(ad)) {
      # create a matrix of read depth per allele
      tps<-rownames(combs);
      tps<-tps[tps %in% gt];
      cts<-cbind(A1=rep(NA, length(gt)), A2=rep(NA, length(gt)));
      for (i in 1:length(tps)) {
        ind<-combs[tps[i],]+1;
        d<-ad[gt==tps[i]]; # allele depth of calls with given genotype
        if (ind[1]==ind[2]) { # homozygous
          ct<-cbind(sapply(d, function(d) d[ind[1]]), rep(0, length(d)));
          cts[gt==tps[i], ]<-ct;
        } else cts[gt==tps[i], ]<-t(sapply(d, function(d) d[ind])); # heterozygous calls
      }
      sm<-sapply(ad, sum);
      cts<-cbind(cts, Others=sm-rowSums(cts)); # read counts by alleles
      
      # preparing filtering
      tps.homo<-rownames(combs)[combs[,1]==combs[,2]]; # homozygous genotypes
      is.homo<-gt %in% tps.homo; # whether is homozygous call
      ttl<-cts[,1]+cts[,2]; # total reads contributing to the call
      rt<-cts[,2]/cts[,1]; # ratio of reads of two allels, should be 0 if the call is homozygous
      oth<-cts[,3]/sm; # percentage of reads of alleles not contributing to the call
      
      f1<-(!is.na(gq) & gq>=min(score) & gq<=max(score)) | (is.na(gq) & keep.no.geno.qual); # genotype quality score within range?
      f2<-(ttl>=min(depth) & ttl<=max(depth)) | (is.na(ttl) & keep.no.geno.qual); # total read counts within range?
      f3<-(is.homo | (cts[,1]>=min(depth.alleles) & cts[,2]>=min(depth.alleles)&cts[,1]<=max(depth.alleles)&cts[,2]<=max(depth.alleles))) | (is.na(cts[,1]) & is.na(cts[,2]) & keep.no.geno.qual); # both allele reads within range? (heterozygous only)
      f4<-(is.homo | (!is.na(rt) & rt>=min(hetero.alleles.ratio) & rt<=max(hetero.alleles.ratio))) | (is.na(rt) & keep.no.geno.qual); # ratio of reads from two alleles within range? (heteronzygous only)
      f5<-(!is.na(oth) & oth>=min(other.alleles.pct) & oth<=max(other.alleles.pct)) | (is.na(oth) & keep.no.geno.qual); # percentage of reads from other alleles within range?
      flg<-f1 & f2 & f3 & f4 & f5; # keep the genotype call when flag=TRUE
    } else flg<-c();
    
    # generate output matrix
    nocode<-setdiff(unique(gt), names(code));
    if (length(nocode) > 0) {c1<-rep(NA, length(nocode)); names(c1)<-nocode; code<-c(code, c1)}
    g<-code[gt]; # convert genotype characters to integers
    if (length(flg) == length(g)) g[!is.na(g) & !flg]<-missing.default;
    g<-matrix(as.integer(g), nr=nrow(vcf));
    rownames(g)<-rownames(vcf);
    colnames(g)<-colnames(vcf);
        
    if (outputs==0) g else {
      list(genotype=g, annotation=rows, alleles=lll);
    }
  }
}
