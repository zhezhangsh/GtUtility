work.home<-"/nas/is1/zhangz/hts/projects/pcgc/2014-09_gVCF/VCF2R";
function.home<-"/nas/is1/zhangz/hts/projects/pcgc/code/functions"; # home of custom functions

path.out<-"/nas/is1/zhangz/hts/projects/pcgc/2014-09_gVCF/VCF2R";

info.fields<-c('CADD_SCALED', 'SIFT', 'AMR_AF', 'EUR_AF', 'AFR_AF', 'ASN_AF', '1kg_AF', 'CHOP_AF', 'dbNSFP_ESP6500_AA_AF', 'dbNSFP_ESP6500_EA_AF'); # fields to be included in the INFO of VCF fileslib
num.threads<-3; 

setwd(work.home);

fn<-readRDS('vcf_files.rds'); # named character vector, full path to VCF files. Use names as output file names
#fn<-c('test_235'="/nas/is1/PCGC/Ariellas_work/gVCF/test/annotate/wide-20140831-6703cc.235.targets.a.f.snpEff.1kg.evs.msigdb.dbnsfp.chop.cadd.vcf");

doOneSample<-function(fnm, fn, function.home, path.out, info.fields) {
#for (i in 1:length(fn)) {
# fnm<-names(fn)[i];
  print(fnm); 
  
  if (is.null(fnm) | fnm=='') fnm<-paste('Vcf', i, sep='');
  
  # load custom functions
  f<-dir(function.home);
  f<-f[grep('.r$', f, ignore.case=TRUE)];
  f<-paste(function.home, f, sep='/');
  sapply(f, source)->x;
  chr<-readRDS(paste(function.home, '/chromosomes.rds', sep='')); # GRanges object with information about chromosome length
  chr.mapping<-readRDS(paste(function.home, '/chromosome_synonym.rds', sep='')); # GRanges object with information about chromosome length
  
  pth<-paste(path.out, fnm, sep='/'); # folder name of output files
  print(pth);
  
  if (!file.exists(pth)) dir.create(pth);
  
  fs<-loadVcf(fn[fnm], paste(pth, fnm, sep='/'), 'hg19', chr, TRUE, TRUE); # load VCF, split chromosomes
  
  pth0<-paste(pth, '/genotype', sep='');
  if(!file.exists(pth0)) dir.create(pth0); # sub-folder to save genotype calls
  
  f.out<-sapply(fs, function(f) {
    vcf<-eval(parse(text=load(f)));
    f0<-sub(pth, pth0, f);
    gt<-parseVcf(vcf, info.fields=info.fields, depth=c(20, Inf), depth.alleles=c(5, Inf), hetero.alleles.ratio=c(0.25, 4)); 
    #gt<-parseVcf(vcf, info.fields=info.fields); 
    save(gt, file=f0);
    cd<-annotateCoding(gt[[2]], chr.mapping);
    gt[[1]]<-recodeCalls(gt[[1]], cd$coded); 
    save(gt, file=sub('.rdata', '_recoded.rdata', f0));
    f0
  });
  
 
  f.out;
}

library(snow);

if (num.threads>1) { 
  cl<-makeCluster(num.threads, type='SOCK');
  f.out<-clusterApplyLB(cl, names(fn), doOneSample, fn=fn, function.home=function.home, path.out=path.out, info.fields=info.fields);  
  try(stopCluster(cl));
} else {
  f.out<-lapply(names(fn), function(fnm) doOneSample(fnm, fn=fn, function.home=function.home, path.out=path.out, info.fields=info.fields)); 
}

####################################################################################
## merge genotype matrices
## (Optional) merging by two steps, chromosomes first and then all chromosomes
# gt<-mergeRegions(fn); 

## filter variants to remove those without genotype calls, or monopholic
# gt<-filterVariants(gt)

## check duplicated IDs

## annotate variants by mapping to genes
# anno<-map2gene(rows, chromosome.mapping, species='human', targeted=targets, cadd.cnm=cadd)
# gt$map2gene<-anno;

## add full annotation to variants
# fields<-"id\tchr\tpos\ttargeted\tgene\tnum_gene\tCADD\tGATK_QUAL\tGATK_FILTER\tNum_Allele\tREF\tALT\tALT2\tALT3\tALT4\tALT5\tALT6\tintergenic\tpromoter\tfiveUTR\tcoding\tspliceSite\tintron\tthreeUTR\tsynonymous\tnonsynonymous\tnonsense\tframeshift\tCount\tMAF\tHET\tcount_0_0\tcount_0_1\tcount_1_1\tcount_0_2\tcount_1_2\tcount_2_2\tcount_0_3\tcount_1_3\tcount_2_3\tcount_3_3\tcount_0_4\tcount_1_4\tcount_2_4\tcount_3_4\tcount_4_4\tcount_0_5\tcount_1_5\tcount_2_5\tcount_3_5\tcount_4_5\tcount_5_5\tcount_0_6\tcount_1_6\tcount_2_6\tcount_3_6\tcount_4_6\tcount_5_6\tcount_6_6"

## genotype counting function
#fct<-function(g) apply(g, 1, function(g) {x<-g[!is.na(g) & g!=-1]; sapply(0:27, function(i) length(x[x==i]))})