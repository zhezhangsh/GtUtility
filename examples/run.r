# Run description: PCGC cases (total and subgroups) vs. Autism siblings, pathways include SHH/90Gene interactomes and their subsets

# Script description: this script performs one or multiple variant-, gene-, and pathway-level comparisons of sample subgroups

# The script works on Isolon within this R instance: /usr/local/bin/R

##############
### INPUT ####
##############

### Comment out the next two lines to re-install package from GitHub repo    
#library(devtools);
#install_github(host="github.research.chop.edu/api/v3",repo="BiG/GtUtility");
library(GtUtility);

############################################################################################
### OPTIONAL PARAMETERS ###
prms<-list( 
  num.threads=4, # use value between 1-12
  skat.weight=1, # all variants have the same weight if 1; SKAT default if 0
  ar.only=FALSE, # use autosomal recessive variants only
  sex.chromosomes=FALSE, # use variants on X/Y chromosomes
  limit.variant=TRUE, # limit variants to those mapped to tested genes and pathways
  limit.gene=TRUE # limit genes to those mapped to tested pathways
);
############################################################################################


###########################################################################################################
### Location of the output files 
path<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/outputs/mouse_pathways_SHHI_f2f_analysis/";
if(!file.exists(path)) dir.create(path);

############################################################################################################################################
### Location of genetype call matrix ###
# it's an integer matrix using 0,1,2 to represent AA, AB, BB, no calls should be NA
genotype.matrix<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/inputs/variants/genotype_all.rdata";

#############################################################################################################################################
### Location of variant annotation
# the rows of annotation data.frame should match exactly to the rows of genotyping matrix
# it includes pre-defined columns as specified in /nas/is1/PCGC/gVCF6703/README.pdf
variant.annotation<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/inputs/variants/annotation_all.rdata";

#######################################################################################################################################
### Gene to variant mapping
# a named list, each element includes all variants mapped to the same gene
gene2variant<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/inputs/genes/all_genes.rdata";

#######################################################################################################################################
### Pathway to gene mapping
# a named list, each element includes all genes mapped to the same pathway
pathway2gene<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/inputs/pathways/mouse_SHHI.rdata";

#####################################################################################################################################################
# The comparisons to perform
# A nested list, each element is a list of 2 named elements; with column/sample IDs of controls and cases
sample.groups<-"/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/inputs/comparisons/AUT-PCGC_Child_Subgroups.rdata"; 
# !THIS TEMPLATE PERFORMS PAIRWISE COMPARISONS BETWEEN DIFFERENT SAMPLE SUBGROUPS ON THE SAME VARIANTS, GENES, AND PATHWAYS #

##############################################################
### Whether to reload all data
# save time, but risky if FALSE
# rm(list=ls()) before each run to clean up workspace if FALSE
reload.all<-TRUE; 
##############################################################

#########################################################################################################################################################



    


# END OF USER-SPECIPIC PARAMETERS

    
    
    
    
    
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
########################################################################################################################################################################
    
# save log script
run.time<-format(Sys.time());
file.copy("/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/runs/run.r", 
          paste("/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/runs/run_log/", run.time, '.r', sep=''));
    
    ####################
    ### LOADING DATA ###
    ####################
    
if (!file.exists(path)) dir.create(path); 
    
cat('Loading data ...\n');
if (exists('geno')) {if (class(geno)!='matrix' | reload.all) geno<-eval(parse(text=load(genotype.matrix)))} else geno<-eval(parse(text=load(genotype.matrix))); 
if (!exists('anno') | reload.all) anno<-eval(parse(text=load(variant.annotation))); 
g2v<-eval(parse(text=load(gene2variant)));
p2g<-eval(parse(text=load(pathway2gene)));
grps<-eval(parse(text=load(sample.groups)));

# Automatically assign a name to an unnamed pathway
pth.names<-names(p2g);
if (length(pth.names)==0) pth.names<-paste('Pathway', 1:length(p2g), sep='_') else {
  ind<-which(is.na(pth.names) | pth.names=='');
  if (length(ind) > 0) pth.names[ind]<-paste('Pathway', ind, sep='_');
}
names(p2g)<-pth.names;
    
# limit genes to tested pathways
gn<-sort(unique(unlist(p2g, use.names=FALSE)));
if (prms$limit.gene) g2v<-g2v[names(g2v) %in% gn];

# limit variants to tested genes
v<-sort(unique(unlist(g2v, use.names=FALSE)));
if (prms$limit.variant) {
  geno<-geno[rownames(anno) %in% v, ];
  anno<-anno[rownames(anno) %in% v, ];
}

# no variants on chromosome X/Y
if (!prms$sex.chromosomes) {
  geno<-geno[anno$chr %in% as.character(1:22), ];
  anno<-anno[anno$chr %in% as.character(1:22), ];
}

cat('Analyzing', nrow(geno), 'variants.\n');
cat('Analyzing', length(g2v), 'genes.\n');
cat('Analyzing', length(p2g), 'pathways.\n');

    
    #####################################################
    ### A FUNCTION TO PERFORM ONE SUBGROUP COMPARISON ###
    #####################################################
    
#########################################################################################################################################################
#########################################################################################################################################################
doComp<-function(smps, path, prms, geno, anno, g2v, p2g) {
  num.threads<-prms$num.threads;
  skat.weight<-prms$skat.weight;
  ar.only<-prms$ar.only;
  
  smps<-lapply(smps, function(s) s[s %in% colnames(geno)]);
  smpn<-sapply(smps, length);
  
  # create a subfolder for results
  path.out<-paste(path, '/', names(smps)[1], '_vs_', names(smps)[2], sep='');
  if (!file.exists(path.out)) dir.create(path.out);
  cat("Results will be written to", path.out, '\n');
  
  
  # data for comparison
  g<-geno[, c(smps[[1]], smps[[2]])];
  g[g>2 | g<0]<-NA;
  
  #########################################################
  ### count allele frequency
  cat('Calculating allele frequency of both groups ...\n');
  if (num.threads>1 & nrow(g)>100000) {
    ct0<-countGenotypeCl(g, smps[[1]], 2, num.threads);
    ct1<-countGenotypeCl(g, smps[[2]], 2, num.threads);
  } else {
    #####################################
    ct0<-countGenotype(g, smps[[1]], 2);#
    ct1<-countGenotype(g, smps[[2]], 2);#
    #####################################
  }

  
  ########################################################################################
  ### swap genotype if alternative allele is more common in controls
  ind<-which(ct0[,3]>ct0[,1]);
  if (length(ind)>0) { # if found such variants, flip the genotype
    g[ind, ]<-2-g[ind, ];
    ct0[ind, ]<-ct0[ind, 3:1];
    ct1[ind, ]<-ct1[ind, 3:1];  
  }
  a<-anno;
  freq<-cbind(ct0, ct1);
  colnames(freq)<-paste(rep(names(smps), each=3) , rep(c('AA', 'AB', 'BB'), 2), sep='_')
  save(freq, file=paste(path.out, 'allele_frequency_all_variants.rdata', sep='/'));
  
  
  #####################################################################################################
  ### remove variants without AB nor BB calls
  n<-rowSums(freq[, c(2, 3, 5, 6)], na.rm=TRUE);
  ind<-which(n==0);
  if (length(ind) > 0) {
    g<-g[-ind, , drop=FALSE];
    a<-anno[rownames(anno) %in% rownames(g), , drop=FALSE][rownames(g), ,drop=FALSE];
    freq<-freq[-ind, , drop=FALSE];
    cat("Removed", length(ind), 'variants having no alternative allele calls in all tested samples.\n');
  } 
  if (ar.only) { # use autosomal recessive variants only
    ar<-freq[, 3]+freq[,6];
    g<-g[ar>0, ];
    g[!is.na(g) & g==1]<-0;
    a<-a[ar>0, ];
    freq<-freq[ar>0, ];
    cat("Removed", length(ar[ar==0]), 'non-autosomal recessive variants.\n');
  }
  maf<-50*cbind((freq[,2]+freq[,3]*2)/rowSums(freq[,1:3]), (freq[,5]+freq[,6]*2)/rowSums(freq[,4:6]), (freq[,2]+freq[,3]*2+freq[,5]+freq[,6]*2)/rowSums(freq[, 1:6]));
  colnames(maf)<-paste('MAF_', c(names(smps), 'All'), sep='');
  freq<-cbind(freq, maf);
  
  
  ############################################################################################
  ### variant-level trend test
  cat('Running Trend test on variants ...\n');
  
  #####################################################################################  
  stat<-runTrendCl(g, as.factor(rep(0:1, smpn)), values=0:2, num.threads=num.threads);#
  #####################################################################################
  
  stat<-cbind(freq[names(stat), ], P_Trend=stat);
  save(stat, file=paste(path.out, 'stat_trend.rdata', sep='/'));
  write.csv(cbind(stat, a[rownames(stat), ]), paste(path.out, 'stat_trend.csv', sep='/'));
  
  
  ###################################################################################################
  #### variant-level fisher test
  cat('Running Fisher test on variants ...\n');
  
  #############################################################
  stat.fsh<-runFisherCl(freq[, 1:6], num.threads=num.threads);#
  #############################################################
  
  ### allele model
  nm<-rep(names(smps), each=2);
  nm<-c(paste(nm, c('A', 'B', 'A', 'B'), sep='_'), 'Odds_Ratio', 'P_Fisher');
  stat<-stat.fsh[, grep('Allelic', colnames(stat.fsh))];
  colnames(stat)<-nm;
  save(stat, file=paste(path.out, 'stat_fisher_allele.rdata', sep='/'));
  write.csv(cbind(stat, a[rownames(stat), ]), paste(path.out, 'stat_fisher_allele.csv', sep='/'));
  
  ### recessive model
  stat<-stat.fsh[, grep('Recessive', colnames(stat.fsh))];
  colnames(stat)<-nm;
  save(stat, file=paste(path.out, 'stat_fisher_recessive.rdata', sep='/'));
  write.csv(cbind(stat, a[rownames(stat), ]), paste(path.out, 'stat_fisher_recessive.csv', sep='/'));
  
  ### dominant model
  stat<-stat.fsh[, grep('Dominant', colnames(stat.fsh))];
  colnames(stat)<-nm;
  save(stat, file=paste(path.out, 'stat_fisher_dominant.rdata', sep='/'));
  write.csv(cbind(stat, a[rownames(stat), ]), paste(path.out, 'stat_fisher_dominant.csv', sep='/'));
  
  
  #######################################################################################################################################################
  ### gene-level SKAT test
  cat('Running gene-level SKAT test ...\n');  
  gs<-rep(names(g2v), sapply(g2v, length)); 
  vs<-unlist(g2v, use.names=FALSE); 
  g2v<-split(vs[vs %in% rownames(g)], gs[vs %in% rownames(g)]); 
  #g2v<-g2v[sapply(g2v, length)>1]; 
  
  ####################################################################################################################################################
  stat<-runSkatCl(ids=g2v, gt=g, samples0=smps[[1]], samples1=smps[[2]], num.threads=num.threads, trimming=FALSE, weights=skat.weight, verbose=TRUE);#
  ####################################################################################################################################################
  
  save(stat, file=paste(path.out, 'stat_gene_skat.rdata', sep='/')); 
  write.csv(stat, paste(path.out, 'stat_gene_skat.csv', sep='/'));
  
  
  ################################################################################################################################################
  ### pathway-level SKAT test
  cat('Running pathway-level SKAT test ...\n');
  p2v<-lapply(p2g, function(g) unlist(g2v[g[g %in% names(g2v)]], use.names=FALSE));
  p2v<-lapply(p2v, unique);
  p2v<-lapply(p2v, sort);
  
  ######################################################################################################################################
  stat<-runSkatCl(ids=p2v, gt=g, samples0=smps[[1]], samples1=smps[[2]], num.threads=num.threads, trimming=FALSE, weights=skat.weight);#
  ######################################################################################################################################
  
  save(stat, file=paste(path.out, 'stat_pathway_skat.rdata', sep='/'));
  write.csv(stat, paste(path.out, 'stat_pathway_skat.csv', sep='/'));
  
  inputs<-list(samples=smps, parameters=prms, genes=g2v, pathways=p2g);
  save(inputs, file=paste(path.out, 'inputs.rdata', sep='/'));
  
  path.out;
}
#########################################################################################################################################################
#########################################################################################################################################################

    
    ######################
    ### DO COMPARISONS ###
    ######################

####################################################################################
fld<-sapply(grps, function(grps) doComp(grps, path, prms, geno, anno, g2v, p2g));###
####################################################################################

write.table(fld, row=FALSE, col=FALSE, sep='\t', qu=FALSE, 
          paste("/nas/is1/PCGC/gVCF6703/minimal_requirement_variants/subgroup_comparisons/runs/run_log/", run.time, '.txt', sep=''));
    
