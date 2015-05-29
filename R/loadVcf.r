# Import all variants or a subset of them from a VCF filecomplete or a subset of a VCF file and save them in an R object
loadVcf<-function(fn.in, fn.out, genome.name, regions=NA, split.regions=FALSE, verbose=FALSE) {
# fn.in         Full path and name of the VCF file 
# fn.out        Prefix of the output file(s)
# genome.name   Genome name, such as 'hg19' and 'mm10'
# regions       If not NA, only import variants within given regions defined by a GRange object; make sure the chromosome names match the names in VCF file
# split.region  If TRUE, import variants of each region and save them in a separate file
  
  library(Rsamtools);
  library(VariantAnnotation);
  
  # compress VCF file into bgzip format if not done yet
  if (!grepl('gz$', fn.in)) {
    f<-paste(fn.in, '.bgz', sep='');
    if (file.exists(f)) fn.in<-f else fn.in<-bgzip(fn.in);
  }
  
  
  # Index VCF file for fast access
  if (!file.exists(paste(fn.in, '.tbi', sep=''))) {
    if (verbose) print('Indexing VCF file...');
    indexTabix(fn.in, format='vcf');
  } 
  
  chr.names<-headerTabix(fn.in)$seqnames; # seqnames in VCF file
  if (!identical(regions, NA)) regions<-regions[names(regions) %in% chr.names]; # remove chr names not in VCF file
  if (length(regions)==0) regions<-NA; # if none of the regions are in VCF file, import all file
  
  tab<-TabixFile(fn.in);
  if (identical(regions, NA)) { # read in all VCF and save variants in one R file
    if (verbose) print('Reading in all variants ...')
    vcf<-readVcf(tab, genome.name);
    fnm<-paste(fn.out, '.rdata', sep='');
    save(vcf, file=fnm);
    fnm;
  } else {
    if (!split.regions) { # read in all variant in given regions and save all variant in one R file
      if (verbose) print('Reading in variants in given regions ...');
      vcf<-readVcf(tab, genome.name, param=regions);
      fnm<-paste(fn.out, '.rdata', sep='');
      save(vcf, file=fnm);
      fnm;
    } else { # 
      nm<-paste(names(regions), '_', start(regions), '-', end(regions), sep='');
      fnm<-paste(fn.out, '_', nm, '.rdata', sep='');
      names(regions)<-nm;
      lapply(nm, function(nm) {
        if (verbose) print(paste('Reading variants within', nm, '...'));
        vcf<-readVcf(tab, genome.name, param=regions[nm]);
        save(vcf, file=paste(fn.out, '_', nm, '.rdata', sep=''));
      });
      fnm;
    }
  }
  
  
}
