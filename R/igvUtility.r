# Generate the XML file for 
generateIgvXml<-function(fn.yaml) {
  yaml<-yaml::yaml.load_file(fn.yaml);
  
  fn<-yaml$xml.path;
  bams<-yaml$bams;
  
  igvHeaderLine(yaml$project.name, fn);
  ln<-sapply(names(bams), function(nm) {
    igvCategoryLine(nm, bams[[nm]], file.name=fn);    
  });
  igvFooterLine(fn);
  
  fn;
}

generateIgvYaml<-function(project.name, url.path, xml.path, category2library, fn.yaml) {
  ln<-paste('project.name', project.name, sep=': ');
  ln<-c(ln, paste('url.path', url.path, sep=': '));
  ln<-c(ln, paste('xml.path', xml.path, sep=': '));
  ln<-c(ln, 'bams:');
  
  url.path<-sub('/$', '', url.path);
  
  lns<-lapply(names(category2library), function(nm) {
    l<-paste('    ', names(category2library[[nm]]), ': ', paste(url.path, category2library[[nm]], sep='/'), sep='');
    c(paste('  ', nm, ':', sep=''), l);
  });
  
  writeLines(c(ln, unlist(lns)), fn.yaml);
}

igvHeaderLine<-function(project.name, file.name=NA) {
  file.name<-file.name[1];
  ln<-'<?xml version="1.0" encoding="UTF-8"?>';
  ln[2]<-paste('<Global name="', project.name, '"  infolink="http://www.broadinstitute.org/igv/" version="1">', sep='');
  
  if (!identical(NA, file.name)) writeLines(ln, file.name);
  
  ln;
}

igvCategoryLine<-function(category.name, samples, indent=1, file.name=NA) {
  indent<-max(1, indent);
  file.name<-file.name[1];
  
  ln<-paste(rep('\t', indent), collapse='');
  ln<-paste(ln, '<Category name="', category.name, '">', sep='');
  if (!identical(NA, file.name)) write(ln, file.name, append=TRUE);
  
  lns<-sapply(names(samples), function(nm) igvSampleLine(nm, samples[nm], file.name=file.name));
  if (!identical(NA, file.name)) write(paste(paste(rep('\t', indent), collapse=''), '</Category>', sep=''), file.name, append=TRUE);
  
  lns;
}

igvSampleLine<-function(sample.name, bam.path, indent=2, file.name=NA) {
  indent<-max(1, indent);
  file.name<-file.name[1];
  
  ln<-paste(rep('\t', indent), collapse='');
  ln<-paste(ln, '<Resource name="', sample.name, '" path="', bam.path, '" />', sep='');
  
  if (!identical(NA, file.name)) write(ln, file.name, append=TRUE);
  
  ln;
}

igvFooterLine<-function(file.name=NA) {
  file.name<-file.name[1];
  
  ln<-'</Global>';
  if (!identical(NA, file.name)) write(ln, file.name, append=TRUE);
  
  ln;
}


