#### Utility functions for processing and analysis of genotype data

#### To load library within your R console
```
library(devtools);
install_github(host="github.research.chop.edu/api/v3",repo="BiG/GtUtility");
library(GtUtility);
```


#### To run the runSmpComp.r script
```
# Copy both 'runSmpComp.r' and 'runSmpComp.yaml' files into your working folder
# Edit the 'runSmpComp.yaml' to specify inputs and output path
### !!! avoid putting any space into the empty lines !!!
### !!! do not change the yaml file name, a copy of it will be automatically saved to output folder !!!
# Run the following shell command
### specify the full path if the files are not within your working directory
### change the path to Rscript executible if it's installed elsewhere
/usr/local/bin/Rscript ./runSmpComp.r ./runSmpComp.yaml
```
