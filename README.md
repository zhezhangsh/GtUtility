#### Utility functions for processing and analysis of genotype data

#### To load library within your R console
```
library(devtools);
install_github("zhezhangsh/GtUtility");
library(GtUtility);
```


#### To run the runSmpComp.r script
```
# Copy both 'runSmpComp.r' and 'runSmpComp.yaml' files into your working folder
# Edit the 'runSmpComp.yaml' to specify inputs and output path
### avoid putting any space into the empty lines 
### do not change the yaml file name, a copy of it will be automatically saved to output folder
# Never change the 'runSmpComp.r' file
# Run the following shell command
### replace "./" with the full path of both files
### replace "/usr/local/bin/Rscript" with the full path where Rscript was installed

/usr/local/bin/Rscript ./runSmpComp.r ./runSmpComp.yaml

```
