# Overview

This repository maintains the implementation of our meta-analysis that integrate public protemic data from human tissues and perform a exploratory data analysis through the discovery of abundant proteins in tumoral tissue, in order to enhance the identification genes as molecular targes for cancer at the protein level.

## How to use?

This project is also publicly available in CodeOcean here. Hence, you can execute our code and obtain the results. Our main script invokes five scripts developed for the steps  for processing the data were cleaning, integrating and transforming data. The scripts developed for each of the steps can be seen below: 

1. **Cleaning data:** 
- processing_protein_groups.R: removes proteins considered as contaminants, as well as proteins with the missing gene name and with the duplicate gene name. Besides processes the peptides number and the genes number.
3. **Integrating data:**
- integrating_protein_groups.R: changes the format of the protein group tables, adds descriptive columns about the biological samples and integrates the protein group tables.
4. **Transformating data:**
- normalizating_integrated_data.R: transforms (using log2) and then normalizes (using Z-score) the intensity measurements.
* imputating_integrated_data.R: replaces missing values with imputed values.
+ categorizing_integrated_data.R: calculates, for each protein, the median of the intensity values for all samples of the same tissue and converts the median into three types of categories (low, medium and high).

In house packages [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html), [truncnorm](https://cran.r-project.org/web/packages/truncnorm/), [arules](https://cran.r-project.org/web/packages/arules/) and [this.path](https://cran.r-project.org/web/packages/this.path/index.html) installation are required.
