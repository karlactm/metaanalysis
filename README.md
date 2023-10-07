# Overview

This repository maintains the implementation of our meta-analysis that integrates public protemic data from human tissues and performs exploratory data analysis for the discovery of abundant proteins in tumoral tissue, in order to enhance the identification of relevant genes as molecular targets for cancer also detected at the protein level.

## How to use?

This project is also publicly available in CodeOcean here. Hence, you can execute our code and obtain the results. Our main script invokes five other scripts developed for processing, cleaning, integrating and transforming data. The scripts developed for each of the steps can be seen below:

1. **Data Cleaning:** 
- processing_protein_groups.R: removes proteins considered as contaminants, as well as proteins with missing gene name and with duplicate gene name. Besides, it processes the peptide number and the gene number.
3. **Data Integration:**
- integrating_protein_groups.R: changes the format of the protein group tables, adds descriptive columns about the biological samples and integrates the protein group tables.
4. **Data Transformation:**
- normalizing_integrated_data.R: transforms (using log2) and then normalizes (using Z-score) the intensity measurements.
* imputing_integrated_data.R: replaces missing values with imputed values using a Random Tail Imputation (RTI).
+ converting_integrated_data.R: calculates, for each protein, the median of the intensity values for all samples of the same tissue and converts the median into three types of categories (low, medium and high).

CRAN packages [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html), [truncnorm](https://cran.r-project.org/web/packages/truncnorm/), [arules](https://cran.r-project.org/web/packages/arules/) and [this.path](https://cran.r-project.org/web/packages/this.path/index.html) installation are required.
