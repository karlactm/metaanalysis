# Overview

This repository maintains the implementation of our meta-analysis that integrates and analyzes independent proteomics datasets from tumor and healthy human tissues to verify protein levels of CTAs predicted by transcriptomics, and evaluate their potential as biomarker candidates.

## How to use?

This project is also publicly available in CodeOcean here. Hence, you can execute our code and obtain the results. Our main script invokes five other scripts developed for processing, cleaning, integrating and transforming data. The scripts developed for each of the steps can be seen below:

1. **Data Cleaning:** 
- processing_protein_groups.R: removes proteins considered as contaminants and false-positives, as well as proteins with missing gene names. In addition, for protein groups containing more than one protein and one
or more genes,  selects the protein entry number and its coding gene from the group member with the most number of supporting unique peptides.
3. **Data Integration:**
- integrating_protein_groups.R: changes proteingroups.txt file structure, classify and group individual samples according to their tissue of origin, adds descriptive attributes for each sample and integrates proteingroups.txt files.
4. **Data Transformation:**
- normalizing_integrated_data.R: transforms (using log2) and then normalizes (using Z-score) the intensity measurements.
* imputing_integrated_data.R: replaces missing values with imputed values using a Random Tail Imputation (RTI).
+ converting_integrated_data.R: categorizes the values of median protein abundance as low, medium, high, or absent.

CRAN packages [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html), [truncnorm](https://cran.r-project.org/web/packages/truncnorm/), [arules](https://cran.r-project.org/web/packages/arules/) and [this.path](https://cran.r-project.org/web/packages/this.path/index.html) installation are required.
