load_required_packages <- function(x) {
  for(i in x){
    if(!require(i, character.only = TRUE)){
      install.packages(i, dependencies = TRUE)
      require(i, character.only = TRUE)
    }
  }
}

# first, we need to install/load all required packages
load_required_packages(c("this.path","tidyverse","truncnorm","arules"))
#load_required_packages(c("this.path","tidyverse"))

current_dir <- this.path::this.dir()

# loading needed underlying functions
source(file.path(current_dir, "processing_protein_groups.R"))

source(file.path(current_dir, "integrating_protein_groups.R"))

source(file.path(current_dir, "normalizating_integrated_data.R"))

source(file.path(current_dir, "imputating_integrated_data.R"))

source(file.path(current_dir, "categorizing_integrated_data.R"))
