get_dataset_name <- function(name) {
  dataset_name <- str_replace(str_replace(str_split(
    name, "_", simplify = T)[, 2], ".txt", ""), " ", "")
  return(dataset_name)
}

# reading the dataset from its standard location and using its standard name
sample_data <- read_tsv("../data/samples_data_processed.tsv")
print("> Loaded importation data: samples_data_processed.tsv")

quantification_dt <- sample_data %>% distinct(Dataset, `Label-free / SILAC`)
quantification_dt <- quantification_dt %>%
  mutate(`Label-free / SILAC` = ifelse(
    `Label-free / SILAC` == "Label-free", "FREELABEL", `Label-free / SILAC` ))

file_list <- list.files(file.path(current_dir,"../results/"))
file_list <- tibble("File name" = file_list)
file_list <- file_list %>% filter(grepl("proteinGroups", `File name`)) 
file_list <- file_list %>%
  mutate(Dataset = str_to_upper(str_replace(str_replace(str_split(
    `File name`, "_", simplify = T)[, 2], ".txt", ""), " ", ""))) 
file_list <- file_list %>%
  inner_join(quantification_dt, by = "Dataset")

file_list_free <- file_list %>% filter(`Label-free / SILAC` == "FREELABEL") %>% select(`File name`) 
file_list_free <- as.vector(file_list_free$`File name`)
file_list_silac <- file_list %>% filter(`Label-free / SILAC` == "SILAC") %>% select(`File name`) 
file_list_silac <- as.vector(file_list_silac$`File name`)
freelabel_data <- tibble()

# processing label-free quantitative data
for (i in 1:length(file_list_free)){
  
  #getting dataset name
  dataset_name <- get_dataset_name(file_list_free[i])
  
  print(paste0("> Loaded importation data: ",file_list_free[i]))
  data <- read_tsv(file.path("../results/",file_list_free[i]))
  
  data <- data %>% select(`Gene names`, starts_with("Intensity "))
  names(data) <- gsub(x = names(data), pattern = "Intensity.", replacement = "")
  
  data_longer <- data %>%
    pivot_longer(!`Gene names`, names_to = "Sample", values_to = "Intensity") %>%
    mutate(`Dataset` = toupper(dataset_name)) 
  
  data_longer <- inner_join(data_longer, sample_data, by = c("Sample", "Dataset")) 
  
  data_longer <- data_longer %>%
    relocate(Dataset, Description, Tissue, `Id Tissue`, Sample, `Sample Name`, `Label-free / SILAC`)
  
  # integrating data
  freelabel_data <- bind_rows(freelabel_data, data_longer)
  
}
# removing Library 
freelabel_data <- freelabel_data %>% filter(Sample != "Library")

# removing Testis
freelabel_data <- freelabel_data %>% filter(Sample != "Testis")

# renaming Immune T
freelabel_data <- freelabel_data %>%
  mutate(Tissue = if_else(str_detect(Tissue, "Immune T"), "Healthy Immune T", Tissue))

#file_list_silac <- list.files(file.path(current_dir,"/result/silac/"))
#file_list_silac <- file_list_silac[!is.na(file_list_silac)]
silac_data <- tibble()

# processing silac quantitative data
for (i in 1:length(file_list_silac)){
  
  dataset_name <- get_dataset_name(file_list_silac[i])
  
  data <-  read_tsv(file.path("../results/",file_list_silac[i]))
  print(paste0("> Loaded importation data: ",file_list_silac[i]))
  
  data <- data %>% select(`Gene names`, starts_with("Intensity "))
  
  names(data) <- gsub(x = names(data), pattern = "Intensity ", replacement = "")
  
  data_longer <- data %>%
    pivot_longer(!`Gene names`, names_to = "Sample", values_to = "Intensity") %>%
    mutate(`Dataset` = toupper(dataset_name)) 
  
  data_longer <- data_longer %>%  
    mutate(Sample = toupper(Sample)) %>%
    mutate(Sample = gsub(Sample, pattern = "\\L ", replacement = "")) %>%
    mutate(Sample = gsub(Sample, pattern = "\\H [[:alnum:]]+", replacement = "Intensity H")) 
  
  data_longer <- inner_join(data_longer, sample_data, by = c("Sample", "Dataset")) 
  
  data_longer <- data_longer %>%
    relocate(Dataset, Description, Tissue, `Id Tissue`, Sample, `Sample Name`, `Label-free / SILAC`)
  
  # integrating data
  silac_data <- bind_rows(silac_data, data_longer)
  
}

# integrating data
integrated_data <- bind_rows(freelabel_data, silac_data)

# writing the raw results as a tsv file
write_tsv(integrated_data, "../results/integrated_data.tsv")
print("> Exportation data: integrated_data.tsv")



