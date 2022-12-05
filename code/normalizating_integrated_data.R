# reading the input datasets 
integrated_data <- read_tsv("../results/integrated_data.tsv")
print("> Loaded importation data: integrated_data.tsv")
sample_data <- read_tsv("../data/samples_data_processed.tsv")
print("> Loaded importation data: samples_data_processed.tsv")

# removing tissues with less than 10 samples
list_tissue <- integrated_data %>% group_by(Tissue,Dataset,Sample) %>%
  distinct(Sample)
list_tissue <- list_tissue %>% group_by(Tissue) %>% count(Tissue)

cancer_tissue <- list_tissue %>%
  filter(!grepl("Healthy", Tissue)) %>%
  filter( n >= 10)

healthy_tissue <- list_tissue %>%
  filter(grepl("Healthy", Tissue))

integrated_data <- integrated_data %>% 
  filter(Tissue %in% c(cancer_tissue$Tissue,healthy_tissue$Tissue))

# replacing 0 with NA
integrated_data <- integrated_data %>%
  mutate(Intensity = ifelse(Intensity == 0, NA, Intensity))
hist(integrated_data$Intensity, breaks = 20)

# log 2
normalized_data <- integrated_data %>%
  mutate(Intensity = log2(Intensity))
hist(normalized_data$Intensity, breaks = 20)

# z-score normalization 
normalized_data <- normalized_data %>% 
  group_by(Dataset, Sample) %>% 
  mutate(Mean_group = mean(Intensity, na.rm=TRUE)) %>%
  mutate(Sd_group = sd(Intensity, na.rm=TRUE)) %>%
  mutate(Z_score_group = (Intensity - mean(Intensity, na.rm= TRUE)) / sd(Intensity, na.rm=TRUE)) 
summary(normalized_data)

normalized_data <- normalized_data %>% 
  select(-c(Mean_group, Sd_group))
summary(normalized_data)

# writing the results as a tsv file
write_tsv(normalized_data,"../results/normalized_data.tsv")
print("> Exportation data: normalized_data.tsv")
