# reading the input datasets 
imputed_data <- read_tsv("../results/imputed_data.tsv")
print("> Loaded importation data: imputed_data.tsv")
ctas_transcriptomic_level <- read_tsv("../data/ct_antigens.tsv")
print("> Loaded importation data: ct_antigens.tsv")

# selecting interest columns
imputed_data <- imputed_data %>% select(c(Dataset, Sample, Tissue, `Gene names`, Z_score_group))

# calculating median
median_data_grouped  <-  imputed_data %>% group_by(`Gene names`, Tissue) %>%
  summarise("Median Z-score" = median(Z_score_group, na.rm= TRUE), .groups = "keep")

# applying equal frequency discretization
table(discretize(median_data_grouped $`Median Z-score`, method = "frequency", breaks = 3))
hist(median_data_grouped $`Median Z-score`, breaks = 20, main = "Equal Frequency")
abline(v = discretize(median_data_grouped $`Median Z-score`, method = "frequency", breaks = 3,
                      onlycuts = TRUE), col = "red")

discretize_data <- discretizeDF(median_data_grouped ,
                                            default = list(method = "frequency",
                                                           breaks = 3,labels = c("low", "medium","high")))

discretize_data <- discretize_data %>%
  rename(Class = `Median Z-score`) %>% 
  left_join(intensity_z_score_data, by = c("Gene names", "Tissue")) %>%
  relocate(`Gene names`, Tissue, `Median Z-score`, Class)

# ctas observed at the protein level
ctas_proteomic_level <- intensity_z_score_data_disc %>% 
  inner_join(ctas_transcriptomic_level, by = "Gene names") 

# writing the results as a tsv file  
write_tsv(ctas_proteomic_level, "../results/ctas_proteomic_level.tsv")
print("> Exportation data: ctas_proteomic_level.tsv")
