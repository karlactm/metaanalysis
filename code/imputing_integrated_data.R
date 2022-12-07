# reading the dataset from its standard location and using its standard name
normalized_data <- read_tsv("../results/normalized_data.tsv")
print("> Loaded importation data: normalized_data.tsv")

# grouping all healthy samples 
healthy_tissue <- normalized_data %>% filter(str_detect(Tissue, "Healthy")) %>% distinct(`Id Tissue`)
healthy_tissue <- pull(healthy_tissue)
normalized_data <- normalized_data %>%
  mutate(Tissue = ifelse(`Id Tissue` %in% healthy_tissue, "Healthy Tissue", Tissue))

# calculating min, mean and sd
min_dt <- min(normalized_data$Z_score_group, na.rm = T)
mean_dt <- mean(normalized_data$Z_score_group, na.rm = T)
sd_dt <- sd(normalized_data$Z_score_group, na.rm = T)

# producing the table with the statitics results
statistics <- normalized_data %>% 
  group_by(`Gene names`) %>%
  dplyr::summarise(
    n = n(),
    min =  min(Z_score_group, na.rm = T),
    mean =  mean(Z_score_group, na.rm = T),
    median =  median(Z_score_group, na.rm = T),
    sd =  sd(Z_score_group, na.rm = T),
    nas = sum(is.na(Z_score_group)),
    notnas = (n - nas)
  )

# imputating values
imputed_data <- normalized_data %>% 
  group_by(`Gene names`) %>%
  inner_join(statistics, by = "Gene names") %>%
  mutate(Z_score_group =  ifelse(is.na(Z_score_group),
                                 ifelse((notnas) > (ceiling((5*n)/100)),
                                        rtruncnorm(sum(is.na(Z_score_group)), 
                                                   a=-Inf, 
                                                   b=min, 
                                                   mean=median, 
                                                   sd=sd),
                                        rtruncnorm(sum(is.na(Z_score_group)), 
                                                   a=-Inf, 
                                                   b=min_dt, 
                                                   mean=(mean_dt - (1.8*sd_dt)), 
                                                   sd=(sd_dt * 0.3))),
                                 Z_score_group)) %>%
  select(-c(n, mean, min, sd, nas, notnas))

# writing the results as a tsv file
write_tsv(imputed_data, "../results/imputed_data.tsv")
print("> Exportation data: imputed_data.tsv")
