file_list_free <- list.files(file.path(current_dir,"../data/free_label/"))
file_list_silac <- list.files(file.path(current_dir,"../data/silac/"))

process_data_free = function(data, output_name) {
  
  # removing potential contaminant and Reverse and select interest columns
  protein_groups <- data %>%  
    mutate(Reverse = replace_na(Reverse, "-"),
           `Potential contaminant` = replace_na(
             `Potential contaminant`, "-")) %>%
    filter(Reverse != "+" & `Potential contaminant` != "+") %>%
    select (`Protein IDs`, `Protein names`,`Gene names`, `Fasta headers`,
            `Peptide counts (unique)`, `Sequence coverage [%]`, 
            Score, starts_with("Intensity "))
  
  protein_groups_processed <-  protein_groups
  
  # identifying the gene name
  protein_groups_processed$`Gene names` <- mapply(function(x,y) {
    fasta <- unlist(str_split(y, "\\|"))[3]
    if (is.na(x)) {
      if (!str_detect(fasta, "Immunoglobulin")) {
        gene_name <- unlist(str_split(str_extract(fasta, "\\w*GN=\\w+"), "GN="))[2]
        return(gene_name)
      } else {
        return(NA)
      }
    } else {
      return(x)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Fasta headers`)
  
  # removing missing values 
  protein_groups_processed <- protein_groups_processed %>% drop_na(`Gene names`)
  
  # processing the peptides number
  protein_groups_processed$`Peptide counts (unique)` <- mapply(function(x, y) {
    if (str_detect(y, "([0-9]+;[0-9])+")) {
      gene_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      counts <- table(peptides_count)
      peptides_count <- as.integer(peptides_count)
      value <- counts[names(counts) == max(peptides_count)]
      peptides_count <- peptides_count[1:value]
      peptides_count <- paste0(peptides_count, collapse = ";")
      return(peptides_count) 
    } else {
      return(y)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Peptide counts (unique)`)
  
  
  # processing the genes number
  protein_groups_processed$`Gene names` <- mapply(function(x, y) {
    if (str_detect(x, "([a-zA-Z0-9]+;[a-zA-Z0-9])+")) {
      gene_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      length(gene_name) <- length(peptides_count) 
      gene_name <- paste0(gene_name, collapse = ";")
      return(gene_name)
    } else {
      return(x)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Peptide counts (unique)`)
  
  # processing the proteins number
  protein_groups_processed$`Protein IDs` <- mapply(function(x, y) {
    if (str_detect(x, "([a-zA-Z0-9]+;[a-zA-Z0-9])+")) {
      protein_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      length(protein_name) <- length(peptides_count) 
      protein_name <- paste0(protein_name, collapse = ";")
      return(protein_name)
    } else {
      return(x)
    }
  }, protein_groups_processed$`Protein IDs`, protein_groups_processed$`Peptide counts (unique)`)
  
  # removing missing values 
  protein_groups_processed <- protein_groups_processed %>% drop_na(`Gene names`)
  
  # separating a gene and a protein by line
  protein_groups_processed <- protein_groups_processed %>% 
    separate_rows( `Protein IDs`,`Gene names`, `Peptide counts (unique)`, sep=";")
  
  # sorting by peptides number and removing lines with duplicate genes  
  protein_groups_processed <- protein_groups_processed %>% 
    filter(`Peptide counts (unique)` != 0) %>%
    arrange(desc(as.numeric(`Peptide counts (unique)`))) %>%
    distinct(`Gene names`, .keep_all = TRUE)
  
  # removing missing values 
  protein_groups_processed <- protein_groups_processed %>% filter(`Gene names` != "NA")
  
  # writing the results as a tsv file
  write_tsv(protein_groups_processed, file.path(current_dir,"../results/", output_name))
  
}

process_data_silac = function(name, data, output_name) {
  
  col_name <- paste0("Intensity H ", name)
  
  # removing potential contaminant and Reverse and select interest columns
  protein_groups <- data %>%
    mutate(Reverse = replace_na(Reverse, "-"), 
           `Potential contaminant` = replace_na(
             `Potential contaminant`, "-")) %>%
    filter(Reverse != "+" & `Potential contaminant` != "+") %>%
    select (`Protein IDs`, `Protein names`,`Gene names`, `Fasta headers`, 
            `Peptide counts (unique)`, `Sequence coverage [%]`, 
            Score, `Intensity H`,starts_with("Intensity L ")) %>%
    rename(!!quo_name(col_name) := `Intensity H`)
  
  protein_groups_processed <-  protein_groups
  
  # identifying the gene name
  protein_groups_processed$`Gene names` <- mapply(function(x,y) {
    fasta <- unlist(str_split(y, "\\|"))[3]
    if (is.na(x)) {
      if (!str_detect(fasta, "Immunoglobulin")) {
        gene_name = unlist(str_split(str_extract(fasta, "\\w*GN=\\w+"), "GN="))[2]
        return(gene_name)
      } else {
        return(NA)
      }
    } else {
      return(x)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Fasta headers`)
  
  # removing missing values
  protein_groups_processed <- protein_groups_processed %>% drop_na(`Gene names`)
  
  # processing the peptides number
  protein_groups_processed$`Peptide counts (unique)` <- mapply(function(x, y) {
    if (str_detect(y, "([0-9]+;[0-9])+")) {
      gene_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      # Conta o número de elementos igual ao maior valor de peptides_count
      counts <- table(peptides_count)
      peptides_count <- as.integer(peptides_count)
      value <-counts[names(counts) == max(peptides_count)]
      peptides_count <- peptides_count[1:value]
      peptides_count <- paste0(peptides_count, collapse = ";")
      return(peptides_count) 
    } else {
      return(y)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Peptide counts (unique)`)
  
  # processing the genes number
  protein_groups_processed$`Gene names` <- mapply(function(x, y) {
    if (str_detect(x, "([a-zA-Z0-9]+;[a-zA-Z0-9])+")) {
      gene_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      length(gene_name) <- length(peptides_count) 
      gene_name <- paste0(gene_name, collapse = ";")
      return(gene_name)
    } else {
      return(x)
    }
  }, protein_groups_processed$`Gene names`, protein_groups_processed$`Peptide counts (unique)`)
  
  # processing the proteins number
  protein_groups_processed$`Protein IDs` <- mapply(function(x, y) {
    if (str_detect(x, "([a-zA-Z0-9]+;[a-zA-Z0-9]+)")) {
      protein_name <- unlist(str_split(x, ";"))
      peptides_count <- unlist(str_split(y, ";"))
      length(protein_name) <- length(peptides_count) 
      protein_name <- paste0(protein_name, collapse = ";")
      return(protein_name)
    } else {
      return(x)
    }
  }, protein_groups_processed$`Protein IDs`, protein_groups_processed$`Peptide counts (unique)`)
  
  # removing missing values
  protein_groups_processed <- protein_groups_processed %>% drop_na(`Gene names`)
  
  # Separating a gene and a protein by line 
  protein_groups_processed <- protein_groups_processed %>% 
    separate_rows( `Protein IDs`,`Gene names`, `Peptide counts (unique)`, sep=";")
  
  # sorting by peptides number and removing lines with duplicate genes  
  protein_groups_processed <- protein_groups_processed %>%
    arrange(desc(as.numeric(`Peptide counts (unique)`))) %>%
    distinct(`Gene names`, .keep_all = TRUE)
  
  # removing missing values
  protein_groups_processed <- protein_groups_processed %>% filter(`Gene names` != "NA")
  
  # writing the results as a tsv file
  write_tsv(protein_groups_processed, file.path(current_dir,"../results/", output_name))
  
}

for (i in 1:length(file_list_free)){
  protein_groups <- read_tsv(file.path("../data/free_label/", file_list_free[i]))
  print(paste0("> Loaded importation data: ", file_list_free[i]))
  process_data_free(protein_groups,  file_list_free[i])
}

for (i in 1:length(file_list_silac)){
  protein_groups <- read_tsv(file.path("../data/silac/", file_list_silac[i]))
  print(paste0("> Loaded importation data: ",file_list_silac[i]))
  name <- gsub(file_list_silac[i], pattern = "\\s+", replacement = "")
  name <- unlist(str_split(unlist(str_split(name, "\\."))[1], "s_"))[2]
  if (str_detect(name, "\\_")) {
    name <- gsub(name, pattern = "_", replacement = "") 
  } 
  process_data_silac(name, protein_groups, file_list_silac[i])
}


