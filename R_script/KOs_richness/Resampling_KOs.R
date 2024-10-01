library(iNEXT)
library(reshape2)
library("plyr")
library("dplyr")
library("tidyverse")
library("readr")
library("tidyr")

# Function to process and analyze abundance data
process_and_analyze <- function(file_path, datatype = "abundance", knots = 3) {

  
  data <- fread(file_path)
  
  data_wide <- reshape2::dcast(data, Ages ~ KEGG_ko, value.var = "SUM_CPM")

  
  numeric_columns <- sapply(data_wide, is.numeric)

  
  data_wide[numeric_columns] <- lapply(data_wide[numeric_columns], round, digits = 0)

  
  data_transposed <- t(data_wide)

  
  mydf <- as.data.frame(data_transposed)

  
  mydf[is.na(mydf)] <- 0

  
  endpoint <- min(colSums(data_transposed))

  
  result <- iNEXT(mydf[-1,], q = 0, datatype = datatype, endpoint = min(colSums(mydf[-1, 1:42])), knots = knots)

  return(result)
}

# input file
file_paths <- c("prokgap_raw_count_eggnog_abundance_prokaryotes_KOs_abundance.csv",
                "eukgap_raw_count_eggnog_abundance_prokaryotes_KOs_abundance.csv",
                "preclass_prokgap_raw_count_eggnog_abundance_prokaryotes_KOs_abundance.csv",
                "preclass_eukgap_raw_count_eggnog_abundance_prokaryotes_KOs_abundance.csv")
results_list <- lapply(file_paths, process_and_analyze)

# Save results
save(results_list, file = "abundance_analysis_results_prokaryotes.RData")
