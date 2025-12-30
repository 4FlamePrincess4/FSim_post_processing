library(tidyverse)

#Set your working directory to wherever you want the split files to be saved
setwd("C:/Users/laure/OneDrive/Documents/WFSETP/Data/National_Run/National_FSim_Inputs/FOA_1/")
wd <- getwd()

#Assuming you have all the ERC files in subdirectories of a parent directory, use recursive=TRUE
input_dir <- "C:/Users/laure/OneDrive/Documents/WFSETP/Data/National_Run/National_FSim_Inputs/"
season_erc_total <- list.files(input_dir, pattern = "FOA.*_SeasonERC\\.csv$", full.names = TRUE,
                               recursive=TRUE)

for (file in season_erc_total) {
  df <- read_csv(file, show_col_types = FALSE)
  
  #Get the base filename to then extract the FOA label
  base_name <- str_extract(basename(file), "FOA\\d+")
  
  #Split into 8 equal parts of 2500 rows
  split_dfs <- df %>%
    mutate(part = ceiling(row_number() / 2500)) %>%
    group_split(part)%>%
    lapply(select, -part)
  
  #Save each split
  for (i in seq_along(split_dfs)) {
    out_file <- file.path(wd, paste0(base_name, "_pt", i, "_SeasonERC.csv"))
    write_csv(split_dfs[[i]], out_file)
  }
}

