library(tidyverse)
library(glue)

setwd("C:/WFSETP/TTN_paper/FSim_Outputs")
#setwd("C:/WFSETP/TTN_paper/test_cmdx")
wd <- getwd()

#Read in the template cmdx
template <- read_file("C:/WFSETP/R_scripts/cmdx_template.txt")

#Create dataframe of scenarios
#list the names of the run directories
run_dirs <- list.dirs(wd, recursive = FALSE, full.names = TRUE)

#Function to extract input files for one area
get_input_files <- function(foa_dir){
  inputs_dir <- file.path(foa_dir, "_inputs")
  
  filename_parts <- str_split(basename(foa_dir), "_", simplify = TRUE) %>% as.vector()
  
  study_area <- filename_parts[1]
  timepoint  <- filename_parts[length(filename_parts)]
  foa_run    <- paste(filename_parts[2:3], collapse = "_")   
  scenario   <- paste(filename_parts[4:(length(filename_parts)-1)], collapse = "_")
  
  subdirs <- c("adj", "erc", "fdist", "frisk", "idg", "lcp")
  
  erc_files <- list.files(file.path(inputs_dir, "erc"), full.names = FALSE)
  
  map_dfr(erc_files, function(erc_file){
    # extract part number from filename
    part_num <- str_match(erc_file, "ERC(\\d+)\\.csv$")[,2]
    
    base_tbl <- tibble(
      study_area, foa_run, scenario, timepoint,
      part = part_num,
      !!!set_names(
        map(subdirs, ~ {
          files <- list.files(file.path(inputs_dir, .x), full.names = FALSE)
          
          if(.x == "erc") {
            # use the current erc file for this row
            files <- erc_file
          }
          if(.x == "lcp"){
            # only keep the .tif, not the .tif.aux.xml
            files <- files[str_detect(files, "(?i)\\.tif$")]
          }
          if(length(files) == 0) NA_character_ else files[1]
        }),
        subdirs
      )
    )
    
    # handle fms separately
    fms_files <- list.files(file.path(inputs_dir, "fms"), full.names = FALSE)
    
    fms_tbl <- tibble(
      fms80 = fms_files[str_detect(fms_files, "80")],
      fms90 = fms_files[str_detect(fms_files, "90")],
      fms97 = fms_files[str_detect(fms_files, "97")]
    )
    
    bind_cols(base_tbl, fms_tbl)
  })
}

#Use this if you want to automate the process over several run directories
scenario_df <- map_dfr(run_dirs, get_input_files)
view(scenario_df)

for (i in seq_len(nrow(scenario_df))) {
  
  row <- scenario_df[i, ]
  
  output_dir <- file.path(
    wd,
    paste0(row$study_area, "_", row$foa_run, "_", row$scenario, "_", row$timepoint)
  )
  
  # fill placeholders
  contents <- glue(template, .envir = list(
    study_area = row$study_area,
    foa_run = row$foa_run,
    scenario = row$scenario,
    timepoint = row$timepoint,
    part = row$part,
    idg = row$idg,
    adj = row$adj,
    lcp = row$lcp,
    fdist = row$fdist,
    frisk = row$frisk,
    erc = row$erc
  ))
  
  # output filename
  cmdx_file <- file.path(
    output_dir,
    paste0(row$foa_run, "_pt", row$part, "_", row$scenario, "_", row$timepoint, ".cmdx")
  )
  
    writeLines(contents, cmdx_file)
  
  message("Generated ", row$part, " cmdx files for ",
          row$study_area, " ", row$foa_run, ", scenario ", row$scenario, ", run timepoint ", row$timepoint, ".")
}
