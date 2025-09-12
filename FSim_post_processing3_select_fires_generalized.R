library(tidyverse)
library(furrr)
library(terra)
library(tidyterra)

#######################################################################################
# NOTE: To run this code, you need to make sure the following FSim outputs are in the #
#       working directory: the FireSizeList.csv files, the All.tif files, the         #
#       ArrivalDays tifs, the FlameLengths tifs, and the Perimeters.sqlite files.     #
#######################################################################################

#First, we need to crop the ignitions to our 15 km buffer around okawen and filter the merged firesizelists by those fireIDs
okawen_15km_buf <- vect("../Shapefiles/okawen_only_buffer_15km/okawen_only_Buffer_15km.shp")
ex_ig <- vect("C:/WFSETP/TTN_paper/FSim_Outputs/okawen_foa1c_r16_LF2022_recoff3_time2/ignitions_foa1c_r16_LF2022_recoff3_time2/ignitions_foa1c_r16_LF2022_recoff3_time2.shp")
ex_firelists <- read_csv("C:/WFSETP/TTN_paper/FSim_Outputs/okawen_foa1c_r16_LF2022_recoff3_time2/foa1c_r16_LF2022_recoff3_time2_merged_firelists.csv")

ex_ig <- ex_ig %>%
  rename(Season_FireID = Season_Fir)
head(ex_ig$Season_FireID)
head(ex_firelists$Season_FireID)

summarize_and_sample_fires <- function(study_area,
                                       foas, runs, scenarios, timepoints,
                                       n_seasons = 20000,
                                       n_samples = 100,
                                       k_years = 5,
                                       out_dir = "C:/WFSETP/TTN_paper/FSim_Outputs",
                                       workers = parallel::detectCores() - 1,
                                       seed = 42,
                                       save_per_foa = FALSE) {


  params_df <- tibble(foa = foas, run = runs) %>%
    crossing(tibble(scenario = scenarios, timepoint = timepoints))

  
  # Setup parallel plan for reading
  plan(multisession, workers = workers)
  
  # Read all merged firelists (one per param row). Each read returns a tibble with FOA/run metadata.
  # We try to be robust if files are missing.
  all_firelists_list <- future_pmap(params_df, function(foa, run, scenario, timepoint) {
    folder <- file.path(out_dir, paste0(study_area, "_", run, "_", scenario, "_", timepoint))
    file   <- file.path(folder, paste0(run, "_", scenario, "_", timepoint, "_merged_firelists.csv"))
    
    # Read firelist
    firelist <- read_csv(file, show_col_types = FALSE)
    
    # Ignitions shapefile path (assumes naming convention matches folder)
    ignitions_file <- file.path(folder, paste0("ignitions_", run, "_", scenario, "_", timepoint, "/", "ignitions_", run, "_", scenario, "_", timepoint, ".shp"))
    
    if (!file.exists(ignitions_file)) {
      warning("Ignitions file missing for ", ignitions_file, " — skipping filtering")
      firelist_filt <- firelist
    } else {
      ignitions <- vect(ignitions_file)
      okawen_15km_buf <- vect("C:/WFSETP/TTN_paper/Shapefiles/okawen_only_buffer_15km/okawen_only_Buffer_15km.shp")
      
      # Crop ignitions to 15 km buffer
      ignitions_crop <- ignitions[okawen_15km_buf, ]
      print(head(ignitions_crop))
      
            # Filter firelist based on ignitions inside buffer
      firelist_filt <- firelist %>%
        dplyr::filter(Season_FireID %in% ignitions_crop$Season_Fir)
    }
    
    # Tag with FOA/run/scenario/timepoint info
    return(firelist_filt %>%
      mutate(FOA = foa, run = run, scenario = scenario, timepoint = timepoint))
  }, .options = furrr_options(seed = NULL)) # reading doesn't need RNG
  
  # Combine into a single table (this is the full study-area firelist across FOAs)
  all_firelists <- bind_rows(all_firelists_list)
  
  if (nrow(all_firelists) == 0) {
    stop("No firelist rows were read. Check file paths and input arguments.")
  }
  
  total_seasons <- seq_len(n_seasons)
  
  # Create sampling sets on the main process for reproducibility
  set.seed(seed)
  sampled_sets <- replicate(n_samples, sample(total_seasons, k_years, replace = FALSE), simplify = FALSE)
  
  # For each sampled 5-year set, compute total area burned across ALL FOAs (use all_firelists)
  # This step can be parallelized if n_samples large; use future_map_dbl for safety
  plan(multisession, workers = max(1, min(workers, n_samples))) # switch plan to map over sets
  total_area_burned <- future_map_dbl(sampled_sets, function(season_set) {
    all_firelists %>%
      filter(Season %in% season_set) %>%
      summarise(total_area = sum(Acres, na.rm = TRUE)) %>%
      pull(total_area)
  }, .options = furrr_options(seed = NULL))
  
  # Also create long version with one row per season
  distribution_long_df <- tibble(
    sample_id = rep(seq_along(total_area_burned), each = k_years),
    Season = unlist(sampled_sets),
    total_area_burned = rep(total_area_burned, each = k_years)
  )
  hist(distribution_long_df$total_area_burned)
  
  write_csv(distribution_long_df, file.path(out_dir, paste0(study_area, "_five_year_set_distribution_long.csv")))
  
  # Find the 20th, 60th, 90th percentile values and the closest sampled sets
  distribution_df <- tibble(sample_id = seq_along(total_area_burned),
                            total_area_burned = total_area_burned)
  targets <- quantile(distribution_df$total_area_burned, probs = c(0.2, 0.6, 0.9), na.rm = TRUE, names = FALSE)
  idxs <- sapply(targets, function(t) which.min(abs(distribution_df$total_area_burned - t)))
  names(idxs) <- c("p20_idx", "p60_idx", "p90_idx")
  
  decile_list <- list(
    p20 = list(sample_id = idxs["p20_idx"],
               total_area = distribution_df$total_area_burned[idxs["p20_idx"]],
               seasons = sampled_sets[[idxs["p20_idx"]]]),
    p60 = list(sample_id = idxs["p60_idx"],
               total_area = distribution_df$total_area_burned[idxs["p60_idx"]],
               seasons = sampled_sets[[idxs["p60_idx"]]]),
    p90 = list(sample_id = idxs["p90_idx"],
               total_area = distribution_df$total_area_burned[idxs["p90_idx"]],
               seasons = sampled_sets[[idxs["p90_idx"]]])
  )
  

  # Long version for deciles
  decile_long_df <- tibble(
    decile = rep(c("20th", "60th", "90th"), each = k_years),
    sample_id = rep(c(decile_list$p20$sample_id, decile_list$p60$sample_id, decile_list$p90$sample_id), each = k_years),
    Season = c(decile_list$p20$seasons, decile_list$p60$seasons, decile_list$p90$seasons),
    total_area = rep(c(decile_list$p20$total_area, decile_list$p60$total_area, decile_list$p90$total_area), each = k_years)
  )
  
  write_csv(decile_long_df, file.path(out_dir, paste0(study_area, "_five_year_set_deciles_long.csv")))
  
  # Optionally save per-FOA files that contain the fire rows for each chosen decile set
  if (isTRUE(save_per_foa)) {
    # For every unique FOA / run / scenario / timepoint combination present in the read-in lists,
    # save a CSV containing rows belonging to each decile set, labeled with decile.
    per_params <- all_firelists %>%
      select(FOA, run, scenario, timepoint) %>%
      distinct()
    
    for (i in seq_len(nrow(per_params))) {
      this_foa <- per_params$FOA[i]
      this_run <- per_params$run[i]
      this_scn <- per_params$scenario[i]
      this_tp  <- per_params$timepoint[i]
      
      subset_df <- all_firelists %>%
        filter(FOA == this_foa, run == this_run, scenario == this_scn, timepoint == this_tp)
      
      # for each decile label, grab rows matching the seasons
      df_p20 <- subset_df %>% filter(Season %in% decile_list$p20$seasons) %>% mutate(decile = "20th", sample_id = decile_list$p20$sample_id)
      df_p60 <- subset_df %>% filter(Season %in% decile_list$p60$seasons) %>% mutate(decile = "60th", sample_id = decile_list$p60$sample_id)
      df_p90 <- subset_df %>% filter(Season %in% decile_list$p90$seasons) %>% mutate(decile = "90th", sample_id = decile_list$p90$sample_id)
      
      out_df <- bind_rows(df_p20, df_p60, df_p90) %>%
        arrange(decile, Season)
      
      # output file name
      out_folder <- file.path(out_dir, paste0(study_area, "_", this_run, "_", this_scn, "_", this_tp))
      if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)
      out_file <- file.path(out_folder, paste0(this_run, "_", this_scn, "_", this_tp, "_selected_5yrsets_deciles.csv"))
      
      write_csv(out_df, out_file)
      message("Wrote per-FOA decile file: ", out_file)
    } # end per-FOA loop
  }
  
  # Return a list invisibly for programmatic inspection
  invisible(list(
    distribution = distribution_long_df,
    deciles = decile_long_df,
    decile_sets = decile_list,
    sampled_sets = sampled_sets
  ))
}

study_area <- "okawen"
foas <- c("foa1", "foa2", "foa3")
runs <- c("foa1c_r16", "foa2d_r5", "foa3d_r8")
scenarios <- c("LF2022_recoff3")
timepoints <- c("time2")

summarize_and_sample_fires(study_area, foas, runs, scenarios, timepoints,
                           n_seasons = 20000, n_samples = 100, k_years = 5,
                           out_dir = "C:/WFSETP/TTN_paper/FSim_Outputs",
                           workers = 28,
                           seed = 42,
                           save_per_foa = FALSE)


