#' Analyze Water Data
#'
#' This function analyzes water data provided as a data frame or a CSV file.
#'
#' @param water_data Water data, either a data frame or a CSV file.
#' @param project_code Project code for naming the output file.
#'
#' @return Processed data frame with water analysis results.
#' @export
#'
#' @import tidyverse
#' @import data.table
#'

analyzeWaterData <- function(water_data, predicted_pcv_data, project_code) {
  library(tidyverse)
  library(data.table)

  # Check if the input is a data frame
  if (is.data.frame(water_data)) {
    water_df <- water_data
  } else if (is.character(water_data)) {
    # Check if the input is a CSV file
    if (tools::file_ext(water_data) == "csv") {
      water_df <- fread(water_data, header = TRUE,sep=",",data.table=FALSE)
      water_df$genotype <- substr(water_df$plantbarcode, 3,5)
      water_df$treatment <- substr(water_df$plantbarcode, 6, 7)
      water_df$treatment <- ifelse(water_df$treatment == "AA", "100",
                                   ifelse(water_df$treatment == "AB", "30", "0"))
    } else {
      message("Error: Unsupported file format")
      return(NULL)
    }
  } else {
    message("Error: Invalid input")
    return(NULL)
  }

  # Check for outliers in weight.before
  check_outlier <- water_df %>%
    filter(weight.before < 0)

  if (dim(check_outlier)[1] == 0) {
    print("No outliers in weight.before.")
  } else {
    print("Outliers in weight.before removed.")
  }

  # Check for outliers in weight.after
  check_outlier <- water_df %>%
    filter(weight.after < 0)

  if (dim(check_outlier)[1] == 0) {
    print("No outliers in weight.after.")
  } else {
    print("Outliers in weight.after removed.")
  }

  # Calculate water amount lost
  water_df <- water_df %>% filter(weight.before > -1) %>%
    filter(weight.after > -1) %>%
    group_by(plantbarcode) %>%
    arrange(plantbarcode, day, timestamp) %>%
    mutate(weight.after.lag1 = lag(weight.after,1)) %>%
    mutate(water.amount.plus = weight.after.lag1 -  weight.before) %>%
    filter(!is.na(water.amount.plus)) %>%
    group_by(plantbarcode) %>%
    mutate(time.water.loss = as_datetime(timestamp) -
             as_datetime(lag(timestamp, 1))) %>%
    group_by(plantbarcode, day) %>% arrange(timestamp) %>%
    mutate(watering.job.num = 1:n()) %>%
    arrange(plantbarcode, day, watering.job.num) %>%
    group_by(plantbarcode, day, genotype, treatment, replicate) %>%
    summarize(water.amount.plus.dap = sum(water.amount.plus), .groups = "drop")

  water.join <- water_df %>% select(plantbarcode, genotype, day,
                                    treatment, replicate, water.amount.plus.dap)
  sv2 <- predicted_pcv_data %>% mutate(day = 12 + DAP) %>%
    select(-c(genotype, treatment, DAP))
  joined1 <- left_join(water.join, sv2, by = c("plantbarcode",
                                               "day"))
  joined1 <- distinct(joined1)
  joined1 <- joined1 %>% filter(day != 12)
  joined1 <- joined1 %>% filter(water.amount.plus.dap > 0)
  joined1 <- drop_na(joined1)
  joined1$dailyWUE <- joined1$daily_area_diff_pred/joined1$water.amount.plus.dap
  num_na <- sum(is.na(joined1$dailyWUE))
  if (num_na > 0) {
    cat("There are", num_na, "NA values in joined1$dailyWUE.\n")
  }
  else {
    print("No NA values in joined1$dailyWUE.")
  }
  wue.fit.models <- joined1 %>% group_by(day) %>% do(mod = lm(daily_area_diff_pred ~
                                                                water.amount.plus.dap, data = .))
  m <- left_join(joined1, wue.fit.models, by = c("day"))
  x <- m %>% group_by(day) %>%
    do(add_predictions(., first(.$mod), var = "WUE.fit")) %>%
    do(add_residuals(., first(.$mod), var = "WUE.resid")) %>%
    select(-c(mod))
  output_file <- paste(project_code, "large_plot_file.csv",
                       sep = "_")
  write.csv(x, file = output_file, row.names = FALSE)
  return(x)
}
