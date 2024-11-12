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

clean_water <- function (water_data, project_code) {

  if (is.data.frame(water_data)) {
    water_df <- water_data
  }
  else if (is.character(water_data)) {
    if (tools::file_ext(water_data) == "csv") {
      water_df <- fread(water_data, header = TRUE, sep = ",",
                        data.table = FALSE)
      water_df$treatment <- substr(water_df$plantbarcode,
                                   6, 7)
      water_df$treatment <- ifelse(water_df$treatment ==
                                     "AA", "100", ifelse(water_df$treatment == "AB",
                                                         "30", "0"))
    }
    else {
      message("Error: Unsupported file format")
      return(NULL)
    }
  }
  else {
    message("Error: Invalid input")
    return(NULL)
  }
  check_outlier <- water_df %>% filter(weight.before < 0)
  if (dim(check_outlier)[1] == 0) {
    print("No outliers in weight.before.")
  }
  else {
    print("Outliers in weight.before removed.")
  }
  check_outlier <- water_df %>% filter(weight.after < 0)
  if (dim(check_outlier)[1] == 0) {
    print("No outliers in weight.after.")
  }
  else {
    print("Outliers in weight.after removed.")
  }
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


  return(water.join)
}

