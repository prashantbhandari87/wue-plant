#' Merge pcv_data and water_data based on plantbarcode and day
#'
#' @param pcv_data Data frame containing pcv data
#' @param water_data Data frame containing water data
#' @param project_code Code for the project
#'
#' @return Merged data frame with additional calculated columns
#'
#' @export

merge_pcv_water <- function(pcv_data,water_data,project_code){

  overlapping_columns <- intersect(names(pcv_data), names(water_data))
  cat("Overlapping columns:", paste(overlapping_columns, collapse = ", "), "\n")

  pcv_data <- pcv_data %>% mutate(day = 12 + DAP) %>%
    select(-c(genotype, treatment, DAP))

  joined1 <- left_join(water_data, pcv_data, by = c("plantbarcode",
                                               "day"))

  output_file <- paste(project_code, "large_plot_file.csv",
                       sep = "_")
  write.csv(joined1, file = output_file, row.names = FALSE)


  joined1 <- joined1 %>%
    group_by(genotype,treatment) %>%
    filter(!is.na(area.pixels) & !is.na(water.amount.plus.dap)) %>%
    do({
      mod_area = loess(area.pixels ~ day, data = ., span = 0.75, degree = 2)
      mod_water = loess(water.amount.plus.dap ~ day, data = ., span = 0.75, degree = 2)

      residuals_area = resid(mod_area)
      residuals_water = resid(mod_water)

      fits_area = predict(mod_area)
      fits_water = predict(mod_water)

      data.frame(
        .,
        residuals_area = residuals_area,
        residuals_water = residuals_water,
        fits_area = fits_area,
        fits_water = fits_water
      )
    }) %>%
    group_by(plantbarcode, day) %>%
    #select(genotype, treatment, day, replicate,area.pixels, water.amount.plus.dap,residuals_area,
   #        residuals_water,fits_area,fits_water) %>%
    distinct() %>%
    group_by(genotype,treatment) %>%
    arrange(genotype,treatment, day) %>%
    mutate(daily_area_diff_pred = fits_area - lag(fits_area, default = 0)) %>%
    mutate(daily_area_diff_pred = if_else(daily_area_diff_pred < 0, NA_real_, daily_area_diff_pred)) %>%
    mutate(dailyWUE = daily_area_diff_pred/fits_water)


  wue.fit.models <- joined1 %>% group_by(day,treatment) %>% do(mod = lm(daily_area_diff_pred ~
                                                             fits_water, data = .))

  m <- left_join(joined1, wue.fit.models, by = c("day","treatment"))
  x <- m %>% group_by(day,treatment) %>%
    do(add_predictions(., first(.$mod), var = "WUE.fit")) %>%
    do(add_residuals(., first(.$mod),
                     var = "WUE.resid")) %>% select(-c(mod))

x <- x %>% filter(daily_area_diff_pred>0)

  output_file1 <- paste(project_code, "smooth.csv",
                       sep = "_")

  write.csv(x, file = output_file1, row.names = FALSE)


  return(x)

}
