#' Process Traits
#'
#' This function processes traits data provided as a data frame or a CSV file.
#'
#' @param project_code Project code for naming the output files.
#' @param x Traits data, either a data frame or a CSV file.
#'
#' @return Processed data frames with traits analysis results.
#' @export
#'
#' @import tidyverse
#' @import data.table
#'
#' @examples
#' processTraits(project_code = "project001", x = my_traits_data)
#' processTraits(project_code = "project001", x = "path/to/traits_data.csv")
processTraits <- function(project_code, x) {
  library(tidyverse)
  library(data.table)

  # Check if the input is a data frame
  if (is.data.frame(x)) {
    x_df <- x
  } else if (is.character(x)) {
    # Check if the input is a CSV file
    if (tools::file_ext(x) == "csv") {
      x_df <- fread(x, data.table = FALSE, sep = ",", header = TRUE)
    } else {
      message("Error: Unsupported file format")
      return(NULL)
    }
  } else {
    message("Error: Invalid input")
    return(NULL)
  }

  # Remove day 12 from data frame x
  x_df <- x_df %>% filter(day != 12)

  # Remove plantbarcode column from data frame x
  x_df <- x_df %>% select(-c(plantbarcode))

  # Cast the data frame x using dcast function
  castedLT <- dcast(melt(x_df, id.vars = c("genotype", "treatment", "day", "replicate")),
                    ... ~ day + variable)

  # Filter castedLT for treatment 30 and treatment 100
  castedLT_30 <- castedLT %>% filter(treatment == 30)
  castedLT_100 <- castedLT %>% filter(treatment == 100)

  # Write the results to CSV files
  write.csv(castedLT, file = paste(project_code, "preblup_traits", sep = "_"), row.names = FALSE)
  write.csv(castedLT_30, file = paste(project_code, "T30_preblup_traits", sep = "_"), row.names = FALSE)
  write.csv(castedLT_100, file = paste(project_code, "T100_preblup_traits", sep = "_"), row.names = FALSE)

  # Calculate median values
  median.x <- aggregate(. ~ genotype + day + treatment, x_df, median, na.rm = TRUE)
  median.x <- median.x %>% select(-c(replicate))

  # Cast the median.x data frame
  median.x.LT <- dcast(melt(median.x, id.vars = c("genotype", "treatment", "day")),
                       ... ~ day + variable)

  # Filter median.x.LT for treatment 30 and treatment 100
  median.x.LT.30 <- median.x.LT %>% filter(treatment == 30)
  median.x.LT.100 <- median.x.LT %>% filter(treatment == 100)

  # Write the results to CSV files
  write.csv(median.x.LT, file = paste(project_code, "median_traits", sep = "_"), row.names = FALSE)
  write.csv(median.x.LT.30, file = paste(project_code, "T30_median_traits", sep = "_"), row.names = FALSE)
  write.csv(median.x.LT.100, file = paste(project_code, "T100_median_traits", sep = "_"), row.names = FALSE)
}
