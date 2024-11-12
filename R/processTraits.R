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

processTraits <- function(project_code, x) {


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
  x_df <- x_df %>% filter(day > 12)

  # Calculate median values
  median.x <- x_df %>%
    group_by(genotype, day, treatment) %>%
    summarise(across(where(is.numeric), median))


  median.x <- median.x %>% select(-c(replicate))

  # Cast the median.x data frame
  median.x.LT <- dcast(melt(median.x, id.vars = c("genotype", "treatment", "day")),
                       ... ~ day + variable)

  # Filter median.x.LT for treatment 30 and treatment 100
  median.x.LT.30 <- median.x.LT %>% filter(treatment == 30)
  median.x.LT.100 <- median.x.LT %>% filter(treatment == 100)

  # Write the results to CSV files
  write.csv(median.x.LT, file = paste(project_code, "median_traits.csv", sep = "_"), row.names = FALSE)
  write.csv(median.x.LT.30, file = paste(project_code, "T30_median_traits.csv", sep = "_"), row.names = FALSE)
  write.csv(median.x.LT.100, file = paste(project_code, "T100_median_traits.csv", sep = "_"), row.names = FALSE)
  print("Files written out!")

return(median.x.LT)
}
