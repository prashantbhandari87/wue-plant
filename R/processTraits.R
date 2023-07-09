processTraits <- function(project_code, x) {
  library(tidyverse)

  # Remove day 12 from data frame x
  x <- x %>% filter(day != 12)

  # Remove plantbarcode column from data frame x
  x <- x %>% select(-c(plantbarcode))

  # Cast the data frame x using dcast function
  castedLT <- dcast(melt(x, id.vars = c("genotype", "treatment", "day", "replicate")),
                    ... ~ day + variable)

  # Filter castedLT for treatment 30 and treatment 100
  castedLT_30 <- castedLT %>% filter(treatment == 30)
  castedLT_100 <- castedLT %>% filter(treatment == 100)

  # Write the results to CSV files
  write.csv(castedLT, file = paste(project_code, "preblup_traits", sep = "_"), row.names = FALSE)
  write.csv(castedLT_30, file = paste(project_code, "T30_preblup_traits", sep = "_"), row.names = FALSE)
  write.csv(castedLT_100, file = paste(project_code, "T100_preblup_traits", sep = "_"), row.names = FALSE)

  # Calculate median values
  median.x <- aggregate(. ~ genotype + day + treatment, x, median, na.rm = TRUE)
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
