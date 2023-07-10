runPredictionAnalysis <- function(data) {
  library(tidyverse)
  library(data.table)

  # Check if the input is a data frame
  if (is.data.frame(data)) {
    sv2 <- data
  } else if (is.character(data)) {
    # Check if the input is a CSV file
    if (tools::file_ext(data) == "csv") {
      sv2 <- fread(data, sep = ",", header = TRUE, data.table = FALSE)
    } else {
      message("Error: Unsupported file format")
      return(NULL)
    }
  } else {
    message("Error: Invalid input")
    return(NULL)
  }

  # Empty lists to store predictions and problematic genotypes
  loess_fits <- list()
  predictions <- list()
  problematic_genotypes <- character()

  # Function to subset data by genotype and fit loess model
  genotype_subset <- function(geno, data) {
    subset_data <- subset(data, genotype == geno)
    subset_data$area.pixels <- as.numeric(subset_data$area.pixels)

    tryCatch({
      loess_fit <- loess(area.pixels ~ DAP, data = subset_data, span = 1)
      predicted_values <- predict(loess_fit, newdata = subset_data)

      return(predicted_values)
    }, error = function(e) {
      problematic_genotypes <<- c(problematic_genotypes, geno)
      return(NULL)
    })
  }

  # Get unique genotypes
  genotypes <- unique(sv2$genotype)

  # Apply the genotype_subset function for each genotype
  predictions <- lapply(genotypes, function(geno) {
    genotype_subset(geno, data = sv2)
  })

  # Print problematic genotypes
  if (length(problematic_genotypes) > 0) {
    cat("The following genotypes caused an error:\n")
    print(problematic_genotypes)

    # Remove problematic genotypes from sv2
    sv2 <- sv2 %>% filter(!genotype %in% problematic_genotypes)

    # Remove corresponding predictions
    predictions <- predictions[!genotypes %in% problematic_genotypes]
  } else {
    cat("No problematic genotypes found.\n")
  }

  # Add predicted_area_pixels as a new column with NA as default value
  sv2 <- sv2 %>% mutate(predicted_area_pixels = NA_real_)

  # Update sv2 with predictions
  for (i in seq_along(predictions)) {
    if (!is.null(predictions[[i]])) {
      sv2 <- sv2 %>% mutate(predicted_area_pixels = if_else(genotype == genotypes[i], predictions[[i]], predicted_area_pixels))
    }
  }

  # Calculate daily area differences and filter out non-shrinking areas
  predictions <- sv2 %>%
    mutate(daily_area_diff_pred = predicted_area_pixels - lag(predicted_area_pixels, default = 0)) %>%
    mutate(daily_area_diff_pred = if_else(daily_area_diff_pred < 0, NA_real_, daily_area_diff_pred)) %>%
    filter(DAP != 0)

  # Check if there are areas not shrinking
  if (sum(is.na(predictions$daily_area_diff_pred)) == 0) {
    print("Areas not shrinking")
  } else {
    cat("The number of cases with shrinking areas is:", sum(is.na(predictions$daily_area_diff_pred)), "\n")
  }

  # Return the final data frame
  return(sv2)
}
