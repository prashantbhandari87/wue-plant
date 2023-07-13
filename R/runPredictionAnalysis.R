#' Run Prediction Analysis
#'
#' This function performs prediction analysis on the provided data using the loess model.
#'
#' @param data A dataframe or the path to a CSV file containing the data.
#' @return The processed dataframe after prediction analysis.
#' @export
runPredictionAnalysis <- function(data) {
  library(tidyverse)
  ########
  # Load data if a CSV file path is provided
  if (is.character(data) && tolower(file_ext(data)) == ".csv") {
    if (!file.exists(data)) {
      stop("The specified CSV file does not exist.")
    }
    data <- read.csv(data)
  }

  # Convert 'area.pixels' to numeric
  if (!"area.pixels" %in% colnames(data)) {
    stop("The 'area.pixels' column does not exist in the provided data.")
  }
  data$area.pixels <- as.numeric(data$area.pixels)

  # Filter data for treatment == 30
  sv2.30 <- data %>% filter(treatment == 30)

  # Filter data for treatment == 100
  sv2.100 <- data %>% filter(treatment == 100)

  # Get unique genotypes
  genotypes <- unique(data$genotype)

  # List to store loess fits
  loess_fits <- list()

  # List to store predictions
  predictions <- list()

  # Subset the data by genotype and fit loess models
  genotype_subset <- function(geno, data) {
    subset_data <- subset(data, genotype == geno)
    subset_data$area.pixels <- as.numeric(subset_data$area.pixels)

    loess_fit <- loess(area.pixels ~ DAP, data = subset_data, span = 1)
    predicted_values <- predict(loess_fit, newdata = subset_data)

    return(predicted_values)
  }

  # Empty list to store problematic genotypes
  problematic_genotypes <- c()

  # Apply the function for each genotype using lapply for treatment == 30
  predictions.30 <- lapply(genotypes, function(geno) {
    tryCatch(
      {
        genotype_subset(geno, data = sv2.30)
      },
      error = function(e) {
        problematic_genotypes <<- c(problematic_genotypes, geno)
        NULL
      }
    )
  })

  # Apply the function for each genotype using lapply for treatment == 100
  predictions.100 <- lapply(genotypes, function(geno) {
    tryCatch(
      {
        genotype_subset(geno, data = sv2.100)
      },
      error = function(e) {
        problematic_genotypes <<- c(problematic_genotypes, geno)
        NULL
      }
    )
  })

  # Print the problematic genotypes
  if (length(problematic_genotypes) > 0) {
    cat("The following genotypes caused an error:\n")
    print(problematic_genotypes)
  } else {
    cat("No problematic genotypes found.\n")
  }

  data <- data %>% filter(!genotype %in% problematic_genotypes)
  # Convert 'area.pixels' to numeric
  if (!"area.pixels" %in% colnames(data)) {
    stop("The 'area.pixels' column does not exist in the provided data.")
  }
  data$area.pixels <- as.numeric(data$area.pixels)

  # Filter data for treatment == 30
  sv2.30 <- data %>% filter(treatment == 30)

  # Filter data for treatment == 100
  sv2.100 <- data %>% filter(treatment == 100)

  # Get unique genotypes
  genotypes <- unique(data$genotype)

  # List to store loess fits
  loess_fits <- list()

  # List to store predictions
  predictions <- list()

  # Subset the data by genotype and fit loess models
  genotype_subset <- function(geno, data) {
    subset_data <- subset(data, genotype == geno)
    subset_data$area.pixels <- as.numeric(subset_data$area.pixels)

    loess_fit <- loess(area.pixels ~ DAP, data = subset_data, span = 1)
    predicted_values <- predict(loess_fit, newdata = subset_data)

    return(predicted_values)
  }

  # Empty list to store problematic genotypes
  problematic_genotypes <- c()

  # Apply the function for each genotype using lapply for treatment == 30
  predictions.30 <- lapply(genotypes, function(geno) {
    tryCatch(
      {
        genotype_subset(geno, data = sv2.30)
      },
      error = function(e) {
        problematic_genotypes <<- c(problematic_genotypes, geno)
        NULL
      }
    )
  })

  # Apply the function for each genotype using lapply for treatment == 100
  predictions.100 <- lapply(genotypes, function(geno) {
    tryCatch(
      {
        genotype_subset(geno, data = sv2.100)
      },
      error = function(e) {
        problematic_genotypes <<- c(problematic_genotypes, geno)
        NULL
      }
    )
  })

  # Print the problematic genotypes
  if (length(problematic_genotypes) > 0) {
    cat("The following genotypes caused an error:\n")
    print(problematic_genotypes)
  } else {
    cat("No problematic genotypes found.\n")
  }


  # Create dataframes from predictions and set column names dynamically
  # Combine data frames using bind_cols() instead of cbind()
  predictions.30 <- bind_cols(sv2.30, as.data.frame(unlist(predictions.30)))
  colnames(predictions.30)[ncol(predictions.30)] <- "predicted_area_pixels"

  # Combine data frames using bind_cols() instead of cbind()
  predictions.100 <- bind_cols(sv2.100, as.data.frame(unlist(predictions.100)))

  colnames(predictions.100)[ncol(predictions.100)] <- "predicted_area_pixels"

  # Filter out problematic genotypes from sv2.30 and sv2.100
  predictions.30 <- predictions.30 %>% filter(!genotype %in% problematic_genotypes)
  predictions.100 <- predictions.100 %>% filter(!genotype %in% problematic_genotypes)

  # Calculate daily area difference for predictions.30
  predictions.30 <- predictions.30 %>%
    mutate(daily_area_diff_pred = predicted_area_pixels - lag(predicted_area_pixels, default = 0)) %>%
    mutate(daily_area_diff_pred = if_else(daily_area_diff_pred < 0, NA_real_, daily_area_diff_pred)) %>%
    mutate(out=predicted_area_pixels-
             lag(predicted_area_pixels,default=0)) %>%
    mutate(daily_area_diff_pred = if_else(out<0,NA_real_,daily_area_diff_pred)) %>%
    select(-c(out))

  # Check if areas are shrinking in predictions.30
  if (predictions.30$daily_area_diff_pred %>% is.na() %>% sum() == 0) {
    print("Areas not shrinking")
  } else {
    cat("The number of cases with shrinking areas is:\n")
    print(predictions.30$daily_area_diff_pred %>% is.na() %>% sum())
  }

  # Calculate daily area difference for predictions.100
  predictions.100 <- predictions.100 %>%
    mutate(daily_area_diff_pred = predicted_area_pixels - lag(predicted_area_pixels, default = 0)) %>%
    mutate(daily_area_diff_pred = if_else(daily_area_diff_pred < 0, NA_real_, daily_area_diff_pred)) %>%
    mutate(out=predicted_area_pixels-
             lag(predicted_area_pixels,default=0)) %>%
    mutate(daily_area_diff_pred = if_else(out<0,NA_real_,daily_area_diff_pred)) %>%
    select(-c(out))

  # Check if areas are shrinking in predictions.100
  if (predictions.100$daily_area_diff_pred %>% is.na() %>% sum() == 0) {
    print("Areas not shrinking")
  } else {
    cat("The number of cases with shrinking areas is:\n")
    print(predictions.100$daily_area_diff_pred %>% is.na() %>% sum())
  }

  # Combine predictions for treatments 30 and 100
  data <- as.data.frame(rbind(predictions.100, predictions.30))
  data <- drop_na(data)

  # Return the final processed dataframe
  return(data)

}

