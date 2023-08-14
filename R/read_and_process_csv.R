#' Read and Process CSV Files
#'
#' This function reads and processes CSV files from a specified directory that
#' match a given pattern. It extracts information from the filenames and adds
#' it as new columns to the data frames. It then combines all the data frames
#' into a single data frame.
#'
#' @param csv_dir Character string. The directory where the CSV files are located.
#' @param pattern Character string. The pattern to match the CSV files.
#'
#' @return A data frame containing the combined data from the CSV files.
#'
#' @import data.table
#' @import stringr
#' @export
#'

read_and_process_csv <- function(csv_dir, pattern) {

  # Get a list of matching CSV files in the directory
  csv_files <- list.files(csv_dir, pattern = pattern,
                          full.names = TRUE)

  # Create an empty list to store the data frames
  data_list <- list()

  # Read and process each CSV file
  for (file in csv_files) {
    # Read the CSV file
    data <- fread(file, sep = ",", data.table = FALSE, header = TRUE)

    # Check if the data frame has any rows (non-empty)
    if (nrow(data) > 0) {
      # Add a new column with the filename
      data$filename <- basename(file)
      data$projectcode <- str_extract(data$filename, "^([^_]+)")
      data$treatment <- str_extract(data$filename, "(?<=_)[^_]+(?=X)")
      data$trait <- str_extract(data$filename, "(?<=X)[^_]+_.*(?=_gwas)")
      data$day <- str_extract(data$trait, "\\d+")
      data$Gtrait <- str_extract(data$trait, "(?<=\\d_).*")

      # Add the data frame to the list
      data_list <- c(data_list, list(data))
    } else {
      # If the data frame has no rows, show a warning message
      warning(paste("Empty data frame in file:", basename(file)))
    }
  }

  # Check if any valid data frames were read
  if (length(data_list) == 0) {
    stop("No valid data found in the CSV files.")
  }

  # Combine all the data frames into a single data frame
  combined_data <- do.call(rbind, data_list)

  return(combined_data)
}
