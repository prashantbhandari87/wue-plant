#' Process Arrow files and create a summary CSV file
#'
#' This function reads Arrow files, filters data based on a p-value threshold,
#' and creates a summary CSV file.
#'
#' @param csv_dir The directory containing the Arrow files.
#' @param pattern The pattern to match Arrow files.
#' @param outname The name of the output CSV file.
#' @param pvalue_threshold The threshold for filtering data based on p-value.
#'
#' @return NULL if there is no data matching the specified criteria, otherwise,
#' it returns the name of the output CSV file.
#'
#' @import arrow
#' @import tidyverse
#' @import parallel
#' @import doParallel
#'
#' @examples
#' \dontrun{
#' process_arrow_files_for_summary(pattern = ".arrow", outname = "summary.csv", csv_dir = "/path/to/your/directory", pvalue_threshold = 4)
#' }
#'
#' @export
process_arrow_files_for_summary <- function(csv_dir = getwd(), pattern = ".arrow", outname = "s1_terra_height.pixels.csv", pvalue_threshold = 4) {
  tryCatch({
    setwd(csv_dir)  # Set working directory
    
    library(arrow)
    library(tidyverse)
    library(parallel)
    library(doParallel)
    
    arrow_files <- list.files(pattern = pattern, full.names = TRUE)
    
    if (length(arrow_files) == 0) {
      stop("No files found matching the specified pattern.")
    }
    
    print("Pattern:"); print(pattern)
    print("Reading:"); print(arrow_files)
    
    data_list <- list()
    for (file in arrow_files) {
      dataset <- arrow::open_dataset(file)
      dataset <- dataset %>%
        filter(-pvalue > pvalue_threshold) %>%
        as.data.frame() %>%
        mutate(snp = id) %>%
        select(-c(id, score)) %>%
        arrange(POS, CHR, estim, std.err, trait, snp) %>%
        collect()
      data_list <- c(data_list, list(dataset))
    }
    
    if (length(data_list) == 0) {
      warning("No data matching the specified criteria.")
      return(NULL)
    }
    
    random_subset <- data_list %>% bind_rows() %>% arrange(snp)
    write.csv(random_subset, outname, row.names = FALSE)
    cat("Process completed successfully. Output written to:", outname, "\n")
    return(outname)
  }, error = function(e) {
    cat("An error occurred:", conditionMessage(e), "\n")
  })
}

