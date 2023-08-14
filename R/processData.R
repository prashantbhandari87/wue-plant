#' Process SV File
#'
#' This function processes an SV file and performs various data manipulations and calculations.
#'
#' @param sv_file Path to the SV file.
#' @param project_code Code for the project.
#'
#' @return A processed data frame containing the calculated values.
#' @export
#'
#' @import pcvr
#' @import data.table
#' @import tidyverse
#' @import lubridate
#' @import modelr
#' @import patchwork
#' @import naniar
#' @importFrom memoise memoise

.processData <- function(sv_file, project_code) {


  tryCatch({
    sv <- read.pcv(sv_file, mode = "wide", singleValueOnly = TRUE, reader = "fread")
    sv <- sv %>%
      rename(plantbarcode = barcode) %>%
      separate(timestamp, into = c("date1", "date2"), sep = "_") %>%
      mutate(across(c(date1, date2), ~ ymd_hms(.))) %>%
      rename(timestamp = date1) %>%
      separate(plantbarcode, into = c("first_part", "second_part"), sep = "_") %>%
      rename(plantbarcode = first_part) %>%
      filter(!grepl("^Cc", plantbarcode)) %>%
      select(-c(date2, second_part))

    sv$genotype <- substr(sv$plantbarcode, 3, 5)
    sv$treatment <- substr(sv$plantbarcode, 6, 7)
    sv$treatment <- ifelse(sv$treatment == "AA", "100",
                           ifelse(sv$treatment == "AB", "30", "0"))

    sv$rotation <- as.integer(ifelse(sv$rotation == "90_90", 90,
                                     ifelse(sv$rotation == "0_0", 0,
                                            ifelse(sv$rotation == "0", 0,
                                                   ifelse(sv$rotation == "90", "90", NA)))))

    sv <- bw.time(sv, plantingDelay = 0, phenotype = "area.pixels", cutoff = 10,
                  timeCol = "timestamp", group = c("plantbarcode", "rotation"), plot = FALSE)


    sv <- bw.outliers(sv, phenotype = "area.pixels", group = c("DAP", "genotype", "treatment"),
                      plotgroup = c("plantbarcode", "rotation"))

    phenotypes <- c('area.pixels', 'convex_hull_area.pixels', 'convex_hull_vertices',
                    'ellipse_angle.degrees', 'ellipse_eccentricity',
                    'ellipse_major_axis.pixels', 'ellipse_minor_axis.pixels',
                    'height.pixels', 'hue_circular_mean.degrees',
                    'hue_circular_std.degrees', 'hue_median.degrees',
                    'longest_path.pixels', 'perimeter.pixels', 'solidity',
                    'width.pixels')

    phenoForm <- paste0("cbind(", paste0(phenotypes, collapse = ", "), ")")
    groupForm <- "DAP + plantbarcode + genotype + treatment"
    form <- as.formula(paste0(phenoForm, "~", groupForm))

    sv <- aggregate(form, data = sv, mean, na.rm = TRUE)


    return(sv)
  }, error = function(e) {
    message("An error occurred: ", conditionMessage(e))
    return(NULL)
  })
}

processData <- memoise::memoise(.processData)
