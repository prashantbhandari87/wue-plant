#' Process Data
#'
#' This function performs the processing steps on the provided SV data.
#'
#' @param sv_file The file path of the SV data.
#' @param project_code The project code.
#' @return Processed SV data.
#'
process_data <- function(sv_file, project_code) {
  library(pcvr)
  library(data.table)
  library(tidyverse)
  library(lubridate)
  library(modelr)
  library(patchwork)
  library(naniar)

  sv <- read.pcv(sv_file, mode = "wide", singleValueOnly = TRUE, reader = "fread")

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
}
