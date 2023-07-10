#' Create project folders
#'
#' This function creates the project directories and subdirectories if they don't already exist.
#'
#' @param projectDir Character string specifying the project directory path.
#'
#' @return None
#' @export
#'
#' @examples
#' createProjectFolders("path/to/project")
createProjectFolders <- function(projectDir) {
  # Create project directories if they don't already exist
  createDirectory <- function(directory) {
    if (!dir.exists(directory)) {
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      cat("Folder created:", directory, "\n")
    } else {
      cat("Folder already exists:", directory, "\n")
    }
  }

  # Create project subdirectories
  pcvrDir <- file.path(projectDir, "pcvr")
  createDirectory(pcvrDir)

  pcvrRawdataDir <- file.path(pcvrDir, "rawdata")
  createDirectory(pcvrRawdataDir)

  pcvrCleandataDir <- file.path(pcvrDir, "cleandata")
  createDirectory(pcvrCleandataDir)

  pcvrPlotsDir <- file.path(pcvrDir, "plots")
  createDirectory(pcvrPlotsDir)

  pcvrReportsDir <- file.path(pcvrDir, "reports")
  createDirectory(pcvrReportsDir)

  pcvrMiscDir <- file.path(pcvrDir, "misc")
  createDirectory(pcvrMiscDir)

  pcvrScriptsDir <- file.path(pcvrDir, "scripts")
  createDirectory(pcvrScriptsDir)

  mashDir <- file.path(projectDir, "mash")
  createDirectory(mashDir)

  mashRawdataDir <- file.path(mashDir, "rawdata")
  createDirectory(mashRawdataDir)

  mashCleandataDir <- file.path(mashDir, "cleandata")
  createDirectory(mashCleandataDir)

  mashPlotsDir <- file.path(mashDir, "plots")
  createDirectory(mashPlotsDir)

  mashReportsDir <- file.path(mashDir, "reports")
  createDirectory(mashReportsDir)

  mashMiscDir <- file.path(mashDir, "misc")
  createDirectory(mashMiscDir)

  mashScriptsDir <- file.path(mashDir, "scripts")
  createDirectory(mashScriptsDir)

  gwasDir <- file.path(projectDir, "gwas")
  createDirectory(gwasDir)

  gwasRawdataDir <- file.path(gwasDir, "rawdata")
  createDirectory(gwasRawdataDir)

  gwasCleandataDir <- file.path(gwasDir, "cleandata")
  createDirectory(gwasCleandataDir)

  gwasPlotsDir <- file.path(gwasDir, "plots")
  createDirectory(gwasPlotsDir)

  gwasReportsDir <- file.path(gwasDir, "reports")
  createDirectory(gwasReportsDir)

  gwasMiscDir <- file.path(gwasDir, "misc")
  createDirectory(gwasMiscDir)

  gwasScriptsDir <- file.path(gwasDir, "scripts")
  createDirectory(gwasScriptsDir)

  # Print completion message
  cat("Project folders created successfully.\n")
}
