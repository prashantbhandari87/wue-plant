#' Calculate Heritability for a Treatment
#'
#' This function calculates the heritability for a given treatment using linear mixed-effects modeling.
#'
#' @param obj The data frame containing the phenotype data.
#' @param pheno.col.no The column number of the phenotype of interest in the data frame.
#'
#' @return The heritability value for the treatment or NA if an error occurs.
#'
#' @import lme4
#' @import tidyverse
#' @import data.table
#' @import foreach
#' @import doParallel
#' @import ggsci
#' @import bigsnpr
#'
#' @export
get_h2_treatment <- function(obj, pheno.col.no) {
  library(lme4)
  library(tidyverse)
  library(data.table)
  library(foreach)
  library(doParallel)
  library(ggsci)
  library(bigsnpr)

  # Wrap the code inside tryCatch
  tryCatch({
    m <- lmer(obj[, pheno.col.no] ~ (1|genotype), data = obj)
    summary(m)
    re <- as.numeric(VarCorr(m))
    res <- attr(VarCorr(m), "sc")^2
    geno.var <- VarCorr(m)$genotype[1]
    tot.var <- sum(sum(re), res)
    h <- geno.var / tot.var
    return(h)
  }, error = function(e) {
    # If an error occurs (e.g., heritability cannot be calculated),
    # return NA and print the error message
    print(paste("Error:", e))
    return(NA)
  })
}
