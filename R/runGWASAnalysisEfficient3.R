#' Run GWAS analysis efficiently. this is for aws cluster
#'
#' This function reads genotype data from a bed/bim/fam file, performs GWAS
#' analysis for each phenotype column, and writes the results to CSV files.
#'
#' @param bedfile The path to the bed/bim/fam file.
#' @param famFile The path to the fam file.
#' @param phenoFile The path to the phenotype file.
#' @param phenoCols The names of the phenotype columns to analyze.
#' @param projectCode The project code to use for the output file names.
#'
#' @return None.
#' @export
#'
#' @import bigsnpr
#' @import tidyverse
#'
runGWASAnalysisEfficient3 <- function (bedfile, famFile, phenoFile, phenoCols, projectCode) {
  library(bigsnpr)
  library(tidyverse)
  library(arrow)
  library(here)

  obj.bigSNP <- snp_attach(list.files(here("gwas/rawdata/"),
                                      pattern = "\\.rds$",
                                      full.names = TRUE))


  G <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  ind.excl <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)
  ind.keep <- snp_clumping(G, infos.chr = CHR, exclude = ind.excl,
                           ncores = nb_cores())
  obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(),
                           ind.col = ind.keep, ncores = nb_cores())
  PC <- predict(obj.svd)
  covar <- PC[, 1:10]
  pheno <- read.csv(phenoFile, header = TRUE)
  colnames(pheno)[1] <- "genotype"
  pheno <- pheno %>% select_if(~n_distinct(.) >= 50)
  read_tab_or_space_sep <- function(file_path) {
    con <- file(file_path, "r")
    first_line <- readLines(con, n = 1)
    close(con)
    if (grepl("\t", first_line)) {
      data <- read.table(file_path, sep = "\t", header = FALSE,
                         fill = TRUE)
    }
    else if (grepl(" ", first_line)) {
      data <- read.table(file_path, sep = " ", header = FALSE,
                         fill = TRUE)
    }
    else {
      stop("File is not separated by either tabs or spaces.")
    }
    return(data)
  }
  fam <- read_tab_or_space_sep(famFile)
  colnames(fam)[2] <- "genotype"
  fam <- fam %>% select(genotype)
  pheno <- left_join(fam, pheno, by = "genotype")
  ord <- match(obj.bigSNP$fam$sample.ID, pheno$genotype)
  pheno2 <- pheno[ord, ]
  obj.bigSNP$fam <- cbind(obj.bigSNP$fam, pheno[-c(1, 2)])
  for (col in phenoCols) {
    y <- obj.bigSNP$fam[[col]]
    ind.gwas <- which(!is.na(y) & complete.cases(covar))
    pcs_to_include = c(1, 2, 3, 4)
    for (j in 5:10) {
      if (cor.test(y, covar[, j])[3] < 0.05) {
        pcs_to_include = c(pcs_to_include, j)
      }
    }

    print(paste0(length(pcs_to_include),"PCS added for column number",col, sep=" "))

    #covar <- PC[, pcs_to_include]
    gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas,
                           covar.train = PC[, pcs_to_include][ind.gwas, ],
                           ncores = nb_cores())
    obj.gwas.gc <- snp_gc(gwas)
    res_file <- paste0(projectCode, colnames(pheno2)[col],
                       "_gwas.csv")
    rds_file <- paste0(projectCode, colnames(pheno2)[col],
                       "_gwas.arrow")

    out <- as.data.frame(cbind(CHR, POS, obj.gwas.gc))
    out$id <- paste0(CHR,":",POS)
    out$trait <- colnames(pheno2)[col]
    order1 <- order(CHR, POS)
    out$pvalue <- stats::predict(obj.gwas.gc)[order1]
    out1 <- arrow_table(out)
    write_dataset(out1, rds_file, partitioning = c("CHR"))


  }


  print("Files written out!")
}
