#' Run GWAS analysis efficiently
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
runGWASAnalysisEfficient2 <- function (bedfile, famFile, phenoFile, phenoCols, projectCode) {
  library(bigsnpr)
  library(tidyverse)
  tmpfile <- tempfile()
  snp_readBed(bedfile, backingfile = tmpfile)
  obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))


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
    #covar <- PC[, pcs_to_include]
    gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas,
                           covar.train = PC[, pcs_to_include][ind.gwas, ],
                           ncores = nb_cores())
    obj.gwas.gc <- snp_gc(gwas)
    res_file <- paste0(projectCode, colnames(pheno2)[col],
                       "_gwas.csv")
    rds_file <- paste0(projectCode, colnames(pheno2)[col],
                       "_gwas.rds.gz")
    saveRDS(obj.gwas.gc, file = rds_file,compress ="xz",xz=9 )
    out <- as.data.frame(cbind(CHR, POS, obj.gwas.gc))
    out <- out %>% filter(-log10(score) > 5)
    write.csv(out, file=gzfile(res_file))

  }
  # Free up memory and remove unnecessary objects
  unlink(tmpfile)
  rm(tmpfile, obj.bigSNP, G, CHR, POS, ind.excl, ind.keep, obj.svd, PC, covar, pheno, fam, pheno2, gwas, obj.gwas.gc, out)

  print("Files written out!")
}
