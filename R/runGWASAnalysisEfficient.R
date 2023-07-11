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
runGWASAnalysisEfficient <- function(bedfile, famFile, phenoFile, phenoCols, projectCode) {
  library(bigsnpr)
  library(tidyverse)
  # Read from bed/bim/fam
  tmpfile <- tempfile()
  snp_readBed(bedfile, backingfile = tmpfile)

  # Attach the "bigSNP" object in R session
  obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))

  # Get aliases for useful slots
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos

  # Exclude Long-Range Linkage Disequilibrium Regions of the human genome
  ind.excl <- snp_indLRLDR(infos.chr = CHR, infos.pos = POS)

  # Use clumping to keep SNPs weakly correlated with each other
  ind.keep <- snp_clumping(G, infos.chr = CHR, exclude = ind.excl, ncores = nb_cores())

  # Get the specified number of PCs, corresponding to pruned SNPs
  obj.svd <- big_randomSVD(G, fun.scaling = snp_scaleBinom(), ind.col = ind.keep, ncores = nb_cores())
  PC <- predict(obj.svd)
  covar <- PC[, 1:6]
  # Read phenotype data
  pheno <- read.csv(phenoFile, header = TRUE)
  colnames(pheno)[1] <- "genotype"
  fam <- read.table(famFile, sep = "\t", header = FALSE)
  colnames(fam)[2] <- "genotype"
  fam <- fam %>% select(genotype)
  pheno <- left_join(fam, pheno, by = "genotype")

  # Order the phenotypes based on sample IDs
  ord <- match(obj.bigSNP$fam$sample.ID, pheno$genotype)
  pheno2 <- pheno[ord, ]
  obj.bigSNP$fam <- cbind(obj.bigSNP$fam, pheno[-c(1, 2)])


  # Perform GWAS analysis for each phenotype column
  for (col in phenoCols) {
    y <- obj.bigSNP$fam[[col]]
    ind.gwas <- which(!is.na(y) & complete.cases(covar))

    gwas <- big_univLinReg(G, y[ind.gwas], ind.train = ind.gwas, covar.train = covar[ind.gwas, ], ncores = nb_cores())


    obj.gwas.gc <- snp_gc(gwas)

    res_file <- paste0(projectCode, colnames(pheno2)[col],"_gwas.csv")
    out <- as.data.frame(cbind(CHR,POS,obj.gwas.gc))
    out <- out %>% filter(-log10(score)>5)
    write.csv(out,res_file,row.names = FALSE)


  }
  print("Files written out!")

}
