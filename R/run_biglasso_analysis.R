#awk '$1~/^Chr02/ {print}' TM021.vcf > tmp.vcf
#cat header.vcf tmp.vcf > chr02_tm021.vcf
#create header.vcf
#awk '$1~/^#/{print}' TM027_028.vcf > header.vcf
#parallel -j 4 --dry-run "awk '\$1~/^{}/ {{print}}' TM021.vcf > tmp_{}.vcf; cat header.vcf tmp_{}.vcf > chr{}_tm021.vcf" ::: Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
#parallel -j 4 "awk '\$1~/^{}/ {{print}}' TM021.vcf > tmp_{}.vcf; cat header.vcf tmp_{}.vcf > chr{}_tm021.vcf" ::: Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10

#parallel -j 4  "awk '\$1~/^{}/ {{print}}' TM027_028.vcf > tmp_{}.vcf; cat header.vcf tmp_{}.vcf > chr{}_tm027_028.vcf" ::: Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10

#Rscript workflow.R > output.txt 2> errors.txt &
#
#Rscript lasso.R > lasso_output.txt 2> lasso_errors.txt &

library(biglasso)
library(data.table)
library(tidyverse)
library(vcfR)
library(parallel)
library(furrr)

run_biglasso_analysis <- function(vcf_file, pheno_file, string, outname) {
  tryCatch(
    {
      # Read VCF file
      print("Reading VCF file...")
      tm021_vcf <- read.vcfR(vcf_file, convertNA = TRUE)
      
      # Extract genotype information
      print("Extracting genotype information...")
      tm021_vcf_num <- extract.gt(tm021_vcf,
                                  element = "GT",
                                  IDtoRowNames  = TRUE,
                                  as.numeric = TRUE,
                                  convertNA = TRUE)
      
      tm021_vcf_num_t <- t(tm021_vcf_num) %>% as.data.frame()
      tm021_vcf_num_t$genotype <- rownames(tm021_vcf_num_t)
      
      # Read phenotype file and merge with genotype information
      print("Reading phenotype file and merging with genotype information...")
      pheno <- fread(pheno_file, data.table = FALSE)
      
      merged1 <- pheno %>% 
        left_join(tm021_vcf_num_t, by = "genotype")
      
      merged1 <- na.omit(merged1)
      
      # Select columns containing the specified string
      print("Selecting columns based on the specified string...")
      merged2 <- merged1 %>% select(contains(string))
      
      # Select columns starting with "Chr" for predictors
      print("Selecting predictor columns starting with 'Chr'...")
      X <- merged1 %>% select(starts_with("Chr"))
      
      X.bm <- as.big.matrix(X)
      
      # Loop over response columns and perform biglasso analysis
      for (response_col in colnames(merged2)) {
        print(paste("Running biglasso analysis for response variable:", response_col))
        
        y <- merged2[[response_col]]
        
        fit <- biglasso(X.bm, y)
        
        # 10-fold cross-validation in parallel
        print("Running cross-validation...")
        cvfit <- tryCatch(
          {
            cv.biglasso(X.bm, y, seed = 1234, nfolds = 10, ncores = detectCores() - 1)
          },
          error = function(cond) {
            message("An error occurred during cross-validation. Trying with fewer cores.")
            cv.biglasso(X.bm, y, seed = 1234, nfolds = 10, ncores = detectCores() - 2)
          }
        )
        
        # Save cross-validation results
        print("Saving cross-validation results...")
        saveRDS(cvfit, paste0(outname, "_cvfit_", response_col), compress = TRUE)
        
        # Save non-zero coefficients
        non_zero <- coef(cvfit)[which(coef(cvfit) != 0),]
        print("Saving non-zero coefficients...")
        saveRDS(non_zero, paste0(outname, "_non_zero_", response_col), compress = TRUE)
      }
    },
    error = function(cond) {
      message("An error occurred. Please check your inputs.")
    }
  )
}

# Example usage:
# run_biglasso_analysis(
#   vcf_file = "gwas/rawdata/chrChr01_tm021.vcf",
#   pheno_file = "TM021_nov_T100_median_traits.csv",
#   string = "height",  # Replace with the actual string
#   outname = "tm021_dec28_chsr01"
# )



# Plan for parallelization (e.g., using multicore)

# Create a list of inputs for each run
inputs <- lapply(1:2, function(i) {
  list(
    vcf_file = sprintf("gwas/rawdata/chrChr%02d_tm021.vcf", i),
    pheno_file = "TM021_nov_T100_median_traits.csv",
    string = "height",  # Replace with the actual string
    outname = sprintf("tm021_dec28_chr%02d", i)
  )
})

# Run the function in parallel using furrr
future_map(inputs, ~ run_biglasso_analysis(.x$vcf_file, .x$pheno_file, .x$string, .x$outname))
