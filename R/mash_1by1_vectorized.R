#' Perform MASH analysis for multiple columns simultaneously
#'
#' This function performs MASH analysis using the ashr package for multiple columns
#' of Bhat and Shat matrices simultaneously. It calculates posterior mean,
#' posterior standard deviation, local false sign rate (lfsr), and log likelihood.
#'
#' @param data A list containing Bhat and Shat matrices.
#' @param alpha The tuning parameter controlling the amount of sparsity. Default is 0.
#'
#' @return A list containing posterior matrices and log likelihood.
#' @export
#'


mash_1by1_vectorized <- function(data, alpha = 0, ...) {
  Bhat <- data$Bhat
  Shat <- data$Shat
  post_mean = post_sd = lfsr = matrix(nrow = nrow(Bhat), ncol= ncol(Bhat))
  loglik = 0

  # Validate input dimensions to ensure Bhat and Shat have the same number of rows
  if (nrow(Bhat) != nrow(Shat)) {
    stop("Error: The number of rows in Bhat and Shat must be the same.")
  }

  nr <- nrow(Bhat)
  ncores <- detectCores()-1

  # Apply the 'ash' function to each column of 'Bhat' and 'Shat' simultaneously
  ash_results <- mclapply(1:ncol(Bhat), function(i) {
    ashres <- ashr::ash(Bhat[, i], Shat[, i], mixcompdist = "normal", alpha = alpha,...)
    list(post_mean = ashr::get_pm(ashres),
         post_sd = ashr::get_psd(ashres),
         lfsr = ashr::get_lfsr(ashres),
         loglik = ashr::get_loglik(ashres))
  },mc.cores =ncores)

  # Extract the posterior matrices and log likelihood
  post_mean <- do.call(cbind, mclapply(ash_results, function(res) res$post_mean,mc.cores =ncores))
  post_sd <- do.call(cbind, mclapply(ash_results, function(res) res$post_sd,mc.cores =ncores))
  lfsr <- do.call(cbind, mclapply(ash_results, function(res) res$lfsr,mc.cores =ncores))
  loglik <- sum(sapply(ash_results, function(res) res$loglik))

  # Set column and row names for the posterior matrices
  colnames(post_mean) <- colnames(post_sd) <- colnames(lfsr) <- colnames(Bhat)
  rownames(post_mean) <- rownames(post_sd) <- rownames(lfsr) <- rownames(Bhat)

  # Create a list with the posterior matrices and log likelihood
  posterior_matrices <- list(
    PosteriorMean = post_mean,
    PosteriorSD = post_sd,
    lfsr = lfsr
  )

  # Create the output list
  output_list <- list(result = posterior_matrices, loglik = loglik)

  # Set the class attribute for the output list
  class(output_list) <- "mash_1by1"

  return(output_list)
}

