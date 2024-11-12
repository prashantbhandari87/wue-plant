#' Perform clustering and extract data for each cluster
#'
#' This function takes a hierarchical clustering object, a pivoted genotype data frame,
#' and a joined data frame, and performs clustering to create k clusters. It then extracts
#' the corresponding data for each cluster and returns a list of data frames.
#'
#' @param hc_sbd Hierarchical clustering object.
#' @param dat_pivoted_geno Pivoted genotype data frame.
#' @param joined_data Joined data frame containing information about plants.
#' @param k Number of clusters to create.
#'
#' @return A list of data frames, where each data frame contains information about plants
#'         belonging to a specific cluster.
#'
#' @export
#'
#'
#' @seealso
#' \code{\link{cutree}}, \code{\link{future_map}}, \code{\link{filter}}
#'
#' @import dplyr
#' @import furrr
performCluster <- function(hc_sbd, dat_pivoted_geno, joined_data, k = 10) {

  k_clusters <- cutree(hc_sbd, k = k)

  # Generate cluster indices in parallel
  cluster_indices <- future_map(1:k, ~which(k_clusters == .x))

  # Generate filtered data frames in parallel
  dat_clusters <- future_map(1:k, ~dat_pivoted_geno[cluster_indices[[.x]], ])

  cluster_sv <- future_map(1:k, ~joined_data %>%
                             filter(plantbarcode %in% dat_clusters[[.x]]$plantbarcode))

  return(cluster_sv)

}
