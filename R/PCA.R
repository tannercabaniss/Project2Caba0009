#' Principle Component Analysis
#'
#' This function performs a principle component analysis of the provided data.
#'
#' @param data a multivariate data set in matrix or data frame form
#'
#' @return the coefficients of the principle components, the variances of the components, and the cumulative percentage of total variance
#' @export
#' @importFrom stats cov
#'
#' @examples
#' data <- T8_5
#' PCA(data)
PCA <- function (data) {

  if (!dim(data)[2] > 1) {
    rlang::abort("The data provided is not multivariate")
    return(NULL)
  }

  if (!all(apply(data, 2, function(x) all(is.numeric(x))))) {
    rlang::abort("The data provided is not numeric")
    return(NULL)
  }

  data <- as.data.frame(data)

  cov_mat <- cov(data)
  eigens_cov <- eigen(cov_mat)
  eigenvalues <- eigens_cov$values

  cumulative_var <- (cumsum(eigenvalues) / sum(eigenvalues)) * 100

  num_prinp_comps <- 1:length(eigens_cov$values)

  scree_data <- data.frame(num_prinp_comps, cumulative_var)

  scree_plot <- ggplot2::ggplot(scree_data, ggplot2::aes(x = num_prinp_comps, y = cumulative_var)) +
    ggplot2::geom_bar(stat = "identity", fill = "blue") +
    ggplot2::geom_line(ggplot2::aes(x = num_prinp_comps, y = cumulative_var), color = "red") +
    ggplot2::labs(title = "Scree Plot",
         x = "Principal Component",
         y = "Cumulative Variance") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  results_list <- list(coefficients = round(eigens_cov$vectors, 3),
                       variances = round(eigens_cov$values, 2),
                       cumulative = round(cumulative_var, 1),
                       plot = scree_plot)

  return(results_list)
}
