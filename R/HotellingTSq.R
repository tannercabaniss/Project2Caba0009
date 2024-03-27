#' Hotelling T Squared Test
#'
#' This function tests a provided dataset against a hypothesized mean vector for conformity.
#'
#' @param data the multivariate input data in matrix or dataframe form to be evaluated
#' @param H0 the hypothesized mean vector
#' @param alpha the alpha value used in the determination of the critical value
#'
#' @return True if the null hypothesis is accepted and false if the null hypothesis is rejected
#' @export
#' @importFrom stats qf
#'
#' @examples
#' vec1 <- rnorm(100)
#' vec2 <- rnorm(100)
#' vec3 <- rnorm(100)
#' data <- data.frame(x=vec1, y=vec2, z=vec3)
#' H0 <- c(0, 0, 0)
#' alpha <- 0.05
#' HotellingTSq(data, H0, alpha)
HotellingTSq <- function(data, H0, alpha=0.05) {

  if (!dim(data)[2] > 1) {
    rlang::abort("The data provided is not multivariate")
    return(NULL)
  }

  if (!all(apply(data, 2, function(x) all(is.numeric(x))))) {
    rlang::abort("The data provided is not numeric")
    return(NULL)
  }

  data_mat <- as.matrix(data)

  xbar <- colMeans(data_mat)
  xbar_minus_H0 <- xbar - H0
  xbar_minus_H0_t <- t(xbar) - H0
  cov_mat <- cov(data_mat)
  inv_cov_mat <- solve(cov_mat)
  n <- dim(data_mat)[1]
  p <- dim(data_mat)[2]

  T_sq <- n * xbar_minus_H0_t %*% inv_cov_mat %*% xbar_minus_H0

  critical_value = ((n-1)*p/(n-p)) * qf(1 - alpha, p, n - p)

  if (T_sq > critical_value) {
    conclusion <- "Reject null hypothesis"
  } else {
    conclusion <- "Fail to reject null hypothesis"
  }

  confidence_plots <- list()
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      plot <- ggplot2::ggplot(data, ggplot2::aes_string(x = colnames(data)[i], y = colnames(data)[j])) +
        ggplot2::geom_point() +
        ggplot2::stat_ellipse(level = 1-alpha, type = "t", linetype = "dashed", color = "red") +
        ggplot2::labs(title = paste("1-alpha confidence ellipse:", colnames(data)[i], "and", colnames(data)[j])) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

      confidence_plots[[length(confidence_plots) + 1]] <- plot
    }
  }

  ncols_conf <- floor(sqrt(dim(data)[2]))
  nrows_conf <- ceiling(dim(data)[2] / ncols_conf)

  conf_plots_arranged <- gridExtra::grid.arrange(
    grobs = confidence_plots,
    ncol = ncols_conf,
    nrow = nrows_conf,
    widths = rep(1/ncols_conf, ncols_conf),
    heights = rep(1, nrows_conf)
  )

  result <- list(T_squared = T_sq, critical_value = critical_value, conclusion = conclusion)
  return(result)
}
