#' Multivariate Normality Check
#'
#' This function checks a multivariate numerical data set for normality by utilizing the squared generalized distances
#' and their relationship to specific chi-square values.
#'
#' @param data 2D or greater input matrix or data frame containing multivariate numerical data
#'
#' @return chi-sq plots for the original and transformed data, the optimal lambda values used for transformation, and the resulting transformed data
#' @export
#' @importFrom stats cov qchisq optimize
#'
#' @examples
#' vec1 <- rexp(100)
#' vec2 <- rexp(100)
#' vec3 <- rexp(100)
#' data <- data.frame(x=vec1, y=vec2, z=vec3)
#' mvNormCheck(data)
mvNormCheck <- function(data) {
  data_mat <- as.matrix(data)

  if (!dim(data)[2] > 1) {
    rlang::abort("The data provided is not multivariate")
    return(NULL)
  }

  if (!all(apply(data, 2, function(x) all(is.numeric(x))))) {
    rlang::abort("The data provided is not numeric")
  }

  if (!all(data >0)) {
    rlang::abort("The data provided is not all positive")
  }

  xbar <- colMeans(data_mat)
  x_minus_xbar <- sweep(data_mat, 2, xbar, "-")
  cov_data <- cov(data_mat)
  inv_cov_data <- solve(cov_data)
  sq_gen_dist <- vector(mode="numeric", length=dim(data)[1])

  for (i in 1:dim(data_mat)[1]) {
    x_minus_xbar_i <- x_minus_xbar[i,]
    sq_gen_dist[i] <- t(x_minus_xbar_i) %*% inv_cov_data %*% x_minus_xbar_i
  }
  sorted_sq_gen_dist <- sort(sq_gen_dist)

  deg_free <- dim(data_mat)[2]
  n <- dim(data_mat)[1]
  quantiles <- vector(mode="numeric", length=dim(data)[1])
  for (i in 1:length(quantiles)) {
    quantiles[i] <- qchisq((n - i + 0.5) / n, df = deg_free)
  }
  sorted_quantiles <- sort(quantiles)

  xlimit_max <- max(sorted_quantiles) + 0.5
  xlimit_min <- 0

  ylimit_max <- max(sorted_sq_gen_dist) + 0.5
  ylimit_min <- 0

  df_orig <- data.frame(x = sorted_quantiles, y = sorted_sq_gen_dist)
  orig_plot <- ggplot2::ggplot(df_orig, ggplot2::aes(x = sorted_quantiles, y = sorted_sq_gen_dist)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Quantiles",
                  y = "Squared Generalized Distance",
                  title = "Chi-Square Plot (Original)") +
    ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") +
    ggplot2::xlim(xlimit_min, xlimit_max) +
    ggplot2::ylim(ylimit_min, ylimit_max) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  x_lambda_4_34 <- function(data_column, lambda) {
    augmented_data <- numeric(length(data_column))
    for (i in seq_along(data_column)) {
      if (lambda != 0) {
        augmented_data[i] <- ((data_column[i]^lambda) - 1)/(lambda)
      }
      else if (lambda == 0) {
        augmented_data[i] <- log(data_column[i])
      }
    }
    return(augmented_data)
  }

  xbar_lambda_4_36 <- function(data_column, lambda) {
    x_lambdas <- x_lambda_4_34(data_column, lambda)
    xbar_lambda <- sum(x_lambdas)/length(x_lambdas)
    return(xbar_lambda)
  }

  likelihood_4_35 <- function(lambda, data_column) {
    x_lambdas <- x_lambda_4_34(data_column, lambda)
    xbar_lambda <- xbar_lambda_4_36(data_column, lambda)
    n <- length(x_lambdas)
    likelihood_func <- ((-n/2) * log((1/n) * sum((x_lambdas - xbar_lambda)^2))) + (lambda - 1) * sum(log(data_column))
    return(likelihood_func)
  }

  data <- as.data.frame(data)
  column_names <- colnames(data)

  augmented_data <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
  colnames(augmented_data) <- paste0(colnames(data), "_aug")

  optimal_lambda_list <- list()
  for (col_name in column_names) {
    data_column <- data[[col_name]]
    maximum <- optimize(f = function(lambda) likelihood_4_35(lambda, data_column),
                        interval = c(-5, 5), maximum = TRUE)
    augmented_data[[paste0(col_name, "_aug")]] <- x_lambda_4_34(data_column, maximum$maximum)
    optimal_lambda_list[[col_name]] <- round(maximum$maximum, 3)
  }

  likelihood_data <- data.frame(matrix(ncol = length(column_names), nrow = 1000))
  colnames(likelihood_data) <- paste0(column_names, "_l(lambda)")
  lambda_values <- seq(-5, 5, length.out = 1000)

  for (col_name in column_names) {
    data_column <- data[[col_name]]
    likelihood_values <- sapply(lambda_values, function(lambda) likelihood_4_35(lambda, data_column))
    new_col_name <- paste0(col_name, "_l(lambda)")
    likelihood_data[[new_col_name]] <- likelihood_values
  }

  # lambda_plots <- list()
  # for (i in 1:ncol(likelihood_data)) {
  #   col_name <- colnames(likelihood_data)[i]
  #
  #   max_index <- which.max(likelihood_data[[i]])
  #   optimal_lambda <- lambda_values[max_index]
  #
  #   percentile_05 <- quantile(likelihood_data[, i], probs = 0.05)
  #
  #   p <- ggplot2::ggplot(likelihood_data, ggplot2::aes(x = lambda_values, y = likelihood_data[,i])) +
  #     ggplot2::geom_line() +
  #     ggplot2::geom_vline(xintercept = optimal_lambda, linetype = "dashed", color = "red") +
  #     ggplot2::annotate("text", x = optimal_lambda - 0.125, y = percentile_05,
  #                       label = paste("Optimal value of lambda =", round(optimal_lambda, 3)),
  #                       color = "black", angle = 90, hjust = -0.2) +
  #     ggplot2::labs(title = col_name, x = "lambda", y = "l(lambda)") +
  #     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  #
  #   lambda_plots[[i]] <- p
  # }
  #
  # # Arrange plots in a grid
  # ncols_lambda <- floor(sqrt(ncol(likelihood_data)))
  # nrows_lambda <- ceiling(ncol(likelihood_data) / ncols_lambda)
  #
  # lambda_plots_arranged <- gridExtra::grid.arrange(
  #   grobs = lambda_plots,
  #   ncol = ncols_lambda,
  #   nrow = nrows_lambda,
  #   widths = rep(1/ncols_lambda, ncols_lambda),
  #   heights = rep(1, nrows_lambda),
  #   padding = unit(0.01, "line")
  # )


  aug_data_mat <- as.matrix(augmented_data)

  aug_xbar <- colMeans(aug_data_mat)
  aug_x_minus_xbar <- sweep(aug_data_mat, 2, aug_xbar, "-")
  aug_cov_data <- cov(aug_data_mat)
  aug_inv_cov_data <- solve(aug_cov_data)
  aug_sq_gen_dist <- vector(mode="numeric", length=dim(aug_data_mat)[1])

  for (i in 1:dim(aug_data_mat)[1]) {
    aug_x_minus_xbar_i <- aug_x_minus_xbar[i,]
    aug_sq_gen_dist[i] <- t(aug_x_minus_xbar_i) %*% aug_inv_cov_data %*% aug_x_minus_xbar_i
  }
  aug_sorted_sq_gen_dist <- sort(aug_sq_gen_dist)

  deg_free <- dim(aug_data_mat)[2]
  n <- dim(aug_data_mat)[1]
  aug_quantiles <- vector(mode="numeric", length=dim(data)[1])
  for (i in 1:length(aug_quantiles)) {
    aug_quantiles[i] <- qchisq((n - i + 0.5) / n, df = deg_free)
  }
  aug_sorted_quantiles <- sort(aug_quantiles)

  xlimit_max <- max(aug_sorted_quantiles) + 0.5
  xlimit_min <- 0

  ylimit_max <- max(aug_sorted_sq_gen_dist) + 0.5
  ylimit_min <- 0

  df_aug <- data.frame(x = aug_sorted_quantiles, y = aug_sorted_sq_gen_dist)
  aug_plot <- ggplot2::ggplot(df_aug, ggplot2::aes(x = aug_sorted_quantiles, y = aug_sorted_sq_gen_dist)) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Quantiles",
         y = "Squared Generalized Distance",
         title = "Chi-Square Plot (Transformed)") +
    ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") +
    ggplot2::xlim(xlimit_min, xlimit_max) +
    ggplot2::ylim(ylimit_min, ylimit_max) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  chisq_plots <- list()
  chisq_plots[[1]] <- orig_plot
  chisq_plots[[2]] <- aug_plot


  ncols_chisq <- 2
  nrows_chisq <- 1

  chisq_plots_arranged <- gridExtra::grid.arrange(
    grobs = chisq_plots,
    ncol = ncols_chisq,
    nrow = nrows_chisq,
    widths = rep(1/ncols_chisq, ncols_chisq),
    heights = rep(0.25, nrows_chisq)
  )

  results_list <- list(
    optimal_lambdas = optimal_lambda_list,
    transformed_data = augmented_data
  )

  return(results_list)
}
