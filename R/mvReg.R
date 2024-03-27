#' Multivariate Regression
#'
#' This function uses a least squares approach to fit a multivariate linear model to the input data.
#'
#' @param data multivariate input data including a target varialbe and related independent variables
#' @param target_index the column index of the target variable
#' @param ind_index an array of the independent variables by their respective column indices
#'
#' @return the calculated multivariate linear regression  model and a plot showing the fitment of the model with the data
#' @export
#' @importFrom graphics abline
#' @importFrom stats as.formula fitted lm
#'
#' @examples
#' data <- T7_1
#' mvReg(data, 3, c(1,2))
mvReg <- function(data, target_index=NULL, ind_index=NULL) {

  if (!dim(data)[2] > 1) {
    rlang::abort("The data provided is not multivariate")
    return(NULL)
  }

  if (!all(apply(data, 2, function(x) all(is.numeric(x))))) {
    rlang::abort("The data provided is not numeric")
    return(NULL)
  }

  if (is.null(target_index)) {
    rlang::abort("No target index provided")
    return(NULL)
  }

  if (floor(target_index) != target_index) {
    rlang::abort("Target index is not an integer")
    return(NULL)
  }

  if (target_index > dim(data)[2] || target_index < 1) {
    rlang::abort("Target index is not within provided data")
    return(NULL)
  }

  if (is.null(ind_index)) {
    rlang::abort("No independent indices provided")
    return(NULL)
  }

  if (!all(floor(ind_index) == ind_index)) {
    rlang::abort("Independent indices are not integers")
    return(NULL)
  }

  if (any(ind_index) > dim(data)[2] || any(ind_index) < 1) {
    rlang::abort("Independent indices are not within provided data")
    return(NULL)
  }

  data <- as.data.frame(data)

  formula_str <- paste("data[,target_index]", "~", paste("data[,",ind_index, "]",collapse = " + "))
  formula <- as.formula(formula_str)

  lm_model <- lm(formula, data=data)

  coeff <- summary(lm_model)$coefficients[0:length(ind_index)+1]

  stds <- summary(lm_model)$coefficients[(length(ind_index)+2):(2*length(ind_index)+2)]

  r_squared <- summary(lm_model)$r.squared

  plot <- plot(data[, target_index], fitted(lm_model), xlab = "Observed Values", ylab = "Fitted Values", main = "Fitted vs. Observed") +
  abline(0, 1, col = "red")

  results_list <- list(fitted_model = lm_model, coefficients = coeff, standardDev = stds, rSquared = r_squared, fitPlot = plot)

  return(results_list)
}
