#' @title Evaluate smooths from a fitted gdam model at covariate values
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is basically a wrapper around the [gratia::smooth_estimates()] function for valuating a smooth at a grid of evenly spaced value over the range of the covariate, but designed to work with `gdam` objects. Please see the associated help file in the `gratia` package for more details.
#'
#' @param object A fitted `gdam` object.
#' @param select A character vector of smooths to evaluate. If `NULL`, all smooths are evaluated.
#' @param n The number of points at which to evaluate the smooths.
#' @param n_3d The number of points over the range of the last covariate must in a 3D plot; see [gratia::smooth_estimates()] for details.
#' @param n_4d The number of points over the range of the last covariate in a 4D plot; see [gratia::smooth_estimates()] for details.
#' @param data A data frame containing the values of the covariates at which to evaluate the smooths. If `NULL`, the original data is used.
#' @param dist Numeric; if greater than zero, this is used to determine when a location is too far from data to be plotted when plotting 2D smooths; see [gratia::smooth_estimates()] for details.
#' @param partial_match Logical. In the case of character `select`, should `select` match partially against `smooths`? If `TRUE` then `select` only be a single string.
#' @param ... Currently not used.
#'
#' @return
#' A data frame which is of class `smooth_estimates`.
#'
#' @note
#' Acknowledgments to Gavin Simpson and the `gratia` package for the original [gratia::smooth_estimates()] function.
#'
#' @author
#' Francis K.C. Hui <fhui28@gmail.com>
#'
#' @examples
#' \donttest{
#' # See the main `gdam` help file for examples, particularly for constructing partial smooth plots.
#' }
#'
#' @export
#' @md
smooth_estimates <- function(object, ...) {
    UseMethod("smooth_estimates")
    }

#' @export
#' @rdname smooth_estimates
#' @importFrom gratia smooth_estimates
#' @md
smooth_estimates.gdam <- function(object,
                                  select = NULL,
                                  n = 100,
                                  n_3d = 16,
                                  n_4d = 4,
                                  data = NULL,
                                  dist = NULL,
                                  partial_match = FALSE,
                                  ...) {
    if(!inherits(object, "gdam"))
        stop("model must be of class \"gdam\".")

    ##-----------------------
    #' ## Duplicate GAM object but remove as many things as you do not need (mainly done so that gratia does not take in the wrong information)
    ##-----------------------
    new_gam <- object$gamObject
    new_gam$aic <- NULL
    new_gam$assign <- NULL
    new_gam$boundary <- NULL
    new_gam$call <- NULL
    new_gam$control <- NULL
    new_gam$converged <- NULL
    new_gam$db.drho <- NULL
    new_gam$deviance <- NULL
    new_gam$df.null <- NULL
    new_gam$df.residual <- NULL
    new_gam$edf <- NULL
    new_gam$edf1 <- NULL
    new_gam$edf2 <- NULL
    new_gam$full.sp <- NULL
    new_gam$F <- NULL
    new_gam$gcv.ubre <- NULL
    new_gam$hat <- NULL
    new_gam$iter <- NULL
    new_gam$method <- NULL
    new_gam$mgcv.conv <- NULL
    new_gam$min.edf <- NULL
    new_gam$nsdf <- NULL
    new_gam$null.deviance <- NULL
    new_gam$optimizer <- NULL
    new_gam$outer.info <- NULL
    new_gam$paraPen <- NULL
    new_gam$pterms <- NULL
    new_gam$R <- NULL
    new_gam$rank <- NULL
    new_gam$reml.scale <- NULL
    new_gam$residuals <- NULL
    new_gam$rV <- NULL
    new_gam$var.summary <- NULL
    new_gam$Ve <- NULL
    new_gam$Vp <- NULL
    new_gam$Vc <- NULL
    new_gam$weights <- NULL


    ##-----------------------
    #' ## Modify relevant items and pass to gratia::smooth_estimates
    ##-----------------------
    new_gam$coefficients <- object$coefficients
    new_gam$sp <- object$sp
    new_gam$scale <- object$sigma2
    new_gam$sig2 <- object$sigma2
    new_gam$fitted.values <- object$fitted_values
    new_gam$linear.predictors <- object$linear_predictor
    new_gam$Vp <- object$covariance_matrix

    out <- gratia::smooth_estimates(object = new_gam,
                                    select = select,
                                    n = n,
                                    n_3d = n_3d,
                                    n_4d = n_4d,
                                    data = data,
                                    overall_uncertainty = TRUE,
                                    dist = dist,
                                    partial_match = partial_match)
    return(out)
    }

