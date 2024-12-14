#' @title Model diagnostics plots for gdam objects
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is basically a wrapper around the [gratia::appraise()] function, but designed to work with `gdam` objects. The main difference is that observations (residuals, quantiles etc...) in the diagnostic plots are colored by the final weights from the `gdam` object.
#'
#' @param model A fitted `gdam` object.
#' @param use_worm Logical. Should the worm plot be drawn in place of the QQ plot? Default is `FALSE`.
#' @param n_uniform Integer. The number of times to randomize uniform quantiles in the direct computation method or QQ plots.
#' @param n_simulate Integer. The number of times to simulate from the estimated model in the simulation-based method for QQ plots.
#' @param method Character. The method to use for to generate theoretical quantiles. Options are "uniform" or "simulate".
#' @param type Character. The type of residuals to use. Options are "deviance" or "pearson" residuals.
#' @param n_bins Character. The number of bins to use in the histogram of residuals. See [gratia::appraise()] for more details.
#' @param ncol Integer. The number of columns in the grid of plots. If `NULL`, the number of columns is determined automatically.
#' @param nrow Integer. The number of rows in the grid of plots. If `NULL`, the number of rows is determined automatically.
#' @param guides Character. The guide to use for the color scale. See [gratia::appraise()] for more details
#' @param level Numeric. The confidence level for the reference intervals in the QQ plot. Only used if `method = "simulate"`.
#' @param ci_col Character. The color of the reference interval lines in the QQ plot. Only used if `method = "simulate"`.
#' @param ci_alpha Numeric. The alpha transparency level of the reference interval lines in the QQ plot. Only used if `method = "simulate"`.
#' @param line_col Character. The color of the reference 1:1 line in the QQ plot, and the reference line in the residuals vs. linear predictor plot.
#' @param ... Currently not used.
#'
#'
#' @note
#' Acknowledgments to Gavin Simpson and the `gratia` package for the original [gratia::appraise()] function.
#'
#' @author
#' Francis K.C. Hui <fhui28@gmail.com>
#'
#' @examples
#' \dontrun{
#' See the main `gdam` help file for examples.
#' }
#'
#' @export
#' @importFrom gratia appraise
#' @importFrom colorspace scale_color_continuous_sequential
#' @importFrom ggplot2 geom_point theme_bw aes_string
#' @md

appraise.gdam <- function(model,
                          use_worm = FALSE,
                          n_uniform = 10,
                          n_simulate = 50,
                          method = "uniform",
                          type = "deviance",
                          n_bins = "sturges",
                          ncol = NULL,
                          nrow = NULL,
                          guides = "keep",
                          level = 0.9,
                          ci_col = "black",
                          ci_alpha = 0.2,
                          line_col = "red",
                          ...) {

    if(!inherits(model, "gdam"))
        stop("model must be of class \"gdam\".")

    ##-----------------------
    #' ## Duplicate GAM object but remove as many things as you do not need (mainly done so that gratia does not take in the wrong information)
    ##-----------------------
    new_gam <- model$gamObject
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
    new_gam$smooth <- NULL
    new_gam$var.summary <- NULL
    new_gam$Ve <- NULL
    new_gam$Vp <- NULL
    new_gam$Vc <- NULL
    new_gam$weights <- NULL


    ##-----------------------
    #' ## Modify relevant items and pass to gratia::appraise
    ##-----------------------
    new_gam$coefficients <- model$coefficients
    new_gam$sp <- model$sp
    new_gam$scale <- model$sigma2
    new_gam$sig2 <- model$sigma2
    new_gam$fitted.values <- model$fitted_values
    new_gam$linear.predictors <- model$linear_predictor
    new_gam$Vp <- model$covariance_matrix

    make_plots <- gratia::appraise(new_gam,
                                   method = method,
                                   type = type,
                                   use_worm = use_worm,
                                   n_uniform = n_uniform,
                                   n_simulate = n_simulate,
                                   n_bins = n_bins,
                                   ncol = ncol,
                                   nrow = nrow,
                                   guides = guides,
                                   level = level,
                                   ci_col = ci_col,
                                   ci_alpha = ci_alpha,
                                   line_col = line_col
    )

    make_plots[[1]]$data$weights <- model$final_weights
    make_plots[[1]]$layers[[2]] <- NULL
    make_plots[[1]] <- make_plots[[1]] +
        ggplot2::geom_point(aes_string(x = 'theoretical', y = 'residuals', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    make_plots[[2]]$data$weights <- model$final_weights
    make_plots[[2]]$layers[[2]] <- NULL
    make_plots[[2]] <- make_plots[[2]] +
        ggplot2::geom_point(aes_string(x = 'eta', y = 'residuals', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    make_plots[[3]] <- make_plots[[3]] +
        ggplot2::theme_bw()

    make_plots[[4]]$data$weights <- model$final_weights
    make_plots[[4]]$layers[[2]] <- NULL
    make_plots[[4]] <- make_plots[[4]] +
        ggplot2::geom_point(aes_string(x = 'fitted', y = 'observed', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    print(make_plots)
    }


