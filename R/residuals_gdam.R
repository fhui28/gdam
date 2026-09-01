#' @title Residuals from a fitted gdam model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Computes residuals from a fitted `gdam` object. The residuals can be of different types: "response" corresponds to the raw residuals, "deviance" corresponds to the deviance residuals (note these are slightly ad-hoc as the concept of deviance itself is not clearly defined in the case of gamma divergence; here the deviance residual is computed assuming normally distributed errors), "pearson" corresponds to the Pearson residuals, "scaled_pearson" corresponds to the Pearson residuals scaled by the square root of the estimated dispersion parameter, and "working" corresponds to the working residuals (i.e., the residuals from a post-hoc GAM fit, though in this case they should correspond to the raw response residuals).
#'
#' @param object A fitted `gdam` object.
#' @param type The type of residuals to compute. Options are "response", "deviance", "pearson", "scaled_pearson", or "working".
#' @param ... Currently not used.
#'
#' @return A vector of residuals.
#'
#' @author
#' Francis K.C. Hui <fhui28@gmail.com>
#'
#' @examples
#' \donttest{
#' # See the main `gdam` help file for examples.
#' }
#'
#' @export
#' @importFrom mgcv gam residuals.gam
#' @importFrom stats naresid
#' @md

residuals.gdam <- function(object, type = "deviance", ...) {
    if(!inherits(object, "gdam"))
        stop("model must be of class \"gdam\".")

    type <- match.arg(type, choices = c("response", "deviance", "pearson", "scaled_pearson", "working"))

    mu <- object$fitted.values
    wts <- object$prior.weights
    if(type == "response")
        res <- object$gamObject$y - object$fitted_values
    if(type == "deviance") {
        res <- object$gamObject$family$dev.resids(y = object$gamObject$y,
                                                  mu = object$fitted_values,
                                                  wt = rep(1, length(object$fitted_values)))
        s <- attr(res, "sign")
        if(is.null(s))
            s <- sign(object$gamObject$y - object$fitted_values)
        res <- sqrt(pmax(res, 0)) * s
        }
    if(type %in% c("pearson", "scaled_pearson")) {
        var <- object$gamObject$family$variance
        res <- (object$gamObject$y - object$fitted_values) / sqrt(var(object$fitted_values))
        if(type == "scaled_pearson")
            res <- res / sqrt(object$sigma2)
        }
    if(type %in% c("working")) { #' In additive model case, should be the same as raw residuals
        posthoc_mgcv_fit_fn <- function(object, new_weight) {
            dofit <- mgcv::gam(object$gamObject$formula,
                               data = cbind(object$gamObject$model, wts = new_weight),
                               weights = wts,
                               sp = object$sp,
                               scale = object$sigma2,
                               method = "REML")
            return(dofit)
            }
        posthoc_mgcv_fit <- posthoc_mgcv_fit_fn(object = object,
                                                new_weight = object$final_weights)
        res <- mgcv::residuals.gam(posthoc_mgcv_fit, type = "working")
        }

    res <- naresid(object$gamObject$na.action, res)
    return(res)
    }
