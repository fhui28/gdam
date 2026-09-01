#' @title Predictions from a fitted gdam model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted gam object produced by [gdam()] and produces predictions given a new set of values for the model covariates or the original values used for the model fit. Predictions can be accompanied by standard errors and uncertainty intervals based on a large sample normality assumption for the distribution of the model coefficients. The routine can optionally return the matrix by which the model coefficients must be premultiplied in order to yield the values of the linear predictor.
#'
#' @param object A fitted `gdam` object.
#' @param newdata A data frame containing the values of the covariates at which predictions are required. If `newdata` is omitted, predictions are made at the original data used to fit the model.
#' @param type The type of prediction required. The default is on the scale of the linear predictor; the alternative `response` is on the scale of the response variable...which in the case of the additive models i.e., identity link function, yields the same predictions; the alternative `lpmatrix` returns the matrix by which the model coefficients must be pre-multiplied in order to yield the values of the linear predictor (in this case `se.fit` is ignored).
#' @param se.fit Logical. Should standard errors be returned? Default is `FALSE`.
#' @param coverage Numeric. The coverage of the uncertainty intervals. Default is 0.95.
#' @param na.action A function which indicates what should happen when `newdata` contain `NA` values. The default is `na.pass`; see [mgcv::predict.gam()] for more details.].
#' @param ... Currently not used.
#'
#' @return A vector of predictions, or a data frame containing the predictions, standard errors, and uncertainty intervals.
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
#' @importFrom mgcv predict.gam
#' @importFrom stats qnorm na.pass
#' @md

predict.gdam <- function(object,
                         newdata,
                         type = "response",
                         se.fit = FALSE,
                         coverage = 0.95,
                         na.action = na.pass,
                         ...) {

    if(!inherits(object, "gdam"))
        stop("object must be of class \"gdam\".")

    type <- match.arg(type, choices = c("link", "response", "lpmatrix"))

    get_MM <- mgcv::predict.gam(object = object$gamObject,
                                newdata = newdata,
                                type = "lpmatrix",
                                na.action = na.action)

    if(type == "lpmatrix")
        return(get_MM)

    point_prediction <- as.vector(get_MM %*% object$coefficients)
    names(point_prediction) <- rownames(get_MM)
    if(!se.fit)
        return(point_prediction)

    if(se.fit) {
        getCov <- object$covariance_matrix[1:ncol(get_MM), 1:ncol(get_MM)]
        getCov <- 0.5*(getCov + t(getCov))
        v <- rowSums((get_MM %*% getCov) * get_MM)
        names(v) <- rownames(get_MM)
        out <- data.frame(fit = point_prediction, se.fit = sqrt(v))
        #crit <- abs(qnorm(0.5*(1-coverage)))

        #' #' Construct simultaneous confidence intervals using the approach of (https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/)
        #' BUdiff <- mgcv::rmvn(n = nsim, mu = rep(0, ncol(get_MM)), V = getCov)
        #' simDev <- tcrossprod(get_MM, BUdiff)
        #' absDev <- abs(sweep(simDev, 1, out$se.fit, FUN = "/"))
        #' masd <- apply(absDev, 2, max)
        #' crit <- quantile(masd, prob = 0.95)
        #' rm(BUdiff, simDev, absDev, masd)
        crit <- qnorm(0.975)

        out$lower <- out$fit - crit * out$se.fit
        out$upper <- out$fit + crit * out$se.fit
        return(out)
        }
    }



