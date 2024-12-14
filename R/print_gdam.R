#' @title Print a fitted gdam model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function prints out some details of the fitted `gdam` object in a human-readable format, similar to that of [mgcv::print.gam()]. This includes the tuning parameter used, the formula, the estimated degrees of freedom, and the gamma divergence. Note the **estimated degrees of freedom should be taken with a grain of salt**, as they are estimated from a post-hoc GAM fit (we conjecture that they are likely to be marginally underestimated).
#'
#' @param x A fitted `gdam` object.
#' @param ... Currently not used.
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
#' @importFrom mgcv gam
#' @method print gdam
#' @md

print.gdam <- function(x, ...) {
    if(!inherits(x, "gdam"))
        stop("model must be of class \"gdam\".")

    posthoc_mgcv_fit_fn <- function(object, new_weight) {
        dofit <- mgcv::gam(object$gamObject$formula,
                           data = cbind(object$gamObject$model, wts = new_weight),
                           weights = wts,
                           sp = object$sp,
                           scale = object$sigma2,
                           method = "REML")
        return(dofit)
        }
    posthoc_mgcv_fit <- posthoc_mgcv_fit_fn(object = x,
                                            new_weight = x$final_weights)

    cat("Additive model fitting using gamma divergence with gamma =", x$gamma_tuning, "\n\n")

    cat("Formula:\n")
    if(is.list(posthoc_mgcv_fit$formula))
        for(i in 1:length(posthoc_mgcv_fit$formula)) print(posthoc_mgcv_fit$formula[[i]])
    else print(posthoc_mgcv_fit$formula)

    n.smooth <- length(posthoc_mgcv_fit$smooth)
    if (n.smooth == 0)
        cat("Total model degrees of freedom", sum(posthoc_mgcv_fit$edf), "\n")
    else {
        edf <- 0
        cat("\nEstimated degrees of freedom (estimated from a post-hoc GAM; please take results with a gain of salt!):\n")
        for (i in 1:n.smooth) edf[i] <- sum(posthoc_mgcv_fit$edf[posthoc_mgcv_fit$smooth[[i]]$first.para:posthoc_mgcv_fit$smooth[[i]]$last.para])
        edf.str <- format(round(edf, digits = 4), digits = 3,
                          scientific = FALSE)
        for (i in 1:n.smooth) {
            cat(edf.str[i], " ", sep = "")
            if (i%%7 == 0)
                cat("\n")
            }
        cat(" total =", round(sum(posthoc_mgcv_fit$edf), digits = 2), "\n")
        }

    cat("\n")
    cat("Gamma divergence: ", as.numeric(x$gamma_divergence))
    cat("\n")
    invisible(x)
    }
