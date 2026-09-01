#' @title Summary for a fitted gdam model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted `gdam` object and produces some useful *post-hoc* summaries of the fitted model. Note all results **should be taken with a grain of salt**, as they are estimated from a post-hoc GAM fit using [mgcv::gam()] (we conjecture quantities like degrees of freedom are likely to be marginally underestimated, and hence P-values are likely to be marginally smaller than they should be). Please see [mgcv::summary.gam()] for details on the output.
#'
#' @param object A fitted `gdam` object.
#' @param digits The number of significant digits to print.
#' @param signif.stars If `TRUE`, `p`-values are printed with stars indicating significance.
#' @param ... Currently not used.
#'
#' @return A summary of the fitted `gdam` object, similar to that of [mgcv::summary.gam()]. Please see the help file for that function for more details.
#'
#' @note
#' **Warning**
#'
#' The p-values are approximate and neglect smoothing parameter uncertainty. They are likely to be somewhat too low when smoothing parameter estimates are highly uncertain; please see the [mgcv::summary.gam()] help file for more details.
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
#' @importFrom mgcv gam summary.gam
#' @importFrom stats printCoefmat
#' @method summary gdam
#' @md

summary.gdam <- function(object,
                         digits = max(3, getOption("digits") - 3),
                         signif.stars = getOption("show.signif.stars"),
                         ...) {
    if(!inherits(object, "gdam"))
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
    posthoc_mgcv_fit <- posthoc_mgcv_fit_fn(object = object,
                                            new_weight = object$final_weights)
    posthoc_mgcv_fit$Vp <- object$covariance_matrix[-nrow(object$covariance_matrix),-nrow(object$covariance_matrix)]

    cat("Additive model fitting using gamma divergence with gamma = ", object$gamma_tuning, "\n")
    cat("Post-hoc GAM summary: **Please take all results below with a gain of salt!**\n\n")
    cat("#--------------------------#\n")
    .print2(x = summary(posthoc_mgcv_fit),
            digits = digits,
            signif.stars = signif.stars)
    cat("#--------------------------#")
    }



#' @noRd
#' @noMd
.print2 <- function(x,
                    digits,
                    signif.stars,
                    ...) {
    cat("Formula:\n")

    if (is.list(x$formula))
        for (i in 1:length(x$formula)) print(x$formula[[i]]) else print(x$formula)

    if (length(x$p.coeff) > 0) {
        cat("\nParametric coefficients:\n")
        stats::printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
        }

    cat("\n")
    if(x$m > 0) {
        cat("Approximate significance of smooth terms:\n")
        stats::printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
        }

    cat("\n")

    #if(!is.null(x$rank) && x$rank< x$np) cat("Rank: ",x$rank,"/",x$np,"\n",sep="")
    #if(!is.null(x$r.sq)) cat("R-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5),"  ")
    #if (length(x$dev.expl)>0) cat("Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%",sep="")
    # cat("\n")
    # if !is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))
    #     cat(x$method," = ",formatC(x$sp.criterion,digits=5),sep="")

    cat("Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
    invisible(x)
    }
