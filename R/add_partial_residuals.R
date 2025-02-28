#' @title Add partial residuals from a fitted gdam model
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is basically a wrapper around the [gratia::add_partial_residuals()] function for adding partial residuals, but designed to work with `gdam` objects. Please see the associated help file in the `gratia` package for more details.
#'
#' @param data A data frame containing values for the variables used to fit the `gdam` model.
#' @param model A fitted `gdam` object.
#' @param select Character, logical, or numeric; which smooths to plot; see [gratia::add_partial_residuals()] for details.
#' @param partial_match Should smooths be selected by partial matches with select? If `TRUE`, select can only be a single string to match against.
#' @param ... Arguments passed to other methods.
#'
#' @return
#' A data frame of the same dimension as `data` with partial residuals added on, while the final column also contains the final weights from the MM algorithm used to fit the `gdam` model.
#'
#' @note
#' Acknowledgments to Gavin Simpson and the `gratia` package for the original [gratia::add_partial_residuals()] function.
#'
#' @author
#' Francis K.C. Hui <fhui28@gmail.com>
#'
#' @export
#' @rdname add_partial_residuals
#' @importFrom gratia add_partial_residuals smooths
#' @md
add_partial_residuals.gdam <- function(data,
                                       model,
                                       select = NULL,
                                       partial_match = FALSE,
                                       ...) {

    if(!inherits(model, "gdam"))
        stop("model must be of class \"gdam\".")

    sms <- gratia::smooths(model$gamObject)
    ## which were selected; select = NULL -> all selected
    take <- gratia:::check_user_select_smooths(sms,
                                              select = select,
                                              partial_match = partial_match)
    if (!any(take)) {
        stop("No smooth label matched 'select'. Try 'partial_match = TRUE'?", call. = FALSE)
        }
    sms <- sms[take] # subset to selected smooths


    ## compute partial resids
    p_resids <- .compute_partial_residuals2(model,
                                            terms = sms,
                                            data = data)

    ## bind partial residuals to data
    data <- cbind(data, as.data.frame(p_resids), data.frame(weights = model$final_weights))
    return(data)
    }



#' @title Hidden function for computing partial residuals -- acknowledgements goes to Gavin Simpson for this!
#' @noRd
#' @noMd
.compute_partial_residuals2 <- function(object,
                                        terms = NULL,
                                        data = NULL) {

    w <- object$gamObject$weights
    ## need as.numeric for gamm() objects
    w_resid <- as.numeric(residuals.gdam(object, "working")) * sqrt(w)

    ## if data is null, just grab the $model out of object
    if (is.null(data)) {
        data <- object$gamObject[["model"]]
        }
    else {
        ## check size of data
        if (nrow(data) != length(w_resid)) {
            stop("Length of model residuals not equal to number of rows in 'data'", call. = FALSE)
        }
    }

    ## get the contributions for each selected smooth
    new_gam <- object$gamObject
    new_gam$sp <- object$sp
    new_gam$scale <- object$sigma2
    new_gam$sig2 <- object$sigma2
    new_gam$fitted.values <- object$fitted_values
    new_gam$linear.predictors <- object$linear_predictor
    new_gam$Vp <- object$covariance_matrix

    p_terms <- if (is.null(terms)) {
        mgcv::predict.gam(new_gam,
                          type = "terms",
                          newdata = data)
        }
    else {
        mgcv::predict.gam(new_gam,
                      type = "terms",
                      terms = terms,
                      newdata = data)
        }
    attr(p_terms, "constant") <- NULL # remove intercept attribute

    ## and compute partial residuals
    p_resids <- p_terms + w_resid

    return(as.data.frame(p_resids))
    }
