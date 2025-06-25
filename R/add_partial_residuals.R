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
    take <- .check_user_select_smooths2(sms,
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



#' @title Hidden function for Select smooths based on user's choices -- acknowledgements goes to Gavin Simpson for this!
#' @noRd
#' @noMd
.check_user_select_smooths2 <- function(smooths, select = NULL,
                                       partial_match = FALSE,
                                       model_name = NULL) {
    lenSmo <- length(smooths)
    select <- if (!is.null(select)) {
        lenSel <- length(select)
        if (is.numeric(select)) {
            if (lenSmo < lenSel) {
                stop("Trying to select more smooths than are in the model.")
                }
            if (any(select > lenSmo)) {
                stop("One or more indices in 'select' > than the number of smooths in the model.")
                }
            l <- rep(FALSE, lenSmo)
            l[select] <- TRUE
            l
            }
        else if (is.character(select)) {
            take <- if (isTRUE(partial_match)) {
                if (length(select) != 1L) {
                    stop("When 'partial_match' is 'TRUE', 'select' must be a single string")
                    }
                grepl(select, smooths, fixed = TRUE)
                } else {
                smooths %in% select
                }
            # did we fail to match?
            if (sum(take) < length(select)) {
                # must have failed to match at least one of `smooth`
                if (all(!take)) {
                    stop("Failed to match any smooths in model",
                         ifelse(is.null(model_name), "",
                                paste0(" ", model_name)
                         ),
                         ".\nTry with 'partial_match = TRUE'?",
                         call. = FALSE
                    )
                    } else {
                    stop("Some smooths in 'select' were not found in model ",
                         ifelse(is.null(model_name), "", model_name),
                         ":\n\t",
                         paste(select[!select %in% smooths], collapse = ", "),
                         call. = FALSE)
                    }
                }
            take
            }
        else if (is.logical(select)) {
            if (lenSmo != lenSel) {
                stop("When 'select' is a logical vector, 'length(select)' must equal\nthe number of smooths in the model.")
            }
            select
            } else {
            stop("'select' is not numeric, character, or logical.")
            }
        }
    else {
        rep(TRUE, lenSmo)
        }

    select
    }
