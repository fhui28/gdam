#' @title Fast robust additive models using gamma divergence.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function fits a robust additive model, assuming an identity link function and a normal distribution for the errors, using a gamma divergence adaption of the restricted maximum likelihood function (REML) as the loss function (hence the name gdam, which stands for gamma divergence additive models).
#'
#' The functions works by taking an object of class [mgcv::gam()] and uses at as a starting point for optimization of the gamma divergence, which itself is performed via a minorization-maximization (MM) algorithm.
#'
#'
#' @param gamObject An object of class [mgcv::gam()] fitted using the `mgcv` package. It is **strongly recommended that the additive model is fitted using `method = REML`.**
#' @param gamma_tuning A positive tuning parameter for the gamma divergence, and typically less than 0.5 for almost most additive models. The closer this value is to zero, the more the gamma divergence behaves like standard restricted maximum likelihood estimation i.e., the same as what is done in `mgcv`.
#' @param control A list of control parameters for the MM algorithm.
#' \itemize{
#' \item \code{tol} The tolerance for convergence of the MM algorithm. Default is 1e-6.
#' \item \code{max_iteration} The maximum number of iterations for the MM algorithm. Default is 1000.
#' \item \code{rinse_and_repeat} An integer value for the number of times to perform the entire MM algorithm, but with a refitted GAM using the final weights from the previous application of the MM algorithm. This tends to be useful for more complex smooths like tensor product smooths, where the MM algorithm can sometimes get stuck in local maxima. Default is 0.
#' }
#'
#' @details
#' See manuscript for more details (currently in preparation).
#'
#' @return An object of class `gdam` with the following components:
#' \itemize{
#' \item \code{call} The matched call.
#' \item \code{coefficients} The estimated (smoothing) coefficients of the additive model.
#' \item \code{sp} The estimated smoothing parameters of the additive model.
#' \item \code{sigma2} The estimated residual variance of the additive model.
#' \item \code{linear_predictors} The linear predictors of the additive model.
#' \item \code{fitted_values} The fitted values of the additive model.
#' \item \code{residuals} The raw residuals of the additive model.
#' \item \code{final_weights} The final weights from the MM algorithm. The weights sum to the total number of observations in the data. Weights comparably less than 1 indicate that the observation is down-weighted in the estimation process i.e., is considered more likely to be an outlying observation
#' \item \code{gamma_tuning} The gamma tuning parameter used in the estimation process.
#' \item \code{gamma_divergence} The value of the gamma divergence at convergence.
#' \item \code{gamObject} The original `gamObject` used to fit the additive model.
#' \item \code{covariance_matrix} The estimated covariance matrix of the parameters, including of the residual variance (but not the smoothing parameters).
#' \item \code{Hscore} The conditional Hyvarinen score (H-score) of the additive model. This is a measure that can be used to select the tuning parameter `gamma_tuning`, when comparing across difference `gdam` fits. The H-score calculation is based on the proposal of [Sugusawa and Yonekura, 2021](https://doi.org/10.3390/e23091147).
#' \item \code{aic} The conditional Akaike Information Criterion (AIC) of the additive model. This is calculated as \eqn{-2 \times \text{gamma_divergence} + 2 \times \text{edf}}, where `edf` is the effective degrees of freedom of the additive model. Effectively adapts the idea of  [Kurara, 2024](https://doi.org/10.1080/03610926.2022.2155788); see also [mgcv::logLik.gam()].
#' }
#'
#' @note
#' **Warning**
#'
#' Not all smooth specifications available in `mgcv` are supported let alone being properly tested e.g., we don't think smooths linked using the "id" argument, linear functional terms, or varying coefficient models as available using the "by" argument are supported; see [mgcv::gam.models()] for specification of GAMs.
#' Apologies in advance if you encounter any issues with these, but do email the author if encounter any issues or have a feature request.
#'
#' @author
#' Francis K.C. Hui <fhui28@gmail.com>
#'
#' @examples
#' library(mgcv)
#' library(tidyverse)
#' library(gratia)
#'
#' # Simulate a additive model and add some outliers.
#' set.seed(122024)
#' simdat <- gamSim(eg = 7,
#' n = 200,
#' dist = "normal",
#' scale = 0.5,
#' verbose = TRUE)
#'
#' create_outliers <- rbinom(nrow(simdat), size = 1, prob = 0.05)
#' # Could also the probability homogeneous e.g., vary with one or more covariate values.
#' table(create_outliers)
#' simdat$contaminated_y <- simdat$y
#' simdat$contaminated_y[create_outliers == 1] <- rnorm(sum(create_outliers == 1),
#' mean = simdat$f[create_outliers==1],
#' sd = 9)
#'
#'
#' # Fit the model using mgcv and then fit the robust additive model using gdam
#' fit_mgcv <- gam(contaminated_y ~ s(x0) + s(x1) + s(x2) + s(x3),
#' data = simdat,
#' family = gaussian(),
#' method = "REML")
#'
#' fit_gdam <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.1)
#' fit_gdam
#'
#' summary(fit_gdam)
#'
#' # Basic illustration of functions
#' residuals(object = fit_gdam)
#'
#' appraise(fit_gdam)
#'
#' ## Automated partial smooth plots
#' gdam::smooth_estimates(fit_gdam) %>%
#' gratia::draw() & theme_bw()
#'
#' ## A bit more manual but customized partial smooth plots -- all covariates
#' gdam::smooth_estimates(fit_gdam) %>%
#' gratia::add_confint() %>%
#' pivot_longer(x0:x3, values_to = "x", names_to = "covariate") %>%
#' filter(!is.na(x)) %>%
#' ggplot() +
#' geom_line(aes(x = x, y = .estimate), lwd = 1.2) +
#' geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x), alpha = 0.2) +
#' geom_rug(aes(x = x),
#' data = simdat %>% pivot_longer(x0:x3, values_to = "x", names_to = "covariate"),
#' sides = "b") +
#' facet_wrap(. ~ covariate, scale = "free", nrow = 2) +
#' theme_bw()
#'
#' ## A bit more manual but customized partial smooth plots -- one covariate
#' ## but adding partial residuals
#' simdat_expanded <- add_partial_residuals(data = simdat, model = fit_gdam, select = NULL)
#' gdam::smooth_estimates(fit_gdam) %>%
#' gratia::add_confint() %>%
#' filter(.smooth == s(x0)) %>%
#' ggplot() +
#' geom_point(aes(x = x0, y = `s(x0)`, color = weights), data = simdat_expanded, cex = 1.5) +
#' geom_line(aes(x = x0, y = .estimate), lwd = 1.2) +
#' geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x0), alpha = 0.2) +
#' scale_color_viridis_c() +
#' theme_bw()
#'
#'
#' predict(fit_gdam, se.fit = TRUE)
#'
#' # Simple check of the in-sample performance.
#' # Please see the associated manuscript for a more systematic evaluation of
#' # gdam including against other robust methods for additive models
#' cbind(simdat$f, fit_mgcv$fitted.values, fit_gdam$fitted_values) %>%
#' t %>%
#' dist()
#'
#'
#' \dontrun{
#' # Checking various selection criteria -- mainly for choosing the gamma tuning parameter
#' # Personally we prefer using the H-score for this as it is slightly more
#' # principled in its construction for gamma divergence
#' fit_gdam00 <- gdam(gamObject = fit_mgcv, gamma_tuning = 1e-6)
#' fit_gdam001 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.01)
#' fit_gdam01 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.1)
#' fit_gdam02 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.2)
#' fit_gdam03 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.3)
#' fit_gdam04 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.4)
#' fit_gdam05 <- gdam(gamObject = fit_mgcv, gamma_tuning = 0.5)
#'
#' fit_gdam00$Hscore
#' fit_gdam001$Hscore
#' fit_gdam01$Hscore
#' fit_gdam02$Hscore
#' fit_gdam03$Hscore
#' fit_gdam04$Hscore
#' fit_gdam05$Hscore
#' }
#'
#' @export gdam
#' @importFrom mgcv gam model.matrix.gam notExp notLog
#' @importFrom stats dnorm gaussian nlminb
#' @import Matrix
#' @md

gdam <- function(gamObject,
                 gamma_tuning = 0.3,
                 control = list(tol = 1e-6,
                                max_iteration = 1000,
                                rinse_and_repeat = 0)) {

    ##-----------------------------
    #' # Checks and balances
    ##-----------------------------
    if(length(gamObject$sp) == 0)
        stop("gdam is designed for robust generalized additive modeling, meaning at least one smooth term should be included in the model.")

    if(class(gamObject)[1] != "gam")
        stop("gamObject must be a object of class gam, fitted using the mgcv package.")
    if(gamObject$family$family[1] != "gaussian")
        stop("gamObject must be have assumed family is \"gaussian\".")
    if(gamma_tuning <= 0)
        stop("gamma_tuning should be a positive tuning parameter for the gamma divergence. The closer this value is to zero, the more the gamma divergence behaves like standard maximum likelihood estimation.")

    if(is.null(control$tol))
        control$tol <- 1e-6
    if(is.null(control$max_iteration))
        control$max_iteration <- 1000
    if(is.null(control$rinse_and_repeat))
        control$rinse_and_repeat <- 0

    MM <- mgcv::model.matrix.gam(gamObject)

    ##-----------------------------
    #' # Initial values
    ##-----------------------------
    cw_betas <- gamObject$coefficients
    cw_lambda <- gamObject$sp
    cw_sigma2  <- gamObject$sig2
    cw_gamma_divergence <- -Inf
    diff_gamma_divergence <- 1
    counter <- 0
    calc_weights <- dnorm(gamObject$y, mean = MM %*% cw_betas, sd = sqrt(cw_sigma2))^gamma_tuning
    calc_weights <- calc_weights/mean(calc_weights)


    ##-----------------------------
    #' # Run MM algorithm
    ##-----------------------------
    while(diff_gamma_divergence > control$tol & counter < control$max_iteration) {

        ##----------------------
        #' ### Construct minorization and its gradient, and perform maximization
        ##----------------------
        update_gdam_fn <- function(x, return_more = FALSE) {
            get_sigma2 <- notExp(x[1])
            get_lambda <- notExp(x[-1])

            P <- .get_bigS(gamObject = gamObject, smoothing_parameters = get_lambda)
            new_betas <- crossprod(MM, Diagonal(x = calc_weights)) %*% MM + P
            new_betas <- as.vector(solve(new_betas) %*% crossprod(MM, Diagonal(x = calc_weights)) %*% gamObject$y)
            new_eta <- as.vector(MM %*% new_betas) + gamObject$offset
            new_resid <- gamObject$y - new_eta

            out <- nrow(MM) * gamma_tuning * log(get_sigma2) / (2 * (1 + gamma_tuning)) - 0.5 * sum(calc_weights * new_resid^2) / get_sigma2 - 0.5 * as.vector(crossprod(new_betas, P / get_sigma2) %*% new_betas)
            e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
            out <- out - 0.5 * determinant((crossprod(MM) + P), logarithm = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps])) # Improved calculation of logdet(Sigma) where Sigma = I + Z P^{-1} Z^T, which avoids inverse of P; (https://en.wikipedia.org/wiki/Determinant#Sylvester's_determinant_theorem) part b and see (https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00749.x)
            out <- out - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) * log(get_sigma2)

            if(!return_more)
                return(as.vector(-out))
            if(return_more) {
                new_gamma_divergence <- nrow(MM)/gamma_tuning * log(mean(dnorm(gamObject$y, mean = MM %*% new_betas, sd = sqrt(get_sigma2))^gamma_tuning)) + nrow(MM) * (1 + 2*gamma_tuning) /(2*(1+gamma_tuning)) * log(get_sigma2)
                new_gamma_divergence <- new_gamma_divergence - 0.5 * as.vector(crossprod(new_betas, P/ get_sigma2) %*% new_betas)
                new_gamma_divergence <- new_gamma_divergence - 0.5 * determinant((crossprod(MM) + P), logarithm = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps]))
                new_gamma_divergence <- new_gamma_divergence - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) * log(get_sigma2)

                return(list(neg_minorizer = as.vector(out),
                            new_gamma_divergence = new_gamma_divergence,
                            betas = new_betas,
                            sigma2 = get_sigma2,
                            lambda = get_lambda))
            }
        }

        update_gdam_grad <- function(x) {
            get_sigma2 <- notExp(x[1])
            get_lambda <- notExp(x[-1])

            Slists <- .get_Slists(gamObject = gamObject)
            P <- .get_bigS(gamObject = gamObject, smoothing_parameters = get_lambda)
            e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
            ZtWZplusPinv <- solve(crossprod(MM, Diagonal(x = calc_weights)) %*% MM + P)
            ZtWy <- crossprod(MM, Diagonal(x = calc_weights)) %*% gamObject$y
            new_betas <- as.vector(ZtWZplusPinv %*% ZtWy)
            new_eta <- as.vector(MM %*% new_betas) + gamObject$offset
            new_resid <- gamObject$y - new_eta

            #' ## Differentiating all but the one log determinant term in the minorized gamma divergence
            dbeta_dsigma2lambda <- NULL
            for(j1 in 1:length(get_lambda))
                dbeta_dsigma2lambda <- cbind(dbeta_dsigma2lambda, - (ZtWZplusPinv %*% Slists[[j1]] %*% ZtWZplusPinv %*% ZtWy))

            out <- c(nrow(MM) * gamma_tuning / (2 * (1 + gamma_tuning) * get_sigma2), numeric(length(get_lambda)))
            out[1] <- out[1] + 0.5 / get_sigma2^2 * sum(calc_weights * new_resid^2) + 0.5 / get_sigma2^2 * as.vector(crossprod(new_betas, P) %*% new_betas) - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) / get_sigma2

            for(j1 in 1:length(get_lambda))
                out[j1+1] <- out[j1+1] + 1 / get_sigma2 * as.vector(dbeta_dsigma2lambda[,j1] %*% crossprod(MM, calc_weights * new_resid)) - as.vector(crossprod(P %*% new_betas, dbeta_dsigma2lambda[,j1]) / get_sigma2) - 0.5 * sum(diag(Slists[[j1]] %*% tcrossprod(new_betas))) / get_sigma2

            #' ## Differentiating log determinant term in the minorized gamma divergence
            #' Since P is not full rank, then apply a transformation as per Appendix B in (https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00749.x)
            e2 <- eigen(Reduce("+", sapply(Slists, function(x) { x/norm(x, type = "F") })), symmetric = TRUE)
            Stilde_lists <- sapply(Slists, function(x) crossprod(e2$vectors[, e2$values > .Machine$double.eps], x) %*% e2$vectors[, e2$values > .Machine$double.eps])
            P_transformed <- Reduce("+", lapply(1:length(get_lambda), function(x) get_lambda[x] * Stilde_lists[[x]]))
            P_transformedinv <- solve(P_transformed)
            rm(e2)

            out_logdet <- numeric(length(x))
            for(j1 in 1:length(get_lambda))
                out_logdet[j1+1] <- -0.5 * sum(diag(solve(crossprod(MM) + P) %*% Slists[[j1]])) + 0.5 * sum(diag(P_transformedinv %*% Stilde_lists[[j1]]))


            #' ## Combine all gradients together
            return(-as.vector(out + out_logdet) * .notExp_grad(x))
            }

        #update_gdam_fn(x = notLog(c(cw_sigma2, cw_lambda)))
        update_gdam <- nlminb(start = notLog(c(cw_sigma2, cw_lambda)),
                                objective = update_gdam_fn,
                                gradient = update_gdam_grad,
                                lower = -1e5,
                                upper = 1e5,
                                control = list(trace = 0))

        ##----------------------
        #' ### Finish update of parameters
        ##----------------------
        new_stuff <- update_gdam_fn(x = update_gdam$par, return_more = TRUE)
        cw_betas <- new_stuff$betas
        cw_lambda <- new_stuff$lambda
        cw_sigma2  <- new_stuff$sigma2
        diff_gamma_divergence <- new_stuff$new_gamma_divergence - cw_gamma_divergence
        cw_gamma_divergence <- new_stuff$new_gamma_divergence

        ##----------------------
        #' ### Recalculate weights and complete iteration
        ##----------------------
        calc_weights <- dnorm(gamObject$y, mean = MM %*% cw_betas, sd = sqrt(cw_sigma2))^gamma_tuning
        calc_weights <- calc_weights/mean(calc_weights)

        message("Iteration: ", counter, "\t Gamma divergence: ", round(as.vector(new_stuff$new_gamma_divergence), 6))
        counter <- counter + 1
        }

    ##-----------------------------
    #' # An optional refitting process where. based on the final weights from the fitted gdam above. Fit a new gam using mgcv and then rerun the gdam again.
    #' This can sometimes be useful for more complex multivariate smooths, to get things must closer to global maxima of the gamma divergence likelihood
    #' From empirical experience, we find that often one extra R&R is sufficient =D
    ##-----------------------------
    if(control$rinse_and_repeat > 0) {
        rinse_counter <- 1
        while(rinse_counter <= control$rinse_and_repeat) {
            message("Doing a rinse and repeat for the ", rinse_counter, "-th time...")

            new_gamObject_fn <- function(object, new_weight) {
                dofit <- mgcv::gam(object$formula,
                                   data = cbind(object$model, wts = new_weight),
                                   weights = wts,
                                   method = "REML",
                                   family = gaussian())
                return(dofit)
                }
            new_gamObject <- new_gamObject_fn(object = gamObject,
                                              new_weight = calc_weights)

            ##-----------------------------
            #' # Initial values
            ##-----------------------------
            cw_betas <- new_gamObject$coefficients
            cw_lambda <- new_gamObject$sp
            cw_sigma2  <- new_gamObject$sig2
            cw_gamma_divergence <- -Inf
            diff_gamma_divergence <- 1
            counter <- 0
            calc_weights <- dnorm(gamObject$y, mean = MM %*% cw_betas, sd = sqrt(cw_sigma2))^gamma_tuning
            calc_weights <- calc_weights/mean(calc_weights)

            ##-----------------------------
            #' # Run MM algorithm
            ##-----------------------------
            while(diff_gamma_divergence > control$tol & counter < control$max_iteration) {

                ##----------------------
                #' ### Construct minorization and its gradient, and perform maximization
                ##----------------------
                update_gdam_fn <- function(x, return_more = FALSE) {
                    get_sigma2 <- notExp(x[1])
                    get_lambda <- notExp(x[-1])

                    P <- .get_bigS(gamObject = gamObject, smoothing_parameters = get_lambda)
                    new_betas <- crossprod(MM, Diagonal(x = calc_weights)) %*% MM + P
                    new_betas <- as.vector(solve(new_betas) %*% crossprod(MM, Diagonal(x = calc_weights)) %*% gamObject$y)
                    new_eta <- as.vector(MM %*% new_betas) + gamObject$offset
                    new_resid <- gamObject$y - new_eta

                    out <- nrow(MM) * gamma_tuning * log(get_sigma2) / (2 * (1 + gamma_tuning)) - 0.5 * sum(calc_weights * new_resid^2) / get_sigma2 - 0.5 * as.vector(crossprod(new_betas, P / get_sigma2) %*% new_betas)
                    e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
                    out <- out - 0.5 * determinant((crossprod(MM) + P), logarithm = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps])) # Improved calculation of logdet(Sigma) where Sigma = I + Z P^{-1} Z^T, which avoids inverse of P; (https://en.wikipedia.org/wiki/Determinant#Sylvester's_determinant_theorem) part b and see (https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00749.x)
                    out <- out - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) * log(get_sigma2)

                    if(!return_more)
                        return(as.vector(-out))
                    if(return_more) {
                        new_gamma_divergence <- nrow(MM)/gamma_tuning * log(mean(dnorm(gamObject$y, mean = MM %*% new_betas, sd = sqrt(get_sigma2))^gamma_tuning)) + nrow(MM) * (1 + 2*gamma_tuning) /(2*(1+gamma_tuning)) * log(get_sigma2)
                        new_gamma_divergence <- new_gamma_divergence - 0.5 * as.vector(crossprod(new_betas, P/ get_sigma2) %*% new_betas)
                        new_gamma_divergence <- new_gamma_divergence - 0.5 * determinant((crossprod(MM) + P), logarithm = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps]))
                        new_gamma_divergence <- new_gamma_divergence - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) * log(get_sigma2)

                        return(list(neg_minorizer = as.vector(out),
                                    new_gamma_divergence = new_gamma_divergence,
                                    betas = new_betas,
                                    sigma2 = get_sigma2,
                                    lambda = get_lambda))
                        }
                    }

                update_gdam_grad <- function(x) {
                    get_sigma2 <- notExp(x[1])
                    get_lambda <- notExp(x[-1])

                    Slists <- .get_Slists(gamObject = gamObject)
                    P <- .get_bigS(gamObject = gamObject, smoothing_parameters = get_lambda)
                    e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
                    ZtWZplusPinv <- solve(crossprod(MM, Diagonal(x = calc_weights)) %*% MM + P)
                    ZtWy <- crossprod(MM, Diagonal(x = calc_weights)) %*% gamObject$y
                    new_betas <- as.vector(ZtWZplusPinv %*% ZtWy)
                    new_eta <- as.vector(MM %*% new_betas) + gamObject$offset
                    new_resid <- gamObject$y - new_eta

                    #' ## Differentiating all but the one log determinant term in the minorized gamma divergence
                    dbeta_dsigma2lambda <- NULL
                    for(j1 in 1:length(get_lambda))
                        dbeta_dsigma2lambda <- cbind(dbeta_dsigma2lambda, - (ZtWZplusPinv %*% Slists[[j1]] %*% ZtWZplusPinv %*% ZtWy))

                    out <- c(nrow(MM) * gamma_tuning / (2 * (1 + gamma_tuning) * get_sigma2), numeric(length(get_lambda)))
                    out[1] <- out[1] + 0.5 / get_sigma2^2 * sum(calc_weights * new_resid^2) + 0.5 / get_sigma2^2 * as.vector(crossprod(new_betas, P) %*% new_betas) - 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) / get_sigma2

                    for(j1 in 1:length(get_lambda))
                        out[j1+1] <- out[j1+1] + 1 / get_sigma2 * as.vector(dbeta_dsigma2lambda[,j1] %*% crossprod(MM, calc_weights * new_resid)) - as.vector(crossprod(P %*% new_betas, dbeta_dsigma2lambda[,j1]) / get_sigma2) - 0.5 * sum(diag(Slists[[j1]] %*% tcrossprod(new_betas))) / get_sigma2

                    #' ## Differentiating log determinant term in the minorized gamma divergence
                    #' Since P is not full rank, then apply a transformation as per Appendix B in (https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00749.x)
                    e2 <- eigen(Reduce("+", sapply(Slists, function(x) { x/norm(x, type = "F") })), symmetric = TRUE)
                    Stilde_lists <- sapply(Slists, function(x) crossprod(e2$vectors[, e2$values > .Machine$double.eps], x) %*% e2$vectors[, e2$values > .Machine$double.eps])
                    P_transformed <- Reduce("+", lapply(1:length(get_lambda), function(x) get_lambda[x] * Stilde_lists[[x]]))
                    P_transformedinv <- solve(P_transformed)
                    rm(e2)

                    out_logdet <- numeric(length(x))
                    for(j1 in 1:length(get_lambda))
                        out_logdet[j1+1] <- -0.5 * sum(diag(solve(crossprod(MM) + P) %*% Slists[[j1]])) + 0.5 * sum(diag(P_transformedinv %*% Stilde_lists[[j1]]))


                    #' ## Combine all gradients together
                    return(-as.vector(out + out_logdet) * .notExp_grad(x))
                    }

                #update_gdam_fn(x = notLog(c(cw_sigma2, cw_lambda)))
                update_gdam <- nlminb(start = notLog(c(cw_sigma2, cw_lambda)),
                                        objective = update_gdam_fn,
                                        gradient = update_gdam_grad,
                                        lower = -1e5,
                                        upper = 1e5,
                                        control = list(trace = 0))

                ##----------------------
                #' ### Finish update of parameters
                ##----------------------
                new_stuff <- update_gdam_fn(x = update_gdam$par, return_more = TRUE)
                cw_betas <- new_stuff$betas
                cw_lambda <- new_stuff$lambda
                cw_sigma2  <- new_stuff$sigma2
                diff_gamma_divergence <- new_stuff$new_gamma_divergence - cw_gamma_divergence
                cw_gamma_divergence <- new_stuff$new_gamma_divergence

                ##----------------------
                #' ### Recalculate weights and complete iteration
                ##----------------------
                calc_weights <- dnorm(gamObject$y, mean = MM %*% cw_betas, sd = sqrt(cw_sigma2))^gamma_tuning
                calc_weights <- calc_weights/mean(calc_weights)

                message("Iteration: ", counter, "\t Gamma divergence: ", round(as.vector(new_stuff$new_gamma_divergence), 6))
                counter <- counter + 1
                }

            rinse_counter <- rinse_counter + 1
            }
        }

    ##-----------------------------
    #' # Complete fitting and attempt covariance matrix calculations
    ##-----------------------------
    out <- list(call = match.call(),
                coefficients = cw_betas,
                sp = cw_lambda,
                sigma2 = cw_sigma2,
                linear_predictors = as.vector(MM %*% cw_betas) + gamObject$offset,
                fitted_values = as.vector(MM %*% cw_betas) + gamObject$offset,
                residuals = gamObject$y - as.vector(MM %*% cw_betas) - gamObject$offset,
                final_weights = calc_weights,
                gamma_tuning = gamma_tuning,
                gamma_divergence = as.vector(unlist(cw_gamma_divergence)),
                gamObject = gamObject)
    names(out$coefficients) <- names(gamObject$coefficients)


    #' ## Hessian of minorizing function
    P <- .get_bigS(gamObject = gamObject, smoothing_parameters = cw_lambda)
    e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
    Hessminorizer_11 <- -1/cw_sigma2 * (crossprod(MM * sqrt(calc_weights)) + P)
    Hessminorizer_12 <- -1/cw_sigma2^2 * crossprod(MM*calc_weights, out$residuals) + 1/cw_sigma2^2 * (P %*% cw_betas)
    Hessminorizer_22 <- -1/cw_sigma2^3 * sum(calc_weights * out$residuals^2) - nrow(MM) * gamma_tuning/(2 *(1 + gamma_tuning) * cw_sigma2^2) - 1/cw_sigma2^3 * crossprod(cw_betas, P) %*% cw_betas + 0.5 * (nrow(MM) - sum(e < .Machine$double.eps)) / cw_sigma2^2
    Hessminorizer <- rbind(cbind(Hessminorizer_11, Hessminorizer_12),
                           cbind(t(Hessminorizer_12), Hessminorizer_22))
    rm(Hessminorizer_11, Hessminorizer_12, Hessminorizer_22)

    #' ## Derivative of weights wrt to parameters
    r_i <- cbind(1/cw_sigma2 * (MM * out$residuals), -0.5/cw_sigma2 + 0.5/cw_sigma2^2 * out$residuals^2)
    calc_weights_derivative <- calc_weights * gamma_tuning * (r_i - mean(calc_weights * r_i))
    rm(r_i)

    #' ## Cross derivative of minorizing function wrt to parameters from current and previous iteration
    #' Based on testing with the Hessian derived using numerical differentiation, construct the cross derivative of the minorizer in a more manual way (the last column of CrossDerivminorizer_1 appears to be more correct than the values of CrossDerivminorizer_2. Maybe double check the maths here?
    CrossDerivminorizer_1 <- 1/cw_sigma2 * crossprod(MM * out$residuals, calc_weights_derivative) # 1/cw_sigma2 * (P %*% cw_betas)
    CrossDerivminorizer_2 <- 0.5/cw_sigma2^2 * crossprod(out$residuals^2, calc_weights_derivative) #- nrow(MM)/(2 *(1 + gamma_tuning) * cw_sigma2) + 0.5/cw_sigma2^2 * crossprod(cw_betas, P) %*% cw_betas
    CrossDerivminorizer <- rbind(CrossDerivminorizer_1, rep(0, ncol(CrossDerivminorizer_1)))
    CrossDerivminorizer[nrow(CrossDerivminorizer),] <- c(CrossDerivminorizer_1[,ncol(CrossDerivminorizer_1)], CrossDerivminorizer_2[length(CrossDerivminorizer_2)])
    rm(CrossDerivminorizer_1, CrossDerivminorizer_2)
    CrossDerivminorizer <- 0.5*(CrossDerivminorizer + t(CrossDerivminorizer)) #' Symmetrize

    negHess_gamma_divergence <- -Hessminorizer - CrossDerivminorizer

    out$covariance_matrix <- as.matrix(solve(negHess_gamma_divergence))
    rownames(out$covariance_matrix) <- c(names(gamObject$coefficients), "sigma2")
    colnames(out$covariance_matrix) <- c(names(gamObject$coefficients), "sigma2")

    ##-----------------------------
    #' # Selection criterion
    ##-----------------------------
    class(out) <- "gdam"
    out$Hscore <- .calc_Hscore(object = out)
    s <- mgcv::gam(out$gamObject$formula,
                   data = cbind(out$gamObject$model, wts = out$final_weights),
                   weights = wts,
                   sp = out$sp,
                   scale = out$sigma2,
                   method = "REML")
    out$aic <- (-2*out$gamma_divergence + 2*sum(s$edf)) #' AIC: adopts the ideas of (https://doi.org/10.1080/03610926.2022.2155788)
    # {
    #     sink("/dev/null")
    #     s <- invisible(summary.gdam(object = out))
    #     out$aic <- (-2*out$gamma_divergence + 2*sum(s$edf)) #' AIC: adopts the ideas of (https://doi.org/10.1080/03610926.2022.2155788)
    #     sink(NULL)
    # }

    return(out)
    }



#' @title Hidden function to compute the conditional Hyvarinen score (H-score)
#' @description Adapts the idea of (https://doi.org/10.3390/e23091147) to the gamma divergence of the conditional log-likelihood of an additive model
#' @noRd
#' @noMd

.calc_Hscore <- function(object) {
    Cgammasigma2 <- ((1 + object$gamma_tuning)^(-0.5) * (2*pi*object$sigma2)^(-0.5*object$gamma_tuning))^(object$gamma_tuning / (1+object$gamma_tuning))
    Hscore1 <- 2*(object$gamma_tuning * object$residuals^2 - object$sigma2) / (object$sigma2^2 * Cgammasigma2) * dnorm(object$gamObject$y, mean = object$linear_predictors, sd = sqrt(object$sigma2))^object$gamma_tuning
    Hscore2 <- object$residuals^2 / (object$sigma2^2 * Cgammasigma2^2) * dnorm(object$gamObject$y, mean = object$linear_predictors, sd = sqrt(object$sigma2))^(2*object$gamma_tuning)

    out <- sum(Hscore1 + Hscore2)
    return(out)
    }


# gamdiv_fn <- function(x) {
#     new_betas <- x[1:ncol(MM)]
#     get_sigma2 <- (x[1 + ncol(MM)])
#     get_lambda <- out$sp
#
#     P <- .get_bigS(gamObject = gamObject, smoothing_parameters = get_lambda)
#     new_eta <- as.vector(MM %*% new_betas) + gamObject$offset
#     new_resid <- gamObject$y - new_eta
#
#     out <- - nrow(MM) * log(get_sigma2) / (2 * (1 + gamma_tuning)) - 0.5 * sum(calc_weights * new_resid^2) / get_sigma2 - 0.5 * as.vector(crossprod(new_betas, P / get_sigma2) %*% new_betas)
#     e <- eigen(P, only.values = TRUE, symmetric = TRUE)$values
#     out <- out - 0.5 * determinant((crossprod(MM) + P), log = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps])) # Improved calculation of logdet(Sigma) where Sigma = I + Z P^{-1} Z^T, which avoids inverse of P; see (https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2010.00749.x)
#
#     new_gamma_divergence <- nrow(MM)/gamma_tuning * log(mean(dnorm(gamObject$y, mean = MM %*% new_betas, sd = sqrt(get_sigma2))^gamma_tuning)) + nrow(MM) * gamma_tuning /(2*(1+gamma_tuning)) * log(get_sigma2)
#     new_gamma_divergence <- new_gamma_divergence - 0.5 * as.vector(crossprod(new_betas, P/ get_sigma2) %*% new_betas)
#     new_gamma_divergence <- new_gamma_divergence - 0.5 * determinant((crossprod(MM) + P), log = TRUE)$mod + 0.5 * sum(log(e[e > .Machine$double.eps]))
#
#     return(new_gamma_divergence)
#     }
#
#
# gamdiv_fn(x = c(out$coefficients, cw_sigma2))
#
# goldstandardHess <- numDeriv::hessian(gamdiv_fn, x = c(out$coefficients, cw_sigma2))
#
# goldstandard_covariancematrix <- solve(-goldstandardHess)
#
#
# plot(goldstandard_covariancematrix %>% diag, out$covariance_matrix %>% diag, log = "xy")
# abline(0,1)


