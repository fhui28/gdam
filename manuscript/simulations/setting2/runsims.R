##-----------------------
#' # Template code for replicating simulation setting II in the main paper
#' NOTE: The code below assumes the simulations are run on a machine that has multiple cores available for parallelizing the simulations.
##-----------------------
rm(list = ls())
library(tidyverse)
library(colorspace)
library(Matrix)
library(mgcv)
library(DoubleRobGam) #remotes::install_github("ilapros/DoubleRobGam")
library(GJRM)
here::i_am("setting2/runsims.R")
library(here)

source(here("Robust-and-efficient-estimation-of-nonparametric-GLMs/Main function.R")) #https://github.com/ioanniskalogridis/Robust-and-efficient-estimation-of-nonparametric-GLMs
source(here("utils.R"))

library(gdam) #devtools::install_github("fhui28/gdam")

library(foreach)
library(doParallel)
registerDoParallel(cores = detectCores() - 3)


##---------------
#' #' Run simulations
##---------------
sim_fn <- function(seed = NULL,
                   N = 200,
                   prob_outlier = 0.05,
                   k = 30, # Through personal communication with Prof. Simon Wood, it was suggested that a larger number of basis functions than the default number (around 9 to 10) was needed to overcome the approximation bias of the penalized splines to the true smoothing function. After some experimentation, it was found that k = 30 was a good choice for this.
                   true_sigma2 = 0.5) {
    message("Onto seed ", seed)

    ##-----------------------------
    #' # Simulate data and add outliers
    ##-----------------------------
    set.seed(seed)
    simdat <- gamSim(eg = 7,
                     n = N,
                     dist = "normal",
                     scale = sqrt(true_sigma2),
                     verbose = FALSE)

    create_outliers <- rbinom(nrow(simdat), size = 1, prob = prob_outlier)
    table(create_outliers)
    simdat$contaminated_y <- simdat$y
    simdat$contaminated_y[create_outliers == 1] <- rnorm(sum(create_outliers == 1), mean = simdat$f[create_outliers==1], sd = 9)

    simdat_test <- gamSim(eg = 7,
                          n = 1000,
                          dist = "normal",
                          scale = sqrt(true_sigma2),
                          verbose = FALSE)


    ##-----------------------------
    #' # Fit existing methods to the literature, along with new gamma divergence
    #' Robust local polynomial estimator of Azadeh and Salibian-Barrera (2011), denoted by RGAM and in principle available in the robustgam/rgam packages, is *not* available since out of the box it can not do Gaussian responses.
    ##-----------------------------
    #' ## Standard non-robust GAMs using mgcv
    message("Fitting mgcv GAM")
    tic <- proc.time()
    fit_mgcv <- gam(contaminated_y ~ s(x0, k = k) + s(x1, k = k) + s(x2, k = k) + s(x3, k = k),
                data = simdat,
                select = FALSE,
                family = gaussian(),
                method = "REML")
    toc <- proc.time()
    fit_mgcv$time_taken <- toc - tic

    #' ## Implements the monotone Huber bounded influence function method of (https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1541-0420.2011.01630.x) coupled with P-spline approach on a quasi-likelihood function
    message("Fitting DoubleRobGam")
    tic <- proc.time()
    choose_knots <- k
    make_formula <- paste0("contaminated_y ~ bsp(x0, nknots = ",  choose_knots, ") + bsp(x1, nknots = ", choose_knots, ") + bsp(x2, nknots = ", choose_knots, ") + bsp(x3, nknots = ", choose_knots, ")") %>%
        as.formula()
    fit_RobGam <- DoubleRobGam(formulaM = make_formula,
                                data = simdat,
                                family = "gaussian")
    toc <- proc.time()
    fit_RobGam$time_taken <- toc - tic


    #' ## DPD method of (https://doi.org/10.1007/s11749-023-00866-x). Note they did not set up the method and code for more than one covariate, and so was not actually implemented in this setting!

    #' ## The Psi-divergence approach of (https://link.springer.com/article/10.1007/s11222-020-09979-x), which wraps log-likelihood contributions with a smooth convex function that down-weights small log-likelihood contributions
    message("Fitting Psi divergence")
    tic <- proc.time()
    fl <- list(contaminated_y ~ s(x0, k = k) + s(x1, k = k) + s(x2, k = k) + s(x3, k = k))
    fit_psidiv <- try(GJRM::gamlss(fl,
                                   data = simdat,
                                   family = "N",
                                   robust = TRUE,
                                   sp.method = "efs"),
                      silent = TRUE)
    if(inherits(fit_psidiv, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    rm(fl)
    toc <- proc.time()
    fit_psidiv$time_taken <- toc - tic


    #' ## GDAMs
    message("Fitting gdams")
    k <<- k
    tic <- proc.time()
    fit_gdam0 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 1e-6), silent = TRUE)
    if(inherits(fit_gdam0, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam0$time_taken <- toc - tic

    fit_gdam001 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.01), silent = TRUE)
    if(inherits(fit_gdam001, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam001$time_taken <- toc - tic

    tic <- proc.time()
    fit_gdam01 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.1), silent = TRUE)
    if(inherits(fit_gdam01, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam01$time_taken <- toc - tic

    tic <- proc.time()
    fit_gdam02 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.2), silent = TRUE)
    if(inherits(fit_gdam02, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam02$time_taken <- toc - tic

    tic <- proc.time()
    fit_gdam03 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.3), silent = TRUE)
    if(inherits(fit_gdam03, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam03$time_taken <- toc - tic

    tic <- proc.time()
    fit_gdam04 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.4), silent = TRUE)
    if(inherits(fit_gdam04, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam04$time_taken <- toc - tic

    tic <- proc.time()
    fit_gdam05 <- try(gdam(gamObject = fit_mgcv, gamma_tuning = 0.5), silent = TRUE)
    if(inherits(fit_gdam05, "try-error")) {
        return(list(ME = NA, interval_performance = NA, simdat_test = NA, time_taken = NA))
        }
    toc <- proc.time()
    fit_gdam05$time_taken <- toc - tic

    select_hscore <- which.min(c(fit_gdam0$Hscore, fit_gdam001$Hscore, fit_gdam01$Hscore, fit_gdam02$Hscore, fit_gdam03$Hscore, fit_gdam04$Hscore, fit_gdam05$Hscore))
    fit_gdamhscore <- list(fit_gdam0, fit_gdam001, fit_gdam01, fit_gdam02, fit_gdam03, fit_gdam04, fit_gdam05)[[select_hscore]]
    fit_gdamhscore$time_taken <- fit_gdam0$time_taken + fit_gdam001$time_taken + fit_gdam01$time_taken + fit_gdam02$time_taken + fit_gdam03$time_taken + fit_gdam04$time_taken + fit_gdam05$time_taken
    rm(select_hscore)

    ##-----------------------------
    #' # In and out of sample performance
    ##-----------------------------
    message("Evaluating performance...")

    all_time_taken <- c(fit_mgcv$time_taken[3],
                        fit_RobGam$time_taken[3],
                        fit_psidiv$time_taken[3],
                        fit_gdam01$time_taken[3],
                        fit_gdam05$time_taken[3],
                        fit_gdamhscore$time_taken[3]
                        )
    names(all_time_taken) <- c("gamfit", "RobGam", "psidivergence", "gdam01", "gdam05", "gdamhscore")


    mgcv_predict_MM <- predict(fit_mgcv, newdata = simdat_test, type = "lpmatrix")
    simdat_test <- simdat_test %>%
        mutate(gamfit = predict(fit_mgcv, newdata = simdat_test, se = FALSE, type = "response"),
               RobGam = predict_RobGam(fit_RobGam, newdata = simdat_test %>% dplyr::select(x0:x3), choose_knots = choose_knots),
               psidivergence = as.vector(mgcv_predict_MM %*% fit_psidiv$coefficients[1:length(fit_mgcv$coefficients)]),
               gdam01 = predict(fit_gdam01, newdata = simdat_test, se.fit = FALSE),
               gdam05 = predict(fit_gdam05, newdata = simdat_test, se.fit = FALSE),
               gdamhscore = predict(fit_gdamhscore, newdata = simdat_test, se.fit = FALSE)
               )


    ME <- simdat_test %>%
        mutate(across(gamfit:last_col(), ~ (.x - f)^2)) %>%
        dplyr::select(gamfit:last_col()) %>%
        colMeans


    ##-----------------------------
    #' # Coverage probability
    ##-----------------------------
    base_predict <- predict(fit_mgcv, newdata = simdat_test, se.fit = TRUE, unconditional = FALSE)
    mgcv_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                         lower = base_predict$fit - qnorm(0.975)*base_predict$se.fit,
                                         upper = base_predict$fit + qnorm(0.975)*base_predict$se.fit,
                                         method = "mgcv")

    RobGam_testMM <- cbind(1,
                           bsp(simdat_test$x0 %>% na.omit(), nknots = choose_knots)$B,
                           bsp(simdat_test$x1 %>% na.omit(), nknots = choose_knots)$B,
                           bsp(simdat_test$x2 %>% na.omit(), nknots = choose_knots)$B,
                           bsp(simdat_test$x3 %>% na.omit(), nknots = choose_knots)$B)
    RobGam_predict_intervals <- data.frame(fit = as.vector(RobGam_testMM %*% fit_RobGam$coefficients),
                                           se.fit = sqrt(diag(RobGam_testMM %*% tcrossprod(fit_RobGam$cov.coef, RobGam_testMM))))
    RobGam_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                           lower = RobGam_predict_intervals$fit - qnorm(0.975)*RobGam_predict_intervals$se.fit,
                                           upper = RobGam_predict_intervals$fit + qnorm(0.975)*RobGam_predict_intervals$se.fit,
                                           method = "RobGam")


    psidiv_testMM <- predict(fit_mgcv, newdata = simdat_test, type = "lpmatrix")
    psidiv_predict_intervals <- data.frame(fit = as.vector(psidiv_testMM %*% fit_psidiv$coefficients[1:length(fit_mgcv$coefficients)]),
                                           se.fit = sqrt(diag(psidiv_testMM %*% tcrossprod(fit_psidiv$Vb[1:length(fit_mgcv$coefficients),1:length(fit_mgcv$coefficients)], psidiv_testMM))))
    psidiv_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                           lower = psidiv_predict_intervals$fit - qnorm(0.975)*psidiv_predict_intervals$se.fit,
                                           upper = psidiv_predict_intervals$fit + qnorm(0.975)*psidiv_predict_intervals$se.fit,
                                           method = "psidiv")


    gdam01_predict_intervals <- predict(fit_gdam01, newdata = simdat_test, se.fit = TRUE)
    gdam01_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                             lower = gdam01_predict_intervals$lower,
                                             upper = gdam01_predict_intervals$upper,
                                             method = "gdam01")

    gdam05_predict_intervals <- predict(fit_gdam05, newdata = simdat_test, se.fit = TRUE)
    gdam05_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                             lower = gdam05_predict_intervals$lower,
                                             upper = gdam05_predict_intervals$upper,
                                             method = "gdam05")

    gdamhscore_predict_intervals <- predict(fit_gdamhscore, newdata = simdat_test, se.fit = TRUE)
    gdamhscore_predict_intervals <- data.frame(trueeta = simdat_test$f,
                                           lower = gdamhscore_predict_intervals$lower,
                                           upper = gdamhscore_predict_intervals$upper,
                                           method = "gdamhscore")

    all_predict_intervals <- bind_rows(mgcv_predict_intervals,
                                       RobGam_predict_intervals,
                                       psidiv_predict_intervals,
                                       gdam01_predict_intervals,
                                       gdam05_predict_intervals,
                                       gdamhscore_predict_intervals
                                       )

    rm(mgcv_predict_intervals,
       RobGam_predict_intervals,
       psidiv_predict_intervals,
       gdam01_predict_intervals,
       gdam05_predict_intervals,
       gdamhscore_predict_intervals)

    interval_performance <- all_predict_intervals %>%
        mutate(inout = (lower < trueeta & upper > trueeta)) %>%
        mutate(method = fct_inorder(method)) %>%
        group_by(method) %>%
        summarise(coverage = mean(inout),
                  interval_widths = mean(upper - lower))


    set.seed(NULL)
    return(list(ME = ME,
                interval_performance = interval_performance,
                simdat_test = simdat_test,
                time_taken = all_time_taken))
    }


cw_results <- foreach(l0 = 1:num_datasets) %dopar% sim_fn(seed = l0, N = 400, prob_outlier = 0.1, true_sigma2 = 0.5)
save(cw_results, file = "setting2_N400_sigma205_prob01.RData")

## ---------------
sessioninfo::session_info()
## ---------------
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.3 (2024-02-29)
# os       Linux Mint 22.1
# system   x86_64, linux-gnu
# ui       RStudio
# language en_AU:en
# collate  en_AU.UTF-8
# ctype    en_AU.UTF-8
# tz       Australia/Sydney
# date     2025-06-25
# rstudio  2024.12.0+467 Kousa Dogwood (desktop)
# pandoc   3.2 @ /usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
#
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package      * version    date (UTC) lib source
# abind          1.4-8      2024-09-12 [1] CRAN (R 4.3.3)
# ADGofTest      0.3        2011-12-28 [1] CRAN (R 4.3.3)
# bitops         1.0-9      2024-10-03 [1] CRAN (R 4.3.3)
# cachem         1.0.8      2023-05-01 [3] CRAN (R 4.3.0)
# callr          3.7.5      2024-02-19 [3] CRAN (R 4.3.3)
# cli            3.6.4      2025-02-13 [1] CRAN (R 4.3.3)
# cluster        2.1.6      2023-12-01 [4] CRAN (R 4.3.3)
# codetools      0.2-19     2023-02-01 [4] CRAN (R 4.2.2)
# colorspace   * 2.1-1      2024-07-26 [1] CRAN (R 4.3.3)
# copula         1.1-6      2025-03-18 [1] CRAN (R 4.3.3)
# crayon         1.5.2      2022-09-29 [3] CRAN (R 4.2.2)
# curl           6.2.0      2025-01-23 [1] CRAN (R 4.3.3)
# DBI            1.2.2      2024-02-16 [3] CRAN (R 4.3.2)
# DEoptimR       1.1-3-1    2024-11-23 [1] CRAN (R 4.3.3)
# desc           1.4.3      2023-12-10 [3] CRAN (R 4.3.2)
# deSolve      * 1.40       2023-11-27 [3] CRAN (R 4.3.2)
# devtools       2.4.5      2022-10-11 [3] CRAN (R 4.2.2)
# digest         0.6.34     2024-01-11 [3] CRAN (R 4.3.2)
# distr          2.9.7      2025-01-12 [1] CRAN (R 4.3.3)
# distrEx        2.9.6      2025-01-13 [1] CRAN (R 4.3.3)
# doParallel   * 1.0.17     2022-02-07 [1] CRAN (R 4.3.3)
# DoubleRobGam * 0.1        2025-06-25 [1] Github (ilapros/DoubleRobGam@2808246)
# dplyr        * 1.1.4      2023-11-17 [1] CRAN (R 4.3.3)
# ellipsis       0.3.2      2021-04-29 [3] CRAN (R 4.1.1)
# evaluate       0.23       2023-11-01 [3] CRAN (R 4.3.2)
# evd            2.3-6.1    2022-07-04 [3] CRAN (R 4.2.1)
# farver         2.1.2      2024-05-13 [1] CRAN (R 4.3.3)
# fastmap        1.1.1      2023-02-24 [3] CRAN (R 4.2.2)
# fda          * 6.2.0      2024-09-17 [1] CRAN (R 4.3.3)
# fds          * 1.8        2018-10-31 [1] CRAN (R 4.3.3)
# forcats      * 1.0.0      2023-01-29 [3] CRAN (R 4.2.2)
# foreach      * 1.5.2      2022-02-02 [1] CRAN (R 4.3.3)
# fs             1.6.5      2024-10-30 [1] CRAN (R 4.3.3)
# gamlss.dist    6.1-1      2023-08-23 [1] CRAN (R 4.3.3)
# gdam         * 0.1        2025-02-28 [1] local
# generics       0.1.3      2022-07-05 [1] CRAN (R 4.3.3)
# ggokabeito     0.1.0      2021-10-18 [1] CRAN (R 4.3.3)
# ggplot2      * 3.5.1      2024-04-23 [1] CRAN (R 4.3.3)
# GJRM         * 0.2-6.8    2025-06-23 [1] CRAN (R 4.3.3)
# glue           1.8.0      2024-09-30 [1] CRAN (R 4.3.3)
# gmp            0.7-4      2024-01-15 [3] CRAN (R 4.3.2)
# gratia         0.10.0     2024-12-19 [1] CRAN (R 4.3.3)
# gsl            2.1-8      2023-01-24 [3] CRAN (R 4.2.2)
# gtable         0.3.6      2024-10-25 [1] CRAN (R 4.3.3)
# hdrcde         3.4        2021-01-18 [1] CRAN (R 4.3.3)
# here         * 1.0.1      2020-12-13 [3] CRAN (R 4.1.1)
# hms            1.1.3      2023-03-21 [3] CRAN (R 4.3.1)
# htmltools      0.5.8.1    2024-04-04 [1] CRAN (R 4.3.3)
# htmlwidgets    1.6.4      2023-12-06 [3] CRAN (R 4.3.2)
# httpuv         1.6.14     2024-01-26 [3] CRAN (R 4.3.3)
# ismev          1.42       2018-05-10 [1] CRAN (R 4.3.3)
# iterators    * 1.0.14     2022-02-05 [1] CRAN (R 4.3.3)
# kableExtra   * 1.4.0      2024-01-24 [3] CRAN (R 4.3.2)
# KernSmooth     2.23-22    2023-07-10 [4] CRAN (R 4.3.3)
# knitr          1.45       2023-10-30 [3] CRAN (R 4.3.2)
# ks             1.14.2     2024-01-15 [3] CRAN (R 4.3.2)
# later          1.3.2      2023-12-06 [3] CRAN (R 4.3.2)
# lattice        0.22-5     2023-10-24 [4] CRAN (R 4.3.3)
# lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.3.3)
# lubridate    * 1.9.3      2023-09-27 [3] CRAN (R 4.3.1)
# magic          1.6-1      2022-11-16 [3] CRAN (R 4.2.2)
# magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.3.3)
# MASS         * 7.3-60.0.1 2024-01-13 [4] CRAN (R 4.3.2)
# Matrix       * 1.6-5      2024-01-11 [4] CRAN (R 4.3.2)
# matrixStats    1.2.0      2023-12-11 [3] CRAN (R 4.3.2)
# maxLik       * 1.5-2      2021-07-26 [3] CRAN (R 4.1.1)
# mclust         6.0.1      2023-11-15 [3] CRAN (R 4.3.2)
# memoise        2.0.1      2021-11-26 [3] CRAN (R 4.1.2)
# mgcv         * 1.9-1      2023-12-21 [1] CRAN (R 4.3.3)
# mime           0.12       2021-09-28 [3] CRAN (R 4.1.1)
# miniUI         0.1.1.1    2018-05-18 [3] CRAN (R 4.0.1)
# miscTools    * 0.6-28     2023-05-03 [3] CRAN (R 4.3.1)
# mitools        2.4        2019-04-26 [3] CRAN (R 4.0.1)
# mnormt         2.1.1      2022-09-26 [3] CRAN (R 4.2.2)
# munsell        0.5.1      2024-04-01 [1] CRAN (R 4.3.3)
# mvnfast        0.2.8      2023-02-23 [1] CRAN (R 4.3.3)
# mvtnorm        1.3-3      2025-01-10 [1] CRAN (R 4.3.3)
# nlme         * 3.1-164    2023-11-27 [4] CRAN (R 4.3.3)
# numDeriv       2016.8-1.1 2019-06-06 [1] CRAN (R 4.3.3)
# patchwork      1.3.0      2024-09-16 [1] CRAN (R 4.3.3)
# pcaPP        * 2.0-4      2023-12-07 [3] CRAN (R 4.3.2)
# pillar         1.10.1     2025-01-07 [1] CRAN (R 4.3.3)
# pkgbuild       1.4.7      2025-03-24 [1] CRAN (R 4.3.3)
# pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.3.3)
# pkgload        1.3.4      2024-01-16 [3] CRAN (R 4.3.2)
# pracma         2.4.4      2023-11-10 [3] CRAN (R 4.3.2)
# processx       3.8.5      2025-01-08 [1] CRAN (R 4.3.3)
# profvis        0.3.8      2023-05-02 [3] CRAN (R 4.3.1)
# promises       1.2.1      2023-08-10 [3] CRAN (R 4.3.1)
# ps             1.8.1      2024-10-28 [1] CRAN (R 4.3.3)
# pspline        1.0-21     2024-12-11 [1] CRAN (R 4.3.3)
# psych          2.4.1      2024-01-18 [3] CRAN (R 4.3.2)
# purrr        * 1.0.4      2025-02-05 [1] CRAN (R 4.3.3)
# R6             2.6.1      2025-02-15 [1] CRAN (R 4.3.3)
# rainbow      * 3.8        2024-01-23 [1] CRAN (R 4.3.3)
# Rcpp           1.0.14     2025-01-12 [1] CRAN (R 4.3.3)
# RCurl        * 1.98-1.14  2024-01-09 [3] CRAN (R 4.3.3)
# readr        * 2.1.5      2024-01-10 [3] CRAN (R 4.3.2)
# remotes      * 2.4.2.1    2023-07-18 [3] CRAN (R 4.3.1)
# RhpcBLASctl    0.23-42    2023-02-11 [3] CRAN (R 4.3.2)
# rlang          1.1.5      2025-01-17 [1] CRAN (R 4.3.3)
# rmarkdown      2.25       2023-09-18 [3] CRAN (R 4.3.2)
# Rmpfr          0.9-5      2024-01-21 [3] CRAN (R 4.3.2)
# RobStatTM    * 1.0.11     2024-10-17 [1] CRAN (R 4.3.3)
# robustbase   * 0.99-4-1   2024-09-27 [1] CRAN (R 4.3.3)
# rprojroot      2.0.4      2023-11-05 [3] CRAN (R 4.3.2)
# rstudioapi     0.15.0     2023-07-07 [3] CRAN (R 4.3.1)
# sandwich       3.1-0      2023-12-11 [3] CRAN (R 4.3.2)
# scales         1.3.0      2023-11-28 [1] CRAN (R 4.3.3)
# scam           1.2-19     2025-05-26 [1] CRAN (R 4.3.3)
# sessioninfo    1.2.2      2021-12-06 [3] CRAN (R 4.1.2)
# sfsmisc        1.1-16     2023-08-10 [3] CRAN (R 4.3.1)
# shiny          1.8.0      2023-11-17 [3] CRAN (R 4.3.2)
# SparseM        1.84-2     2024-07-17 [1] CRAN (R 4.3.3)
# stabledist     0.7-1      2016-09-12 [3] CRAN (R 4.0.1)
# startupmsg     1.0.0      2025-01-12 [1] CRAN (R 4.3.3)
# stringi        1.8.4      2024-05-06 [1] CRAN (R 4.3.3)
# stringr      * 1.5.1      2023-11-14 [1] CRAN (R 4.3.3)
# survey         4.2-1      2023-05-03 [3] CRAN (R 4.3.1)
# survival       3.5-8      2024-02-14 [4] CRAN (R 4.3.2)
# svglite        2.1.3      2023-12-08 [3] CRAN (R 4.3.3)
# systemfonts    1.0.5      2023-10-09 [3] CRAN (R 4.3.1)
# tibble       * 3.2.1      2023-03-20 [1] CRAN (R 4.3.3)
# tidyr        * 1.3.1      2024-01-24 [1] CRAN (R 4.3.3)
# tidyselect     1.2.1      2024-03-11 [1] CRAN (R 4.3.3)
# tidyverse    * 2.0.0      2023-02-22 [3] CRAN (R 4.3.1)
# timechange     0.3.0      2024-01-18 [3] CRAN (R 4.3.2)
# trust          0.1-8      2020-01-10 [1] CRAN (R 4.3.3)
# tzdb           0.4.0      2023-05-12 [3] CRAN (R 4.3.1)
# urlchecker     1.0.1      2021-11-30 [3] CRAN (R 4.2.1)
# usethis        2.2.2      2023-07-06 [3] CRAN (R 4.3.1)
# utf8           1.2.4      2023-10-22 [1] CRAN (R 4.3.3)
# vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.3.3)
# VGAM           1.1-9      2023-09-19 [3] CRAN (R 4.3.1)
# VineCopula     2.6.1      2025-03-24 [1] CRAN (R 4.3.3)
# viridisLite    0.4.2      2023-05-02 [1] CRAN (R 4.3.3)
# withr          3.0.2      2024-10-28 [1] CRAN (R 4.3.3)
# xfun           0.52       2025-04-02 [1] CRAN (R 4.3.3)
# xml2           1.3.6      2023-12-04 [3] CRAN (R 4.3.2)
# xtable         1.8-4      2019-04-21 [3] CRAN (R 4.0.1)
# zoo            1.8-12     2023-04-13 [3] CRAN (R 4.3.1)
#
# [1] /home/fkch/R/x86_64-pc-linux-gnu-library/4.3
# [2] /usr/local/lib/R/site-library
# [3] /usr/lib/R/site-library
# [4] /usr/lib/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
