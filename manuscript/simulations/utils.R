## HIDDEN: Function to trick mgcv and subsequently gratia so that the right standard errors are obtained, along with everything else, when applying gratia::smooth_estimates
.calc_smooth_estimates <- function(gamObject, new_betas, new_Vp = NULL) {
    if(length(gamObject$smooth) == 0) 
        out <- NULL
    if(length(gamObject$smooth) > 0) {
        gamObject$coefficients <- new_betas
        if(!is.null(new_Vp))
            gamObject$Vp <- new_Vp
        out <- suppressMessages(gratia::smooth_estimates(object = gamObject, 
                                                         overall_uncertainty = TRUE))
        }
    
    return(out)
    }


## Extra function to construct predictions for DoubleRobGam model
predict_RobGam <- function(RobGam, newdata, choose_knots = 15) {
    newdata <- as.matrix(newdata)
    if(ncol(newdata) != (length(RobGam$dims)-1))
        stop("Number of columns in newdata must match the number of smooths in RobGam.")
    
    newMM <- do.call("cbind", lapply(1:ncol(newdata), function(k0) bsp(newdata[,k0], nknots = choose_knots)$B))
    newMM <- cbind(1, newMM)
    
    out <- as.vector(newMM %*% RobGam$coefficients)
    return(out)
}


## Extra function to construct predictions for list of univariate DPD GAMs. Not the intercept is not available/added here, and so the practitioner has to adjust this themselves downstream e.g., using the mean of the responses in the remaining data 
predict_dpd <- function(dpdGam, newdata) {
    if(ncol(newdata) != length(dpdGam))
        stop("Number of columns in newdata must match the length of dpdGam.")
    
    newMM <- do.call("cbind", lapply(1:ncol(newdata), function(k0) eval.basis(dpdGam[[k0]]$bsb, newdata[,k0])))
    
    out <- as.vector(newMM %*% as.vector(sapply(dpdGam, function(x) x$coefficients)))
    return(out)
    }


rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
    }


.siminterval_critvalue <- function(V, MM, se, level = 0.95) {
    BUdiff <- rmvn(n = 5000, mu = numeric(nrow(V)), sig = V)
    simDev <- tcrossprod(MM, BUdiff)
    absDev <- abs(sweep(simDev, 1L, se, FUN = "/"))
    masd <- apply(absDev, 2L, max)
    unname(quantile(masd, probs = level, type = 8))
    }

