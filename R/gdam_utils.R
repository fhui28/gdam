#' @title Get the full S matrix from GAMs. Relies on the fact that gam always move the parametric terms first
#' @noRd
#' @noMd

.get_bigS <- function(gamObject, smoothing_parameters) {
    num_MM <- ncol(mgcv::model.matrix.gam(gamObject))

    bigS <- Matrix::Matrix(0, num_MM, num_MM, sparse = TRUE)
    num_smooth_terms <- length(gamObject$smooth)
    if(num_smooth_terms == 0)
        return(bigS)

    num_Smatrices_per_smooth <- lapply(gamObject$smooth, function(x) length(x$S)) # The sum of this should equal length(gamObject$sp)
    sp_index <- split(1:length(gamObject$sp), rep(1:num_smooth_terms, num_Smatrices_per_smooth)) # Because fs, te, and ti smooths have multiple S and smoothing parameters, then this tells you how many and indexes the S/sp's within each smooth term. This is very similar to extracting first.sp and last.sp from each smooth
    num_smooth_cols <- sum(sapply(gamObject$smooth, function(x) x$last.para - x$first.para + 1))

    subS <- lapply(1:num_smooth_terms, function(j) {
        out <- smoothing_parameters[sp_index[[j]][1]] * gamObject$smooth[[j]]$S[[1]]
        if(length(sp_index[[j]]) > 1) { # To deal with smooths that have multiple S matrices and smoothing parameters
            for(l0 in 2:length(sp_index[[j]]))
                out <- out + smoothing_parameters[sp_index[[j]][l0]] * gamObject$smooth[[j]]$S[[l0]]
        }
        return(out)
    })

    subS <- Matrix::bdiag(subS)
    bigS[-(1:gamObject$nsdf), -(1:gamObject$nsdf)] <- subS

    return(bigS)
    }


#' @title Get the full list of S matrices from GAMs, converting them to be of the same dimension.
#' @description Assumption is made that there are always smooths in the model!
#' @noRd
#' @noMd

.get_Slists <- function(gamObject) {
    num_MM <- ncol(mgcv::model.matrix.gam(gamObject))
    num_smooth_terms <- length(gamObject$smooth)
    num_Smatrices_per_smooth <- lapply(gamObject$smooth, function(x) length(x$S)) # The sum of this should equal length(gamObject$sp)
    sp_index <- split(1:length(gamObject$sp), rep(1:num_smooth_terms, num_Smatrices_per_smooth)) # Because fs, te, and ti smooths have multiple S and smoothing parameters, then this tells you how many and indexes the S/sp's within each smooth term. This is very similar to extracting first.sp and last.sp from each smooth


    subS <- lapply(1:num_smooth_terms, function(j) {
        select_indices <- gamObject$smooth[[j]]$first.para:gamObject$smooth[[j]]$last.para
        if(length(sp_index[[j]]) == 1) {
            out <- Matrix::Matrix(0, num_MM, num_MM, sparse = TRUE)
            out[select_indices, select_indices] <- gamObject$smooth[[j]]$S[[1]]
            }
        if(length(sp_index[[j]]) > 1) { # To deal with smooths that have multiple S matrices
            out <- vector("list", length(sp_index[[j]]))
            for(l0 in 1:length(sp_index[[j]])) {
                out[[l0]] <- Matrix::Matrix(0, num_MM, num_MM, sparse = TRUE)
                out[[l0]][select_indices, select_indices] <- gamObject$smooth[[j]]$S[[l0]]
                }
            }
        return(out)
        })

    subS <- do.call(c, subS)
    return(subS)
    }


#' @noRd
#' @noMd
.notExp_grad <- function(x) {
    f <- x
    ind <- which(x > 1)
    f[ind] <- exp(1) * x[ind]
    ind <- which((x <= 1) & (x > -1))
    f[ind] <- exp(x[ind])
    ind <- which(x <= -1)
    f[ind] <- -x[ind] / exp(1) * ((x[ind]^2 + 1)/2)^-2

    return(f)
    }


#' @title A local pseudo-inverse function -- straight from [mgcv::summary.gam()] package. Full credit goes to Simon Wood for this!
#' @noRd
#' @noMd

.pinv <- function(V, M, rank.tol = 1e-6) {
    if(missing(M))
        M <- ncol(V)
    D <- eigen(V, symmetric = TRUE)
    M1 <- length(D$values[D$values > rank.tol * D$values[1]])
    if(M > M1)
        M<-M1 # avoid problems with zero eigen-values
    if(M+1 <= length(D$values))
        D$values[(M+1):length(D$values)]<-1
    D$values<- 1/D$values
    if(M+1 <= length(D$values))
        D$values[(M+1):length(D$values)]<-0
    res <- D$vectors %*% diag(x = D$values, nrow = length(D$values)) %*% t(D$vectors)

    return(res)
    }


