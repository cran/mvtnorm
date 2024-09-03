
# mvnorm

### allow more than one distribution
mvnorm <- function(mean, chol, invchol) {

  # mvnorm chol invchol
  
  if (missing(chol) && missing(invchol))
      chol <- as.chol(ltMatrices(1, diag = TRUE))
  stopifnot(xor(missing(chol), missing(invchol)))

  if (!missing(chol)) {
      if (!is.ltMatrices(chol))
          chol <- as.ltMatrices(chol)
      scale <- as.chol(chol)
  }

  if (!missing(invchol)) {
      if (!is.ltMatrices(invchol))
          invchol <- as.ltMatrices(invchol)
      scale <- as.invchol(invchol)
  }
  ret <- list(scale = scale)
  
  # mvnorm mean
  
  if (!missing(mean)) {
      stopifnot(is.numeric(mean))
      stopifnot(NROW(mean) == dim(scale)[2L])
      if (!is.matrix(mean)) {
          mean <- matrix(mean, nrow = NROW(mean))
          rownames(mean) <- names(mean)
      }
      nm <- dimnames(scale)[[2L]]
      if (is.null(rownames(mean)))
          rownames(mean) <- nm
      if (!isTRUE(all.equal(rownames(mean), nm)))
          stop("rownames of mean do not match")        
      nm <- dimnames(scale)[[1L]]
      if (!is.null(nm) && dim(scale)[[2L]] == ncol(mean)) {
          if (is.null(colnames(mean)))
              colnames(mean) <- nm
          if (!isTRUE(all.equal(colnames(mean), nm)))
              stop("colnames of mean do not match")        
      }
      ret$mean <- mean
  }
  
  class(ret) <- "mvnorm"
  return(ret)
}

# mvnorm methods

names.mvnorm <- function(x)
    dimnames(x$scale)[[2L]]

aperm.mvnorm <- function(a, perm, ...) {

    ret <- list(scale = aperm(a$scale, perm = perm, ...))
    if (!is.null(a$mean))
        ret$mean <- a$mean[perm,,drop = FALSE]
    class(ret) <- "mvnorm"
    ret
}

# mvnorm simulate

simulate.mvnorm <- function(object, nsim = dim(object$scale)[1L], seed = NULL, 
                            standardize = FALSE, as.data.frame = FALSE, ...) {

    J <- dim(object$scale)[2L]
    N <- dim(object$scale)[1L]
    if (N > 1)
        stopifnot(nsim == N)
    if (standardize) {
        if (is.chol(object$scale)) {
            object$scale <- standardize(chol = object$scale)
        } else {
            object$scale <- standardize(invchol = object$scale)
        }
    }
    Z <- matrix(rnorm(nsim * J), nrow = J)
    if (is.chol(object$scale)) {
        Y <- Mult(object$scale, Z)
    } else {
        Y <- solve(object$scale, Z)
    }
    ret <- Y
    if (!is.null(object$mean))
        ret <- ret + c(object$mean)
    rownames(ret) <- dimnames(object$scale)[[2L]]
    if (!as.data.frame)
        return(ret)
    return(as.data.frame(t(ret)))
}

# mvnorm margDist

margDist <- function(object, which, ...)
    UseMethod("margDist")

margDist.mvnorm <- function(object, which, ...) {

    if (is.chol(object$scale)) {
        ret <- list(scale = as.chol(marg_mvnorm(chol = object$scale, 
                                                which = which)$chol))
    } else {
        ret <- list(scale = as.invchol(marg_mvnorm(invchol = object$scale, 
                                                   which = which)$invchol))
    }
    if (!is.null(object$mean))
        ret$mean <- object$mean[which,,drop = FALSE]
    class(ret) <- "mvnorm"
    return(ret)
}

# mvnorm condDist

condDist <- function(object, which_given, given, ...)
    UseMethod("condDist")

condDist.mvnorm <- function(object, which_given = 1L, given, ...) {

    if (is.chol(object$scale)) {
        ret <- cond_mvnorm(chol = object$scale, which_given = which_given, 
                           given = given, ...)
        ret$scale <- as.chol(ret$chol)
        ret$chol <- NULL
    } else {
        ret <- cond_mvnorm(invchol = object$scale, which_given = which_given, 
                           given = given, ...)
        ret$invchol <- as.chol(ret$invchol)
        ret$invchol <- NULL
    }
    if (!is.null(object$mean)) {
        if (is.character(which_given)) 
            which_given <- match(which_given, dimnames(object$scale)[[2L]])
        if (ncol(object$mean) > 1L && ncol(ret$mean) > 1)
            stop("dimensions do not match")
        if (ncol(object$mean) == 1L && ncol(ret$mean) > 1L) {
            ret$mean <- object$mean[-which_given,,drop = TRUE] + ret$mean
        } else {
            ret$mean <- object$mean[-which_given,,drop = FALSE] + c(ret$mean)
        }
        
    }
    class(ret) <- "mvnorm"
    return(ret)
}

# mvnorm logLik

logLik.mvnorm <- function(object, obs, lower, upper, standardize = FALSE, 
                          ...) {
    # argchecks
    
    args <- c(object, list(...))
    nargs <- missing(obs) + missing(lower) + missing(upper)
    stopifnot(nargs < 3L)

    nmobs <- NULL
    if (!missing(obs)) {
        if (!is.null(obs)) {
            stopifnot(is.matrix(obs))
            nmobs <- rownames(obs)
        }
    }
    nmlower <- nmupper <- nmlu <- NULL
    if (!missing(lower)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(lower))
            nmlu <- nmlower <- rownames(lower)
        }
    }
    if (!missing(upper)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(upper))
            nmupper <- rownames(upper)
            if (!missing(lower)) {
                stopifnot(isTRUE(all.equal(nmlower, nmupper)))
            } else {
                nmlu <- nmupper
            }
        }
    }

    nm <- c(nmobs, nmlu)
    no <- names(object)
    stopifnot(nm %in% no)
    perm <- NULL
    if (!isTRUE(all.equal(nm, no)))
        perm <- c(nm, no[!no %in% nm])

    if (!missing(obs)) args$obs <- obs
    if (!missing(lower)) args$lower <- lower
    if (!missing(upper)) args$upper <- upper
    
    if (is.chol(object$scale)) {
        # logLik chol
        
        names(args)[names(args) == "scale"] <- "chol"
        if (standardize)
            args$chol <- standardize(chol = args$chol)
        if (!is.null(perm)) {
            args$chol <- aperm(as.chol(args$chol), perm = perm)
            if (length(nm) < length(no))
                args$chol <- marg_mvnorm(chol = args$chol, which = nm)$chol
            args$mean <- args$mean[nm,,drop = FALSE]
        }
        return(do.call("ldpmvnorm", args))
        
    }
    # logLik invchol
    
    names(args)[names(args) == "scale"] <- "invchol"
    if (standardize)
        args$invchol <- standardize(invchol = args$invchol)
    if (!is.null(perm)) {
        args$invchol <- aperm(as.invchol(args$invchol), perm = perm)
        if (length(nm) < length(no))
            args$invchol <- marg_mvnorm(invchol = args$invchol, 
                                        which = nm)$invchol
        args$mean <- args$mean[nm,,drop = FALSE]
    }
    return(do.call("ldpmvnorm", args))
    
}

# mvnorm lLgrad

lLgrad <- function(object, ...)
    UseMethod("lLgrad")

lLgrad.mvnorm <- function(object, obs, lower, upper, standardize = FALSE, 
                          ...) {
    # argchecks
    
    args <- c(object, list(...))
    nargs <- missing(obs) + missing(lower) + missing(upper)
    stopifnot(nargs < 3L)

    nmobs <- NULL
    if (!missing(obs)) {
        if (!is.null(obs)) {
            stopifnot(is.matrix(obs))
            nmobs <- rownames(obs)
        }
    }
    nmlower <- nmupper <- nmlu <- NULL
    if (!missing(lower)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(lower))
            nmlu <- nmlower <- rownames(lower)
        }
    }
    if (!missing(upper)) {
        if (!is.null(lower)) {
            stopifnot(is.matrix(upper))
            nmupper <- rownames(upper)
            if (!missing(lower)) {
                stopifnot(isTRUE(all.equal(nmlower, nmupper)))
            } else {
                nmlu <- nmupper
            }
        }
    }

    nm <- c(nmobs, nmlu)
    no <- names(object)
    stopifnot(nm %in% no)
    perm <- NULL
    if (!isTRUE(all.equal(nm, no)))
        perm <- c(nm, no[!no %in% nm])

    if (!missing(obs)) args$obs <- obs
    if (!missing(lower)) args$lower <- lower
    if (!missing(upper)) args$upper <- upper
    
    if (is.chol(object$scale)) {
        # lLgrad chol
        
        names(args)[names(args) == "scale"] <- "chol"
        sc <- args$chol
        if (standardize)
            args$chol <- sc <- standardize(chol = args$chol)
        if (!is.null(perm)) {
            if (!attr(args$chol, "diag")) {
                diagonals(args$chol) <- 1
                sc <- args$chol
            }
            args$chol <- pc <- aperm(as.chol(args$chol), perm = perm)
            if (length(nm) < length(no))
                args$chol <- marg_mvnorm(chol = args$chol, which = nm)$chol
            args$mean <- args$mean[nm,,drop = FALSE]
        }
        ret <- do.call("sldpmvnorm", args)
        # lLgrad mean

        ### sldmvnorm returns mean score as -obs
        if (is.null(ret$mean)) ret$mean <- - ret$obs
        
        # lLgrad marginalisation

        om <- length(no) - length(nm)
        if (om > 0) {
            am <- matrix(0, nrow = om, ncol = ncol(ret$mean))
            rownames(am) <- no[!no %in% nm]
            ret$mean <- rbind(ret$mean, am)
            Jo <- dim(object$scale)[[2L]]
            pJ <- dim(args$invchol)[[2L]]
            am <- matrix(0, nrow = Jo * (Jo + 1) / 2 - pJ * (pJ + 1) / 2, 
                         ncol = dim(ret$invchol)[1L])
            byrow_orig <- attr(ret$chol, "byrow")
            ret$chol <- ltMatrices(ret$chol, byrow = TRUE)
            ### rbind only works for byrow = TRUE
            ret$chol <- ltMatrices(rbind(unclass(ret$chol), am), 
                                   byrow = TRUE, 
                                   diag = TRUE,
                                   names = perm)
            ret$chol <- ltMatrices(ret$chol, byrow = byrow_orig)
        }
        
        # lLgrad deperma

        if (!is.null(perm))
            ret$chol <- deperma(chol = sc, permuted_chol = pc, 
                                perm = match(perm, no), 
                                score_schol = ret$chol)
        
        # lLgrad destandarized

        if (standardize)
            ret$chol <- destandardize(chol = object$scale, 
                                      score_schol = ret$chol)
        
        # lLgrad diagonals

        if (!attr(sc, "diag"))
            ret$chol <- ltMatrices(Lower_tri(ret$chol, diag = FALSE),
                                   diag = FALSE, 
                                   byrow = attr(ret$chol, "byrow"), 
                                   names = dimnames(ret$chol)[[2L]])
        
        # lLgrad return

        ret$scale <- ret$chol
        ret$chol <- NULL
        ret$mean <- ret$mean[no,,drop = FALSE]
        return(ret)
        
        
    }
    # lLgrad invchol
    
    names(args)[names(args) == "scale"] <- "invchol"
    si <- args$invchol
    if (standardize)
        args$invchol <- si <- standardize(invchol = args$invchol)
    if (!is.null(perm)) {
        if (!attr(args$invchol, "diag")) {
            diagonals(args$invchol) <- 1
            si <- args$invchol
        }
        args$invchol <- pi <- aperm(as.invchol(args$invchol), perm = perm)
        if (length(nm) < length(no))
            args$invchol <- marg_mvnorm(invchol = args$invchol,
                                        which = nm)$invchol
        args$mean <- args$mean[nm,,drop = FALSE]
    }
    ret <- do.call("sldpmvnorm", args)
    ### sldmvnorm returns mean score as -obs
    if (is.null(ret$mean)) ret$mean <- - ret$obs
    om <- length(no) - length(nm)
    if (om > 0) {
        am <- matrix(0, nrow = om, ncol = ncol(ret$mean))
        rownames(am) <- no[!no %in% nm]
        ret$mean <- rbind(ret$mean, am)
        Jo <- dim(object$scale)[[2L]]
        pJ <- dim(args$invchol)[[2L]]
        am <- matrix(0, nrow = Jo * (Jo + 1) / 2 - pJ * (pJ + 1) / 2, 
                     ncol = dim(ret$invchol)[1L])
        byrow_orig <- attr(ret$invchol, "byrow")
        ret$invchol <- ltMatrices(ret$invchol, byrow = TRUE)
        ### rbind only works for byrow = TRUE
        ret$invchol <- ltMatrices(rbind(unclass(ret$invchol), am), 
                                  byrow = TRUE,
                                  diag = TRUE,
                                  names = perm)
        ret$invchol <- ltMatrices(ret$invchol, byrow = byrow_orig)
    }
    if (!is.null(perm))
        ret$invchol <- deperma(invchol = si, permuted_invchol = pi, 
                               perm = match(perm, no), 
                               score_schol = -vectrick(pi, ret$invchol))
    if (standardize)
        ret$invchol <- destandardize(invchol = object$scale, 
                                     score_schol = -vectrick(si, ret$invchol))
    if (!attr(si, "diag"))
        ret$invchol <- ltMatrices(Lower_tri(ret$invchol, diag = FALSE),
                                  diag = FALSE, 
                                  byrow = attr(ret$invchol, "byrow"), 
                                  names = dimnames(ret$invchol)[[2L]])
    ret$scale <- ret$invchol
    ret$invchol <- NULL
    ret$mean <- ret$mean[no,,drop = FALSE]
    return(ret)
    
}

