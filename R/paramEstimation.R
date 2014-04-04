# file : paramEstimation.R author: Gavin Ha
# <gha@bccrc.ca> Dept of Molecular Oncolgy British
# Columbia Cancer Agency University of British
# Columbia date : March 3, 2014

## wrapper that makes subroutine calls to
## updateParameters based on the normal estimate
## method chosen
estimateClonalCNParamsMap <- function(x, N, l, rho, 
    n_0, s_0, var_0, phi_0, gParams, nParams, sParams, 
    pParams, maxiter = 500, normalEstimateMethod = "fixed", 
    estimateS = TRUE, estimatePloidy = TRUE, verbose = TRUE) {
    K <- dim(rho)[1]
    Z <- dim(rho)[2]
    T <- dim(rho)[3]
    # assign some variables to simpify the update
    # equations
    a <- matrix(0, K, Z)
    b <- matrix(0, K, Z)
    c <- matrix(0, K, Z)
    d <- matrix(0, K, Z)
    e <- matrix(0, K, Z)
    g <- matrix(0, K, Z)
    for (z in 1:Z) {
        for (k in 1:K) {
            a[k, z] <- x %*% rho[k, z, ]
            b[k, z] <- (N - x) %*% rho[k, z, ]
            c[k, z] <- l %*% rho[k, z, ]
            e[k, z] <- sum(rho[k, z, ])
            g[k, z] <- (l^2) %*% rho[k, z, ]
            d[k, z] <- lchoose(N, x) %*% rho[k, z, 
                ]
        }
    }
    
    estimateOut <- updateParameters(n_0, s_0, var_0, 
        phi_0, a, b, c, d, e, g, gParams, nParams, 
        sParams, pParams, normalEstimateMethod = normalEstimateMethod, 
        estimateS = estimateS, estimatePloidy = estimatePloidy, 
        maxiter = maxiter, verbose = verbose)
    
    if (verbose == TRUE) {
        message("Using Coordinate Descent iteration ", 
            estimateOut$maxFind, " with Fval=", format(estimateOut$maxF, 
                digits = 4), " and n=", format(estimateOut$n, 
                digits = 4), " (", normalEstimateMethod, 
            ")", ", s=[", paste(format(estimateOut$s, 
                digits = 2), collapse = ","), "], phi=", 
            format(estimateOut$phi, digits = 4))
    }
    return(estimateOut)
}


## handles normal contamination estimation using 1)
## map and 2) fixed approaches
updateParameters <- function(n_0, s_0, var_0, phi_0, 
    a, b, c, d, e, g, gParams, nParams, sParams, pParams, 
    normalEstimateMethod, estimateS, estimatePloidy, 
    maxiter = 500, miniter = 10, verbose = TRUE) {
    
    K <- length(var_0)
    Z <- length(s_0)
    s <- rep(0, Z)
    var = rep(0, K)
    phi = 0
    n = 0
    
    intervalS <- c(.Machine$double.eps, 1 - .Machine$double.eps)
    intervalVar <- c(.Machine$double.eps, 5)
    intervalPhi <- c(.Machine$double.eps, 10)
    intervalN <- c(.Machine$double.eps, 1 - .Machine$double.eps)
    converged <- 0
    i <- 1
    maxFind <- 1
    objfun <- rep(0, maxiter)
    n_prev <- n_0
    s_prev <- s_0
    var_prev <- var_0
    phi_prev <- phi_0
    F <- matrix(0, K + Z + 2, maxiter)
    F[, i] <- c(n_prev, s_prev, var_prev, phi_prev)
    F0 <- rep(0, maxiter)  ## keep track of likelihood function (objective function)
    F0[i] <- likelihoodFunc(F[, i], a, b, c, d, e, 
        g, gParams$rt, gParams$rn, gParams$ct, sParams$alphaSHyper, 
        sParams$betaSHyper, nParams$alphaNHyper, nParams$betaNHyper, 
        gParams$alphaKHyper, gParams$betaKHyper, pParams$alphaPHyper, 
        pParams$betaPHyper, normalEstimateMethod, estimateS, 
        estimatePloidy)
    objfun[i] <- F0[i]
    if (verbose == TRUE) {
        message("Coord Descent[", i, "]=", format(objfun[i], 
            digits = 4), " n=", format(n_prev, digits = 4), 
            " s=[", paste(format(s_prev, digits = 2), 
                collapse = ","), "]", " phi=", format(phi_prev, 
                digits = 4))
    }
    
    ### begin coordinate descent ###
    while ((!converged && (i < maxiter)) || (i <= miniter)) {
        i <- i + 1
        if (normalEstimateMethod == "map") {
            funcN <- function(n) {
                clonalCNDerivativeNUpdateEqn(n, s_prev, 
                  var_prev, phi_prev, a, b, c, e, gParams$rt, 
                  gParams$rn, gParams$ct, nParams$alphaNHyper, 
                  nParams$betaNHyper)
            }
            n <- uniroot(funcN, intervalN, tol = 1e-15)$root
        } else if (normalEstimateMethod == "fixed") {
            n <- n_prev
        }
        
        if (normalEstimateMethod != "optim" && estimateS) {
            # estimate S independently since not using optim
            for (z in 1:Z) {
                funcS <- function(s) {
                  clonalCNDerivativeSUpdateEqn(s, n, 
                    var_prev, phi_prev, a[, z], b[, 
                      z], c[, z], e[, z], gParams$rt, 
                    gParams$rn, gParams$ct, sParams$alphaSHyper[z], 
                    sParams$betaSHyper[z])
                }
                s[z] <- uniroot(funcS, intervalS, tol = 1e-15)$root
            }
        } else {
            s <- s_prev
            
        }
        
        ## Estimate Gaussian variance for each copy number
        ## level, rather than each genotype state ##
        cnLevel <- unique(gParams$ct)
        for (ck in 1:length(cnLevel)) {
            cnInd <- which(gParams$ct == cnLevel[ck])
            varTmp <- clonalCNDerivativeVarExactUpdateEqn(n, 
                s, phi_prev, c, e, g, gParams$rt, gParams$rn, 
                gParams$ct, cnLevel[ck], gParams$alphaKHyper[cnInd[1]], 
                gParams$betaKHyper[cnInd[1]])
            var[cnInd] <- varTmp
        }
        
        
        if (estimatePloidy) {
            funcP <- function(phi) {
                clonalCNDerivativePloidyUpdateEqn(phi, 
                  n_prev, s, var, c, e, gParams$rt, 
                  gParams$rn, gParams$ct, pParams$alphaPHyper, 
                  pParams$betaPHyper)
            }
            phi <- uniroot(funcP, intervalPhi, tol = 1e-15)$root
        } else {
            # estimatePloidy==FALSE
            phi <- phi_prev
        }
        
        
        
        F[, i] <- c(n, s, var, phi)
        F0[i] <- likelihoodFunc(F[, i], a, b, c, d, 
            e, g, gParams$rt, gParams$rn, gParams$ct, 
            sParams$alphaSHyper, sParams$betaSHyper, 
            nParams$alphaNHyper, nParams$betaNHyper, 
            gParams$alphaKHyper, gParams$betaKHyper, 
            pParams$alphaPHyper, pParams$betaPHyper, 
            normalEstimateMethod, estimateS, estimatePloidy)
        objfun[i] <- F0[i]
        if (verbose == TRUE) {
            message("Coord Descent[", i, "]=", format(objfun[i], 
                digits = 4), " n=", format(n, digits = 4), 
                " s=[", paste(format(s, digits = 2), 
                  collapse = ","), "]", " phi=", format(phi, 
                  digits = 4))
        }
        n_prev <- n
        s_prev <- s
        var_prev <- var
        phi_prev <- phi
        
        if (objfun[maxFind] >= objfun[i]) {
            message("Coordinate descent decreases likelihood.")
        } else {
            maxFind <- i
        }
        if ((abs(objfun[i] - objfun[i - 1])/abs(objfun[i])) <= 
            0.001) {
            converged <- 1
        }
    }
    # return estimated parameters
    output <- vector("list", 0)
    output$maxFind <- maxFind
    output$maxF <- objfun[maxFind]
    output$n <- F[1, maxFind]
    output$s <- F[2:(Z + 1), maxFind]
    output$var <- F[(Z + 2):(dim(F)[1] - 1), maxFind]
    output$phi <- F[dim(F)[1], maxFind]
    return(output)
}

clonalNSUpdateEqn <- function(unk, var, phi, a, b, 
    c, e, gParams, nParams, sParams) {
    K <- length(gParams$alphaKHyper)
    Z <- length(sParams$alphaSHyper)
    F <- rep(0, Z + 1)
    
    n <- unk[1]
    s <- unk[2:(Z + 1)]
    
    F[1] <- clonalCNDerivativeNUpdateEqn(n, s, var, 
        phi, a, b, c, e, gParams$rt, gParams$rn, gParams$ct, 
        nParams$alphaNHyper, nParams$betaNHyper)
    
    for (z in 1:Z) {
        F[z + 1] <- clonalCNDerivativeSUpdateEqn(s[z], 
            n, var, phi, a[, z], b[, z], c[, z], e[, 
                z], gParams$rt, gParams$rn, gParams$ct, 
            sParams$alphaSHyper[z], sParams$betaSHyper[z])
    }
    # return(F)
    return(sum(abs(F)))
}

clonalCNSVarUpdateEqn <- function(unk, a, b, c, d, 
    e, g, gParams, nParams, sParams, pParams, normalEstimateMethod, 
    estimateS, estimatePloidy) {
    K <- length(gParams$alphaKHyper)
    Z <- length(sParams$alphaSHyper)
    F <- rep(0, K + Z + 2)
    
    n <- unk[1]
    s <- unk[2:(Z + 1)]
    var <- unk[(Z + 2):(length(unk) - 1)]
    phi <- unk[length(unk)]
    
    if (normalEstimateMethod == "fixed") {
        F[1] <- 0
    } else if (normalEstimateMethod == "map") {
        F[1] <- clonalCNDerivativeNUpdateEqn(n, s, 
            var, phi, a, b, c, e, gParams$rt, gParams$rn, 
            gParams$ct, nParams$alphaNHyper, nParams$betaNHyper)
    }
    
    for (z in 1:Z) {
        if (estimateS) {
            F[z + 1] <- clonalCNDerivativeSUpdateEqn(s[z], 
                n, var, phi, a[, z], b[, z], c[, z], 
                e[, z], gParams$rt, gParams$rn, gParams$ct, 
                sParams$alphaSHyper[z], sParams$betaSHyper[z])
        } else {
            # estimateS ==FALSE
            F[z + 1] <- 0
        }
    }
    
    cnLevel <- unique(gParams$ct)
    for (ck in 1:length(cnLevel)) {
        cnInd <- which(gParams$ct == cnLevel[ck])
        varTmp <- 0
        for (k in cnInd) {
            varTmp <- varTmp + clonalCNDerivativeVarUpdateEqn(var[k], 
                n, s, phi, c[k, ], e[k, ], g[k, ], 
                gParams$rt[k], gParams$rn, cnLevel[ck], 
                gParams$alphaKHyper[k], gParams$betaKHyper[k])
        }
        F[(Z + 1) + cnInd] <- varTmp
    }
    
    if (estimatePloidy) {
        F[K + Z + 2] <- clonalCNDerivativePloidyUpdateEqn(phi, 
            n, s, var, c, e, gParams$rt, gParams$rn, 
            gParams$ct, pParams$alphaPHyper, pParams$betaPHyper)
    } else {
        # estimatePloidy==FALSE
        F[K + Z + 2] <- 0
    }
    return(F)
}


## computed once for all parameters ##
likelihoodFunc <- function(unk, a, b, c, d, e, g, rt, 
    rn, ct, alphaS, betaS, alphaN, betaN, alphaK, betaK, 
    alphaP, betaP, normalEstimateMethod, estimateS, 
    estimatePloidy) {
    K <- length(rt)
    Z <- length(alphaS)
    
    n <- unk[1]
    s <- unk[2:(Z + 1)]
    var <- unk[(Z + 2):(length(unk) - 1)]
    phi <- unk[length(unk)]
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    lik = 0
    ## find copy number level indices - short list of
    ## the full genotype indices
    cnLevel <- unique(ct)
    cnInd <- sapply(cnLevel, function(x) {
        which(ct == x)[1]
    })
    
    for (z in 1:Z) {
        for (k in 1:K) {
            lik <- lik + d[k, z] + a[k, z] * log(mus$R[k, 
                z]) + b[k, z] * log(1 - mus$R[k, z]) + 
                e[k, z] * log(1/sqrt(2 * pi * var[k])) - 
                (g[k, z] - 2 * c[k, z] * mus$C[k, z] + 
                  e[k, z] * mus$C[k, z]^2)/(2 * var[k])
        }
    }
    
    ## prior for S ##
    prior_beta_s <- 0
    if (estimateS) {
        for (z in 1:Z) {
            prior_beta_s <- prior_beta_s + log(1/beta(alphaS[z], 
                betaS[z])) + (alphaS[z] - 1) * log(s[z]) + 
                (betaS[z] - 1) * log(1 - s[z])
            # prior_beta_s <- prior_beta_s +
            # dbeta(s[z],alphaS[z],betaS[z],log=TRUE)
        }
    } else {
        prior_beta_s <- 0
    }
    ## prior for n ##
    if (normalEstimateMethod == "map") {
        prior_beta_n <- log(1/beta(alphaN, betaN)) + 
            (alphaN - 1) * log(n) + (betaN - 1) * log(1 - 
            n)
        # prior_beta_n <- dbeta(n,alphaN,betaN,log=TRUE)
    } else {
        prior_beta_n <- 0
    }
    ## prior for variance var ##
    prior_gamma_var <- 0
    for (k in cnInd) {
        prior_gamma_var <- prior_gamma_var + alphaK[k] * 
            log(betaK[k]) - lgamma(alphaK[k]) + (-alphaK[k] - 
            1) * log(var[k]) - betaK[k]/var[k]
    }
    ## prior for phi ##
    if (estimatePloidy) {
        prior_gamma_phi <- alphaP * log(betaP) - lgamma(alphaP) + 
            (-alphaP - 1) * log(phi) - betaP/phi
    } else {
        prior_gamma_phi <- 0
    }
    
    return(lik + prior_beta_s + prior_beta_n + prior_gamma_var + 
        prior_gamma_phi)
}


clonalCNDerivativeNUpdateEqn <- function(n, s, var, 
    phi, a, b, c, e, rt, rn, ct, alphaN, betaN) {
    K <- length(ct)
    Z <- length(s)
    cn <- 2
    # data likelihood derivative wrt to n
    dlik_dn <- 0
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    
    for (z in 1:Z) {
        dmuC_dn <- (cn - s[z] * cn - (1 - s[z]) * ct)/(n * 
            cn + (1 - n) * s[z] * cn + (1 - n) * (1 - 
            s[z]) * ct) - (cn - phi)/(n * cn + (1 - 
            n) * phi)
        dmuR_dn <- ((1 - s[z]) * cn * ct * (rn - rt))/((n * 
            cn + (1 - n) * s[z] * cn + (1 - n) * (1 - 
            s[z]) * ct)^2)
        for (k in 1:K) {
            term1 <- a[k, z] * dmuR_dn[k]/mus$R[k, 
                z]
            term2 <- b[k, z] * dmuR_dn[k]/(1 - mus$R[k, 
                z])
            term3 <- (c[k, z] - e[k, z] * mus$C[k, 
                z]) * dmuC_dn[k]/var[k]
            dlik_dn <- dlik_dn + (term1 - term2 + term3)
        }
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dbeta_dn <- (alphaN - 1)/n - (betaN - 1)/(1 - n)
    f <- dlik_dn + dbeta_dn
    return(f)
}

clonalCNDerivativeSUpdateEqn <- function(s, n, var, 
    phi, a, b, c, e, rt, rn, ct, alphaS, betaS) {
    K <- length(ct)
    cn <- 2
    # data likelihood derivative wrt to s
    dlik_ds <- 0
    
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    dmuC_ds <- ((1 - n) * cn - (1 - n) * ct)/(n * cn + 
        (1 - n) * s * cn + (1 - n) * (1 - s) * ct)
    dmuR_ds <- ((1 - n) * cn * ct * (rn - rt))/((n * 
        cn + (1 - n) * s * cn + (1 - n) * (1 - s) * 
        ct)^2)
    
    for (k in 1:K) {
        term1 <- a[k] * dmuR_ds[k]/mus$R[k]
        term2 <- b[k] * dmuR_ds[k]/(1 - mus$R[k])
        term3 <- (c[k] - e[k] * mus$C[k]) * dmuC_ds[k]/var[k]
        dlik_ds <- dlik_ds + (term1 - term2 + term3)
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dbeta_ds <- (alphaS - 1)/s - (betaS - 1)/(1 - s)
    f <- dlik_ds + dbeta_ds
    return(f)
}

# ck is the copy number level that we want to
# marginalize over ct is vector of copy number;
# each corresponding to genotype state
clonalCNDerivativeVarExactUpdateEqn <- function(n, 
    s, phi, c, e, g, rt, rn, ct, ck, alphaK, betaK) {
    Z <- length(s)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    cnInd <- which(ct == ck)
    var <- 0
    term1 <- 0
    term2 <- 0
    for (k in cnInd) {
        for (z in 1:Z) {
            term1 = term1 + -sum(g[k, z] - 2 * c[k, 
                z] * mus$C[k, z] + e[k, z] * (mus$C[k, 
                z])^2)
            term2 = term2 + -sum(e[k, z])
        }
    }
    term1 = term1 - 2 * betaK
    term2 = term2 - 2 * (alphaK + 1)
    var = term1/term2
    return(var)
}

clonalCNDerivativeVarUpdateEqn <- function(var, n, 
    s, phi, c, e, g, rt, rn, ct, alphaK, betaK) {
    Z <- length(s)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    
    # data likelihood derivative wrt to sigma
    dlik_dvar <- 0
    
    for (z in 1:Z) {
        term1 = -e[z]/(2 * var)
        term2 = (g[z] + e[z] * (mus$C[z]^2))/(2 * var^2) - 
            c[z] * mus$C[z]/var^2
        dlik_dvar = dlik_dvar + term1 + term2
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dinvgamma_dvar = (-alphaK - 1)/var + betaK/var^2
    f = dlik_dvar + dinvgamma_dvar
    return(f)
}

clonalCNDerivativePloidyUpdateEqn <- function(phi, 
    n, s, var, c, e, rt, rn, ct, alphaP, betaP) {
    K <- length(var)
    Z <- length(s)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, 
        ct, phi)
    cn = 2
    dlik_dphi = 0
    
    for (z in 1:Z) {
        dmuC_dphi <- -(1 - n)/(n * cn + (1 - n) * phi)
        for (k in 1:K) {
            term1 <- c[k, z]/var[k]
            term2 <- e[k, z] * mus$C[k, z]/var[k]
            dlik_dphi <- dlik_dphi + (term1 - term2) * 
                dmuC_dphi
        }
    }
    
    # beta prior likelihood (of phi) derivative wrt to
    # phi
    dinvgamma_dphi = (-alphaP - 1)/phi + betaP/phi^2
    f = dlik_dphi + dinvgamma_dphi
    return(f)
}

estimateClonalMixWeightsParamMap <- function(rho, kappa) {
    K <- dim(rho)[1]
    pi = (rowSums(rho) + kappa - 1)/(sum(rowSums(rho)) + 
        sum(kappa) - K)
    return(pi)
}

estimateGenotypeMixWeightsParamMap <- function(rho, 
    kappa) {
    K <- dim(rho)[1]
    pi <- (rho[, 1] + kappa - 1)/(sum(rho[, 1]) + sum(kappa) - 
        K)
    return(pi)
} 
