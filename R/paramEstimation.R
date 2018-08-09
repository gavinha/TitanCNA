# author: Gavin Ha 
# 		  Dana-Farber Cancer Institute
#		  Broad Institute
# contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
# date:	  November 13, 2014

## wrapper that makes subroutine calls to
## updateParameters based on the normal estimate
## method chosen
estimateClonalCNParamsMap <- function(x, N, l, rho, rhoG, rhoZ,
    n_0, s_0, phi_0, var_0, varR_0, corRho_0, piG_0, piZ_0, chrPos_0,
    gParams, nParams, sParams, pParams, 
    maxiter = 500, normalEstimateMethod = "fixed", 
    estimateS = TRUE, estimatePloidy = TRUE, verbose = TRUE) {
    K <- dim(rho)[1]
    Z <- dim(rho)[2]
    T <- dim(rho)[3]
    # assign some variables to simpify the update equations
    a <- matrix(0, K, Z)
    b <- matrix(0, K, Z)
    c <- matrix(0, K, Z)
    d <- matrix(0, K, Z)
    e <- matrix(0, K, Z)
    f <- matrix(0, K, Z)
    g <- matrix(0, K, Z)
    h <- matrix(0, K, Z)
    for (z in 1:Z) {
        for (k in 1:K) {
            if (gParams$alleleEmissionModel == "Gaussian"){
              a[k, z] <- (x / N) %*% rho[k, z, ]
              f[k, z] <- (x / N)^2 %*% rho[k, z, ]
              h[k, z] <- ((x / N) * l) %*% rho[k, z, ]
            }else{
              a[k, z] <- x %*% rho[k, z, ]
              b[k, z] <- (N - x) %*% rho[k, z, ]
              d[k, z] <- lchoose(N, x) %*% rho[k, z, ]
            }
            c[k, z] <- l %*% rho[k, z, ]
            e[k, z] <- sum(rho[k, z, ])
            g[k, z] <- (l^2) %*% rho[k, z, ]
        }
    }
    likInitG <- sum(t(log(piG_0)) %*% rhoG[, chrPos_0])
    likInitZ <- sum(t(log(piZ_0)) %*% rhoZ[, chrPos_0])
    
    estimateOut <- updateParameters(n_0, s_0, var_0, varR_0,
        phi_0, corRho_0, chrPos_0, piG_0, piZ_0, likInitG, likInitZ, 
        a, b, c, d, e, f, g, h, 
        gParams, nParams, sParams, pParams, 
        normalEstimateMethod = normalEstimateMethod, 
        estimateS = estimateS, estimatePloidy = estimatePloidy,
        maxiter = maxiter, miniter = 10, verbose = verbose)
    ## update pi's ##
    piZ <- estimateClonalMixWeightsParamMap(rhoZ, sParams$kappaZHyper, chrPos_0)
    piG <- estimateGenotypeMixWeightsParamMap(rhoG, gParams$kappaGHyper, chrPos_0)
    estimateOut$piG <- piG
    estimateOut$piZ <- piZ
    rm(a, b, c, d, e, f, g, h)
    gc(verbose = FALSE, reset = TRUE)
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
updateParameters <- function(n_0, s_0, var_0, varR_0, phi_0, corRho_0, chrPos_0,
                             piG, piZ, likInitG, likInitZ, a, b, c, d, e, f, g, h, 
                             gParams, nParams, sParams, pParams, 
                             normalEstimateMethod, estimateS, estimatePloidy,
    maxiter = 500, miniter = 10, verbose = TRUE) {
    
    K <- length(var_0)
    Z <- length(s_0)
    s <- rep(0, Z)
    var <- rep(0, K)
    varR <- rep(0, K)
    phi <- 0
    n <- 0
    
    intervalS <- c(.Machine$double.eps, 1 - .Machine$double.eps)
    intervalVar <- c(.Machine$double.eps, 5)
    intervalVarR <- c(.Machine$double.eps, 5)
    intervalPhi <- c(.Machine$double.eps, 10)
    intervalN <- c(.Machine$double.eps, 1 - .Machine$double.eps)
    converged <- 0
    i <- 1
    maxFind <- 1
    objfun <- rep(0, maxiter)
    n_prev <- n_0
    s_prev <- s_0
    var_prev <- var_0
    varR_prev <- varR_0
    phi_prev <- phi_0
    numParameters <- K + Z + K + 2 # logR var + s + allele varR + n + phi
    F <- matrix(0, numParameters, maxiter)
    F[, i] <- c(n_prev, s_prev, var_prev, varR_prev, phi_prev)
    F0 <- rep(0, maxiter)  ## keep track of likelihood function (objective function)
    F0[i] <- likelihoodFunc(F[, i], a, b, c, d, e, f, g, h, gParams,
                            gParams$rt, gParams$rn, gParams$ct, corRho_0, piG, piZ, likInitG, likInitZ, chrPos_0,
                            sParams$alphaSHyper, sParams$betaSHyper, nParams$alphaNHyper, nParams$betaNHyper, 
                            gParams$alphaKHyper, gParams$betaKHyper, gParams$alphaRHyper,
                            gParams$betaRHyper, pParams$alphaPHyper, pParams$betaPHyper, 
                            gParams$kappaGHyper, sParams$kappaZHyper,
                            normalEstimateMethod, estimateS, estimatePloidy)
    objfun[i] <- F0[i]
    if (verbose == TRUE) {
        #message("Coord Descent[", i, "]=", format(objfun[i], 
        #    digits = 4), " n=", format(n_prev, digits = 4), 
        #    " s=[", paste(format(s_prev, digits = 2), 
        #        collapse = ","), "]", " phi=", format(phi_prev, 
        #        digits = 4))
    }
    
    ### begin coordinate descent ###
    while ((!converged && (i < maxiter)) || (i <= miniter)) {
        i <- i + 1
        if (normalEstimateMethod == "map") {
            funcN <- function(n) {
                clonalCNDerivativeNUpdateEqn(n, s_prev, 
                  var_prev, varR_prev, phi_prev, corRho_0, a, b, c, e, 
                  gParams$alleleEmissionModel, gParams$rt, 
                  gParams$rn, gParams$ct, nParams$alphaNHyper, 
                  nParams$betaNHyper)
            }
            n <- tryCatch({
        		uniroot(funcN, intervalN, tol = 1e-15)$root
		  	}, error = function(x){ 
		  		message("TitanCNA updateParameters: Issue maximizing n, using previous iteration (", n_prev, ")")
				return(n_prev) 
		  	})
        } else if (normalEstimateMethod == "fixed") {
            n <- n_prev
        }
        
        if (estimateS) {
            # estimate S independently since not using optim
            for (z in 1:Z) {
                funcS <- function(s) {
                  clonalCNDerivativeSUpdateEqn(s, n, var_prev, 
                      varR_prev, phi_prev, corRho_0, a[, z], b[, z], c[, z], e[, z], 
                      gParams$alleleEmissionModel,
                      gParams$rt, gParams$rn, gParams$ct, 
                      sParams$alphaSHyper[z], sParams$betaSHyper[z])
                }
                s[z] <- tryCatch({
        			uniroot(funcS, intervalS, tol = 1e-15)$root
				}, error = function(x){ 
					message("TitanCNA updateParameters: Issue maximizing s[,", z, "], using previous iteration (", s_prev[z], ")")
					return(s_prev[z]) 
				})
            }
        } else {
            s <- s_prev
        }
        cnLevel <- unique(gParams$ct)
        if (gParams$alleleEmissionModel == "Gaussian"){
          ## Estimate Gaussian variance for each copy number
          ## level, rather than each genotype state ##
          
          for (ck in 1:length(cnLevel)) {
            cnInd <- which(gParams$ct == cnLevel[ck])
            funcVar <- function(var){
              clonalCNDerivativeVarUpdateEqn(var, n, s, varR_prev[cnInd], phi_prev, corRho_0, 
                a[cnInd, , drop=FALSE], c[cnInd, , drop=FALSE], e[cnInd, , drop=FALSE], 
                g[cnInd, , drop=FALSE], h[cnInd, , drop=FALSE],
                gParams$rt[cnInd], gParams$rn, gParams$ct[cnInd], 
                gParams$alphaKHyper[cnInd], gParams$betaKHyper[cnInd], type = "logRatio")
            }
            var[cnInd] <- tryCatch({
            	uniroot(funcVar, intervalVar, tol = 1e-15)$root
            }, error = function(x){
            	message("TitanCNA updateParameters: Issue maximizing var[,", cnInd, "], using previous iteration (", var_prev[cnInd], ")")
            	return(var_prev[cnInd])
            })
          }
          # Estimate Gaussian variance for each allelic state 
          for (k in 1:K){
            funcVarR <- function(varR){
              clonalCNDerivativeVarUpdateEqn(varR, n, s, var[k], phi_prev, corRho_0, 
                    c[k, , drop=FALSE], a[k, , drop=FALSE], e[k, , drop=FALSE], 
                    f[k, , drop=FALSE], h[k, , drop=FALSE],
                    gParams$rt[k], gParams$rn, gParams$ct[k], 
                    gParams$alphaRHyper[k], gParams$betaRHyper[k], type = "allelicRatio")
            }
            varR[k] <- tryCatch({
            	uniroot(funcVarR, intervalVarR, tol = 1e-15)$root
            }, error = function(x){
            	message("TitanCNA updateParameters: Issue maximizing s[,", z, "], using previous iteration (", varR_prev[k], ")")
            	return(varR_prev[k])
            })
          }
        }else{ # alleleEmissionModel == "binomial"
          for (ck in 1:length(cnLevel)) {
            cnInd <- which(gParams$ct == cnLevel[ck])
            varTmp <- clonalCNDerivativeVarExactUpdateEqn(n, 
                s, phi_prev, c, e, g, 
                gParams$rt, gParams$rn, gParams$ct, cnLevel[ck],
                gParams$alphaKHyper[cnInd[1]], 
                gParams$betaKHyper[cnInd[1]])
            var[cnInd] <- varTmp
          }
          # allelic emission variance not estiamted
          varR <- varR_prev
        }
        
        if (estimatePloidy) {
            funcP <- function(phi) {
                clonalCNDerivativePloidyUpdateEqn(phi, 
                  n, s, var, varR, corRho_0, a, c, e, gParams$alleleEmissionModel,
                  gParams$rt, gParams$rn, gParams$ct, pParams$alphaPHyper, 
                  pParams$betaPHyper)
            }
            phi <- tryCatch({
        		uniroot(funcP, intervalPhi, tol = 1e-15)$root
		  	}, error = function(x){ 
		  		message("TitanCNA updateParameters: Issue maximizing phi, using previous iteration (", phi_prev, ")")
				return(phi_prev) 
		  	})
        } else {
            # estimatePloidy==FALSE
            phi <- phi_prev
        }
        
        
        F[, i] <- c(n, s, var, varR, phi)
        F0[i] <- likelihoodFunc(F[, i], a, b, c, d, e, f, g, h,
            gParams, 
            gParams$rt, gParams$rn, gParams$ct, corRho_0, 
            piG, piZ, likInitG, likInitZ, chrPos_0,
            sParams$alphaSHyper, sParams$betaSHyper, 
            nParams$alphaNHyper, nParams$betaNHyper, 
            gParams$alphaKHyper, gParams$betaKHyper, 
            gParams$alphaRHyper, gParams$betaRHyper,
            pParams$alphaPHyper, pParams$betaPHyper, 
            gParams$kappaGHyper, sParams$kappaZHyper,
            normalEstimateMethod, estimateS, estimatePloidy)
        objfun[i] <- F0[i]
        if (verbose == TRUE) {
            #message("Coord Descent[", i, "]=", format(objfun[i], 
            #    digits = 4), " n=", format(n, digits = 4), 
            #    " s=[", paste(format(s, digits = 2), 
            #      collapse = ","), "]", " phi=", format(phi, 
            #      digits = 4))
        }
        n_prev <- n
        s_prev <- s
        var_prev <- var
        varR_prev <- varR
        phi_prev <- phi
        
        if (objfun[maxFind] >= objfun[i]) {
          #if (verbose)
            #message("Coordinate descent decreases likelihood.")
        } else {
            maxFind <- i
        }
        if ((abs(objfun[i] - objfun[i - 1])/abs(objfun[i])) <= 0.001) {
            converged <- 1
        }
    }
    
    maxFind <- i
    # return estimated parameters
    output <- vector("list", 0)
    output$maxFind <- maxFind
    output$maxF <- objfun[maxFind]
    output$n <- F[1, maxFind]; sInd <- 2
    output$s <- F[sInd:(sInd + Z - 1), maxFind]; varInd <- sInd + Z
    output$var <- F[varInd:(varInd + K - 1), maxFind]; varRInd <- varInd + K
    output$varR <- F[varRInd:(varRInd + K - 1), maxFind]
    output$phi <- F[nrow(F), maxFind]
    return(output)
}


## computed once for all parameters ##
likelihoodFunc <- function(unk, a, b, c, d, e, f, g, h, gParams, 
                           rt, rn, ct, corRho, piG, piZ, likInitG, likInitZ, chrPos_0, 
                           alphaS, betaS, alphaN, betaN, 
                           alphaK, betaK, alphaR, betaR, alphaP, betaP, 
                           kappaG, kappaZ,
                           normalEstimateMethod, estimateS, estimatePloidy) {
    K <- length(rt)
    Z <- length(alphaS)

    # c(n_prev, s_prev, var_prev, varR_prev, phi_prev)
    n <- unk[1]; sInd <- 2
    s <- unk[sInd:(sInd + Z - 1)]; varInd <- sInd + Z
    var <- unk[varInd:(varInd + K - 1)]; varRInd <- varInd + K
    varR <- unk[varRInd:(varRInd + K - 1)]
    phi <- unk[length(unk)]
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    likObs <- 0
    
    ## find copy number level indices - short list of
    ## the full genotype indices
    cnLevel <- unique(ct)
    cnInd <- sapply(cnLevel, function(x) {
      which(ct == x)[1]
    })
    
    ## emission / observed likelihood
    if (gParams$alleleEmissionModel == "Gaussian"){
      for (z in 1:Z) {
        for (k in 1:K) {
            #e[k, z] * log(1/sqrt(2 * pi * varR[k])) - 
            #(f[k, z] - 2 * a[k, z] * mus$R[k, z] + e[k, z] * mus$R[k, z]^2) / (2 * varR[k]) + # allelic
            #e[k, z] * log(1/sqrt(2 * pi * var[k])) - 
            #(g[k, z] - 2 * c[k, z] * mus$C[k, z] + e[k, z] * mus$C[k, z]^2) / (2 * var[k]) # logR
            ## bivariate normal
            term1 <- e[k, z] * log(1 / 2 * pi * sqrt(var[k] * varR[k] * (1 - corRho^2))) 
            term2 <- (1 / (2 * (1 - corRho^2))) 
            term3 <- (f[k, z] - 2 * a[k, z] * mus$R[k, z] + e[k, z] * mus$R[k, z]^2) / (2 * varR[k]) 
            term4 <- (g[k, z] - 2 * c[k, z] * mus$C[k, z] + e[k, z] * mus$C[k, z]^2) / (2 * var[k]) 
            term5 <- 2 * corRho * (h[k, z] - a[k, z] * mus$C[k, z] - c[k, z] * mus$R[k, z] + 
                                     e[k, z] * mus$C[k, z] * mus$R[k, z]) / sqrt(var[k] * varR[k])
            likObs <- likObs + term1 - term2 * (term3 + term4 - term5)
        }
      }
    }else{
      for (z in 1:Z) {
        for (k in 1:K) {
          likObs <- likObs + 
            d[k, z] + a[k, z] * log(mus$R[k, z]) + b[k, z] * log(1 - mus$R[k, z]) + # binomial, allelic
            e[k, z] * log(1/sqrt(2 * pi * var[k])) - 
            (g[k, z] - 2 * c[k, z] * mus$C[k, z] + e[k, z] * mus$C[k, z]^2) / (2 * var[k]) # gaussian, logR
        }
      }
    }
    
    ## initial state likelihood ##
    #likInitG <- sum(log(piG) %*% rhoG[, chrPos_0])
    #likInitZ <- sum(log(piZ) %*% rhoZ[, chrPos_0])
    likTxn <- 0
    #for (ks in 1:KS){
    #  likTxn <- likTxn + Zcounts[ks, ] %*% log(A[ks, ])
    #}    
    

    prior <- priorProbs(n, s, phi, var, varR, piG, piZ, gParams, 
                        normalEstimateMethod, estimateS, estimatePloidy,
                        alphaN, betaN, alphaS, betaS, alphaP, betaP,
                        alphaK, betaK, alphaR, betaR, kappaG, kappaZ)
    f <-  likObs + prior$prior + likInitG + likInitZ + likTxn 
    return(f)
}


clonalCNDerivativeNUpdateEqn <- function(n, s, var, varR, phi, corRho,
    a, b, c, e, alleleEmissionModel, rt, rn, ct, alphaN, betaN) {
    K <- length(ct)
    Z <- length(s)
    cn <- 2
    # data likelihood derivative wrt to n
    dlik_dn <- 0
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    lik <- 0
    for (z in 1:Z) {
        dmuC_dn <- (cn - s[z] * cn - (1 - s[z]) * ct) / 
          (n * cn + (1 - n) * s[z] * cn + (1 - n) * (1 - s[z]) * ct) - 
          (cn - phi)/(n * cn + (1 - n) * phi)
        dmuR_dn <- ((1 - s[z]) * cn * ct * (rn - rt)) / 
            ((n * cn + (1 - n) * s[z] * cn + (1 - n) * (1 - s[z]) * ct)^2)
        for (k in 1:K) {
          # allelic
          if (alleleEmissionModel == "Gaussian"){
            term0 <- 1 / (1 - corRho^2)
            term1 <- (a[k, z] - e[k, z] * mus$R[k, z]) * dmuR_dn[k] / varR[k]
            term2 <- (c[k, z] - e[k, z] * mus$C[k, z]) * dmuC_dn[k] / var[k] 
            term3 <- corRho * (((c[k, z] - e[k, z] * mus$C[k, z]) * dmuR_dn[k]) +
                               ((a[k, z] - e[k, z] * mus$R[k, z]) * dmuC_dn[k])) / 
                     sqrt(varR[k] * var[k])
            lik <- term0 * (term1 + term2 - term3)
          }else{ # alleleEmissionModel == "binomial"
            term1 <- a[k, z] * dmuR_dn[k]/mus$R[k, z] -
              b[k, z] * dmuR_dn[k]/(1 - mus$R[k, z]) #binomial
            term2 <- (c[k, z] - e[k, z] * mus$C[k, z]) * dmuC_dn[k]/var[k] # logR
            lik <- term1 + term2
          }
         
          dlik_dn <- dlik_dn + lik
        }
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dbeta_dn <- (alphaN - 1)/n - (betaN - 1)/(1 - n)
    f <- dlik_dn + dbeta_dn
    return(f)
}

clonalCNDerivativeSUpdateEqn <- function(s, n, var, varR, phi, corRho,
    a, b, c, e, alleleEmissionModel, rt, rn, ct, alphaS, betaS) {
    K <- length(ct)
    cn <- 2
    # data likelihood derivative wrt to s
    dlik_ds <- 0
    lik <- 0
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    dmuC_ds <- ((1 - n) * cn - (1 - n) * ct) / 
      (n * cn + (1 - n) * s * cn + (1 - n) * (1 - s) * ct)
    dmuR_ds <- ((1 - n) * cn * ct * (rn - rt)) / 
      ((n * cn + (1 - n) * s * cn + (1 - n) * (1 - s) * ct)^2)
    
    for (k in 1:K) {
      if (alleleEmissionModel == "Gaussian"){
        term0 <- 1 / (1 - corRho^2)
        term1 <- (a[k] - e[k] * mus$R[k]) * dmuR_ds[k] / varR[k]
        term2 <- (c[k] - e[k] * mus$C[k]) * dmuC_ds[k] / var[k] 
        term3 <- corRho * (((c[k] - e[k] * mus$C[k]) * dmuR_ds[k]) +
                           ((a[k] - e[k] * mus$R[k]) * dmuC_ds[k])) / 
                 sqrt(varR[k] * var[k])
        lik <- term0 * (term1 + term2 - term3)
      }else{ # binomial
        term1 <- a[k] * dmuR_ds[k]/mus$R[k] -
          b[k] * dmuR_ds[k]/(1 - mus$R[k])
        term2 <- (c[k] - e[k] * mus$C[k]) * dmuC_ds[k] / var[k]
        lik <- term1 + term2
      }
      
      dlik_ds <- dlik_ds + lik
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dbeta_ds <- (alphaS - 1)/s - (betaS - 1)/(1 - s)
    f <- dlik_ds + dbeta_ds
    return(f)
}

# ck is the copy number level that we want to
# marginalize over ct is vector of copy number;
# each corresponding to genotype state
clonalCNDerivativeVarExactUpdateEqn <- function(n, s, phi, c, e, g, rt, rn, ct, ck, 
                                                alphaK, betaK) {
    Z <- length(s)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    cnInd <- which(ct == ck)
    var <- 0
    term1 <- 0
    term2 <- 0
    for (k in cnInd) {
      for (z in 1:Z) {
        term1 = term1 + -sum(g[k, z] - 2 * c[k, z] * 
                               mus$C[k, z] + e[k, z] * (mus$C[k, z])^2)
        term2 = term2 + -sum(e[k, z])
      }
    }
    term1 = term1 - 2 * betaK
    term2 = term2 - 2 * (alphaK + 1)
    var = term1/term2
    return(var)
}

clonalCNDerivativeVarUpdateEqn <- function(var, n, s, var2, phi, corRho, 
    a, c, e, g, h, rt, rn, ct, alphaK, betaK, type) {
    Z <- length(s)
    K <- length(alphaK)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    if (type == "logRatio"){
      mus1 <- mus$C
      mus2 <- mus$R
    }else{  #type=="allelicRatio"
      mus1 <- mus$R
      mus2 <- mus$C
    }
    # data likelihood derivative wrt to var
    dlik_dvar <- 0
    lik <- 0
    for (z in 1:Z) {
      for (k in 1:K){
        #term1 <- e[z] * var2 * (1 - corRho^2) / (2 * sqrt(var * var2 * (1 - corRho)))
        term1 <- e[k, z] / (2 * var)  
        term2 <- 1 / (2 * (1 - corRho^2))
        term3 <- (g[k, z] - 2 * c[k, z] * mus1[k, z] + e[k, z] * mus1[k, z]^2) / (var^2)
        #term4 <- corRho * (a[z] - e[z] * mus2[z]) * (c[z] - e[z] * mus1[z]) / (var2^(1/2) * var^(3/2))
        #term4 <- corRho * (h[z] - c[z] * mus2[z] - a[z] * mus1[z] + e[z] * mus1[z] * mus2[z]) / (var2^(1/2) * var^(3/2))
        term4 <- corRho * (h[k, z] - c[k, z] * mus2[k, z] - a[k, z] * mus1[k, z] + e[k, z] * mus1[k, z] * mus2[k, z]) / 
          (var * sqrt(var * var2[k]))
        dlik_dvar = dlik_dvar + -term1 + (term2 * (term3 - term4))
      }
    }
    
    # beta prior likelihood (of s) derivative wrt to s
    dinvgamma_dvar = (-alphaK - 1)/var + betaK/var^2
    f = sum(dlik_dvar) + sum(dinvgamma_dvar)
    #message("var=", var, " f=",f)
    return(f)
}

clonalCNDerivativePloidyUpdateEqn <- function(phi, n, s, var, varR, corRho,
    a, c, e, alleleEmissionModel, rt, rn, ct, alphaP, betaP) {
    K <- length(var)
    Z <- length(s)
    mus <- clonalTwoComponentMixtureCN(rt, rn, n, s, ct, phi)
    cn = 2
    dlik_dphi = 0
    lik <- 0
    dmuC_dphi <- -(1 - n)/(n * cn + (1 - n) * phi)
    for (z in 1:Z) {
        for (k in 1:K) {
          if (alleleEmissionModel == "Gaussian"){
            term0 <- 1 / (1 - corRho^2)
            term1 <- (c[k] - e[k] * mus$C[k]) * dmuC_dphi / var[k]
            term2 <- corRho * ((a[k] - e[k] * mus$R[k]) * dmuC_dphi) / 
                     sqrt(varR[k] * var[k])
            lik <- term0 * (term1 - term2)
          }else{ #alleleEmissionModel == "binomial"
            term1 <- c[k, z]/var[k]
            term2 <- e[k, z] * mus$C[k, z]/var[k]
            lik <- (term1 - term2) * dmuC_dphi
          }
          dlik_dphi <- dlik_dphi + lik
        }
    }
    
    # beta prior likelihood (of phi) derivative wrt to
    # phi
    dinvgamma_dphi = (-alphaP - 1)/phi + betaP/phi^2
    f = dlik_dphi + dinvgamma_dphi
    return(f)
}


estimateClonalMixWeightsParamMap <- function(rho, kappa, chrPos_0) {
    K <- nrow(rho)
    #pi = (rowSums(rho[, chrPos_0]) + kappa - 1) / 
    #  (sum(rowSums(rho[, chrPos_0])) + sum(kappa) - K)
    pi <- (sum(rho[, chrPos_0]) + kappa - 1) / 
      (sum(rho[, chrPos_0]) + sum(kappa) - K)
    return(pi)
}

estimateGenotypeMixWeightsParamMap <- function(rho, kappa, chrPos_0) {
    K <-nrow(rho)
    pi <- (sum(rho[, chrPos_0]) + kappa - 1) / 
      (sum(rho[, chrPos_0]) + sum(kappa) - K)
    return(pi)
} 
