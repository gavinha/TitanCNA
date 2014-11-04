# file : hmmClonal.R author: Gavin Ha
# <gha@bccrc.ca> Dept of Molecular Oncolgy British
# Columbia Cancer Agency University of British
# Columbia date : March 3, 2014

#### EM (FWD-BACK) Algorithm ####
runEMclonalCN <- function(data, gParams, nParams, pParams, 
    sParams, txnExpLen = 1e+15, txnZstrength = 5e+5, 
    maxiter = 15, maxiterUpdate = 1500, pseudoCounts = 1e-300, 
    normalEstimateMethod = "map", estimateS = TRUE, 
    estimatePloidy = TRUE, useOutlierState = FALSE, 
    verbose = TRUE) {
    ## check that arguments contain necessary list
    ## elements
    if (length(data$chr) == 0 || length(data$posn) == 
        0 || length(data$ref) == 0 || length(data$tumDepth) == 
        0 || length(data$logR) == 0) {
        stop("data must contain named list elements: chr, posn, ref, tumDepth, 
             logR. See loadDataFromFile()")
    }
    if (length(gParams$rt) == 0 || length(gParams$ct) == 
        0 || length(gParams$rn) == 0 || length(gParams$ZS) == 
        0 || length(gParams$var_0) == 0 || length(gParams$alphaKHyper) == 
        0 || length(gParams$betaKHyper) == 0 || length(gParams$kappaGHyper) == 
        0) {
        stop("genotypeParams must contain named list elements: rt, rn, ct, ZS, 
             var_0, alphaKHyper, betaKHyper, kappaGHyper.\n         
             See loadDefaultParameters().")
    }
    if (length(pParams$phi_0) == 0 || length(pParams$alphaPHyper) == 
        0 || length(pParams$betaPHyper) == 0) {
        stop("ploidyParams must contain named list elements: phi_0, alphaPHyper, 
             betaPHyper. See loadDefaultParameters()")
    }
    if (length(nParams$n_0) == 0 || length(nParams$alphaNHyper) == 
        0 || length(nParams$betaNHyper) == 0) {
        stop("normalParams must contain named list elements: n_0, alphaNHyper, 
             betaNHyper. See loadDefaultParameters()")
    }
    if (length(sParams$s_0) == 0 || length(sParams$kappaZHyper) == 
        0 || length(sParams$alphaSHyper) == 0 || length(sParams$betaSHyper) == 
        0) {
        stop("sParams (clonal parameters) must contain named list elements: s_0, 
             kappaSHyper, alphaSHyper, betaSHyper. \n         
             See setupClonalParameters().")
    }
    
    ## track the total time
    ticTotalId <- proc.time()
    if (verbose == TRUE) {
        message("titan: Running HMM...")
    }
    
    ## requirements for parallelization require(foreach)
    
    ## Allocate memory for our variables
    Z <- length(sParams$s_0)  # number of clonal clusters
    K <- length(gParams$rt)  # number of genotype states + outlier state
    N <- length(data$ref)  # number of data points
    if (useOutlierState) {
        O <- 1
        gParams$rt <- c(0.5, gParams$rt)
        gParams$ct <- c(0, gParams$ct)
        gParams$ZS <- c(-1, gParams$ZS)
        gParams$var_0 <- c(gParams$outlierVar, gParams$var_0)
        gParams$kappaGHyper <- c(2, gParams$kappaGHyper)
        gParams$alphaKHyper <- c(0, gParams$alphaKHyper)
        gParams$betaKHyper <- c(0, gParams$betaKHyper)
        K <- length(gParams$rt)
        gNoOUTStateParams <- excludeGarbageState(gParams, 
            K)
        KnoOutlier <- K - 1
        kRange <- 2:K
        Ktotal <- KnoOutlier * Z + 1
        KtotalRange <- 2:Ktotal
    } else {
        O <- 0
        kRange <- 1:K
        KnoOutlier <- K
        Ktotal <- KnoOutlier * Z
        KtotalRange <- 1:Ktotal
        gNoOUTStateParams <- gParams
    }
    maxiter <- maxiter + 1
    mus <- vector("list", 0)
    mus$R <- array(0, dim = c(KnoOutlier, Z, maxiter))  # state Binomial means: KxZ; don't include outlier
    mus$C <- array(0, dim = c(KnoOutlier, Z, maxiter))  # state Binomial means: KxZ; don't include outlier
    piG <- matrix(0, K, maxiter)  # initial state distribution of genotypes + outlier state
    piZ <- matrix(0, Z, maxiter)  # initial state assignments for clonal clusters
    s <- matrix(0, Z, maxiter)  # clonal frequency parameter
    var <- matrix(0, K, maxiter)  # clonal frequency parameter + outlier state
    phi <- rep(0, maxiter)  # global ploidy parameter
    n <- rep(0, maxiter)  # global normal contamination parameter
    loglik <- rep(0, maxiter)  #log likelihood
    rho <- NULL
    rhoZ <- NULL
    rhoG <- NULL
    outRhoRow <- NULL
    fwdBackPar <- NULL  #marginal probs or responsibilities
    
    ## Set up ## set up the chromosome indicies and make
    ## cell array of chromosome indicies
    chrs <- unique(data$chr)
    posn <- data$posn #assign positions to a new variable
    numChrs <- length(chrs)
    chrsI <- vector("list", numChrs)
    piZi <- vector("list", numChrs)
    piGi <- vector("list", numChrs)
    # initialise the chromosome index and the init
    # state distributions
    for (i in 1:numChrs) {
        chrsI[[i]] <- which(data$chr == chrs[i])
    }
    
    ## INITIALIZATION ##
    piG_0 <- estimateDirichletParamsMap(gParams$kappaGHyper)  #add the outlier state
    piZ_0 <- estimateDirichletParamsMap(sParams$kappaZHyper)
    i <- 1
    loglik[i] <- -Inf
    converged <- 0  # flag for convergence
    s[, i] <- sParams$s_0
    var[, i] <- gParams$var_0
    # var[1,i] <- gParams$outlierVar
    phi[i] <- pParams$phi_0
    n[i] <- nParams$n_0
    piG[, i] <- piG_0
    if (useOutlierState) 
        {
            piG[1, i] <- piG_0[1] * piZ_0[1]
        }  #initialize outlier state 
    piZ[, i] <- piZ_0
    musTmp <- clonalTwoComponentMixtureCN(gNoOUTStateParams$rt, 
        gParams$rn, nParams$n_0, sParams$s_0, gNoOUTStateParams$ct, 
        pParams$phi_0)
    mus$R[, , i] <- musTmp$R
    mus$C[, , i] <- musTmp$C
    
    # compute likelihood conditional on k and z
    # (K-by-Z-by-N)
    pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
        matrix(mus$R[, , i], KnoOutlier, Z))
    pyC <- computeNormalObslik(data$logR, matrix(mus$C[, 
        , i], KnoOutlier, Z), gNoOUTStateParams$var_0)
    if (useOutlierState == 1) {
        pyO <- outlierObslik(data$ref, data$tumDepth, 
            data$logR, gParams$outlierVar)
        py <- exp(rbind(pyO$R, pyR) 
        			+ rbind(pyO$C, pyC)) + pseudoCounts  #add the outlier state
    } else {
        #joint likelihood between Binomial and Gaussian
        py <- pyR + pyC + pseudoCounts
    }
        
    ## EXPECTATION MAXIMIZATION
    while (!converged && (i < (maxiter))) {
        # clear the previous iteration of rho and garbage
        # collect
        rm(fwdBackPar, rhoZ, rhoG, musTmp, pyR, pyC)
        gc(verbose = FALSE, reset = TRUE)
        ticId <- proc.time()
        i <- i + 1
        ## E-STEP: FORWARDS-BACKWARDS ALGORITHM;
        ## PARALLELIZATION
        loglik[i] <- 0
        if (verbose == TRUE) {
            message("fwdBack: Iteration ", i - 1, " chr: ", 
                appendLF = FALSE)
        }
        # cFun <- getLoadedDLLs()[['TitanCNA']][[2]]
        
        piGiZi <- matrix(0, numChrs, Ktotal)
        for (c in 1:numChrs) {
            piZi[[c]] <- piZ[, i - 1, drop = FALSE]
            piGi[[c]] <- piG[kRange, i - 1, drop = FALSE]
            piGiZi[c, KtotalRange] <- as.vector(piGi[[c]] %*% 
                t(piZi[[c]]))  #inner product then flatten to 1-by-K*Z
            if (useOutlierState) 
                {
                  piGiZi[c, ] <- c(piG[1, i - 1], piGiZi[c, 
                    KtotalRange])
                }  #add outlier state
        }
        gc(verbose = FALSE, reset = TRUE)
        ## PARALLELIZATION
        fwdBackPar <- foreach(c = 1:numChrs, .combine = rbind, .noexport = c("bigPy","data")) %dopar% 
            {
                ## Fwd-back returns rho (responsibilities) and
                ## loglik (log-likelihood, p(Data|Params)) include outlier state
                if (verbose == TRUE) {
                  message(c, " ", appendLF = FALSE)
                }
                .Call("fwd_backC_clonalCN", 
                  log(piGiZi[c, ]), py[, chrsI[[c]]], gNoOUTStateParams$ct, 
                  gNoOUTStateParams$ZS, Z, posn[chrsI[[c]]], 
                  txnZstrength * txnExpLen, txnExpLen, O)
            }
        if (verbose == TRUE) {
            message("")
        }
        if (numChrs > 1) {
            loglik[i] <- sum(do.call(rbind, fwdBackPar[, 
                2]))  #combine loglik
            # marginal probs or responsibilities,
            # p(G_t,Z_t|Data,Params)
            rho <- do.call(cbind, fwdBackPar[, 1])  #combine rho from parallel runs  
        } else {
            loglik[i] <- sum(do.call(rbind, fwdBackPar[2]))
            rho <- do.call(cbind, fwdBackPar[1])
        }
        
        outRhoRow <- rho[1, ]
        rho <- exp(array(rho[KtotalRange, ], dim = c(KnoOutlier, Z, N))) 
        ##only use states that we will need to for estimation later
        # marginalize to get rhoZ and rhoG from rho
        rhoZ <- colSums(rho)
        rhoG <- colSums(aperm(rho, c(2, 1, 3)))
        if (useOutlierState) {
            rhoG <- rbind(outRhoRow, rhoG)
        }
        
        ## M-STEP: COORDINATE DESCENT UPDATE OF PARAMETERS
        estimateOut <- estimateClonalCNParamsMap(data$ref, data$tumDepth, 
            data$logR, rho, n[i - 1], 
            s[, i - 1], var[kRange, i - 1], phi[i - 1], 
            gNoOUTStateParams, nParams, sParams, pParams, 
            maxiter = maxiterUpdate, 
            normalEstimateMethod = normalEstimateMethod, 
            estimateS = estimateS, estimatePloidy = estimatePloidy, 
            verbose = verbose)
        n[i] <- estimateOut$n
        s[, i] <- estimateOut$s
        var[kRange, i] <- estimateOut$var
        phi[i] <- estimateOut$phi
        if (useOutlierState) 
            {
                var[1, i] <- gParams$outlierVar
            }  #garbage state 
        piZ[, i] <- estimateClonalMixWeightsParamMap(rhoZ, 
            sParams$kappaZHyper)
        piG[, i] <- estimateGenotypeMixWeightsParamMap(rhoG, 
            gParams$kappaGHyper)
        rm(rho, estimateOut, outRhoRow)  #clear rho and estimateOut
        ## Recompute the likelihood conditional on k and z
        musTmp <- clonalTwoComponentMixtureCN(gNoOUTStateParams$rt, 
            gParams$rn, n[i], s[, i], gNoOUTStateParams$ct, 
            phi[i])
        mus$R[, , i] <- musTmp$R
        mus$C[, , i] <- musTmp$C
        pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
            matrix(mus$R[, , i], KnoOutlier, Z))
        pyC <- computeNormalObslik(data$logR, matrix(mus$C[, 
            , i], KnoOutlier, Z), var[kRange, i])
        if (useOutlierState == 1) {
            pyO <- outlierObslik(data$ref, data$tumDepth, 
                data$logR, gParams$outlierVar)
            py <- exp(rbind(pyO$R, pyR) 
            		+ rbind(pyO$C, pyC)) + pseudoCounts  #add the outlier state
        } else {
            #py <- exp(pyR + pyC) + pseudoCounts  #joint likelihood between Binomial and Gaussian
            py <- pyR + pyC + pseudoCounts
        }
        
        ## Compute log-likelihood and check converge
        priorS <- rep(0, Z)
        if (estimateS) {
            for (z in 1:Z) {
                priorS[z] <- betapdflog(s[z, i], sParams$alphaSHyper[z], 
                  sParams$betaSHyper[z])
            }
        } else {
            priorS[1:Z] <- 0
        }
        cnLevel <- unique(gParams$ct)
        priorVar <- rep(0, length(cnLevel))
        cnInd <- sapply(cnLevel, function(x) {
            which(gParams$ct == x)[1]
        })
        for (k in cnInd) {
            priorVar[k] <- invgammapdflog(var[k, i], 
                gParams$alphaKHyper[k], gParams$betaKHyper[k])
        }
        if (estimatePloidy) {
            priorPhi <- invgammapdflog(phi[i], pParams$alphaPHyper, 
                pParams$betaPHyper)
        } else {
            priorPhi <- 0
        }
        if (normalEstimateMethod != "fixed") {
            priorN <- betapdflog(n[i], nParams$alphaNHyper, 
                nParams$betaNHyper)
        } else {
            priorN <- 0
        }
        priorPiG <- dirichletpdflog(piG[, i], gParams$kappaGHyper)
        priorPiZ <- dirichletpdflog(piZ[, i], sParams$kappaZHyper)
        if (verbose == TRUE) {
            message("fwdBack: loglik=", sprintf("%0.4f", 
                loglik[i]))
        }
        loglik[i] <- loglik[i] + priorPiG + priorPiZ + 
            sum(priorS, na.rm = TRUE) + sum(priorVar, 
            na.rm = TRUE) + priorPhi + priorN
        if (verbose == TRUE) {
            message("fwdBack: priorN=", sprintf("%0.4f", 
                priorN))
            message("fwdBack: priorS=", sprintf("%0.4f", 
                sum(priorS, na.rm = TRUE)))
            message("fwdBack: priorVar=", sprintf("%0.4f", 
                sum(priorVar, na.rm = TRUE)))
            message("fwdBack: priorPhi=", sprintf("%0.4f", 
                priorPhi))
            message("fwdBack: priorPiG=", sprintf("%0.4f", 
                priorPiG))
            message("fwdBack: priorPiZ=", sprintf("%0.4f", 
                priorPiZ))
            message("fwdBack: EM iteration ", i - 1, 
                " complete loglik=", sprintf("%0.4f", 
                  loglik[i]))
        }
        
        if (((abs(loglik[i] - loglik[i - 1])/abs(loglik[i])) <= 
            0.001) && (loglik[i] >= loglik[i - 1])) {
            converged <- 1
        } else if (loglik[i] < loglik[i - 1]) {
            # stop('Failed EM!')
            converged <- 1
            message("fwdBack: Optimization during update decreased complete 
                    likelihood.  Stopping EM...")
        }
        
        elapsedTime <- (proc.time() - ticId)/60
        if (verbose == TRUE) {
            message("fwdBack: Elapsed time for iteration ", 
                i - 1, ": ", sprintf("%0.4f", elapsedTime[3]), 
                "m")
        }
        gc(verbose = FALSE, reset = TRUE)
    }
    
    ## print time elapsed to screen if(verbose==TRUE){
    ## message('s=\t',sprintf('%0.4f ',s[,dim(s)[2]]))
    ## }
    elapsedTimeTotal <- (proc.time() - ticTotalId)/60
    if (verbose == TRUE) {
        message("fwdBack: Total elapsed time: ", sprintf("%0.4f", 
            elapsedTimeTotal[3]), "m")
    }
    
    output <- vector("list", 0)
    output$n <- n[1:i]
    output$s <- s[, 1:i, drop = FALSE]
    output$var <- var[, 1:i, drop = FALSE]
    output$phi <- phi[1:i]
    output$piG <- piG[, 1:i, drop = FALSE]
    output$piZ <- piZ[, 1:i, drop = FALSE]
    output$muR <- mus$R[, , 1:i, drop = FALSE]
    output$muC <- mus$C[, , 1:i, drop = FALSE]
    output$loglik <- loglik[1:i]
    # output$rho <- rho
    output$rhoG <- rhoG
    output$rhoZ <- rhoZ
    # store parameters for use with viterbi later
    output$txnExpLen <- txnExpLen
    output$txnZstrength <- txnZstrength
    output$useOutlierState <- useOutlierState
    output$genotypeParams <- gNoOUTStateParams
    output$ploidyParams <- pParams
    output$normalParams <- nParams
    output$clonalParams <- sParams
    output$symmetric <- gParams$symmetric
    return(output)
}

viterbiClonalCN <- function(data, convergeParams, genotypeParams = NULL) {
    ## requirements for parallelization require(foreach)
    ## use genotypeParams found in convergeParams unless
    ## genotypeParams given
    if (is.null(genotypeParams)) {
        if (length(convergeParams$genotypeParams) > 
            0) {
            genotypeParams <- convergeParams$genotypeParams
        } else {
            stop("convergeParams does not contain list element: genotypeParams. 
                 Please specify one to viterbiClonalCN() or add the element to 
                 convergeParams.")
        }
    }
    # if convergeParams contains txnExpLen,
    # txnZstrength, useOutlier, then use that unless
    # given
    if (length(convergeParams$txnExpLen) == 0 || length(convergeParams$txnZstrength) == 
        0 || length(convergeParams$useOutlierState) == 
        0) {
        stop("convergeParams does not contain list elements: txnExpLen, 
             txnZstrength and/or useOutlierState. ")
    }
    txnExpLen <- convergeParams$txnExpLen
    txnZstrength <- convergeParams$txnZstrength
    useOutlierState <- convergeParams$useOutlierState
    
    K <- dim(convergeParams$var)[1]
    Z <- dim(convergeParams$s)[1]
    T <- length(data$chr)
    i <- dim(convergeParams$s)[2]
    
    if (useOutlierState) {
        O <- 1
        KnoOutlier <- K - 1
        kRange <- 2:K
        Ktotal <- KnoOutlier * Z + 1
        KtotalRange <- 2:Ktotal
        if (dim(convergeParams$var)[1] == length(genotypeParams$var_0)) {
            stop("convergeParams does not contain the outlier state but 
                 useOutlierState is set to TRUE.")
        }
    } else {
        O <- 0
        kRange <- 1:K
        KnoOutlier <- K
        Ktotal <- KnoOutlier * Z
        KtotalRange <- 1:Ktotal
        if (dim(convergeParams$var)[1] == (length(genotypeParams$var_0) + 
            1)) {
            stop("convergeParams contains the outlier state but useOutlierState 
                 is set to FALSE.")
        }
    }
    
    chrs <- unique(data$chr)
    numChrs <- length(chrs)
    chrsI <- vector("list", numChrs)
    piZi <- vector("list", numChrs)
    piGi <- vector("list", numChrs)
    # initialise the chromosome index
    for (j in 1:numChrs) {
        chrsI[[j]] <- which(data$chr == chrs[j])
    }
    
    ## compute likelihood
    pseudoCounts <- 1e-200
    pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
        matrix(convergeParams$muR[, , i], KnoOutlier, 
            Z))
    pyC <- computeNormalObslik(data$logR, matrix(convergeParams$muC[, 
        , i], KnoOutlier, Z), convergeParams$var[kRange, 
        i])
    if (useOutlierState) {
        pyO <- outlierObslik(data$ref, data$tumDepth, 
            data$logR, genotypeParams$outlierVar)
        py <- rbind(pyO$R, pyR) 
        		+ rbind(pyO$C, pyC) + pseudoCounts  #add the outlier state
    } else {
        py <- pyR + pyC + pseudoCounts  #joint likelihood between Binomial and Gaussian
    }
    
    piGiZi <- matrix(0, numChrs, Ktotal)
    for (c in 1:numChrs) {
        piZi[[c]] <- convergeParams$piZ[, i - 1, drop = FALSE]
        piGi[[c]] <- convergeParams$piG[kRange, i - 
            1, drop = FALSE]
        piGiZi[c, KtotalRange] <- as.vector(piGi[[c]] %*% 
            t(piZi[[c]]))  #inner product then flatten to 1-by-K*Z
        if (useOutlierState) 
            {
                piGiZi[c, ] <- c(convergeParams$piG[1, 
                  i - 1], piGiZi[c, KtotalRange])
            }  #add outlier state
    }
    G <- foreach(c = 1:numChrs, .combine = c) %dopar% 
        {
            # dyn.load(cFun)
            viterbiOut <- .Call("viterbiC_clonalCN", 
                log(piGiZi[c, ]), py[, chrsI[[c]]], 
                genotypeParams$ct, genotypeParams$ZS, 
                Z, data$posn[chrsI[[c]]], 
                txnZstrength * txnExpLen, txnExpLen, O)
        }
    return(G)
}


# Computes the binomial observation likelihood
# given discrete reference read counts and total
# depth.  mus is K-by-Z
computeBinomialObslik <- function(ref, depth, mus) {
    N <- length(ref)
    K <- dim(mus)[1]
    Z <- dim(mus)[2]
    py <- matrix(0, K * Z, N)  #local evidence is 2-D matrix: K*ZxN
    for (z in 1:Z) {
        for (k in 1:K) {
            # if (k==1){ py[k+(z-1)*K,] <-
            # log(dunif(depth,min=0,max=depth)) }else{
            py[k + (z - 1) * K, ] <- binomialpdflog(ref, 
                depth, mus[k, z])
            # }
        }
    }
    # py <- t(matrix(py,N,K*Z))
    return(py)
}

# Computes the Gaussian observation likelihood
# given log ratio data.  mus is K-by-Z
computeNormalObslik <- function(x, mus, var) {
    N <- length(x)
    K <- dim(mus)[1]
    Z <- dim(mus)[2]
    py <- matrix(0, K * Z, N)  #local evidence is 2-D matrix: K*ZxN
    for (z in 1:Z) {
        for (k in 1:K) {
            py[k + (z - 1) * K, ] <- normalpdflog(x, 
                mus[k, z], var[k])
        }
    }
    # py <- t(matrix(py,N,K*Z))
    return(py)
}

# 2-component mixtures for the Binomial mean given
# stromal contamination parameter (s) and Gaussian
# mean given s and ploidy (phi)
clonalTwoComponentMixtureCN <- function(rt, rn, n, 
    s, ct, phi) {
    K <- length(rt)
    Z <- length(s)
    musR <- matrix(NA, K, Z)
    musC <- matrix(NA, K, Z)
    cn <- 2
    
    for (k in 1:K) {
        for (z in 1:Z) {
            numRefAlleles <- n * rn * cn + (1 - n) * 
                s[z] * rn * cn + (1 - n) * (1 - s[z]) * 
                rt[k] * ct[k]  #number of reference alleles based on copy number
            totalAlleles <- n * cn + (1 - n) * s[z] * 
                cn + (1 - n) * (1 - s[z]) * ct[k]  #total number of alleles based on copy number
            totalPloidy <- n * cn + (1 - n) * phi  #total ploidy of tumour sample (includes normal + tumour)
            musR[k, z] <- numRefAlleles/totalAlleles
            musC[k, z] <- log(totalAlleles/totalPloidy)
        }
    }
    output <- vector("list", 0)
    output$R <- musR
    output$C <- musC
    return(output)
}

outlierObslik <- function(ref, depth, logR, var) {
    outR <- log(dunif(depth, min = 0, max = depth))
    # outR <- betabinomialpdflog(ref,depth,0.5,3)
    outC <- normalpdflog(logR, 0, var)
    output <- vector("list", 0)
    output$R <- outR
    output$C <- outC
    return(output)
}

# Binomial prob density function computed and
# returned in log scale
binomialpdflog <- function(k, N, mu) {
    c <- lgamma(N + 1) - lgamma(k + 1) - lgamma(N - 
        k + 1)  #normalizing constant  
    l <- k * log(mu) + (N - k) * log(1 - mu)  #likelihood
    y <- c + l  #together
    # y <- exp(y)
    return(y)
}

# Gaussian prob density function computed and
# returned in log scale
normalpdflog <- function(x, mu, var) {
    c <- log(1/(sqrt(var) * sqrt(2 * pi)))  #normalizing constant
    l <- -((x - mu)^2)/(2 * var)  #likelihood
    y <- c + l  #together 
    # y <- exp(y)
    return(y)
}

betabinomialpdflog <- function(k, N, mu, M) {
    theta <- k/N
    c <- lgamma(M)/(lgamma(mu * M) * lgamma(1 - mu))
    l <- (M * mu - 1) * log(theta) + (M * (1 - mu) - 
        1) * log(1 - theta)
    l[is.infinite(l)] <- -1 * .Machine$double.xmin
    y <- c + l
    return(y)
}

betapdflog <- function(x, a, b) {
    y = -lbeta(a, b) + (a - 1) * log(x) + (b - 1) * 
        log(1 - x)
    return(y)
}

invgammapdflog <- function(x, a, b) {
    c <- a * log(b) - lgamma(a)  # normalizing constant  
    l <- (-a - 1) * log(x) + (-b/x)  #likelihood  
    y <- c + l
    return(y)
}

dirichletpdflog <- function(x, k) {
    c <- lgamma(sum(k, na.rm = TRUE)) - sum(lgamma(k), 
        na.rm = TRUE)  #normalizing constant
    l <- sum((k - 1) * log(x), na.rm = TRUE)  #likelihood
    y <- c + l
    return(y)
}

estimateDirichletParamsMap <- function(kappa) {
    K <- length(kappa)
    pi <- (kappa - 1)/(sum(kappa, na.rm = TRUE) - K)
    return(pi)
}

estimateInvGammaParamsMap <- function(alpha, beta) {
    mu <- beta/(alpha + 1)
    return(mu)
}

estimateBetaParamsMap <- function(alpha, beta) {
    mu <- (alpha - 1)/(alpha + beta - 2)
    return(mu)
} 
