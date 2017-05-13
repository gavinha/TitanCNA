# author: Gavin Ha 
# 		  Dana-Farber Cancer Institute
#		  Broad Institute
# contact: <gavinha@gmail.com> or <gavinha@broadinstitute.org>
# date:	  March 17, 2017

#### EM (FWD-BACK) Algorithm ####
runEMclonalCN <- function(data, params,
                          #gParams = NULL, nParams = NULL, pParams = NULL, sParams = NULL, 
                          txnExpLen = 1e+15, txnZstrength = 5e+5, 
                          maxiter = 15, maxiterUpdate = 1500, pseudoCounts = 1e-300, 
                          normalEstimateMethod = "map", estimateS = TRUE, estimatePloidy = TRUE, 
                          useOutlierState = FALSE, likChangeThreshold = 0.001, verbose = TRUE) {
    ## check that arguments contain necessary list elements
    if (is.null(params$genotypeParams) || is.null(params$ploidyParams) ||
        is.null(params$normalParams) || is.null(params$cellPrevParams)){
      stop("params must include genotypeParams, ploidyParams, normalParams, cellPrevParams.")
    }else{
      gParams <- params$genotypeParams
      pParams <- params$ploidyParams
      sParams <- params$cellPrevParams
      nParams <- params$normalParams
    }
  
    if (is.null(data$chr) || is.null(data$posn) || 
        is.null(data$ref) || is.null(data$tumDepth) ||
        is.null(data$logR)) {
        stop("data must contain named list elements: chr, posn, ref, tumDepth, 
             logR. See loadDataFromFile()")
    }
    if (is.null(gParams$rt) || is.null(gParams$ct) || 
        is.null(gParams$rn)  || is.null(gParams$ZS) || 
        is.null(gParams$var_0) || is.null(gParams$alphaKHyper) || 
        is.null(gParams$betaKHyper) || is.null(gParams$kappaGHyper)) {
        stop("genotypeParams must contain named list elements: rt, rn, ct, ZS, 
             var_0, alphaKHyper, betaKHyper, kappaGHyper.\n         
             See loadDefaultParameters().")
    }
    if (is.null(pParams$phi_0) || is.null(pParams$alphaPHyper) ||
        is.null(pParams$betaPHyper)) {
        stop("ploidyParams must contain named list elements: phi_0, alphaPHyper, 
             betaPHyper. See loadDefaultParameters()")
    }
    if (is.null(nParams$n_0) || is.null(nParams$alphaNHyper) || 
        is.null(nParams$betaNHyper)) {
        stop("normalParams must contain named list elements: n_0, alphaNHyper, 
             betaNHyper. See loadDefaultParameters()")
    }
    if (is.null(sParams$s_0) || is.null(sParams$kappaZHyper) || 
        is.null(sParams$alphaSHyper) || is.null(sParams$betaSHyper)) {
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
        gNoOUTStateParams <- excludeGarbageState(gParams, K)
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
    var <- matrix(0, K, maxiter)  # variance parameter (logR) + outlier state
    varR <- matrix(0, K, maxiter) # variance parameter (allelic) + outlier state
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
    chrPos_0 <- unlist(lapply(chrsI, head, 1)) # index of first datapoint for each chromosome
    ## INITIALIZATION ##
    piG_0 <- estimateDirichletParamsMap(gParams$kappaGHyper)  #add the outlier state
    piZ_0 <- estimateDirichletParamsMap(sParams$kappaZHyper)
    i <- 1
    loglik[i] <- -Inf
    converged <- 0  # flag for convergence
    s[, i] <- sParams$s_0
    var[, i] <- gParams$var_0
    varR[, i] <- gParams$varR_0
    # var[1,i] <- gParams$outlierVar
    phi[i] <- pParams$phi_0
    n[i] <- nParams$n_0
    piG[, i] <- gParams$piG_0
    if (useOutlierState){ #initialize outlier state 
      piG[1, i] <- gParams$piG_0[1] * gParams$piZ_0[1]
    }  
    piZ[, i] <- sParams$piZ_0
    musTmp <- clonalTwoComponentMixtureCN(gNoOUTStateParams$rt, 
        gParams$rn, nParams$n_0, sParams$s_0, gNoOUTStateParams$ct, 
        pParams$phi_0)
    mus$R[, , i] <- musTmp$R
    mus$C[, , i] <- musTmp$C
    
    # compute likelihood conditional on k and z
    # (K-by-Z-by-N)
    if (gParams$alleleEmissionModel == "Gaussian"){
    #  pyR <- computeNormalObslik(data$ref / data$tumDepth, 
    #                             matrix(mus$R[, , i], KnoOutlier, Z),
    #                             gNoOUTStateParams$varR_0)
      py <- computeBivariateNormalObslik(data$ref / data$tumDepth, data$logR, 
                                        matrix(mus$R[, , i], KnoOutlier, Z), 
                                        matrix(mus$C[, , i], KnoOutlier, Z),
                                        gNoOUTStateParams$varR_0, gNoOUTStateParams$var_0,
                                        gNoOUTStateParams$corRho_0)
      py <- py + pseudoCounts
    }else{ # if (gParams$alleleEmissionModel == "binomial"){
      pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
                                   matrix(mus$R[, , i], KnoOutlier, Z))
      pyC <- computeNormalObslik(data$logR, matrix(mus$C[, , i], KnoOutlier, Z), 
                               gNoOUTStateParams$var_0)
      if (useOutlierState == 1) {
        pyO <- outlierObslik(data$ref, data$tumDepth, 
            data$logR, gParams$outlierVar)
        py <- exp(rbind(pyO$R, pyR) 
        			+ rbind(pyO$C, pyC)) + pseudoCounts  #add the outlier state
      } else {
          #joint likelihood between Binomial and Gaussian
          py <- pyR + pyC + pseudoCounts
          rm(pyR, pyC)
      }
    }
    ## EXPECTATION MAXIMIZATION
    while (!converged && (i < (maxiter))) {
        # clear the previous iteration of rho and garbage
        # collect
        rm(fwdBackPar, rhoZ, rhoG, musTmp)
        gc(verbose = FALSE, reset = TRUE)
        ticId <- proc.time()
        i <- i + 1
        ## E-STEP: FORWARDS-BACKWARDS ALGORITHM;
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
            #inner product then flatten to 1-by-K*Z
            piGiZi[c, KtotalRange] <- as.vector(piGi[[c]] %*% t(piZi[[c]]))  
            if (useOutlierState) { #add outlier state
                  piGiZi[c, ] <- c(piG[1, i - 1], piGiZi[c, KtotalRange])
            }  
        }
        gc(verbose = FALSE, reset = TRUE)
        ## PARALLELIZATION
        fwdBackPar <- foreach(c = 1:numChrs, .combine = rbind, .noexport = c("data")) %dopar% {
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
            loglik[i] <- sum(do.call(rbind, fwdBackPar[, 2]))  #combine loglik
            # marginal probs or responsibilities,
            # p(G_t,Z_t|Data,Params)
            rho <- do.call(cbind, fwdBackPar[, 1])  #combine rho from parallel runs  
            #Zcounts <- do.call(t(colSums(aperm(output$xi, c(3, 2, 1)))))
        } else {
            loglik[i] <- sum(do.call(rbind, fwdBackPar[2]))
            rho <- do.call(cbind, fwdBackPar[1])
            #Zcounts <- do.call(cbind, )
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
            data$logR, rho, rhoG, rhoZ, n[i - 1], s[, i - 1], phi[i - 1], 
            var[kRange, i - 1], varR[kRange, i - 1], gNoOUTStateParams$corRho_0, 
            piG[, i - 1, drop = FALSE], piZ[, i - 1, drop = FALSE], chrPos_0,
            gNoOUTStateParams, nParams, sParams, pParams, 
            maxiter = maxiterUpdate, 
            normalEstimateMethod = normalEstimateMethod, 
            estimateS = estimateS, estimatePloidy = estimatePloidy, 
            verbose = verbose)
        n[i] <- estimateOut$n
        s[, i] <- estimateOut$s
        var[kRange, i] <- estimateOut$var
        varR[kRange, i] <- estimateOut$varR
        phi[i] <- estimateOut$phi
        piZ[, i] <- estimateOut$piZ #estimateClonalMixWeightsParamMap(rhoZ, sParams$kappaZHyper)
        piG[, i] <- estimateOut$piG #estimateGenotypeMixWeightsParamMap(rhoG, gParams$kappaGHyper)
        if (useOutlierState) { #garbage state
                var[1, i] <- gParams$outlierVar
        }   
        rm(rho, estimateOut, outRhoRow)  #clear rho and estimateOut
        ## Recompute the likelihood conditional on k and z
        musTmp <- clonalTwoComponentMixtureCN(gNoOUTStateParams$rt, 
            gParams$rn, n[i], s[, i], gNoOUTStateParams$ct, phi[i])
        mus$R[, , i] <- musTmp$R
        mus$C[, , i] <- musTmp$C
        if (gParams$alleleEmissionModel == "Gaussian"){
        #  pyR <- computeNormalObslik(data$ref / data$tumDepth, 
        #                             matrix(mus$R[, , i], KnoOutlier, Z),
        #                             varR[kRange, i])
          py <- computeBivariateNormalObslik(data$ref / data$tumDepth, data$logR, 
                                             matrix(mus$R[, , i], KnoOutlier, Z), 
                                             matrix(mus$C[, , i], KnoOutlier, Z),
                                             varR[kRange, i], var[kRange, i],
                                             gNoOUTStateParams$corRho_0)
          py <- py + pseudoCounts
        }else{ # if (gParams$alleleEmissionModel == "binomial"){
          pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
                                       matrix(mus$R[, , i], KnoOutlier, Z))
        
          pyC <- computeNormalObslik(data$logR, 
                                   matrix(mus$C[, , i], KnoOutlier, Z), var[kRange, i])
            if (useOutlierState == 1) {
                pyO <- outlierObslik(data$ref, data$tumDepth, 
                    data$logR, gParams$outlierVar)
                py <- exp(rbind(pyO$R, pyR) 
                		+ rbind(pyO$C, pyC)) + pseudoCounts  #add the outlier state
            } else {
                #py <- exp(pyR + pyC) + pseudoCounts  #joint likelihood between Binomial and Gaussian
                py <- pyR + pyC + pseudoCounts
                rm(pyR, pyC)
            }
        }
        prior <- priorProbs(n[i], s[, i], phi[i], var[, i], varR[, i], piG[, i], piZ[, i], gParams, 
                            normalEstimateMethod, estimateS, estimatePloidy,
                            nParams$alphaNHyper, nParams$betaNHyper, sParams$alphaSHyper, 
                            sParams$betaSHyper, pParams$alphaPHyper, pParams$betaPHyper,
                            gNoOUTStateParams$alphaKHyper, gNoOUTStateParams$betaKHyper, 
                            gNoOUTStateParams$alphaRHyper, gNoOUTStateParams$betaRHyper, 
                            gNoOUTStateParams$kappaGHyper, sParams$kappaZHyper)
        ## Compute log-likelihood and check converge
        
        if (verbose == TRUE) {
            message("fwdBack: loglik=", sprintf("%0.4f", loglik[i]))
            message("fwdBack: priorN=", sprintf("%0.4f", prior$n))
            message("fwdBack: priorS=", sprintf("%0.4f", prior$s))
            message("fwdBack: priorVar=", sprintf("%0.4f", prior$var))
            message("fwdBack: priorVarR=", sprintf("%0.4f", prior$varR))
            message("fwdBack: priorPhi=", sprintf("%0.4f", prior$phi))
            message("fwdBack: priorPiG=", sprintf("%0.4f", prior$piG))
            message("fwdBack: priorPiZ=", sprintf("%0.4f", prior$piZ))
        }
        loglik[i] <- loglik[i] + prior$prior
        if (verbose == TRUE){
          message("fwdBack: EM iteration ", i - 1, 
                " complete loglik=", sprintf("%0.4f", loglik[i]))
        }
        if (((abs(loglik[i] - loglik[i - 1])/abs(loglik[i])) <= likChangeThreshold) 
            && (loglik[i] >= loglik[i - 1])) {
            converged <- 1
        } else if (loglik[i] < loglik[i - 1]) {
            # stop('Failed EM!')
            converged <- 1
            i <- i - 1
            message("fwdBack: Optimization during update decreased complete 
                    likelihood.  Stopping EM...")
        }
        
        elapsedTime <- (proc.time() - ticId)/60
        if (verbose == TRUE) {
            message("fwdBack: Elapsed time for iteration ", i - 1, 
                    ": ", sprintf("%0.4f", elapsedTime[3]), "m")
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
    output$varR <- varR[, 1:i, drop = FALSE]
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
    output$cellPrevParams <- sParams
    output$symmetric <- gParams$symmetric
    return(output)
}

viterbiClonalCN <- function(data, convergeParams, genotypeParams = NULL) {
    ## requirements for parallelization require(foreach)
    ## use genotypeParams found in convergeParams unless
    ## genotypeParams given
    if (is.null(genotypeParams)) {
        if (length(convergeParams$genotypeParams) > 0) {
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
    if (is.null(genotypeParams$alleleEmissionModel)){
      genotypeParams$alleleEmissionModel <- "binomial" #default
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
    pseudoCounts <- 1e-300
    if (genotypeParams$alleleEmissionModel == "Gaussian"){
      #pyR <- computeNormalObslik(data$ref / data$tumDepth, 
      #                             matrix(convergeParams$muR[, , i], KnoOutlier, Z),
      #                             convergeParams$varR[kRange, i])
      py <- computeBivariateNormalObslik(data$ref / data$tumDepth, data$logR, 
                                  matrix(convergeParams$muR[, , i], KnoOutlier, Z), 
                                  matrix(convergeParams$muC[, , i], KnoOutlier, Z),
                                  convergeParams$varR[kRange, i], convergeParams$var[kRange, i],
                                  convergeParams$corRho_0)
    }else{ # if (gParams$alleleEmissionModel == "binomial"){
      pyR <- computeBinomialObslik(data$ref, data$tumDepth, 
                                  matrix(convergeParams$muR[, , i], KnoOutlier, Z))
      pyC <- computeNormalObslik(data$logR, matrix(convergeParams$muC[, , i], KnoOutlier, Z), 
                               convergeParams$var[kRange, i])
      if (useOutlierState) {
        pyO <- outlierObslik(data$ref, data$tumDepth, 
            data$logR, genotypeParams$outlierVar)
        py <- rbind(pyO$R, pyR) 
        		+ rbind(pyO$C, pyC) + pseudoCounts  #add the outlier state
      } else {
          py <- pyR + pyC + pseudoCounts  #joint likelihood between Binomial and Gaussian
          rm(pyR, pyC)
      }
    }
    
    piGiZi <- matrix(0, numChrs, Ktotal)
    for (c in 1:numChrs) {
        piZi[[c]] <- convergeParams$piZ[, i - 1, drop = FALSE]
        piGi[[c]] <- convergeParams$piG[kRange, i - 1, drop = FALSE]
        piGiZi[c, KtotalRange] <- as.vector(piGi[[c]] %*% 
            t(piZi[[c]]))  #inner product then flatten to 1-by-K*Z
        if (useOutlierState) 
            {
                piGiZi[c, ] <- c(convergeParams$piG[1, i - 1], piGiZi[c, KtotalRange])
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
            py[k + (z - 1) * K, ] <- normalpdflog(x, mus[k, z], var[k])
        }
    }
    # py <- t(matrix(py,N,K*Z))
    return(py)
}

computeBivariateNormalObslik <- function(x1, x2, mu1, mu2, var1 = NULL, var2 = NULL, corRho = NULL){
  N <- length(x1)
  K <- nrow(mu1)
  Z <- ncol(mu1)
  coVar <- getBivariateCovariance(x1, x2)
  if (is.null(var1) || is.null(var2)){
    coVar <- getBivariateCovariance(x1, x2)
    var1 <- rep(coVar$covar[1, 1], K)
    var2 <- rep(coVar$covar[2, 2], K)
    corRho <- rep(coVar$cor, K)
  }
  if (is.null(corRho)){
    corRho <- coVar$cor
  }
  py <- matrix(0, K * Z, N)  #local evidence is 2-D matrix: K*ZxN
  for (k in 1:K) {
    for (z in 1:Z) {
      py[k + (z - 1) * K, ] <- bivariateNormalpdflog(x1, x2, mu1[k, z], mu2[k, z], 
                            var1[k], var2[k], corRho)
    }
  }
  return(py)
}

# Bivariate Gaussian density function, return in log scale #
bivariateNormalpdflog <- function(x1, x2, mu1, mu2, var1, var2, corRho){
  c <- log(2 * pi * sqrt(var1 * var2 * (1 - corRho ^ 2)))
  l_0 <- 1 / (2 * (1 - corRho ^ 2))
  l_1 <- ((x1 - mu1) ^ 2) / var1
  l_2 <- ((x2 - mu2) ^ 2) / var2 
  l_3 <-  2 * corRho * (x1 - mu1) * (x2 - mu2) / sqrt(var1 * var2)
  y <- -c - l_0 * (l_1 + l_2 - l_3)
  y[is.na(y)] <- 1
  return(y)
}

getBivariateCovariance <- function(x1, x2){
  # sample covariance 
  covar <- cov(cbind(x1, x2), use = "pairwise.complete.obs")
  # correlation
  cor <- covar[1, 2] / sqrt(covar[1, 1] * covar[2, 2])
  return(list(covar = covar, cor = cor))
}

asCovarianceMatrix <- function(varA, varB, rho){
  covMat <- matrix(c(varA, rho * sqrt(varA * varB), rho * sqrt(varA * varB), varB), 
                   ncol = 2, nrow = 2, byrow = TRUE)
  return(covMat)
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

# likelihood function for covariance #
covarLikelihoodFun <- function(covar, mu, params, D, rho){
  #print(covar)
  covar <- asCovarianceMatrix(covar[1], covar[2], params$cor_0)
  lik <- rho %*% bivariateNormalpdflog(D, mu, covar, params$cor_0)
  prior <- invWishartpdflog(covar, params$psi, params$nu)
  f <- lik + prior
  return(f)
}



# Multivariate gamma function 
multiGammaFunlog <- function(p, a){
  if (p == 1){
    return(lgamma(a))
  }else{
    return(log(pi ^ ((p - 1) / 2) + multiGammaFunlog(p - 1, a) + multiGammaFunlog(1, a + (1 - p) / 2)))
  }
}

# Inverse Wishart density function in log scale 
invWishartpdflog <- function(covar, psi, nu){
  p <- ncol(covar)
  const <- log(det(psi)) * (nu / 2) - log(2) * (nu * p / 2) + multiGammaFunlog(p, nu / 2)
  lik <- log(det(covar)) * -((nu + p + 1) / 2) - 0.5 * sum(diag(psi %*% solve(covar)))
  y <- const + lik
  return(y)
}



priorProbs <- function(n, s, phi, var, varR, piG, piZ, gParams,
                       normalEstimateMethod, estimateS, estimatePloidy,
                       alphaN, betaN, alphaS, betaS, alphaP, betaP,
                       alphaK, betaK, alphaR, betaR, kappaG, kappaZ){
  Z <- length(s)
  K <- length(var)
  ## prior for pi's ##
  priorPiG <- dirichletpdflog(piG, kappaG)
  priorPiZ <- dirichletpdflog(piZ, kappaZ)
  ## prior for S ##
  prior_beta_s <- 0
  if (estimateS) {
    for (z in 1:Z) {
      #prior_beta_s <- prior_beta_s + log(1/beta(alphaS[z], 
      #    betaS[z])) + (alphaS[z] - 1) * log(s[z]) + 
      #    (betaS[z] - 1) * log(1 - s[z])
      prior_beta_s <- prior_beta_s + 
        betapdflog(s[z], alphaS[z], betaS[z])
    }
  } else {
    prior_beta_s <- 0
  }
  ## prior for n ##
  if (normalEstimateMethod == "map") {
    #prior_beta_n <- log(1/beta(alphaN, betaN)) + 
    (alphaN - 1) * log(n) + (betaN - 1) * log(1 - n)
    prior_beta_n <- betapdflog(n, alphaN, betaN)
  } else {
    prior_beta_n <- 0
  }
  ## prior for variance var ##
  cnLevel <- unique(gParams$ct)
  cnInd <- sapply(cnLevel, function(x) {
    which(gParams$ct == x)[1]
  })
  prior_gamma_var <- 0
  for (k in cnInd) {
    #prior_gamma_var <- prior_gamma_var + 
    #  alphaK[k] * log(betaK[k]) - lgamma(alphaK[k]) + # constant term
    #  (-alphaK[k] - 1) * log(var[k]) - betaK[k]/var[k] # beta likelihood
    prior_gamma_var <- prior_gamma_var + invgammapdflog(var[k], alphaK[k], betaK[k])
  }
  
  ## prior for phi ##
  if (estimatePloidy) {
    #prior_gamma_phi <- alphaP * log(betaP) - lgamma(alphaP) + 
    #    (-alphaP - 1) * log(phi) - betaP/phi
    prior_gamma_phi <- invgammapdflog(phi, alphaP, betaP)
  } else {
    prior_gamma_phi <- 0
  }
  
  prior_gamma_varR <- 0
  if (gParams$alleleEmissionModel == "Gaussian"){
    for (k in 1:K){
      #prior_gamma_varR <- prior_gamma_varR + 
      #  alphaR[k] * log(betaR[k]) - lgamma(alphaR[k]) + #constant term
      #  (-alphaR[k] - 1) * log(varR[k]) - betaR[k]/varR[k]
      prior_gamma_varR <- prior_gamma_varR + invgammapdflog(varR[k], alphaR[k], betaR[k])
    }
  }
  
  prior <- prior_beta_s + prior_beta_n + prior_gamma_phi + 
    prior_gamma_var + prior_gamma_varR + priorPiG + priorPiZ
  return(list(prior = prior, s = prior_beta_s, n = prior_beta_n, 
              phi = prior_gamma_phi,
              var = prior_gamma_var, varR = prior_gamma_varR, 
              piG = priorPiG, piZ = priorPiZ))
}

