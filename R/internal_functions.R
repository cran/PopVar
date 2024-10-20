#' Internal functions
#' 
#' @name internal
#' 
#' @description 
#' Internal functions generally not to be called by the user.
#' 
#' @param crossing.table The crossing table.
#' @param par.entries The parent entries.
#' @param crossing.mat The crossing matrix.
#' @param GEBVs The genomic estimated breeding values.
#' @param tail.p The proportion from the population to select.
#' @param G The marker genotypes
#' @param y.CV The phenotypes for cross-validation.
#' @param G.CV The marker genotypes for cross-validation.
#' @param models.CV The models for cross-validation.
#' @param frac.train.CV The fraction of data to use as training data in cross-validation.
#' @param nCV.iter.CV The number of iterations of cross-validation.
#' @param burnIn.CV The burn-in number for cross-validation.
#' @param nIter.CV The number of iterations for Bayesian models in cross-validation.
#' @param nFold.CV The number of folds in k-fold cross-validation.
#' @param nFold.CV.reps The number of replications of k-fold cross-validation.
#' @param M The marker matrix.
#' @param y.df The phenotype data.
#' @param models The models to use.
#' @param nIter The number of iterations.
#' @param burnIn The burn-in rate.
#' 
#' @importFrom stats cor sd
#' @importFrom BGLR BGLR
#' 

#' 
#' @rdname internal
#' 
par_position <- function(crossing.table, par.entries) { # Used when a crossing table is defined
  
  par.pos <- matrix(nrow = nrow(crossing.table), ncol = 2)
  crosses.possible <- matrix(nrow = nrow(crossing.table), ncol = 2)
  for(i in 1:nrow(crossing.table)){
    par1 <- as.character(crossing.table[i,1])
    par2 <- as.character(crossing.table[i,2])
    
    if(par1 %in% par.entries & par2 %in% par.entries){
      par.pos[i,1] <- which(par.entries == par1)
      par.pos[i,2] <- which(par.entries == par2)
      crosses.possible[i,1] <- par1
      crosses.possible[i,2] <- par2  
    }
  }
  
  par.pos <- par.pos[which(!is.na(par.pos[,1])), ]
  crosses.possible <- crosses.possible[which(!is.na(crosses.possible[,1])), ]
  
  for.dup <- paste(par.pos[,1], ".", par.pos[,2], sep=""); duplicated <- which(duplicated(for.dup))
  if(length(duplicated) > 0){
    par.pos <- par.pos[-duplicated, ]
    crosses.possible <- crosses.possible[-duplicated, ]
  }
  
  return(list(parent.positions=par.pos, crosses.possible=crosses.possible))
}

#'
#' @rdname internal
#' 
par_name <- function(crossing.mat, par.entries){ ## Used when all combinations of parents are crossed
  crosses.possible <- matrix(nrow = nrow(crossing.mat), ncol = 2)
  for(i in 1:nrow(crossing.mat)){
    crosses.possible[i,1] <- par.entries[crossing.mat[i,1]]
    crosses.possible[i,2] <- par.entries[crossing.mat[i,2]]
  }
  return(crosses.possible)
}

#'
#' @rdname internal
#' 
tails <- function(GEBVs, tail.p){ #Calculates means of tails; set tail.p to the proportion of the populaiton you want to take the mean of, default is 10%
  u.top <- mean(GEBVs[which(GEBVs >= quantile(GEBVs, 1-tail.p))], na.rm=T)
  u.bot <- mean(GEBVs[which(GEBVs <= quantile(GEBVs, tail.p))], na.rm=T)
  
  return(rbind(u.top, u.bot))
}

#'
#' @rdname internal
#' 
maf_filt <- function(G){
  G_noNA <- sum(!is.na(G), na.rm = TRUE)
  freq1 <- sum(G == 1, na.rm = TRUE) / G_noNA #+ .5*(sum(G == 0, na.rm = TRUE) / G_noNA)
  min(freq1, 1 - freq1)
}


#test.map <- qtl::sim.map(len = c(50,50), n.mar = c(100,100), anchor.tel = FALSE, include.x = FALSE, sex.sp = TRUE)
#View(qtl::map2table(test.map))


### FldTrial.lm -- this fucntion is actually better, replace the other one with this
#Xlm <- function(y, X){
#  n <- nrow(X)
#  k <- ncol(X)
#  
#  inv.tX.X <- solve(crossprod(X))
#  #if(!exists("inv.tX.X")) inv.tX.X <- MASS::ginv(t(X) %*% X) ## If chol fails then use general inverse of MASS package
#  
#  beta <- matrix(inv.tX.X %*% t(X) %*% y, ncol = 1, dimnames = list(colnames(X), NULL))
#  
#  yhat <- X %*% beta
#  e <- as.vector(y - yhat)
#  
#  sigma.samp <- (t(e) %*% e) / (n-k)
#  varB <- sigma.samp %x% inv.tX.X
#  se.B <- sqrt(sigma.samp) %x% sqrt(diag(inv.tX.X))
#  rownames(varB) <- colnames(varB) <- rownames(beta)
#  
#  return(list(betahat=beta, resids=e, varbetahat= diag(varB), stderrbetahat=se.B, dfErr=(n-k)))
#}




## G.CV must have individuals as rows and markers as columns... so 100 individuals with 500 markers would be a 100x500 matrix
## y.CV must be a numeric vector representing a continuous variable


#'
#' @rdname internal
#' 
XValidate_nonInd <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, frac.train.CV=NULL, nCV.iter.CV=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
  gc(verbose = F) ## Close unused connections
  #con.path <- getwd() ## BGLR will write temp files to the wd
  
  non.BGLR <- models.CV[models.CV %in% c("rrBLUP")]
  BGLR <- models.CV[models.CV %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")]
  
  for(i in 1:nCV.iter.CV){
    if(i==1){
      rrBLUP.cv <- c()
      BayesA.cv <- c()
      BayesB.cv <- c()
      BayesC.cv <- c()
      BL.cv <- c()
      BRR.cv <- c()
    }
    
    TP.sample <- sample(x=1:length(y.CV), size=round(frac.train.CV*length(y.CV)), replace=F)
    VP.sample <- setdiff(1:length(y.CV), TP.sample)
    
    TP.G <- G.CV[TP.sample, ]
    TP.y <- y.CV[TP.sample]
    
    VP.G <- G.CV[VP.sample, ]
    VP.y <- y.CV[VP.sample]
    
    ### non-BGLR models below
    if("rrBLUP" %in% non.BGLR) {RR.pred <- rrBLUP::kinship.BLUP(y=TP.y, G.train=TP.G, G.pred=VP.G, K.method="RR"); rrBLUP.cv[i] <- cor(VP.y, RR.pred$g.pred, use="pairwise.complete.obs")}
    
    ### models from BGLR package
    ## Bayesian A
    if("BayesA" %in% BGLR){
      tryCatch({
        BayesA.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesA")), verbose=F, nIter=1500, burnIn=500)
      }, warning=function(w){
        gc(verbose = F)
        BayesA.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BayesA.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesA.fit)){
        mkr.effs <- as.numeric(BayesA.fit$ETA[[1]]$b); BayesA.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesA.cv[i] <- NA}
    }
    
    ## Bayesian B 
    if("BayesB" %in% BGLR){
      tryCatch({
        BayesB.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesB")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BayesB.fit <- NULL
      }, error=function(e){
        gc(verbose = )
        BayesB.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesB.fit)){
        mkr.effs <- as.numeric(BayesB.fit$ETA[[1]]$b)
        BayesB.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesB.cv[i] <- NA}
    }
    
    ## Bayesian C
    if("BayesC" %in% BGLR){
      tryCatch({
        BayesC.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesC")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BayesC.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BayesC.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BayesC.fit)){
        mkr.effs <- as.numeric(BayesC.fit$ETA[[1]]$b)
        BayesC.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BayesC.cv[i] <- NA}
    }
    
    ## Bayesian LASSO
    if("BL" %in% BGLR){
      tryCatch({
        BL.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BL")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BL.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BL.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BL.fit)){
        mkr.effs <- as.numeric(BL.fit$ETA[[1]]$b)
        BL.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BL.cv[i] <- NA}
    }
    
    ### Bayesian Ridge Reg.
    if("BRR" %in% BGLR){
      tryCatch({
        BRR.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BRR")), verbose=F, nIter=1500, burnIn=1000)
      }, warning=function(w){
        gc(verbose = F)
        BRR.fit <- NULL
      }, error=function(e){
        gc(verbose = F)
        BRR.fit <- NULL
      }
      ); gc(verbose = F)
      
      if(!is.null(BRR.fit)){
        mkr.effs <- as.numeric(BRR.fit$ETA[[1]]$b)
        BRR.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
      } else{BRR.cv[i] <- NA}
    }
    
  }; gc(verbose = F) ## End of ith iteration
  
  cv.list <- paste(models.CV, "cv", sep=".")
  
  for(l in 1:length(cv.list)){ ## Tried using sapply() across cv.list using the "get" function, but did not work within PopVar function
    toget <- cv.list[l]
    if(l == 1) cvs <- get(toget)
    if(l > 1) cvs <- cbind(cvs, get(toget))
  }
  
  if(length(models.CV) == 1){
    bad.models <- F
    if(length(which(is.na(cvs))) > 0.025*nCV.iter.CV) bad.models <- T
  }
  
  if(length(models.CV) > 1){
    bad.models <- apply(cvs, 2, function(X){if(length(which(is.na(X))) > 0.025*nCV.iter.CV){ ## if more than 2.5% of the iterations resulted in an error 
      return(T)
    }else{return(F)}
    })
  }
  
  if(length(models.CV) > 1) {CV.results <- data.frame(Model=models.CV, r_avg=apply(cvs, 2, mean, na.rm=T), r_sd=apply(cvs, 2, sd, na.rm=T)) ; rownames(CV.results) <- NULL}
  if(length(models.CV) == 1) {CV.results <- data.frame(Model=models.CV, r_avg=mean(cvs), r_sd=sd(cvs)) ; rownames(CV.results) <- NULL}
  
  if(length(which(bad.models)) > 0){
    if(length(models.CV) == 1 | (length(which(bad.models)) == length(models.CV))) stop("All model(s) tested was/were removed due to excessive negative values of nu being returned by BGLR::BGLR")
    CV.results <- CV.results[-which(bad.models), ]
    warning(paste("Model(s)", models.CV[which(bad.models)], "was/were removed due to excessive negative values of nu being returned by BGLR::BGLR."))
  }
  
  #CV.lists <- as.data.frame(t(rbind(as.character(models.CV), matrix(c(rrBLUP.cv, BayesA.cv, BayesB.cv, BayesC.cv, BL.cv, BRR.cv), ncol=length(models.CV)))))
  
  return(list(CV.summary = CV.results)) #, iter.CV = CV.lists))
  
} ## End of XValidate function



## G.CV must have individuals as rows and markers as columns... so 100 individuals with 500 markers would be a 100x500 matrix
## y.CV must be a numeric vector representing a continuous variable

#'
#' @rdname internal
#' 
XValidate_Ind <- function(y.CV=NULL, G.CV=NULL, models.CV=NULL, nFold.CV=NULL, nFold.CV.reps=NULL, burnIn.CV=NULL, nIter.CV=NULL){
  
  gc(verbose = F) ## Close unused connections
  #con.path <- getwd() ## BGLR will write temp files to the wd
  
  non.BGLR <- models.CV[models.CV %in% c("rrBLUP")]
  BGLR <- models.CV[models.CV %in% c("BayesA", "BayesB", "BayesC", "BL", "BRR")]
  if(nFold.CV >= length(y.CV)) stop("nFold too large given the TP size")
  
  for(j in 1:nFold.CV.reps){
    for(i in 1:nFold.CV){
      
      if(i==1){
        rrBLUP.cv <- c()
        BayesA.cv <- c()
        BayesB.cv <- c()
        BayesC.cv <- c()
        BL.cv <- c()
        BRR.cv <- c()
        
        for.VP.mat <- 1:length(y.CV)
        while(length(for.VP.mat) %% nFold.CV != 0) for.VP.mat <- c(for.VP.mat, NA)
        VP.mat <- matrix(for.VP.mat[order(sample(1:length(for.VP.mat), replace = F))], nrow = nFold.CV, byrow = T)
      }
      
      VP.sample <- as.numeric(VP.mat[i, ]); VP.sample <- VP.sample[which(!is.na(VP.sample))]
      TP.sample <- setdiff(1:length(y.CV), VP.sample)
      
      TP.G <- G.CV[TP.sample, ]
      TP.y <- y.CV[TP.sample]
      
      VP.G <- G.CV[VP.sample, ]
      VP.y <- y.CV[VP.sample]
      
      ### non-BGLR models below
      if("rrBLUP" %in% non.BGLR) {RR.pred <- rrBLUP::kinship.BLUP(y=TP.y, G.train=TP.G, G.pred=VP.G, K.method="RR"); rrBLUP.cv[i] <- cor(VP.y, RR.pred$g.pred, use="pairwise.complete.obs")}
      
      ### models from BGLR package
      ## Bayesian A
      if("BayesA" %in% BGLR){
        tryCatch({
          BayesA.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesA")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesA.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesA.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesA.fit)){
          mkr.effs <- as.numeric(BayesA.fit$ETA[[1]]$b); BayesA.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesA.cv[i] <- NA}
      }
      
      
      ## Bayesian B 
      if("BayesB" %in% BGLR){
        tryCatch({
          BayesB.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesB")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesB.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesB.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesB.fit)){
          mkr.effs <- as.numeric(BayesB.fit$ETA[[1]]$b)
          BayesB.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesB.cv[i] <- NA}
      }
      
      
      ## Bayesian C
      if("BayesC" %in% BGLR){
        tryCatch({
          BayesC.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BayesC")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BayesC.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BayesC.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BayesC.fit)){
          mkr.effs <- as.numeric(BayesC.fit$ETA[[1]]$b)
          BayesC.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BayesC.cv[i] <- NA}
      }
      
      
      ## Bayesian LASSO
      if("BL" %in% BGLR){
        tryCatch({
          BL.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BL")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BL.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BL.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BL.fit)){
          mkr.effs <- as.numeric(BL.fit$ETA[[1]]$b)
          BL.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BL.cv[i] <- NA}        
      }
      
      
      ### Bayesian Ridge Reg.
      if("BRR" %in% BGLR){
        tryCatch({
          BRR.fit <- BGLR::BGLR(y=TP.y, ETA=list(list(X=TP.G, model="BRR")), verbose=F, nIter=1500, burnIn=1000)
        }, error=function(e){
          gc(verbose=F)
          BRR.fit <- NULL
        }, warning=function(w){
          gc(verbose=F)
          BRR.fit <- NULL
        }); gc(verbose=F)
        
        if(!is.null(BRR.fit)){
          mkr.effs <- as.numeric(BRR.fit$ETA[[1]]$b)
          BRR.cv[i] <- cor((VP.G %*% mkr.effs), VP.y, use="pairwise.complete.obs")
        } else{BRR.cv[i] <- NA}
      }
      
    }; gc(verbose=F) ## End of i-fold iteration
    
    if(j==1) cv.list <- paste(models.CV, "cv", sep=".")
    
    for(l in 1:length(cv.list)){ ## Tried using sapply() across cv.list using the "get" function, but did not work within PopVar function
      toget <- cv.list[l]
      if(l == 1) cvs <- get(toget)
      if(l > 1) cvs <- cbind(cvs, get(toget))
    }
    
    if(j == 1) cvs.all <- cvs
    if(j > 1) cvs.all <- rbind(cvs.all, cvs)
    
  }
  
  if(length(models.CV) == 1){
    bad.models <- F
    if(length(which(is.na(cvs.all))) > 0.025*(nFold.CV*nFold.CV.reps)) bad.models <- T
  }
  
  if(length(models.CV) > 1){
    bad.models <- apply(cvs.all, 2, function(X){if(length(which(is.na(X))) > 0.025*(nFold.CV*nFold.CV.reps)){ ## if more than 2.5% of the iterations resulted in an error 
      return(T)
    }else{return(F)}
    })
  }
  
  if(length(models.CV) > 1) {CV.results <- data.frame(Model=models.CV, r_avg=apply(cvs.all, 2, mean, na.rm=T), r_sd=apply(cvs.all, 2, sd, na.rm=T)) ; rownames(CV.results) <- NULL}
  if(length(models.CV) == 1) {CV.results <- data.frame(Model=models.CV, r_avg=mean(cvs.all), r_sd=sd(cvs.all)) ; rownames(CV.results) <- NULL}
  
  if(length(which(bad.models)) > 0){
    if(length(models.CV) == 1 | (length(which(bad.models)) == length(models.CV))) stop("All model(s) tested was/were removed due to excessive negative values of nu being returned by BGLR::BGLR")
    CV.results <- CV.results[-which(bad.models), ]
    warning(paste("Model(s)", models.CV[which(bad.models)], "was/were removed due to excessive negative values of nu being returned by BGLR::BGLR."))
  }
  
  #CV.lists <- as.data.frame(t(rbind(as.character(models.CV), matrix(c(rrBLUP.cv, BayesA.cv, BayesB.cv, BayesC.cv, BL.cv, BRR.cv), ncol=length(models.CV)))))
  return(list(CV.summary = CV.results))#, iter.CV = CV.lists))
  
} ## End of XValidate function





# Function to calculate marker effects
# 
# Allows other arguments to be passed
# 

#'
#' @rdname internal
#' 
calc_marker_effects <- function(M, y.df, models = c("rrBLUP", "BayesA", "BayesB","BayesC", "BL", "BRR"), nIter, burnIn) {
  
  models <- match.arg(models)
  
  # Error out if model != rrBLUP and nIter & burnIn are missing
  if (models != "rrBLUP") {
    if (missing(nIter) | missing(burnIn)) stop ("You must provide the arguments 'nIter' and 'burnIn' for that model choice.")
    
  }
  
  # Determine the function to use for marker effect estimation
  if (models == "rrBLUP") {
    cme <- function(M, y) {
      fit <- mixed.solve(y = y, Z = M, method = "REML")
      # Return marker effects and the grand mean
      list(effects = as.matrix(fit$u), grand_mean = fit$beta)
    }
  } else {
    cme <- function(M, y) {
      suppressWarnings(fit <- BGLR(y = y, ETA = list(list(X = M, model = models)), verbose = FALSE, nIter = nIter, burnIn = burnIn))
      list(effects = as.matrix(fit$ETA[[1]]$b), grand_mean = fit$mu)
    }
  }
    
  ## Calculate marker effects for each trait
  me_out <- lapply(X = y.df, FUN = cme, M = M)
  
  # Clean up files if models != rrBLUP
  if (models != "rrBLUP") {
    invisible(file.remove(list.files(path = ".", pattern = paste0("mu|varE|", models), full.names = TRUE)))
  }
 
  # Return me_out
  return(me_out)
   
}





