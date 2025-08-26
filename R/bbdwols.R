#' Bayesian Bootstrap dWOLS
#' @param outcome.mod A formula for the outcome model
#' @param trt.mod A Formula for the treatment model
#' @param dat the data.frame with all data, must have variables for outcome and treatment model
#' @param trt.name A string giving the name of the outcome in dat
#' @param L The number of bootstrap samples to take
#' 
#' @return A data.frame with columns corresponding to parameters in the outcome model and each row is an iteration of the bbdwols

bbdwols <- function(outcome.mod,
                    trt.mod,
                    dat, L,
                    trt.name, maxit = 500) {
  
  n <- nrow(dat)
  
  for(i in 1:L) {
    
    dat$omega <- t(rdirichlet(1, rep(1,n)))
    
    # Get treatment model weights 
    tm <- multinom(trt.mod,
                   data = dat,
                   maxit = maxit,
                   weights = omega,
                   trace = F)
    if(tm$convergence !=0) {
      
      warning(paste0("Iteration l = ", i, " nnet did not converge"))
      
    }
    
    pis <- predict(tm, type = "probs") # conditional trt probs for all trts
    ix <- match(dat[,trt.name], colnames(pis)) # Index of treatment they actually received
    
    dat$iptws <- 1/pis[cbind(1:n, ix)] 
    
    # Estimate blip parameters
    bmod <- lm(outcome.mod,
               weights = omega*iptws,
               data = dat)
    
    if(i == 1) { # set up matrix after first iteration so we know how many parameters there are
      
      bpars <- matrix(nrow = L, ncol = length(coef(bmod)))
      colnames(bpars) <- names(coef(bmod))
      
    }
    
    bpars[i,] <- coef(bmod)
    
  }
  
  return(as.data.frame(bpars))
  
}


# Assess convergence

mc_check <- function(obj, ...) {
  
  # Remove parameters that are unestimable so all NA
  notna <- apply(obj, 2, function(x) !anyNA(x))
  
  obj <- obj[,notna]
  
  L <- nrow(obj)
  
  npar <- ncol(obj)
  
  inc <- ceiling(L/1000)
  
  # set up data.frame
  calc <- seq(2, L, by = inc)
  mc_df <- matrix(nrow = length(calc), ncol = npar+1)
  
  for(i in 1:length(calc)) {
    
    mc_df[i,1] <- calc[i]
    
    mc_df[i,2:ncol(mc_df)] <- apply(obj[1:calc[i],], 2, sd)
    
  }
  
  # Make plots
  xlab <- "MC Iteration"
  ylabs <- paste0(colnames(obj), "\nMC SD")
  par(mfrow = c(4, ceiling(npar/4)), mar = c(4, 6, 2, 1), ...)
  
  for(k in 1:npar) {
    
    p <- plot(mc_df[,1], mc_df[,k+1], type = "l", col = "orchid4", xlab = xlab, ylab = ylabs[k])
    invisible(p)
    
  }
  
  # return(mc_df)
  
  
}