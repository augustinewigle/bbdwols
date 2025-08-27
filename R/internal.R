# Internal functions

# Function that executes bbdwols

i.bbdwols <- function(outcome.mod,
                      trt.mod,
                      dat, L,
                      trt.name, maxit = 500) {

  n <- nrow(dat)
  dat <- as.data.frame(dat)

  for(i in 1:L) {

    dat$omega <- t(rdirichlet(1, rep(1,n)))

    # Get treatment model weights
    tm <- nnet::multinom(trt.mod,
                   data = dat,
                   maxit = maxit,
                   weights = omega,
                   trace = F, model = TRUE)
    if(tm$convergence !=0) {

      warning(paste0("Iteration l = ", i, " nnet did not converge"))

    }

    pis <- predict(tm, type = "probs") # conditional trt probs for all trts

    if(length(unique(dat[,trt.name])) == 2) { # binary treatment, need to process probabilities into matrix

      pis <- cbind(pis, 1-pis)
      colnames(pis) <- c(tm$lev[1], tm$lev[2])

    }

    ix <- match(dat[,trt.name], colnames(pis)) # Index of treatment they actually received

    dat$iptws <- 1/pis[cbind(1:n, ix)]

    # return(list(ix = ix, n = n, dat = dat, pis = pis, tm = tm)) # TODO remove

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
