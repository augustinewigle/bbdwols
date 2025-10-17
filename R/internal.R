# Internal functions

# Function that executes bbdwols

i.bbdwols <- function(outcome.mod,
                      trt.mod,
                      dat, L,
                      trt.name, maxit = 500,
                      missingoutcome = FALSE,
                      outmiss.mod = NULL
                      ) {

  dat <- as.data.frame(dat)

  outcome.name <- paste(f_lhs(outcome.mod))

  # check for missing values in the outcome
  dat$outmiss <- is.na(dat[,outcome.name])*1

  if(sum(dat$outmiss) == 0) { # no missing outcomes

    if(missingoutcome) {

      warning("No missing outcomes detected")
      missingoutcome <- FALSE

    }

  } else { # missing outcomes

    if (!missingoutcome) {

      warning("Observations with missing outcomes will be removed from the analysis. To adjust for missing outcomes through weighting, set `missingoutcome = TRUE`")

      dat <- dat[(dat$outmiss==0),]

    } else if(is.null(outmiss.mod)) {

      stop("The outcome missingness model must be specified through `outmiss.mod` if `missingoutcome = TRUE`.")

    }

  }

  n <- nrow(dat)

  # Start the Bayesian Bootstrap
  for(i in 1:L) {

    dat$omega <- t(rdirichlet(1, rep(1,n)))
    dat$omegan <- dat$omega*n # improves convergence rate

    if(missingoutcome) {

      # estimate outcome missingness propensity

      outmiss.mod <- reformulate(paste(rlang::f_rhs(outmiss.mod)), "outmiss")

      outmissm <- glm(outmiss.mod, family = quasibinomial(),
                       data = dat,
                       weights = omegan)

      dat$prob.miss <- outmissm$fitted.values
      # dat$omegamiss <- dat$omega/dat$prob.miss

      # Get treatment model weights

      # scale numeric data columns so they are between 0 and 1 for faster nnet convergence
      tmod.dat <- dat[,colnames(dat)!=outcome.name]
      ind <- sapply(tmod.dat, is.numeric)
      tmod.dat[ind] <- lapply(tmod.dat[ind], function(x) (x-min(x))/(max(x)-min(x)))

      # Add missing probability with dirichlet weights
      tmod.dat$omegamiss <- dat$omega/dat$prob.miss

      tm <- nnet::multinom(trt.mod,
                           data = tmod.dat,
                           maxit = maxit,
                           weights = omegamiss, # weight for prob missing
                           trace = F, model = TRUE)


    } else {

      # scale numeric data columns so they are between 0 and 1 for faster convergence in nnet
      tmod.dat <- dat[,colnames(dat)!=outcome.name]
      ind <- sapply(tmod.dat, is.numeric)
      tmod.dat[ind] <- lapply(tmod.dat[ind], function(x) (x-min(x))/(max(x)-min(x)))

      # Add missing probability with dirichlet weights
      tmod.dat$omegan <- dat$omegan

      # Get treatment model weights

      tm <- nnet::multinom(trt.mod,
                           data = tmod.dat,
                           maxit = maxit,
                           weights = omegan,
                           trace = F, model = TRUE)

    }


    if(tm$convergence !=0) {

      warning(paste0("Iteration l = ", i, " nnet did not converge"))

    }

    pis <- predict(tm, type = "probs") # conditional trt probs for all trts

    if(length(unique(dat[,trt.name])) == 2) { # binary treatment, need to process probabilities into matrix

      pis <- cbind(pis, 1-pis)
      colnames(pis) <- c(tm$lev[1], tm$lev[2])

    }

    ix <- match(dat[,trt.name], colnames(pis)) # Index of treatment they actually received

    dat$iptws <- 1/pis[cbind(1:nrow(dat), ix)]

    # Estimate blip parameters
    bmod <- lm(outcome.mod,
               weights = omega*iptws,
               data = dat,
               subset = outmiss==0)

    if(i == 1) { # set up matrix after first iteration so we know how many parameters there are

      bpars <- matrix(nrow = L, ncol = length(coef(bmod)))
      colnames(bpars) <- names(coef(bmod))

    }

    bpars[i,] <- coef(bmod)

  }

  return(as.data.frame(bpars))

}
