#' Function that executes 1 iteration of bbdwols
#' @export
one.bbdwols <- function(data,
                      w,
                      outcome.mod,
                      trt.mod,
                      maxit = 500,
                      missingoutcome = FALSE,
                      outmiss.mod = NULL,
                      scale_vars,
                      trim_weights = 1,
                      mode = "parameter") {

  trt.name <- paste(f_lhs(trt.mod))
  outcome.name <- paste(f_lhs(outcome.mod))

  data <- as.data.frame(data)

  data$outmiss <- is.na(data[,outcome.name])

  if(missingoutcome&(!anyNA(data[,outcome.name]))) { # no missing outcomes in data, don't do missing outcome routine

    missingoutcome <- FALSE

  }

  if(missingoutcome) {

    # estimate outcome missingness propensity

    # missind.name <- paste(f_lhs(outmiss.mod))
    outmiss.mod <- reformulate(paste(rlang::f_rhs(outmiss.mod)), "outmiss")

    outmissm <- glm(outmiss.mod, family = quasibinomial(),
                    data = data,
                    weights = w)

    data$prob.miss <- outmissm$fitted.values
    # data$omegamiss <- data$omega/data$prob.miss

    # Get treatment model weights

    tmod.data <- data[,colnames(data)!=outcome.name]

    if(scale_vars) {

      # scale numeric data columns so they are between 0 and 1 for faster nnet convergence
      ind <- sapply(tmod.data, is.numeric)
      tmod.data[ind] <- lapply(tmod.data[ind], function(x) (x-min(x))/(max(x)-min(x)))

    }

    # Add missing probability with dirichlet weights
    tmod.data$omegamiss <- w/data$prob.miss

    tm <- nnet::multinom(trt.mod,
                         data = tmod.data,
                         maxit = maxit,
                         weights = omegamiss, # weight for prob missing
                         trace = F, Hess = FALSE,
                         model = TRUE)

  } else {

    tmod.data <- data[,colnames(data)!=outcome.name]

    if(scale_vars) {

      # scale numeric data columns so they are between 0 and 1 for faster convergence in nnet

      ind <- sapply(tmod.data, is.numeric)
      tmod.data[ind] <- lapply(tmod.data[ind], function(x) (x-min(x))/(max(x)-min(x)))

    }

    # Add missing probability with dirichlet weights
    tmod.data$w <- w

    # Get treatment model weights

    tm <- nnet::multinom(trt.mod,
                         data = tmod.data,
                         maxit = maxit,
                         weights = w,
                         trace = F, model = TRUE)

  }


  if(tm$convergence !=0) {

    warning(paste0("Treatment model did not converge"))

  }

  pis <- predict(tm, type = "probs") # conditional trt probs for all trts

  if(length(unique(data[,trt.name])) == 2) { # binary treatment, need to process probabilities into matrix

    pis <- cbind(pis, 1-pis)
    colnames(pis) <- c(tm$lev[1], tm$lev[2])

  }

  ix <- match(data[,trt.name], colnames(pis)) # Index of treatment they actually received

  iptws <- 1/pis[cbind(1:nrow(data), ix)]


  # Weight trimming
  cap <- quantile(iptws, trim_weights)
  iptws[iptws > cap] <- cap

  if(mode == "iptw") {

    return(iptws)

  } else if(mode == "ipmw") {

    return(1/data$prob.miss)

  } else {

    data$wts <- iptws*w

    # if(scale_iptw) {
    #
    #   data$wts <- data$wts/sum(data$wts)*length(data$wts)
    #
    # }

    # Estimate blip parameters
    bmod <- lm(outcome.mod,
               weights = wts,
               data = data,
               subset = outmiss==0)

    res <- coef(bmod)

    return(res[which(stringr::str_detect(names(res), trt.name))])

  }

}
