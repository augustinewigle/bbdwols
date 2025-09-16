#' Bayesian Bootstrap dWOLS
#' @param outcome.mod A formula for the outcome model
#' @param trt.mod A Formula for the treatment model
#' @param dat the data.frame with all data, must have variables for outcome and treatment model
#' @param trt.name A string giving the name of the outcome in dat
#' @param L The number of bootstrap samples to take
#' @param diag Logical indicating if diagnostic plots to assess convergence of bootstrap algorithm should be made, default TRUE
#' @param full_posterior Logical indicating if full posterior samples should be returned, or just point estimates. Default is FALSE
#' @importFrom MCMCpack rdirichlet
#' @importFrom nnet multinom
#' @importFrom dplyr contains select
#' @export
#' @return A data.frame with columns corresponding to parameters in the outcome model and each row is an iteration of the bbdwols

bbdwols <- function(outcome.mod,
                    trt.mod,
                    dat, L,
                    trt.name,
                    maxit = 500,
                    diag = TRUE,
                    full_posterior = FALSE,
                    ...) {

  res <- i.bbdwols(outcome.mod = outcome.mod,
                   trt.mod = trt.mod,
                   dat = dat, L = L,
                   trt.name = trt.name,
                   maxit = maxit)

  # Remove parameters that aren't in the blip function

  blip <- dplyr::select(res, dplyr::contains(trt.name))

  # Process results to give point estimates and standard errors

  estimates <- apply(blip, 2, mean)
  names(estimates) <- colnames(blip)

  # ses <- apply(blip, 2, sd)
  # names(ses) <- colnames(blip)

  # df <- matrix(c(estimates, ses), ncol = 2, dimnames = list(names(estimates), c("Estimate", "Std. Error")))

  if(diag) {

    parix <- ifelse(ncol(res) < 12, ncol(res), 12)

    convcheck(res[,sample(1:ncol(res), parix)], ...)

    convergence <- res

  } else {

    convergence <- NULL

  }

  varcovar <- cov(blip, blip)

  if(full_posterior) {

    return(list(coefficients = estimates,
                varcovar = varcovar,
                full_posterior = blip))

  } else {

    return(list(coefficients = estimates,
                varcovar = varcovar,
                full_posterior = FALSE))

  }



}
