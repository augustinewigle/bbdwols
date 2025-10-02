#' Bayesian Bootstrap dWOLS
#' @param outcome.mod A formula for the outcome model
#' @param trt.mod A Formula for the treatment model
#' @param dat the data.frame with all data, must have variables for outcome and treatment model
#' @param missingoutcome Logical indicating if probability of missingness weighting should be used to account for missing outcome data. If FALSE, observations with missing outcomes are removed for the analysis.
#' @param outmiss.mod A formula describing the explanatory variables to include in the outcome missingness model, of the format ~ (explanatory variables)
#' @param L The number of bootstrap samples to take
#' @param diag Logical indicating if diagnostic plots to assess convergence of bootstrap algorithm should be made, default TRUE
#' @param full_posterior Logical indicating if full posterior samples should be returned, or just point estimates. Default is FALSE
#' @param alpha 1-confidence level for percentile-based credible intervals. Default to 0.05, ie 95% credible intervals
#' @importFrom MCMCpack rdirichlet
#' @importFrom nnet multinom
#' @importFrom dplyr contains select
#' @importFrom rlang f_lhs f_rhs
#' @importFrom survey svydesign svyglm
#' @export
#' @return A data.frame with columns corresponding to parameters in the outcome model and each row is an iteration of the bbdwols

bbdwols <- function(outcome.mod,
                    trt.mod,
                    dat, L,
                    missingoutcome,
                    outmiss.mod = NULL,
                    maxit = 500,
                    diag = TRUE,
                    full_posterior = FALSE,
                    alpha = 0.05,
                    ...) {

  trt.name <- paste(f_lhs(trt.mod))

  res <- i.bbdwols(outcome.mod = outcome.mod,
                   trt.mod = trt.mod,
                   dat = dat, L = L,
                   trt.name = trt.name,
                   missingoutcome = missingoutcome,
                   outmiss.mod = outmiss.mod,
                   maxit = maxit)

  # Remove parameters that aren't in the blip function

  blip <- dplyr::select(res, dplyr::contains(trt.name))

  # Process results to give point estimates and standard errors

  estimates <- apply(blip, 2, mean)
  names(estimates) <- colnames(blip)

  # Compute CIs
  lower <- apply(blip, 2, function(x) quantile(x, alpha/2))
  upper <- apply(blip, 2, function(x) quantile(x, 1-alpha/2))

  ests <- data.frame(parname = names(estimates),
                     mean = estimates,
                     lower = lower,
                     upper = upper)


  if(diag) {

    parix <- ifelse(ncol(res) < 12, ncol(res), 12)

    convcheck(res[,sample(1:ncol(res), parix)], ...)

    convergence <- res

  } else {

    convergence <- NULL

  }

  varcovar <- cov(blip, blip)

  if(full_posterior) {

    return(list(results = ests,
                varcovar = varcovar,
                full_posterior = blip))

  } else {

    return(list(results = ests,
                varcovar = varcovar,
                full_posterior = FALSE))

  }



}
