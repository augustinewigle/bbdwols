#' Execute 1 Bayesian Bootstrap iteration of DWOLS given data and Dirichlet weights w
#' @description
#' Do 1 DWOLS with Bayesian Bootstrap weights. Designed for use with the fwb R package.
#' @param data Data.frame of study data
#' @param w vector of weights, one for every row in data, drawn from a Dirichlet (1,...,1) distribution. Supplied by fwb
#' @param outcome.mod A formula for the outcome model using variable names in data. Outcome name extracted from LHS
#' @param trt.mod A formula for modelling treatment assignment using variable names in data. Treatment assignment extracted from LHS
#' @param maxit Maximum number of iterations to allow for fitting outcome model using nnet
#' @param missingoutcome Logical for whether weights that consider the missing outcomes should be used. Otherwise, just uses standard DWOLS IPT weights.
#' @param outmiss.mod A formula for modelling missing outcomes using variables in data. Missing outcomes indicated with NA in response column
#' @param scale_vars Logical for if variables should be scaled to approximately between 0,1 when fitting treatment assignment model to improve convergence of nnet
#' @param trim_weights Percentile at which to trim weights. Default is 1 (no trimming).
#' @param mode Character string indicating what to return. Default is "parameter" which returns blip parameter posterior samples.
#' Other options useful for debugging or assessing suitability of models are "iptw" which returns inverse probability of treatment, "ipmw" which returns inverse probability of being missing,
#' and "allweights" which returns the weights used in finding the blip parameter estimates.
#' @export
#' @importFrom rlang f_lhs
#' @importFrom rlang f_rhs
#' @importFrom stats reformulate
#' @importFrom nnet multinom
#' @importFrom stringr str_detect
#' @import stats
#' @examples
#'
#' if(require("fwb")) {
#' # Create example data
#' example_dat <- data.frame(X1 = runif(200),
#'                           X2 = rnorm(200, 8),
#'                           trt = sample(c("1", "2", "3"), 200, replace = TRUE))
#' example_dat$outcome <- 5 + 0.1*example_dat$X1 + 0.2*example_dat$X2 +
#'                          as.numeric(example_dat$trt) +
#'                          as.numeric(example_dat$trt)*example_dat$X1
#'
#' # BBDWOLS
#' result <- fwb(data = example_dat,
#'               statistic = do_one_bbdwols,
#'               outcome.mod = outcome ~ X2 + trt*X1,
#'               trt.mod = trt ~ 1,
#'               R = 199,
#'               missingoutcome = FALSE,
#'               mode = "parameter")
#' convcheck(result, tix = c(1,2))
#'
#' }
#'
do_one_bbdwols <- function(data,
                        w,
                        outcome.mod,
                        trt.mod,
                        maxit = 500,
                        missingoutcome = FALSE,
                        outmiss.mod = NULL,
                        scale_vars = TRUE,
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

    outmiss.mod <- reformulate(paste(rlang::f_rhs(outmiss.mod)), "outmiss")

    outmissm <- glm(outmiss.mod, family = quasibinomial(),
                    data = data,
                    weights = w)

    data$prob.miss <- outmissm$fitted.values

  } else {

    data$prob.miss <- 1

  }

  tmod.data <- data[,colnames(data)!=outcome.name]

  if(scale_vars) {

    # scale numeric data columns so they are between 0 and 1 for faster convergence in nnet

    ind <- sapply(tmod.data, is.numeric)
    tmod.data[ind] <- lapply(tmod.data[ind], function(x) (x-min(x))/(max(x)-min(x)))

  }

  tmod.data$w <- w

  # Get treatment model weights

  tm <- nnet::multinom(trt.mod,
                       data = tmod.data,
                       maxit = maxit,
                       weights = w,
                       trace = F, model = TRUE)


  if(tm$convergence !=0) {

    warning(paste0("Treatment model did not converge"))

  }

  pis <- predict(tm, type = "probs") # conditional trt probs for all trts

  if(length(unique(data[,trt.name])) == 2) { # binary treatment, need to process probabilities into matrix

    pis <- cbind(pis, 1-pis)
    colnames(pis) <- c(tm$lev[2], tm$lev[1])

  }

  ix <- match(data[,trt.name], colnames(pis)) # Index of treatment they actually received

  iptws <- 1/pis[cbind(1:nrow(data), ix)]

  if(mode == "iptw") {

    return(iptws)

  } else if(mode == "ipmw") {

    return(1/data$prob.miss)

  } else {

    data$rawwts <- iptws/data$prob.miss

    # Weight trimming
    cap <- quantile(data$rawwts, trim_weights)
    data$rawwts[data$rawwts > cap] <- cap

    data$wts <- data$rawwts*w


    if(mode == "allweights") {

      return(data$rawwts/sum(data$rawwts))

    }

    # Estimate blip parameters
    bmod <- lm(outcome.mod,
               weights = wts,
               data = data,
               subset = outmiss==0)

    res <- coef(bmod)

    return(res[which(stringr::str_detect(names(res), trt.name))])

  }

}
