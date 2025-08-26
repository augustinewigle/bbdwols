#' Assess convergence of Bayesian Bootstrap for a bbdwols object
#' @param obj A bbdwols object
#' @param ... additional arguments to `par`
#' @returns Empty
#' @export

convcheck <- function(obj, ...) {

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
