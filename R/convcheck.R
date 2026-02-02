#' Assess convergence of Bayesian Bootstrap from an fwb object object
#' @param fwbobj Result of running fwb with do_one_bbdwols function
#' @param tix index of blip parameters to plot. Either numeric vector,
#' "random" which will randomly choose 1/4 of parameters to check, or
#' "all" which will plot all parameters
#' @returns Empty
#' @export

convcheck <- function(fwbobj, tix = "random") {

  objtemp <- fwbobj$t

  nblip <- ncol(objtemp)

  if(tix[1] == "random") {

    # if nblip < 6, just plot them all
    if(nblip < 5) {

      message("Less than 5 parameters, plotting all...")
      tix <- 1:nblip

    } else {

      message("Randomly selecting one quarter of parameters...")
      ncheck <- ceiling(nblip*0.25)
      tix <- sample(1:nblip, ncheck)

    }


  } else if (tix[1] == "all") {

    tix <- 1:nblip

  }

  obj <- objtemp[,tix]

  if(length(tix)==1) {

    obj <- matrix(obj, ncol = 1, dimnames = list(NULL, colnames(objtemp)[tix]))

  }

  runningsd <- matrix(nrow = nrow(obj)-1, ncol = ncol(obj))

  for(i in 1:ncol(obj)) {

    for(j in 2:nrow(obj)) {

      runningsd[j-1,i] <- sd(obj[1:j, i])

    }

  }

  runningsd <- data.frame(runningsd, iter = 2:nrow(obj))

  colnames(runningsd)[1:ncol(obj)] <- colnames(obj)

  plotlist <- vector(length = length(tix), mode = "list")

  for(j in 1:ncol(obj)) {

    plotlist[[j]] <- plot(x = runningsd$iter, y = runningsd[,j], xlab = "Iters", ylab = "SD over MC iters",
               main = paste0("Parameter: ", colnames(obj)[j]))

  }

  return(invisible(plotlist))


}
