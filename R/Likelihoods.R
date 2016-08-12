#
## Utilities for computation of likelihoods and test statistics
#


#' @export
#' @title LogLike pseudo
#' @param hobs
LogLike.PsedoExp <- function(vobs, data_exp) {



    v <- apply(data_exp, 1, LogLike.Binned, vobs = vobs)
}



#' @export
#' @title LogLike of binned data
#' @description compute negative loglike
#' @param hobs: binned observed data
#' @param hexp: binned expected data
LogLike.binned <- function(vobs, vexp, opt="") {


  LL <- vobs - vexp - vobs * log(vobs/vexp)

  LL[!is.finite(LL)] <- 0

  LL <- sum(LL, na.rm = TRUE)

  invisible(-2 * LL)
}


#' @export
#' @title LogLike of two histogram
#' @description compute negative loglike
#' @param hobs: binned observed data
#' @param hexp: binned expected data
LogLike.histo <- function(hobs, hexp, opt="") {

  LL <- LogLike.binned(hobs$counts,hexp$counts,opt)

  invisible(LL)

}


#' @export
#' @title Chi2 of binned data
#' @param hobs: binned observed data
#' @param herr: error of observed data bin
#' @param hexp: binned expected data
Chi2.binned <- function(vobs, verr, vexp) {

  # rm zero err entries
  chi2 <- (vobs - vexp)^2/(verr^2)

  chi2[is.infinite(chi2)] <- 0

  chi2 <- 0.5 * sum(chi2, na.rm=TRUE )

}



#' @export
#' @title Chi2 of two histograms
#' @param hobs: binned observed data
#' @param herr: error of observed data bin
#' @param hexp: binned expected data
Chi2.histo <- function(hobs, herr, hexp) {

    chi2 <- Chi2.binned(hobs$counts, herr$counts, hexp$counts)

    invisible(chi2)

}




