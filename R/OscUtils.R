#
#
#
#

#' @export
#' @title Oscillation Probability Envelope
#' @param e: vector or energy points
#' @param func: oscillation probability function
#' @param ...: list of arguments to func
#' @description calculate for a given range of values of one of the oscillation parameters the
#' envelope containing the max and min oscillation probability, as a function of the energy
#' @return a list containing \code{x} and \code{y} coordinates of the envelope polygon
prob.env <- function(e, func, ...) {

  func <- match.fun(func)

  m <- sapply(e, function(x) func(x,...))

  #upper and lower envelope line for each value of "e"
  pmin <- apply(m,2, min)
  pmax <- apply(m,2, max)

  emax <- e
  emin <- e

  #invert order of lower envelope line
  emin <- rev(emin)
  pmin <- rev(pmin)

  #join points
  x <- c(emax,emin)
  y <- c(pmax,pmin)

  #add a last point equal to the first one in order to close the envelope
  x <- c(x,x[1])
  y <- c(y,y[1])


  l <- list(x=x,y=y)

}



#' @export
#' @title Histogram Reweight with Oscillation Probability
#' @description weight histogram counts with value of function \code{func} evaluated at the center of the bins (i.e.:\code{mids} vector in \link{hist})
#' @param \code{histo}\cr
#' @param \code{func} weight function. \code{histo$mids} vector is passed as first argument
#' @param \code{...} additional arguments to \code{func}
#' @param opt: other configuration options
#' @return weighted histogram. \code{histo} obj contains a vector with the used wieghts
hist.osc <- function(histo, func, ..., opt="") {

  func <- match.fun(func)

  e <- histo$mids

  w <- func(e,...)

  histo$counts <- histo$counts*w

  histo$oscw <- w

  invisible(histo)
}


#' @export
#' @title transition \code{plotmath} expression
#' @description quote the \code{plotmath} expression for the neutrino oscillation transition \code{nu_in} -> \code{nu_out}
#' @param Standard flavour codes of the initial and final states  (e:1,mu:2,tau:3; negative for antiparticles)
nu.trans.exp <- function(nu_in, nu_out){


  nu_in <- as.integer(nu_in)
  nu_out <- as.integer(nu_out)

  if(nu_in==1 & nu_out==1) return( quote(nu[e]%->%nu[e]) )
  if(nu_in==1 & nu_out==2) return( quote(nu[e]%->%nu[mu]) )
  if(nu_in==1 & nu_out==3) return( quote(nu[e]%->%nu[tau]) )

  if(nu_in==2 & nu_out==1) return( quote(nu[mu]%->%nu[e]) )
  if(nu_in==2 & nu_out==2) return( quote(nu[mu]%->%nu[mu]) )
  if(nu_in==2 & nu_out==3) return( quote(nu[mu]%->%nu[tau]) )

  if(nu_in==3 & nu_out==1) return( quote(nu[tau]%->%nu[e]) )
  if(nu_in==3 & nu_out==2) return( quote(nu[tau]%->%nu[mu]) )
  if(nu_in==3 & nu_out==3) return( quote(nu[tau]%->%nu[tau]) )


  if(nu_in==-1 & nu_out==-1) return( quote(bar(nu)[e]%->%bar(nu)[e]) )
  if(nu_in==-1 & nu_out==-2) return( quote(bar(nu)[e]%->%bar(nu)[mu]) )
  if(nu_in==-1 & nu_out==-3) return( quote(bar(nu)[e]%->%bar(nu)[tau]) )

  if(nu_in==-2 & nu_out==-1) return( quote(bar(nu)[mu]%->%bar(nu)[e]) )
  if(nu_in==-2 & nu_out==-2) return( quote(bar(nu)[mu]%->%bar(nu)[mu]) )
  if(nu_in==-2 & nu_out==-3) return( quote(bar(nu)[mu]%->%bar(nu)[tau]) )

  if(nu_in==-3 & nu_out==-1) return( quote(bar(nu)[tau]%->%bar(nu)[e]) )
  if(nu_in==-3 & nu_out==-2) return( quote(bar(nu)[tau]%->%bar(nu)[mu]) )
  if(nu_in==-3 & nu_out==-3) return( quote(bar(nu)[tau]%->%bar(nu)[tau]) )

  return( quote(bar(nu)[alpha]%->%bar(nu)[beta]) )
}


