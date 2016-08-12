#' Stat utils
#'
#'example stopifnot(!missing(x), is.matrix(x), dim(x)==c(2,2), x>0)


#' @export
#' @title List Statistical Descriptive Quantities
#' @param x: input array of values \cr
#' quant: "quantile around the median", see return values for a description.
#' @return
#' Mean, Median,SD (standard deviation), Min, Max \cr
#' Qr/Ql: value which limit quant/2 portion of the sample on the right/left of the median \cr
#' it means that in the range [Qr,Ql] it is contained quant fraction of the sample \cr
#' while each of the ranges [Ql,Median] and [Median,Qr] contains quant/2 fraction of the sample\cr
#' Qb/Qe: 5\% and 95\% quantile
stat.band <- function(x, quant=0.68) {

    Mean <- mean(x, na.rm = TRUE)
    Median <- median(x, na.rm = TRUE)
    SD <- sd(x, na.rm = TRUE)

    Min <- min(x, na.rm = TRUE)
    Max <- max(x, na.rm = TRUE)

    xx <- x[x > Median]
    Qr <- quantile(xx, quant)
    xx <- x[x <= Median]
    Ql <- quantile(xx, 1 - quant)

    #change default suffix to quantile values to reflect the actual meaning
    names(Ql)[1] <- paste(as.character(quant*100),"%",sep="")
    names(Qr)[1] <- paste(as.character(quant*100),"%",sep="")

    Qb <- quantile(x,0.05)
    Qe <- quantile(x,0.95)

    l <- c(Mean = Mean, SD = SD, Min = Min, Max = Max, Median = Median, Ql = Ql, Qr = Qr, Qb=Qb, Qe=Qe)
}


#' @export
#' @title Statistical Fluctuation
#' @description apply a Poisson (default) or Gaussian Statistical Fluctuation \cr
#' to the input number n.\cr
#' For Gauss the variation is sampled from a gaussian with sigma=sqrt(n) and mean 0\cr
#' For Poisson the variation is sampled from a poissonian with lambda=sqrt(n)
#' @param value to be flucuated
#' @param option ("Poisson", "Gauss")
#' @return fluctuated value
stat.fluct <- function(n, opt = "Pois") {

    if (n < 0) {
        return(0)
    }


    if (opt == "Gaus") {

        dn <- rnorm(1, mean = 0, sd = sqrt(n))
        n <- n + dn

    }else{

        n <- rpois(1, lambda = n)*1.0
    }


    if (n > 0) return(n)
    else return(0)
}


#' @export
#' @title Histogram Statistical Fluctuation Matrix
#' @description generate bin-by-bin statistical fluctuations
#' @param histo: input histogram, statistical fluctuations are calculate around each bin count
#' @param n: number of repetitions
#' @param opt: type of fluctuation, es. "Gauss" or "Pois", check \link{stat.fluct} function for details
#' @return a matrix with a colum for each histogram bin, each column as "n" values corresponding to the
#' "n" throws
hist.fluct.m <- function(histo, n=1000, opt="") {

  y <- histo$counts

  m <- replicate(n,sapply(y, stat.fluct, opt = opt))

  #transpose in order to have a column for each bin
  m <- t(m)

}



#' @export
#' @title Histogram Statistical Fluctuation Boxplot
#' @description gnerate bin-by-bin statistical fluctuations calculated from \link{hist.fluct.m} function
#' @return a boxplot-like object with desciptive
#' quantities for each bin.
#' The fivenum-style descriptive values are the one returned from stat.band function: (N.B.: they are not the usual boxplot ones):\cr
#' Qb,Ql,Meadia,Qr,Qe.
hist.fluct.b <- function(histo, n=1000, opt="") {

  m <- hist.fluct.m(histo,n,opt)

  sm <- apply(m, 2, stat.band)

  smm <- apply(sm,2, function(x) {
                           l <- as.list(x);
                           c(Qb=l$Qb, Ql=l$Ql, Median=l$Median, Qr=l$Qr, Qe=l$Qe);
                          }
                          )


  b <- boxplot(m,plot=F)

  b$stats <- smm

  #rm unnecessary structures
  b$out <- NULL
  b$group <- NULL
  b$conf <- NULL

  return(b)
}



#' @export
#' @title Histogram Statistical Fluctuation Envelop
#' @description gnerate bin-by-bin statistical fluctuations calculated from \link{hist.fluct.m} function
#' @return points of a close polygon deliminting the envelope
#' @details upper an lower points of the envelope are given by the value Qr and Ql defined in the \link{stat.band} function
hist.fluct.env <- function(histo, n=1000, opt="", ...) {


  b <- hist.fluct.b(histo,n,opt,...)

  #extracts lower/upper bounds (Ql and Qr)
  lo <- apply(b$stats,2,function(x) as.list(x)$Ql)
  up <- apply(b$stats,2,function(x) as.list(x)$Qr)

  breaks <- histo$breaks

  #create histogram object
  hlo <- list(counts=lo,density=lo,breaks=breaks)
  hup <- list(counts=up,density=up,breaks=breaks)

  #calculate upper and lower lineshapes
  loxy <- hist.shxy(hlo)
  upxy <- hist.shxy(hup)

  #remove first and last point
  loxy$x <- loxy$x[c(-1,-length(loxy$x))]
  loxy$y <- loxy$y[c(-1,-length(loxy$y))]

  upxy$x <- upxy$x[c(-1,-length(upxy$x))]
  upxy$y <- upxy$y[c(-1,-length(upxy$y))]

  #reverse order for lower line
  loxy$x <- rev(loxy$x)
  loxy$y <- rev(loxy$y)

  #join upper and lower
  x <- c(upxy$x,loxy$x)
  y <- c(upxy$y,loxy$y)

  #add last-point=first-point to close the envelope
  x <- c(x,x[1])
  y <- c(y,y[1])

  l <- list(x=x,y=y)

}







