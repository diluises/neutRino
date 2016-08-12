#
#
#' @import fields
#' @import weights
#' @import ggplot2
#' @import compiler

#' @export
#' @title 1D Histogram
#' @description General 1D (weighted) histogram structure
#' @param \code{x}: list of entries
#' @param \code{w}: list of weights
#' @param \code{breaks}: bins edges
#' @param \code{show}: plot the histo
#' @param \code{main}: histo title
#' @param \code{xlab,ylab}: axis titles
#' @param \code{...}: other graphic options
#' @return obj elements:\cr
#' @return \code{nentries}
#' @return \code{outlayers}
#' @return \code{counts}
#' @return \code{breaks}
#' @return \code{nbins}
#' @return \code{binw}: bin width
#' @return \code{errs}: bin error
#' @return \code{equidist}
#' @return \code{main,xlab,ylab}
#' @return \code{opts}: arguments in \code{...}
#' @return \code{call}
histo1d <- function (x=NULL, weight=NULL, breaks, na.rm = TRUE,
                     show = TRUE,
                     main, xlab, ylab, ...)
{


  opts <- list(...)

  nbins <- length(breaks) - 1

  nentries <- length(x)

  if( is.null(x) ) {
    x <- rep(0,nbins)
  }


  if( na.rm ){
    x <- x[!is.na(x)]
  }

  # associate entries to bin
  idx <- cut(x, breaks, include.lowest = TRUE)

  #count outlayers
  outlayers <- length(idx[is.na(idx)])



  v <- vector()

  if( is.null(weight) ){
    v <- tapply(x,idx,length)
  }else{
    v <- tapply(weight,idx,sum)
  }

  v[is.na(v)] <- 0

  #
  #
  #

  if( missing(main) ){
    main <- "histogram"
  }

  if( missing(xlab) ){
    xlab <- "x"
  }

  if( missing(ylab) ){
    ylab <- "Entries"
  }


  #
  ## Build Object
  #
  obj <- list()

  obj$nentries  <- nentries
  obj$outlayers <- outlayers

  obj$counts <- as.vector(v)
  obj$breaks <- breaks
  obj$nbins  <- length(breaks) -1
  obj$binw   <- diff(breaks)
  obj$errs   <- histo1d.errs(as.vector(v))

  obj$equidist <- are.equidist(breaks)

  obj$main <- main
  obj$xlab <- xlab
  obj$ylab <- ylab

  obj$opts <- opts
  obj$call <- match.call()

  class(obj) <- "histo1d"
  #
  ##
  #

  if( show ){

    plot.histo1d(obj,main,xlab,ylab,...)
  }


  #
  ##
  #
  invisible(obj)

}#histo1d




#' @export
#' @title plot \code{histo1d}
#' @description Plot \code{histo1d} by means of the \code{plot.hist} function
#' @param \code{histo1d}: a histo1d obj
#' @param \code{main,xlab,ylab}: arguments of \code{plot}
#' @param \code{...}: other graphics options
plot.histo1d <- function(histo1d,main,xlab,ylab,...){

  hist <- to.hist(histo1d)

  if(missing(main)) main <- histo1d$main
  if(missing(xlab)) xlab <- histo1d$xlab
  if(missing(ylab)) ylab <- histo1d$ylab

  plot(hist, main=main, xlab=xlab, ylab=ylab, ...)

}#plot.histo1d



#' @export
#' @title convert \code{histo1d} to \code{histogram}
#' @description build a \code{histogram} obj from a \code{histo1d} one
#' @return a \code{hist} obj
to.hist <- function(histo1d) {

  h <- list()

  h$breaks  <- histo1d$breaks

  h$counts  <- as.vector(histo1d$counts)
  h$density <- as.vector(histo1d$counts)

  h$mids  <- histo1d$mids

  h$xname <- histo1d$xlab

  h$equidist <- histo1d$equidist

  class(h) <- "histogram"

  invisible(h)

}#to.hist




#' @export
#' @title Compute  1D histogram bin erros
#' @description compute error of each bin content, default is to take the \code{sqrt} of the
#' bin content
#' @param vector of bin contents
#' @return vector of bin errors
histo1d.errs <- function(x, opt=NULL) {

  errs <- sapply(x, function(x) {
                       if(!is.finite(x) | x<=0) return (0)
                       else return (sqrt(x))
                       }  )

  invisible(errs)

}#histo1d.errs




#' @export
#' @title 2D Histogram
#' @description General 2D (weighted) histogram structure
#' @param \code{x}: list of entries
#' @param \code{w}: list of weights
#' @param \code{breaks}: bins edges
#' @param \code{show}: plot the histo
#' @param \code{main}: histo title
#' @param \code{xlab,ylab}: axis titles
#' @param \code{col}
#' @param \code{...}: other graphic options
#' @return obj elements:\cr
#' @return \code{nentries}
#' @return \code{outlayers}
#' @return \code{counts}
#' @return \code{breaks}
#' @return \code{nbins}
#' @return \code{binw}: bin width
#' @return \code{errs}: bin error
#' @return \code{equidist}
#' @return \code{main,xlab,ylab}
#' @return \code{opts}: arguments in \code{...}
#' @return \code{call}
histo2d <- function (x, y = NULL, weight=NULL, breaks.x, breaks.y, na.rm = TRUE
                     ,main,xlab, ylab,
                     show = TRUE, legend=TRUE, col = c("white", heat.colors(12))
                     ,...)
{


  opts <- list(...)

  if (is.null(y)) {
    if (ncol(x) != 2) stop(" if y is ommitted, x must be a 2 column matrix")
    y <- x[,2]
    x <- x[,1]
  }


  if ( na.rm ) {
    nas <- is.na(x) | is.na(y)
    x <- x[!nas]
    y <- y[!nas]
  }

  #associate entries to bins
  idx <- cut(x, breaks.x, include.lowest = TRUE)
  idy <- cut(y, breaks.y, include.lowest = TRUE)


  #compute bin counts
  m <- matrix()

  if(is.null(weight)){

    m <- tapply(x, list(idx, idy), length)

  }else{

    m <- tapply(weight, list(idx, idy), sum)

  }

  m[is.na(m)] <- 0


  #  #if (identical(FUN, base::length))

  if (missing(main))
    main <- deparse(substitute(main))
  if (missing(xlab))
    xlab <- deparse(substitute(xlab))
  if (missing(ylab))
    ylab <- deparse(substitute(ylab))


  midpoints <- function(x) (x[-1] + x[-length(x)])/2


  #
  ## Build object
  #
  obj <- list()

  obj$nentries = length(x)
  obj$outlayers = obj$nentries - sum(m)

  obj$counts <- m
  obj$breaks.x <- breaks.x
  obj$breaks.y <- breaks.y
  obj$breaks <- c(breaks.x,breaks.y)

  obj$errs <- histo2d.errs(m)

  obj$mids.x <- midpoints(breaks.x)
  obj$mids.y <- midpoints(breaks.y)
  obj$mids <- c(midpoints(breaks.x), midpoints(breaks.y) )

  obj$nbins.x <- length(breaks.x)-1
  obj$nbins.y <- length(breaks.y)-1
  obj$nbins <- c(length(breaks.x)-1,length(breaks.y)-1)

  obj$main <- main
  obj$xlab <- xlab
  obj$ylab <- ylab

  obj$opts <- opts
  obj$call <- match.call()

  class(obj) <- "histo2d"

  #
  ##
  #
  if (show){
    plot.histo2d(histo2d=obj,main=main,legend=legend,col=col,xlab=xlab,ylab=ylab,...)
  }


  invisible(obj)

}#histo2d



#' @export
#' @title plot \code{histo2d}
#' @description Plot \code{histo2d} by means of the \code{image} function
#' @param \code{histo2d}: a histo2d obj
#' @param \code{legend}: add color gradied scale for z axis (if \code{TRUE} \code{image.plot} is used)
#' @param \code{col}: heat color map for the \code{legend}
#' @param \code{xlab,ylab}: axis titles
#' @param \code{...}: graphics options
plot.histo2d <- function(histo2d, legend=T, col=c("white", heat.colors(12)), main,xlab, ylab,...) {

    breaks.x <- histo2d$breaks.x
    breaks.y <- histo2d$breaks.y
    counts   <- histo2d$counts


    if(missing(main)) main <- histo2d$main
    if(missing(xlab)) xlab <- histo2d$xlab
    if(missing(ylab)) ylab <- histo2d$ylab

    if(legend){
      image.plot(breaks.x, breaks.y, counts, col = col, xlab = xlab, ylab = ylab,
                 ...)
    }else{
      image(breaks.x, breaks.y, counts , col = col, xlab = xlab, ylab = ylab,
            ...)
    }

    #add main text above the plot region (pos=3)
    text(x=0.5*(breaks.x[1]+breaks.x[length(breaks.x)]),y=max(breaks.y),pos=3,adj=0.5,label=main)

}#plot.histo2d



#' @export
#' @title Compute  2D histogram bin erros
#' @description compute error of each bin content, default is to take the \code{sqrt} of the
#' bin content
#' @param matrix of bin contents
#' @return matrix of bin errors
histo2d.errs <- function(x, opt=NULL) {

  geterr <- function(x){
    if(!is.finite(x) | x<=0) return (0)
    else return (sqrt(x))
  }

  errs <- apply(x, 2, function(x) { sapply(x,geterr)  })

  invisible(errs)

}#histo2d.errs



#' @export
#' @title Histogram Bin Content Text
#' @description show bin content as text on top of each bin \cr
#' @param \code{histo}: \code{histo1d} obj \cr
#' @param \code{opt}= "v" labels are vertical, "o": labels are orizontal
#' @param \code{round}: number of decimals for rounding
#' @param \code{cex}: font size
#' @param \code{srt}: angle
#' @param \code{...}: all the other options of the "text" function
#' @details check also parameter \code{labels} in \code{hist}
bintext.histo1d <- function(histo, opt, angle="v",round=3, cex=0.8, srt=90, ...) {

  x <- histo$mids

  y <- as.numeric(histo$counts)
  l <- as.numeric(histo$counts)

  l <- round(l,round)

  adj <- 1
  pos <- NULL

  if(angle=="v"){
    adj <- -0.4
  }else if(angle=="o") {
    srt <- 0
    pos <- 3
  }

  text(x=x, y=y, labels=l
       ,srt=srt, adj=adj, pos=pos, cex=cex
       ,...)

}#bintext.histo1d



#' @export
#' @title 1D Histo line shape
#' @description calculate line defining the shape of the histogram
#' @param \code{histo1d}
#' @return list with x and y vectors defining the hist shape
lshape.hist1d <- function(histo1d, opt=NULL) {

  x <- histo$breaks
  x <- sapply(x,function(x) c(x,x))
  x <- as.vector(x)


  y <- as.numeric()
  if(opt=="c") y <- histo$counts
  else if(opt=="d") y <- histo$density

  y <- sapply(y,function(x) c(x,x))
  y <- as.vector(y)
  y <- c(0,y,0)

  list(x=x,y=y)

}#lshape.hist1d



#################################################################################################
#
#
#
#################################################################################################



# histo2d.v0
# local version
histo2d.v0 <- function (x, y = NULL, weight=NULL, nbins=200, same.scale = FALSE, na.rm = TRUE,
                        show = TRUE, col = c("white", heat.colors(12)), FUN = base::length,
                        main,xlab, ylab,x.breaks,y.breaks,legend=T, ...)
{


  if (is.null(y)) {
    if (ncol(x) != 2)
      stop("If y is ommitted, x must be a 2 column matrix")
    y <- x[, 2]
    x <- x[, 1]
  }
  if (length(nbins) == 1)
    nbins <- rep(nbins, 2)
  nas <- is.na(x) | is.na(y)
  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
  }
  else stop("missinig values not permitted if na.rm=FALSE")

  if(same.scale){
    x.cuts = x.breaks;
    y.cuts = x.breaks;
  }else{
    x.cuts <- x.breaks
    y.cuts <- y.breaks
  }

  index.x <- cut(x, x.cuts, include.lowest = TRUE)
  index.y <- cut(y, y.cuts, include.lowest = TRUE)


  m <- matrix()

  if(is.null(weight)){

    m <- tapply(x, list(index.x, index.y), FUN)

  }else{

    m <- tapply(weight, list(index.x, index.y), sum)

  }


  #
  #if (identical(FUN, base::length))
  #  m[is.na(m)] <- 0

  #do it anyhow
  m[is.na(m)] <- 0
  #
  #

  if (missing(main))
    main <- deparse(substitute(main))
  if (missing(xlab))
    xlab <- deparse(substitute(xlab))
  if (missing(ylab))
    ylab <- deparse(substitute(ylab))

  if (show){
    if(legend){
      image.plot(x.cuts, y.cuts, m, col = col, xlab = xlab, ylab = ylab,
                 ...)
    }else{
      image(x.cuts, y.cuts, m, col = col, xlab = xlab, ylab = ylab,
            ...)
    }
  }


  midpoints <- function(x) (x[-1] + x[-length(x)])/2


  retval <- list()

  retval$nentries = length(x)
  retval$outlayers = retval$nentries - sum(m)

  retval$counts <- m
  retval$afreq <- m/retval$nentries
  retval$breaks.x = x.cuts
  retval$breaks.y = y.cuts
  retval$mids.x = midpoints(x.cuts)
  retval$mids.y = midpoints(y.cuts)
  retval$mids = c(midpoints(x.cuts),midpoints(y.cuts))


  retval$nbins.x = length(x.cuts)
  retval$nbins.y = length(y.cuts)
  retval$nbins = c(length(x.cuts),length(y.cuts))


  retval$main = main
  retval$xlab = xlab
  retval$ylab = ylab

  retval$call <- match.call()

  class(retval) <- "histo2d"

  invisible(retval)

}#histo2d.v0

