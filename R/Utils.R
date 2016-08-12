



#' @export
#' @title  Replace \code{NULLs}
#' @details Replace \code{NULL} with a default value
#' @param a is the value to be checked, b is the substitute value
`%on.null%` <- function(a, b) if (is.null(a)) b else a


#' @export
#' @title Replace \code{NAs}
#' @details Replace \code{NA} with a default value
#' @param a is the value to be checked, b is the substitute value
`%on.na%` <- function(a, b) if (is.na(a)) b else a

#' @export
#' @title Replace \code{NAs} or \code{NLLs}
#' @details Replace \code{NA} or \code{NULL} with a default value
#' @param a is the value to be checked, b is the substitute value
`%on.nn%` <- function(a, b) if (is.null(a) | is.na(a)) b else a


#' @export
#' @title rm \code{NULL}
#' @description Remove \code{NULLs} from a \code{list}
#' @param modified \code{list}
rm.null <- function(x) {
    x[!vapply(x, is.null, logical(1))]
}

#' @export
#' @title rm edges values
#' @description first and last element of a \code{list}
#' @param subsetted \code{list}
rm.edges <- function(x) {
  x[c(-1,length(x))]
}


#' @export
#' @title remove duplicated list-element names
#' @description remove duplicates element names in a list: only first occurrence is retained
#' @param \code{list}
#' @return \code{list} with removed duplicates
rm.duplnames <- function( list ) {

  isd <- duplicated(names(list))

  list[!isd]
}


#' @export
#' @title compute array mid-points
#' @description compute n-1 mid-points of an array of n elements
#' @param \code{x} \code{numeric} list of values
#' @return list of n-1 elements
mid.points <- function(x) {

  invisible( (x[-1]+x[-lentgh(x)])/2  )

}


#' @export
#' @title check equidistat values
#' @description check if values in the array are equally spaced
#' @param numeric vector
#' @return \code{logical}
are.equidist <- function(x) {

  b <- length( unique(diff(x)) ) == 1

  invisible(b)
}


#' @export
#' @title new window
#' @description open a new window for graphics.
#' For OSX it calls quartz function
#' @param \code{title}: window title
new.win <- function(title, opt="") {

  pkg <- .Platform$pkgType

  if( any(grepl("mac.", pkg)) )   return( quartz(title) )
  else stop("OS system not found")

  #use Sys.info alternatively
}



#' @export
#' @title log10-spaced values
#' @description compute a vector of \code{n} log10-spaced values in a given range
#' @details used when you need finer sampling of functions which are varying more rapidly
#' in the lower part of the sampling range
#' @param \code{n}: number of points
#' @param \code{min,max}: interval range
#' @param \code{force.first}: put first element equal to \code{min} (last element is always equal to \code{max})
logspaced <- function(n, min, max, force.first=F){

  if(min<=0 | max<=0 | max<=min | n<=1) stop(call.=T," invalid inputs")

  s <- log10(max/min)/n

  x <- s*(1:n)

  e <- min*10^x

  if(force.first) e[1]<-min

  invisible(e)
}



#' @export
#' @title Check if function is compiled
#' @description check if the function hase been byte-compiled (Returns \code{TRUE} if it finds the text "bytecode:")
#' @param \code{func} function to check
#' @return \code{logical}
is.compiled <- function(func)
{
  # this function lets us know if a function has been byte-coded or not
  #If you have a better idea for how to do this - please let me know...
  if(class(func) != "function") stop("You need to enter a function")
  last_2_lines <- tail(capture.output(func),2)
  any(grepl("bytecode:", last_2_lines)) # returns TRUE if it finds the text "bytecode:" in any of the last two lines of the function's print

}



#' @export
#' @title expand grid
#' @description just a wrapper around \code{\link{expand.grid}} for naming consistency
#' performs combination of all the elements in the set of lists provided via \code{...}
#' @param \code{...}: set of \code{lists} to be combined (as in \code{\link{expand.grid}})
#' @return a \code{data.frame}
grid.df <- function(..., opt=NULL) {

  len <- length(list(...))

  l <- list(NULL)

  if(len ==1) l <- as.list(...)
  else l <- list(...)

  expand.grid(l)

}



#' @export
#' @title \code{\link{expand.grid}} as \code{list}
#' @description performs combination of all the elements in the set of lists provided via \code{...}
#' and store each combiantion into a \code{list} (\code{\link{expand.grid}} produce a \code{data.frame} instead, check also \code{\link{grid.df}})
#' if \code{pre} and \code{post} \code{list}s are provided they are prepended or appended to each combination \code{list}
#' @return a \code{list} or \code{list}s
#' @param \code{...}: set of \code{lists} to be combined (as in \code{\link{expand.grid}})
#' @param \code{pre}: prepend this \code{list} to each combination
#' @param \code{post}: append this \code{list} to each combination
#' @param \code{print}: dump all the sets of combined parameter (N.B.: params in \code{pre} and \code{post} are not shown)
grid.list <- function(pre, post, ..., print=F) {

  l <- list(...)



  df <- grid.df(l)

  if( print ){
    print(df)
  }

  dl <- df.to.list(df, level=2)


  if( !missing(pre) ){

    if( !is.list(pre) ) stop(call.=T, " argument pre is not a list")

    dl <- lapply( dl, function(l) c(pre, l)  )

  }

  if( !missing(post) ){

    if( !is.list(post) ) stop(call.=T, " argument post is not a list")

    dl <- lapply( dl, function(l) c(l, post) )
  }

  invisible(dl)

}


#' @export
#' @title data frame to by-row-list
#' @description convert a \code{data.frame} (or similar) into a \code{list} of rows, prior each row is converted \code{list} \cr
#' For default value of \code{level} (1) From an input \code{data.frame} each element of the \code{list} is a one-row \code{data.frame}.\cr
#' If \code{level}=2, the \code{data.frame} row itself is split into a \code{list} of elements, one of each column.
#' @param \code{d}: \code{data.frame} or matrix-like data structure
#' @return a \code{list} of row-\code{list}s
df.to.list <- function( d , level=1) {

  dl <- lapply(as.list(1:dim(d)[1]), function(x) d[x[1],])

  if( level==2 ){
   #temporary fix
   dl <- lapply( dl, function(l) c(l,list())  )
  }

  invisible(dl)

  # alt.
  #nrow <- nrow(d)
  #
  # l <- list()
  #
  # for(r in 1:nrow) {
  #
  #   l <- c(l, list(d[r,]) )
  #
  # }
  #
  # invisible(l)

}



#' @export
#' @title from a \code{list} of rows to a \code{\link{data.frame}}
#' @description each element of the input \code{list} can be one-row \code{data.frame}
#' or a \code{list} of element each representing a column of the final \code{data.frame}
#' @return a \code{data.frame}
list.to.df <- function( list ) {

  if(length(list) == 0) {
    stop(call.=T, " emply list")
  }

  if( is.data.frame(list[[1]]) ) {

    dl <- plyr:::list_to_dataframe(list, attr(list, "split_labels")  )
    invisible(dl)

  } else {

    dl <- lapply(list, as.data.frame)
    dl <- plyr:::list_to_dataframe(dl, attr(dl, "split_labels")  )
    invisible(dl)

  }

}#list.to.df



#' @export
#' @title function grid-scan
#' @description Perform grid-scan of the function for a subset of its parameters. Vectorizing the input function it's possible to compare the
#' input function shape on range of the x (y) axis for a different combination of the parameters.
#' @param \code{grid.pars}: list of list of values to be scanned for the selected parameters \cr
#' @param \code{fun}: function to be evaluated. Function should be vectorized in order to evaluate the function on each point of
#' the grid (\code{grid.pars}) and for each element of the additional parameter (list \code{...}) passed as a vector.
#' @param \code{...}: other function parameters. A parameter can be passed as a vector (\code{fun} should be vectorized), in this case
#' an output is produced for each value of the input vactor (check the \code{return} section)
#' @param \code{progress}: argument of \code{alply} function in \code{plyr} package
#' @param \code{parallel}: run on multiple cores. Check \code{alply} function in \code{plyr} package.
#' @return a \code{data.frame} with the input parameters and the corresponding function values. A separate line is produced for each
#' point of the \code{grid.pars} and for each value of the parameters passed as vector via the \code{...} list.
fun.grid.scan <- function( grid.pars, fun=NULL, ... , .progress="none", .parallel=FALSE ){


   df <- do.call("expand.grid",grid.pars)

   if ( ! all(names(list(...)) %in% names(formals(fun))) ) {

      #temporary
      stop(call.=T, " cannot pass additional parameters nested in a named list")

   }else{
      # Fix this
   }


   f.df <- function(v,...)
   {

     plist <- c(v,list(...))

     #data.frame(c(v,list(...)),value=do.call(fun,c(v,list(...))))

     data.frame( plist, value=do.call(fun,plist) )
   }

   res <- alply(df,1,f.df, ..., .progress=.progress, .parallel=.parallel)

   df <- list.to.df(res)

}#



#
#
edply <- function(vars, fun, ...)
{
  elply <- function(vars,fun,...,.progress="none",.parallel=FALSE)
  {
    df <- do.call("expand.grid",vars)
    if (all(names(vars) %in% names(formals(fun))))
    {
      #We assume that fun takes the variables in vars as named arguments
      funt <- function(v,...)
      {
        do.call(fun,c(v,list(...)))
      }
      res <- alply(df,1,funt,...,.progress=.progress,.parallel=.parallel)
    }
    else
    {
      #We assume that fun takes a named list as first argument
      res <- alply(df,1,fun,...,.progress=.progress,.parallel=.parallel)
    }
    res
  }

  res <- elply(vars, fun, ...)
  plyr:::list_to_dataframe(res,attr(res, "split_labels"))
}


#' @export
#' @title Read as data frame
#' @param file: data file name
#' @param \code{check.firstCol}: remove first col if all entries are \code{NA}
#' @param "...": inputs of \code{\link{read.delim}} othern than file name
#' @description wrapper around read.delim to add the option of removing
#' the first column if all entries are \code{NA}
#' may be used to fix txt files written out with a blanck at the
#' beginning of each row (like files wittern via C++ cout)
dataRead.asFrame <- function(file, check.firstCol = T, ...) {

  data <- read.delim(file = file, ...)

  if (check.firstCol) {

    # remove first colums if all NA's
    nevents <- nrow(data)

    col.1 <- data[,1]

    nna <- sum(is.na(col.1))

    if (nna == nevents) {
      #data <- data[, 2:ncol(data)]
      data[,1] <- NULL
    }
  }

  invisible(data)
}



