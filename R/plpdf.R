
# basic functions for class plpdf

#' Check or object can be a PL-PDF
#'
#' Check or object can be a PL-PDF
#' @param  obj   to be tested object
plpdfCheck <- function(obj) {

   # check obj
   if (! is.list(obj)) {
      warning("Object is of wrong type, should be a list or data frame.")
      return(FALSE)
   }
   if (! all(c("x","y") %in% names(obj))) {
      warning("Object has not enough data columns, should have at least x and y.")
      return(FALSE)
   }
   if (length(obj$x) < 2) {
      warning("PDF should consists of at least 2 data points")
      return(FALSE)
   }
   if (length(obj$x) != length(obj$y)) {
      warning("Lengths of x and y are not equal")
      return(FALSE)
   }
   if (any(is.na(obj$x)) | any(is.na(obj$y))) {
      warning("NA values not allowed")
      return(FALSE)
   }
   if (any(diff(obj$x) < 0)) {
      # x values must increase
      warning("Values of x not increasing")
      return(FALSE)
   }
   if (any(obj$y < 0)) {
      # y values must be >= 0
      warning("Values of x not increasing")
      return(FALSE)
   }

   # return
   return(TRUE)
}


#' Convert a list into a PL-PDF
#'
#' Convert a list into a Piecewise Linear Probability Density Function
#' @param  obj   list; object to be converted into a \code{plpdf}. The list must contain at least the data \code{(x,y)}.
#' @param  normalize  logical; normalize the function so it integrates to 1.
#' @details A piecewise pinear PDF (\code{plpdf}) consists of at least an \code{x} and \code{y} array. Here \code{x} is the value and \code{y} is the probability density. In addition, the array \code{Y} can be defined. This array contains the cumulative probability.
#' @seealso \code{link{is.plpdf}}
#' @author Aris Lourens
#' @export
as.plpdf <- function(obj,normalize=TRUE) {

   # init
   out = list()

   # check object
   ret = plpdfCheck(obj)
   if (! ret) {
      return(NULL)
   }

   # create out
   out$x = obj$x
   out$y = obj$y
#   out$Y = obj$Y

   # normalize
   if (normalize) {
      cum   = sum(diff(out$x)*(head(out$y,-1)+tail(out$y,-1)))/2.
      out$y = out$y/cum
   }

   # assign a class
   class(out) = "plpdf"

   # return
   return(out)
}


#' Test object being of class plpdf
#'
#' Test object being of class plpdf. 
#' @param  obj   to be tested object
#' @author Aris Lourens
#' @seealso \code{link{as.plpdf}}
#' @export
is.plpdf <- function(obj) {

   # init
   out = FALSE

   # test
   if ("plpdf" %in% class(obj)) {
      out = TRUE
   } else {
      # check content of obj
      out = plpdfCheck(obj)
   }

   # return
   return(out)
}


