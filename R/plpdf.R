
#    plPDF is a library to perform calcualtions with piecewise linear random variables.
#
#    Copyright (C) 2016, 2021  Aris Lourens
#
#    This file is part of plPDF.
#
#    plPDF is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    plPDF is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' plPDF: A package for computations with piecewise linear PDFs
#'
#'
#' @docType package
#' @name plPDF
#' @useDynLib plPDF
NULL
#> NULL

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
   } else if (is.null(obj)) {
      out = FALSE
   } else if (! is.list(obj)) {
      out = FALSE
   } else {
      # check content of obj
      out = plpdfCheck(obj)
   }

   # return
   return(out)
}


#' Determine Number of bins
#'
#' Determine Number of bins based on number of bins of PDFs
#' @param  nbin   if nbin has a finite value, this value is returned. If not finte then the returnded value is based on the given PDFs and a default number of bins.
#' @param  pdf1   if object is of class plpdf then its number of bins is used in the evaluation (optional)
#' @param  pdf2   as \code{pdf1}
#' @param  pdf3   as \code{pdf1}
plpdfNbin <- function(nbin,pdf1=NULL,pdf2=NULL,pdf3=NULL) {

   # init
   nbindef = 51   # default number of bins

   # check or nbin is defined
   if (is.finite(nbin)) {
      if (nbin > 0) {
         # to be sure it is an integer value
         nbin = as.integer(ceiling(nbin))
         # return
         return(nbin)
      }
   }

   # determine nbin
   n = 0
   for (p in list(pdf1,pdf2,pdf3)) {
      if (is.plpdf(p)) {n = max(n,length(p$x)-1)}
   }
   nbin = as.integer(ceiling((n+nbindef)/2))

   # return
   return(nbin)
}
