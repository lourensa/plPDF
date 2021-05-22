
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

# methods for generic functions for class plpdf objects

#' Method for is.na function for PL-PDFs
#'
#' Method for generic \code{is.na} function for PL-PDFs of class \code{plpdf}. The density array is tested for \code{NA} values. In fact, \code{NA} values are not allowed in a PL-PDF.
#' @param pdf  object of class \code{plpdf}.
#' @export
is.na.plpdf <- function(pdf) {
   out = is.na(pdf$y)
   # return
   return(out)
}


#' Method for binary operations on PL-PDFs
#'
#' Method for binary operations on PL-PDFs. Implemented operations are \code{+,-,*,/}.
#' @param pdf1  object of class \code{plpdf} or scalar variable
#' @param pdf2  object of class \code{plpdf} or scalar variable
#' @export
Ops.plpdf <- function(pdf1,pdf2) {

   # init

   # calc
   if (.Generic %in% c("+","-","*","/")) {
      # determine number of bins
      nbin = plpdfNbin(nbin=NA,pdf1,pdf2)
      #
      if     (.Generic == "+") {ret = plpdfSum2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "-") {ret = plpdfSub2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "*") {ret = plpdfMul2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "/") {ret = plpdfDiv2(pdf1,pdf2,nbin=nbin)}
#   } else if (.Generic %in% c("==","!=","<","<=",">=",">")) {
   } else {
      stop("Binary operation ",.Generic," not implemented")
   }

   # return
   return(ret)
}

#' Plot method for PL-PDF objects
#'
#' Plot method for PL-PDF objects
#' @param pdf    object of class \code{plpdf}
#' @param xlab   a title for the x axis. See \code{plot}
#' @param ylab   a title for the y axis. See \code{plot}
#' @param type   what type of plot should be drawn, See \code{plot}.
#' @param ...    arguments to be passed to graphical functions
#' @export
plot.plpdf <- function(pdf,xlab="Value",ylab="Probability density",add=FALSE,type="l",...) {

   # init

   if (! add) {
      # initialize graph
      plot(x=range(pdf$x),y=range(c(0,pdf$y)),xlab=xlab,ylab=ylab,type="n",...)
   }

   # draw line
   lines(x=pdf$x,y=pdf$y,type=type,...)

   # return
   return(invisible())
}


#' Math method for PL-PDF objects
#'
#' Math method for PL-PDF objects.
#' The implemented generic functions are \code{log}, \code{exp}, \code{sqrt}, \code{sq}.
#' @export
Math.plpdf <- function(pdf,...) {

   if (.Generic %in% c("log","ln")) {
      out = plpdfFunLn2(pdf)
   } else if (.Generic == "exp") {
      out = plpdfFunExp2(pdf)
   } else if (.Generic == "sq") {
      out = plpdfFunSq2(pdf)
   } else if (.Generic == "sqrt") {
      out = plpdfFunSqrt2(pdf)
   } else {
      stop("Funtion ", .Generic, " not implemented")
   }
   # return
   return(out)
}

