
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
#' Method for binary operations on PL-PDFs. Implemented operations are \code{+,-,*,/,^}.
#' @param pdf1  object of class \code{plpdf} or scalar variable
#' @param pdf2  object of class \code{plpdf} or scalar variable
#' @export
Ops.plpdf <- function(pdf1,pdf2) {

   # init

   # calc
   if (.Generic %in% c("+","-","*","/","^")) {
      # determine number of bins
      nbin = plpdfNbin(nbin=NA,pdf1,pdf2)
      #
      if     (.Generic == "+") {ret = plpdfSum2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "-") {ret = plpdfSub2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "*") {ret = plpdfMul2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "/") {ret = plpdfDiv2(pdf1,pdf2,nbin=nbin)}
      else if(.Generic == "^") {ret = plpdfPow2(pdf1,pdf2,nbin=nbin)}
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
#' @param xlab   a title for the x axis. See \code{\link[base]{plot}}
#' @param ylab   a title for the y axis. See \code{\link[base]{plot}}
#' @param add    logical; If `TRUE`, add the graph to the current plot, if `FALSE`, start a new plot.
#' @param type   what type of plot should be drawn, See \code{\link[base]{plot}}.
#' @param derivative integer; the order derivative which should be plotted. Implemented: 0=PDF (default), -1=CDF.
#' @param prob   numeric; probability interval to limit the extent x-axis. See details.
#' @param ...    arguments to be passed to graphical functions
#' @details `prob` defines a probability interval. The accompanying quantiles of `pdf` are used to limit the exent of the x-axis. 
#'    The interval may be defined by one or two values. If more than two values are given then all but the first two are ignored.
#'    If one value is given this is the symmetric probability interval `[0.5-prob/2,0.5+prob/2]`. 
#'    If two values are given they are the lower and upper probability value.
#'    `NA` values are allowed in which case the respective extent of `pdf` is used.
#'    If `NULL` it is not used. This argument is only used if `add=FALSE`.
#' @export
plot.plpdf <- function(pdf,xlab="Value",ylab=NA,add=FALSE,type="l",derivative=0,prob=NULL,...) {

   # init
   ylab_pdf = "Probability density"
   ylab_cdf = "Cumulative probability"

   # derivative
   if (derivative == -1) {
      # CDF
      # add Y to pdf
      y = pplpdf(pdf$x,pdf)
      if(is.na(ylab)) {ylab=ylab_cdf}
      yrange = range(y)
   } else if (derivative == 0) {
      # PDF
      y = pdf$y
      if(is.na(ylab)) {ylab=ylab_pdf}
      yrange = range(c(0,y))
   } else {
      stop("Unknown derivative value ",derivative)
   }

   # initialize graph
   if (! add) {
      # xrange
      xrange = range(pdf$x)
      if (! is.null(prob)) {
         if (! all(is.na(prob))) {
            if (length(prob) == 1) {
               tprob = c(0.5-prob/2,0.5+prob/2)
            } else {
               tprob = prob[1:2]
            }
            sel = is.na(tprob)
            tprob[sel] = c(0,1)[sel]
            txrange = qplpdf(tprob,pdf)
            xrange[! sel] = txrange[! sel]
         }
      }
      # init graph
      plot(x=xrange,y=yrange,xlab=xlab,ylab=ylab,type="n",...)
   }

   # draw line
   lines(x=pdf$x,y=y,type=type,...)

   # return
   return(invisible())
}


#' Math method for PL-PDF objects
#'
#' Math method for PL-PDF objects.
#' The implemented generic functions are \code{log}, \code{exp}, \code{sqrt}, \code{sq}.
#' @param pdf    object of class \code{plpdf}
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

#' Summary method for PL-PDF objects
#'
#' Summary method for PL-PDF objects.
#' @param pdf    object of class \code{plpdf}
#' @details The implemented generic functions are \code{sum}, \code{min}, \code{max}, \code{range}.
#'
#' \code{sum} returns the integral of the PL-PDF function. If multiple pdf's are given then the total sum is returned.
#'
#' \code{min}, \code{max} and \code{range} return the respective values of the domain (\code{x}-value) of the pdf's.
#' If multiple PL-PDFs are given then the minimum or maximum af all domain values is given.
#' @export
Summary.plpdf <- function(...,na.rm=FALSE) {

   if (.Generic %in% c("sum")) {
      out = 0
      for (pdf in list(...)) {
         s   = sum(diff(pdf$x)*(head(pdf$y,-1)+tail(pdf$y,-1)))/2.
         out = out + s
      }
   } else if (.Generic %in% c("min")) {
      out = NULL
      for (pdf in list(...)) {
         out = min(out,pdf$x,na.rm=na.rm)
      }
   } else if (.Generic %in% c("max")) {
      out = NULL
      for (pdf in list(...)) {
         out = max(out,pdf$x,na.rm=na.rm)
      }
   } else if (.Generic %in% c("range")) {
      out = NULL
      for (pdf in list(...)) {
         out = range(out,pdf$x,na.rm=na.rm)
      }
   } else {
      stop("Funtion ", .Generic, " not implemented")
   }
   # return
   return(out)
}


#' Calculate the mean value of a PL-PDF
#'
#' Calculate the mean value of a PL-PDF
#' @param pdf    object of class \code{plpdf}
#' @export
mean.plpdf <- function(pdf) {

   # init
   n   = length(pdf$x)
   ii  = 1:(n-1)

   # calc mean
   out = sum(diff(pdf$x)*((pdf$y[ii]  *(   pdf$x[ii+1]+  2.*pdf$x[ii] ))   
                         +(pdf$y[ii+1]*(2.*pdf$x[ii+1]+     pdf$x[ii] )) ) )
   out = out/6.

   # devide by total probability, to be sure (should be 1)
   prb = 0.5*sum((pdf$y[ii]+pdf$y[ii+1])*diff(pdf$x))
   out = out/prb

   # return
   return(out)
}


# This function length has an unwanted side effect: e.g. str() does not work anymore on class plpdf objects.
# So it must be removed.
##' Get the length of a PL-PDF
##'
##' Get the length of the arrays in a PL-PDF. 
##' @param pdf    object of class \code{plpdf}
##' @return The number of discretization points of the PL-PDF, which is equal to \code{nbin}+1. 
##'         If the object \code{pdf} is not of class \code{plpdf} then \code{NULL} is returned.
##' @export
#length.plpdf <- function(pdf) {
#
#   # init
#   out = NULL
#
#   # get length
#   if (is.plpdf(pdf)) {
#      out = length(pdf$x)
#   }
#
#   # return
#   return(out)
#}


