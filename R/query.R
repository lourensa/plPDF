
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


#' Query functions for PL-PDF's
#'
#' These functions provide information about the plpdf distribution. \code{dplpdf} gives the density, \code{pplpdf} gives the distribution function, \code{qplpdf} gives the quantile function and \code{rplpdf} generates random deviates.
#' @describeIn dplpdf
#' Get the probability densities of a PDF given the values
#' @param pdf   object of class \code{plpdf}
#' @param x     vector of values. Non-finite values are ignored.
#' @export
dplpdf <- function(x,pdf) {

   out = plpdfFunValDen(pdf,x=x)

   # return
   return(out$y)
}

# Get probabilities given the quantiles
#
# Get the probabilities of a PDF given the quantiles
# @inheritParams dplpdf

#' @describeIn dplpdf
#' Get the probabilities of a PDF given the quantiles
#' @param q     vector of quantiles. Non-finite values are ignored.
#' @export
pplpdf <- function(q,pdf) {

   out = plpdfFunValCum(pdf,x=q)

   # return
   return(out$Y)
}

# Get values at cumulative probabilities 
# 
# @inheritParams dplpdf

#' @describeIn dplpdf
#' Get values at cumulative probability values
#' @param p     vector of probabilities. Non-finite values are ignored.
#' @export
qplpdf <- function(p,pdf) {


   out = plpdfFunCumVal(pdf,Y=p)

   # return
   return(out$x)
}

# Draw a number of values from a PDF
#
# @inheritParams dplpdf

#' @describeIn dplpdf
#' Draw a number of values from a PDF
#' @param n     number of observations. If length(n) > 1, the length is taken to be the number required.
#' @export
rplpdf <- function(n,pdf) {
   # init
   if (length(n) > 1) {n = length(n)}

   # generate n probabilities of a standard uniform distribution
   p = runif(n)

   # get quantiles
   q = qplpdf(p,pdf)

   # return
   return(q)
}


