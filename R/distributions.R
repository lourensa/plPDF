
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

# functions to generate PL-PDFs of standard distributions

#' Create a Piecewise Linear uniform Distribution
#'
#' Create a Piecewise Linear uniform Distribution
#' @param min  lower limit of the distribution
#' @param max  upper limit of the distribution
#' @export
plUnif <- function(min=0,max=1) {

   # init
   n = 2   # number of discretization points

   # create distribution
   pdf = list(x=seq(from=min,to=max,length.out=n),
              y=rep(1/(max-min),n))

   # assign class
   pdf = as.plpdf(pdf)

   # return
   return(pdf)
}

#' Create a Piecewise Linear trapezium shaped Distribution
#'
#' Create a Piecewise Linear trapezium shaped Distribution.
#' @param x  array of ordinate values of length \code{n}, with \code{n >= 3}
#' @details The densities of the outermost values of \code{x} will be 0, the density of the other values will be equal and greater than 0. If \code{x[1]==x[2]} then \code{x[1]} is dropped, if \code{x[n]==x[n-1]} then \code{x[n]} is dropped.
#' @export
plTrap <- function(x=c(0,1,2)) {

   # check
   n = length(x)
   if (n < 3) {
      stop("length of x should be at least 3")
   }

   # create distribution
   y    = rep(1,n)   # density=1, normalization is performed in as.plpdf()
   y[1] = y[n] = 0
   if (x[n] == x[n-1]) {
      x = x[-n]
      y = y[-n]
   }
   if (x[1] == x[2]) {
      x = x[-1]
      y = y[-1]
   }
   pdf = list(x=x,y=y)

   # assign class
   pdf = as.plpdf(pdf)

   # return
   return(pdf)
}

#' Create a Piecewise Linear normal Distribution
#'
#' Create a Piecewise Linear normal (Gaussian) Distribution
#' @param mean mean value of the distribution
#' @param sd   standard deviation of the ditribution
#' @export
plNorm <- function(mean=0,sd=1) {

   # init
   n = 51   # number of discretization points
   nsigma = 5

   # create distribution
   x = seq(-sd*nsigma,sd*nsigma,length.out=n) + mean
   pdf = list(x=x)
   pdf$y = dnorm(x,mean=mean,sd=sd)
   pdf$Y = pnorm(x,mean=mean,sd=sd)

   # assign class
   pdf = as.plpdf(pdf)

   # return
   return(pdf)
}

#' Create a Piecewise Linear log-normal Distribution
#'
#' Create a Piecewise Linear log-normal Distribution
#' @param meanlog mean value of the distribution on the log scale
#' @param sdlog   standard deviation of the ditribution on the log scale
#' @export
plLnorm <- function(meanlog=0,sdlog=1) {

   # use plnorm
   pdf = plNorm(mean=meanlog,sd=sdlog)
   pdf = exp(pdf)

   # assign class
   pdf = as.plpdf(pdf)

   # return
   return(pdf)
}


