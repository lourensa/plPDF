
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

# define stub routines for Fortran source of plpdfFun-routines
# =========================================================


#' Get values at cumulative probability values
#'
#' Get values at cumulative probability values.
#' @param pdf    object of type \code{plpdf}.
#' @param Y      cumulative probability values to find an ordinate value for.
plpdfFunCumVal <- function(pdf,Y=NA) {

   # get length of arrays
   nx  = as.integer(length(pdf$x))
   sel = is.finite(Y)
   nv  = as.integer(sum(sel))
   x   = rep(as.double(NA),length(Y))

   exitcode = as.integer(9999)

   # call to Fortran routine
   #    arguments: xval,xden,nx,zval,zden,zcum,nz,exitcode
   if (any(sel)) {
      out <- .Fortran("rvfuncumval",PACKAGE="plPDF",
                                    xval=as.double(pdf$x),
                                    xden=as.double(pdf$y),
                                    nx=nx,
                                    vval=double(nv),
                                    vcum=as.double(Y[sel]),
                                    nv=nv,
                                    exitcode=exitcode)
      # return value
      x[sel]    = out$vval
      exitcode  = out$exitcode
   }

   # output
   ret = list(x=x,Y=Y)
   attr(ret,"exitcode") = exitcode

   # return
   return(ret)
}

#' Get densities at values of PL-PDF
#'
#' Get densities at values of PL-PDF.
#' @param pdf    object of type \code{plpdf}.
#' @param x      ordinate values to find density values for.
plpdfFunValDen <- function(pdf,x=NA) {
   # get probability density of pdf at x
   # -----------------------------------
   # return value: v, structure calculated densities

   # get length of arrays
   nx  = as.integer(length(pdf$x))
   sel = is.finite(x)
   nv  = as.integer(sum(sel))
   y   = rep(as.double(NA),length(x))

   exitcode = as.integer(9999)

   # call to Fortran routine
   if (any(sel)) {
      out <- .Fortran("rvfunvalden",PACKAGE="plPDF",
                                    xval=as.double(pdf$x),
                                    xden=as.double(pdf$y),
                                    nx=nx,
                                    vval=as.double(x[sel]),
                                    vden=double(nv),
                                    nv=nv,
                                    exitcode=exitcode)
      # return value
      y[sel]    = out$vden
      exitcode  = out$exitcode
   }

   # output
   ret = list(x=x,y=y)
   attr(ret,"exitcode") = exitcode

   # return
   return(ret)
}


#' Get cumulative probabilities at values of PL-PDF
#'
#' Get cumulative probabilities at values of PL-PDF.
#' @param pdf    object of type \code{plpdf}.
#' @param x      ordinate values to find density values for.
plpdfFunValCum <- function(pdf,x=NA) {
   # get cumulative probability of x at val
   # --------------------------------------------------------------
   # return value: z

   # get length of arrays
   nx  = as.integer(length(pdf$x))
   sel = is.finite(x)
   nv  = as.integer(sum(sel))
   Y   = rep(as.double(NA),length(x))

   exitcode = as.integer(9999)

   # call to Fortran routine
   if (any(sel)) {
      out <- .Fortran("rvfunvalcum",PACKAGE="plPDF",
                                    xval=as.double(pdf$x),
                                    xden=as.double(pdf$y),
                                    nx=nx,
                                    vval=as.double(x[sel]),
                                    vcum=double(nv),
                                    nv=nv,
                                    exitcode=exitcode)
      # return value
      Y[sel]    = out$vcum
      exitcode  = out$exitcode
   }

   # output
   ret = list(x=x,Y=Y)
   attr(ret,"exitcode") = exitcode

   # return
   return(ret)
}

#' Get cumulative probabilities and densities at values of PL-PDF
#'
#' Get cumulative probabilities and densities at values of PL-PDF.
#' @param pdf    object of type \code{plpdf}.
#' @param x      ordinate values to find density values for.
plpdfFunValDenCum <- function(pdf,x=NA) {
   # get cumulative probability and densities of x at val
   # --------------------------------------------------------------
   # return value: z

   # get length of arrays
   nx  = as.integer(length(pdf$x))
   sel = is.finite(x)
   nv  = as.integer(sum(sel))
   y   = rep(as.double(NA),length(x))
   Y   = rep(as.double(NA),length(x))

   exitcode = as.integer(9999)

   # call to Fortran routine
   if (any(sel)) {
      out <- .Fortran("rvfunvaldencum",PACKAGE="plPDF",
                                    xval=as.double(pdf$x),
                                    xden=as.double(pdf$y),
                                    nx=nx,
                                    vval=as.double(x[sel]),
                                    vden=double(nv),
                                    vcum=double(nv),
                                    nv=nv,
                                    exitcode=exitcode)
      # return value
      y[sel]    = out$vden
      Y[sel]    = out$vcum
      exitcode  = out$exitcode
   }

   # output
   ret = list(x=x,y=y,Y=Y)
   attr(ret,"exitcode") = exitcode

   # return
   return(ret)
}

