
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

# Stubs to Fortran routines for math functions

plpdfFunMath2 <- function(x,nbin,routine) {
   # general function for Math functions
   #    all math functions have the same function format, the only
   #    difference is the to be called Fortran routine

   # determine nbin
   nbin = plpdfNbin(nbin,x)

   # copy data
   xval = as.double(x$x)
   xden = as.double(x$y)
   nx   = as.integer(length(xval))
   nz   = nbin + 1L

   # arguments: xval,xden,nx,zval,zden,zcum,nz,exitcode
   out <- .Fortran(routine,PACKAGE="plPDF",
                               xval=xval,xden=xden,nx=nx,
                               zval=double(nz),zden=double(nz),zcum=double(nz),nz=nz,
                               exitcode=as.integer(9999))

   # output
   z = list(x=out$zval,y=out$zden,Y=out$zcum)
   z = as.plpdf(z)
   attr(z,"exitcode") = out$exitcode

   # return
   return(z)
}

#plpdfInv <- function(x,nbin=NA){return(plpdfFunInv2(x,nbin))}  # alias for rvFunInv2()
plpdfFunInv2 <- function(x,nbin=NA) {
   # calc the inverse (reciprocal) of X
   # ----------------------------------
   # perform     : z = 1/x
   # return value: z

   # call function
   z = plpdfFunMath2(x=x,nbin=nbin,routine="rvfuninv2")

   # return
   return(z)
  }

#plpdfExp <- function(x,nbin=NA){return(plpdfFunExp2(x,nbin))}  # alias
plpdfFunExp2 <- function(x,nbin=NA) {
   # calc the inverse (reciprocal) of X
   # ----------------------------------
   # perform     : z = 1/x
   # return value: z

   # call function
   z = plpdfFunMath2(x=x,nbin=nbin,routine="rvfunexp2")

   # return
   return(z)
}

#plpdfLn <- function(x,nbin=NA){return(plpdfFunLn2(x,nbin))}  # alias
plpdfFunLn2 <- function(x,nbin=NA) {
   # calculate the natural logarithm of RV x
   # ---------------------------------------
   # perform     : z = ln(x)
   # return value: z

   # call function
   z = plpdfFunMath2(x=x,nbin=nbin,routine="rvfunln2")

   # return
   return(z)
}

#plpdfSq <- function(x,nbin=NA){return(plpdfFunSq2(x,nbin))}  # alias
plpdfFunSq2 <- function(x,nbin=NA) {
   # calculate the square of PDF x
   # -----------------------------
   # perform     : z = x^2
   # return value: z

   # call function
   z = plpdfFunMath2(x=x,nbin=nbin,routine="rvfunsq2")

   # return
   return(z)
}

#plpdfSqrt <- function(x,nbin=NA){return(plpdfFunSqrt2(x,nbin))}  # alias
plpdfFunSqrt2 <- function(x,nbin=NA)
  {# calc the squareroot of X
   # ----------------------------------
   # perform     : z = sqrt(x)
   # return value: z

   # call function
   z = plpdfFunMath2(x=x,nbin=nbin,routine="rvfunsqrt2")

   # return
   return(z)
  }

