
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


# define stub routines for Fortran source
# =======================================



# define functions
# ----------------
# Stub functions
# ..............

plpdfBop2 <- function(x,y,nbin=NA,routine,normalize=TRUE) {
   # general function for binary operations
   #    all binary operations have the same function format, the only
   #    difference is the to be called Fortran routine

   # check
   if (missing(routine)) {
      stop("routine name missing")
   }

   # aray lengths
   nx   = length(x$x)
   ny   = length(y$x)
   nbin = plpdfNbin(nbin=nbin,x,y)   
#   if(is.na(nbin)) {
#      nz = max(nx,ny)
#   } else {
      nz = nbin+1
#   }

   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
   out <- .Fortran(routine,PACKAGE="plPDF",
              xval=as.double(x$x),xden=as.double(x$y),nx=as.integer(nx),
              yval=as.double(y$x),yden=as.double(y$y),ny=as.integer(ny),
              zval=as.double(rep(0.,nz)),
              zden=as.double(rep(0.,nz)),
              zcum=as.double(rep(0.,nz)),
              as.integer(nz),
              exitcode=as.integer(0))

   # output
   z = list(x=out$zval,y=out$zden,Y=out$zcum)
   z = as.plpdf(z,normalize=normalize)
   attr(z,"exitcode") = out$exitcode

   # return
   return(z)

}

#' Binary operations of PL-PDF objects
#' 
#' Perform binary operations on two PL-PDF's. Especially useful when the default values of \code{nbin} and \code{normalize} are not suitable. Otherwise the generic binary operation (+,-,*,/,^) can be used.
#' @describeIn plpdfSum
#' Calculate the sum of two PL-PDF's.
#' @param x     object of class \code{plpdf}
#' @param y     object of class \code{plpdf}
#' @param nbin  integer; number of bins in output object. If \code{NA} this number is automatically defined.
#' @param normalize logical; If \code{TRUE} (default) the resulting PDF is normalized. See \code{\link{as.plpdf}}.
#' @export
plpdfSum  <- function(x,y,nbin=NA,normalize=TRUE){
   # alias for plpdfSum2()
   return(plpdfSum2(x,y,nbin,normalize=normalize))
}

#' @describeIn plpdfSum
#' Subtract two PL-PDF's.
#' @export
plpdfSub  <- function(x,y,nbin=NA,normalize=TRUE){
   # alias for plpdfSub2()
   return(plpdfSub2(x,y,nbin,normalize=normalize))
}

#' @describeIn plpdfSum
#' Multiply two PL-PDF's.
#' @export
plpdfMul  <- function(x,y,nbin=NA,normalize=TRUE){
   # alias for plpdfMul2()
   return(plpdfMul2(x,y,nbin,normalize=normalize))
}

#' @describeIn plpdfSum
#' Divide two PL-PDF's.
#' @export
plpdfDiv  <- function(x,y,nbin=NA,normalize=TRUE){
   # alias for plpdfDiv2()
   return(plpdfDiv2(x,y,nbin,normalize=normalize))
}

#' @describeIn plpdfSum
#' Calculate the power of two PL-PDF's
#' @export
plpdfPow <- function(x,y,nbin=NA,normalize=TRUE) {
   # alias for plpdfPow2()
   return(plpdfDiv2(x,y,nbin,normalize=normalize))
}


plpdfSum2 <- function(x,y,nbin=NA,normalize=TRUE) {
   # summarize PL-PDFs of independent variables x and y
   # --------------------------------------------------
   # perform     : z = x + y
   # return value: z

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopsum2",normalize=normalize)
   } else {
      if (lx) {
         z = x
         if (is.numeric(y)) {
            z$x = z$x + y
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = y
            z$x = z$x + x
         } else if (is.numeric(y)) {
            z = x + y  # should not occur
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

# -------------------------------

plpdfSum2R <- function(x,y,nclass=NA,normalize=TRUE) {
   # summarize PDFs of independent variables x and y
   # -----------------------------------------------
   # perform     : z = x + y
   # return value: z

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopsum2r",normalize=normalize)
   } else {
      if (lx) {
         z = x
         if (is.numeric(y)) {
            z$x = z$x + y
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = y
            z$x = z$x + x
         } else if (is.numeric(y)) {
            z = x + y  # should not occur
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

# -------------------------------

plpdfSub2 <- function(x,y,nbin=NA,normalize=TRUE) {
   # summarize PDFs of independent variables x and y
   # -----------------------------------------------
   # perform     : z = x - y
   # return value: z

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopsub2",normalize=normalize)
   } else {
      if (lx) {
         z = x
         if (is.numeric(y)) {
            z$x = z$x - y
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = plpdfMul2(y,-1.)
            z$x = z$x + x
         } else if (is.numeric(y)) {
            z = x - y  # should not occur
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

# -------------------------------

plpdfSub2R <- function(x,y,nbin=NA,normalize=TRUE) {
   # Subtract PDFs of independent variables x and y, RAW version
   # TEMPORARY FUNCTION, TO BE REPLACED BY FORTRAN EQUIVALENCE
   # --------------------------------------------
   # perform     : z = x - y
   # return value: z

   iy = list(x=rev(-y$x),y=rev(y$y))
   z  = plpdfSum2R(x,iy,nbin,normalize=normalize)

   # return
   return(z)
}

# -------------------------------

plpdfMul2 <- function(x,y,nbin=NA,normalize=TRUE) {
   # Multiply PDFs of independent variables x and y
   # ----------------------------------------------
   # perform     : z = x * y
   # return value: z

   locRev <- function(v) {
      # reverse order of arrays
      v   = unclass(v)
      attr(v,"exitcode") = NULL
      v$x = rev(v$x)
      v$y = rev(v$y)
      if ("Y" %in% names(v)) {
         v$Y = 1-rev(v$Y)
      }
      return(v)
   }

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopmul2",normalize=normalize)
   } else {
      if (lx) {
         z = x
         if (is.numeric(y)) {
            if (y < 0) {z=locRev(z)}
            z$x = z$x * y
            z$y = z$y / abs(y)
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = y
            if (x < 0) {z=locRev(z)}
            z$x = z$x * x
            z$y = z$y / abs(x)
         } else if (is.numeric(y)) {
            z = x / y
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

# -------------------------------

plpdfMul2R <- function(x,y,nbin=NA,normalize=TRUE) {
   # Multiply PDFs of independent variables x and y
   # ----------------------------------------------
   # perform     : z = x * y
   # return value: z

   locRev <- function(v) {
      # reverse order of arrays
      v   = unclass(v)
      attr(v,"exitcode") = NULL
      v$x = rev(v$x)
      v$y = rev(v$y)
      if ("Y" %in% names(v)) {
         v$Y = 1-rev(v$Y)
      }
      return(v)
   }

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopmul2r",normalize=normalize)
   } else {
      if (lx) {
         z = x
         if (is.numeric(y)) {
            if (y < 0) {z=locRev(z)}
            z$x = z$x * y
            z$y = z$y / abs(y)
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = y
            if (x < 0) {z=locRev(z)}
            z$x = z$x * x
            z$y = z$y / abs(x)
         } else if (is.numeric(y)) {
            z = x / y
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

# -------------------------------

plpdfDiv2 <- function(x,y,nbin=NA,normalize=TRUE) {
   # Divide PDFs of independent variables x and y
   # --------------------------------------------
   # perform     : z = x / y
   # return value: z

   locRev <- function(v) {
      # reverse order of arrays
      v   = unclass(v)
      attr(v,"exitcode") = NULL
      v$x = rev(v$x)
      v$y = rev(v$y)
      if ("Y" %in% names(v)) {
         v$Y = 1-rev(v$Y)
      }
      return(v)
   }

   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (lx & ly) {
      z = plpdfBop2(x=x,y=y,nbin=nbin,routine="rvbopdiv2",normalize=normalize)
   } else {
      if (lx) {
         if (is.numeric(y)) {
            z = plpdfMul2(x,1/y,nbin)
         } else {
            stop("Wrong type of argument y")
         }
      } else if (is.numeric(x)) {
         if (ly) {
            z = plpdfFunInv2(y,nbin=nbin)
            z = plpdfMul2(x,z,nbin)
         } else if (is.numeric(y)) {
            z = x / y
         } else {
            stop("Wrong type of argument y")
         }
      } else {
         stop("Wrong type of argument x")
      }
      z = as.plpdf(z,normalize=normalize)
   }

   # return
   return(z)
}

#' Calucalate the power of a PDF
#'
#' Interim function for calculation of the power of a PL-PDF, \code{x^y}. Fortran implementation coming later.
#' The function is currently implemented as: z=exp(y*log(x))
plpdfPow2 <- function(x,y,nbin=NA,normalize=TRUE) {

   # nbin
   nbin = plpdfNbin(nbin=nbin,x,y)   
   
   # check type of arguments
   lx = is.plpdf(x)
   ly = is.plpdf(y)
   if (! lx) {if (! is.numeric(x)) {stop("Argument x should be a PL-PDF or numeric")}}
   if (! ly) {if (! is.numeric(y)) {stop("Argument y should be a PL-PDF or numeric")}}

   # check lx value
   if (lx) {
      if (any(x$x < 0)) {
         stop("Power function not valid for x-values less than 0")
      } else if (is.numeric(x)) {
         if (x < 0) {
            stop("Power function not valid for x less than 0")
         }
      }
   }

   # calc
   if (lx) {
      lnx  = plpdfFunLn2(x,nbin=nbin)
   } else {
      lnx = log(x)
   }
   ylnx = plpdfMul2(y,lnx,nbin=nbin)
   z    = plpdfFunExp2(ylnx,nbin=nbin)

   # return
   return(z)
}

###xxx###rvBopDiv2R <- function(x,y,nclass=NA) {
###xxx###   # Divide PDFs of independent variables x and y, RAW VERSION
###xxx###   # --------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: z
###xxx###
###xxx###   # aray lengths
###xxx###   nx   = length(x$x)
###xxx###   ny   = length(y$x)
###xxx###   if(is.na(nbin))
###xxx###     {nz = max(nx,ny)}
###xxx###   else
###xxx###     {nz = nbin+1}
###xxx###
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopdiv2r",PACKAGE="plPDF",
###xxx###              xval=as.double(x$d.val),xden=as.double(x$d.prob),nx=as.integer(nx),
###xxx###              yval=as.double(y$d.val),yden=as.double(y$d.prob),ny=as.integer(ny),
###xxx###              zval=as.double(rep(0.,nz)),
###xxx###              zden=as.double(rep(0.,nz)),
###xxx###              zcum=as.double(rep(0.,nz)),
###xxx###              as.integer(nz),
###xxx###              exitcode=as.integer(0))
###xxx###
###xxx###   pdfz = list(d.val=out$zval,d.prob=out$zden,c.prob=out$zcum,attr=list(exitcode=out$exitcode))
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx###rvBopDiv2Abs <- function(x,y,nclass=NA) {
###xxx###   # Divide PDFs of independent variables x and y
###xxx###   # --------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: z
###xxx###
###xxx###   # aray lengths
###xxx###   nx   = length(x$x)
###xxx###   ny   = length(y$x)
###xxx###   if(is.na(nbin))
###xxx###     {nz = max(nx,ny)}
###xxx###   else
###xxx###     {nz = nbin+1}
###xxx###
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopdiv2abs",PACKAGE="plPDF",
###xxx###              xval=as.double(x$d.val),xden=as.double(x$d.prob),nx=as.integer(nx),
###xxx###              yval=as.double(y$d.val),yden=as.double(y$d.prob),ny=as.integer(ny),
###xxx###              zval=as.double(rep(0.,nz)),
###xxx###              zden=as.double(rep(0.,nz)),
###xxx###              zcum=as.double(rep(0.,nz)),
###xxx###              as.integer(nz),
###xxx###              exitcode=as.integer(0))
###xxx###
###xxx###   pdfz = list(d.val=out$zval,d.prob=out$zden,c.prob=out$zcum,attr=list(exitcode=out$exitcode))
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopSum3 <- function(x,y,nclass=NA,z=NA,mxden=FALSE)
###xxx###  {# summarize PDFs of independent variables x and y
###xxx###   # -----------------------------------------------
###xxx###   # perform     : z = x + y
###xxx###   # return value: z
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="sum",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopsum3",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopSum3R <- function(x,y,nclass=NA,z=NA,mxden=FALSE)
###xxx###  {# summarize PDFs of independent variables x and y
###xxx###   # -----------------------------------------------
###xxx###   # perform     : z = x + y
###xxx###   # return value: z
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="sum",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopsum3r",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopSub3 <- function(x,y,nclass=NA,z=NA,mxden=FALSE)
###xxx###  {# subtract PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x - y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="sub",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopsub3",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopSub3R <- function(x,y,nclass=NA,z=NA,mxden=FALSE) {
###xxx###   # Subtract PDFs of independent variables x and y, RAW version
###xxx###   # TEMPORARY FUNCTION, TO BE REPLACED BY FORTRAN EQUIVALENCE
###xxx###   # --------------------------------------------
###xxx###   # perform     : z = x - y
###xxx###   # return value: z
###xxx###
###xxx###   iy   = list(d.val=rev(-y$d.val),d.prob=rev(y$d.prob))
###xxx###   pdfz = rvBopSum3R(x,iy,nclass=nclass,z=z,mxden=mxden)
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopSub3tmp <- function(x,y,nclass=NA)
###xxx###  {# subtract PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x - y
###xxx###   # return value: z
###xxx###
###xxx###   # make use of summation function
###xxx###   ix        = (length(y$d.val):1)
###xxx###   yy        = list()
###xxx###   yy$d.val  = -1*y$d.val[ix]
###xxx###   yy$d.prob = y$d.prob[ix]
###xxx###
###xxx###   pdfz = rvBopSum3(x,yy,nclass)
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopMul3 <- function(x,y,nclass=NA,z=NA,mxden=FALSE)
###xxx###  {# Multiply PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x * y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="mul",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopmul3",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopMul3R <- function(x,y,nclass=NA,z=NA,mxden=FALSE)
###xxx###  {# Multiply PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x * y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="mul",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopmul3r",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopMul3Q <- function(x,y,nclass=NA,z=NA,q=1)
###xxx###  {# Multiply PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x * y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="mul",nclass=nclass)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   routine=paste("srvbopmulq",q,sep="")
###xxx###   for(k in 1:nz){
###xxx###      z =zval[k]
###xxx###      zd=as.double(0)
###xxx###      zc=as.double(0)
###xxx###      out <- .Fortran(routine    ,xval,xden,nx,
###xxx###                                  yval,yden,ny,
###xxx###                                  z=z,zd=zd,zc=zc,exitcode=exitcode,PACKAGE="plPDF")
###xxx###      zden[k]=out$zd
###xxx###      zcum[k]=out$zc
###xxx###   exitcode    = out$exitcode
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode,k,z))}
###xxx###     }
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = zden
###xxx###   pdfz$c.prob = zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###  }
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopDiv3 <- function(x,y,nclass=NA,z=NA,mxden=FALSE) {
###xxx###   # Divide PDFs of independent variables x and y
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="div",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopdiv3",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopDiv3R <- function(x,y,nclass=NA,z=NA,mxden=FALSE) {
###xxx###   # Divide PDFs of independent variables x and y, RAW VERSION
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="div",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopdiv3r",xval,xden,nx,
###xxx###                               yval,yden,ny,
###xxx###                               zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###rvBopDiv3Abs <- function(x,y,nclass=NA,z=NA,mxden=FALSE) {
###xxx###   # Divide PDFs of independent variables x and y, RAW VERSION
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: z
###xxx###
###xxx###
###xxx###   # create Z values
###xxx###   if (any(is.na(z))){z = spdfLinZ(x,y,type="div",nclass=nclass,mxden=mxden)}
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z)
###xxx###   nz   = as.integer(length(zval))
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = as.double(rep(0,nz))
###xxx###   zcum     = zden
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran("rvbopdiv3abs",xval,xden,nx,
###xxx###                                  yval,yden,ny,
###xxx###                                  zval,zden=zden,zcum=zcum,nz,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   pdfz        = list()
###xxx###   pdfz$d.val  = zval
###xxx###   pdfz$d.prob = out$zden
###xxx###   pdfz$c.prob = out$zcum
###xxx###   exitcode    = out$exitcode
###xxx###
###xxx###   if(exitcode != 0){print(paste("Exitcode:",exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(pdfz)
###xxx###}
###xxx###
###xxx#### -------------------------------
###xxx###
###xxx###srvBopDivQn <- function(x,y,z,rtn) {
###xxx###   # Multiple quadrant probability routines for division
###xxx###   # ----------------------------------------------
###xxx###   # perform     : z = x / y
###xxx###   # return value: list(zval,zden,zcum)
###xxx###
###xxx###
###xxx###   # routine choice
###xxx###   #   0 subroutine srvBopDivQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #   1 subroutine srvBopDivQ1(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #   2 subroutine srvBopDivQ2(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #   3 subroutine srvBopDivQ3(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #   4 subroutine srvBopDivQ4(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #  10 subroutine srvBopDivQ0Abs(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #  11 subroutine srvBopDivQ1Abs(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #  12 subroutine srvBopDivQ2Abs(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #  13 subroutine srvBopDivQ3Abs(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###   #  14 subroutine srvBopDivQ4Abs(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 
###xxx###
###xxx###   # init
###xxx###   routines = c("srvbopdivq0"   ,"srvbopdivq1"   ,"srvbopdivq2"   ,"srvbopdivq3"   ,"srvbopdivq4",
###xxx###                NA,NA,NA,NA,NA,
###xxx###                "srvbopdivq0abs","srvbopdivq1abs","srvbopdivq2abs","srvbopdivq3abs","srvbopdivq4abs")
###xxx###   routine  = routines[rtn+1]
###xxx###
###xxx###   # copy data
###xxx###   xval = as.double(x$d.val)
###xxx###   xden = as.double(x$d.prob)
###xxx###   nx   = as.integer(length(xval))
###xxx###   yval = as.double(y$d.val)
###xxx###   yden = as.double(y$d.prob)
###xxx###   ny   = as.integer(length(yval))
###xxx###   zval = as.double(z[1])
###xxx###
###xxx###   # create data memory for output variables
###xxx###   zden     = double(1)
###xxx###   zcum     = double(1)
###xxx###   exitcode = as.integer(9999)
###xxx###
###xxx###   # arguments: xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode
###xxx###   out <- .Fortran(routine,xval,xden,nx,
###xxx###                           yval,yden,ny,
###xxx###                           zval=zval,zden=zden,zcum=zcum,exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy returned values
###xxx###   out     = out[c("zval","zden","zcum","exitcode")]
###xxx###   out$rtn = routine
###xxx###
###xxx###   if(out$exitcode != 0){print(paste("Exitcode:",out$exitcode))}
###xxx###
###xxx###   # return
###xxx###   return(out)
###xxx###}
###xxx###
###xxx#### -------------------------------

plpdfMean <- function(pdf1) {
   # calculate mean of PDF
   #
   # Return:
   #   mean

   # calc mu
   out=NULL
   out <- .Fortran("rvstatmean",PACKAGE="plPDF",
                   as.double(pdf1$x),as.double(pdf1$y),as.integer(length(pdf1$x)),
                   mu=as.double(-9e300),
                   exitcode=as.integer(0))
   mu   = out$mu

   # return
   return(mu)
}

# -------------------------------

plpdfStat <- function(pdf1) {
   # calculate mean and standard deviation of PDF
   #
   # Return:
   #   c(mean,standard deviation)

   if (! plpdfCheck(pdf1)) {
      print(pdf1)
      warning("rvStat error in input PDF")
   }

   # calc stat
   out=NULL
   out <- .Fortran("rvstat",PACKAGE="plPDF",
                   as.double(pdf1$x),as.double(pdf1$y),as.integer(length(pdf1$x)),
                   stat=as.double(rep(-9e300,2)),
                   exitcode=as.integer(0))
   out   = out$stat

   # return
   return(out)
}


# -------------------------------

plpdfStatProb <- function(pdf1) {
   # calculate the total probability as described by the densities of the PDF
   #
   # Return:
   #   probability

   # calc mu
   out=NULL
   
   out <- .Fortran("rvstatprob",PACKAGE="plPDF",
                   xval=as.double(pdf1$x),
                   xden=as.double(pdf1$y),
                   nx=as.integer(length(pdf1$x)),
                   prob=as.double(-1.),
                   exitcode=as.integer(0))
   out   = out$prob

   # return
   return(out)
}


# -------------------------------

plpdfClean <- function(pdf1,nsigma=6,minprob=0.999,nscale=1,switch=0) {
   # Clean and truncate a PDF
   #
   #!! switch options:
   #!!
   #!!  1. make 0 <= xcum <= 1
   #!!
   #!! 64. remove empty bins at outside:
   #!!        xcum(1)==xcum(2)  ==0 & xden(1)==xden(2)  ==0
   #!!        xcum(n)==xcum(n-1)==1 & xden(n)==xden(n-1)==0
   #!!
   #!!  2. remove bin borders with no contribution: xcum(i-1)==xcum(i)==xcum(i+1)
   #!!
   #!!  4. truncate ends if PDF width > 2*nsigma*sigma
   #!!        The truncation is performed in such way that the probability densities
   #!!        at both ends are equal (if possible)
   #!!        Further restriction: retained probability must be at least minprob
   #!!                      
   #!!  8. extent end bins
   #!!        if xcum(ixb) > 0 extent bin ixb
   #!!        if xcum(ixe) < 1 extent bin ixe
   #!!
   #!! 16. rescale densities to let probability of density be close to cum-prob
   #!!
   #!! 32. renormalize the PDF
   #!!        set the total probability of the PDF to 1
   #!!        rescale the densities to estimate xcum as close as possible
   #!!
   #!! Multiple switches may be specified through addition: e.g. switch = 1+8+32
   #!! If switch == 0, all options are used               : switch = 1+2+4+8+16+32+64
   #!! If switch < 0 , the specified switched are excluded: switch = 1+2+4+8+16+32+64 + switch
   #!!
   #!! Order of processing: 1,64,2,4,8,16,32
   #!!
   #!! When multiple switch options are given (not summarized) the fortran routine is called
   #!! multiple times. This may be usefull to change the order of processing, or to perform 
   #!! some options multiple times.
   # Return:
   #   cleaned PDF

   # preprocess pdf
   x = as.plpdf(pdf1)
   if(is.null(x)) {
      warning("pdf1 incomplete!")
      print(pdf1)
      return(pdf1)
   }
   if(is.null(x$Y)) {
      nx  = length(x$x)
      cum = diff(x$x)*(x$y[-1]+x$y[-nx])/2.
      x$Y = cumsum(c(0,cum))/sum(cum)
   }

   # values
   xval     = as.double(x$x)
   xden     = as.double(x$y)
   xcum     = as.double(x$Y)
   nx       = as.integer(length(xval))
   nsigma   = as.double(nsigma)
   minprob  = as.double(minprob)
   nscale   = as.integer(nscale)
   switch   = as.integer(switch)
   exitcode = as.integer(0)


   for(sw in switch)
     {out=NULL
      # xval,xden,xcum,nx,nsigma,minprob,switch,exitcode
      out <- .Fortran("rvclean",PACKAGE="plPDF",
                      xval=xval,xden=xden,xcum=xcum,nx=nx,
                      nsigma=nsigma,minprob=minprob,nscale=nscale,switch=sw,
                      exitcode=exitcode)
      nx   = out$nx
      xval = out$xval[1:nx]
      xden = out$xden[1:nx]
      xcum = out$xcum[1:nx]
     }

   # copy result
   x = list(x=xval,y=xden,Y=xcum)
   x = as.plpdf(x)
   attr(x,"exitcode") = out$exitcode

   # return
   return(x)
}

# -------------------------------

###xxx###srvNxtZ <- function(lst)
###xxx###  {# Find next z-value based on probability density and cumulative probality
###xxx###   #
###xxx###   # Arguments
###xxx###   #    lst      list with arrays
###xxx###   #       $val  values
###xxx###   #       $den  probability densities
###xxx###   #       $cum  cumulative probability
###xxx###   # Return
###xxx###   #    lst      list with new arrays for val, den and cum
###xxx###   #       iz    position of inserted value
###xxx###
###xxx###   # copy data and add one extra position
###xxx###   zval = as.double(c(lst$val,0))
###xxx###   zden = as.double(c(lst$den,0))
###xxx###   zcum = as.double(c(lst$cum,0))
###xxx###   nz   = as.integer(length(zval))
###xxx###   iz   = as.integer(-1)
###xxx###   exitcode = as.integer(0)
###xxx###
###xxx###   # call Fortran routine srvNxtZ(zval,zden,zcum,nz,iz,exitcode)
###xxx###   out <- .Fortran("srvnxtz",zval=zval,zden=zden,zcum=zcum,nz=nz,
###xxx###                             iz=iz,
###xxx###                             exitcode=exitcode,PACKAGE="plPDF")
###xxx###
###xxx###   # copy result
###xxx###   lst$val      = out$zval
###xxx###   lst$den      = out$zden
###xxx###   lst$cum      = out$zcum
###xxx###   lst$iz       = out$iz
###xxx###   lst$exitcode = out$exitcode
###xxx###
###xxx###   # return
###xxx###   return(lst)
###xxx###  }

# -------------------------------

plpdfMix2 <- function(x,y,fx=1,fy=1,nbin=NA,normalize=TRUE) {
   # Mix two variables x and y, using weights fx and fy
   # --------------------------------------------------
   # THIS FUNCTION IS TEMPORARILY NOT CODED IN FORTRAN, TO BE DONE
   # perform     : density     f_Z(z) = (fx*f_X(z) + fy*f_Y(z))/(fx+fy)
   #               probability F_Z(z) = (fx*F_X(z) + fy*F_Y(z))/(fx+fy)
   # return value: z

   # init
   f  = fx + fy  # normalize factors
   fx = fx/f
   fy = fy/f

   # aray lengths
   nx   = length(x$x)
   ny   = length(y$x)
   if(is.na(nclass))
     {nz = max(nx,ny,3)}
   else
     {nz = nbin+1}
   # complete pdfs
   tmpAddCum <- function(v) {
      if("Y" %in% names(v)){return(v)}
      nv = length(v$x)
      ii = 1:(nv-1)
      prob = (v$y[ii]+v$y[ii+1])*diff(v$x)
      v$Y = cumsum(c(0,prob))
      return(v)
   }
   x = tmpAddCum(x)
   y = tmpAddCum(y)

   # init z
   zv = range(c(x$x,y$x))
   zv = c(zv[1],mean(zv),zv[2])
   z  = list(val=zv)
   z$den = fx*rvFunValDen(x,zv)$y + fy*rvFunValDen(y,zv)$y
   z$cum = fx*rvFunValCum(x,zv)$Y + fy*rvFunValCum(y,zv)$Y

   # loop
   if(nz > 3) {
      for(i in 4:nz) {
         z  = srvNxtZ(z)
         iz = z$iz
         zv = z$val[iz]
         z$den[iz] = fx*rvFunValDen(x,zv)$y + fy*rvFunValDen(y,zv)$y
         z$cum[iz] = fx*rvFunValCum(x,zv)$Y + fy*rvFunValCum(y,zv)$Y
      }
   }

   pdfz = list(x=z$val,y=z$den,Y=z$cum)
   pdfz = as.plpdf(pdfz,normalize=normalize)

   # return
   return(pdfz)
}


