
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

# set options for the plPDF funtions
# ==================================

#' Set the Default Number of Bins
#'
#' Set the default number of bins for PL-PDFS. This value is stored using the R \code{options} function and can be set through that function too. This function is supplied for convenience. The name of the option is \code{plpdfNbin}.
#' @param   nbin    integer; default number if bins for an output PL-PDF
#' @return  The function invisibly returns the number of bins set.
#' @seealso Use \code{\link{plpdfGetNbin}} to retrieve the option value, and \code{\link{plpdfUnsetNbin}} to remove it.
#' @export
plpdfSetNbin <- function(nbin=51L) {

   # init
   optname = "plpdfNbin"
   nbindef = 50L            # default  number of bins
   nbinmin = 3L

   # check nbin
   if (is.null(nbin)) {
      nbin = nbindef
   } else if (! is.finite(nbin)) {
      nbin = nbindef
   }
   if (nbin < nbinmin) {
      warning("Default number of bins too small, set to",nbinmin)
      nbin = nbinmin
   }
   nbin = nbin[1]
   nbin = as.integer(ceiling(nbin))

   # set option
   opts = list()
   opts[[optname]] = nbin
   options(opts)

   # return
   return(invisible(nbin))
}

#' Remove the default number of bins from the options list
#'
#' Remove the default number of bins from the options list. This value is stored using the R \code{options} function and can be set through that function too. This function is supplied for convenience. The name of the option is \code{plpdfNbin}.
#' @return  \code{NULL}
#' @seealso Use \code{\link{plpdfSetNbin}} to set option value, and \code{\link{plpdfGetNbin}} to retrieve it.
#' @export
plpdfUnsetNbin <- function() {

   # init
#   optname = "plpdfNbin"

   # unset option
   options(plpdfNbin=NULL)

   # return
   return(invisible(NULL))
}


#' Get the Default Number of Bins
#'
#' Get the default number of bins for PL-PDFS. This value is stored using the R \code{options} function and can be accessed through that function too. This function is supplied for convenience. The name of the option is \code{plpdfNbin}. If the option is not defined it will be set.
#' @param   default  Default value if the option is not set. If \code{default=NULL} and the option is currently not set then \code{\link{plpdfSetNbin}} is called and the option is set.
#' @return  This function returns the number of bins set.
#' @seealso Use \code{\link{plpdfSetNbin}} to set the option value, and \code{\link{plpdfUnsetNbin}} to remove it.
#' @export
plpdfGetNbin <- function(default=NULL) {

   # init
   optname = "plpdfNbin"

   # get nbin
   nbin = getOption(x=optname,default=default)
   if (is.null(nbin)) {
      nbin = plpdfSetNbin(NULL)
   }

   # return
   return(nbin)
}

