
#    plPDF is a library to perform calcualtions with piecewise linear random variables.
#
#    Copyright (C) 2016, 2021, 2023  Aris Lourens
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

# statistical functions without a generic equivalent
# ==================================================


#' Variance of a pdf
#'
#' Calculate the variance of a pdf of class `plpdf`. Unfortunately, the R-function `var` is not generic.
#' 
#' @param pdf    object of type `plpdf`.
#' @seealso [`plStd`]
#' @export
plVar <- function(pdf) {
   ret = plpdfStat(pdf)
   return(ret[2]^2)
}


#' Standard deviation of a pdf
#'
#' Calculate the Standard deviation of a pdf of class `plpdf`.
#' 
#' @param pdf    object of type `plpdf`.
#' @seealso [`plVar`]
#' @export
plStd <- function(pdf) {
   ret = plpdfStat(pdf)
   return(ret[2])
}


