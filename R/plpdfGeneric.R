
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
#' Method for binary operations on PL-PDFs
#' @param pdf1  object of class \code{plpdf} or scalar variable
#' @param pdf2  object of class \code{plpdf} or scalar variable
#' @export
Ops.plpdf <- function(pdf1,pdf2) {

   # calc
   if (.Generic %in% c("+","-","*","/")) {
      nbin = NA
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


