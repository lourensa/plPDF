
!> \file 

!> Inquire properties or measures of a piecewise linear PDF
!! All routines rvInq* are subroutines and the routines srvInq* are accompanying functions

!> Interpolate the cumulative probability within bin ibin at value x \n
!! In this routine, xcum is used in addition to xval and xden.
!! Since the linearization of a PDF is an approximation, mostly the probability
!! described by xval,xden will not be equal to the difference in xcum. In this routine
!! xcum is used as well to achieve a weighted interpolated cum prob.
subroutine rvInqIntCum(xval,xden,xcum,nx,ibin,x,cum,exitcode)

! History
! programmer------------date------------version---------------------------------
! loure001              2015-08-19      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
double precision, dimension(nx), intent(in)     :: xcum        !< cumulative probabilities of the CDF of X

integer                        , intent(in)     :: ibin        !< bin number to be used for calculation

double precision               , intent(in)     :: x           !< position to interpolate the cumulative probability. \n
                                                               !! x may lay outside bin ibin, but the sense of this
                                                               !! is up to the user of this routine.

double precision               , intent(out)    :: cum         !< cumulative probability at x

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
double precision  :: f,rx,xi,xi1,pxi,pxi1,c1,c2


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! values
xi   = xval(ibin  )
xi1  = xval(ibin+1)
pxi  = xden(ibin)
pxi1 = xden(ibin+1)
rx   = (pxi1-pxi)/(xi1-xi)


! weight factor
f = (x - xi)/(xi1-xi)


! calc
c1  = xcum(ibin)   + (x   - xi)*(pxi + rx*((x+xi )/2.d0-xi) )
c2  = xcum(ibin+1) - (xi1 - x )*(pxi + rx*((x+xi1)/2.d0-xi) )
cum = (1.d0-f)*c1 + f*c2


! end of program
return
end subroutine rvInqIntCum

! ------------------------------------------------------------------------------

!> Function version of rvInqIntCum
function srvInqIntCum(xval,xden,xcum,nx,ibin,x)

! History
! programmer------------date------------version---------------------------------
! loure001              2015-08-19      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none

! function
double precision                                :: srvInqIntCum !< function

! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
double precision, dimension(nx), intent(in)     :: xcum        !< cumulative probabilities of the CDF of X

integer                        , intent(in)     :: ibin        !< bin number to be used for calculation

double precision               , intent(in)     :: x           !< position to interpolate the cumulative probability. \n
                                                               !! x may lay outside bin ibin, but the sense of this
                                                               !! is up to the user of this routine.

! local variables
double precision  :: cum

integer           :: exitcode

! program section
! ------------------------------------------------------------------------------

! calc value
call rvInqIntCum(xval,xden,xcum,nx,ibin,x,cum,exitcode)

! assign function value
srvInqIntCum = cum

! end of program
return
end function srvInqIntCum

! ******************************************************************************

!> Interpolate the probability density within bin ibin at value x \n
!! When the bin width is <= 0, xden(ibin) is returned.
subroutine rvInqIntDen(xval,xden,nx,ibin,x,den,exitcode)

! History
! programmer------------date------------version---------------------------------
! loure001              2015-08-20      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ibin        !< bin number to be used for calculation

double precision               , intent(in)     :: x           !< position to interpolate the cumulative probability. \n
                                                               !! x may lay outside bin ibin, but the sense of this
                                                               !! is up to the user of this routine.

double precision               , intent(out)    :: den         !< probability density at x

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
double precision  :: dx,rx


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! calc
den  = xden(ibin)
dx   = (xval(ibin+1)-xval(ibin))
if (dx.gt.0.d0) then
   rx  = (xden(ibin+1)-xden(ibin))/dx
   den = den + rx*(x-xval(ibin))
endif


! end of program
return
end subroutine rvInqIntDen
! ------------------------------------------------------------------------------
!> Function version of rvInqIntDen
function srvInqIntDen(xval,xden,nx,ibin,x)

! History
! programmer------------date------------version---------------------------------
! loure001              2015-08-20      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none

! function
double precision                                :: srvInqIntDen !< function


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ibin        !< bin number to be used for calculation

double precision               , intent(in)     :: x           !< position to interpolate the cumulative probability. \n
                                                               !! x may lay outside bin ibin, but the sense of this
                                                               !! is up to the user of this routine.


! local variables
double precision  :: den

integer           :: exitcode


! program section
! ------------------------------------------------------------------------------

! calc value
call rvInqIntDen(xval,xden,nx,ibin,x,den,exitcode)

! assign function value
srvInqIntDen = den

! end of program
return
end function srvInqIntDen
! ******************************************************************************



