
!    plPDF is a library to perform calcualtions with piecewise linear random variables.
!
!    Copyright (C) 2016, 2021  Aris Lourens
!
!    This file is part of plPDF.
!
!    plPDF is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    plPDF is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> \file

!> Random Variable Binary FUNction INVerse level 2 \n
!!    Z = 1/X
!!
!! Implementation of the reciprocal (inversion) of a piecewise linear PDF
!!
!! Algorithm: \n
!!  If all values of X>0 or X<0 then \n
!!     all values of Z are defined: min_Z = 1/max(X), max_Z = 1/min(X) \n
!!  else \n
!!     transformed values of Z for X arround 0 tend to infinity \n
!!     start with a very small probability interval around 0 to approximate \n
!!     some limits for Z \n
subroutine rvFunInv2(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-05      initial version   


! declaration section
! ------------------------------------------------------------------------------


implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
!double precision, dimension(nx), intent(in)     :: xcum        !< cumulative probabilities of the CDF of Z


integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
!   3: nz too small
integer,parameter    :: rtnerrcode = 140000000 +  61*1000  ! routine error code for routine rvFunInv2

! local variables
integer          :: k,iz,ix

double precision :: zmin,zmax,x,cum,den,xc0

double precision, parameter :: epsProb=1.e-3   ! minimum probability to calculate Z-limits when x=0 exists


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
call rvFunValCum(xval,xden,nx,0.d0,xc0,1,exitcode) ! probability at x=0

! check nz
if(nz.lt.3) then
   ! ERROR: nz too small, should be at least 3
   exitcode = rtnerrcode+3
   return
endif

! find start values of z
if(xval(1).le.0.d0 .and. xval(nx).ge.0.d0) then
   ! find bin of X containing x=0
   call srvFndBin(xval,nx,0.d0,ix,.false.)
   ! calc probability density and cumulative probability at x=0
   call rvFunValDenCum(xval,xden,nx,0.d0,den,cum,1,exitcode)

   ! approximate limits
   if(cum.lt.epsProb) then
      ! ignore the negative part of X
      ! zmin is easy
      zmin = 1.d0/xval(nx)
      ! find x for prob(x)=prob(x=0)+epsProb
      call rvFunCumVal(xval,xden,nx,x,xc0+epsProb,1,exitcode)
      ! calc zmax
      zmax = 1.d0/x
   else if(cum.gt.(1.d0-epsProb)) then
      ! ignore the positive part of X
      ! find x for prob(x)=prob(x=0)-epsProb
      call rvFunCumVal(xval,xden,nx,x,xc0-epsProb,1,exitcode)
      ! calc zmin
      zmin = 1.d0/x
      ! zmax is easy
      zmax = 1.d0/xval(1)
   else
      ! find x for prob(x)=prob(x=0)-epsProb
      call rvFunCumVal(xval,xden,nx,x,xc0-epsProb,1,exitcode)
      ! calc zmin
      zmin = 1.d0/x
      ! find x for prob(x)=prob(x=0)+epsProb
      call rvFunCumVal(xval,xden,nx,x,xc0+epsProb,1,exitcode)
      ! calc zmax
      zmax = 1.d0/x
   endif
else
   ! both limits are easy
   zmin = 1.d0/xval(nx)   ! mimimum z-value
   zmax = 1.d0/xval(1)    ! maximum z-value
endif
zval(1) = zmin
zval(2) = (zmin+zmax)/2.d0
zval(3) = zmax

call rvFunInv3(xval,xden,nx,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvFunInv3(xval,xden,nx,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvFunInv2

! ******************************************************************************

!> Calc the inverse values of X for Z \n
!!   z=1/x
!!   
subroutine rvFunInv3(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-05      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
!double precision, dimension(nx), intent(in)     :: xcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(in)     :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z


integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  71*1000  ! routine error code for routine rvFunInv3

! local variables
integer          :: k

double precision :: z,x,xc0,cum,den


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
call rvFunValCum(xval,xden,nx,0.d0,xc0,1,exitcode) ! probability at x=0

! loop over all z values
! ----------------------
do k=1,nz
   z = zval(k)
   if(z.gt.0.d0) then
      ! z>0
      ! get values at x
      x = 1.d0/z
      call rvFunValDenCum(xval,xden,nx,x,den,cum,1,exitcode)
      ! convert
      den = den * x**2 ! den/abs((x**(-2)))
      cum = 1.d0-cum + xc0
   else if (z.lt.0.d0) then
      ! z<0
      ! get values at x
      x = 1.d0/z
      call rvFunValDenCum(xval,xden,nx,x,den,cum,1,exitcode)
      ! convert
      den = den * x**2 ! den/abs((x**(-2)))
      cum = xc0-cum
   else
      ! z=0, only from x=inf or -inf gives z=0
      !      we assume, for good reasons, that the probability density
      !      at x=inf and x=-inf is 0. 
      !      so: zden is 0 too, and
      !          zcum(0)=xcum(0)
      den = 0.d0
      cum = xc0
   endif
   ! store result
   zden(k) = den
   zcum(k) = cum
enddo


! end of program
return
end subroutine rvFunInv3


