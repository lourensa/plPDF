
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

!> Random Variable FUNction SQuare, level 2 \n
!!    Z = X^2
!!
!! Implementation of the logarithmic function
!!
!! Algorithm: \n
!!    prerequisites
!!       nz >= 3
!!    method
!!       z    = x^2
!!       zden = (xden(-x)+xden(x))/(|2x|)
!!       zcum = xcum(x)-xcum(-x)
!!
!!       At z=0 then zden tends to infinity when xden(x=0)>0. To avoid this
!!       the minimum value of z is set to a value >0 with probability epsProb.
!!       If nz>3 then z=0 will afterwards be added with an approximated density
!!       and zcum(z=0)=0.
!!     
subroutine rvFunSq2(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! Aris                  2017-05-17      initial version   


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
integer,parameter    :: rtnerrcode = 140000000 +  65*1000  ! routine error code for routine rvFunSq2


! local variables
integer          :: k,iz,izshift

double precision :: zmin,zmax,x,xc0

double precision, parameter :: epsProb=1.e-3   ! minimum probability to calculate Z-limits when x=0 exists

logical          :: lz0

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! check nz
if(nz.lt.3) then
   ! ERROR: nz too small, should be at least 3
   exitcode = rtnerrcode+3
   return
endif


! find start values of z
lz0 = .false.  ! .true. : calc density of z=0 afterwards
               ! .false.: don't calc density at z=0
izshift = 0
if (xval(1).le.0.d0) then          ! mimimum z-value
   if (xval(nx).lt.0.d0) then
      zmin = xval(nx)**2
   else
      ! find x for prob(x)=prob(x=0)+epsProb
      call rvFunValCum(xval,xden,nx,0.d0,xc0,1,exitcode) ! probability at x=0
      call rvFunCumVal(xval,xden,nx,x,xc0+0.5*epsProb,1,exitcode)
      zmin = x**2
      if (nz.gt.3) then
         ! approximate the density of z=0 separately
         lz0     = .true.
         izshift = 1
      endif
   endif
else
   zmin = xval(1)**2
endif
zmax = max(xval(1)**2,xval(nx)**2)   ! maximum z-value

zval(1+izshift) = zmin
zval(2+izshift) = (zmin+zmax)/2.d0
zval(3+izshift) = zmax

call rvFunSq3(xval,xden,nx,zval(1+izshift),zden(1+izshift),zcum(1+izshift),3,exitcode)
if(exitcode.ne.0) return

! approximate density at z=0
if (lz0) then
   zval(1) = 0.d0
   zcum(1) = 0.d0
   zden(1) = 2.d0*zcum(2)/zval(2) - zden(2)
endif

! loop over all z values
! ----------------------
do k=4+izshift,nz

   ! next z value
   call srvNxtZ(zval(1+izshift),zden(1+izshift),zcum(1+izshift),k-izshift,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   iz = iz + izshift
   call rvFunSq3(xval,xden,nx,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvFunSq2

! ******************************************************************************

!> Calc the square of X for Z \n
!!   z=x^2
!!   
subroutine rvFunSq3(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! Aris                  2017-05-17      initial version   


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
integer,parameter    :: rtnerrcode = 140000000 +  75*1000  ! routine error code for routine rvFunSq3

! local variables
integer          :: k

double precision :: z,x,cum,den
double precision, dimension(2) :: x2,cum2,den2

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all z values
! ----------------------
do k=1,nz
   z = zval(k)
   x = sqrt(z)
   x2(1)=-x
   x2(2)= x
   call rvFunValDenCum(xval,xden,nx,x2,den2,cum2,2,exitcode)
   ! convert
   den = sum(den2) /(2.d0*x)
   cum = cum2(2)-cum2(1)
   ! store result
   zden(k) = den
   zcum(k) = cum
enddo


! end of program
return
end subroutine rvFunSq3

! ******************************************************************************


