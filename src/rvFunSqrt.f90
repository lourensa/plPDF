
!> \file

!> Random Variable FUNction SQuare RooT, level 2 \n
!!    Z = sqrt(X)
!!
!! Implementation of the square root function
!!
!! Algorithm: \n
!!    prerequisites
!!       nz >= 3
!!       x >= 0
!!    method
!!       z    = sqrt(x)
!!       zden = xden/((1/2)*x^(-1/2)) = xden*2*z
!!       zcum = xcum
!!     
subroutine rvFunSqrt2(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

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
!   4: x contains negative values
integer,parameter    :: rtnerrcode = 140000000 +  66*1000  ! routine error code for routine rvFunSqrt2


! local variables
integer          :: k,iz

double precision :: zmin,zmax

!double precision, parameter :: epsProb=1.e-3   ! minimum probability to calculate Z-limits when x=0 exists

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

! check x
if(xval(1).lt.0.d0) then
   ! ERROR: x should be greater than 0 for a valid logarithm
   exitcode = rtnerrcode+4
   return
endif

! find start values of z
zmin = sqrt(xval(1))    ! mimimum z-value
zmax = sqrt(xval(nx))   ! maximum z-value
zval(1) = zmin
zval(2) = (zmin+zmax)/2.d0
zval(3) = zmax

call rvFunSqrt3(xval,xden,nx,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvFunSqrt3(xval,xden,nx,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvFunSqrt2

! ******************************************************************************

!> Calc the natural logarithm of X for Z \n
!!   z=ln(x)
!!   
subroutine rvFunSqrt3(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

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
integer,parameter    :: rtnerrcode = 140000000 +  76*1000  ! routine error code for routine rvFunSqrt3

! local variables
integer          :: k

double precision :: z,x,cum,den

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all z values
! ----------------------
do k=1,nz
   z = zval(k)
   x = z**2
   call rvFunValDenCum(xval,xden,nx,x,den,cum,1,exitcode)
   ! convert
   den = den * 2.d0 * z
   ! store result
   zden(k) = den
   zcum(k) = cum
enddo


! end of program
return
end subroutine rvFunSqrt3

! ******************************************************************************


