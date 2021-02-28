
!> \file

!> Random Variable FUNction EXPonent, level 2 \n
!!    Z = exp(X)
!!
!! Implementation of the exponential function
!!
!! Algorithm: \n
!!    prerequisites
!!       nz >= 3
!!    method
!!       z    = exp(x)
!!       zden = xden/exp(x) = xden/z
!!       zcum = xcum
!!     
subroutine rvFunExp2(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

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
integer,parameter    :: rtnerrcode = 140000000 +  62*1000  ! routine error code for routine rvFunExp2

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

! find start values of z
zmin = exp(xval(1))    ! mimimum z-value
zmax = exp(xval(nx))   ! maximum z-value
zval(1) = zmin
zval(2) = (zmin+zmax)/2.d0
zval(3) = zmax

call rvFunExp3(xval,xden,nx,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvFunExp3(xval,xden,nx,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvFunExp2

! ******************************************************************************

!> Calc the exponential of X for Z \n
!!   z=exp(x)
!!   
subroutine rvFunExp3(xval,xden,nx,zval,zden,zcum,nz,exitcode) 

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
integer,parameter    :: rtnerrcode = 140000000 +  72*1000  ! routine error code for routine rvFunExp3

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
   x = log(z)
   call rvFunValDenCum(xval,xden,nx,x,den,cum,1,exitcode)
   ! convert
   den = den / z
   ! store result
   zden(k) = den
   zcum(k) = cum
enddo


! end of program
return
end subroutine rvFunExp3

! ******************************************************************************


