
!> \file

!> Random Variable Binary OPeration SUMmation level 2 \n
!!    Z = X + Y
!!
!! Implementation of the summation of two independent piecewise linear PDFs
!!
subroutine rvBopSum2(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-27      initial version   


! declaration section
! ------------------------------------------------------------------------------


implicit none

! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
!   3: nz too small
integer,parameter    :: rtnerrcode = 140000000 +  21*1000  ! routine error code for routine rvBopSum2


! local variables
integer          :: k,iz

double precision :: zmin,zmax


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

! start values of z
zmin     = xval(1) +yval(1)   ! mimimum z-value for all joint bins
zmax     = xval(nx)+yval(ny)  ! maximum z-value for all joint bins
zval(1)  = zmin
zval(2)  = (zmin+zmax)/2.
zval(3)  = zmax

call rvBopSum3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvBopSum3(xval,xden,nx,yval,yden,ny,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvBopSum2

! ******************************************************************************

!> Random Variable Binary OPeration SUBtraction level 2 \n
!!    Z = X - Y
!!
!! Implementation of the subtraction of two independent piecewise linear PDFs
!!
subroutine rvBopSub2(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-27      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
!   3: nz too small
integer,parameter    :: rtnerrcode = 140000000 +  22*1000  ! routine error code for routine rvBopSub2


! local variables
integer          :: k,iz

double precision :: zmin,zmax


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

! start values of z
zmin     = xval(1) -yval(ny)  ! mimimum z-value for all joint bins
zmax     = xval(nx)-yval(1)   ! maximum z-value for all joint bins
zval(1)  = zmin
zval(2)  = (zmin+zmax)/2.
zval(3)  = zmax

call rvBopSub3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvBopSub3(xval,xden,nx,yval,yden,ny,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvBopSub2

! ******************************************************************************

!> Random Variable Binary OPeration MULtiplication level 2 \n
!!    Z = X * Y
!!
!! Implementation of the multiplication of two independent piecewise linear PDFs
!!
subroutine rvBopMul2(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-27      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
!   3: nz too small
integer,parameter    :: rtnerrcode = 140000000 +  23*1000  ! routine error code for routine rvBopMul2


! local variables
integer          :: k,iz

double precision :: zmin,zmax,z1,z2,z3,z4,xmean,ymean


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

! extreme z-values available
z1 = xval(1) * yval(1)
z2 = xval(1) * yval(ny)
z3 = xval(nx)* yval(1)
z4 = xval(nx)* yval(ny)

! calc mean values
call rvStatMean(xval,xden,nx,xmean,exitcode)
call rvStatMean(yval,yden,ny,ymean,exitcode)

! start values of z
zmin     = min(z1,z2,z3,z4) ! mimimum z-value for all joint bins
zmax     = max(z1,z2,z3,z4) ! maximum z-value for all joint bins
zval(1)  = zmin
zval(2)  = xmean*ymean !(zmin+zmax)/2.
zval(3)  = zmax

call rvBopMul3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,3,exitcode)
if(exitcode.ne.0) return

! loop over all z values
! ----------------------
do k=4,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvBopMul3(xval,xden,nx,yval,yden,ny,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvBopMul2

! ******************************************************************************

!> Random Variable Binary OPeration DIVision level 2 \n
!!    Z = X / Y
!!
!! Implementation of the division of two independent piecewise linear PDFs
!!
subroutine rvBopDiv2(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-27      vcrc    initial version   
! Aris                  2016-12-13      zmin,zmax calculation improved. Iterative loop added.


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! routine error code
!   3: nz too small
integer,parameter    :: rtnerrcode = 140000000 +  24*1000  ! routine error code for routine rvBopDiv2

! local variables
integer          :: k,iz,n !,iter


double precision, parameter:: minprob=0.999  ! minimum amount of probability of Z within initial z-values
                                             ! This parameter is used when 0 is in the domain of Y

!double precision :: tcum,tcumx,tcumy,tcumz,cumlow,cumhig,dval
double precision :: zlim(2)


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


! get z-value starting limits
call rvDivZLimits(xval,xden,nx,yval,yden,ny,minprob,zlim,exitcode)
n = 3
zval(1) = zlim(1)
zval(2) = 0.5d0*(zlim(1)+zlim(2))
zval(3) = zlim(2)

! calculate initial results for Z
call rvBopDiv3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,n,exitcode)

! loop over all z values
! ----------------------
do k=n+1,nz

   ! next z value
   call srvNxtZ(zval,zden,zcum,k,iz,exitcode)
   if(exitcode.ne.0) return

   ! calc probs for zval(iz)
   call rvBopDiv3(xval,xden,nx,yval,yden,ny,zval(iz),zden(iz),zcum(iz),1,exitcode)
   if(exitcode.ne.0) return

enddo


! end of program
return
end subroutine rvBopDiv2

! ******************************************************************************

!> Random Variable Binary OPeration SUMmation level 3 \n
!!    Z = X + Y
!!
!! Implementation of the summation of two independent piecewise linear PDFs
!!
subroutine rvBopSum3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-16      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(in)     :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  31*1000  ! routine error code for routine rvBopSum3

! local variables
integer          :: i,j,k,iend,istp,ilst,jend,jstp,jlst,ii,jj

double precision :: z,zd,zc,zmin,zmax
double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
zmin     = xval(1) +yval(1)   ! mimimum z-value for all joint bins
zmax     = xval(nx)+yval(ny)  ! maximum z-value for all joint bins

! loop over all z values
! ----------------------
do k=1,nz

   ! init
   z  = zval(k)
   zd = 0.0d0
   zc = 0.0d0

   ! let z be within defined bins
   z  = min(max(zmin,z),zmax)

   ! intersection of function with axis
   xu = z - yval(1)
   xu = max(xval(1),min(xu,xval(nx)))
   yl = z - xu
   yl = max(yval(1),min(yl,yval(ny)))

   ! find first joint bin
   call srvFndBin(xval,nx,xu,i,.true.)
   call srvFndBin(yval,ny,yl,j,.false.)

   ! set counters
   iend = 1
   istp = -1
   ilst = i-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
   jend = ny-1
   jstp = 1
   jlst = 0       ! last bin of y used to calculate yc

   ! loop over all intersectd joint bins
   cont = .true.
   yc   = 0.0d0
   do while (cont)

      ! get function parameters
      xi   = xval(i)
      xi1  = xval(i+1)
      pxi  = xden(i)
      pxi1 = xden(i+1)
      rx   = (pxi1-pxi)/(xi1-xi)
      yj   = yval(j)
      yj1  = yval(j+1)
      pyj  = yden(j)
      pyj1 = yden(j+1)
      ry   = (pyj1-pyj)/(yj1-yj)

      ! calc xl,yu
      xl   = z - yj1
      nxtx = (xl.lt.xi)
      xl   = max(xi,min(xl,xi1))
      yu   = z - xl
      yu   = max(yj,min(yu,yj1))

      dy   = yu - yl
      pxl  = pxi + rx*(xl-xi)
      pxu  = pxi + rx*(xu-xi)
      pyl  = pyj + ry*(yl-yj)
      pyu  = pyj + ry*(yu-yj)
      

      ! calc cumulative probability
      ! ---------------------------
      zc = zc + (1./2.)*pxl*pyu*dy**2                         & 
              - (1./3.)*pxl*ry *dy**3                         &
              + (1./6.)*rx *pyu*dy**3                         &
              - (1./8.)*rx *ry *dy**4                         

      zc = zc + (1./4.)*(pyu+pyl)*(yu-yl)*(pxl+pxi)*(xl-xi)   &        !  |
              + (1./4.)*(pyl+pyj)*(yl-yj)*(pxu+pxi)*(xu-xi)            !  --

      ! calc cumulative probability for y until bin j-1
      do jj=jlst+jstp,j-jstp,jstp
         yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
      enddo
      jlst = j-jstp

      ! add probability for all joint bins (i,<j)
      xc = 0.0d0
      do ii=ilst+istp,i,istp
         xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
      enddo
      zc = zc + xc*yc
      ilst = i

      ! calc probability density
      ! ------------------------
      ! keep xl and yl fixed and let xu and yu be a function of z
      ! First derivative (of area d) with respect to z is:
!      zd = zd + (1./2.)*pxl*ry *  dy**2                       &
!              + (1./2.)*pxl*pyu*2*dy                          &
!              - (1./3.)*pxl*ry *3*dy**2                       &
!              + (1./6.)*rx *ry *  dy**3                       &
!              + (1./6.)*rx *pyu*3*dy**2                       &
!              - (1./8.)*rx *ry *4*dy**3
      zd = zd +         pxl*pyu*dy                            &
              - (1./2.)*pxl*ry *dy**2                         &
              + (1./2.)*rx *pyu*dy**2                         &
              - (1./3.)*rx *ry *dy**3


      ! next joint bin
      ! --------------
      if (nxtx) then
         i = i + istp
         if (istp*(i-iend).gt.0) cont=.false.
      else
         j = j + jstp
         if (jstp*(j-jend).gt.0) cont=.false.
      endif

      xu = xl
      yl = yu
   enddo


   ! complete yc
   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,jend,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = jend

   ! add forgotten cumulative probability
   ! When the line z = x + y intersects the line y=yval(ny) at some bin i>1 then
   ! not all probbaility is added to zc. 
   xc = 0.0d0
   do ii=ilst+istp,iend,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                    ! the calculated value yc is used

   ! store results
!   zden(k) = zd
   zden(k) = max(0.d0,zd)   ! due to rounding errors, zd is sometimes slightly negative
   zcum(k) = zc
enddo


! end of program
return
end subroutine rvBopSum3

! ******************************************************************************

!> Random Variable Binary OPeration SUBtraction level 3 \n
!!    Z = X - Y
!!
!! Implementation of the subtraction of two independent piecewise linear PDFs
!!
subroutine rvBopSub3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-21      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(in)     :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  32*1000  ! routine error code for routine rvBopSub3

! local variables
integer          :: i,j,k,iend,istp,ilst,jend,jstp,jlst,ii,jj

double precision :: z,zd,zc,zmin,zmax
double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
zmin     = xval(1) -yval(ny)  ! mimimum z-value for all joint bins
zmax     = xval(nx)-yval(1)   ! maximum z-value for all joint bins

! loop over all z values
! ----------------------
do k=1,nz

   ! init
   z  = zval(k)
   zd = 0.0d0
   zc = 0.0d0

   ! let z be within defined bins
   z  = min(max(zmin,z),zmax)

   ! intersection of function with axis
   xu = z + yval(ny)
   xu = max(xval(1),min(xu,xval(nx)))
   yu = xu - z
   yu = max(yval(1),min(yu,yval(ny)))

   ! find first joint bin
   call srvFndBin(xval,nx,xu,i,.true.)
   call srvFndBin(yval,ny,yu,j,.true.)

   ! set counters
   iend = 1
   istp = -1
   ilst = i-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
   jend = 1
   jstp = -1
   jlst = ny      ! last bin of y used to calculate yc

   ! loop over all intersectd joint bins
   cont = .true.
   yc   = 0.0d0
   do while (cont)

      ! get function parameters
      xi   = xval(i)
      xi1  = xval(i+1)
      pxi  = xden(i)
      pxi1 = xden(i+1)
      rx   = (pxi1-pxi)/(xi1-xi)
      yj   = yval(j)
      yj1  = yval(j+1)
      pyj  = yden(j)
      pyj1 = yden(j+1)
      ry   = (pyj1-pyj)/(yj1-yj)

      ! calc xl,yl
      xl   = z + yj
      nxtx = (xl.lt.xi)
      xl   = max(xi,min(xl,xi1))
      yl   = xl - z
      yl   = max(yj,min(yl,yj1))

      dy   = yu - yl
      pxl  = pxi + rx*(xl-xi)
      pxu  = pxi + rx*(xu-xi)
      pyl  = pyj + ry*(yl-yj)
      pyu  = pyj + ry*(yu-yj)
      

      ! calc cumulative probability
      ! ---------------------------
      zc = zc + (1./2.)*pxl*pyl*dy**2                         & 
              + (1./3.)*pxl*ry *dy**3                         &
              + (1./6.)*rx *pyl*dy**3                         &
              + (1./8.)*rx *ry *dy**4                         

      zc = zc + (1./4.)*(pyj1+pyu)*(yj1-yu)*(pxu+pxi)*(xu-xi) &        !  --
              + (1./4.)*(pyu +pyl)*(yu -yl)*(pxl+pxi)*(xl-xi)          !  |

      ! calc cumulative probability for y until bin j-1
      do jj=jlst+jstp,j-jstp,jstp
         yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
      enddo
      jlst = j-jstp

      ! add probability for all joint bins (i,<j)
      xc = 0.0d0
      do ii=ilst+istp,i,istp
         xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
      enddo
      zc = zc + xc*yc
      ilst = i

      ! calc probability density
      ! ------------------------
      ! keep xl and yu fixed and let xu and yl be a function of z
      ! First derivative (of area d) with respect to z is:
!      zd = zd - (1./2.)*pxl*ry *  dy**2                       &
!              + (1./2.)*pxl*pyl*2*dy                          &
!              + (1./3.)*pxl*ry *3*dy**2                       &
!              - (1./6.)*rx *ry *  dy**3                       &
!              + (1./6.)*rx *pyl*3*dy**2                       &
!              + (1./8.)*rx *ry *4*dy**3
      zd = zd + (1./2.)*pxl*ry   *dy**2                       &
              +         pxl*pyl  *dy                          &
              + (1./2.)*rx *pyl  *dy**2                       &
              + (1./3.)*rx *ry   *dy**3


      ! next joint bin
      ! --------------
      if (nxtx) then
         i = i + istp
         if (istp*(i-iend).gt.0) cont=.false.
      else
         j = j + jstp
         if (jstp*(j-jend).gt.0) cont=.false.
      endif

      xu = xl
      yu = yl
   enddo


   ! complete yc
   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,jend,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = jend

   ! add forgotten cumulative probability
   ! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
   ! not all probbaility is added to zc. 
   xc = 0.0d0
   do ii=ilst+istp,iend,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                    ! the calculated value yc is used

   ! store results
!   zden(k) = zd
   zden(k) = max(0.d0,zd)   ! due to rounding errors, zd is sometimes slightly negative
   zcum(k) = zc
enddo


! end of program
return
end subroutine rvBopSub3

! ******************************************************************************

!> Random Variable Binary OPeration MULtiplication level 3 \n
!!    Z = X * Y
!!
!! Implementation of the multiplication of two independent piecewise linear PDFs
!!
subroutine rvBopMul3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-21      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(in)     :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  33*1000  ! routine error code for routine rvBopMul3

! local variables
integer          :: k

double precision :: z,zd,zc,zmin,zmax,z1,z2,z3,z4

logical          :: q1,q2,q3,q4

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! extreme z-values available
z1 = xval(1) * yval(1)
z2 = xval(1) * yval(ny)
z3 = xval(nx)* yval(1)
z4 = xval(nx)* yval(ny)

zmin     = min(z1,z2,z3,z4) ! mimimum z-value for all joint bins
zmax     = max(z1,z2,z3,z4) ! maximum z-value for all joint bins

! which quadrants are available
q1 = .false.
q2 = .false.
q3 = .false.
q4 = .false.
if (xval(1) .lt. 0.d0) then
   if (yval(1) .lt. 0.d0) q3 = .true.
   if (yval(ny).gt. 0.d0) q2 = .true.
endif
if (xval(nx).gt. 0.d0) then
   if (yval(1) .lt. 0.d0) q4 = .true.
   if (yval(ny).gt. 0.d0) q1 = .true.
endif


! loop over all z values
! ----------------------
do k=1,nz

   ! init
   z  = zval(k)
   zd = 0.0d0
   zc = 0.0d0

   ! let z be within defined bins
   z  = min(max(zmin,z),zmax)

   ! workaround for z.eq.0
   ! The density at z=0 may be infinite, so a z-value close to 0 will be used.
   if (z.eq.0.d0 .and. zmin.lt.0.d0 .and. zmax.gt.0.d0) then
      if (-zmin.gt.zmax) then
         z=-(zmax-zmin)*1.d-100
      else
         z= (zmax-zmin)*1.d-100
      endif
   endif

   ! calculate all quadrants
   if (z.gt.0.d0) then
      if (q1) call srvBopMulQ1(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      if (q3) call srvBopMulQ3(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
   else if (z.lt.0.d0) then
      if (q2) call srvBopMulQ2(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      if (q4) call srvBopMulQ4(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
   else
    !  call srvBopMulQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
   endif


   ! store results
!   zden(k) = zd
   zden(k) = max(0.d0,zd)   ! due to rounding errors, zd is sometimes slightly negative
   zcum(k) = zc

   if(exitcode.ne.0) return
enddo


! end of program
return
end subroutine rvBopMul3

! ******************************************************************************

!> Random Variable Binary OPeration DIVision level 3 \n
!!    Z = X / Y
!!
!! Implementation of the division of two independent piecewise linear PDFs
!!
subroutine rvBopDiv3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,nz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-24      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

integer                        , intent(in)     :: nz          !< number of elements of variable Z
double precision, dimension(nz), intent(in)     :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  34*1000  ! routine error code for routine rvBopDiv3

! local variables
integer          :: k

double precision :: z,zd,zc !,zmin,zmax

logical          :: q1,q2,q3,q4

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
!zmin     = min(xval(1) *yval(ny),xval(nx)*yval(1) )   ! mimimum z-value for all joint bins
!zmax     = max(xval(1) *yval(1) ,xval(nx)*yval(ny))   ! maximum z-value for all joint bins

! which quadrants are available
q1 = .false.
q2 = .false.
q3 = .false.
q4 = .false.
if (xval(1) .lt. 0.d0) then
   if (yval(1) .lt. 0.d0) q3 = .true.
   if (yval(ny).gt. 0.d0) q2 = .true.
endif
if (xval(nx).gt. 0.d0) then
   if (yval(1) .lt. 0.d0) q4 = .true.
   if (yval(ny).gt. 0.d0) q1 = .true.
endif

   
! loop over all z values
! ----------------------
do k=1,nz

   ! init
   z  = zval(k)
   zd = 0.0d0
   zc = 0.0d0

   ! let z be within defined bins
!   z  = min(max(zmin,z),zmax)

   ! calculate all quadrants
   if (z.gt.0.d0) then
      if (q1) call srvBopDivQ1(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      if (q3) call srvBopDivQ3(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      if (xval(1).ge.0.d0 .or. xval(nx).le.0.d0) then
         call srvBopDivQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      endif
   else if (z.lt.0.d0) then
      if (q2) call srvBopDivQ2(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
      if (q4) call srvBopDivQ4(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
   else
      call srvBopDivQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode)
   endif


   ! store results
!   zden(k) = zd
   zden(k) = max(0.d0,zd)   ! due to rounding errors, zd is sometimes slightly negative
   zcum(k) = zc
enddo


! end of program
return
end subroutine rvBopDiv3

! ******************************************************************************

!> Calculate the probability for each quadrant and add to zc \n
!! Only the probabilities of the selected quadrants are added
!!
subroutine srvBopPrbQ(xval,xden,nx,yval,yden,ny,zc,q1,q2,q3,q4,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(inout)  :: zc          !< summation for cumulative probability

logical                        , intent(in)     :: q1,q2,q3,q4 !< quadrants of which the probability has 
                                                               !! to be added to zc

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: i,j

double precision :: prxlow,prxupp,prylow,pryupp,vall,denl,valu,denu


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! 
! ----------------------

! init marginal probabilities
prxlow = 0.d0    ! probability for x<0
prxupp = 0.d0    ! probability for x>0
prylow = 0.d0    ! probability for y<0
pryupp = 0.d0    ! probability for y>0


! calculate marginal probabilities

! prxlow: q2,q3
if (q2 .or. q3) then
   i = 1
   do while(xval(i).lt.0.d0 .and. i.lt.nx)
   !   prxlow = prxlow + 0.5*(xden(i+1)+xden(i))*(min(xval(i+1),0.d0)-xval(i))
      if (xval(i+1) .le. 0.d0) then
         valu = xval(i+1)
         denu = xden(i+1)
      else
         valu = 0.d0    ! = min(xval(i+1),0.d0)
         denu = xden(i) + (valu-xval(i)) * (xden(i+1)-xden(i))/(xval(i+1)-xval(i))
      endif
      prxlow = prxlow + 0.5d0*(denu+xden(i))*(valu-xval(i))
      i = i + 1
   end do
endif

! prxupp: q1,q4
if (q1 .or. q4) then
   i = nx
   do while(xval(i).gt.0.d0 .and. i.gt.1)
      i = i - 1
   !   prxupp = prxupp + 0.5*(xden(i+1)+xden(i))*(xval(i+1)-max(xval(i),0.d0))
      if (xval(i) .ge. 0.d0) then
         vall = xval(i)
         denl = xden(i)
      else
         vall = 0.d0    ! = min(xval(i+1),0.d0)
         denl = xden(i) + (xval(i+1)-vall) * (xden(i+1)-xden(i))/(xval(i+1)-xval(i))
      endif
      prxupp = prxupp + 0.5d0*(xden(i+1)+denl)*(xval(i+1)-vall)
   end do
endif

! prylow: q3,q4
if (q3 .or. q4) then
   j = 1
   do while(yval(j).lt.0.d0 .and. j.lt.ny)
   !   prylow = prylow + 0.5*(yden(j+1)+yden(j))*(min(yval(j+1),0.d0)-yval(j))
      if (yval(j+1) .le. 0.d0) then
         valu = yval(j+1)
         denu = yden(j+1)
      else
         valu = 0.d0    ! = min(yval(j+1),0.d0)
         denu = yden(j) + (valu-yval(j)) * (yden(j+1)-yden(j))/(yval(j+1)-yval(j))
      endif
      prylow = prylow + 0.5d0*(denu+yden(j))*(valu-yval(j))
      j = j + 1
   end do
endif

! pryupp: q1,q2
if (q1 .or. q2) then
   j = ny
   do while(yval(j).gt.0.d0 .and. j.gt.1)
      j = j - 1
   !   pryupp = pryupp + 0.5*(yden(j+1)+yden(j))*(yval(j+1)-max(yval(j),0.d0))
      if (yval(j) .ge. 0.d0) then
         vall = yval(j)
         denl = yden(j)
      else
         vall = 0.d0    ! = min(yval(j+1),0.d0)
         denl = yden(j) + (yval(j+1)-vall) * (yden(j+1)-yden(j))/(yval(j+1)-yval(j))
      endif
      pryupp = pryupp + 0.5d0*(yden(j+1)+denl)*(yval(j+1)-vall)
   end do
endif


! add probabilities of the selected quadrants
if(q1) zc = zc + prxupp*pryupp
if(q2) zc = zc + prxlow*pryupp
if(q3) zc = zc + prxlow*prylow
if(q4) zc = zc + prxupp*prylow


! end of program
return
end subroutine srvBopPrbQ


! ******************************************************************************

!> Perform multiplcaion for z=0 for rvBopMul3
!!
subroutine srvBopMulQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !> exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 330*1000  ! routine error code for routine srvBopMulQ0

! local variables
integer          :: i,j,i0,j0

double precision :: p0x,rx,p0x0,rx0,xl,xu,xeps,p0y,ry,p0y0,ry0,yeps,yl,yu

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! Cumulative probaility
! ---------------------

! add probabilities of quadrant 2 and quadrant 4
call srvBopPrbQ(xval,xden,nx,yval,yden,ny,zc,.false.,.true.,.false.,.true.,exitcode)



! Probability density
! -------------------

! find x-bin containing 0
if (xval(1).gt.0.d0 .or. xval(nx).lt.0.d0) then
   ! no bin containing 0
   p0x0 = 0.d0
   rx0  = 0.d0
   xeps = 0.d0
   if (xval(1).gt.0.d0) then
      i0 = 1
   else
      i0 = nx - 1
   endif
else
   call srvFndBin(xval,nx,0.d0,i,.true.)
   i0   = i
   rx0  = (xden(i+1)-xden(i))/(xval(i+1)-xval(i))
   p0x0 = xden(i) + rx0*(0.d0-xval(i))

   ! minimum value to allow x close to zero
   xeps = xval(i+1)-xval(i)
   xeps = xeps*10.**(-13)
endif

! find y-bin containing 0
if (yval(1).gt.0.d0 .or. yval(ny).lt.0.d0) then
   ! no bin containing 0
   p0y0 = 0.d0
   ry0  = 0.d0
   yeps = 0.d0
   if (yval(1).gt.0.d0) then
      j0 = 1
   else
      j0 = ny - 1
   endif
else
   call srvFndBin(yval,ny,0.d0,j,.true.)
   j0   = j
   ry0  = (yden(j+1)-yden(j))/(yval(j+1)-yval(j))
   p0y0 = yden(j) + ry0*(0.d0-yval(j))

   ! minimum value to allow y close to zero
   yeps = yval(j+1)-yval(j)
   yeps = yeps*10.**(-13)
endif


do i=1,i0!-1
   ! derived from Q3
   rx  = (xden(i+1)-xden(i))/(xval(i+1)-xval(i))
   p0x = xden(i) + rx0*(0.d0-xval(i))
   xl  = xval(i)
   xu  = xval(i+1)

   xl  = min(xl,-xeps)
   xu  = min(xu,-xeps)

   zd  = zd - p0x*p0y0*log(xu/xl)    &
            - rx *p0y0*(xu-xl)
end do

do j=1,j0!-1
   ! derived from Q3
   ry  = (yden(j+1)-yden(j))/(yval(j+1)-yval(j))
   p0y = yden(j) + ry0*(0.d0-yval(j))
   yl  = yval(j)
   yu  = yval(j+1)

   yl  = min(yl,-yeps)
   yu  = min(yu,-yeps)

   zd  = zd - p0x0*p0y*log(yu/yl)    &
            - p0x0*ry *(yu-yl)
end do


do i=i0+1*0,nx-1
   ! derived from Q1
   rx  = (xden(i+1)-xden(i))/(xval(i+1)-xval(i))
   p0x = xden(i) + rx0*(0.d0-xval(i))
   xl  = xval(i)
   xu  = xval(i+1)

   xl  = max(xl,+xeps)
   xu  = max(xu,+xeps)

   zd  = zd + p0x*p0y0*log(xu/xl)    &
            + rx *p0y0*(xu-xl)
end do

do j=j0+1*0,ny-1
   ! derived from Q1
   ry  = (yden(j+1)-yden(j))/(yval(j+1)-yval(j))
   p0y = yden(j) + ry0*(0.d0-yval(j))
   yl  = yval(j)
   yu  = yval(j+1)

   yl  = max(yl,+yeps)
   yu  = max(yu,+yeps)

   zd  = zd + p0x0*p0y*log(yu/yl)    &
            + p0x0*ry *(yu-yl)
end do




! end of program
return
end subroutine srvBopMulQ0

! ******************************************************************************

!> Perform multiplicaion for quadrant 1 for rvBopMul3
!!
!! If xmin<0 then the probability of quadrant 2 is calculated as well in this routine \n
!! else the probability of quadrant 4 is also calculated
!!
subroutine srvBopMulQ1(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 331*1000  ! routine error code for routine srvBopMulQ1

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,lny

logical          :: addq2

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! which quadrant has to be added to the probability
if (xval(1).lt.0.d0) then
   addq2 = .true.    ! add quadrant 2
else
   addq2 = .false.   ! add quadrant 4
endif

! 
! ----------------------

! intersection of function with axis
yl = z / xval(nx)
yl = max(yval(1),min(yl,yval(ny)))
xu = z / yl
xu = max(xval(1),min(xu,xval(nx)))

! find first joint bin
call srvFndBin(xval,nx,xu,istr,.true.)
call srvFndBin(yval,ny,yl,jstr,.false.)

! set counters
iend =  1
istp = -1
jend = ny-1
jstp =  1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
if (addq2) then
   ! add probability of quadrant 2
   call srvFndBin(yval,ny,0.d0,jlst,.false.)     ! last bin of y used to calculate yc

   ! correction for y
   ! Necessary when y=0.d0 exists and does not coincide with a bin boundary
   if (yval(jlst).lt. 0.d0) then
      yj   = yval(jlst)
      yj1  = yval(jlst+1)
      pyj  = yden(jlst)
      pyj1 = yden(jlst+1)
      ry   = (pyj1-pyj)/(yj1-yj)
      p0y  = pyj + ry*(0.d0-yj)
      yc   = -1.*0.5*(p0y+pyj)*(0.d0-yj)
   else
      yc   = 0.0d0
   endif

   jlst = jlst-jstp
else
   ! add probability of quadrant 4
   jlst = 0     ! last bin of y used to calculate yc

   yc   = 0.d0
endif

! loop over all intersectd joint bins
cont = .true.
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xl,yu
   xl   = z / yj1
   nxtx = (xl.lt.xi)
   xl   = max(xi,min(xl,xi1))
   yu   = z / xl
   yu   = max(yj,min(yu,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   lny  = log(yu/yl)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc + p0x*p0y*(z*lny-xl*dy)                                 & 
           + p0x*ry *(z*dy-(1./2.)*xl*dy*(yu+yl))                  &
!           + (1./2.)*rx *p0y*(-z**2*(1./yu - 1./yl)-xl**2*dy)      &
!           + (1./2.)*rx *p0y*(z*(xu-xl) - xl**2*dy)                &
           + (1./2.)*rx *p0y*(xu*xl*dy - xl**2*dy)                &
           + (1./2.)*rx *ry *(z**2*lny - (1./2.)*xl**2*dy*(yu+yl)) 

   zc = zc + (1./4.)*(pyu+pyl)*(yu-yl)*(pxl+pxi)*(xl-xi)           &   !  |
           + (1./4.)*(pyl+pyj)*(yl-yj)*(pxu+pxi)*(xu-xi)               !  --

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   ! keep xl and yl fixed and let xu and yu be a function of z
   !    yu = z/xl       => d(yu)/d(z) = 1/xl
   !    xu = z/yl       => d(xu)/d(z) = 1/yl
   !    dy = yu - yl    => d(dy)/d(z) = 1/xl
   !    lny= log(yu/yl) => d(lny)/d(z)= 1/z
   ! First derivative (of area d) with respect to z is:
!   zd = zd + p0x*p0y*(lny + z/z - xl/xl)                                 & 
!           + p0x*ry *(dy+z/xl -(1./2.)*xl/xl*(yu+yl) -(1./2.)*xl*dy*(1/xl+yl) )    &
!           + (1./2.)*rx *p0y*(-z +2*z/yl   -xl**2/xl)      &
!           + (1./2.)*rx *ry *(2*z*lny + z  -(1./2.)*xl**2*(1/xl)*(yu+yl) -(1./2.)*xl**2*dy*(1/xl))
   zd = zd + p0x*p0y*(lny)        & 
           + p0x*ry*dy            &
           + rx *p0y*(xu-xl)      &
           + rx *ry *(z*lny)


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xu = xl
   yl = yu

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopMulQ1

! ******************************************************************************

!> description
!! Perform multiplcaion for quadrant 2 for rvBopMul3
!!
subroutine srvBopMulQ2(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 332*1000  ! routine error code for routine srvBopMulQ2

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,lny

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! 
! ----------------------

! intersection of function with axis
xu = z / yval(ny)
xu = max(xval(1),min(xu,xval(nx)))
yu = z / xu
yu = max(yval(1),min(yu,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xu,istr,.true.)
call srvFndBin(yval,ny,yu,jstr,.true.)

! set counters
iend =  1
istp = -1
call srvFndBin(yval,ny,0.d0,jend,.false.)
jstp = -1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
jlst = ny         ! last bin of y used to calculate yc


! loop over all intersectd joint bins
cont = .true.
yc   = 0.d0
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xl,yj
   yl   = z / xi
   nxtx = (yl.gt.yj)
   yl   = max(yj,min(yl,yj1))
   xl   = z / yl
   xl   = max(xi,min(xl,xi1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   lny  = log(yu/yl)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc + p0x*p0y*(z*lny-xl*dy)                                 & 
           + p0x*ry *(z*dy-(1./2.)*xl*dy*(yu+yl))                  &
!           + (1./2.)*rx *p0y*(-z**2*(1./yu - 1./yl)-xl**2*dy)      &
!           + (1./2.)*rx *p0y*(-z*(xu-xl) - xl**2*dy)                &
           + (1./2.)*rx *p0y*(xu*xl*dy - xl**2*dy)                &
           + (1./2.)*rx *ry *(z**2*lny - (1./2.)*xl**2*dy*(yu+yl)) 

   zc = zc + (1./4.)*(pyj1+pyu)*(yj1-yu)*(pxu+pxi)*(xu-xi)         &   !  --
           + (1./4.)*(pyu +pyl)*(yu -yl)*(pxl+pxi)*(xl-xi)             !  |

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   zd = zd + p0x*p0y*(lny)        & 
           + p0x*ry*dy            &
           - rx *p0y*(xu-xl)      &
           + rx *ry *(z*lny)


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xu = xl
   yu = yl

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopMulQ2

! ******************************************************************************

!> description
!! Perform multiplcaion for quadrant 3 for rvBopMul3
!! If xmax>0 then the probability of quadrant 4 is calculated as well in this routine
!! else the probability of quadrant 2 is also calculated
!!
subroutine srvBopMulQ3(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 333*1000  ! routine error code for routine srvBopMulQ3

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,lny

logical          :: addq2

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! which quadrant has to be added to the probability
if (xval(nx).gt.0.d0) then
   addq2 = .false.   ! add quadrant 4
else
   addq2 = .true.    ! add quadrant 2
endif

! 
! ----------------------

! intersection of function with axis
yu = z / xval(1)
yu = max(yval(1),min(yu,yval(ny)))
xl = z / yu
xl = max(xval(1),min(xl,xval(nx)))

! find first joint bin
call srvFndBin(xval,nx,xl,istr,.false.)
call srvFndBin(yval,ny,yu,jstr,.true.)

! set counters
iend = nx-1
istp =  1
jend =  1
jstp = -1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z

if (addq2) then
   jlst = ny     ! last bin of y used to calculate yc

   yc   = 0.d0
else
   call srvFndBin(yval,ny,0.d0,jlst,.true.)     ! last bin of y used to calculate yc

   ! correction for y
   ! Necessary when y=0.d0 exists and does not coincide with a bin boundary
   j = jlst
   if (yval(j+1).gt. 0.d0) then
      yj   = yval(j)
      yj1  = yval(j+1)
      pyj  = yden(j)
      pyj1 = yden(j+1)
      ry   = (pyj1-pyj)/(yj1-yj)
      p0y  = pyj + ry*(0.d0-yj)
      yc   = -1.*0.5*(pyj1+p0y)*(yj1-0.d0)
   else
      yc   = 0.0d0
   endif

   jlst = jlst-jstp
endif


! loop over all intersectd joint bins
cont = .true.
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xu,yl
   xu   = z / yj
   nxtx = (xu.gt.xi1)
   xu   = max(xi,min(xu,xi1))
   yl   = z / xu
   yl   = max(yj,min(yl,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   lny  = log(yu/yl)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc - p0x*p0y*(z*lny-xu*dy)                                 & 
           - p0x*ry *(z*dy-(1./2.)*xu*dy*(yu+yl))                  &
!           - (1./2.)*rx *p0y*(-z**2*(1./yu - 1./yl)-xu**2*dy)      &
!           - (1./2.)*rx *p0y*(z*(xu-xl) - xu**2*dy)                &
           - (1./2.)*rx *p0y*(xu*xl*dy - xu**2*dy)                &
           - (1./2.)*rx *ry *(z**2*lny - (1./2.)*xu**2*dy*(yu+yl)) 

   zc = zc + (1./4.)*(pyj1+pyu)*(yj1-yu)*(pxi1+pxl)*(xi1-xl)       &   !  --
           + (1./4.)*(pyu+pyl) *(yu-yl) *(pxi1+pxu)*(xi1-xu)           !   |

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   zd = zd - p0x*p0y*(lny)        & 
           - p0x*ry*dy            &
           - rx *p0y*(xu-xl)      &
           - rx *ry *(z*lny)


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xl = xu
   yu = yl

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopMulQ3

! ******************************************************************************

!> description
!! Perform multiplcaion for quadrant 4 for rvBopMul3
!!
subroutine srvBopMulQ4(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-22      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 334*1000  ! routine error code for routine srvBopMulQ4

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,lny

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! 
! ----------------------

! intersection of function with axis
xl = z / yval(1)
xl = max(xval(1),min(xl,xval(nx)))
yl = z / xl
yl = max(yval(1),min(yl,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xl,istr,.false.)
call srvFndBin(yval,ny,yl,jstr,.false.)

! set counters
iend = nx-1
istp =  1
call srvFndBin(yval,ny,0.d0,jend,.true.)
jstp =  1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
jlst =  0         ! last bin of y used to calculate yc


! loop over all intersectd joint bins
cont = .true.
i    = istr
j    = jstr
yc   = 0.d0
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xu,yu
   yu   = z / xi1
   nxtx = (yu.lt.yj1)
   yu   = max(yj,min(yu,yj1))
   xu   = z / yu
   xu   = max(xi,min(xu,xi1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   lny  = log(yu/yl)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc - p0x*p0y*(z*lny-xu*dy)                                 & 
           - p0x*ry *(z*dy-(1./2.)*xu*dy*(yu+yl))                  &
!           - (1./2.)*rx *p0y*(-z**2*(1./yu - 1./yl)-xu**2*dy)      &
!           - (1./2.)*rx *p0y*(-z*(xu-xl) - xu**2*dy)                &
           - (1./2.)*rx *p0y*(xu*xl*dy - xu**2*dy)                &
           - (1./2.)*rx *ry *(z**2*lny - (1./2.)*xu**2*dy*(yu+yl)) 

   zc = zc + (1./4.)*(pyu+pyl)*(yu-yl)*(pxi1+pxu)*(xi1-xu)         &   !   |
           + (1./4.)*(pyl+pyj)*(yl-yj)*(pxi1+pxl)*(xi1-xl)             !  --

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   zd = zd - p0x*p0y*(lny)        & 
           - p0x*ry*dy            &
           + rx *p0y*(xu-xl)      &
           - rx *ry *(z*lny)


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xl = xu
   yl = yu

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopMulQ4

! ******************************************************************************

!> description
!! Perform division for z=0 for rvBopDiv3
!!
subroutine srvBopDivQ0(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-29      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 340*1000  ! routine error code for routine srvBopDivQ0

! local variables
integer          :: j,i0

double precision :: rx,p0x,yj,yj1,pyj,pyj1,ry,p0y,yl,yu,dy,dy2,dy3


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0



! Cumulative probaility
! ---------------------

! add probabilities of quadrant 2 and quadrant 4
call srvBopPrbQ(xval,xden,nx,yval,yden,ny,zc,.false.,.true.,.false.,.true.,exitcode)



! Probability density
! -------------------

if (xval(1).le.0.d0 .and. xval(nx).ge.0.d0 .and. z.eq.0.d0) then
   ! find x-bin containing 0
   call srvFndBin(xval,nx,0.d0,i0,.true.)

   ! calc rx and p0x
   rx  = (xden(i0+1)-xden(i0))/(xval(i0+1)-xval(i0))
   p0x = xden(i0) + rx*(0.d0-xval(i0))

   ! loop over all bins of y
   j = 1
   do while(yval(j).lt.0.d0 .and. j.lt.ny)
      yj   = yval(j)
      yj1  = yval(j+1)
      pyj  = yden(j)
      pyj1 = yden(j+1)
      ry   = (pyj1-pyj)/(yj1-yj)
      p0y  = pyj + ry*(0.d0-yj) 

      yl   = yj
      yu   = min(0.d0,yj1)

      dy   = (yu - yl)
      dy2  = dy*(yu + yl)                                ! (yu**2 - yl**2)
      dy3  = dy*(3.*yl*yu+dy**2)                         ! (yu**3 - yl**3)

      zd   = zd - 0.5*p0x*p0y*dy2 -(1./3.)*p0x*ry*dy3

      j = j + 1
   end do

   j = ny
   do while(yval(j).gt.0.d0 .and. j.gt.1)
      j = j - 1

      yj   = yval(j)
      yj1  = yval(j+1)
      pyj  = yden(j)
      pyj1 = yden(j+1)
      ry   = (pyj1-pyj)/(yj1-yj)
      p0y  = pyj + ry*(0.d0-yj) 

      yl   = max(0.d0,yj)
      yu   = yj1

      dy   = (yu - yl)
      dy2  = dy*(yu + yl)                                ! (yu**2 - yl**2)
      dy3  = dy*(3.*yl*yu+dy**2)                         ! (yu**3 - yl**3)

      zd   = zd + 0.5*p0x*p0y*dy2 +(1./3.)*p0x*ry*dy3

   end do

endif

! end of program
return
end subroutine srvBopDivQ0

! ******************************************************************************

!> description
!! Perform division for quadrant 1 for rvBopDiv3
!!
subroutine srvBopDivQ1(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-25      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 341*1000  ! routine error code for routine srvBopDivQ1

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,dy2,dy3,dy4

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! intersection of function with axis
xu = z * yval(ny)
xu = max(xval(1),min(xu,xval(nx)))
yu = xu / z
yu = max(yval(1),min(yu,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xu,istr,.true.)
call srvFndBin(yval,ny,yu,jstr,.true.)

! set counters
iend =  1
istp = -1

call srvFndBin(yval,ny,0.d0,jend,.false.)
jstp = -1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
jlst = ny         ! last bin of y used to calculate yc


! loop over all intersectd joint bins
cont = .true.
yc   = 0.0d0
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xl,yl
   xl   = z * max(0.d0,yj)
   nxtx = (xl.lt.xi)
   xl   = max(xi,min(xl,xi1))
   yl   = xl / z
   yl   = max(yj,min(yl,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   dy2  = dy*(yu + yl)                                ! (yu**2 - yl**2)
   dy3  = dy*(3.*yl*yu+dy**2)                         ! (yu**3 - yl**3)
   dy4  = dy*(yl*(yl*(4.*yl+6.*dy)+4.*dy**2)+dy**3)   ! (yu**4 - yl**4)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc + p0x*p0y*((1./2.)*z*dy2 - xl*dy)                       & 
           + p0x*ry *((1./3.)*z*dy3 -(1./2.)*xl*dy2)               &
           + (1./2.)*rx *p0y*((1./3.)*z**2*dy3 - xl**2*dy)         &
           + (1./2.)*rx *ry *((1./4.)*z**2*dy4 - (1./2.)*xl**2*dy2) 

   zc = zc + (1./4.)*(pyj1+pyu)*(yj1-yu)*(pxu+pxi)*(xu-xi)         &   !  --
           + (1./4.)*(pyu +pyl)*(yu -yl)*(pxl+pxi)*(xl-xi)             !  |

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   ! keep xl and yu fixed and let xu and yl be a function of z
   ! First derivative (of area d) with respect to z is:
   zd = zd + (1.d0/2.d0)*p0x*p0y*  dy2        & 
           + (1.d0/3.d0)*p0x*ry *  dy3        &
           + (1.d0/3.d0)*rx *p0y*z*dy3        &
           + (1.d0/4.d0)*rx *ry *z*dy4


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xu = xl
   yu = yl

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yj   = yval(jj)
   yj1  = yval(jj+1)
   pyj  = yden(jj)
   pyj1 = yden(jj+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   yl   = max(0.d0,yj)
   pyl  = pyj + ry*(yl-yj)
!   yl   = yj
!   pyl  = pyj
   yu   = yj1
   pyu  = pyj1
!   yu   = max(0.d0,yj1)
!   pyu  = pyj + ry*(yu-yj)

   yc = yc + 0.5d0*(pyu+pyl)*(yu-yl)
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopDivQ1

! ******************************************************************************

!> description
!! Perform division for quadrant 2 for rvBopDiv3
!!
subroutine srvBopDivQ2(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-25      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 342*1000  ! routine error code for routine srvBopDivQ2

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc

double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,dy2,dy3,dy4

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! intersection of function with axis
xu = z * 0.d0
xu = max(xval(1),min(xu,xval(nx)))
yl = xu / z
yl = max(yval(1),min(yl,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xu,istr,.true.)
call srvFndBin(yval,ny,yl,jstr,.false.)

! set counters
iend =  1
istp = -1
jend = ny-1
jstp =  1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
call srvFndBin(yval,ny,0.d0,jlst,.false.)

! correction for y
! Necessary when y=0.d0 exists and does not coincide with a bin boundary
j = jlst
if (yval(j).lt. 0.d0) then
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)
   p0y  = pyj + ry*(0.d0-yj)
   yc   = -1.d0*0.5d0*(p0y+pyj)*(0.d0-yj)
else
   yc   = 0.0d0
endif

jlst = jlst-jstp


! loop over all intersectd joint bins
cont = .true.
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xl,yu
   xl   = z * yj1
   nxtx = (xl.lt.xi)
   xl   = max(xi,min(xl,xi1))
   yu   = xl / z
   yu   = max(yj,min(yu,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   dy2  = dy*(yu + yl)                                      ! (yu**2 - yl**2)
   dy3  = dy*(3.d0*yl*yu+dy**2)                             ! (yu**3 - yl**3)
   dy4  = dy*(yl*(yl*(4.d0*yl+6.d0*dy)+4.d0*dy**2)+dy**3)   ! (yu**4 - yl**4)


   ! calc cumulative probability
   ! ---------------------------
!   zc = zc + p0x*p0y*((1./2.)*z*dy2 -                    xl   *dy) & 
!           + p0x*ry *((1./3.)*z*dy3            - (1./2.)*xl   *dy2)&
!           + (1./2.)*rx *p0y*((1./3.)*z**2*dy3 -         xl**2*dy) &
!           + (1./2.)*rx *ry *((1./4.)*z**2*dy4 - (1./2.)*xl**2*dy2) 

   zc = zc + p0x*p0y*((1.d0/2.d0)*z*dy2 -                            xl   *dy) & 
           + p0x*ry *((1.d0/3.d0)*z*dy3                - (1.d0/2.d0)*xl   *dy2)&
           + (1.d0/2.d0)*rx *p0y*((1.d0/3.d0)*z**2*dy3 -             xl**2*dy) &
           + (1.d0/2.d0)*rx *ry *((1.d0/4.d0)*z**2*dy4 - (1.d0/2.d0)*xl**2*dy2) 

   zc = zc + (1.d0/4.d0)*(pyu+pyl)*(yu-yl)*(pxl+pxi)*(xl-xi)           &   !  |
           + (1.d0/4.d0)*(pyl+pyj)*(yl-yj)*(pxu+pxi)*(xu-xi)               !  --

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xi   = xval(ii)
      xi1  = xval(ii+1)
      pxi  = xden(ii)
      pxi1 = xden(ii+1)
      rx   = (pxi1-pxi)/(xi1-xi)

   !   xl   = min(0.d0,xi)
   !   pxl  = pxi + rx*(xl-xi)
      xl   = xi
      pxl  = pxi
   !   xu   = xi1
   !   pxu  = pxi1
      xu   = min(0.d0,xi1)
      pxu  = pxi + rx*(xu-xi)

      xc = xc + 0.5d0*(pxu+pxl)*(xu-xl)
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   ! keep xl and yu fixed and let xu and yl be a function of z
   ! First derivative (of area d) with respect to z is:
!   zd = zd + (1./2.)*p0x*p0y*  dy2        & 
!           + (1./3.)*p0x*ry *  dy3        &
!           + (1./3.)*rx *p0y*z*dy3        &
!           + (1./4.)*rx *ry *z*dy4

   zd = zd + (1.d0/2.d0)*p0x*p0y*  dy2    & 
           + (1.d0/3.d0)*p0x*ry *  dy3    &
           + (1.d0/3.d0)*rx *p0y*z*dy3    &
           + (1.d0/4.d0)*rx *ry *z*dy4


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xu = xl
   yl = yu

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopDivQ2

! ******************************************************************************

!> description
!! Perform division for quadrant 3 for rvBopDiv3
!!
subroutine srvBopDivQ3(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-25      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 343*1000  ! routine error code for routine srvBopDivQ3

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,dy2,dy3,dy4

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! intersection of function with axis
xl = z * yval(1)
xl = max(xval(1),min(xl,xval(nx)))
yl = xl / z
yl = max(yval(1),min(yl,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xl,istr,.false.)
call srvFndBin(yval,ny,yl,jstr,.false.)

! set counters
iend = nx-1
istp =  1
call srvFndBin(yval,ny,0.d0,jend,.true.)
jstp =  1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
jlst =  0         ! last bin of y used to calculate yc

! loop over all intersectd joint bins
cont = .true.
yc   = 0.0d0
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xu,yu
   xu   = z * min(0.d0,yj1)
   nxtx = (xu.gt.xi1)
   xu   = max(xi,min(xu,xi1))
   yu   = xu / z
   yu   = max(yj,min(yu,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   dy2  = dy*(yu + yl)                                ! (yu**2 - yl**2)
   dy3  = dy*(3.*yl*yu+dy**2)                         ! (yu**3 - yl**3)
   dy4  = dy*(yl*(yl*(4.*yl+6.*dy)+4.*dy**2)+dy**3)   ! (yu**4 - yl**4)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc - p0x*p0y*((1./2.)*z*dy2 - xu*dy)                       & 
           - p0x*ry *((1./3.)*z*dy3 -(1./2.)*xu*dy2)               &
           - (1./2.)*rx *p0y*((1./3.)*z**2*dy3 - xu**2*dy)         &
           - (1./2.)*rx *ry *((1./4.)*z**2*dy4 - (1./2.)*xu**2*dy2) 

   zc = zc + (1./4.)*(pyu+pyl)*(yu-yl)*(pxi1+pxu)*(xi1-xu)         &   !   |
           + (1./4.)*(pyl+pyj)*(yl-yj)*(pxi1+pxl)*(xi1-xl)             !  --

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   ! keep xl and yu fixed and let xu and yl be a function of z
   ! First derivative (of area d) with respect to z is:
   zd = zd - (1.d0/2.d0)*p0x*p0y*  dy2        & 
           - (1.d0/3.d0)*p0x*ry *  dy3        &
           - (1.d0/3.d0)*rx *p0y*z*dy3        &
           - (1.d0/4.d0)*rx *ry *z*dy4


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xl = xu
   yl = yu

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yj   = yval(jj)
   yj1  = yval(jj+1)
   pyj  = yden(jj)
   pyj1 = yden(jj+1)
   ry   = (pyj1-pyj)/(yj1-yj)

!   yl   = min(0.d0,yj)
!   pyl  = pyj + ry*(yl-yj)
   yl   = yj
   pyl  = pyj
!   yu   = yj1
!   pyu  = pyj1
   yu   = min(0.d0,yj1)
   pyu  = pyj + ry*(yu-yj)

   yc = yc + 0.5d0*(pyu+pyl)*(yu-yl)
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopDivQ3

! ******************************************************************************

!> description
!! Perform division for quadrant 4 for rvBopDiv3
!!
subroutine srvBopDivQ4(xval,xden,nx,yval,yden,ny,z,zd,zc,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-25      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: ny          !< number of elements of variable Y
double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
                                                               !! \n yval and yden together form the PDF of Y

double precision               , intent(in)     :: z           !< z-value
double precision               , intent(inout)  :: zd          !< summation for probability density
double precision               , intent(inout)  :: zc          !< summation for cumulative probability

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! routine error code
integer,parameter    :: rtnerrcode = 140000000 + 344*1000  ! routine error code for routine srvBopDivQ4

! local variables
integer          :: i,j,istr,iend,istp,ilst,jstr,jend,jstp,jlst,ii,jj

double precision :: xi,xi1,pxi,pxi1,rx, xl,xu,pxl,pxu, xc
double precision :: yj,yj1,pyj,pyj1,ry, yl,yu,pyl,pyu, yc
double precision :: dy
double precision :: p0x,p0y,dy2,dy3,dy4

logical          :: cont

logical          :: nxtx     ! .true. : shift to next x bin (i=i+istp)
                             ! .false.: shift to next y bin (j=j+jstp)

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! intersection of function with axis
xl = z * 0.d0
xl = max(xval(1),min(xl,xval(nx)))
yu = xl / z
yu = max(yval(1),min(yu,yval(ny)))

! find first joint bin
call srvFndBin(xval,nx,xl,istr,.false.)
call srvFndBin(yval,ny,yu,jstr,.true.)

! set counters
iend = nx-1
istp =  1
jend =  1
jstp = -1

ilst = istr-istp  ! last bin of x used to calculate probabiility of not intersected
                  ! joint  with max-z < z
call srvFndBin(yval,ny,0.d0,jlst,.true.)  ! last bin of y used to calculate yc

! correction for y
! Necessary when y=0.d0 exists and does not coincide with a bin boundary
j = jlst
if (yval(j+1).gt. 0.d0) then
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)
   p0y  = pyj + ry*(0.d0-yj)
   yc   = -1.*0.5*(pyj1+p0y)*(yj1-0.d0)
else
   yc   = 0.0d0
endif

jlst = jlst-jstp


! loop over all intersectd joint bins
cont = .true.
i    = istr
j    = jstr
do while (cont)

   ! get function parameters
   xi   = xval(i)
   xi1  = xval(i+1)
   pxi  = xden(i)
   pxi1 = xden(i+1)
   rx   = (pxi1-pxi)/(xi1-xi)
   yj   = yval(j)
   yj1  = yval(j+1)
   pyj  = yden(j)
   pyj1 = yden(j+1)
   ry   = (pyj1-pyj)/(yj1-yj)

   ! calc xu,yl
   xu   = z * yj
   nxtx = (xu.gt.xi1)
   xu   = max(xi,min(xu,xi1))
   yl   = xu / z
   yl   = max(yj,min(yl,yj1))

   dy   = yu - yl
   pxl  = pxi + rx*(xl-xi)
   pxu  = pxi + rx*(xu-xi)
   pyl  = pyj + ry*(yl-yj)
   pyu  = pyj + ry*(yu-yj)

   p0x  = pxi + rx*(0.d0-xi)
   p0y  = pyj + ry*(0.d0-yj)
   dy2  = dy*(yu + yl)                                ! (yu**2 - yl**2)
   dy3  = dy*(3.*yl*yu+dy**2)                         ! (yu**3 - yl**3)
   dy4  = dy*(yl*(yl*(4.*yl+6.*dy)+4.*dy**2)+dy**3)   ! (yu**4 - yl**4)


   ! calc cumulative probability
   ! ---------------------------
   zc = zc - p0x*p0y*((1./2.)*z*dy2 - xu*dy)                       & 
           - p0x*ry *((1./3.)*z*dy3 -(1./2.)*xu*dy2)               &
           - (1./2.)*rx *p0y*((1./3.)*z**2*dy3 - xu**2*dy)         &
           - (1./2.)*rx *ry *((1./4.)*z**2*dy4 - (1./2.)*xu**2*dy2) 

   zc = zc + (1./4.)*(pyj1+pyu)*(yj1-yu)*(pxi1+pxl)*(xi1-xl)       &   !  --
           + (1./4.)*(pyu+pyl) *(yu-yl) *(pxi1+pxu)*(xi1-xu)           !   |

   ! calc cumulative probability for y until bin j-1
   do jj=jlst+jstp,j-jstp,jstp
      yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
   enddo
   jlst = j-jstp

   ! add probability for all joint bins (i,<j)
   xc = 0.0d0
   do ii=ilst+istp,i,istp
      xi   = xval(ii)
      xi1  = xval(ii+1)
      pxi  = xden(ii)
      pxi1 = xden(ii+1)
      rx   = (pxi1-pxi)/(xi1-xi)

      xl   = max(0.d0,xi)
      pxl  = pxi + rx*(xl-xi)
   !   xl   = xi
   !   pxl  = pxi
      xu   = xi1
      pxu  = pxi1
   !   xu   = min(0.d0,xi1)
   !   pxu  = pxi + rx*(xu-xi)

      xc = xc + 0.5d0*(pxu+pxl)*(xu-xl)
   enddo
   zc = zc + xc*yc
   ilst = i

   ! calc probability density
   ! ------------------------
   ! keep xl and yu fixed and let xu and yl be a function of z
   ! First derivative (of area d) with respect to z is:
   zd = zd - (1.d0/2.d0)*p0x*p0y*  dy2        & 
           - (1.d0/3.d0)*p0x*ry *  dy3        &
           - (1.d0/3.d0)*rx *p0y*z*dy3        &
           - (1.d0/4.d0)*rx *ry *z*dy4


   ! next joint bin
   ! --------------
   if (nxtx) then
      i = i + istp
      if (istp*(i-iend).gt.0) cont=.false.
   else
      j = j + jstp
      if (jstp*(j-jend).gt.0) cont=.false.
   endif

   xl = xu
   yu = yl

enddo


! complete yc
! calc cumulative probability for y until bin j-1
do jj=jlst+jstp,jend,jstp
   yc = yc + 0.5d0*(yden(jj+1)+yden(jj))*(yval(jj+1)-yval(jj))
enddo
jlst = jend

! add forgotten cumulative probability
! When the line z = x - y intersects the line y=yval(1) at some bin i>1 then
! not all probbaility is added to zc. 
xc = 0.0d0
do ii=ilst+istp,iend,istp
   xc = xc + 0.5d0*(xden(ii+1)+xden(ii))*(xval(ii+1)-xval(ii))
enddo
zc = zc + xc*yc  ! yc should be 1. by now, but in case of roundoff errors
                 ! the calculated value yc is used


! end of program
return
end subroutine srvBopDivQ4


