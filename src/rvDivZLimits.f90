
! routines to get the Z-limits in case of division: Z = X / Y
! If the variable Y contains 0 in its domain and the probability
! density is not equal 0 in the neighbourhood then the limits of Z are
! +/- inifinite. The next routines make a good guess for these limits.
!
! Next routines are available:
! rvDivZLimit   calculate the limits for well defined PDFs,
!               i.e. non-zero densities and the CDF sums up to 1


!> Calculate the Z-limits for division of two Random Variables \n
!!    Z = X / Y
!!
!! Implementation of the division of two independent piecewise linear PDFs
!!
subroutine rvDivZLimits(xval,xden,nx,yval,yden,ny,minprob,zlim,exitcode) 

! History
! programmer------------date------------version---------------------------------
! Aris                  2019-11-22      initial version, copied from rvBopDiv2.


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

double precision,               intent(in)      :: minprob     !< minimum amount of probability (fraction) to 
                                                               !! be retained between the z-limits

double precision, dimension(2), intent(out)     :: zlim        !< found z-limits

!integer                        , intent(in)     :: nz          !< number of elements of variable Z
!double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
!double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
!double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

! local variables
integer          :: n,iter

double precision :: zmin,zmax,z1,z2,z3,z4,yeps

double precision :: ymin,ymax

!double precision, parameter:: minprob=0.999  ! minimum amount of probability of Z within initial z-values
!                                             ! This parameter is used when 0 is in the domain of Y

double precision :: tcum,tcumx,tcumy,tcumz,cumlow,cumhig,dval

logical          :: zfinite,cont

double precision, dimension(2) :: zval,zden,zcum

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! start values of z
yeps = (yval(ny)-yval(1))*(10.**(-3))

if (yval(1).gt.0.d0 .or. yval(ny).lt.0.d0) then
   ! y<0 OR y>0
   ! zmin and zmax will be finite
   zfinite = .true.
   ymin = yval(1)
   ymax = yval(ny)
else
   ! zmin and/or zmax is infinite
   zfinite = .false.
   if (yval(1).lt.0.d0) then
      if (yval(ny).gt.0.d0) then
         ! y(1)<0 AND y(ny)>0
         ymax =  yeps
         ymin = -yeps
      else
         ! y<=0
         ymax = min(yval(ny),-yeps)
         ymin = min(yval(1) ,-yeps)
      endif
   else
      ! y>=0
      ymax = max(yval(ny),yeps)
      ymin = max(yval(1) ,yeps)
   endif
endif

! zmin and zmax available
z1 = xval(nx)/ ymax
z2 = xval(1) / ymax
z3 = xval(1) / ymin
z4 = xval(nx)/ ymin

zmin     = min(z1,z2,z3,z4) ! mimimum z-value for all joint bins
zmax     = max(z1,z2,z3,z4) ! maximum z-value for all joint bins

n = 2
zval(1)  = zmin
!zval(2)  = (zmin+zmax)/2.
zval(2)  = zmax

if (.not.zfinite) then
   ! calculate the minimum probability to be declared by Z
   ! Not all PDFs do integrate exactly to 1, so the integral has to be calculated first
   call rvStatProb(xval,xden,nx,tcumx,exitcode)
   call rvStatProb(yval,yden,ny,tcumy,exitcode)
   tcum   = tcumx*tcumy
   cumlow = 0.5*(1.d0-minprob)*tcum
   cumhig = tcum - cumlow
endif

! calculate density and probability for initial values of Z
! If necessary, zmin and zmax are extended
cont = .true.
iter = 1
do while(cont)
   ! calc zden and zcum
   call rvBopDiv3(xval,xden,nx,yval,yden,ny,zval,zden,zcum,n,exitcode)
   if(exitcode.ne.0) return

   ! check or iteration is needed
   if (zfinite) then
      cont = .false.
   else
      iter = iter + 1
      if (iter .lt. 100) then   ! maximum number of iterations to prevent a deadlock
         tcumz = zcum(n)-zcum(1)
         if (tcumz.gt.(tcum*minprob)) then
            cont = .false.
         else
            ! check which boundary should be extended
            ! Extension is calculated as width (dval) of a bin with
            ! density zden and probability zcum:
            !    dval*zden = cum
            !
            cont = .false.
            ! lower boundary
            if (zcum(1).gt.cumlow) then
               if (zden(1).gt.0.d0) then
                  dval = zcum(1)/zden(1)
                  zval(1) = zval(1) - dval
                  cont = .true.
               endif
            endif
            ! upper boundary
            if (zcum(n).lt.cumhig) then
               if (zden(n).gt.0.d0) then
                  dval = (tcum-zcum(n))/zden(n)
                  zval(n) = zval(n) + dval
                  cont = .true.
               endif
            endif
!            ! new middle value
!            if (cont) then
!               zval(2) = (zval(1)+zval(n))/2.d0
!            endif
         endif
      else
         cont = .false.
      endif
   endif

enddo

! assign results
zlim(1) = zval(1)
zlim(2) = zval(n)

! end of program
return
end subroutine rvDivZLimits

! ******************************************************************************

!!!###!!!!> Calculate the Z-limits for division of two Random Variables \n
!!!###!!!!!    Z = X / Y
!!!###!!!!!
!!!###!!!!! Implementation of the division of two independent piecewise linear PDFs
!!!###!!!!! This is the Raw version of rvDivZLimits.
!!!###!!!!!
!!!###!!!subroutine rvDivZLimitsR(xval,xden,nx,yval,yden,ny,minprob,zlim,exitcode) 
!!!###!!!
!!!###!!!! History
!!!###!!!! programmer------------date------------version---------------------------------
!!!###!!!! Aris                  2019-11-22      initial version, copied from rvBopDiv2.
!!!###!!!
!!!###!!!
!!!###!!!! declaration section
!!!###!!!! ------------------------------------------------------------------------------
!!!###!!!
!!!###!!!
!!!###!!!implicit none
!!!###!!!
!!!###!!!
!!!###!!!! arguments
!!!###!!!integer                        , intent(in)     :: nx          !< number of elements of variable X
!!!###!!!double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
!!!###!!!double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
!!!###!!!                                                               !! \n xval and xden together form the PDF of X
!!!###!!!
!!!###!!!integer                        , intent(in)     :: ny          !< number of elements of variable Y
!!!###!!!double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
!!!###!!!double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
!!!###!!!                                                               !! \n yval and yden together form the PDF of Y
!!!###!!!
!!!###!!!double precision,               intent(in)      :: minprob     !< minimum amount of probability (fraction) to 
!!!###!!!                                                               !! be retained between the z-limits
!!!###!!!
!!!###!!!double precision, dimension(2), intent(out)     :: zlim        !< found z-limits
!!!###!!!
!!!###!!!!integer                        , intent(in)     :: nz          !< number of elements of variable Z
!!!###!!!!double precision, dimension(nz), intent(out)    :: zval        !< values of the PDF of variable Z
!!!###!!!!double precision, dimension(nz), intent(out)    :: zden        !< probability densities of the PDF of variable Z
!!!###!!!!double precision, dimension(nz), intent(out)    :: zcum        !< cumulative probabilities of the CDF of Z
!!!###!!!
!!!###!!!integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK
!!!###!!!
!!!###!!!! local variables
!!!###!!!integer          :: n,bin0,ibin
!!!###!!!
!!!###!!!double precision :: zmin,zmax,z1,z2,z3,z4,yeps
!!!###!!!
!!!###!!!double precision :: ymin,ymax
!!!###!!!
!!!###!!!!double precision, parameter:: minprob=0.999  ! minimum amount of probability of Z within initial z-values
!!!###!!!!                                             ! This parameter is used when 0 is in the domain of Y
!!!###!!!double precision, parameter:: minyeps=1.d-100 ! minimum valkue for yeps
!!!###!!!
!!!###!!!double precision :: tcumy,cumlow
!!!###!!!double precision :: p0,r
!!!###!!!
!!!###!!!logical          :: zfinite
!!!###!!!
!!!###!!!double precision, dimension(2) :: zval,zden,zcum
!!!###!!!
!!!###!!!! program section
!!!###!!!! ------------------------------------------------------------------------------
!!!###!!!
!!!###!!!! init
!!!###!!!exitcode = 0
!!!###!!!
!!!###!!!! start values of z
!!!###!!!yeps = (yval(ny)-yval(1))*(10.**(-3))
!!!###!!!
!!!###!!!if (yval(1).gt.0.d0 .or. yval(ny).lt.0.d0) then
!!!###!!!   ! y<0 OR y>0
!!!###!!!   ! zmin and zmax will be finite
!!!###!!!   zfinite = .true.
!!!###!!!   ymin = yval(1)
!!!###!!!   ymax = yval(ny)
!!!###!!!else
!!!###!!!   ! zmin and/or zmax is infinite
!!!###!!!   zfinite = .false.
!!!###!!!
!!!###!!!   ! calculate the minimum probability to be around 0 for yeps
!!!###!!!   call rvStatAbsArea(yval,yden,ny,tcumy,exitcode)
!!!###!!!   cumlow = (1.d0-minprob)*tcumy
!!!###!!!
!!!###!!!   ! find bin containing 0.0
!!!###!!!   call srvFndBin(yval,ny,0.d0,bin0,.false.)
!!!###!!!
!!!###!!!   ! try two bins
!!!###!!!   do ibin=max(1,bin0-1),bin0
!!!###!!!      if (.not. yval(ibin+1).lt.0.d0) then
!!!###!!!         r  = (abs(yden(ibin+1))-abs(yden(ibin)))/(yval(ibin+1)-yval(ibin))  ! slope of abs densities
!!!###!!!         p0 = yden(ibin) + r*(0.d0-yval(ibin))
!!!###!!!         r  = abs(r)  ! absolute value gives the smallest value of yeps and makes the calculations easyer
!!!###!!!         if (r.gt.0.d0 .or. p0.gt.0.d0) then
!!!###!!!            if (r.gt.0.d0) then
!!!###!!!               ! solve: cumlow = p0*yeps + 0.5*r*yeps^2
!!!###!!!               yeps = min(yeps,max(minyeps,(sqrt(p0**2 + 2.d0*r*cumlow)-p0)/r))
!!!###!!!            else
!!!###!!!               ! solve: cumlow = p0*yeps
!!!###!!!               yeps = min(yeps,max(minyeps,cumlow/p0))
!!!###!!!            endif
!!!###!!!         else
!!!###!!!            ! use the bin-boundary as yeps
!!!###!!!            if (yval(ibin).lt.0.d0) then
!!!###!!!               ! bin [<0,?]
!!!###!!!               if (yval(ibin+1).gt.0.d0) then
!!!###!!!                  ! bin [<0,>0]
!!!###!!!                  yeps = min(yeps,max(minyeps,abs(yval(ibin))))
!!!###!!!                  yeps = min(yeps,max(minyeps,abs(yval(ibin+1))))
!!!###!!!               else
!!!###!!!                  ! bin [<0,0]
!!!###!!!                  yeps = min(yeps,max(minyeps,abs(yval(ibin))))
!!!###!!!               endif
!!!###!!!            else
!!!###!!!               ! bin [0,>0]
!!!###!!!               yeps = min(yeps,max(minyeps,abs(yval(ibin+1))))
!!!###!!!            endif
!!!###!!!         endif
!!!###!!!      endif
!!!###!!!   enddo
!!!###!!!
!!!###!!!   !
!!!###!!!   if (yval(1).lt.0.d0) then
!!!###!!!      if (yval(ny).gt.0.d0) then
!!!###!!!         ! y(1)<0 AND y(ny)>0
!!!###!!!         ymax =  yeps
!!!###!!!         ymin = -yeps
!!!###!!!      else
!!!###!!!         ! y<=0
!!!###!!!         ymax = min(yval(ny),-yeps)
!!!###!!!         ymin = min(yval(1) ,-yeps)
!!!###!!!      endif
!!!###!!!   else
!!!###!!!      ! y>=0
!!!###!!!      ymax = max(yval(ny),yeps)
!!!###!!!      ymin = max(yval(1) ,yeps)
!!!###!!!   endif
!!!###!!!   ymax =  yeps
!!!###!!!   ymin = -yeps
!!!###!!!endif
!!!###!!!
!!!###!!!! zmin and zmax available
!!!###!!!z1 = xval(nx)/ ymax
!!!###!!!z2 = xval(1) / ymax
!!!###!!!z3 = xval(1) / ymin
!!!###!!!z4 = xval(nx)/ ymin
!!!###!!!
!!!###!!!zmin     = min(z1,z2,z3,z4) ! mimimum z-value for all joint bins
!!!###!!!zmax     = max(z1,z2,z3,z4) ! maximum z-value for all joint bins
!!!###!!!
!!!###!!!
!!!###!!!n = 2
!!!###!!!zval(1)  = zmin
!!!###!!!!zval(2)  = (zmin+zmax)/2.
!!!###!!!zval(2)  = zmax
!!!###!!!
!!!###!!!! calc zden and zcum
!!!###!!!call rvBopDiv3R(xval,xden,nx,yval,yden,ny,zval,zden,zcum,n,exitcode)
!!!###!!!if(exitcode.ne.0) return
!!!###!!!
!!!###!!!
!!!###!!!! assign results
!!!###!!!zlim(1) = zval(1)
!!!###!!!zlim(2) = zval(n)
!!!###!!!
!!!###!!!! end of program
!!!###!!!return
!!!###!!!end subroutine rvDivZLimitsR
!!!###!!!
!!!###!!!! ******************************************************************************
!!!###!!!
!!!###!!!!> Calculate the Z-limits for division of two Random Variables \n
!!!###!!!!!    Z = X / Y
!!!###!!!!!
!!!###!!!!! Implementation of the division of two independent piecewise linear PDFs
!!!###!!!!!
!!!###!!!subroutine rvDivZLimitsAbs(xval,xden,nx,yval,yden,ny,minprob,zlim,exitcode) 
!!!###!!!
!!!###!!!! History
!!!###!!!! programmer------------date------------version---------------------------------
!!!###!!!! Aris                  2019-11-24      adapted copy of rvDivZLimits
!!!###!!!
!!!###!!!
!!!###!!!! declaration section
!!!###!!!! ------------------------------------------------------------------------------
!!!###!!!
!!!###!!!
!!!###!!!implicit none
!!!###!!!
!!!###!!!
!!!###!!!! arguments
!!!###!!!integer                        , intent(in)     :: nx          !< number of elements of variable X
!!!###!!!double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
!!!###!!!double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
!!!###!!!                                                               !! \n xval and xden together form the PDF of X
!!!###!!!
!!!###!!!integer                        , intent(in)     :: ny          !< number of elements of variable Y
!!!###!!!double precision, dimension(ny), intent(in)     :: yval        !< values of the PDF of variable Y
!!!###!!!double precision, dimension(ny), intent(in)     :: yden        !< probability densities of the PDF of variable Y
!!!###!!!                                                               !! \n yval and yden together form the PDF of Y
!!!###!!!
!!!###!!!double precision,               intent(in)      :: minprob     !< minimum amount of probability (fraction) to 
!!!###!!!                                                               !! be retained between the z-limits
!!!###!!!
!!!###!!!double precision, dimension(2), intent(out)     :: zlim        !< found z-limits
!!!###!!!
!!!###!!!
!!!###!!!integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK
!!!###!!!
!!!###!!!! local variables
!!!###!!!integer          :: n,iter
!!!###!!!
!!!###!!!double precision :: zmin,zmax,z1,z2,z3,z4,yeps
!!!###!!!
!!!###!!!double precision :: ymin,ymax
!!!###!!!
!!!###!!!double precision :: tcum,tcumx,tcumy,tcumz,cumlow,cumhig,dval
!!!###!!!
!!!###!!!logical          :: zfinite,cont
!!!###!!!
!!!###!!!double precision, dimension(2) :: zval,zden,zcum
!!!###!!!
!!!###!!!! program section
!!!###!!!! ------------------------------------------------------------------------------
!!!###!!!
!!!###!!!! init
!!!###!!!exitcode = 0
!!!###!!!
!!!###!!!! start values of z
!!!###!!!yeps = (yval(ny)-yval(1))*(10.**(-3))
!!!###!!!
!!!###!!!if (yval(1).gt.0.d0 .or. yval(ny).lt.0.d0) then
!!!###!!!   ! y<0 OR y>0
!!!###!!!   ! zmin and zmax will be finite
!!!###!!!   zfinite = .true.
!!!###!!!   ymin = yval(1)
!!!###!!!   ymax = yval(ny)
!!!###!!!else
!!!###!!!   ! zmin and/or zmax is infinite
!!!###!!!   zfinite = .false.
!!!###!!!   if (yval(1).lt.0.d0) then
!!!###!!!      if (yval(ny).gt.0.d0) then
!!!###!!!         ! y(1)<0 AND y(ny)>0
!!!###!!!         ymax =  yeps
!!!###!!!         ymin = -yeps
!!!###!!!      else
!!!###!!!         ! y<=0
!!!###!!!         ymax = min(yval(ny),-yeps)
!!!###!!!         ymin = min(yval(1) ,-yeps)
!!!###!!!      endif
!!!###!!!   else
!!!###!!!      ! y>=0
!!!###!!!      ymax = max(yval(ny),yeps)
!!!###!!!      ymin = max(yval(1) ,yeps)
!!!###!!!   endif
!!!###!!!endif
!!!###!!!
!!!###!!!! zmin and zmax available
!!!###!!!z1 = xval(nx)/ ymax
!!!###!!!z2 = xval(1) / ymax
!!!###!!!z3 = xval(1) / ymin
!!!###!!!z4 = xval(nx)/ ymin
!!!###!!!
!!!###!!!zmin     = min(z1,z2,z3,z4) ! mimimum z-value for all joint bins
!!!###!!!zmax     = max(z1,z2,z3,z4) ! maximum z-value for all joint bins
!!!###!!!
!!!###!!!n = 2
!!!###!!!zval(1)  = zmin
!!!###!!!!zval(2)  = (zmin+zmax)/2.
!!!###!!!zval(2)  = zmax
!!!###!!!
!!!###!!!if (.not.zfinite) then
!!!###!!!   ! calculate the minimum probability to be declared by Z
!!!###!!!   ! Not all PDFs do integrate exactly to 1, so the integral has to be calculated first
!!!###!!!   call rvStatAbsArea(xval,xden,nx,tcumx,exitcode)
!!!###!!!   call rvStatAbsArea(yval,yden,ny,tcumy,exitcode)
!!!###!!!   tcum   = tcumx*tcumy
!!!###!!!   cumlow = 0.5*(1.d0-minprob)*tcum
!!!###!!!   cumhig = tcum - cumlow
!!!###!!!endif
!!!###!!!
!!!###!!!! calculate density and probability for initial values of Z
!!!###!!!! If necessary, zmin and zmax are extended
!!!###!!!cont = .true.
!!!###!!!iter = 1
!!!###!!!do while(cont)
!!!###!!!   ! calc zden and zcum
!!!###!!!   call rvBopDiv3Abs(xval,xden,nx,yval,yden,ny,zval,zden,zcum,n,exitcode)
!!!###!!!   if(exitcode.ne.0) return
!!!###!!!
!!!###!!!   ! check or iteration is needed
!!!###!!!   if (zfinite) then
!!!###!!!      cont = .false.
!!!###!!!   else
!!!###!!!      iter = iter + 1
!!!###!!!      if (iter .lt. 100) then   ! maximum number of iterations to prevent a deadlock
!!!###!!!         tcumz = zcum(n)-zcum(1)
!!!###!!!         if (tcumz.gt.(tcum*minprob)) then
!!!###!!!            cont = .false.
!!!###!!!         else
!!!###!!!            ! check which boundary should be extended
!!!###!!!            ! Extension is calculated as width (dval) of a bin with
!!!###!!!            ! density zden and probability zcum:
!!!###!!!            !    dval*zden = cum
!!!###!!!            !
!!!###!!!            cont = .false.
!!!###!!!            ! lower boundary
!!!###!!!            if (zcum(1).gt.cumlow) then
!!!###!!!               if (zden(1).gt.0.d0) then
!!!###!!!                  dval = zcum(1)/zden(1)
!!!###!!!                  zval(1) = zval(1) - dval
!!!###!!!                  cont = .true.
!!!###!!!               endif
!!!###!!!            endif
!!!###!!!            ! upper boundary
!!!###!!!            if (zcum(n).lt.cumhig) then
!!!###!!!               if (zden(n).gt.0.d0) then
!!!###!!!                  dval = (tcum-zcum(n))/zden(n)
!!!###!!!                  zval(n) = zval(n) + dval
!!!###!!!                  cont = .true.
!!!###!!!               endif
!!!###!!!            endif
!!!###!!!!            ! new middle value
!!!###!!!!            if (cont) then
!!!###!!!!               zval(2) = (zval(1)+zval(n))/2.d0
!!!###!!!!            endif
!!!###!!!         endif
!!!###!!!      else
!!!###!!!         cont = .false.
!!!###!!!      endif
!!!###!!!   endif
!!!###!!!
!!!###!!!enddo
!!!###!!!
!!!###!!!! assign results
!!!###!!!zlim(1) = zval(1)
!!!###!!!zlim(2) = zval(n)
!!!###!!!
!!!###!!!! end of program
!!!###!!!return
!!!###!!!end subroutine rvDivZLimitsAbs

