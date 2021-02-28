
!> \file

!> Find next z-value based on probability density and cumulative probality
subroutine srvNxtZ(zval,zden,zcum,nz,iz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-27      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nz          !< number of elements of variable Z   \n
                                                               !! On entry: number of known z-values \n
                                                               !!           is nz-1
double precision, dimension(nz), intent(inout)  :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(inout)  :: zden        !< probability densities of the PDF of variable Z
double precision, dimension(nz), intent(inout)  :: zcum        !< cumulative probabilities of the CDF of Z

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

integer                        , intent(out)    :: iz          !< position of the inserted z-value


! routine error code
!   : 
integer,parameter    :: rtnerrcode = 140000000 +  82*1000  ! routine error code for routine srvNxtZ

! local variables
integer          :: k

double precision :: dmax,cumd,cumc,diff,ldif,fct

logical          :: lden

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! which type of objective used
if (nz.gt.5) then
   if (mod(nz,4) .eq.0) then
      lden = .true.    ! find largest step in probability density
   else
      lden = .false.   ! find largest difference in cumulative probaility
   endif
else
   lden = .false.
endif


if (lden) then
   ! find bin with maximum difference between successive priobability denisties
   dmax = 0.d0
   iz   = 1
   do k=1,nz-2
      diff = zden(k+1) - zden(k)
      if (abs(diff).gt.dmax) then
         iz   = k
         dmax = abs(diff)
      endif
   enddo
   k = iz
   cumd = 0.5*(zden(k+1)+zden(k))*(zval(k+1)-zval(k))
   cumc = zcum(k+1)-zcum(k)
   ldif = cumd-cumc
else
   ! find bin with maximum difference between cumulative probabilities
   dmax = 0.d0
   iz   = 1
   do k=1,nz-2
      cumd = 0.5*(zden(k+1)+zden(k))*(zval(k+1)-zval(k))
      cumc = zcum(k+1)-zcum(k)
      diff = cumd-cumc
      if (abs(diff).gt.dmax) then
         iz   = k
         dmax = abs(diff)
         ldif = diff
      endif
   enddo
endif


! create empty value at position iz+1
do k=nz-1,iz+1,-1
   zval(k+1)=zval(k)
   zden(k+1)=zden(k)
   zcum(k+1)=zcum(k)
enddo
iz = iz + 1


! calculate new z-value at centre of mass at bin iz
if (ldif.gt.0.d0) then
   ! density overestimates the probability
   fct = zden(iz+1)+zden(iz-1)
   if (fct.ne.0.d0) then
      fct = (2.*zden(iz+1)+zden(iz-1))/(3.d0*fct)
      fct = max(1./3.,min(fct,2./3.))
   else
      fct = 1./2.
   endif
else
   fct = 1./2.
endif
zval(iz) = zval(iz-1) + fct*(zval(iz+1)-zval(iz-1))
zden(iz) = 0.d0
zcum(iz) = 0.d0


! end of program
return
end subroutine srvNxtZ

! ******************************************************************************

!> Find next z-value based on probability density only
subroutine srvNxtZden(zval,zden,nz,iz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-07-15      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nz          !< number of elements of variable Z   \n
                                                               !! On entry: number of known z-values \n
                                                               !!           is nz-1
double precision, dimension(nz), intent(inout)  :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(inout)  :: zden        !< probability densities of the PDF of variable Z
                                                               !! \n zval and zden together form the PDF of Z


integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

integer                        , intent(out)    :: iz          !< position of the inserted z-value

! routine error code
!   : 
integer,parameter    :: rtnerrcode = 140000000 +  83*1000  ! routine error code for routine srvNxtZden

! local variables
integer           :: k,izd,izv

double precision  :: diff,dmax,fct,ldif,f,fd,fv,dthr,vthr

integer, parameter:: method=4   ! 1: bin with max abs density difference  dzden(k)
                                ! 2:          max diff density fraction   dzden(k)/dzden(k+1)
                                ! 3:          max diff value fraction     dzval(k)/dzval(k+1)
                                ! 4: combination of 2 and 3

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
iz       = 1
izd      = -1
izv      = -1
fd       = -1.d0
fv       = -1.d0


! find bin with maximum difference between successive priobability denisties
if (method.eq.1) then
   dmax = 0.d0
   do k=1,nz-2
      diff = zden(k+1) - zden(k)
      if (abs(diff).gt.dmax) then
         iz   = k
         dmax = abs(diff)
      endif
   enddo
endif


! find bin with max diff density fraction   dzden(k)/dzden(k+1)
if (method.eq.2 .or. method.eq.4) then
   ! threshold value den
   dthr = 0.d0
   do k=1,nz-2
      dthr = max(dthr,abs(zden(k+1) - zden(k)))
   enddo
   dthr = 0.5*dthr

   izd  = 1
   fd   = 0.d0
   ldif = abs(zden(2) - zden(1))
   do k=2,nz-2
      diff = abs(zden(k+1) - zden(k))
      ! backward: dzden(k)/dzden(k-1)
      if (ldif.gt.0.d0) then
         f = diff/max(ldif,dthr)
         if (f.gt.fd) then
            izd  = k
            fd   = f
         endif
      endif
      ! forward : dzden(k-1)/dzden(k)
      if (diff.gt.0.d0) then
         f = ldif/max(diff,dthr)
         if (f.gt.fd) then
            izd  = k-1
            fd   = f
         endif
      endif

      ! save diff
      ldif = diff
   enddo
endif


! find bin with max diff value fraction   dzval(k)/dzval(k+1)
if (method.eq.3 .or. method.eq.4) then
   ! threshold value val
   vthr = 0.d0
   do k=1,nz-2
      vthr = max(vthr,abs(zval(k+1) - zval(k)))
   enddo
   vthr = 0.5*vthr

   izv  = 1
   fv   = 0.d0
   ldif = abs(zval(2) - zval(1))
   do k=2,nz-2
      diff = abs(zval(k+1) - zval(k))
      ! backward: dzden(k)/dzden(k-1)
      if (ldif.gt.0.d0) then
         f = diff/max(ldif,vthr)
         if (f.gt.fv) then
            izv  = k
            fv   = f
         endif
      endif
      ! forward : dzden(k-1)/dzden(k)
      if (diff.gt.0.d0) then
         f = ldif/max(diff,vthr)
         if (f.gt.fv) then
            izv  = k-1
            fv   = f
         endif
      endif

      ! save diff
      ldif = diff
   enddo
endif


! combine methods
if (method.gt.1) then
   if (fv.gt.fd) then
      iz = izv
   else
      iz = izd
   endif
endif


! create empty value at position iz+1
do k=nz-1,iz+1,-1
   zval(k+1)=zval(k)
   zden(k+1)=zden(k)
enddo
iz = iz + 1


! calculate new z-value at centre of mass at bin iz
fct = zden(iz+1)+zden(iz-1)
if (fct.ne.0.d0) then
   fct = (2.*zden(iz+1)+zden(iz-1))/(3.d0*fct)
   fct = max(1./3.,min(fct,2./3.))
else
   fct = 1./2.
endif

zval(iz) = zval(iz-1) + fct*(zval(iz+1)-zval(iz-1))
zden(iz) = 0.d0


! end of program
return
end subroutine srvNxtZden

! ******************************************************************************

!> Find next z-value using srvNxtZden \n
!! Insert empty space in zx and zy
subroutine srvNxtZml(zval,zml,zx,zy,nz,iz,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-07-15      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nz          !< number of elements of variable Z   \n
                                                               !! On entry: number of known z-values \n
                                                               !!           is nz-1
double precision, dimension(nz), intent(inout)  :: zval        !< values of the PDF of variable Z
double precision, dimension(nz), intent(inout)  :: zml         !< maximum likelihood
double precision, dimension(nz), intent(inout)  :: zx          !< value for zml of marginal distribution X
double precision, dimension(nz), intent(inout)  :: zy          !< value for zml of marginal distribution Y

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK

integer                        , intent(out)    :: iz          !< position of the inserted z-value

! routine error code
integer,parameter    :: rtnerrcode = 140000000 +  84*1000  ! routine error code for routine srvNxtZml

! local variables
integer           :: k


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! find new value for z
call srvNxtZden(zval,zml,nz,iz,exitcode)


! create empty value at position iz+1
do k=nz-1,iz,-1
   zx(k+1)=zx(k)
   zy(k+1)=zy(k)
enddo


! end of program
return
end subroutine srvNxtZml

