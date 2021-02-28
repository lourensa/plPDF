
! calculate cumulative probabilities and probability densities for certain values

!> \file

!> Return the probability density of X at locations val
subroutine rvFunValDen(xval,xden,nx,val,vden,nv,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-07      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: nv          !< number of elements of val
double precision, dimension(nv), intent(in)     :: val         !< x-values of the PDF of X to return cum
double precision, dimension(nv), intent(out)    :: vden        !< probability densities of X at the locations val


integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: k,ix

double precision :: x,r,den

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all values of val
! ---------------------------
do k=1,nv
   x = val(k)
   
   if (x.lt.xval(1)) then
      den = 0.d0
   else if (x.gt.xval(nx)) then
      den = 0.d0
   else
      ! find bin of X containing x
      call srvFndBin(xval,nx,x,ix,.false.)
      ! calc den
      r   = (xden(ix+1)-xden(ix))/(xval(ix+1)-xval(ix))
      den = xden(ix) + r*(x-xval(ix))
   endif
   ! store result
   vden(k) = den
enddo


! end of program
return
end subroutine rvFunValDen

! ******************************************************************************

!> Return the cumulative distribution of X at locations val
subroutine rvFunValCum(xval,xden,nx,val,vcum,nv,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-07      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: nv          !< number of elements of val
double precision, dimension(nv), intent(in)     :: val         !< x-values of the PDF of X to return cum
double precision, dimension(nv), intent(out)    :: vcum        !< cumulative probabilities of X at locations val


integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: i,k,ix

double precision :: x,r,cum,den

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all values of val
! ---------------------------
do k=1,nv
   x = val(k)
   
   if (x.lt.xval(1)) then
      den = 0.d0
      cum = 0.d0
   else if (x.gt.xval(nx)) then
      den = 0.d0
      cum = 1.d0
   else
      ! find bin of X containing x
      call srvFndBin(xval,nx,x,ix,.false.)
      ! calc den
      r   = (xden(ix+1)-xden(ix))/(xval(ix+1)-xval(ix))
      den = xden(ix) + r*(x-xval(ix))
      ! calc cum
      cum = 0.d0
      do i=1,(ix-1)
         cum = cum + 0.5d0*(xden(i+1)+xden(i))*(xval(i+1)-xval(i))
      enddo
      cum = cum + 0.5d0*(den+xden(ix))*(x-xval(ix))
   endif
   ! store result
   vcum(k) = cum
enddo


! end of program
return
end subroutine rvFunValCum

! ******************************************************************************

!> Return the probability density and the cumulative distribution of X at locations val
subroutine rvFunValDenCum(xval,xden,nx,val,vden,vcum,nv,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-07      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: nv          !< number of elements of val
double precision, dimension(nv), intent(in)     :: val         !< x-values of the PDF of X to return cum
double precision, dimension(nv), intent(out)    :: vden        !< probability densities of X at the locations val
double precision, dimension(nv), intent(out)    :: vcum        !< cumulative probabilities of X at locations val


integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: i,k,ix

double precision :: x,r,cum,den

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all values of val
! ---------------------------
do k=1,nv
   x = val(k)
   
   if (x.lt.xval(1)) then
      den = 0.d0
      cum = 0.d0
   else if (x.gt.xval(nx)) then
      den = 0.d0
      cum = 1.d0
   else
      ! find bin of X containing x
      call srvFndBin(xval,nx,x,ix,.false.)
      ! calc den
      r   = (xden(ix+1)-xden(ix))/(xval(ix+1)-xval(ix))
      den = xden(ix) + r*(x-xval(ix))
      ! calc cum
      cum = 0.d0
      do i=1,(ix-1)
         cum = cum + 0.5d0*(xden(i+1)+xden(i))*(xval(i+1)-xval(i))
      enddo
      cum = cum + 0.5d0*(den+xden(ix))*(x-xval(ix))
   endif
   ! store result
   vden(k) = den
   vcum(k) = cum
enddo


! end of program
return
end subroutine rvFunValDenCum

! ******************************************************************************

!> Return x at the cumulative distribution of X given by vcum
subroutine rvFunCumVal(xval,xden,nx,val,vcum,nv,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2016-09-07      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

integer                        , intent(in)     :: nv          !< number of elements of val
double precision, dimension(nv), intent(out)    :: val         !< x-values of the PDF of X to return cum
double precision, dimension(nv), intent(in)     :: vcum        !< cumulative probabilities of X at locations val

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: i,k

double precision :: x,r,cum,bcum,tcum,dcum

logical          :: cont

double precision, parameter :: epsR   =1.e-3   ! minimum slope fraction of r to ignore it

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! loop over all values of val
! ---------------------------
do k=1,nv
   cum = vcum(k)
   
   if (cum.le.0.d0) then
      x = xval(1)
   else
      cont = .true.
      tcum = 0.d0
      i    = 1
      do while(cont)
         bcum = 0.5d0*(xden(i+1)+xden(i))*(xval(i+1)-xval(i))
         if (cum.gt.(tcum+bcum)) then
            ! cum not found, try next bin
            i = i + 1
            if (i.ge.nx) then
               ! no next bin anymore, use last value
               x = xval(nx)
               cont = .false.
            else
               ! next bin
               tcum = tcum + bcum
            endif
         else
            ! find x in the current bin i
            dcum = cum - tcum
            r   = (xden(i+1)-xden(i))
            ! decide or r is (nearly) 0 and may be ignored
            if(abs(r/(xden(i+1)+xden(i))) .lt. epsR) then
               ! ignore r, approximate x by a linear equation
               x = 2.d0*dcum/(xden(i+1)+xden(i))
               x = xval(i) + x
            else
               ! solve a quadratic equation: g(x) = 0.5*r*x^2 + xden[i]*x - dcum = 0
               r = r / (xval(i+1)-xval(i))
               ! for r<0: the lowest  value of x is the right solution
               ! for r>0: the highest value of x is the right solution
               !    This yields for both types of r the same formula.
               x = (-xden(i) + sqrt(xden(i)**2 + 2.d0*r*dcum)) / r
               x = xval(i) + x
            endif
            ! correct roundoff errors
            x = max(xval(i),min(x,xval(i+1)))
            ! ready
            cont = .false.
         endif
      enddo
   endif
   ! store result
   val(k)  = x
enddo


! end of program
return
end subroutine rvFunCumVal

