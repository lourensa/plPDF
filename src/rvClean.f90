
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

!> Clean an RV. Remove small errors probably caused by numerical inaccuracies
!!
!! switch options:
!!
!!  1. make 0 <= xcum <= 1
!!
!! 64. remove empty bins at outside
!!        xcum(1)==xcum(2)  ==0 & xden(1)==xden(2)  ==0
!!        xcum(n)==xcum(n-1)==1 & xden(n)==xden(n-1)==0
!!
!!  2. remove bin borders with no contribution: xcum(i-1)==xcum(i)==xcum(i+1)
!!
!!  4. truncate ends if PDF width > 2*nsigma*sigma
!!        The truncation is performed in such a way that the probability densities
!!        at both ends are equal (if possible)
!!        Further restriction: retaining probability must be at least minprob
!!                      
!!  8. extent end bins
!!        if xcum(ixb) > 0 extent bin ixb
!!        if xcum(ixe) < 1 extent bin ixe
!!
!! 16. rescale densities to let probability of density be close to cum-prob
!!
!! 32. renormalize the PDF
!!        set the total probability of the PDF to 1
!!        rescale the densities to estimate xcum as close as possible
!!     
!! Multiple switches may be specified through addition: e.g. switch = 1+8+32
!! If switch == 0, all options are used               : switch = 1+2+4+8+16+32+64
!! If switch < 0 , the specified switched are excluded: switch = 1+2+4+8+16+32+64 + switch
!!
subroutine rvClean(xval,xden,xcum,nx,nsigma,minprob,nscale,switch,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2015-07-07      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(inout)  :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(inout)  :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(inout)  :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
double precision, dimension(nx), intent(inout)  :: xcum        !< cumulative probabilities of the CDF of X

double precision               , intent(in)     :: nsigma      !< width of PDF to retain after truncation

double precision               , intent(in)     :: minprob     !< minimum probability to be retained
                                                               !! after truncation

integer                        , intent(in)     :: nscale      !< number of iterations to run the scale algorithm

integer                        , intent(in)     :: switch      !< Select which parts of the clean process have to be
                                                               !! executed.
                                                               !! The value is a sum of the wanted switch options.
                                                               !! If 0 , all options are executed.
                                                               !! If <0, the specified switches are excluded

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer           :: i,j,ixb,ixe,n
integer           :: usw,uns

logical           :: optXcum,optNoCont,optTrunc,optExtent,optRescale,optRenorm,optRmEmpty
logical           :: keep,shiftright,shiftleft

double precision  :: width,prob,v1,vn
double precision  :: cp,cw,lp,lw,f,sumden

double precision, dimension(2) :: stat


! program section
! ------------------------------------------------------------------------------


! init
exitcode = 0
ixb = 1     ! ixb begin position of data, may vary when bins are removed
ixe = nx    ! ixe end   position of data, may vary when bins are removed
uns = nscale
if (uns.lt.0) uns = 1
usw = switch


! find options
if (usw.le.0) usw = 1+2+4+8+16+32+64 + usw
optXcum    = (mod(int(usw/ 1),2).eq.1)   !  1.
optRmEmpty = (mod(int(usw/64),2).eq.1)   ! 64.
optNoCont  = (mod(int(usw/ 2),2).eq.1)   !  2.
optTrunc   = (mod(int(usw/ 4),2).eq.1)   !  4.
optExtent  = (mod(int(usw/ 8),2).eq.1)   !  8.
optRescale = (mod(int(usw/16),2).eq.1)   ! 16.
optRenorm  = (mod(int(usw/32),2).eq.1)   ! 32.


! calc statistics
call rvStat(xval,xden,nx,stat,exitcode)

! save first and last value of xval
v1 = xval(1)
vn = xval(nx)

!  1. make xcum >= 0 and <= 1
if (optXcum) then
   do i=ixb,ixe
      xcum(i) = max(0.d0,min(xcum(i),1.d0))
   enddo
endif


! 64. remove empty bins at outside
if (optRmEmpty) then
   do while(ixb.lt.ixe .and. xcum(ixb).eq.0.d0 .and. xcum(ixb+1).eq.0.d0  &
                       .and. xden(ixb).eq.0.d0 .and. xden(ixb+1).eq.0.d0)
      ixb = ixb + 1
   enddo
   do while(ixb.lt.ixe .and. xcum(ixe).eq.1.d0 .and. xcum(ixe-1).eq.1.d0  &
                       .and. xden(ixe).eq.0.d0 .and. xden(ixe-1).eq.0.d0)
      ixe = ixe - 1
   enddo
endif


!  2. remove bin borders with no contribution
if (optNoCont) then
   j = ixb
   do i=ixb+1,ixe-1
      if (xcum(j).ne.xcum(i+1)) then
         keep = .true.
      else
         if (xcum(j).ne.xcum(i)) then
            keep = .true.
         else
            keep = .false.
         endif
      endif
      if (keep) then
         ! keep bin border i
         j = j + 1
         xval(j) = xval(i)
         xden(j) = xden(i)
         xcum(j) = xcum(i)
      endif
   enddo
   j = j + 1
   xval(j) = xval(ixe)
   xden(j) = xden(ixe)
   xcum(j) = xcum(ixe)
   ixe = j
endif


!  4. truncate
if (optTrunc) then
   width = 2.0d0*nsigma*stat(2)
   call srvTrunc(xval,xden,xcum,nx,width,minprob,ixb,ixe,exitcode)
endif


!  8. extent outer bins
if (optExtent) then

   ! try to shift data to create an empty space at the ends where necessary
   n = nx - (ixe-ixb+1)  ! number of unused values
   shiftright = .false.
   shiftleft  = .false.
   if (xcum(ixb).gt.0.d0) then
      if (xcum(ixe).lt.1.d0) then
         ! extend at both sides
         if (n.ge.2) then
            ! enough space
            if (ixb.eq.1) then
               ! shift right
               shiftright = .true.
            else if (ixe.eq.nx) then
               ! shift left
               shiftleft  = .true.
            endif
!         else
!            ! not enough space to extend both sides
!            ! leave unchanged
         endif
      else
         ! extend only lower side
         if (n.ge.1) then
            if (ixb.eq.1) shiftright = .true.
         endif
      endif
   else
      if (xcum(ixe).lt.1.d0) then
         ! extent only upper sides
         if (n.ge.1) then
            if (ixe.eq.nx) shiftleft = .true.
         endif
!      else
!         ! extend none
!         ! leave unchanged
      endif
   endif
   ! shift
   if (shiftright) then
      ! shift one position to the right
      j = ixe + 1
      do i=ixe,ixb,-1
         xval(j) = xval(i)
         xden(j) = xden(i)
         xcum(j) = xcum(i)
         j = j - 1
      enddo
   else if (shiftleft) then
      ! shift one position to the left
      j = ixb - 1
      do i=ixb,ixe
         xval(j) = xval(i)
         xden(j) = xden(i)
         xcum(j) = xcum(i)
         j = j + 1
      enddo
   endif

   ! first bin
   if (xcum(ixb).gt.0.d0) then
      if (ixb.eq.1) then
         ! no space to create extra bin, increase probability of current bin
         ixb = ixb + 1
      endif

      ! create extra bin
      if (xval(ixb).gt.v1) then
         ! use original least value (v1) as end of PDF
         ! Unless v1 is too low
         if (xden(ixb).gt.0.d0) then
            xval(ixb-1) = xval(ixb) - 2.d0*xcum(ixb)/xden(ixb)
            xval(ixb-1) = max(xval(ixb-1),v1)
         else
            xval(ixb-1) = v1
         endif
         xden(ixb-1) = 2.d0*xcum(ixb)/(xval(ixb)-xval(ixb-1)) - xden(ixb)
         xcum(ixb-1) = 0.d0
         ixb = ixb - 1
      else
         if (xden(ixb).gt.0.d0) then
            ! create new value
            xval(ixb-1) = xval(ixb) - 2.d0*xcum(ixb)/xden(ixb)
            xden(ixb-1) = 0.d0
            xcum(ixb-1) = 0.d0
            ixb = ixb - 1
         else
            ! very arbitrary solution for xval(ixb-1)
            ! * last bin same width as second last bin
            xval(ixb-1) = 2.d0*xval(ixb) - xval(ixb+1)
            xden(ixb-1) = 2.d0*xcum(ixb)/(xval(ixb)-xval(ixb-1)) - xden(ixb)
            xcum(ixb-1) = 0.d0
            ixb = ixb - 1
         endif
      endif
   endif

   ! last bin
   if (xcum(ixe).lt.1.d0) then
      if (ixe.eq.nx) then
         ! no space to create extra bin, increase probability of current bin
         ixe = ixe - 1
      endif

      ! create extra bin
      if (xval(ixe).lt.vn) then
         ! use original least value (vn) as end of PDF
         ! Unless v1 is too low
         if (xden(ixe).gt.0.d0) then
            xval(ixe+1) = xval(ixe) + 2.d0*(1.d0-xcum(ixe))/xden(ixe)
            xval(ixe+1) = min(xval(ixe+1),vn)
         else
            xval(ixe+1) = vn
         endif
         xden(ixe+1) = 2.d0*(1.d0-xcum(ixe))/(xval(ixe+1)-xval(ixe)) - xden(ixe)
         xcum(ixe+1) = 1.d0
         ixe = ixe + 1
      else
         if (xden(ixe).gt.0.d0) then
            ! create new value
            xval(ixe+1) = xval(ixe) + 2.d0*(1.d0-xcum(ixe))/xden(ixe)
            xden(ixe+1) = 0.d0
            xcum(ixe+1) = 1.d0
            ixe = ixe + 1
         else
            ! very arbitrary solution for xval(ixe+1)
            ! * last bin same width as second last bin
            xval(ixe+1) = 2.d0*xval(ixe) - xval(ixe-1)
            xden(ixe+1) = 2.d0*(1.d0-xcum(ixe))/(xval(ixe+1)-xval(ixe)) - xden(ixe)
            xcum(ixe+1) = 1.d0
            ixe = ixe + 1
         endif
      endif
   endif

endif


! 16. rescale densities
if (optRescale) then
   ! remove discrepancy in xcum and xden probability for each bin
   do j=1,uns
      lw = 0.d0 ! width of former bin
      lp = 0.d0 ! probability discrepancy of former bin, to be processed in current bin
      do i=ixb,ixe-1
         ! current bin width
         cw = xval(i+1)-xval(i)
         ! error between density and cumulative probability
         sumden = xden(i+1)+xden(i)
         cp = xcum(i+1)-xcum(i) - cw*0.5d0*sumden
         if (sumden.gt.0.d0) then
            ! error is devided proportional to both border densities
            f = xden(i)/sumden
         else
            ! densities are 0, arbitrary choice of 0.5
            f = 0.5d0
         endif

         ! correct in left bin border
         lp = lp + f*cp    ! amount of error in probability for left border
         xden(i) = xden(i) + lp/(lw+cw)
         ! to avoid negative densities
         xden(i) = max(0.d0,xden(i))

         ! save values for the right border
         lw = cw
         lp = (1.d0 - f)*cp
      enddo

      ! correct last bin border
      i = ixe
      xden(i) = xden(i) + lp/(lw)
      ! to avoid negative densities
      xden(i) = max(0.d0,xden(i))
   enddo

endif


! 32. renormalize probability of PDF to 1
if (optRenorm) then
   ! calculate total probability of PDF
   call rvStatProb(xval(ixb),xden(ixb),ixe-ixb+1,prob,exitcode)

   ! renormalize densities
   prob = 1.d0/prob  ! multiplication is cheaper than division
   do i=ixb,ixe
      xden(i) = xden(i)*prob
!      xcum(i) = xcum(i)*prob
   enddo
endif


! move data to begin position
if (ixb.gt.1) then
   j = 1
   do i=ixb,ixe
      xval(j) = xval(i)
      xden(j) = xden(i)
      xcum(j) = xcum(i)
      j = j + 1
   enddo
endif
nx = ixe - ixb + 1


! end of program
return
end subroutine rvClean

! *********************
! *********************
!> version 3 of truncation algorithm
!! The algorithm is separated into nine main parts based on the slopes of the PDFs: rxb and rxe
!! These values can be: >0, ==0, <0
!! 
subroutine srvTrunc(xval,xden,xcum,nx,width,minprob,ixb,ixe,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2015-08-21      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(inout)  :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(inout)  :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(inout)  :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X
double precision, dimension(nx), intent(inout)  :: xcum        !< cumulative probabilities of the CDF of X

double precision               , intent(in)     :: width       !< width of PDF to retain after truncation

double precision               , intent(in)     :: minprob     !< minimum probability to be retained
                                                               !! after truncation

integer                        , intent(inout)  :: ixb         !< begin position of remaining PDF
integer                        , intent(inout)  :: ixe         !< end   position of remaining PDF

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer           :: ib,ie,truncType

logical           :: cont,trunc,parallel,pdown

double precision  :: prob,xb,xe,rxb,rxe
double precision  :: a,b,c,D,den,denb,dene,dxb,dxe,dx,eps,f

! functions
double precision  :: srvInqIntCum,srvInqIntDen

! program section
! ------------------------------------------------------------------------------


! init
exitcode = 0


! truncate
cont  = ((xval(ixe)-xval(ixb)).gt.width)                  ! current width > wanted width
cont  = (cont .and. ((ixe-ixb).gt.1))                     ! at least 2 bins (3 bin borders)
cont  = (cont .and. ((xcum(ixe) - xcum(ixb)).gt.minprob) )! at least enough probability in current PDF
do while(cont)
   ! entring this loop means
   ! - enough probability in the remaining PDF  (xcum(ixe)-xcum(ixb)) > minprob
   ! - enough width in the remaining PDF        () > width
   ! convinient indices
   ib = ixb
   ie = ixe - 1

   ! slopes of the PL-PDFs
   rxb  = (xden(ib+1)-xden(ib))/(xval(ib+1)-xval(ib)) ! assume non-zero bin width
   rxe  = (xden(ie+1)-xden(ie))/(xval(ie+1)-xval(ie)) ! assume non-zero bin width

   ! test or lines are parallel
   parallel = (.not. abs(rxb-rxe).gt.0.d0)


   ! perform for each of the nine combinations
   ! calc densities
   denb= xden(ib) + rxb*(xb-xval(ib))
   dene= xden(ie) + rxe*(xe-xval(ie))
   trunc = .false.
   pdown = .true.   ! extend the width (xe-xb) in the direction of downwards density
   if (rxb.gt.0.d0) then
      if (rxe.gt.0.d0) then
         ! + +
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie+1)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib+1).lt.xden(ie)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      else if (rxe.lt.0.d0) then
         ! + -
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib+1).lt.xden(ie+1)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      else
         ! + 0
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib+1).lt.xden(ie+1)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      endif
   else if (rxb.lt.0.d0) then
      if (rxe.gt.0.d0) then
         ! - +
         !=====
         ! find truncation type
         if (xden(ib).lt.xden(ie)) then
            ! 1. only trunc bin ie (higher bin!)
            truncType = 1
         else if (xden(ib+1).gt.xden(ie+1)) then
            ! 2. only trunc bin ib (higher bin!)
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
            pdown = .false.
         endif
      else if (rxe.lt.0.d0) then
         ! - -
         !=====
         ! find truncation type
         if (xden(ib+1).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib).lt.xden(ie+1)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      else
         ! - 0
         !=====
         ! find truncation type
         if (xden(ib+1).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib).lt.xden(ie)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      endif
   else
      if (rxe.gt.0.d0) then
         ! 0 +
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie+1)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib).lt.xden(ie)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      else if (rxe.lt.0.d0) then
         ! 0 -
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib).lt.xden(ie+1)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 3. trunc both
            truncType = 3
         endif
      else
         ! 0 0
         !=====
         ! find truncation type
         if (xden(ib).gt.xden(ie)) then
            ! 1. only trunc bin ie
            truncType = 1
         else if (xden(ib).lt.xden(ie)) then
            ! 2. only trunc bin ib
            truncType = 2
         else
            ! 4. trunc both (special situation, horizontal lines at the same level)
            truncType = 4
         endif
      endif
   endif

   ! perform truncation
   trunc = .false.
   if (truncType.eq.1) then
!exitcode=exitcode+32**0
      ! 1. only trunc bin ie
      xb = xval(ib)
      xe = xb + width
      xe = max(xval(ie),min(xe,xval(ie+1)))

      ! test prob
      prob = srvInqIntCum(xval,xden,xcum,nx,ie,xe) - &
             srvInqIntCum(xval,xden,xcum,nx,ib,xb)
      prob = minprob - prob   ! remaining probability
      if (prob.gt.0.d0) then
         ! not enough probability, extent xe
         den = xden(ie) + rxe*(xe-xval(ie))
         if (abs(rxe).gt.0.d0) then
            D   = den**2 + 2.d0*rxe*prob
            D   = max(0.d0,D)  ! avoid roundof errors
            dx  = (-den+sqrt(D))/rxe
            ! select right root
            if (dx.gt.(xval(ie+1)-xe) .or. dx.lt.0.d0) then
               dx  = (-den-sqrt(D))/rxe
            endif
         else
            dx = prob/den
         endif
         ! new xe
         xe = xe + dx
         trunc = .true.
      else
         ! check width only
         if (xe.eq.xval(ie)) then
            ! still too wide
            ixe = ixe - 1
         else
            ! truncate
            trunc = .true.
         endif
      endif
   else if (truncType.eq.2) then
!exitcode=exitcode+32**1
      ! 2. only trunc bin ib
      xe = xval(ie+1)
      xb = xe - width
      xb = max(xval(ib),min(xb,xval(ib+1)))

      ! test prob
      prob = srvInqIntCum(xval,xden,xcum,nx,ie,xe) - &
             srvInqIntCum(xval,xden,xcum,nx,ib,xb)
      prob = minprob - prob   ! remaining probability
      if (prob.gt.0.d0) then
         ! not enough probability, extent xb
         den = xden(ib) + rxb*(xb-xval(ib))
         if (abs(rxe).gt.0.d0) then
            D   = den**2 - 2.d0*rxb*prob
            D   = max(0.d0,D)  ! avoid roundof errors
            dx  = (-den+sqrt(D))/(-rxb)
            ! select right root
            if (dx.gt.(xb-xval(ib)) .or. dx.lt.0.d0) then
               dx  = (-den-sqrt(D))/(-rxb)
            endif
         else
            dx = prob/den
         endif
         ! new xb
         xb = xb - dx
         trunc = .true.
      else
         ! check width only
         if (xb.eq.xval(ib+1)) then
            ! still too wide
            ixb = ixb + 1
         else
            ! truncate
            trunc = .true.
         endif
      endif
   else if (truncType.eq.3) then
!exitcode=exitcode+32**2
      ! 3. trunc both

      ! calc xb and xe for width, only for non-parallel lines
      if (.not.parallel) then
         ! find xb and xe for (xe-xb)==width and p(xb)==p(xe)
         ! The found values are not necessarily within their respective bin borders
         xb = (xden(ib) - rxb*xval(ib) - xden(ie) - rxe*(width - xval(ie))) / (rxe-rxb)
      else
         ! set selection to highest probability
         xb = xval(ib)
         xe = xb + width
         ! calc densities
         denb= xden(ib) + rxb*(xb-xval(ib))
         dene= xden(ie) + rxe*(xe-xval(ie))
         ! get largest probability xb value
         if (rxb.gt.0.d0) then
            if (denb.gt.dene) then
               xb = xval(ib)
            else
               xb = xval(ib+1)
            endif
         else
            if (denb.lt.dene) then
               xb = xval(ib)
            else
               xb = xval(ib+1)
            endif
         endif
      endif
      ! set x values within their bin borders
      ! hereafter, xe-xb >= width
      xb = max(xval(ib),min(xb,xval(ib+1)))
      xe = xb + width
      xe = max(xval(ie),min(xe,xval(ie+1)))
      xb = xe - width
      eps= 4.*abs((xe-xb)-width)             ! roundoff errors
      xb = max(xval(ib),min(xb,xval(ib+1)))

      ! check probability
      prob = srvInqIntCum(xval,xden,xcum,nx,ie,xe) - &
             srvInqIntCum(xval,xden,xcum,nx,ib,xb)
      prob = minprob - prob   ! remaining probability

      if (prob.gt.0.d0) then
         ! extend area at both sides to meet the probability constraint too

         ! calc densities
         denb = xden(ib) + rxb*(xb-xval(ib))
         dene = xden(ie) + rxe*(xe-xval(ie))

         ! calc dx to let one PDF meet the density of the other
         ! dxb is negative
         if (abs(rxb).gt.0.d0) then
            dxb = (dene - denb)/rxb
            dxb = max(dxb,xval(ib)-xb)
            dxb = min(0.d0,dxb)
         else
            dxb = 0.d0
         endif
         if (dxb.lt.0.d0) then
            ! check dxb for prob 
            a = 0.5d0*rxb
            b = denb
            c = prob
            D = b**2 - 4.d0*a*c
            if (D.ge.0.d0) then
               D   = sqrt(D)
               dxb = max(dxb,min((-b-D)/(2.d0*a),(-b+D)/(2.d0*a)))
            endif

            ! new xb
            prob = prob + dxb*(denb+0.5d0*rxb*dxb)
            xb   = xb + dxb
            denb = xden(ib) + rxb*(xb-xval(ib))
!exitcode=exitcode+32**4
         endif

         if (prob.gt.0.d0) then
            ! dxe is positive
            if (abs(rxe).gt.0.d0) then
               dxe = (denb - dene)/rxe
               dxe = min(dxe,xval(ie+1)-xe)
               dxe = max(0.d0,dxe)
            else
               dxe = 0.d0
            endif
            if (dxe.gt.0.d0) then
               ! check dxe for prob 
               a = 0.5d0*rxe
               b = dene
               c = -prob
               D = b**2 - 4.d0*a*c
               if (D.ge.0.d0) then
                  D   = sqrt(D)
                  dxe = min(dxe,max((-b-D)/(2.d0*a),(-b+D)/(2.d0*a)))
               endif

               ! new xe
               prob = prob - dxe*(dene+0.5d0*rxe*dxe)
               xe   = xe + dxe
               dene = xden(ie) + rxe*(xe-xval(ie))
!exitcode=exitcode+32**5
            endif
         endif

         if (prob.gt.0.d0) then
            ! still some probability left
            ! extend both ends simultaniously
            ! The extent will be done proportional to the remaining distances at both ends.
            ! This is an arbitrary choice.

            dxb = max(0.d0,xb - xval(ib))
            dxe = max(0.d0,xval(ie+1) - xe)

            ! calc f
            ! solve: prob = f*dxb*(denb-0.5*rxb*f*dxb) + f*dxe*(dene+0.5*rxe*f*dxe)
            a = 0.5d0*(rxe*dxe**2 - rxb*dxb**2)
            b = dxb*denb + dxe*dene
            c = -prob
            if (abs(a).gt.0.d0) then
               D = b**2-4.d0*a*c
               if (D.gt.0.d0) then
                  D = sqrt(D)
                  f = (-b-D)/(2.d0*a)
                  if (f.gt.1.d0 .or. f.lt.0.d0) then
                     ! wrong solution, try the other one
                     f = (-b+D)/(2.d0*a)
                     if (f.gt.1.d0 .or. f.lt.0.d0) then
                        ! no solution
                        f = -1.d0
                     endif
                  endif
               else
                  ! no solution
                  f = -1.d0
               endif
            else
               if (abs(b).gt.0.d0) then
                  f = -c/b
               else
                  ! no solution
                  f = -1.d0
               endif
            endif
            if (f.gt.0.d0) then
               xb = xb - f*dxb
               xe = xe + f*dxe
               xb = max(xval(ib),min(xb,xval(ib+1)))
               xe = max(xval(ie),min(xe,xval(ie+1)))
!exitcode=exitcode+32**6
            endif
         endif

         trunc = .true.
      else
         ! remove a bin?
         if (abs(xe-xb-width).gt.eps) then
            ! remove one bin, save bin with largest remaining probability at the end
            if (xcum(ib+1).lt.(1.d0-xcum(ie))) then
               ixb = ixb + 1
            else
               ixe = ixe - 1
            endif
         else
            ! trunc
            trunc = .true.
         endif
      endif

   else
!exitcode=exitcode+32**3
      ! 4. trunc both (special situation, horizontal lines at the same level)

      ! truncate the same distance at both sides
      dx = (xval(ie+1)-xval(ib) - width)*0.5d0
      xb = xval(ib)+dx

      xb = max(xval(ib),min(xb,xval(ib+1)))
      xe = xb + width
      xe = max(xval(ie),min(xe,xval(ie+1)))
      xb = xe - width
      eps= 4.*abs((xe-xb)-width)             ! roundoff errors
      xb = max(xval(ib),min(xb,xval(ib+1)))

      ! check probaility
      prob = srvInqIntCum(xval,xden,xcum,nx,ie,xe) - &
             srvInqIntCum(xval,xden,xcum,nx,ib,xb)
      prob = minprob - prob   ! remaining probability

      if (prob.gt.0.d0) then
         if (xden(ib).gt.0.d0) then
            ! devide probability equaly over the two bins
            dx = 0.5d0*prob/xden(ib)
            xb = max(xval(ib)  ,xb - dx)
            xe = min(xval(ie+1),xe + dx)
            trunc = .true.
         else
            ! strange situation, probably impossible
            ! - density zero, probability too low
            trunc = .true.
         endif
      else
         ! remove a bin?
         if (abs(xe-xb-width).gt.eps) then
            ! remove one bin, save bin with largest remaining probability at the end
            if (xcum(ib+1).lt.(1.d0-xcum(ie))) then
               ixb = ixb + 1
            else
               ixe = ixe - 1
            endif
         else
            ! trunc
            trunc = .true.
         endif
      endif

   endif


   ! truncate bin ixb and ixe
   if (trunc) then
      ! fill new boundary values
      ! if xb == xval(ib) nothing has to be changed
      if(xb.gt.xval(ixb)) then
         if (xb.lt.xval(ixb+1)) then
            ! create new lower boundary
            xcum(ixb) = srvInqIntCum(xval,xden,xcum,nx,ixb,xb)
            xcum(ixb) = max(0.d0,min(xcum(ixb),xcum(ixb+1)))
            xden(ixb) = srvInqIntDen(xval,xden,nx,ixb,xb)
            xval(ixb) = xb
         else
            ! change to upper boundary
            ixb = ixb+1
         endif
      endif

      ! if xe == xval(ie+1) nothing has to be changed
      if(xe.lt.xval(ixe)) then
         if (xe.gt.xval(ixe-1)) then
            ! create new upper boundary
            xcum(ixe) = srvInqIntCum(xval,xden,xcum,nx,ixe-1,xe)
            xcum(ixe) = max(xcum(ixe-1),min(xcum(ixe),1.d0))
            xden(ixe) = srvInqIntDen(xval,xden,nx,ixe-1,xe)
            xval(ixe) = xe
         else
            ! change to lower boundary
            ixe = ixe-1
         endif
      endif

      cont = .false.
   endif

   ! test for next iteration
   cont  = (cont .and. ((xval(ixe)-xval(ixb)).gt.width))     ! current width > wanted width
   cont  = (cont .and. ((ixe-ixb).gt.1))                     ! at least 2 bins (3 bin borders)
   cont  = (cont .and. ((xcum(ixe) - xcum(ixb)).gt.minprob) )! at least enough probability in current PDF
enddo


! end of program
return
end subroutine srvTrunc



