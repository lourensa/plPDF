
!> \file 

!> Calculate the statistics of a piecewise linear PDF (mean value)
subroutine rvStatMean(xval,xden,nx,mean,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-06-04      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X
                                                               !! \n xval and xden together form the PDF of X

double precision               , intent(out)    :: mean        !< mean value of PDF

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer        :: i


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0
mean     = 0.d0

! calc
do i=1,nx-1
   mean = mean + (xval(i+1)-xval(i)) * ( (xden(i)  *(     xval(i+1)+ 2.d0*xval(i) ))   &
                                        +(xden(i+1)*(2.d0*xval(i+1)+      xval(i) )) )
enddo
mean=mean/6.d0


! end of program
return
end subroutine rvStatMean

! ******************************************************************************

!> Calculate the statistics of a piecewise linear PDF (mean and standard deviation)
subroutine rvStat(xval,xden,nx,stat,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2014-06-04      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X

double precision, dimension(2) , intent(out)    :: stat        !< mean and standard deviation of PDF

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer        :: i

double precision :: mean,variance,dx,dd,x1,x2,d1,d2,f


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0

! calc mean
call rvStatMean(xval,xden,nx,mean,exitcode)


! calc variance
variance = 0.d0

f = 1.d0/4.d0
variance = variance + (xden(nx)*(xval(nx)-mean)**3)*f
variance = variance - (xden(1) *(xval(1) -mean)**3)*f

f = 1.d0/12.d0
do i=1,nx-1
   x1 = xval(i)
   x2 = xval(i+1)
   d1 = xden(i)
   d2 = xden(i+1)

   dx = x2 - x1
   dd = d2 - d1

   variance = variance + (d1*dx + dd*(mean-x1))*(dx**2 + 3.d0*(x2-mean)*(x1-mean))*f

enddo


! set values
stat(1) = mean
stat(2) = sqrt(variance)


! end of program
return
end subroutine rvStat

! ******************************************************************************

!> Find the index value of the largest density of a piecewise linear PDF
subroutine rvStatImax(xval,xden,nx,imax,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-07-15      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X

integer                        , intent(out)    :: imax        !< position of maximum density value

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer        :: i


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! calc
imax = 1
do i=2,nx
   if (xden(i).gt.xden(imax)) imax=i
enddo


! end of program
return
end subroutine rvStatImax

! ******************************************************************************

!> Calculate the total probability of a piecewise linear PDF
subroutine rvStatProb(xval,xden,nx,prob,exitcode) 

! History
! programmer------------date------------version---------------------------------
! loure001              2014-08-13      initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X

double precision               , intent(out)    :: prob        !< total probability of a PDF

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer        :: i


! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! calc
prob = 0.d0
do i=2,nx
   prob = prob + (xden(i-1)+xden(i))*(xval(i)-xval(i-1))
enddo
prob = prob/2.d0


! end of program
return
end subroutine rvStatProb

! ******************************************************************************

!> Calculate the total area of the absolute value of a piecewise linear function.
!! This routine honours intersections at the horizontal axis.
subroutine rvStatAbsArea(xval,xden,nx,prob,exitcode) 

! History
! programmer------------date------------version---------------------------------
! Aris                  2019-11-16      initial version, started as a copy of rvStatProb


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)     :: nx          !< number of elements of variable X
double precision, dimension(nx), intent(in)     :: xval        !< values of the PDF of variable X
double precision, dimension(nx), intent(in)     :: xden        !< probability densities of the PDF of variable X

double precision               , intent(out)    :: prob        !< total probability of a PDF

integer                        , intent(out)    :: exitcode    !< exit status of routine, 0=OK


! local variables
integer          :: i
double precision :: d1,d2,xs,x1,x2

! program section
! ------------------------------------------------------------------------------

! init
exitcode = 0


! calc
prob = 0.d0
do i=2,nx
   d1 = xden(i-1)
   d2 = xden(i)
   x1 = xval(i-1)
   x2 = xval(i)
   ! if d1*d1 is negative then the horizontal axis is crossed
   if(.not. ((d1*d2).lt.0.d0)) then
      ! both densities are:
      !   - positive or zero
      !   - negative or zero
      prob = prob + abs((d1+d2))*(x2-x1)
   else
      ! one density positive and one density negative
      ! find x-value for intersection
      xs = (abs(d1)*x2 + abs(d2)*x1)/(abs(d1)+abs(d2))
      ! area
      prob = prob + abs(d1)*(xs-x1) + abs(d2)*(x2-xs)
   endif
enddo
prob = prob/2.d0


! end of program
return
end subroutine rvStatAbsArea

