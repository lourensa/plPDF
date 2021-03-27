
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

!> Find bin containing value v \n
!! Return the bin number ibin
subroutine srvFndBin(val,n,v,ibin,low) 

! History
! programmer------------date------------version---------------------------------
! loure001              2013-05-16      vcrc    initial version   


! declaration section
! ------------------------------------------------------------------------------

implicit none


! arguments
integer                        , intent(in)  :: n               !< number of values
double precision, dimension(n) , intent(in)  :: val             !< values


double precision               , intent(in)  :: v               !< search value

integer                        , intent(out) :: ibin            !< bin number for v      \n
                                                                !! result: 1 <= ibin < n

logical                        , intent(in)  :: low             !< when v coincides with a value of val: \n
                                                                !!   TRUE : the lower bin is chosen:  v in (lower,v    ]  \n
                                                                !!   FALSE: the upper bin is chosen:  v in [v    ,upper)

! local variables
integer          :: i,ilow,iupp


! program section
! ------------------------------------------------------------------------------

! binary search
ilow = 1
iupp = n

!    find ilow for: v in [val(ilow),val(ilow+1))  (left closed, right open)
do while ( (iupp-ilow).gt.1)
   i = int((iupp+ilow)/2)
   if (val(i).gt.v) then
      iupp = i
   else
      ilow = i
   endif
enddo

! check coincidence
! -----------------
!                       ilow      iupp
!                        |         |
! position of v:    1    2    3    4    5
! low=TRUE         -------0000000000+++++++
! low=FALSE        ------0000000000++++++++
!
!  -: ibin=ilow-1
!  0: ibin=ilow
!  +: ibin=ilow+1
!
if (v.gt.val(ilow)) then
   if (v.lt.val(iupp)) then
      ! position 3,                     v within found interval
      ibin=ilow                         ! v in (lower,upper)
   else
      if (v.gt.val(iupp)) then
         ! position 5,                  v greater than found interval
         ibin=ilow+1
      else
         ! position 4                   v equal to upper boundary
         if (low) then
            ibin=ilow                   ! v in (lower,v    ]
         else
            ibin=ilow+1                 ! v in [v    ,upper)
         endif
      endif
   endif
else
   if (v.lt.val(ilow)) then
      ! position 1,                     v less than found interval
      ibin=ilow-1
   else
      ! position 2                      v equal to lower boundary
      if (low) then
         ibin=ilow-1                    ! v in (lower,v    ]
      else
         ibin=ilow                      ! v in [v    ,upper)
      endif
   endif
endif


! 1 <= ibin < n
! -------------
ibin=min(max(1,ibin),n-1)


! end of program
return
end subroutine srvFndBin

