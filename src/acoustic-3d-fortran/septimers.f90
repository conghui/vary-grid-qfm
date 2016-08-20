module sep_timer_type_mod

type septimer
logical :: active
character(len=128) :: tag
real*8 :: tot_time
real*8 :: last_time
integer :: ncalls
end type
end module

module sep_timers_mod
use sep_timer_type_mod
type(septimer) :: timers(100)
integer,private,save :: ilast



contains

logical function init_sep_timers()
integer :: i
init_sep_timers=.false.
do i=1,100
  call clear_time(timers(i))
end do
init_sep_timers=.true.
ilast=0
end  function

subroutine print_timers()
integer :: i

do i=1,100
 call print_timer(timers(i))
end do

end subroutine

subroutine print_timer(st)
type(septimer) :: st
if(st%active) then
write(0,*) trim(st%tag), "=tag ncalls=",st%ncalls," tot_time=",st%tot_time
end if

end subroutine

subroutine stop_timer_num(inum)
call stop_timer(timers(inum))
end subroutine

subroutine stop_timer(st)
type(septimer) :: st
st%tot_time=st%tot_time+get_time()-st%last_time
st%ncalls=st%ncalls+1
end subroutine

subroutine start_timer_num(inum)
call start_timer(timers(inum))
end subroutine

subroutine start_timer(st)
type(septimer)::st
st%last_time=get_time()
end subroutine

logical function setup_cycle_timer(tag,inum)
character(len=*) :: tag
integer :: inum
setup_cycle_timer=setup_next_timer(tag,inum)
end function
logical function setup_next_timer(tag,inum)
integer :: inum
character(len=*) :: tag
setup_next_timer=.false.
if(ilast >99) then
  write(0,*) "out of timers"
  return
end if
ilast=ilast+1
call setup_timer_num(ilast,tag)
setup_next_timer=.true.
inum=ilast
end  function

subroutine setup_timer_num(inum,tag)
integer :: inum
character(len=*) :: tag
call setup_timer(timers(inum),tag)
end subroutine


subroutine setup_timer(st,tag)
type(septimer) :: st
character(len=*) :: tag
st%tag=tag
st%tot_time=0.
st%active=.true.
st%ncalls=0
end subroutine

subroutine clear_time(st)
type(septimer) :: st
st%tag=""
st%active=.false.
st%tot_time=0.
end subroutine



!From U of H (AGL Lab)
!From U of H (AGL Lab)
real*8 function get_time()
IMPLICIT NONE


CHARACTER (LEN=8) 	  :: date
CHARACTER (LEN=10) 	  :: time
CHARACTER (LEN=5) 	  :: zone
INTEGER, DIMENSION(8) 	  :: value
!_____________________________________________________________
! value(1)=year
! value(2)=month
! value(3)=day
! value(4)=difference from UTC in minutes
! value(5)=hour
! value(6)=minute
! value(7)=second
! value(8)=millisecond
!_____________________________________________________________
! calculate the current time on the cpu clock in s.
! use Fortran90 intrinsic function 'date_and_time'
!_____________________________________________________________
CALL DATE_AND_TIME(date,time,zone,value)

get_time=86400.*value(3)+3600.*value(5)+60.*value(6)+value(7)+.001*value(8)


END function

end module
