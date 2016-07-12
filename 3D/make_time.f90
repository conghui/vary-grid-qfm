module times_mod
use iso_c_binding
interface
   subroutine fastmarch(order,s1,s2,s3,b1,b2,b3,n1,n2,n3,d1,d2,d3,s0,slow,ttime) bind(c,name="fastmarch")
     import
     real(c_float),intent(in),value:: s1,s2,s3,d1,d2,d3,s0
     integer(c_int), value,intent(in) :: order,b1,b2,b3,n1,n2,n3
     real(c_float),dimension(*) :: slow,ttime
    end subroutine
 end interface
end module

program ttimes
use times_mod
use sep
implicit none
integer :: n1,n2
real :: o1,o2
real :: d1,d2
integer :: ierr
integer :: j2,i2,n3out
character(len=128) :: junk
real,allocatable :: times(:,:),vel(:,:)
real :: s1,s2,s0

call sep_init()
call from_param("j2",j2,1)
ierr=sep_get_data_axis_par("in",1,n1,o1,d1,junk)
ierr=sep_get_data_axis_par("in",2,n2,o2,d2,junk)
n3out=ceiling(real(n2)/real(j2))
ierr=sep_put_data_axis_par("out",3,n3out,o2,d2*j2,"x position")
allocate(vel(n1,n2),times(n1,n2))
ierr=sreed("in",vel,size(vel)*4)
vel=1./vel
s1=0
do i2=1,n2,j2
  write(0,*) "working on shot ",(i2-1)/j2+1," of ",n3out
  s2=(i2-1)*d2
  s0=vel(1,i2)
  call fastmarch(2,s1,s2,0,1,1,1,n1,n2,1,d1,d2,d2,s0,vel,times)
  ierr=srite("out",times,size(times)*4)
end do
end program
