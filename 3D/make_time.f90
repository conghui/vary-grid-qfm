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
integer :: n1,n2,n3
real :: o1,o2,o3
real :: d1,d2,d3
integer :: ierr
integer :: j2,j3,i3,i2,n4out,n5out
character(len=128) :: junk
real,allocatable :: times(:,:,:),vel(:,:,:)
real :: s1,s2,s3,s0

call sep_init()
call from_param("j2",j2,1)
call from_param("j3",j3,1)
ierr=sep_get_data_axis_par("in",1,n1,o1,d1,junk)
ierr=sep_get_data_axis_par("in",2,n2,o2,d2,junk)
ierr=sep_get_data_axis_par("in",3,n3,o3,d3,junk)
n4out=ceiling(real(n2)/real(j2))
n5out=ceiling(real(n3)/real(j3))
ierr=sep_put_data_axis_par("out",4,n4out,o2,d2*j2,"x position")
ierr=sep_put_data_axis_par("out",5,n5out,o3,d3*j3,"y position")
allocate(vel(n1,n2,n3),times(n1,n2,n3))
ierr=sreed("in",vel,size(vel)*4)
vel=1./vel
s1=0
do i3=1,n3,j3
  do i2=1,n2,j2
    write(0,*) "working on shot ",(i2-1)/j2+1," of ",n4out
    s2=(i2-1)*d2
    s3=(i3-1)*d3
    s0=vel(1,i2,i3)
    call fastmarch(3,s1,s2,s3,1,1,1,n1,n2,n3,d1,d2,d3,s0,vel,times)
    ierr=srite("out",times,size(times)*4)
  end do
end do
end program
