program ab
use sep
integer :: n1,n2
real, allocatable :: in(:,:),out(:,:)
integer :: i1,ierr

call sep_init()
call from_history("n2",n2)
call from_history("n1",n1)
allocate(in(n1,n2),out(n1,n2))
ierr=sreed("in",in,size(in)*4)
out(1:n1-1,:)=in(2:,:)-in(1:n1-1,:)
out(n1,:)=out(n1-1,:)


ierr=srite("out",out,size(in)*4)

end program
