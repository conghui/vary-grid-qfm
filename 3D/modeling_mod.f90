module modeling_mod
implicit none
type modelingT
  real,pointer :: sincT(:,:) ! ASK: what's sincT mean
  integer :: oversamp,oversampS
  real :: o1,o2,o3
  integer :: n1,n2,n3
  real :: d1,d2,d3
  integer :: bnd
  integer,allocatable :: b(:,:), e(:,:)
  integer :: b1,b2,b3,nsect
  double precision :: d1a(5),d2a(5),d3a(5)
  real :: d0,dt,dtextra
  integer :: ntblock
  integer :: jimage
end type


contains

logical function checkSame(a,b) result(res)
  type(modelingT) :: a,b
  real :: d

  res=.false.
  if(a%n1/=b%n1 .or. a%n2/=b%n2 .or. a%n3/=b%n3) then
     return
  end if
  if(abs(a%d1/b%d1-1.)>.01) then
    return
  end if
  if(abs(a%d2/b%d2-1.)>.01) return
  if(abs(a%d3/b%d3-1.)>.01) return

  res=.true.
end function

subroutine clean_modeling(modelingS)
  type(modelingT) :: modelingS
   deallocate(modelingS%b,modelingS%e)
end subroutine


end module
