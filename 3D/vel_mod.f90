module vel_mod
use resample_mod
use modeling_mod
use sep_mod
!use randomMod
use sep
implicit none
real,private :: w0,qfact,gamma,pi
type(modelingT) :: initialS
type velT
  integer :: n1, n2, n3
  real    :: o1, o2, o3
  real    :: d1, d2, d3
  real, allocatable, dimension(:,:,:) :: dat,vgamma
end type
contains

subroutine readVel(fle,vel,full)
  character(len=*)  :: fle
  type(velT) :: vel
  type(modelingT) :: full
  character(len=128) :: junk
  integer :: ierr

  ierr= sep_get_data_axis_par(fle,1,vel%n1,vel%o1,vel%d1,junk)
  ierr= sep_get_data_axis_par(fle,2,vel%n2,vel%o2,vel%d2,junk)
  ierr= sep_get_data_axis_par(fle,3,vel%n3,vel%o3,vel%d3,junk)


  call from_param("qfact",qfact,50.)
  call from_param("w0",w0,60.)
  write(0, *), "vel (n1, n2, n3): (", vel%n1, vel%n2, vel%n3, ")"
  write(0, *), "qfact: ", qfact
  write(0, *), "w0: ", w0

  pi=atan(1.)*4.
  gamma=1./pi*atan(2*pi/qfact)
  allocate(vel%dat(vel%n1,vel%n2, vel%n3))
  ierr=sreed(fle,vel%dat,size(vel%dat)*4)
  !vel%dat=1500

  initialS%n1=vel%n1; initialS%o1=vel%o1; initialS%d1=vel%d1
  initialS%n2=vel%n2; initialS%o2=vel%o2; initialS%d2=vel%d2
  initialS%n3=vel%n3; initialS%o3=vel%o3; initialS%d3=vel%d3
  full%n1=vel%n1; full%o1=vel%o1; full%d1=vel%d1
  full%n2=vel%n2; full%o2=vel%o2; full%d2=vel%d2
  full%n3=vel%n3; full%o3=vel%o3; full%d3=vel%d3
end subroutine


!subroutine cleanVel(vel)
  !type(velT) :: vel
  !deallocate(vel%vgamma,vel%dat)
!end subroutine

!subroutine resampleVel(vel,cur,vuse)
  !type(velT) :: vel,vuse
  !type(modelingT) :: cur

  !deallocate(vuse%dat)

 !! call setVelSize(vel,cur,vuse)
  !deallocate(vuse%vgamma)
  !allocate(vuse%dat(cur%n1,cur%n2))
   !allocate(vuse%vgamma(cur%n1,cur%n2))

  !call interpField(initialS,cur,vel%dat,vuse%dat,.true.)

!!  call addRandom(vuse%dat,cur%bnd)
  !vuse%vgamma=-vuse%dat**(2*gamma-1)*w0**(2.*gamma)*sin(pi*gamma)/cur%dt
  !vuse%dat=vuse%dat*vuse%dat

!end subroutine

!subroutine writeFull(tag,cur,full,field,extend)
  !character(len=*) :: tag
  !type(modelingT) :: full,cur
  !logical :: extend
  !real :: field(:,:)
  !real, allocatable :: ff(:,:)
  !integer :: ierr
  !allocate(ff(full%n1,full%n2))

  !call interpField(cur,full,field,ff,extend)

  !ierr=srite(tag,ff,size(ff)*4)
  !deallocate(ff)

!end subroutine

!subroutine setVelSize(vel,cur,vuse)
  !type(modelingT) :: cur
  !type(velT) :: vuse,vel

  !vuse%o1=cur%o1-cur%bnd; vel%d1=cur%d1; vel%n1=cur%n1+2*cur%bnd
  !vuse%o2=cur%o2-cur%bnd; vel%d2=cur%d2; vel%n1=cur%n2+2*cur%bnd

!end subroutine


end module
