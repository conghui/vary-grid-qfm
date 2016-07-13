module source_mod
  use sep_mod
  use sep
  use modeling_mod
  implicit none
  type sourceT
    real :: x,z,y
    real  :: dt
    real,allocatable :: dat(:)
  end type
  type Wavefield
    real,allocatable :: ar(:,:,:)
  end type
  type sourceW
     type(Wavefield), dimension(:),allocatable :: W
  end type
  contains


  !subroutine  initSourceSave(source,nb)
    !type(sourceW) :: source
    !integer :: nb
    !allocate(source%W(nb))
 !end subroutine

 !subroutine cleanSourceSave(source)
   !type(sourceW) :: source
   !integer :: i
   !do i=1,size(source%W)
     !deallocate(source%w(i)%ar)
   !end do
   !deallocate(source%w)
  !end subroutine

  subroutine setSourceLocation(source,loc)
     type(sourceT) :: source
     real :: loc(:)

     source%x=loc(1)
     source%y = loc(2)
  end subroutine


  subroutine readSource(tag,src)
   character(len=*) :: tag
   type(sourceT) :: src
   integer :: nt
   real :: dt, ot
   real, allocatable :: tmp(:)
   character(len=128) :: tmpS
   integer :: ierr


   ierr= sep_get_data_axis_par(tag,1,nt,ot,dt,tmpS)

   call from_param("src_z",src%z)
   call from_param("src_y", src%y)

   src%dt=dt

   allocate(tmp(nt))
   ierr=sreed(tag,tmp,size(tmp)*4)
   allocate(src%dat(8+nt))
   src%dat(1:4)=0
   src%dat(5:4+nt)=tmp
   src%dat(nt-3:)=0
   deallocate(tmp)

 end subroutine


 !subroutine addSource(src,tm,rnew,cur,scale)
   !type(sourceT) :: src
   !real :: tm
   !integer :: it,ix,iy,iz,ibase,isinc
   !real :: rnew(:,:)
   !real :: scale
   !type(modelingT) :: cur
   !real,pointer :: sincT(:)
   !it=tm/src%dt
   !isinc=(tm-it*src%dt)/src%dt*10000+1
   !sincT=>cur%sincT(:,isinc)
   !ibase=it+4


  !if(ibase+7>size(src%dat)) return
     !!for now nearest neighbor
     !iz=(src%z-cur%o1)/cur%d1+1.5
     !ix=(src%x-cur%o2)/cur%d2+1.5
          !rnew(iz,ix)=rnew(iz,ix)+scale*sum(src%dat(ibase:ibase+7)*sincT)/cur%d1/cur%d2


 !end subroutine




end module
