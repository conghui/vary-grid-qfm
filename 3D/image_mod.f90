module image_mod
  !use resampleMod
  use modeling_mod
  use vel_mod
  implicit none
  real, allocatable,private :: imageTot(:,:,:)
  type(modelingT),private :: myD
  type imageT
    integer :: n1,n2,n3,n4
    real    :: o1,o2,o3,o4
    real    :: d1,d2,d3,d4
    real, allocatable, dimension(:,:,:) :: dat
  end type

contains

  subroutine calcImage(vel,tag)
    type(velT) :: vel
    character(len=*) tag
    integer :: ierr
    ierr= sep_put_data_axis_par(tag,1,vel%n1,vel%o1,vel%d1,"Z")
    ierr= sep_put_data_axis_par(tag,2,vel%n2,vel%o2,vel%d2,"X")
    ierr= sep_put_data_axis_par(tag,3,vel%n2,vel%o2,vel%d2,"Y")

    allocate(imageTot(vel%n1,vel%n2,vel%n3))
    imageTot=0
    myD%n1=vel%n1; myd%o1=vel%o1; myd%d1=vel%d1
    myD%n2=vel%n2; myd%o2=vel%o2; myd%d2=vel%d2
    myD%n3=vel%n3; myd%o3=vel%o3; myd%d3=vel%d3


  end subroutine

  subroutine cleanIuse(iuse)
    type(imageT) :: iuse
    deallocate(iuse%dat)
  end subroutine

  subroutine createIuse(cur,iuse)
  type(modelingT) :: cur
  type(imageT) :: iuse

  iuse%n1=cur%n1; iuse%o1=cur%o1; iuse%d1=cur%d1
  iuse%n2=cur%n2; iuse%o2=cur%o2; iuse%d2=cur%d2
  iuse%n3=cur%n3; iuse%o3=cur%o3; iuse%d3=cur%d3

  allocate(iuse%dat(iuse%n1,iuse%n2,iuse%n3))
  iuse%dat=0
  end subroutine

  subroutine resampleI(iuse,old,cur)
  type(modelingT) :: old,cur
  type(imageT) :: iuse
  real, allocatable :: tmp(:,:,:)
  integer :: ierr
  allocate(tmp(myd%n1,myd%n2,myd%n3))
  call InterpField(old,myD,iuse%dat,tmp,.false.)
  imageTot=imageTot+tmp/cur%ntblock
  iuse%n1=cur%n1; iuse%o1=cur%o1; iuse%d1=cur%d1
  iuse%n2=cur%n2; iuse%o2=cur%o2; iuse%d2=cur%d2
  iuse%n3=cur%n3; iuse%o3=cur%o3; iuse%d3=cur%d3
  deallocate(tmp)
  deallocate(iuse%dat)
  allocate(iuse%dat(cur%n1,cur%n2,cur%n3))
  iuse%dat=0
  end subroutine

  subroutine updateImage(cur,iuse)
    type(modelingT) :: cur
    real, allocatable :: tmp(:,:,:)
    type(imageT) :: iuse
    integer :: ierr
    allocate(tmp(myd%n1,myd%n2,myd%n3))
    call InterpField(cur,myD,iuse%dat,tmp,.false.)
    imageTot=imageTot+tmp/cur%ntblock
    deallocate(tmp)
  end subroutine

  subroutine imageIt(pnew,rnew,iuse,cur)
    integer ::   i1,i2,i3
    type(imageT) :: iuse
    real :: pnew(:,:,:),rnew(:,:,:)
    type(modelingT) :: cur

    !$OMP PARAlLEL DO private(i1,i2,i3)
    do i3=1,size(pnew,3)
      do i2=1,size(pnew,2)
        do i1=1,size(pnew,1)
          iuse%dat(i1,i2,i3)=iuse%dat(i1,i2,i3)+pnew(i1,i2,i3)*&
            rnew(i1,i2,i3)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine

  subroutine writeImage(tag,cur,iuse)
    type(modelingT) :: cur
    real, allocatable :: tmp(:,:,:)
    type(imageT) :: iuse
    integer :: ierr
    character(len=*) :: tag

    allocate(tmp(myd%n1,myd%n2,myd%n3))
    call InterpField(cur,myD,iuse%dat,tmp,.false.)
    imageTot=imageTot+tmp
    ierr=srite(tag,imageTot,size(imageTot)*4)
    deallocate(tmp);
  end subroutine

  !subroutine writeFullImage(tag)
  !integer :: ierr
  !character(len=*) :: tag
  !integer :: ng(2),nw(2),fw(2),jw(2)
  !ng=(/size(imageTot,1),size(imageTot,2)/)
  !nw=ng
  !fw=0
  !jw=1
  !ierr=srite_window(tag,2,ng,nw,fw,jw,4,imageTot);
  !end subroutine
end module

