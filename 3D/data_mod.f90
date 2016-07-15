module data_mod
  use modeling_mod
  use sep_mod
  use sep
  implicit none

  type dataT
    integer :: n1,n2,n3,n4,n5
    real :: o1,o2,o3,o4,o5
    real :: d1,d2,d3,d4,d5
    real :: depth
    real,allocatable :: dat(:,:,:)
    real,allocatable :: floc(:,:,:)
    character(len=128) :: tag
  end type

contains

  subroutine initData(tag,dat)
    character(len=*) :: tag
    type(dataT) :: dat
    integer :: ierr
    character(len=120) :: tmpS

    ierr= sep_get_data_axis_par(tag,1,dat%n1,dat%o1,dat%d1,tmpS)
    ierr= sep_get_data_axis_par(tag,2,dat%n2,dat%o2,dat%d2,tmpS)
    ierr= sep_get_data_axis_par(tag,3,dat%n3,dat%o3,dat%d3,tmpS)
    ierr= sep_get_data_axis_par(tag,4,dat%n4,dat%o4,dat%d4,tmpS)
    ierr= sep_get_data_axis_par(tag,5,dat%n5,dat%o5,dat%d5,tmpS)
    dat%tag=tag

    call from_param("depth",dat%depth,0.)
    allocate(dat%floc(3,dat%n2,dat%n3)) !ASK: right?
    allocate(dat%dat(dat%n1+10,dat%n2,dat%n3)) ! ASK: why n1+10, not both n1,n2

  end subroutine


  function getShotLocation(dat,ishot) result(res)
    type(dataT) :: dat
    integer :: ishot(:)
    integer :: sx, sy
    real :: res(2)

    sx = dat%o4 + dat%d4*(ishot(1)-1)
    sy = dat%o5 + dat%d5*(ishot(2)-1)

    res = (/sx, sy/)
    !res=dat%o3+dat%d3*(ishot-1)

  end function

  subroutine readShot(dat,ishot)
    type(dataT) :: dat
    integer :: ishot(:)
    integer :: ixshot, iyshot
    integer :: ierr,i2, i3
    real :: x, y
    real,allocatable :: tmp(:,:,:)
    integer :: ng(5),nw(5),fw(5),jw(5)

    ixshot = ishot(1)
    iyshot = ishot(2)
    jw=1
    fw=0
    ng=(/dat%n1,dat%n2,dat%n3,dat%n4,dat%n5/)
    nw=ng
    nw(4)=1
    nw(5)=1
    fw(4)=ixshot-1
    fw(5)=iyshot-1


    !ierr=sseek_block(dat%tag,ishot,dat%n2*dat%n1*4,0)

    do i3 = 1, dat%n3
      do i2 = 1, dat%n2
        x=dat%o2+dat%d2*(i2-1)+dat%o4+dat%d4*(ixshot)
        y=dat%o3+dat%d3*(i3-1)+dat%o5+dat%d5*(iyshot)
        dat%floc(:,i2,i3)=(/dat%depth,x,y/)
      end do
    end do

    allocate(tmp(dat%n1,dat%n2,dat%n3))
    ierr=sreed_window(dat%tag,5,ng,nw,fw,jw,4,tmp)
    dat%dat(1:4,:,:)=0
    dat%dat(5:4+dat%n1,:,:)=tmp
    dat%dat(dat%n1-3:,:,:)=0
    deallocate(tmp)


  end subroutine


  !subroutine cleanData(dat)
    !type(dataT) :: dat
    !deallocate(dat%dat)
    !deallocate(dat%floc)
  !end subroutine


  subroutine addReceiver(dat,tm,rnew,cur)
    type(dataT) :: dat
    integer :: it,ix,iy,iz,ibase,i2,i3,isinc
    real :: rnew(:,:,:)
    type(modelingT) :: cur
    real,pointer :: sincT(:)
    real :: tm

    it=tm/dat%d1
    isinc=(tm-it*dat%d1)/dat%d1*10000+1
    sincT=>cur%sincT(:,isinc)
    ibase=it+4

    do i3=1,size(dat%dat,3)
      do i2=1,size(dat%dat,2)
        !for now nearest neighbor
        iz=(dat%floc(1,i2,i3)-cur%o1)/cur%d1+1.5
        ix=(dat%floc(2,i2,i3)-cur%o2)/cur%d2+1.5
        iy=(dat%floc(3,i2,i3)-cur%o3)/cur%d3+1.5
        rnew(iz,ix,iy)=rnew(iz,ix,iy)+sum(dat%dat(ibase:ibase+7,i2,i3)*sincT)/cur%d1/cur%d2/cur%d3
      end do
    end do
  end subroutine

end module
