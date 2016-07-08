module data_mod
 use modeling_mod
 use sep_mod
 use sep
 implicit none

 type dataT
   integer :: n1,n2,n3
   real :: o1,o2,o3
   real :: d1,d2,d3
   real :: depth
   real,allocatable :: dat(:,:)
   real,allocatable :: floc(:,:)
   character(len=128) :: tag
 end type
 
 contains
 
 subroutine initData(tag,dat)
   character(len=*) :: tag
   type(dataT) :: dat
   integer :: nx,ny,i2,nt,ns
   real :: ox,oy,dx,os,ds
   integer :: ierr,i
   real, allocatable :: tmp(:,:)
   real :: ot,dt,depth,x,y
   character(len=128) :: tmpS
 
   ierr= sep_get_data_axis_par(tag,1,dat%n1,dat%o1,dat%d1,tmpS)
   ierr= sep_get_data_axis_par(tag,2,dat%n2,dat%o2,dat%d2,tmpS)
   ierr= sep_get_data_axis_par(tag,3,dat%n3,dat%o3,dat%d3,tmpS)
   dat%tag=tag
   
   
   call from_param("depth",dat%depth,0.)
   allocate(dat%floc(2,dat%n2))
   allocate(dat%dat(dat%n1+10,dat%n2))
 
 end subroutine


 real function getShotLocation(dat,ishot) result(res)
   type(dataT) :: dat
   integer :: ishot
   
   res=dat%o3+dat%d3*(ishot-1)
 
 end function

 subroutine readShot(dat,ishot)
   type(dataT) :: dat
   integer :: ishot
   integer :: ierr,i2
   real :: x
   real,allocatable :: tmp(:,:)
   integer :: ng(3),nw(3),fw(3),jw(3)
   jw=1
   fw=0
   ng=(/dat%n1,dat%n2,dat%n3/)
   nw=ng
   nw(3)=1
   fw(3)=ishot-1   
   

   !ierr=sseek_block(dat%tag,ishot,dat%n2*dat%n1*4,0)


     do i2=1,dat%n2
       x=dat%o2+dat%d2*(i2-1)+dat%o3+dat%d3*(ishot)
       dat%floc(:,i2)=(/dat%depth,x/)
     end do

   allocate(tmp(dat%n1,dat%n2))
   ierr=sreed_window(dat%tag,3,ng,nw,fw,jw,4,tmp)
   dat%dat(1:4,:)=0
   dat%dat(5:4+dat%n1,:)=tmp
   dat%dat(dat%n1-3:,:)=0
   deallocate(tmp)
 
 
 end subroutine
 
   
  subroutine cleanData(dat)
    type(dataT) :: dat
    deallocate(dat%dat)
    deallocate(dat%floc)
  end subroutine
  
 
 subroutine addReceiver(dat,tm,rnew,cur)
   type(dataT) :: dat
   integer :: it,ix,iy,iz,ibase,i2,isinc
   real :: rnew(:,:)
   type(modelingT) :: cur
   real,pointer :: sincT(:)
   real :: tm
   
   
   it=tm/dat%d1
   isinc=(tm-it*dat%d1)/dat%d1*10000+1
   sincT=>cur%sincT(:,isinc)
   ibase=it+4
   do i2=1,size(dat%dat,2)
     !for now nearest neighbor
     iz=(dat%floc(1,i2)-cur%o1)/cur%d1+1.5
     ix=(dat%floc(2,i2)-cur%o2)/cur%d2+1.5
      rnew(iz,ix)=rnew(iz,ix)+sum(dat%dat(ibase:ibase+7,i2)*sincT)/cur%d1/cur%d2
   end do
 end subroutine

end module
