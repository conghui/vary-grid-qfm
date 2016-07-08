module prop_mod
  use modeling_mod
  use resampleMod
  use vel_mod
  use sep
  implicit none
  real, private :: factor,gamma
  real, private :: b_b(20)
  type waveT
    real,allocatable :: dat(:,:)
  end type
   
  interface

 end interface
  
  contains
  
subroutine initWavefields(cur,p1,p2,p3)
  type(modelingT) :: cur
  type(waveT) :: p1,p2,p3
    allocate(p1%dat(cur%n1,cur%n2))
    allocate(p2%dat(cur%n1,cur%n2))
    allocate(p3%dat(cur%n1,cur%n2))
end subroutine

subroutine cleanWavefields(p1,p2,p3)
  type(waveT) :: p1,p2,p3
  deallocate(p1%dat,p2%dat,p3%dat)
end subroutine
  
subroutine propWave(prev,cur,new,vsq,d1a,d2a,d0,b,e,b1,e1)
  real  :: prev(:,:),cur(:,:),new(:,:),vsq(:,:)
  double precision  :: d1a(:),d2a(:)
  real   ::d0,gamma,factor
  integer :: b(:),e(:),b1,e1
  real ::d11,d12,d13,d14,d15,d21,d22,d23,d24,d25
  integer :: i1,i2
  real :: tau,taucur,tauprev
  
  d11=d1a(1); d12=d1a(2); d13=d1a(3); d14=d1a(4); d15=d1a(5);
  d21=d2a(1); d22=d2a(2); d23=d2a(3); d24=d2a(4); d25=d2a(5);
  
!  write(0,*) "CHECK D11",d11,d21
 ! call seperr("")
  do i2=b(2),e(2)
   do i1=b(1),e(1)
     new(i1,i2)=vsq(i1,i2)*(&
       d11*((cur(i1-1,i2)+cur(i1+1,i2)))+&
       d12*((cur(i1-2,i2)+cur(i1+2,i2)))+&
       d13*((cur(i1-3,i2)+cur(i1+3,i2)))+&
       d14*((cur(i1-4,i2)+cur(i1+4,i2)))+&
       d15*((cur(i1-5,i2)+cur(i1+5,i2)))+&
       d21*((cur(i1,i2-1)+cur(i1,i2+1)))+&
       d22*((cur(i1,i2-2)+cur(i1,i2+2)))+&
       d23*((cur(i1,i2-3)+cur(i1,i2+3)))+&
       d24*((cur(i1,i2-4)+cur(i1,i2+4)))+&
       d25*((cur(i1,i2-5)+cur(i1,i2+5)))+&
       d0*( (cur(i1,i2))))+&
       2*cur(i1,i2)-prev(i1,i2)
    end do
  end do
end subroutine
subroutine propWaveQ(prev,cur,new,vsq,vgamma,d1a,d2a,d0,b,e)
  real  :: prev(:,:),cur(:,:),new(:,:),vsq(:,:),vgamma(:,:)
  double precision  :: d1a(:),d2a(:)
  real   ::d0
  integer :: b(:),e(:)
  real ::d11,d12,d13,d14,d15,d21,d22,d23,d24,d25
  integer :: i1,i2
  real :: tau,taucur,tauprev
  
  d11=d1a(1); d12=d1a(2); d13=d1a(3); d14=d1a(4); d15=d1a(5);
  d21=d2a(1); d22=d2a(2); d23=d2a(3); d24=d2a(4); d25=d2a(5);

  do i2=b(2),e(2)
   do i1=b(1),e(1)
     tau=vgamma(i1,i2)
     taucur=1-tau
     tauprev=tau

     new(i1,i2)=vsq(i1,i2)*(&
       d11*(taucur*(cur(i1-1,i2)+cur(i1+1,i2))+tauprev*(prev(i1-1,i2)+prev(i1+1,i2)))+&
       d12*(taucur*(cur(i1-2,i2)+cur(i1+2,i2))+tauprev*(prev(i1-2,i2)+prev(i1+2,i2)))+&
       d13*(taucur*(cur(i1-3,i2)+cur(i1+3,i2))+tauprev*(prev(i1-3,i2)+prev(i1+3,i2)))+&
       d14*(taucur*(cur(i1-4,i2)+cur(i1+4,i2))+tauprev*(prev(i1-4,i2)+prev(i1+4,i2)))+&
       d15*(taucur*(cur(i1-5,i2)+cur(i1+5,i2))+tauprev*(prev(i1-5,i2)+prev(i1+5,i2)))+&
       d21*(taucur*(cur(i1,i2-1)+cur(i1,i2+1))+tauprev*(prev(i1,i2-1)+prev(i1,i2+1)))+&
       d22*(taucur*(cur(i1,i2-2)+cur(i1,i2+2))+tauprev*(prev(i1,i2-2)+prev(i1,i2+2)))+&
       d23*(taucur*(cur(i1,i2-3)+cur(i1,i2+3))+tauprev*(prev(i1,i2-3)+prev(i1,i2+3)))+&
       d24*(taucur*(cur(i1,i2-4)+cur(i1,i2+4))+tauprev*(prev(i1,i2-4)+prev(i1,i2+4)))+&
       d25*(taucur*(cur(i1,i2-5)+cur(i1,i2+5))+tauprev*(prev(i1,i2-5)+prev(i1,i2+5)))+&
       d0*( taucur*(cur(i1,i2))+tauprev*(prev(i1,i2))))+&
       2*cur(i1,i2)-prev(i1,i2)
    end do
  end do
end subroutine

subroutine resampleP(fld,old,cur)
  type(waveT) :: fld
  type(modelingT) :: cur,old
  real, allocatable :: tmp(:,:)
  
  allocate(tmp(old%n1,old%n2))
  tmp=fld%dat
  
  deallocate(fld%dat)
  allocate(fld%dat(cur%n1,cur%n2))
  call interpField(old,cur,tmp,fld%dat,.false.)
  deallocate(tmp)

end subroutine
subroutine setPropSize(vel,cur,vuse)
  type(modelingT) :: cur
  type(velT) :: vuse,vel
  real :: qfact,pi
  integer :: i1



  do i1=1,20
     b_b(i1)=exp(-(.3/20.*real(21-i1))**2)

  end do
end subroutine
subroutine advanceWavefieldQ(old,cur,new,vsq,s,dt)
  type(waveT) :: old,cur,new
  type(modelingT) :: s
  type(velT) :: vsq
  integer :: isect
  real :: dt
  !$OMP PARALLEL DO private(isect)
  do isect=1,s%nsect
    call propwaveQ(old%dat,cur%dat,new%dat,vsq%dat,vsq%vgamma,s%d1a*dt*dt,s%d2a*dt*dt,&
     s%d0*dt*dt,s%b(:,isect),s%e(:,isect))  
  end do
  call boundCond(cur%dat,new%dat)
end subroutine
subroutine advanceWavefield(old,cur,new,vsq,s,dt)
  type(waveT) :: old,cur,new
  type(modelingT) :: s
  type(velT) :: vsq
  integer :: isect
  real :: dt
 ! write(0,*) "OOOH",s%d1a,dt
  !$OMP PARALLEL DO private(isect)
  do isect=1,s%nsect
    call propwave(old%dat,cur%dat,new%dat,vsq%dat,s%d1a*dt*dt,s%d2a*dt*dt,&
     s%d0*dt*dt,s%b(:,isect),s%e(:,isect),s%b1,s%b2)  
  end do
  !$OMP END PARALLEL DO
  call boundCond(cur%dat,new%dat)
end subroutine
subroutine boundCond(cur,new)
  real :: cur(:,:),new(:,:)
  
  integer :: i1,i2
  !$OMP PARALLEL DO private(i1,i2)
  do i2=1,size(cur,2)
    do i1=1,20
      cur(i1,i2)=cur(i1,i2)*b_b(i1)
      cur(size(cur,1)-i1+1,i2)=cur(size(cur,1)-i1+1,i2)*b_b(i1)
      new(i1,i2)=new(i1,i2)*b_b(i1)
      new(size(cur,1)-i1+1,i2)=new(size(cur,1)-i1+1,i2)*b_b(i1)
    end do
  end do
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO private(i1,i2)
  do i2=1,20
    do i1=1,size(cur,1)
      cur(i1,i2)=cur(i1,i2)*b_b(i2)
      cur(i1,size(cur,2)-i2+1)=cur(i1,size(cur,2)-i2+1)*b_b(i2)
      new(i1,i2)=new(i1,i2)*b_b(i2)
      new(i1,size(cur,2)-i2+1)=new(i1,size(cur,2)-i2+1)*b_b(i2)
    end do
  end do
  !$OMP END PARALLEL DO
end subroutine

end module