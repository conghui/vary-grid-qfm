module randomMod
implicit none

contains
subroutine addRandom(vel,bnd)
  real :: vel(:,:)
  integer :: bnd,ie1,ie2
  real,allocatable :: my_rand(:)
  real :: vmin,f
  real :: idist,mys
  integer :: nrand,nlen
  integer :: i1,i2,irand,i,id2,id1
  integer :: n1,n2
  real,allocatable :: mymean(:)
  
  nrand=max(size(vel,1),size(vel,2))
  allocate(mymean(2*nrand))
  mymean=0
  vmin=minval(vel)
  nrand=nrand*nrand*6
  allocate(my_rand(nrand))
  call random_number(my_rand)
  irand=1
  do i1=1,size(mymean)
    mymean(i1)=vmin*.1*(real(i1)/real(size(mymean)))
  end do
  n1=size(vel,1)
  n2=size(vel,2)

  


    do i2=1,size(vel,2)
      id2=0
      if(i2 <=bnd) then
        id2=i2-bnd+1
        ie2=bnd+1
      end if
      if(i2> size(vel,2)-bnd) then
        id2=bnd-(size(vel,2)-i2)
        ie2=size(vel,2)-bnd
      end if
      do i1=1,size(vel,1)
        id1=0
        if(i1 <=bnd) then
          id1=i1-bnd+1
          id1=bnd+1
        end if
        if(i1> size(vel,1)-bnd) then
          id1=bnd-(size(vel,1)-i1)
          ie1=size(vel,1)-bnd
        end if
        idist=id1*id1+id2*id2
        if(idist>0) then
          idist=sqrt(idist)
          f=real(idist)/real(size(mymean))
          vel(i1,i2)=.5*vel(ie1,ie2)*my_rand(irand)+&
            (1.-f)*vel(ie1,ie2)+f*mymean(idist)
          irand=irand+1
        end if
    end do
  end do



  do i2=1,size(vel,2)
    !smooth top
    mys=vel(bnd+1,i2)*bnd*2
    nlen=bnd*2+1
    do i1=bnd,1,-1
      mys=mys+vel(i1,i2)
      vel(i1,i2)=mys/nlen
      mys=mys-2*vel(bnd+1,i2)  
      nlen=nlen-1
    end do
    !smooth bottom
    mys=vel(n1-bnd,i2)*bnd*2
    nlen=bnd*2+1
    do i1=size(vel,1)-bnd+1,size(vel,1)
      mys=mys+vel(i1,i2)
      vel(i1,i2)=mys/nlen
      mys=mys-vel(size(vel,1)-bnd,i2)*2
      nlen=nlen-1
    end do
  end do


 do i1=1,size(vel,1)
   !smooth left
    mys=vel(i1,bnd+1)*bnd*2
    nlen=bnd*2+1
    do i2=bnd,1,-1
      mys=mys+vel(i1,i2)
      vel(i1,i2)=mys/nlen
      mys=mys-2*vel(i1,bnd+1)  
      nlen=nlen-1
    end do
    !smooth right
    mys=vel(i1,n2-bnd)*bnd*2
    nlen=bnd*2+1
    do i2=size(vel,2)-bnd+1,size(vel,2)
      mys=mys+vel(i1,i2)
      vel(i1,i2)=mys/nlen
      mys=mys-vel(i1,size(vel,1)-bnd)*2
      nlen=nlen-1
    end do
  end do



  deallocate(my_rand)

end subroutine





end module