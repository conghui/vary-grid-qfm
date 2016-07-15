module resample_mod
  use modeling_mod
  use sinc_mod
  implicit none
  real, allocatable :: sinc_table(:,:)
  integer,private :: ns


contains

  subroutine init_sinc_table(nsinc,npts)
    integer :: nsinc,npts
    real    :: loc
    integer :: n,i

    allocate(sinc_table(nsinc,npts))
    ns=nsinc
    call mksinc_table(npts,nsinc,sinc_table)
  end subroutine


  subroutine interpField(oldS,newS,oldF,newF,extend)
    logical :: extend
    type(modelingT) :: oldS,newS
    real :: oldf(:,:,:), newf(:,:,:)
    integer,allocatable:: t3(:),t2(:),t1(:)
    integer :: i1,i2,i3,b3,e3,b2,e2,b1,e1
    integer ia,ib,ic
    integer,allocatable :: c3(:,:), c2(:,:),c1(:,:)
    real,allocatable :: x1(:),x2(:),x3(:)
    integer,allocatable :: k1(:),k2(:),k3(:)
    allocate(x3(newS%n3), x2(newS%n2),x1(newS%n1))
    b3=newS%n3; b2=newS%n2; b1=newS%n1;
    e3=1; e2=1; e1=1;
    newf=0;

    do i3=1,newS%n3
      x3(i3)=(newS%o3+newS%d3*(i3-1.)-oldS%o3)/oldS%d3+1
      if(x3(i3) > .999 .and. i3 < b3) b3=i3
      if(x3(i3) < oldS%n3+.001 .and. i3 > e3) e3=i3
    end do
    do i2=1,newS%n2
      x2(i2)=(newS%o2+newS%d2*(i2-1.)-oldS%o2)/oldS%d2+1
      if(x2(i2) > .999 .and. i2 < b2) b2=i2
      if(x2(i2) < oldS%n2+.001 .and. i2 > e2) e2=i2
    end do
    do i1=1,news%n1
      x1(i1)=(newS%o1+newS%d1*(i1-1.)-oldS%o1)/oldS%d1+1
      if(x1(i1) > .999 .and. i1< b1) b1=i1
      if(x1(i1) < oldS%n1+.001 .and. i1 > e1) e1=i1
    end do

    allocate(k1(newS%n1),k2(newS%n2),k3(newS%n3))
    !!write(0,*) "OLD 1",oldS%n1,oldS%o1,oldS%d1
    !write(0,*) "NEW 1",newS%n1,newS%o1,newS%d1
    !write(0,*) "MAP 1",minval(x1),maxval(x1)

    !write(0,*) "OLD 2",oldS%n2,oldS%o2,oldS%d2
    !write(0,*) "NEW 2",newS%n2,newS%o2,newS%d2
    !write(0,*) "MAP 2",minval(x2),maxval(x2)
    !write(0,*) "B E ",b1,b2,e1,e2
    !write(0,*) "CMAP",(newS%o2+newS%d2*(news%n2-1)), (newS%o2+newS%d2*(news%n2-1)-oldS%o2)/oldS%d2+1.
    !if(abs(1.-oldS%d1/newS%d1) < .01) then
    if(1==3) then

      !we just are resizing
      k1=x1+.5
      k2=x2+.5
      k3=x3+.5
      where(x1<1)
        x1=1
      end where
      where(x2<1)
        x2=1
      end where
      where (x3 < 1)
        x3 = 1
      end where
      where(x1>size(oldf,1))
        x1=size(oldf,1)
      end where
      where(x2>size(oldf,2))
        x2=size(oldf,2)
      end where
      where (x3 > size(oldF, 3))
        x3 = size(oldF, 3)
      end where
      newf=0
      !$OMP PARALLEL DO private(i1,i2,i3)
      do i3=1, size(newF, 3)
      do i2=1,size(newf,2)
      do i1=1,size(newf,1)
        newf(i1,i2,i3)=oldf(x1(i1),x2(i2),x3(i3))
      end do
      end do
      end do
      !$OMP END PARALLEL DO

    else !We are doing sinc interpolation
      allocate(c3(ns, newS%n3), c2(ns,newS%n2),c1(ns,newS%n1))
      allocate(t3(newS%n3), t2(newS%n2),t1(newS%n1))

      k1=floor(x1)
      k2=floor(x2)
      k3=floor(x3)
      t1=nint((x1-k1)*size(sinc_table,2))+1
      t2=nint((x2-k2)*size(sinc_table,2))+1
      t3=nint((x3-k3)*size(sinc_table,2))+1

      where(t3 > size(sinc_table, 2))
        t3 = size(sinc_table, 2)
      end where
      where(t2 > size(sinc_table,2))
        t2=size(sinc_table,2)
      end where
      where (t1>size(sinc_table,2))
        t1=size(sinc_table,2)
      end where

      do i1=1,ns
        c1(i1,:)=k1-4+i1
        c2(i1,:)=k2-4+i1
        c3(i1,:)=k3-4+i1
      end do
      where(c1 <1)
        c1=1
      end where
      where(c2 <1)
        c2=1
      end where
      where(c3 < 1)
        c3=1
      end where
      where(c1>oldS%n1)
        c1=oldS%n1
      end where
      where(c2>oldS%n2)
        c2=oldS%n2
      end where
      where (c3>oldS%n3)
        c3=oldS%n3
      end where

      if(extend) then
        !$OMP PARALLEL DO private(i1,i2,i3, ic, ib,ia)
        do i3=1,size(newF,3)
          do i2=1,size(newF,2)
            do i1=1,size(newF,1)
              do ic=1,ns
                do ib=1,ns
                  do ia=1,ns
                    newF(i1,i2,i3)=newF(i1,i2,i3)+&
                      oldF(c1(ia,i1),c2(ib,i2),c3(ic,i3))*&
                      sinc_table(ia,t1(i1))*sinc_table(ib,t2(i2))*sinc_table(ic,t3(i3))
                  end do
                end do
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO

      else
        !$OMP PARALLEL DO private(i1,i2,i3, ic, ib,ia)
        do i3=b3,e3
          do i2=b2,e2
            do i1=b1,e1
              do ic=1,ns
                if(c3(ic,i3) >0 .and. c3(ic,i3) <= oldS%n3) then
                  do ib=1,ns
                    if(c2(ib,i2) >0 .and. c2(ib,i2) <=oldS%n2) then
                      do ia=1,ns
                        if(c1(ia,i1) >0 .and. c1(ia,i1) <=oldS%n1) then
                          newF(i1,i2,i3)=newF(i1,i2,i3)+&
                            oldF(c1(ia,i1),c2(ib,i2),c3(ic,i3))*&
                            sinc_table(ia,t1(i1))*sinc_table(ib,t2(i2))*sinc_table(ic,t3(i3))
                        end if
                      end do
                    end if
                  end do
                end if
              end do
            end do
          end do
        end do
        !$OMP END PARALLEL DO
      end if
      deallocate(c1,c2,c3,t1,t2,t3)
    end if
    deallocate(x1,x2,x3,k1,k2,k3)
  end subroutine

end module
