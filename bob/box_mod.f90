module box_mod
  use modeling_mod
  use sep
  use vel_mod
  use data_mod
  use source_mod
  implicit none
  logical, private :: v
  real, private :: vmin,vmax,dmin,dmax
  double precision :: tot1Cells,tot2Cells
  integer :: totImaging
  integer,private :: timeBlocks
  real,private :: dstable
  integer :: slow
  integer, private :: bound,error,jimage
  real,private :: errorFact,downFact,qfact,maxF,dtBig
  real, pointer :: sincT(:,:)

  type timesT
    integer :: n1,n2,n3
    real    :: o1,o2,o3
    real    :: d1,d2,d3
    real,allocatable :: vals(:,:,:)
  end type

  type boxT
    type(modelingT), pointer :: hyper(:)
    integer :: nt
    real    :: dt
    integer :: totImaging
    real :: blockSize

  end type

  interface
    subroutine fastmarch(order,s1,s2,s3,b1,b2,b3,n1,n2,n3,d1,d2,d3,s0,slow,ttime) bind(c,name="fastmarch")
      import
      real(c_float),intent(in),value:: s1,s2,s3,d1,d2,d3,s0
      integer(c_int), value,intent(in) :: order,b1,b2,b3,n1,n2,n3
      real(c_float),dimension(*) :: slow,ttime
    end subroutine
  end interface

contains

  subroutine initTimes(timesS)
    type(timesT) :: timesS
    integer :: ierr
    character(len=128) :: label 
    ierr=sep_get_data_axis_par("times",1,timesS%n1,timesS%o1,timesS%d1,label)
    ierr=sep_get_data_axis_par("times",2,timesS%n2,timesS%o2,timesS%d2,label)
    ierr=sep_get_data_axis_par("times",3,timesS%n3,timesS%o3,timesS%d3,label)
    allocate(timess%vals(timess%n1,timess%n2,timess%n3))
    ierr=sreed("times",timesS%vals,size(timesS%vals)*4)

  end subroutine

  subroutine cleanBoxShot(boxS)
    type(boxT) :: boxS
    integer :: i
    do i=1,size(boxS%hyper)
    call clean_modeling(boxS%hyper(i))
    end do
    deallocate(boxS%hyper)
  end subroutine

  subroutine initBox(dat,source,vel,verbose)
    type(dataT) :: dat
    type(sourceT) :: source
    type(velT) :: vel
    logical :: verbose
    v=verbose
    call from_param("timeBlocks",timeBlocks,40);
    vmin=minval(vel%dat)
    vmax=maxval(vel%dat)
    dmin=min(vel%d1,vel%d2)
    dmax=max(vel%d1,vel%d2)
    dstable=.499*dmin/vmax !maximum sampling in time for stability

    allocate(sincT(8,10000))
    call mksinc_table(10000,8,sincT)
    call from_param("maxF",maxF,80.)
    call from_param("bound",bound,20)
    call from_param("error",error,20)
    call from_param("errorFact",errorFact,1.2)
    call from_param("qfact",qfact,50.)
    call from_param("downfact",downfact,.04)
    call from_param("jimage",jimage,6)
    call from_param("slow",slow,0)
    dtBig=calcGoodSampling(source,dat,dstable)
  end subroutine

  real function calcShotBox(vel,dat,source,times,domain)
    type(velT) :: vel
    type(dataT) :: dat
    type(boxT)  :: domain
    type(timesT) :: times
    type(sourceT) :: source
    real :: dtuse,blockTime
    integer :: i3,i2,i1,ishot,i
    real, allocatable :: minT(:,:)
    real :: timeMin,timeMax

    ! ASK: what i3 means?
    ! ASK: what times means?
    i3=(source%x-times%o3)/times%d3+1.5
    allocate(minT(times%n1,times%n2))

    minT=times%vals(:,:,i3)

    do ishot=1,size(dat%floc,2)
      i3=min(max(nint((dat%floc(2,ishot)-times%o3)/times%d3+1.5),1),times%n3)
      !$OMP PARALLEL DO private (i1,i2)
      do i2=1,times%n2
        do i1=1,times%n1
          minT(i1,i2)=min(minT(i1,i2),times%vals(i1,i2,i3))
        end do
      end do
    end do

    allocate(domain%hyper(timeBlocks))
    blockTime=dat%n1*dat%d1/timeBlocks
    domain%blockSize=blockTIme
    tot1cells=0
    tot2cells=0
    totImaging=0

    do i=1,size(domain%hyper)
      domain%hyper(i)%bnd=bound
      domain%hyper(i)%dt=domain%dt
      timemin=(i-1)*blockTime
      timemax=timemin+blockTime
      domain%hyper(i)%sincT=>sincT
      call createModeling(domain%hyper(i),source,dat,times,timemin,timeMax,&
        minT,vel,bound,error,errorFact,qfact,downfact)
    end do

    calcShotBox=tot2cells/tot1cells

    if(v) write(0,*) size(domain%hyper),"total speedup factor",tot2cells/tot1cells
    domain%totImaging=totImaging
  end function


  real function calcGoodSampling(source,dat,dmax) result(dt)
    integer :: jdiff
    logical :: found
    integer :: oversampS
    type(sourceT) :: source
    type(dataT) :: dat
    real :: dmax

    dt=dmax*.95
  end function




  subroutine createModeling(mod,source,dat,times,timeMin,timeMax,minT,vel,bound,error,errorFact,qfact,downFact)
    type(modelingT) :: mod
    type(timesT) :: times
    type(sourceT) :: source
    type(dataT) :: dat
    real :: timeMin,timeMax
    real :: minT(:,:)
    type(velT) :: vel
    integer :: bound,error
    real :: errorFact,qFact,downFact
    real :: mymin(2),mymax(2)
    real :: cyclesKill,fkill
    integer :: i1,i2,i3,n(2)
    real :: dsamp,o(2),d(2)
    real :: basic(5)
    integer :: bsize
    integer :: bb(100),ee(100),cb(100),ce(100)
    integer :: nb(2)
    integer :: i
    real ::dtmax

    real :: os(2),ds(2),bv(2),ev(2)
    integer :: ns(2)
    double precision :: fullsize,smallsize
    real ::ff

    ns=(/size(mint,1),size(mint,2)/)
    os=(/times%o1,times%o2/)
    ds=(/times%d1,times%d2/)
    call findExts(os,ds,ns,minT,timeMax,mymin,mymax)

    if(timeMin>.0001) then
      ff=1.-2.*3.14159265/qfact
      cyclesKill=log(downfact)/log(ff)
      fkill=cyclesKill/max(timeMin,.001)
      !I don't care about things downfact less than its original ampliute

      dsamp=max(dmax,vmin/min(maxF,fkill)/3.3/errorFact)
      ! if(v) write(0,*) "killing frequencies greater than ",fkill," at ",timeMin,cyclesKill,dsamp
    else
      dsamp=dmax
    end if

    if (slow==1) dsamp=dmax
    fullsize=vel%n1+2*mod%bnd
    fullsize=fullsize*(vel%n2+mod%bnd*2)

    bv=(/vel%o1,vel%o2/)
    ev=bv+((/vel%n1,vel%n2/)-1)*(/vel%d1,vel%d2/)
    do i=1,2
      o(i)=max(bv(i),mymin(i)-dsamp*error)
      mymax(i)=min(mymax(i)+dsamp*error,ev(i))
      n(i)=ceiling((mymax(i)-mymin(i))/dsamp)+1
      !    mymin(i)=bv(i)
      !    o(i)=bv(i)
      ! mymax(i)=ev(i)
      ! n(i)=ceiling((mymax(i)-mymin(i))/dsamp)+1

      if(slow==1)then
        mymin(i)=bv(i)
        o(i)=bv(i)
        mymax(i)=ev(i)
        n(i)=ceiling((mymax(i)-mymin(i))/dsamp)+1
      end if
    end do


    mod%o1=o(1)-dsamp*(mod%bnd)
    mod%o2=o(2)-dsamp*(mod%bnd)

    mod%d1=dsamp; mod%d2=dsamp;
    mod%n1=n(1)+2*(mod%bnd)
    mod%n2=n(2)+2*(mod%bnd)
    dtmax=.49*dsamp/vmax
    mod%dt=calcGoodSampling(source,dat,dtmax) 
    ! mod%dt=.00025
    if(slow==1) mod%dt=.3*dsamp/vmax
    mod%ntblock=(timeMax-timeMin)/mod%dt
    mod%dtextra=(timeMax-timeMin)-mod%dt*mod%ntblock

    if(slow==1) then
      !  mod%dt=mod%dt/2. !Not sure why blowing up right now
    end if

    basic=(/3.333333333,-.4761904762,.0793650794,-.0099206349,.0006349206/)
    basic=(/8./5.,-.2,8./315.,-1./560.,0./)
    !basic=(/1.,0.,0.,0.,0./)
    !basic=2.*basic/sum(basic)

    mod%d1a=basic!*mod%dt*mod%dt!
    mod%d1a=mod%d1a/dsamp/dsamp
    mod%d2a=basic/dsamp/dsamp

    mod%d0=-2.*(sum(mod%d1a)+sum(mod%d2a))

    ! ASK: what breakAxis means?
    nb(1)=breakAxis(mod%n2,50,bb,ee,5)


    mod%nsect=nb(1)
    allocate(mod%b(2,mod%nsect),mod%e(2,mod%nsect))
    i=0
    do i1=1,nb(1)
      i=i+1
      mod%b(:,i)=(/6,bb(i1)/)
      mod%e(:,i)=(/mod%n1-5,ee(i1)/)
    end do
    smallsize=mod%n1*mod%n2
    mod%jimage=jimage


    tot1cells=tot1cells+smallsize*mod%ntblock
    tot2cells=tot2cells+fullsize*(timeMax-timeMin)/dtBig
    totImaging=totImaging+mod%ntblock/jImage
    !  if(v) write(0,*) "new grid",mod%n1,mod%n2
    if(v) write(0,*) "speedup factor",(fullsize*(timeMax-timeMin)/dtBig)/(smallsize*mod%ntblock),mod%dt
  end subroutine

  integer function breakAxis(n,nblock,b,e,edge) 
    integer :: n,b(100),e(100),edge,nblock,ns
    integer bb,nfinish,nleft,ndo,nsect,i
    ns=n-2*edge
    bb=edge+1
    nsect=ceiling(real(ns)/real(nblock))
    nleft=nsect
    nfinish=ns
    do i=1,nsect
      b(i)=bb
      ndo=ceiling(real(nfinish)/real(nleft))
      e(i)=b(i)-1+ndo
      bb=e(i)+1
      nfinish=nfinish-ndo
      nleft=nleft-1
    end do
    breakAxis=nsect
  end function



  ! conghui: find the minimum/maximum location
  subroutine findExts(o,d,n,minT,maxT,mymin,mymax)
    integer :: i1,i2,i3,n(2)
    real :: o(2),d(2),minT(:,:)
    real :: maxT
    real :: mymin(2),mymax(2)
    mymin=o+d*n
    mymax=o
    do i2=1,size(minT,2)
      do i1=1,size(minT,1)
        if(minT(i1,i2) <= maxT) then
          mymin(1)=min(mymin(1),o(1)+d(1)*(i1-1))
          mymax(1)=max(mymax(1),o(1)+d(1)*(i1-1))
          mymin(2)=min(mymin(2),o(2)+d(2)*(i2-1))
          mymax(2)=max(mymax(2),o(2)+d(2)*(i2-1)) 

        end if
      end do
    end do

  end subroutine





end module
