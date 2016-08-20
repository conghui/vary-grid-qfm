module migrate_shot_mod
  use sep
  use vel_mod
  use prop_mod
  use my_timers_mod
  use image_mod
  use box_mod
  use data_mod
  use source_mod
  implicit none
  logical, private :: qtest
  integer,private  :: writeSection,icount
  integer, private :: writeCheckPoint

contains

  real function migrate_shot(ishot,verb,vel,dat,source,times) result(factor)
    integer :: ishot
    logical :: verb
    type(velT) :: vel,vuse
    type(dataT) :: dat
    type(sourceT) :: source
    type(timesT)  :: times
    type(imageT) :: iuse
    type(boxT)  :: domain
    type(waveT), target :: p1,p2,p3
    type(waveT), pointer :: pold,pcur,pnew
    type(modelingT),pointer :: cur,old
    type(modelingT) :: full
    type(sourceW) :: sou
    integer :: iblock,it,iimage,slow,ierr
    logical :: xx
    real :: t0
    full%o1=vel%o1; full%d2=vel%d2 ;full%n2=vel%n2;
    full%o2=vel%o2; full%d1=vel%d1; full%n1=vel%n1
    call start_timer_num(timerIS)
    call readShot(dat,ishot)
    call setSourceLocation(source, getShotLocation(dat,ishot))
    call from_param("qtest",qtest,.false.)
    call from_param("writeSection",writeSection,-1)
    call from_param("writeCheckpoint",writeCheckpoint,-1)
    icount=0
    factor=calcShotBox(vel,dat,source,times,domain)
    cur=>domain%hyper(1); old=>cur
    if(writeCheckPoint>0) then
      ierr= sep_put_data_axis_par("wfield",1,vel%n1,vel%o1,vel%d1,"Depth");
      ierr= sep_put_data_axis_par("wfield",2,vel%n2,vel%o2,vel%d2,"X Position");
    end if

    call initWavefields(cur,p1,p2,p3)
    call initSourceSave(sou,domain%totImaging)
    if(writeSection>0) then
      ierr= sep_put_data_axis_par("wfield",1,vel%n1,vel%o1,vel%d1,"Depth");
      ierr= sep_put_data_axis_par("wfield",2,vel%n2,vel%o2,vel%d2,"X Position");
      ierr= sep_put_data_axis_par("wfield",3,domain%hyper(writeSection)%ntblock*4,0.,domain%hyper(writeSection)%dt,"time")
    end if
    allocate(vuse%dat(1,1),vuse%vgamma(1,1))
    call setPropSize(vel,cur,vuse) !ASK: how to change in the 3D case?
    call resampleVel(vel,old,vuse)

    p1%dat=0; p2%dat=0; p3%dat=0
    pold=>p1; pcur=>p2; pnew=>p3
    call stop_timer_num(timerIS)

    !FIRST RUN SOURCE FORWARD
    iimage=1
    do iblock=1,size(domain%hyper)
      write(0,*) "FORWARD BLOCK",iblock
      cur=>domain%hyper(iblock)
      if(.not. checkSame(old,cur)) then
        call start_timer_num(timerR)
        call resampleVel(vel,cur,vuse)
        call resampleP(p1,old,cur)
        call resampleP(p2,old,cur)
        call resampleP(p3,old,cur)

        call stop_timer_num(timerR)
        old=>cur
      end if


      t0=(iblock-1)*domain%blockSize
      xx=.false.
      if(writeSection<=iblock .and. writeSection+3>=iblock .and. writeSection>0) xx=.true.
      if(writesection>0 .and. iblock > writeSection+4) return


      call advanceBlock(t0,iimage,cur,pold,pcur,pnew,vuse,source,sou,xx,full)
      !  call writeFull("wfield.H",cur,full,p3%dat,.false.)
      if(writeCheckPoint==iblock) then
        call writeFull("wfield",cur,full,pcur%dat,.false.)
        return
      end if
    end do



    iimage=iimage-1
    p1%dat=0; p2%dat=0; p3%dat=0
    call createIuse(cur,iuse)
    do iblock=size(domain%hyper),1,-1
      write(0,*) "ADJOINT BLOCK",iblock

      cur=>domain%hyper(iblock)
      if(.not. checkSame(cur,old)) then
        call resampleVel(vel,cur,vuse)
        call resampleP(p1,old,cur)
        call resampleP(p2,old,cur)
        call resampleP(p3,old,cur)
        call resampleI(iuse,old,cur)
        old=>cur
      end if

      t0=iblock*domain%blockSize
      call backwardBlock(t0,iimage,cur,pold,pcur,pnew,vuse,dat,sou,iuse)

      ! call writeFull("wfield.H",cur,full,p3%dat,.false.)

    end do
    call updateImage(cur,iuse)

    ! call print_timers()
    call cleanIuse(iuse)
    call cleanWavefields(p1,p2,p3)
    call cleanSourceSave(sou)
    call cleanBoxShot(domain)

    call cleanVel(vuse)

  end function

  subroutine advanceBlock(t0,iimage,cur,pold,pcur,pnew,vuse,source,sou,writeIt,full)
    type(modelingT) :: cur,full
    type(waveT), pointer :: pold,pcur,pnew,pt
    type(sourceW) :: sou
    type(sourceT) :: source
    type(velT) :: vuse
    integer :: it,iimage,ierr
    logical :: done,writeIt
    real :: dt,t0,tm


    done=.false.
    it=1
    dt=cur%dt
    tm=t0

    !cur%ntblock=cur%ntblock*3
    do while(.not. done)
      if(mod(it,20)==0) write(0,*) "working on",it,cur%ntblock
      call start_timer_num(timerPA)
      if(slow==0 .or. qtest) then
        call advanceWavefieldQ(pold,pcur,pnew,vuse,cur,dt)
      else
        call advanceWavefield(pold,pcur,pnew,vuse,cur,dt)
      end if
      call stop_timer_num(timerPA)

      call start_timer_num(timerI)
      call addSource(source,tm,pnew%dat,cur,1.)
      call stop_timer_num(timerI)
      if(mod(it,cur%jimage)==0 .and. it <=cur%ntblock) then
        call start_timer_num(timerIm)
        allocate(sou%w(iimage)%ar(cur%n1,cur%n2))
        sou%w(iimage)%ar=pnew%dat
        iimage=iimage+1
        call stop_timer_num(timerIm)
      end if

      it=it+1
      if(it==cur%ntblock+1) then
        if(cur%dtExtra > dt) then
          dt=cur%dtExtra
          vuse%vgamma=vuse%vgamma*cur%dt/dt
        else
          done=.true.
        end if
      else if(it >cur%ntblock+1) then

        done=.true.
      end if
      if(writeIt .and. it<=cur%ntblock+1) then
        icount=icount+1
        call writeFull("wfield",cur,full,pnew%dat,.false.)
      end if
      tm=tm+dt
      pt=>pold
      pold=>pcur
      pcur=>pnew
      pnew=>pt

    end do


  end subroutine


  subroutine backWardBlock(t0,iimage,cur,pold,pcur,pnew,vuse,dat,sou,iuse)
    type(modelingT) :: cur
    type(waveT), pointer :: pold,pcur,pnew,pt
    type(sourceW) :: sou
    type(dataT) :: dat
    type(velT) :: vuse
    type(imageT) :: iuse
    integer :: it,iimage
    logical :: done
    real :: dt,t0,tm

    done=.false.



    tm=t0
    if(cur%dtExtra > cur%dt/100.) then
      dt=cur%dtExtra
      it=cur%ntblock+1
    else
      dt=cur%dt
      it=cur%ntblock
    end if


    do while(it >= 1)
    call start_timer_num(timerP)
    call advanceWavefield(pold,pcur,pnew,vuse,cur,dt)
    call stop_timer_num(timerP)

    call start_timer_num(timerI)
    call addReceiver(dat,tm,pnew%dat,cur)
    call stop_timer_num(timerI)

    if(mod(it,cur%jimage)==0 .and. it <= cur%ntblock) then

      call start_timer_num(timerIm)

      call imageIt(sou%w(iimage)%ar,pnew%dat,iuse,cur)

      iimage=iimage-1
      call stop_timer_num(timerIm)
    end if

    tm=tm-dt
    it=it-1
    dt=cur%dt
    pt=>pold
    pold=>pcur
    pcur=>pnew
    pnew=>pt
    end do
  end subroutine


end module
