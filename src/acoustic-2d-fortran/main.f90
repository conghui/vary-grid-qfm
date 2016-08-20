program RTM_FAST
  use migrate_shot_mod
  use sep
  logical :: verb
  type(velT) :: vel
  type(timesT) :: times
  type(sourceT) :: source
  type(dataT)   :: dat
  type(modelingT) :: full
  integer :: ishot,ii
  real :: factor,rnd
  integer, allocatable :: shot(:)
  logical, allocatable :: shotL(:)
  
  
  call initpar()
  call setup_timers()
  
  !SETUP global migration parameters
  call start_timer_num(timerInit)
  call from_param("verbose",verb,.false.)
  call readVel("vel",vel,full)
  call initData("data",dat)
  call calcImage(vel,"img")
  call readSource("source",source)
  call initTimes(times)
  call initBox(dat,source,vel,verb)
  
  call  init_sinc_table(8,10000)
  call stop_timer_num(timerInit)
  
  
  allocate(shot(dat%n3),shotL(dat%n3))
  shotL=.false.
  ishot=1
  do while(ishot<=dat%n3)
    call random_number(rnd)
    ii=max(1,min(dat%n3,nint(rnd*dat%n3)+1))
    if(.not. shotL(ii)) then
      shot(ishot)=ii
      ishot=ishot+1
    end if
  end do
  
  
  do ishot=1,dat%n3
      factor= migrate_shot(shot(ishot),verb,vel,dat,source,times)
      write(0,*) ishot," MIGRATING SHOT", shot(ishot), " of ",dat%n3, " sped up ",factor
      if(mod(ishot,5)==0)  call writeFullImage("img")
  end do
  call print_timers()
  call writeFullImage("img")
  
end program
  
  
  
  