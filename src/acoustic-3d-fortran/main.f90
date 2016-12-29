program RTM_FAST
  use my_timers_mod
  use vel_mod
  use data_mod
  use image_mod
  use source_mod
  use box_mod
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
  integer, allocatable :: shot(:, :, :) ! first dimension store the real index
  logical, allocatable :: shotL(:, :)

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

  call init_sinc_table(8,10000)
  call stop_timer_num(timerInit)


  allocate(shot(2, dat%n4, dat%n5), shotL(dat%n4, dat%n5))
  shotL=.false.

  write(0,*) 'dat%n4',dat%n4, ' dat%n5', dat%n5

  ! ASK: what's doing in 2D cases?
  iyshot=1
  do while(iyshot<= dat%n5)
    ixshot = 1
    do while(ixshot<= dat%n4)
      call random_number(rnd)
      ii=max(1,min(dat%n4,nint(rnd*dat%n4)+1))
      jj=max(1,min(dat%n5,nint(rnd*dat%n5)+1))
      if(.not. shotL(ii, jj)) then
        shot(:,ixshot, iyshot)=(/ii, jj/)
        ixshot = ixshot + 1
      end if
    end do
    iyshot = iyshot + 1
  end do

  do iyshot = 1, dat%n5
    do ixshot = 1, dat%n4
      factor= migrate_shot(shot(:, ixshot, iyshot),verb,vel,dat,source,times)
      write(0,*) ixshot, " MIGRATING SHOT", shot(1, ixshot, iyshot), " of ",dat%n4, " sped up ",factor
      !if(mod(ishot,5)==0)  call writeFullImage("img")
    end do
  end do
  call print_timers()
  call writeFullImage("img")

end program
