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
  type(velT) :: vel, bigVel
  type(timesT) :: times
  type(sourceT) :: source
  type(dataT)   :: dat
  type(modelingT) :: small,full
  integer :: ishot,ii
  real :: factor,rnd
  integer, allocatable :: shot(:, :, :) ! first dimension store the real index
  logical, allocatable :: shotL(:, :)


  call initpar()
  call setup_timers()

  call readVel("vel", vel, small)
  call readVel("vel2", bigVel, full)
  bigVel%dat = 0
  call init_sinc_table(8, 10000)

  call interpField(small, full, vel%dat, bigVel%dat, .false.)

  ii = sep_put_data_axis_par('c.H', 1, full%n1, full%o1, full%d1, "Z")
  ii = sep_put_data_axis_par('c.H', 2, full%n2, full%o2, full%d2, "X")
  ii = sep_put_data_axis_par('c.H', 3, full%n3, full%o3, full%d3, "Y")
  ii = srite('c.H', bigVel%dat, size(bigVel%dat)*4)
  write(0,*) 'program exit normally'
end program
