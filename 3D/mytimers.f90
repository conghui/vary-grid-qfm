module my_timers_mod
use sep_timers_mod
integer :: timerInit,timerIS,timerCS,timerR,timerPA,timerP,timerI,timerIM

contains

subroutine setup_timers()
logical :: err


  err=init_sep_timers()
  err=setup_next_timer("Setup Migration",timerInit)
  err=setup_next_timer("Setup shot",timerIS)
  err=setup_next_timer("Resample medium ",timerR)
  err=setup_next_timer("Propagate attenuation",timerPA)
  err=setup_next_timer("Propagate no attenuation",timerP)
  err=setup_next_timer("Inject Source",timerI)
  err=setup_next_timer("Imaging condition",timerIm)


end subroutine







end module
