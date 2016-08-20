module sinc_mod

implicit none

interface
  subroutine mksincit(sinc,lsinc,d,space) bind(c,name="mksincit")
    use iso_c_binding
    implicit none
    real(C_FLOAT),dimension(*),intent(out) :: sinc,space
    integer(C_INT),value,intent(in) :: lsinc
    real(C_FLOAT),value,intent(in) :: d
  end subroutine
end interface
contains


subroutine mksinc_table(nt,nsinc,table)
  integer :: nt,nsinc,it;
  real :: table(:,:),dtab,t,space(1000)

  dtab=1.
  dtab=1./nt
  t=0
  do it=1,nt
    call mksincit(table(:,it),nsinc,t,space)
    t=(it-1)*dtab
  end do
end subroutine





end module
