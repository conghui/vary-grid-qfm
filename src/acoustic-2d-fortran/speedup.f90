program  speedup
  use sep
  implicit none
  integer :: n1,n2,n3
  real :: d1,d2,d3,o1,o2,o3
  integer :: nt,it
  real :: tmax,tmin,errorFact
  real :: dt,deltaT,vel,dist
  real :: fmax,extent,qvalue,b(3),e(3),ff,cyclesKill,downfact
  real :: dmax,dts,fkill,dsamp,base,aqx,aqy
  real,allocatable :: cells(:,:)
  integer :: nblocks,ierr,error,bound,i
  logical :: follow

  call sep_init()
  call from_param("nblocks",nblocks,80)
  
  call from_param("nt",nt,2000)
  call from_param("dt",dt,.004)
  call from_param("vel",vel,2000.)
  call from_param("error",error,10)
  call from_param("bound",bound,20)
  call from_param("downfact",downfact,.06)
  call from_param("errorfact",errorfact,1.2)
  
  call from_param("fmax",fmax,70.)
  call from_param("aqx",aqx,0.001)
  call from_param("aqy",aqy,0.)
  call from_param("extentxy",extent,13000.)
  
  
  
  
  call from_param("qvalue",qvalue,50.)
  call from_param("follow",follow,.true.)
  deltaT=nt*dt/nblocks
  tmin=0
  ierr=sep_put_data_axis_par("out",1,nblocks,0.,deltaT,"time")
  call to_history("n2",2)
  b=(/0.,-extent,-extent/)
  e=extent
  allocate(cells(nblocks,2))
  dmax=vel/fmax/2.8
  dts=.49*dmax/vel
  base=product((e-b)/dmax)*deltaT/dts
  
  
  
  
  
  do i=1,nblocks
    tmax=tmin+deltaT
    if(follow) then
      dist=vel*tmax
      b(2:3)=max(-extent,-dist);
      e=(/min(extent,dist),min(dist+aqx,extent),min(dist+aqy,extent)/)
    end if
    ff=1.-2.*3.14159265/qvalue
    cyclesKill=log(downfact)/log(ff)
    fkill=cyclesKill/max(tmin,.001)
    dsamp=max(dmax,vel/min(fmax,fkill)/2.8/errorFact)
    write(0,*) "kiling frequencies greater",fkill, " sample ",dsamp
    cells(i,1)=product((e-b)/dsamp)
    dts=.49*dsamp/vel
    cells(i,1)=cells(i,1)*deltaT/dts
    tmin=tmax
end do
!   base=product((e-b)/dmax)*deltaT/dts
    cells(:,2)=base/cells(:,1)
    write(0,*) "TOTAL SPEEDUP",1./(sum(1./cells(:,2))/size(cells,1))
    ierr =srite("out",cells,nblocks*8)
    
    
    
 



end program

