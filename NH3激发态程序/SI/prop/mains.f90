!==========================================================================
program main
  use nndms
  implicit none
  real*8, parameter :: au2cm=219474d0
  real*8, parameter :: ang2rad=0.0174533
  real*8, dimension(:,:), allocatable :: gref,dsoc
  real*8, dimension(:,:,:), allocatable :: geom
  real*8 :: r,phi,theta
  real*8 :: eshift
  integer :: i,j,k,total,ios,idx

  !initialize NN PES
  call dmsinit

  total=500

  allocate(geom(3,natoms,total))
  allocate(dsoc(3,nstates))
  allocate(gref(3,natoms))

  gref(:,1)=(/0.00000000,0.00000000,0.00000000/)
  gref(:,2)=(/1.91568930,-0.00000000,0.00000000/)
  gref(:,3)=(/-0.54015829,1.84123135,0.00000000/)
  gref(:,4)=(/-1.15020708,-1.53611589,0.00000000/)

  open(200,file='output')
  do i=1,10
    do j=0,30,5
      geom(:,:,i)=gref
      r=norm2(gref(:,4))+0.4*i
      phi=j*ang2rad
      theta=acos(-gref(1,4)/r)
      geom(1,4,i)=-r*dcos(phi)*dcos(theta)
      geom(2,4,i)=-r*dcos(phi)*dsin(theta)
      geom(3,4,i)=r*dsin(phi)
      call Evaluatediasoc(geom(:,:,i),dsoc)
      write(200,"(2f8.3,6f14.8)") r,dble(j),dsoc(:,1),dsoc(:,2)
    end do
  end do

  stop
end
!==========================================================================
