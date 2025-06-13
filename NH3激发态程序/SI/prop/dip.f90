!===================================================================================
!Module for NN DMS
module nndms
  use nn_class
  implicit none
  integer :: natoms
  type(nn), allocatable :: ANN(:) !Feed-Forward Neural Network
  integer :: na ! number of NNs
  integer :: ncoord   !total number of internal coords
  integer :: nstates  !number of electronic states
  character(len=99), dimension(:),  allocatable ::  wfile
end module nndms
!===================================================================================
subroutine dplcoord(natom,ncoords,cgeom,igeom)
  implicit none
  integer, intent(in) :: natom, ncoords
  real*8, intent(in) :: cgeom(3,natom)
  real*8, intent(out) :: igeom(ncoords)
  real*8 :: dx(3),dcoord(12),x(3,4)
  integer :: i,j,k,info

  x=cgeom(1:3,1:4)
  call order(x,info)

  k=0
  do i=1,natom-1
    do j=i+1,natom
      k=k+1
      dx=cgeom(:,i)-cgeom(:,j)
      igeom(k)=1.d0/sqrt(dot_product(dx,dx))
    end do
  end do

  call OOP(cgeom(:,1),cgeom(:,2),cgeom(:,3),cgeom(:,4),igeom(ncoords),dcoord,0)

  return
end
!===============================================================================
!initialize NN DMS
subroutine dmsinit
  use nndms
  implicit none
  integer :: i

  natoms=4
  nstates=2
  ncoord=7
  na=1

  allocate(ANN(na),wfile(na))

  wfile(1)='WBdip.txt'

  do i=1,na
    ANN(i)%structfl='struct.dip'
    call ANN(i)%init()
    call ANN(i)%savenn(wfile(i),0)
  end do

  return
end subroutine dmsinit
!==================================================================================
subroutine EvaluateDpl0(i,igeom,dpl)
  use nndms
  implicit none
  integer, intent(in) :: i
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  integer :: R

  R=ANN(i)%RX
  call ANN(i)%output(igeom(1:R))

  !(1,1)
  dpl(1,1,1)=ANN(i)%y(1)
  dpl(2,1,1)=ANN(i)%y(2)
  dpl(3,1,1)=ANN(i)%y(3)*igeom(ncoord)
  !(1,2)
  dpl(1,1,2)=ANN(i)%y(4)*igeom(ncoord)
  dpl(2,1,2)=ANN(i)%y(5)*igeom(ncoord)
  dpl(3,1,2)=ANN(i)%y(6)
  !(2,1)
  dpl(:,2,1)=dpl(:,1,2)
  !(2,2)
  dpl(1,2,2)=ANN(i)%y(7)
  dpl(2,2,2)=ANN(i)%y(8)
  dpl(3,2,2)=ANN(i)%y(9)*igeom(ncoord)

  return
end subroutine EvaluateDpl0
!==================================================================================
subroutine EvaluateDpl(igeom,dpl)
  use nndms
  implicit none
  real*8, intent(in) :: igeom(ncoord)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  real*8 :: dpl0(3,nstates,nstates)
  integer :: i

  dpl=0.d0
  do i=1,na
    call EvaluateDpl0(i,igeom,dpl0)
    dpl=dpl+dpl0
  end do
  dpl=dpl/dble(na)

  return
end subroutine EvaluateDpl
!==================================================================================
subroutine Evaluatediadip(cgeom,dpl)
  use nndms
  implicit none
  integer, parameter :: lwork=99
  real*8, intent(in) :: cgeom(3,natoms)
  real*8, intent(out) :: dpl(3,nstates,nstates)
  real*8, allocatable :: gmek(:,:),coord(:),work(:)
  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: ckl0
  real*8, dimension(:,:,:), allocatable :: cg, dcg
  real*8 :: T(3,3)
  integer :: i,j,info

  allocate(gmek(3,natoms),coord(ncoord),work(lwork))
  allocate(ckl0(nstates,nstates))
  allocate(e(nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  gmek=cgeom
  call order(gmek,info)
  call orientation(gmek,T,info)
  !diabatic dipole
  call dplcoord(natoms,ncoord,gmek,coord)
  call EvaluateDpl(coord,dpl)

  !rotation
  do i=1,nstates
    do j=1,nstates
      dpl(:,i,j)=matmul(transpose(T),dpl(:,i,j))
    end do
  end do

  return
end
!=================================================================================
function triple_product(a,b,c)
  implicit none
  real*8 :: triple_product
  real*8, intent(in) :: a(3),b(3),c(3)
  triple_product = a(1)*b(2)*c(3) - a(1)*b(3)*c(2) - a(2)*b(1)*c(3) + &
                   a(2)*b(3)*c(1) + a(3)*b(1)*c(2) - a(3)*b(2)*c(1)
  return
end
!=================================================================================
subroutine OOP(a,b,c,d,coord,dcoord,id)
  implicit none
  real*8, intent(in) :: a(3),b(3),c(3),d(3)
  integer, intent(in) :: id
  real*8, intent(out) :: coord, dcoord(12)
  real*8, external :: triple_product
  real*8 :: ab(3), ac(3), ad(3), rab, rac, rad
  real*8 :: abx(12),acx(12),adx(12),dtp(12)

  ab=a-b
  ac=a-c
  ad=a-d
  rab=sqrt(dot_product(ab,ab))
  rac=sqrt(dot_product(ac,ac))
  rad=sqrt(dot_product(ad,ad))
  coord=triple_product(ab,ac,ad)/(rab*rac*rad)

  if(id.ne.0) then
    abx=0.d0
    abx(1:3)=ab/rab
    abx(4:6)=-abx(1:3)
    acx=0.d0
    acx(1:3)=ac/rac
    acx(7:9)=-acx(1:3)
    adx=0.d0
    adx(1:3)=ad/rad
    adx(10:12)=-adx(1:3)
    !tp = ab(1)*ac(2)*ad(3) - ab(1)*ac(3)*ad(2) - ab(2)*ac(1)*ad(3) + &
    !     ab(2)*ac(3)*ad(1) + ab(3)*ac(1)*ad(2) - ab(3)*ac(2)*ad(1)
    dtp(1)=ac(2)*ad(3)-ac(3)*ad(2)-ab(2)*ad(3)+ab(2)*ac(3)+ab(3)*ad(2)-ab(3)*ac(2)
    dtp(2)=ab(1)*ad(3)-ab(1)*ac(3)-ac(1)*ad(3)+ac(3)*ad(1)+ab(3)*ac(1)-ab(3)*ad(1)
    dtp(3)=ab(1)*ac(2)-ab(1)*ad(2)-ab(2)*ac(1)+ab(2)*ad(1)+ac(1)*ad(2)-ac(2)*ad(1)
    dtp(4)=-ac(2)*ad(3)+ac(3)*ad(2)
    dtp(5)=ac(1)*ad(3)-ac(3)*ad(1)
    dtp(6)=-ac(1)*ad(2)+ac(2)*ad(1)
    dtp(7)=ab(2)*ad(3)-ab(3)*ad(2)
    dtp(8)=-ab(1)*ad(3)+ab(3)*ad(1)
    dtp(9)=ab(1)*ad(2)-ab(2)*ad(1)
    dtp(10)=-ab(2)*ac(3)+ab(3)*ac(2)
    dtp(11)=ab(1)*ac(3)-ab(3)*ac(1)
    dtp(12)=-ab(1)*ac(2)+ab(2)*ac(1)
    dcoord=dtp/(rab*rac*rad)-coord*(abx/rab+acx/rac+adx/rad)
  end if

  return
end
!==================================================================================
