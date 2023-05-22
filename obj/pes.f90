module net
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! number of networks for average
integer, parameter :: nfits=3
! number of hidden layers
integer, parameter :: max_nhid=2
! number of the summation of every layer's outputs
integer, parameter :: max_ny=900
! dimension of input data
integer, parameter :: max_ndim=600
! number of weights
integer, parameter :: max_nw=99999
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: init_flag=0
integer :: nhid(nfits)
integer :: nl(0:max_nhid+1,nfits)
integer :: ny(nfits)
integer :: nw(nfits)
integer :: idyf(2,max_ny,nfits),idwf(2,max_ny,nfits)
real(kind=8) :: minmax(2,max_ndim+1,nfits)
real(kind=8) :: w(max_nw,nfits), wprev(max_nw,nfits)
real(kind=8) :: dw(max_nw,nfits)
end module net
subroutine net_read
use net
implicit none
integer :: i,j,k,ii

init_flag=1

! Write a single NN to a single file...
! if there were more, this could be looped
ii=1

idyf=0; idwf=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(120,file="net.1",status="old")
read(120,*) nhid(ii)
if(nhid(ii).gt.max_nhid) stop "max_nhid..."
ny(ii)=0
do i=0,nhid(ii)
  read(120,*) nl(i,ii)
  ny(ii)=ny(ii)+nl(i,ii)+1
end do
if(nl(0,ii).gt.max_ndim) stop "max_ndim..."
nl(nhid(ii)+1,ii)=1; ny(ii)=ny(ii)+1
read(120,*) nw(ii)
if(nw(ii).gt.max_nw) stop "max_nw..."
do i=1,nl(0,ii)+1
  read(120,*) minmax(:,i,ii)
end do
print *, "We see that nfits is: ", nfits
print *, "We see that nl(0,ii) is: ", nl(0,ii)
print *, "We see that nw is: ", nw(ii), "..."
k=nl(0,ii)+1
do i=1,nhid(ii)+1
  do j=1,nl(i,ii)
    k=k+1
    if(k.gt.max_ny) stop "max_ny..."
    read(120,*) idyf(1:2,k,ii),idwf(1:2,k,ii)
  end do
  k=k+1
end do
do i=1,nw(ii)
  read(120,*) w(i,ii)
end do
close(120)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine net_read

subroutine net_write(fileID)
use net
implicit none
integer, intent(in) :: fileID
integer :: i,j,k,ii
character(20) :: intText

! Name each output file by its fileID
write(intText,FMT=*) fileID

! Write a single NN to a single file...
! if there were more, this could be looped
ii=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(120,file="checkpoints/net.1-"//trim(adjustl(intText)))
write(120,*) nhid(ii)
if(nhid(ii).gt.max_nhid) stop "max_nhid..."
ny(ii)=0
do i=0,nhid(ii)
  write(120,*) nl(i,ii)
  ny(ii)=ny(ii)+nl(i,ii)+1
end do
if(nl(0,ii).gt.max_ndim) stop "max_ndim..."
nl(nhid(ii)+1,ii)=1; ny(ii)=ny(ii)+1
write(120,*) nw(ii)
if(nw(ii).gt.max_nw) stop "max_nw..."
do i=1,nl(0,ii)+1
  write(120,*) minmax(:,i,ii)
end do
k=nl(0,ii)+1
do i=1,nhid(ii)+1
  do j=1,nl(i,ii)
    k=k+1
    if(k.gt.max_ny) stop "max_ny..."
    write(120,*) idyf(1:2,k,ii),idwf(1:2,k,ii)
  end do
  k=k+1
end do
do i=1,nw(ii)
  write(120,*) w(i,ii)
end do
close(120)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine net_write

! We can identify and update the upper and lower
! bounds for the polynomials for a new training set
subroutine net_MINMAXupdate(new_minmax)
use net
implicit none
real(kind=8),intent(in) :: new_minmax(2,nl(0,1)+1)
minmax(1:2,1:nl(0,1)+1,1) = &
  new_minmax(1:2,1:nl(0,1)+1)
end subroutine net_MINMAXupdate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine forward(flag,p,y)
use net
implicit none
integer,intent(in) :: flag
real(kind=8),intent(in) :: p(nl(0,1))
real(kind=8),intent(out) :: y(max_ny)
real(kind=8),external :: ddot
real(kind=8) :: v
integer :: i,j,k,nyx
integer :: inc
nyx=0
do i=0,nhid(flag)
  nyx=nyx+nl(i,flag)+1
  y(nyx)=1.d0
end do
do i=1,nl(0,flag)
  y(i)=2.d0*(p(i)-minmax(1,i,flag))/(minmax(2,i,flag)-minmax(1,i,flag))-1.d0
end do
k=nl(0,flag)+1
do i=1,nhid(flag)+1
  inc=nl(i-1,flag)+1
  do j=1,nl(i,flag)
    k=k+1
    v=ddot(inc,w(idwf(1,k,flag):idwf(2,k,flag),flag),1,y(idyf(1,k,flag):idyf(2,k,flag)),1)
    if(k.lt.ny(flag)) then
      y(k)=v/sqrt(1.d0+v**2)
    else if(k.eq.ny(flag)) then
      y(k)=v
    end if
  end do
  k=k+1
end do
y(ny(flag))=(y(ny(flag))+1.d0)*(minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0+minmax(1,nl(0,flag)+1,flag)
end subroutine forward

subroutine forwardANDbackward(flag,p,truey,y,dvdx)
use net
implicit none
integer,intent(in) :: flag
real(kind=8),intent(in) :: p(nl(0,1))
real(kind=8),intent(in) :: truey
real(kind=8),intent(out) :: y(max_ny)
real(kind=8),intent(out) :: dvdx(nl(0,1))
real(kind=8) :: dy(max_ny)
real(kind=8) :: protoF(nhid(flag),maxval(nl(0:nhid(flag)+1,flag)))
real(kind=8) :: denominator
real(kind=8),external :: ddot
real(kind=8) :: tmpdw(max_nw)
real(kind=8) :: v
integer :: i,j,k,nyx
integer :: inc
nyx=0
do i=0,nhid(flag)
  nyx=nyx+nl(i,flag)+1
  y(nyx)=1.d0
end do
do i=1,nl(0,flag)
  y(i)=2.d0*(p(i)-minmax(1,i,flag))/(minmax(2,i,flag)-minmax(1,i,flag))-1.d0
end do

k=nl(0,flag)+1
do i=1,nhid(flag)+1
  inc=nl(i-1,flag)+1
  do j=1,nl(i,flag)
    k=k+1
    v=ddot(inc,w(idwf(1,k,flag):idwf(2,k,flag),flag),1,y(idyf(1,k,flag):idyf(2,k,flag)),1)
    if(k.lt.ny(flag)) then
      protoF(i,j)=1.0d0/sqrt(1.d0+v**2)
      y(k)=v*protoF(i,j)
    else if(k.eq.ny(flag)) then
      y(k)=v
    end if

  end do
  k=k+1
end do
y(ny(flag))=(y(ny(flag))+1.d0)*(minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0+minmax(1,nl(0,flag)+1,flag)

! The final layer's "delta" is just the error
!dy(ny(flag)) = -(y(ny(flag)) - truey)*2

! Account for scaling at the final layer
!dy(ny(flag)) = dy(ny(flag)) * (minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0

dy(ny(flag)) = -(minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0

tmpdw = 0.0d0

!do i=1,nhid(flag)+1
do i=nhid(flag)+1,1,-1
  k=k-1
  dy(idyf(1,k,flag):idyf(2,k,flag)) = 0.0d0

  inc=nl(i-1,flag)+1
  do j=nl(i,flag),1,-1

    ! In this neural network, no activation is applied on the
    ! final layer, so ignore that derivative
    if(k.lt.ny(flag)) then
      dy(k)=dy(k)*(protoF(i,j)**3)
    else if(k.eq.ny(flag)) then
      dy(k)=dy(k)
    end if

    k=k-1
  end do

  k=k+nl(i,flag)
  do j=nl(i,flag),1,-1

    ! In "back propagation" the previous layer's "delta" is related
    ! by a matrix multiplication of the TRANSPOSE of the previous
    ! layer's weights and the current layer's "delta"
    dy(idyf(1,k,flag):idyf(2,k,flag)) = dy(idyf(1,k,flag):idyf(2,k,flag)) + &
        dy(k) * w(idwf(1,k,flag):idwf(2,k,flag),flag)

    ! The change in the parameters (w here has both the weights and
    ! biases) are directly related to the weights and the "delta"
!   dw(idwf(1,k,flag):idwf(2,k,flag),flag) = dw(idwf(1,k,flag):idwf(2,k,flag),flag) + &
!       dy(k) * y(idyf(1,k,flag):idyf(2,k,flag))
    tmpdw(idwf(1,k,flag):idwf(2,k,flag)) = &
        dy(k) * y(idyf(1,k,flag):idyf(2,k,flag))

    k=k-1
  end do

end do

!All weight derivatives must also be scaled by
!the error, and then they can be added to the
!total weight derivative dw 
dw(:,flag) = dw(:,flag) + tmpdw * (y(ny(flag)) - truey)*2

!Account for the scaling at the initial layer
do i=1,nl(0,flag)
  dy(i) = dy(i)*2.d0/(minmax(2,i,flag)-minmax(1,i,flag))
end do

!The derivative with respect to the inputs is
!simply the first layer's "delta"  
dvdx(1:nl(0,1)) = dy(1:nl(0,1))

end subroutine forwardANDbackward
subroutine forwardANDbackwardNOerror(flag,p,y)
use net
implicit none
integer,intent(in) :: flag
real(kind=8),intent(in) :: p(nl(0,1))
real(kind=8),intent(out) :: y(max_ny)
real(kind=8) :: dy(max_ny)
real(kind=8) :: protoF(nhid(flag),maxval(nl(0:nhid(flag)+1,flag)))
real(kind=8) :: denominator
real(kind=8),external :: ddot
real(kind=8) :: v
integer :: i,j,k,nyx
integer :: inc
nyx=0
do i=0,nhid(flag)
  nyx=nyx+nl(i,flag)+1
  y(nyx)=1.d0
end do
do i=1,nl(0,flag)
  y(i)=2.d0*(p(i)-minmax(1,i,flag))/(minmax(2,i,flag)-minmax(1,i,flag))-1.d0
end do

k=nl(0,flag)+1
do i=1,nhid(flag)+1
  inc=nl(i-1,flag)+1
  do j=1,nl(i,flag)
    k=k+1
    v=ddot(inc,w(idwf(1,k,flag):idwf(2,k,flag),flag),1,y(idyf(1,k,flag):idyf(2,k,flag)),1)
    if(k.lt.ny(flag)) then
      protoF(i,j)=1.0d0/sqrt(1.d0+v**2)
      y(k)=v*protoF(i,j)
    else if(k.eq.ny(flag)) then
      y(k)=v
    end if

  end do
  k=k+1
end do
y(ny(flag))=(y(ny(flag))+1.d0)*(minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0+minmax(1,nl(0,flag)+1,flag)

! Account for scaling at the final layer
dy(ny(flag)) = (minmax(2,nl(0,flag)+1,flag)-minmax(1,nl(0,flag)+1,flag))/2.d0

do i=nhid(flag)+1,1,-1
  k=k-1
  dy(idyf(1,k,flag):idyf(2,k,flag)) = 0.0d0

  inc=nl(i-1,flag)+1
  do j=nl(i,flag),1,-1

    ! In this neural network, no activation is applied on the
    ! final layer, so ignore that derivative
    if(k.lt.ny(flag)) then
      dy(k)=dy(k)*(protoF(i,j)**3)
    else if(k.eq.ny(flag)) then
      dy(k)=dy(k)
    end if

    k=k-1
  end do

  k=k+nl(i,flag)
  do j=nl(i,flag),1,-1

    ! In "back propagation" the previous layer's "delta" is related
    ! by a matrix multiplication of the TRANSPOSE of the previous
    ! layer's weights and the current layer's "delta"
    dy(idyf(1,k,flag):idyf(2,k,flag)) = dy(idyf(1,k,flag):idyf(2,k,flag)) + &
        dy(k) * w(idwf(1,k,flag):idwf(2,k,flag),flag)

    ! The change in the parameters (w here has both the weights and
    ! biases) are directly related to the weights and the "delta"
    dw(idwf(1,k,flag):idwf(2,k,flag),flag) = &
        dy(k) * y(idyf(1,k,flag):idyf(2,k,flag))

    k=k-1
  end do

end do

!Account for the scaling at the initial layer
do i=1,nl(0,flag)
  dy(i) = dy(i)*2.d0/(minmax(2,i,flag)-minmax(1,i,flag))
end do

end subroutine forwardANDbackwardNOerror

subroutine energy(x,v)
use net
implicit none
real(kind=8), intent(in) :: x(max_ndim)
real(kind=8), intent(out) :: v
real(kind=8) :: y(max_ny)
integer :: i

! If this is the first call, then
! initialize the NN
if(init_flag.eq.0) then
  call net_read
end if

call forward(1,x,y)
v = y(ny(1))

end subroutine energy

subroutine energyANDderivative(x,vtrue,v,dvdx)
use net
implicit none
real(kind=8), intent(in) :: x(max_ndim)
real(kind=8), intent(in) :: vtrue
real(kind=8), intent(out) :: v
real(kind=8), intent(out) :: dvdx(nl(0,1))
real(kind=8) :: vF, vB
real(kind=8) :: y(max_ny)
real(kind=8) :: deltaW
character(10) :: time
integer :: i

! If this is the first call, then
! initialize the NN
if(init_flag.eq.0) then
  call net_read
end if

v=0.d0
call forwardANDbackward(1,x,vtrue,y,dvdx)
v = y(ny(1))

end subroutine energyANDderivative

subroutine get_r(Natoms,x,r)
implicit none
integer,intent(in) :: Natoms
double precision,dimension(3,Natoms),intent(in) :: x
double precision,dimension((Natoms*(Natoms-1))/2),intent(out) :: r
double precision :: tmpr
integer :: i,j,k, id
id=0
do j=1,Natoms-1
  do k=j+1,Natoms
    id=id+1
    tmpr=dsqrt(sum((x(:,j)-x(:,k))**2))
    r(id)=dexp(-tmpr*1e0)
  enddo
enddo
end subroutine get_r

subroutine get_rANDdr(Natoms,x,r,dr)
implicit none
integer,intent(in) :: Natoms
real*8,dimension(3,Natoms),intent(in) :: x
real*8,dimension((Natoms*(Natoms-1))/2),intent(out) :: r
real*8,dimension((Natoms*(Natoms-1))/2,3,Natoms),intent(out) :: dr
real*8 :: tmpr
integer :: i,j,k, id
dr = 0.0d0
id=0
do j=1,Natoms-1
  do k=j+1,Natoms
    id=id+1
    tmpr=dsqrt(sum((x(:,j)-x(:,k))**2))
    r(id)=dexp(-tmpr*1e0)

    dr(id,1:3,j) = -(r(id)/tmpr)*(x(:,j)-x(:,k))
    dr(id,1:3,k) = -dr(id,1:3,j)
  enddo
enddo
end subroutine get_rANDdr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine net_initW
use net
implicit none
integer :: i,j,k,ii
logical :: lexist

! Check if the restart file exists and
! if so, read the weights from it
inquire(file="net.1.restart",exist=lexist)
if (lexist) then
  ii = 1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(120,file="net.1.restart",status="old")
  read(120,*) nhid(ii)
  if(nhid(ii).gt.max_nhid) stop "max_nhid..."
  ny(ii)=0
  do i=0,nhid(ii)
    read(120,*) nl(i,ii)
    ny(ii)=ny(ii)+nl(i,ii)+1
  end do
  if(nl(0,ii).gt.max_ndim) stop "max_ndim..."
  nl(nhid(ii)+1,ii)=1; ny(ii)=ny(ii)+1
  read(120,*) nw(ii)
  if(nw(ii).gt.max_nw) stop "max_nw..."
  do i=1,nl(0,ii)+1
    read(120,*) minmax(:,i,ii)
  end do
  print *, "We see that nfits is: ", nfits
  print *, "We see that nl(0,ii) is: ", nl(0,ii)
  print *, "We see that nw is: ", nw(ii), "..."
  k=nl(0,ii)+1
  do i=1,nhid(ii)+1
    do j=1,nl(i,ii)
      k=k+1
      if(k.gt.max_ny) stop "max_ny..."
      read(120,*) idyf(1:2,k,ii),idwf(1:2,k,ii)
    end do
    k=k+1
  end do
  do i=1,nw(ii)
    read(120,*) w(i,ii)
  end do
  close(120)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! If no restart file exist, just set the
! weights to something reasonable (I usually
! do something ~ 1/#neurons)
else

  print *, "We're initializing the weights from random numbers!"
  do i=1,nw(1)
    call random_number(w(i,1))
  end do
  w = w * 0.02d0
end if

end subroutine net_initW

subroutine net_initDW
use net
implicit none

dw = 0.0d0

end subroutine net_initDW

subroutine net_printW
use net
implicit none

!print *, " w:",  w(1:5,1), "..."
!print *, "dw:", dw(1:5,1), "..."
print *, " w:",  w(maxloc(dw(1:nw(1),1)),1), "dw:",  maxval(dw(1:nw(1),1)), " (max) ..."

end subroutine net_printW

subroutine net_update(lr)
use net
implicit none
real(kind=8), intent(in) :: lr

! Update the weights
w(1:nw(1),1) = w(1:nw(1),1) - dw(1:nw(1),1) * lr

end subroutine net_update



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ADAM
use net
implicit none

integer :: ADAMinit_flag=0
integer :: prev_t = 0
real(kind=8),parameter :: ADAMbeta1 = 0.999d0
real(kind=8),parameter :: ADAMbeta2 = 0.9999d0
real(kind=8),parameter :: ADAMepsilon = 1.0d-8
real(kind=8) :: m(max_nw,nfits)
real(kind=8) :: v

real(kind=8) :: ADAMbeta1_var = 1.0d0
real(kind=8) :: ADAMbeta2_var = 1.0d0
end module ADAM
subroutine net_ADAMupdate(lr,t)
use ADAM
implicit none
real(kind=8), intent(in) :: lr
integer, intent(in) :: t

! Initialize the "momenta"
if(ADAMinit_flag.eq.0) then
  m(1:nw(1),1) = 0.0d0
  v = 0.0d0
  ADAMinit_flag = 1
end if

! Update some ADAM variables
if (t > prev_t) then
  ADAMbeta1_var = ADAMbeta1_var * ADAMbeta1
  ADAMbeta2_var = ADAMbeta2_var * ADAMbeta2
  prev_t = t
end if

! Update the "momenta"
m(1:nw(1),1) = m(1:nw(1),1) * ADAMbeta1 + &
               dw(1:nw(1),1) * (1.0d0 - ADAMbeta1)
v = v * ADAMbeta2 + &
    sum(dw(1:nw(1),1)**2) * (1.0d0 - ADAMbeta2)

! Update the weights
w(1:nw(1),1) = w(1:nw(1),1) + m(1:nw(1),1) * lr / &
               ((1.0d0 - ADAMbeta1_var)*&
                sqrt(sqrt(v/(1.0d0 - ADAMbeta2_var))+ADAMepsilon))

end subroutine net_ADAMupdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module LM
use net
implicit none

integer :: LMinit_flag=0
real(kind=8) :: LMlambda = 1.0d4
real(kind=8),parameter :: LMlambdaMIN = 1.0d-7
real(kind=8),parameter :: LMlambdaMAX = 1.0d7
real(kind=8),parameter :: LMlambdareduction = 1.0d0/9.0d0
real(kind=8),parameter :: LMlambdaincrease = 11.0d0

real(kind=8),parameter :: LMepsilon = 1.0d-8

real(kind=8),allocatable :: LM_LHS(:,:)
real(kind=8),allocatable :: LM_RHS(:)

real(kind=8) :: hprev(max_nw)
real(kind=8),allocatable :: LM_LHSprev(:,:)
real(kind=8),allocatable :: LM_RHSprev(:)
real(kind=8) :: RMSEprev
end module LM

subroutine get_rho_denominator(n,h,Jdiag,lambda,rho)
use LM
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: h(n)
real(kind=8), intent(in) :: Jdiag(n)
real(kind=8), intent(in) :: lambda
real(kind=8), intent(out) :: rho

rho = 1.0d0 / abs(dot_product(h(1:n),&
      lambda*Jdiag*h(1:n) +&
      LM_RHSprev(1:n)))

end subroutine get_rho_denominator
subroutine get_hLM(n,h)
use LM
implicit none
integer, intent(in) :: n
real(kind=8), intent(out) :: h(n)

real(kind=8),dimension(n,n) :: A
real(kind=8),dimension(n,1) :: b

real(kind=8),dimension(n,n) :: AF
real(kind=8) :: rcond
real(kind=8), dimension(1) :: FERR, BERR
real(kind=8) :: work(n**2)
integer, dimension(n) :: IWORK, IPIV
integer :: LWORK, INFO

! Make copies because LAPACK overwrites
! the input matrices
A = LM_LHSprev(1:n,1:n)
b = reshape(LM_RHSprev(1:n),(/n,1/))

! This solving of Ax=b would be pretty
! hard to do without LAPACK:
LWORK = n * n
call DSYSVX( "N", "U", n, 1, &
             A,&
             n, AF(1:n,1:n), &
             n, IPIV(1:n), &
             b,&
             n, h(1:n), n,&
             RCOND, FERR, BERR, WORK, LWORK,&
             IWORK(1:n), INFO )

end subroutine get_hLM
subroutine net_LMupdate(RMSE)
use LM
implicit none
real(kind=8), intent(in) :: RMSE
real(kind=8) :: h(max_nw)
real(kind=8) :: Jdiag(max_nw)
real(kind=8) :: rho, rho_denominator
character(10) :: time
integer :: i

call date_and_time(TIME=time)
print *, "LM step start     Time: ", time

! Initialize the LM algorithm
if(LMinit_flag.eq.0) then
  LMinit_flag = 1

  LM_LHSprev = LM_LHS
  LM_RHSprev = LM_RHS
  do i = 1, nw(1)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) * (1.0d0+LMlambda)
  end do
  call get_hLM(nw(1),h(1:nw(1)))
  do i = 1, nw(1)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) / (1.0d0+LMlambda)
  end do

  ! Always accept the zero-th step
  rho = LMepsilon + 1.0
else

  ! Calculate "rho", the goodness of fit
  ! of the new step in minimizing the RMSE
  do i = 1, nw(1)
    Jdiag(i) = LM_LHSprev(i,i)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) * (1.0d0+LMlambda)
  end do

  call get_hLM(nw(1),h(1:nw(1)))
  call get_rho_denominator(nw(1),hprev(1:nw(1)),Jdiag(1:nw(1)),LMlambda,&
          rho_denominator)

! print *, "|Ax-b|: ", sum((&
! matmul(LM_LHSprev(1:nw(1),1:nw(1)),&
!   reshape(h(1:nw(1)),(/nw(1),1/)))-&
!   reshape(LM_RHSprev(1:nw(1)),(/nw(1),1/)))**2)

  rho = (RMSEprev - RMSE) * rho_denominator
end if

! If the previous step was good, then just
! accept it
if (rho > LMepsilon) then
  call date_and_time(TIME=time)
  print *, "LM step has ... been accepted!   Time: ", time

  ! Erase history of the previous steps
  wprev = w
  hprev = h
  LM_LHSprev = LM_LHS
  LM_RHSprev = LM_RHS
  RMSEprev = RMSE

  ! And also decrease lambda
  LMlambda = max(LMlambda*LMlambdareduction,LMlambdaMIN)

! Otherwise, go back a step
else
  call date_and_time(TIME=time)
  print *, "LM step has NOT been accepted!    Time: ", time

  ! Retrieve the original LHS matrix so
  ! that a new step h can be calculated
  do i = 1, nw(1)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) / (1.0d0+LMlambda)
  end do

  ! Retrieve the original weights
  w(1:nw(1),1) = wprev(1:nw(1),1)

  ! And also increase lambda
  LMlambda = min(LMlambda*LMlambdaincrease,LMlambdaMAX)

  ! Calculate the new step h with the new
  ! lambda value
  do i = 1, nw(1)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) * (1.0d0+LMlambda)
  end do
  call get_hLM(nw(1),h(1:nw(1)))
  do i = 1, nw(1)
    LM_LHSprev(i,i) = LM_LHSprev(i,i) / (1.0d0+LMlambda)
  end do
end if

! Update the weights
w(1:nw(1),1) = w(1:nw(1),1) + h(1:nw(1))

! Zero out the LHS and RHS
LM_LHS = 0.0d0
LM_RHS = 0.0d0

end subroutine net_LMupdate
subroutine LM_accumulate(x,vtrue,v,Nweights) !,Ntrain,Nweights,dvs)
use net
use LM
implicit none
real(kind=8), intent(in) :: x(max_ndim)
real(kind=8), intent(in) :: vtrue
real(kind=8), intent(out) :: v
!integer, intent(in) :: Ntrain, Nweights
integer, intent(in) :: Nweights
!real(kind=8), intent(inout) :: dvs(Ntrain,Nweights)
real(kind=8) :: y(max_ny)
integer :: i,j
external :: daxpy, dger

! For each training set point, calculate
! the energy and its derivative
! (WITHOUT the (Y - Ytrue) contribution)
call forwardANDbackwardNOerror(1,x,y)
v = y(ny(1))

! Add the point's contribution to the LHS
call dger(Nweights,Nweights,1.0d0,dw(1:Nweights,1),1,dw(1:Nweights,1),1,LM_LHS(1:Nweights,1:Nweights),Nweights)

! Without LAPACK: 
!do i = 1, Nweights
!  do j = i, Nweights
!    LM_LHS(i,j) = LM_LHS(i,j) + dw(i,1)*dw(j,1)
!    LM_LHS(j,i) = LM_LHS(i,j)
!  end do
!end do

! Add the point's contribution to the RHS
call daxpy(Nweights,vtrue-v,dw(1:Nweights,1),1,LM_RHS(1:Nweights),1)

! Without LAPACK: 
!LM_RHS(1:Nweights) = LM_RHS(1:Nweights) +&
!    dw(1:Nweights,1) * (vtrue - v)

end subroutine LM_accumulate
subroutine LM_initNweights(n,Nweights)
use LM
implicit none
integer, intent(in) :: n
integer, intent(out) :: Nweights

Nweights = nw(1)

allocate(LM_LHS(Nweights,Nweights),&
         LM_LHSprev(Nweights,Nweights))
allocate(LM_RHS(Nweights),&
         LM_RHSprev(Nweights))
end subroutine LM_initNweights
