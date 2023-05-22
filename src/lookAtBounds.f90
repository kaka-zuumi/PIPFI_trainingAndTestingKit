program lookAtBounds
use pipvariables
implicit none

! Number of points in training set to use
! (should be <= the actual maximum)
!integer,parameter :: Ntot=1000 !169824
integer,parameter :: Ntot=14982

! The file with the training set
!character(len=*),parameter :: trainingsetfile = "trainingsets/BrCH5.set1.xyz"
character(len=*),parameter :: trainingsetfile = "trainingsets/BrClH2.setB2.xyz"

! Conversion to internal units (eV)
real(kind=8),parameter :: ev=0.04336412 ! Energy is originally in kcal/mol

! Minimum energy to shift all energies by
!real*8,parameter :: vzero = -1639812.67919d0 ! Energy is originally in kcal/mol
real*8,parameter :: vzero = -1903000.871424484765d0 ! Energy is originally in kcal/mol

! The number of points to train with...
! the rest will be used for validation
integer :: Ntrain = (Ntot*3)/4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,parameter :: Nbonds = (Natoms*(Natoms-1))/2

character(len=2) :: atom(Natoms),ctmp
real*8 :: x(3,Natoms),r(Nbonds,Ntot),f(3,Natoms)
real*8 :: stepsize, error, prevError, minError
integer :: batch_size, lr_counter, error_counter
real*8 :: v(Ntot),ene(Ntot),p(Np,Ntot), tmpp(Np), dvdp(Np)
real*8 :: minmaxp(2,Np+1)
integer :: Nweights
integer :: natom,i,j,k,id,itmp
integer,dimension(Ntot) :: training_indices, shuffled_indices
real*8 :: validationMAE, validationRMSE
character(10) :: time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read in the energies and geometries
open(1,file=trainingsetfile,status='old')
do i=1,Ntot
  read(1,*) natom
  read(1,*) ene(i)
  ene(i)=(ene(i)-(vzero))*ev
  do j=1,natom
    read(1,*) atom(j),x(:,j), f(:,j)
  enddo

  call get_r(natom,x,r(:,i))
  call get_p(r(:,i),tmpp)
  call get_fp(Np,tmpp,p(:,i))

  if (modulo(i,Ntot/20)==0) then
    print *, "Percentage loaded: ", i*100.0/ntot, "%"
  end if
enddo
close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize by doing one calculation
! and reading in weights
call random_seed()
call energyANDderivative(p(:,1),ene(1),v(1),dvdp)
call net_initW()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call date_and_time(TIME=time)
print *, ""
print *, "Iteration: ", 0, " Time: ", time
call net_initDW()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update the minimum and maximum values of
! the polynomial inputs
do i = 1, Np
  minmaxp(1,i) = minval(p(i,1:ntot))
  minmaxp(2,i) = maxval(p(i,1:ntot))
end do
minmaxp(1,Np+1) = minval(ene(1:ntot))
minmaxp(2,Np+1) = maxval(ene(1:ntot))
call net_MINMAXupdate(minmaxp(1:2,:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save to a checkpoint the NN but with
! the lower and upper bounds
call net_write(0)

end program lookAtBounds

subroutine scramble(m,array)
implicit none
integer,intent(in) :: m
integer,dimension(m),intent(out) :: array
integer :: i,j,k,n, itemp
real :: u

do i = 1, m
  array(i) = i
end do

n=1
do k=1,2
  do i=1,m
    call random_number(u)
    j = n + floor((m+1-n)*u)
    itemp=array(j)
    array(j)=array(i)
    array(i)=itemp
  enddo
enddo
end subroutine scramble
