program training
use pipvariables
implicit none

! Number of points in training set to use
! (should be <= the actual maximum)
!integer,parameter :: Ntot=1000 !169824
integer,parameter :: Ntot=14982
!integer,parameter :: Ntot=63041

! The file with the training set
!character(len=*),parameter :: trainingsetfile = "trainingsets/BrCH5.set1.xyz"
character(len=*),parameter :: trainingsetfile = "trainingsets/BrClH2.setB2.xyz"
!character(len=*),parameter :: trainingsetfile = "trainingsets/CH5.set1.xyz"

! Conversion to internal units (eV)
real(kind=8),parameter :: ev=0.04336412 ! Energy is originally in kcal/mol
!real(kind=8),parameter :: ev=27.21138505d0 ! Energy is originally in kcal/mol

! Minimum energy to shift all energies by
!real*8,parameter :: vzero = -1639812.67919d0 ! Energy is originally in kcal/mol
real*8,parameter :: vzero = -1903000.871424484765d0 ! Energy is originally in kcal/mol
!real*8,parameter :: vzero = -40.95588433d0 ! Energy is originally in kcal/mol

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

!   Don't read forces in yet, since we cannot train on them
!   read(1,*) atom(j),x(:,j), f(:,j)

    read(1,*) atom(j),x(:,j)
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
! For this set, decide some random number of
! geometries to put in the training set
call scramble(Ntot,training_indices)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For ADAM and regular optimization:
! The "learning rate"
stepsize = (1.0d2 / (Ntot)) ! (1.0d-5 / (Ntot))
stepsize = 0.000002

! For ADAM optimization:
batch_size=Ntot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the LM optimization:
call LM_initNweights(Ntot,Nweights)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

lr_counter = 0
minError = 1.0d9
prevError = 0.0d0
do j = 1, 100

  call date_and_time(TIME=time)
  print *, ""
  print *, ""
  print *, "Iteration: ", j, " Time: ", time
  call net_initDW()

  ! For each batch, shuffle the training set
  ! so that the gradient is a bit more
  ! "stochastic"
  call scramble(Ntot,shuffled_indices)

  ! Keep track of the validation (non-training) set's error separately
  validationMAE  = 0.0d0
  validationRMSE = 0.0d0
  
  do i=1,Ntot
    k = shuffled_indices(i)

    ! For the LM update:
    ! If this point is part of the training set, then
    ! calculate both its energy and its derivative
    if (training_indices(k) <= Ntrain) then

      call LM_accumulate(p(:,k),ene(k),v(k),Nweights)

    ! Otherwise, just calculate its energy
    else
      call energy(p(:,k),v(k))
      validationMAE = validationMAE + abs(v(k)-ene(k))
      validationRMSE = validationRMSE + (v(k)-ene(k))**2
      v(k) = ene(k)
    end if
  
    ! For the regular and ADAM update:
!   call energyANDderivative(p(:,k),ene(k),v(k),dvdp)
!   if (modulo(i,batch_size)==0) then
!     !call net_printW()
!     if (i==batch_size) call net_printW()

!     !call net_update(stepsize)
!     call net_ADAMupdate(stepsize,j)
!     call net_initDW()
!   end if
!   if (modulo(i,batch_size)==0) print *, i
  enddo

  prevError = error
  error = sqrt(sum((v(1:Ntot)-ene(1:Ntot))**2)/Ntrain)
  
  print *, " MAE: ", sum(abs(v(1:Ntot)-ene(1:Ntot)))/Ntrain, " over ", Ntrain, " points"
  print *, "RMSE: ", error, " over ", Ntrain, " points"

  if (Ntot > Ntrain) then
    print *, "VALIDATION  MAE: ", validationMAE/(Ntot-Ntrain), " over ", Ntot-Ntrain, " points"
    print *, "VALIDATION RMSE: ", sqrt(validationRMSE/(Ntot-Ntrain)), " over ", Ntot-Ntrain, " points"
  end if
  
! call net_update(stepsize)
! call net_ADAMupdate(stepsize,j)
  call net_LMupdate(Ntrain*(error**2))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Save to a checkpoint every now and then, identifying
  ! the checkpoint with the step number
  if (modulo(j,5)==0) call net_write(j)

  ! Reduce the learning rate if the error has not gone
  ! down recently
  if (error < minError) then
    minError = error
    error_counter = 0
  else
    error_counter = error_counter + 1
  end if
  if (error_counter > 150 .or. lr_counter > 500) then
    lr_counter = 0
    error_counter = 0
    stepsize = stepsize * 0.50
    print *, "Learning rate decreased to:", stepsize
  else
    lr_counter = lr_counter + 1
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do

end program training

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
