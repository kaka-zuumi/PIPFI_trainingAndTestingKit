program testEnergyAndOrForces
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
real*8 :: x(3,Natoms),r(Nbonds,Ntot),f(Ntot,3,Natoms)
real*8 :: stepsize, error, prevError, minError
integer :: batch_size, lr_counter, error_counter
real*8 :: v(Ntot),ene(Ntot),p(Np,Ntot),dvdx(Ntot,3,Natoms), dv(Ntot,Np,3,Natoms)
real*8 :: tmpr, tmprs(Nbonds), tmpv, tmpp(Np), tmpdp(Np)
real*8 :: dr(Nbonds,3,Natoms), dp(Np,Nbonds), dvdp(Np)
integer :: Nweights
integer :: natom,i,j,k,n,id,itmp
integer,dimension(Ntot) :: training_indices, shuffled_indices
real*8 :: validationMAE, validationRMSE

logical :: calculate_forces = .true.
real*8,parameter :: inverseDistanceIsZeroThreshhold = dexp(-(20.0d0**2))
real(kind=8),external :: ddot
real*8 :: fmae, frmse
character(10) :: time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read in the energies and geometries
open(1,file=trainingsetfile,status='old')
do i=1,Ntot
  read(1,*) natom
  read(1,*) ene(i)
  ene(i)=(ene(i)-(vzero))*ev
  do j=1,Natoms
    read(1,*) atom(j),x(:,j), f(i,:,j)
  enddo

  if (calculate_forces) then
    call get_rANDdr(Natoms,x,r(:,i),dr)
    call get_p(r(:,i),tmpp)
    call get_fpANDdfp(Np,tmpp,p(:,i),tmpdp)

    f(i,1:3,1:Natoms)=f(i,1:3,1:Natoms)*ev
    call get_dpdr(r(:,i),dp(1:Np,:))
    do k=1,Natoms
      do id=1,3
        do j=1,Np
          dv(i,j,id,k) = ddot(Nbonds,dr(:,id,k),1,dp(j,:),1)*tmpdp(j)

          ! New part to address very LARGE distances
          if (tmpp(j) <= inverseDistanceIsZeroThreshhold) then
            dv(i,j,id,k) = 0.0d0
          end if
        end do
      end do
    end do
  else
    call get_r(Natoms,x,r(:,i))
    call get_p(r(:,i),tmpp)
    call get_fp(Np,tmpp,p(:,i))
  end if

  if (modulo(i,Ntot/20)==0) then
    print *, "Percentage loaded: ", i*100.0/Ntot, "%"
  end if
enddo
close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize by doing one calculation
! and reading in weights
call random_seed()
call energyANDderivative(p(:,1),ene(1),v(1),dvdp(1:Np))
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
do j = 1, 1

  call date_and_time(TIME=time)
  print *, ""
  print *, ""
  print *, "Iteration: ", j, " Time: ", time
  call net_initDW()

  ! For each batch, shuffle the training set
  ! so that the gradient is a bit more
  ! "stochastic"
  call scramble(Ntot,shuffled_indices)

  batch_size=Ntot
  do i=1,Ntot
    k = shuffled_indices(i)
    call energyANDderivative(p(:,k),ene(k),v(k),dvdp(1:Np))
    if (calculate_forces) then
      dvdx(k,1:3,1:Natoms) = 0.0d0
      do id=1,Np
        dvdx(k,1:3,1:Natoms) = dvdx(k,1:3,1:Natoms) + &
                         dvdp(id)*dv(k,id,1:3,1:Natoms)
      end do

!     if (.false.) then
      if (modulo(i,10)==0) then
      print *, "geometry", k
      print *, "TRUE F, PREDICTED F, error F"
      do n=1,Natoms
        write(6,FMT="(3(3(F8.4,1x),2x))") f(k,1:3,n), dvdx(k,1:3,n), f(k,1:3,n)-dvdx(k,1:3,n)
      end do

      if (k==Ntot) then
        print *, "numerical forces:"
        call get_r(Natoms,x,tmprs)
        call get_p(tmprs,tmpp)
        call get_fp(Np,tmpp,tmpdp)
        call energy(tmpdp,v(k))
        do n=1,Natoms
          do id=1,3
            x(id,n)=x(id,n)+0.000005d0
            call get_r(Natoms,x,tmprs)
            call get_p(tmprs,tmpp)
            call get_fp(Np,tmpp,tmpdp)
            call energy(tmpdp,tmpv)
            dvdx(k,id,n) = -(tmpv-v(k))/0.000005d0
            x(id,n)=x(id,n)-0.000005d0
          end do
          print *, dvdx(k,1:3,n)
      end do
      end if

      end if
    end if

    if (modulo(i,batch_size)==0) then
      if (i==batch_size) call net_printW()

      !call net_update(stepsize)
      call net_ADAMupdate(stepsize,j)
      call net_initDW()
    end if
  enddo

  prevError = error
  error = sqrt(sum((v(1:Ntot)-ene(1:Ntot))**2)/Ntot)

  write(6,FMT="(A,F12.6,A,I0,A)") &
      "E(eV)    MAE: ", sum(abs(v(1:Ntot)-ene(1:Ntot)))/Ntot, " over ", Ntot, " points"
  write(6,FMT="(A,F12.6,A,I0,A)") &
      "E(eV)   RMSE: ", error, " over ", Ntot, " points"

  ! When measuring the force error, flip the sign of the
  ! energy gradient before taking the difference
  if (calculate_forces) then
    fmae = sum(abs(f(1:Ntot,1:3,1:Natoms)-dvdx(1:Ntot,1:3,1:Natoms)))/Ntot
    frmse = sqrt(sum((f(1:Ntot,1:3,1:Natoms)-dvdx(1:Ntot,1:3,1:Natoms))**2)/Ntot)
    write(6,FMT="(A,F12.6,A,I0,A,F12.6,A)") &
             "F(eV/A)  MAE: ", fmae, " over ", Ntot, " points   (",&
             fmae/(3*Natoms), " per 3N xyzs)"
    write(6,FMT="(A,F12.6,A,I0,A,F12.6,A)") &
             "F(eV/A) RMSE: ", frmse, " over ", Ntot, " points   (",&
             frmse/sqrt(3.0*Natoms), " per 3N xyzs)"
  end if

! call net_update(stepsize)
  call net_ADAMupdate(stepsize,j)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Save to a checkpoint every now and then, identifying
  ! the checkpoint with the step number
  if (modulo(j,50)==0) call net_write(j)

  ! Reduce the learning rate if the error has not gone
  ! down recently
  if (error < minError) then
    minError = error
    error_counter = 0
  else
    error_counter = error_counter + 1
  end if
  if (error_counter > 50 .or. lr_counter > 200) then
    lr_counter = 0
    error_counter = 0
    stepsize = stepsize * 0.50
    print *, "Learning rate decreased to:", stepsize
  else
    lr_counter = lr_counter + 1
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do


end program testEnergyAndOrForces

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
