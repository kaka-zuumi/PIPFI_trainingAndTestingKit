subroutine get_dpdr(r,dp)
implicit none
real(kind=8),intent(in) :: r(6)
real(kind=8),intent(out) :: dp(7,6)
real(kind=8) :: s2(6)
s2=r**2
dp(1,1)=0.0d0
dp(1,2)=0.0d0
dp(1,3)=0.0d0
dp(1,4)=1.0d0
dp(1,5)=1.0d0
dp(1,6)=0.0d0
dp(2,1)=0.0d0
dp(2,2)=1.0d0
dp(2,3)=1.0d0
dp(2,4)=0.0d0
dp(2,5)=0.0d0
dp(2,6)=0.0d0
dp(3,1)=0.0d0
dp(3,2)=0.0d0
dp(3,3)=0.0d0
dp(3,4)=2*r(4)
dp(3,5)=2*r(5)
dp(3,6)=0.0d0
dp(4,1)=0.0d0
dp(4,2)=2*r(2)
dp(4,3)=2*r(3)
dp(4,4)=0.0d0
dp(4,5)=0.0d0
dp(4,6)=0.0d0
dp(5,1)=0.0d0
dp(5,2)=r(4)
dp(5,3)=r(5)
dp(5,4)=r(2)
dp(5,5)=r(3)
dp(5,6)=0.0d0
dp(6,1)=0.0d0
dp(6,2)=0.0d0
dp(6,3)=0.0d0
dp(6,4)=0.0d0
dp(6,5)=0.0d0
dp(6,6)=1.0d0
dp(7,1)=1.0d0
dp(7,2)=0.0d0
dp(7,3)=0.0d0
dp(7,4)=0.0d0
dp(7,5)=0.0d0
dp(7,6)=0.0d0
end subroutine get_dpdr


subroutine get_fpANDdfp(Np,p,fp,dfp)
use net
implicit none
integer,intent(in) :: Np
double precision,dimension(Np),intent(in) :: p
double precision,dimension(Np),intent(out) :: fp
double precision,dimension(Np),intent(out) :: dfp
logical :: flagOutOfBounds
double precision,parameter :: zeroLimit = 1.0d-44 ! ~ e^-100
integer :: i

fp(1:2)   = p(1:2)
fp(3:5)   = p(3:5)**(0.5d0)
fp(6:7)   = p(6:7)

! In MD simulations, check if the polynomial
! is starting to get out of range (0,1)
flagOutOfBounds = .false.
do i = 1,Np
  if (fp(i)>minmax(2,i,1)) then
!   write(6,FMT="(A,I0,A)") "PIP-NN WARNING: fp(", i, ") too large"
    flagOutOfBounds = .true.

  else if (fp(i)<minmax(1,i,1)) then
!   write(6,FMT="(A,I0,A)") "PIP-NN WARNING: fp(", i, ") too small"
    flagOutOfBounds = .true.
  end if

  ! If all is well, just compute the
  ! derivative as normal
  if (p(i)>zeroLimit) then
    dfp(i) = fp(i)/p(i)
  else
    dfp(i) = 0.0d0
  end if

end do

dfp(1:2)  = dfp(1:2)
dfp(3:5)  = dfp(3:5)/2.0d0
dfp(6:7)  = dfp(6:7)
end subroutine get_fpANDdfp

