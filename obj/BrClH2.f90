module pipvariables
implicit none

! Parameters required for interfacing with
! the "pes" and "training" software...
! must match the subroutines below

integer,parameter :: Natoms = 4  ! Number of atoms
integer,parameter :: Np     = 7  ! Number of polynomials

end module pipvariables

! The PIPs/FIs are defined here... switch
! the content of this subroutine with
! whatever PIPs/FIs you want to use
subroutine get_p(r,p)
implicit none
real(kind=8),intent(in) :: r(6)
real(kind=8),intent(out) :: p(7)
p(1)=r(5)+r(4)
p(2)=r(3)+r(2)
p(3)=r(5)**2+r(4)**2
p(4)=r(3)**2+r(2)**2
p(5)=r(5)*r(3)+r(4)*r(2)
p(6)=r(6)
p(7)=r(1)
end subroutine get_p


! For consistent units, for each Nth degree polynomial,
! apply the N-th root (this is less important but
! still must be manually done each time)
subroutine get_fp(Np,p,fp)
implicit none
integer,intent(in) :: Np
double precision,dimension(Np),intent(in) :: p
double precision,dimension(Np),intent(out) :: fp

fp(1:2)     = p(1:2)
fp(3:5)     = p(3:5)**(0.5d0)
fp(6:7)     = p(6:7)
end subroutine get_fp

