module pseudorandoms
! module for pseudo random numbers
 implicit none
 private
 public                :: randint, r4_uniform, get_random_seed
 contains
!
 subroutine get_random_seed()
  implicit none
  integer, allocatable :: seed(:)
  integer :: n
  !
  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)
  write (*, *) seed
end subroutine get_random_seed
! 
 integer function randint(i,j)
  integer,intent(in)    :: i,j
  real                  :: r
  call random_number(r)
  randint = i + floor((j+1-i)*r)
  return
 end function randint
!
 real function r4_uniform(x,y)
  implicit none
  real,intent(in)       :: x,y
  real                  :: r
  call random_number(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module pseudorandoms
