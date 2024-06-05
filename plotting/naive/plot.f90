program plot
  implicit none

  real    :: step = 0.01
  real :: pi = acos(-1.0)
  
  integer :: n_max
  real :: x

  print *, "Enter maximum number of terms: "
  read(*,*) n_max
  
  ! generate datasets
  x = 0
  open(unit=1, file="exact.dat")
  open(unit=2, file="approx.dat")
  
  do while (x < 2*pi)
    write(1,*) x, exact(x)
    write(2,*) x, approx(x, n_max)
    x = x + step
  end do
  
  close(1)
  close(2)

contains

function exact(x) result(y)
  real, intent(in) :: x
  real             :: y
  
  
  if (x<0) then
    y = 0
  else if (x<pi) then
    y = x
  else if (x<2*pi) then
    y = x - 2*pi
  else
    y = 0
  end if

end function

function approx(x, n_max) result(y)
  real, intent(in)    :: x
  integer, intent(in) :: n_max
  real                :: y
  
  real :: term_b
  integer :: n

  ! solve for term b
  term_b = 0
  do n = 1, n_max, 1
    term_b = term_b + b(real(n))*sin(real(n)*x)
  end do

  y = term_b

end function

function b(n) result(b_n)
  real, intent(in) :: n
  real             :: b_n 
    
  b_n = (2 * (-1)**(n+1))/n
end function

end program
