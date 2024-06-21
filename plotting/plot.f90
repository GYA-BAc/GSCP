program plot
  ! This program uses the fourier series to approximate a sawtooth function
  implicit none

  real(8), parameter :: step = 1d-2
  real(8), parameter :: pi = 4.d0*atan(1.d0)
  
  integer :: n_max
  real(8) :: x

  print *, "Enter maximum number of terms (n_max): "
  read(*,*) n_max
  
  ! generate datasets
  x = 0.d0
  open(unit=100, file="exact.dat")
  open(unit=200, file="approx.dat")
  
  do while (x < 2.d0*pi)
    write(100,*) x, exact(x)
    write(200,*) x, approx(x, n_max)
    x = x + step
  end do
  
  close(100)
  close(200)

contains

function exact(x) result(y)
  real(8), intent(in) :: x
  real(8)             :: y
  
  
  if (x<0) then
    y = 0
  else if (x<pi) then
    y = x
  else if (x<2.d0*pi) then
    y = x - 2.d0*pi
  else
    y = 0
  end if

end function

function approx(x, n_max) result(y)
  real(8), intent(in) :: x
  integer, intent(in) :: n_max
  real(8)             :: y
  
  real(8) :: term_b
  integer :: n

  ! solve for term b
  term_b = 0.d0
  do n = 1, n_max, 1
    term_b = term_b + b(real(n, 8))*sin(real(n, 8)*x)
  end do

  y = term_b

end function

function b(n) result(b_n)
  real(8), intent(in) :: n
  real(8)             :: b_n 
    
  b_n = (2.d0 * (-1)**(n+1))/n
end function

end program
