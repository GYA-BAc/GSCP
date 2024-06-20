program heat
  implicit none
 
  integer, parameter :: N = 61
  real(8), parameter :: &
    L = 1.d0, &
    K = 1.d0

  real(8), parameter :: tau = 0.0001d0
  real(8), parameter :: h = L/N

  integer, parameter :: epochs = 5000
  integer :: i
  integer :: x

  integer :: cur
  integer :: prev = 1

  real(8), dimension(2, N) :: T = 0
  T(1, 31) = 1.d0/tau
  
  open(unit=100, file="u_t.dat")

  do i=1, epochs  
    cur = mod(i, 2)+1
    
    do x=2, N-1

      T(cur, x) = (K*tau/h) &
        * (T(prev, x+1) - 2.d0*T(prev, x) + T(prev, x-1)) &
        + T(prev, x)
      
    end do

    write(100,*) T(cur,:)

    prev = cur
  end do

  close(100)

end program
