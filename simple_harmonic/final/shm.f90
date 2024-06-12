program shm
  real(8), parameter :: k = 0.5d0
  real(8), parameter :: m = 5.d0

  real(8), parameter :: step = 0.001d0
  integer, parameter :: max_i = 20000
 
  real(8) :: x_0
  real(8) :: t(max_i), x(max_i), v(max_i), a(max_i), ke(max_i), spe(max_i)
  integer :: i

  print *, "How far to stretch the spring?"
  read(*,*) x_0

  t(1) = 0
  x(1) = x_0

  t(2) = step
  a(2) = get_a(x(1))
  v(2) = a(1)*step
  x(2) = x(1) + step*v(1) + (0.5)*a(2)*(step)**(2)
  
  open(unit=100, file="x_t.dat")
  open(unit=200, file="v_t.dat")
  open(unit=300, file="a_t.dat")
  open(unit=400, file="k_t.dat")
  open(unit=500, file="u_t.dat")
  open(unit=600, file="e_t.dat")

  do i = 2, max_i - 1
    a(i+1) = get_a(x(i))
    x(i+1) = 2.d0*x(i) - x(i-1) + a(i+1)*(step)**(2)
    v(i+1) = (x(i+1) - x(i-1))/(2.d0*step)

    ke(i+1) = (0.5d0)*m*(v(i+1))**(2) 
    spe(i+1) = (0.5d0)*k*(x(i+1))**(2)

    t(i+1) = t(i) + step

    write(100,*) t(i+1), x(i+1)
    write(200,*) t(i+1), v(i+1)
    write(300,*) t(i+1), a(i+1)
    write(400,*) t(i+1), ke(i+1)
    write(500,*) t(i+1), spe(i+1)
    write(600,*) t(i+1), ke(i+1)+spe(i+1)

  end do

  close(100)
  close(200)
  close(300)
  close(400)
  close(500)
  close(600)

contains

function get_a(x) result(a)
  real(8), intent(in) :: x
  real(8)             :: a

  a = -k*x/m
end function
end program
