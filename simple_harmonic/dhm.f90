program dhm
  ! This program demonstrates dampened, driven harmonic motion using the Verlet Algorithm

  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)

  real(8), parameter :: k = 1.d0
  real(8), parameter :: m = 1.d0
 
  real(8) :: b
  
  real(8) :: x_0, F_0, w_d
  real(8), allocatable, dimension(:) :: t, x, v, a, ke, spe

  real(8), parameter :: step = 0.01d0
  integer, parameter :: cycles = 20
  integer :: i, i_max

  print *, "How far to stretch the spring?"
  read(*,*) x_0

  print *, "Choose a dampening constant:"
  read(*,*) b

  print *, "Choose a driving force:"
  read(*,*) F_0

  print *, "Choose a driving frequency:"
  read(*,*) w_d

  
  i_max = floor((cycles-1.d0)*2*pi*sqrt(m/k)/step)
  allocate(t(i_max), x(i_max), v(i_max), a(i_max), ke(i_max), spe(i_max))

  t(1) = 0
  x(1) = x_0

  a(1) = get_a(x(1)) + (F_0)/m
  v(2) = a(1)*step
  x(2) = x(1) + v(2)*step + (0.5)*a(1)*(step)**(2)
  
  t(2) = step
  
  open(unit=100, file="x_t.dat")
  open(unit=200, file="v_t.dat")
  open(unit=300, file="a_t.dat")
  open(unit=400, file="k_t.dat")
  open(unit=500, file="u_t.dat")
  open(unit=600, file="e_t.dat")

  do i = 2, i_max - 1
    a(i) = get_a(x(i)) - b*v(i)/m + (F_0)*cos(w_d*t(i))/m
    x(i+1) = 2.d0*x(i) - x(i-1) + a(i)*(step)**(2)
    v(i+1) = (x(i+1) - x(i-1))/(2.d0*step)

    ke(i+1) = (0.5d0)*m*(v(i+1))**(2) 
    spe(i+1) = (0.5d0)*k*(x(i+1))**(2)

    t(i+1) = t(i) + step

    write(100,*) t(i+1), x(i+1)
    write(200,*) t(i+1), v(i+1)
    write(300,*) t(i+1), a(i)
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
