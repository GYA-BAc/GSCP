program shm
  real(8), parameter :: k = 0.5d0
  real(8), parameter :: m = 5.d0

  real(8), parameter :: step = 0.05d0
  real(8), parameter :: domain = 20
 
  real(8) :: c_t = step
  real(8) :: n_x
  real(8) :: c_x
  real(8) :: p_x
  real(8) :: c_v
  real(8) :: c_a

  real(8) :: ke
  real(8) :: spe

  print *, "How far to stretch the spring?"
  read(*,*) c_x

  c_a = get_a(c_x)
  c_v = c_a*step
  n_x = c_x + c_v + (0.5)*c_a*(step)**(2)

  open(unit=100, file="x_t.dat")
  open(unit=200, file="v_t.dat")
  open(unit=300, file="a_t.dat")
  open(unit=400, file="k_t.dat")
  open(unit=500, file="u_t.dat")
  open(unit=600, file="e_t.dat")

  do while (c_t < domain)
    p_x = n_x

    c_a = get_a(n_x)
    n_x = 2.d0*n_x - c_x + c_a*(step)**(2)
    
    c_v = (n_x - c_x)/(2.d0*step)
    c_x = p_x

    ke = (0.5d0)*m*(c_v)**(2) 
    spe = (0.5d0)*k*(c_x)**(2)

    write(100,*) c_t, n_x
    write(200,*) c_t, c_v
    write(300,*) c_t, get_a(n_x)
    write(400,*) c_t, ke
    write(500,*) c_t, spe
    write(600,*) c_t, ke+spe

    c_t = c_t + step
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
