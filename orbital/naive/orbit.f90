program orbit
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)
  real(8), parameter :: GM = 4.d0*(pi)**(2) ! AU^3/yr^2
  real(8), parameter :: mass = 1.d0

  real(8), dimension(3) ::         &
    c_r = (/1.d0, 0.d0, 0.d0/),    &
    c_v = (/0.d0, 2.d0*pi, 0.d0/)

  real(8), dimension(4, 3) :: k ! for RK4

  real(8) :: ke, gpe

  real(8), parameter :: step = 1.d-4
  real(8) :: c_t = 0

  real(8) :: domain = 1 !yr

  open(unit=100, file="r_t.dat")
  open(unit=110, file="k_t.dat")
  open(unit=120, file="u_t.dat")
  open(unit=130, file="e_t.dat")

  do while (c_t < domain)
    
    k(1,:) = c_v
    k(2,:) = c_v + (step/2.d0) * get_a(c_r+k(1,:)*step/2.d0)
    k(3,:) = c_v + (step/2.d0) * get_a(c_r+k(2,:)*step/2.d0)
    k(4,:) = c_v + step * get_a(c_r+k(3,:)*step)

    c_r = c_r + (step/6.d0) * (k(1,:)+2.d0*k(2,:)+2.d0*k(3,:)+k(4,:))
    
    c_v = c_v + get_a(c_r)*step

    ke = (0.5d0) * mass * norm2(c_v)**2
    gpe = mass * GM * norm2(c_r)

    write(100,*) c_r(1), c_r(2), c_r(3)  
    write(110,*) c_t, ke
    write(120,*) c_t, gpe
    write(130,*) c_t, ke+gpe


    c_t = c_t + step
  end do

  close(100)
  close(110)
  close(120)
  close(130)

contains

function get_a(r) result(a)
  real(8), dimension(3), intent(in) :: r
  real(8), dimension(3)             :: a
  
  a = -GM*r/(norm2(r)**3.d0)

end function

end program
