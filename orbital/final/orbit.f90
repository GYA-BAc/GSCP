program orbit
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)
  real(8), parameter :: GM = 4.d0*(pi)**(2) ! AU^3/yr^2
  real(8), parameter :: mass = 1.d0

  ! x y z vx vy vz
  real(8), dimension(6) :: state = [1.d0, 0.d0, 0.d0, 0.d0, 2.d0*pi, 0.d0]
  
  real(8), dimension(4, 6) :: k ! for RK4
  real(8), dimension(6) :: temp

  real(8) :: ke, u

  real(8), parameter :: step = 1.d-2
  real(8) :: c_t = 0

  real(8) :: domain = 1 !yr

  open(unit=100, file="r_t.dat")
  open(unit=110, file="k_t.dat")
  open(unit=120, file="u_t.dat")
  open(unit=130, file="e_t.dat")

  do while (c_t < domain)
   
    k(1,1:3) = state(4:6)
    k(1,4:6) = get_a(state)

    temp = state + (0.5)*step*k(1,:)
    
    k(2,1:3) = temp(4:6)
    k(2,4:6) = get_a(temp(1:3))

    temp = state + (0.5)*step*k(2,:)
    
    k(3,1:3) = temp(4:6)
    k(3,4:6) = get_a(temp(1:3))

    temp = state + step*k(3,:)
    
    k(4,1:3) = temp(4:6)
    k(4,4:6) = get_a(temp(1:3))

    state = state + (step/6.d0) * (k(1,:)+2.d0*k(2,:)+2.d0*k(3,:)+k(4,:))
    
    ke = (0.5d0) * mass * norm2(state(4:6))**2
    u = - mass * GM / norm2(state(1:3))

    write(100,*) state(1), state(2), state(3)  
    write(110,*) c_t, ke
    write(120,*) c_t, u
    write(130,*) c_t, ke+u


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
  
  a = -GM*mass*r/(norm2(r)**3.d0)

end function

end program
