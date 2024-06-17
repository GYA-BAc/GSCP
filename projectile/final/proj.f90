program proj
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)

  real(8), parameter :: drag = 0.1d0
  real(8), parameter :: density = 1.225d0   !kg/m^2

  real(8), parameter :: b_mass   = 0.140d0   !kg
  real(8), parameter :: b_rad = 0.0365d0 !meters


  real(8), parameter :: r_0(3) = [0.d0, 0.d0, 1.d0]
  real(8), parameter :: v_0(3) = [20.d0*cosd(50.d0), 0.d0, 20.d0*sind(50.d0)] 

  real(8), parameter :: a_g(3) = [0.d0, 0.d0, -9.8d0]

  real(8) :: c_r(3) = r_0(:)
  real(8) :: c_v(3) = [0.d0, 0.d0, 0.d0]
  real(8) :: p_v(3) = v_0(:)
  real(8) :: c_a(3) = [0.d0, 0.d0, 0.d0]

  real(8) :: c_t = 0.d0
  real(8) :: step = 0.01d0
  
  open(unit=100, file="x_t.dat")
  open(unit=110, file="z_t.dat")
  open(unit=120, file="vx_t.dat")
  open(unit=130, file="vz_t.dat")
  open(unit=140, file="ax_t.dat")
  open(unit=150, file="az_t.dat")

  do while (c_r(3) >= 0)
    c_a = a_g + get_drag(p_v)/b_mass
    c_v = p_v+step*c_a
    c_r = c_r + step*(c_v+p_v)/2.d0
    p_v = c_v
    
    write(100,*) c_t, c_r(1) !x
    write(110,*) c_t, c_r(3) !z
    write(120,*) c_t, c_v(1)
    write(130,*) c_t, c_v(3)
    write(140,*) c_t, c_a(1)
    write(150,*) c_t, c_a(3)
    
    c_t = c_t + step
  end do
 
  close(100) 
  close(110) 
  close(120) 
  close(130) 
  close(140) 
  close(150) 

contains

function get_drag(v) result(f_d)
  real(8), dimension(3), intent(in) :: v
  real(8), dimension(3)             :: f_d

  f_d = -(0.5)*density*drag*(pi*(b_rad)**(2))*v*abs(v)
end function

end program
