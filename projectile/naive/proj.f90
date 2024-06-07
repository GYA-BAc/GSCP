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
  open(unit=200, file="z_t.dat")
  open(unit=300, file="v_x.dat")
  open(unit=400, file="v_z.dat")
  open(unit=500, file="a_x.dat")
  open(unit=600, file="a_z.dat")

  do while (c_r(3) >= 0)
    c_a = a_g + get_drag(p_v)/b_mass
    c_v = p_v+step*c_a
    c_r = c_r + step*(c_v+p_v)/2.d0
    p_v = c_v
    
    write(100,*) c_t, c_r(1) !x
    write(200,*) c_t, c_r(3) !z
    write(300,*) c_t, c_v(1)
    write(400,*) c_t, c_v(3)
    write(500,*) c_t, c_a(1)
    write(600,*) c_t, c_a(3)
    
    c_t = c_t + step
  end do
 
  close(100) 
  close(200) 
  close(300) 
  close(400) 
  close(500) 
  close(600) 

contains

function get_drag(v) result(f_d)
  real(8), dimension(3), intent(in) :: v
  real(8), dimension(3)             :: f_d

  f_d = -(0.5)*density*drag*(pi*(b_rad)**(2))*(v)**(2)
end function

end program
