
program fall
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)

  real(8), parameter :: drag = 0.1d0
  real(8), parameter :: density = 1.225d0   !kg/m^2
  
  real(8), parameter :: b_mass   = 0.140d0   !kg
  real(8), parameter :: b_rad = 0.0365d0 !meters

  real(8), parameter :: f_g = 9.8d0 * b_mass

  real(8), parameter :: v_0 = -3.8d0 !m
  
  real(8), parameter :: c = density*drag*(pi*(b_rad)**(2))
  real(8), parameter :: v_t = sqrt(2.d0*f_g/c)

  real(8) :: step = 0.1d0
  real(8) :: c_t = 0.d0 !seconds
  real(8) :: c_h
  real(8) :: c_v
  real(8) :: p_v = v_0
  real(8) :: c_a

  print *, "Terminal Velocity:", v_t
  print *, "Input initial height (meters)"
  read(*,*) c_h

  open(unit=100, file="h_t.dat")
  open(unit=200, file="v_t.dat")
  open(unit=300, file="a_t.dat")

  do while (c_h >= 0.0)

    c_a = -9.8d0 + get_drag(p_v)/b_mass 
    c_v = p_v+step*c_a
    c_h = c_h + step*(c_v+p_v)/2.d0
    p_v = c_v

    c_t = c_t + step
    
    write(100,*) c_t, c_h 
    write(200,*) c_t, c_v 
    write(300,*) c_t, c_a 
  
  end do
  
  print *, c_t, " seconds"

  close(100)
  close(200)
  close(300)

contains

function get_drag(v) result(f_d)
  real(8), intent(in) :: v
  real(8)             :: f_d

  f_d = (0.5)*c*(v)**(2)

end function

end program

