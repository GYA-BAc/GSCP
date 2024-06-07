program proj

  
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)

  real(8), parameter :: drag = 0.1d0
  real(8), parameter :: density = 1.225d0   !kg/m^2

  real(8), parameter :: b_mass   = 0.140d0   !kg
  real(8), parameter :: b_rad = 0.0365d0 !meters


  real(8), parameter :: r_0(3) = [0.d0, 0.d0, 1.d0]
  real(8), parameter :: v_0(3) = [20.d0*cosd(50), 0.d0, 20.d0*sind(50)] 

  real(8), parameter :: c = density*drag*(pi*(b_rad)**(2))

end program
