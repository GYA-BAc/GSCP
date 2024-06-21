program proj
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)

  real(8), parameter :: drag = 0.1d0
  real(8), parameter :: density = 1.225d0   !kg/m^2

  real(8), parameter :: b_mass = 0.06796d0   !kg
  real(8), parameter :: b_rad = 0.02532d0/2.d0 !meters

  real(8), parameter :: a_g(3) = [0.d0, 0.d0, -9.8d0]

  real(8), dimension(3) :: c_r, p_v, c_v, c_a

  real(8) :: h_max = 0

  real(8) :: c_t = 0.d0
  real(8) :: step = 0.01d0

  real(8) :: v_0
  real(8) :: theta

  print *, "Enter v_0: "
  read(*,*) v_0
  print *, "Enter theta: "
  read(*,*) theta
  
  c_r(:) = [0.d0, 0.d0, 1.075d0]
  c_v(:) = [0.d0, 0.d0, 0.d0]
  p_v(:) = [v_0*cosd(theta), 0.d0, v_0*sind(theta)] 
  c_a(:) = [0.d0, 0.d0, 0.d0]


  open(unit=100, file="r.dat")
  open(unit=110, file="x_t.dat")
  open(unit=120, file="z_t.dat")
  open(unit=130, file="vx_t.dat")
  open(unit=140, file="vz_t.dat")
  open(unit=150, file="ax_t.dat")
  open(unit=160, file="az_t.dat")

  do while (c_r(3) >= 0)
    c_a = a_g! + get_drag(p_v)/b_mass
    c_v = p_v+step*c_a
    c_r = c_r + step*(c_v+p_v)/2.d0
    p_v = c_v
    
    if (c_r(3) > h_max) then
      h_max = c_r(3)
    end if

    write(100,*) c_r(1), c_r(2), c_r(3)
    write(110,*) c_t, c_r(1) !x
    write(120,*) c_t, c_r(3) !z
    write(130,*) c_t, c_v(1)
    write(140,*) c_t, c_v(3)
    write(150,*) c_t, c_a(1)
    write(160,*) c_t, c_a(3)
    
    c_t = c_t + step
  end do

  print *, "Duration of flight: ", c_t, " s"
  print *, "X range was:        ", c_r(1), " m"
  print *, "Max height was:     ", h_max, " m" 

  close(100) 
  close(110) 
  close(120) 
  close(130) 
  close(140) 
  close(150) 
  close(160)

contains

function get_drag(v) result(f_d)
  real(8), dimension(3), intent(in) :: v
  real(8), dimension(3)             :: f_d

  f_d = -(0.5)*density*drag*(pi*(b_rad)**(2))*v*norm2(v)
end function

end program
