program orbit
  implicit none

  real(8), parameter :: pi = 4.d0*atan(1.d0)
  real(8), parameter :: GM = 4.d0*(pi)**(2) ! AU^3/yr^2

  real(8), dimension(3) ::         &
    c_r = (/1.d0, 0.d0, 0.d0/),    &
    c_v = (/0.d0, 2.d0*pi, 0.d0/), &
    c_a = (/0.d0, 0.d0, 0.d0/)
  
  real(8) :: k_1, k_2, k_3, k_4 ! for RK4

  real(8), parameter :: step = 0.01d0
  real(8) :: c_t = 0

  real(8) :: domain = 100

  do while (c_t < domain)
    
    k_1 = c_v
    k_2 = c_v + step*get_a(c_r+k_1*step/2.d0)
    k_3 = c_v + step*get_a(c_r+k_2*step/2.d0)
    k_4 = c_v + step*get_a(c_r+k_3*step)

    c_r = c_r + step/6.d0 * (k_1+k_2+k_3+k_4)

    c_v = c_v + get_a(c_x)*step

    c_t = c_t + step
  end do

contains

function get_a(r) result(a)
  real(8), dimension(3), intent(in) :: r
  real(8), dimension(3)             :: a
  
  a = GM/(r**2)

end function

end program
