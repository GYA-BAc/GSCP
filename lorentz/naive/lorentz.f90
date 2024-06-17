program lorentz
  implicit none

  real(8), parameter :: &
    q = 1.d0, & !C 
    m = 1.d0    !kg

  real(8), dimension(3), parameter :: &
    E = (/0.d0, 0.01d0, 0.1d0/), & !N/C
    B = (/0.d0, 0.d0, 0.1d0/)      !T    

  real(8), dimension(3) ::      &
    c_r = (/0.d0, 0.d0, 0.d0/), &
    c_v = (/1.d0, 0.d0, 0.d0/)
 
  real(8), dimension(4, 3) :: k ! for RK4

  real(8), parameter :: domain = 500
  real(8), parameter :: step = 0.1d0
  
  real(8) :: c_t = 0

  open(unit=100, file="r.dat")
  
  do while (c_t < domain)
    
    k(1,:) = c_v
    k(2,:) = c_v + (step/2.d0) * get_a(c_v+k(1,:)*step/2.d0)
    k(3,:) = c_v + (step/2.d0) * get_a(c_v+k(2,:)*step/2.d0)
    k(4,:) = c_v + step * get_a(c_v+k(3,:)*step)

    c_r = c_r + (step/6.d0) * (k(1,:)+2.d0*k(2,:)+2.d0*k(3,:)+k(4,:))
    
    c_v = c_v + get_a(c_v)*step

    write(100,*) c_r(1), c_r(2), c_r(3)

    c_t = c_t + step
  end do

  close(100)

contains

function cross_product(h, j) result(k)
  real(8), dimension(3), intent(in) :: h, j
  real(8), dimension(3)             :: k

  k(1) = h(2) * j(3) - h(3) * j(2)
  k(2) = h(3) * j(1) - h(1) * j(3)
  k(3) = h(1) * j(2) - h(2) * j(1)

end function

function get_a(v) result(a)
  real(8), dimension(3), intent(in) :: v
  real(8), dimension(3)             :: a
  
  a = (q*cross_product(v, B) + E*q)/m 
end function


end program
