program lorentz
  implicit none

  real(8), parameter :: &
    q = 1.d0, & !C 
    m = 1.d0    !kg

  real(8), dimension(3), parameter :: &
    E = [0.d0, 0.01d0, 0.1d0], & !N/C
    B = [0.d0, 0.d0, 0.1d0]      !T    

  ! x y z vx vy vz
  real(8), dimension(6) :: state = [0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0]
  
  real(8), dimension(4, 6) :: k ! for RK4
  real(8), dimension(6) :: temp

  real(8), parameter :: domain = 1000
  real(8), parameter :: step = 0.1d0
  
  real(8) :: c_t = 0

  open(unit=100, file="r.dat")
  
  do while (c_t < domain)
    
    k(1,1:3) = state(4:6)
    k(1,4:6) = get_a(state(4:6))

    temp = state + (0.5)*step*k(1,:)
    
    k(2,1:3) = temp(4:6)
    k(2,4:6) = get_a(temp(4:6))

    temp = state + (0.5)*step*k(2,:)
    
    k(3,1:3) = temp(4:6)
    k(3,4:6) = get_a(temp(4:6))

    temp = state + step*k(3,:)
    
    k(4,1:3) = temp(4:6)
    k(4,4:6) = get_a(temp(4:6))

    state = state + (step/6.d0) * (k(1,:)+2.d0*k(2,:)+2.d0*k(3,:)+k(4,:))

    write(100,*) state(1), state(2), state(3)

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
