program rlc
  implicit none

  integer, parameter :: max_n = 12000
  integer :: i

  real(8) :: q_0, C
  real(8) :: c_q(max_n)

  real(8) :: L
  real(8) :: d_i(max_n) ! di / dt

  real(8) :: i_0, R
  real(8) :: c_i(max_n) ! current current in circuit

  real(8) :: E_0 ! driving voltage
  real(8) :: w_0 ! resonant frequency

  real(8), parameter :: step = 0.01d0
  real(8), parameter :: domain = 25.d0


  real(8) :: c_t = 0.d0
  
  print *, "Enter capacitor q_0:"
  read(*,*) q_0
  print *, "Enter capacitor C:"
  read(*,*) C

  print *, "Enter inductor L:"
  read(*,*) L
  
  print *, "Enter i_0:"
  read(*,*) i_0
  print *, "Enter resistor R:"
  read(*,*) R

  print *, "Enter driving voltage:"
  read(*,*) E_0

  print *, "Enter w_0"
  read(*,*) w_0
  
  open(unit=100, file="q_t.dat")
  open(unit=120, file="i_t.dat")
  open(unit=130, file="Vr_t.dat")
  open(unit=140, file="Vrlc_t.dat")
  open(unit=150, file="Ec_t.dat")
  open(unit=160, file="El_t.dat")
  open(unit=170, file="E_t.dat")

  c_i(1) = i_0
  c_q(1) = q_0
  d_i(1) = get_di(c_t, c_q(1), c_i(1), L, E_0, w_0, R, C)
  c_t = c_t + step
  c_i(2) = d_i(1)*step
  c_q(2) = c_q(1) + step*c_i(1) + (0.5)*d_i(1)*(step)**(2)

  
  do i = 2, max_n - 1 
   
    d_i(i) = get_di(c_t, c_q(i), c_i(i), L, E_0, w_0, R, C)
    c_q(i+1) = 2.d0*c_q(i) - c_q(i-1) + d_i(i)*(step)**(2)
    c_i(i+1) = (c_q(i+1) - c_q(i-1))/(2.d0*step)

    c_t = c_t + step

    write(100,*) c_t, c_q(i+1)
    write(120,*) c_t, c_i(i+1)
    write(130,*) c_t, c_i(i+1)*R
    write(140,*) c_t, c_i(i+1)*R + d_i(i)*L + c_q(i+1)/C
    write(150,*) c_t, (0.5d0)*C*(c_q(i+1)/C)**(2) 
    write(160,*) c_t, (0.5d0)*L*(c_i(i+1))**(2) 
    write(170,*) c_t, (0.5d0)*C*(c_q(i+1)/C)**(2) + (0.5d0)*L*(c_i(i+1))**(2)  

    c_t = c_t + step
  end do

  close(100)
  close(120)
  close(130)
  close(140)
  close(150)
  close(160)
  close(170)

contains

function get_di(t, c_q, c_i, L, E_0, w_0, R, C) result(d_i)
  real(8), intent(in) :: t, c_q, c_i, L, E_0, w_0, R, C
  real(8)             :: d_i

  !d_i = (1/L) * (E_0*sin(w_0*t) - R*c_i - c_q/C)
  d_i = -(1/L) * (R*c_i + c_q/C)

end function

end program rlc
