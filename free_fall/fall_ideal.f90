
program fall
  ! This program simulates basic free fall without air drag

  implicit none

  real(8) :: step = 0.1
  real(8) :: c_t = 0.0 !seconds
  real(8) :: c_h
  real(8) :: c_v
  real(8) :: p_v = -3.8d0 
  real(8) :: c_a = -9.8d0

  print *, "Input initial height (meters)"
  read(*,*) c_h
  
  open(unit=100, file="ideal_h_t.dat")
  open(unit=200, file="ideal_v_t.dat")
  open(unit=300, file="ideal_a_t.dat")

  do while (c_h >= 0.0)
    c_v = p_v + step*c_a
    c_h = c_h + step*(p_v+c_v)/2
     
    write(100,*) c_t, c_h
    write(200,*) c_t, c_v
    write(300,*) c_t, c_a
    
    p_v = c_v
    c_t = c_t + step
  end do
 
  close(100)
  close(200)
  close(300)

  print *, c_t, " seconds"
  print *, c_v, " m/s"

end program

