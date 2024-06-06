
program fall
  implicit none

  real(8) :: step = 0.1
  real(8) :: c_t = 0.0 !seconds
  real(8) :: c_h
  real(8) :: c_v = -3.8d0
  real(8) :: c_a = -9.8d0

  print *, "Input initial height (meters)"
  read(*,*) c_h

  do while (c_h >= 0.0)
    c_v = c_v + c_a*step
    c_h = c_h + c_v*step
     
    c_t = c_t + step
    
  end do
  
  print *, c_t, " seconds"
  print *, c_v, " m/s"

end program

