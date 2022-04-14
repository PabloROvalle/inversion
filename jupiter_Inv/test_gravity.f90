program test_gravity

  use gravity

  implicit none

  real :: grav

  call newgrav(0.,0.,grav)
  write(*,*) grav

end program test_gravity
