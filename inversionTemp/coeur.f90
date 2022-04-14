subroutine coeur(dt,p,kk,error,spectre,flux,nflux)

  use declaration
  implicit none

  integer :: i, j, k
  integer, intent(in) :: nflux

  Real, parameter :: c=0.75, factor=30e-1
  real :: alpha, beta
  real :: trace_temp, trace_haze, trace_error
  real, dimension(nflux), intent(in) :: error, flux, spectre
  real, dimension(nlevel), intent(inout) :: dt
  real, dimension(nlevel), intent(in) :: p
  real, dimension(nflux,nlevel), intent(in) :: kk
  real, dimension(nlevel,nlevel) :: s
  real, dimension(nlevel,nflux) ::  w, aux3
  real, dimension(1,nflux) ::  v
  real, dimension(nflux,nlevel) :: aux
  real, dimension(nflux,nflux) :: aux1, aux2
!------------------------------------------------------------------------------
  do i=1,nlevel
     do j=1,nlevel
        s(i,j) = exp( -alog(p(i)/p(j))**2 / (2*c**2))
     enddo
  enddo

  trace_error = sum((/(error(i)**2,i=1,nflux)/) / sqrt(real(nflux,kind=4)))
      
  aux = matmul(kk,s) !--- Temperature
  aux1 = matmul(aux,transpose(kk))
  trace_temp = sum((/(aux1(i,i),i=1,nflux)/))
  alpha = factor * (trace_error/trace_temp)

  aux1 = alpha*aux1
  do i=1,nflux
     aux1(i,i) = aux1(i,i) + error(i)**2
  enddo

  call matinv(aux1,nflux,nflux,aux2)

  aux3 = matmul(transpose(kk),aux2)
  w = matmul(s,aux3)
  dt = alpha * matmul(w,(flux-spectre))

  return
end subroutine coeur