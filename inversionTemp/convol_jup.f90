subroutine convol_jup(radiance,rayleigh,repmin,repmax,nflux,wavesup,waveinf&
     &,spectre)

  use declaration
  implicit none

  integer :: i, j, k
  integer, intent(in) ::  repmin, repmax, nflux

  real, dimension(:), intent(in) :: radiance
  real, dimension(size(radiance)) :: fnu, aux
  real, intent(in) :: rayleigh
  real, dimension(:), intent(in) :: wavesup, waveinf
  real, dimension(:), intent(inout) :: spectre
!------------------------------------------------------------------------------
  fnu = gnu00 + (/((j-1),j=(repmin-1)*nfreq+1, repmax*nfreq)/) * dgnu
  do i=1, nflux
     where (fnu >= waveinf(i) .and. fnu<= wavesup(i))
        aux = 1.
     elsewhere
        aux = 0.
     end where
     spectre(i) = sum(radiance*aux) / sum(aux)
  enddo

  return
end subroutine convol_jup
