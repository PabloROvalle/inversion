subroutine convol_jwst(nrad,radiance,fwhm,fnu,nflux,wave,synthetic)

  implicit none

  integer :: i, j, k
  integer, intent(in) ::  nflux, nrad

  real, intent(in) :: fwhm
  real, dimension(nrad), intent(in) :: radiance, fnu
  real, dimension(nrad) :: aux
  real, dimension(nflux), intent(in) :: wave
  real, dimension(nflux), intent(inout) :: synthetic
!------------------------------------------------------------------------------
  write(*,*) shape(nrad),shape(radiance),shape(fwhm),shape(fnu),shape(nflux),shape(wave),shape(synthetic)
  write(*,*) 'convol', nrad,nflux
  do i=1, nflux
  write(*,*) wave(i)
     where (abs(fnu - wave(i)) < 3*fwhm)
        aux = 2**(-((fnu-wave(i))/fwhm)**2)
     elsewhere
        aux = 0.
     end where
     write(*,*)sum(radiance*aux), sum(aux)
     synthetic(i) = sum(radiance*aux) / sum(aux)
     write(*,*)synthetic(i)
  enddo

  return
end subroutine convol_jwst
