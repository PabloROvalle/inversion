module convolution

  implicit none
  
contains

  subroutine convol_jwst(nrad,radiance,fwhm,fnu,nflux,wave,synthetic)

    implicit none
    
    integer :: i, j, k
    integer, intent(in) ::  nflux, nrad
    
    real, dimension(nflux), intent(in) :: fwhm
    real, dimension(nrad), intent(in) :: radiance, fnu
    real, dimension(nrad) :: aux
    real, dimension(nflux), intent(in) :: wave
    real, dimension(nflux), intent(out) :: synthetic
!------------------------------------------------------------------------------
    do i=1, nflux
       where (abs(fnu - wave(i)) < 5*fwhm(i))
          aux = 2**(-((fnu-wave(i))/(fwhm(i)/(2*1.)))**2)
       elsewhere
          aux = 0.
       end where
       synthetic(i) = sum(radiance*aux) / sum(aux)
    enddo
    
  end subroutine convol_jwst

end module convolution
