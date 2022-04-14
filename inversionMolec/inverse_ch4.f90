program inverse_ch4

  use declaration, only : nmol, nlevel
  use read_file
  use info_content
  
  implicit none

  character (len=7) :: gas
  character (len=50), parameter :: file_fuel='jupiter_ch4.fuel'
  character (len=50), parameter :: file_pta='jupiter_ch4.pta'
  character (len=50), parameter :: file_spe='jupiter_ch4.spe'
  character (len=50), parameter :: file_res='jupiter_ch4.res'
  character (len=50), parameter :: file_inv='jupiter_ch4.inv'
  character (len=68), dimension(nmol) :: file_opa
  integer :: p_mol, i
  real :: mui, mue,  lat, dens, chi2, ochi2, dr
  real, dimension(:,:), allocatable :: kk
  real, dimension(nmol) :: vmr
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t, tau_haze, z, grav, dvmr, sigma, tau_aer
  real, dimension(:), allocatable :: wave, spec_obs, spec_syn, error, fwhm
  real, dimension(nlevel,nlevel) :: A

!==============================================================================
  call read_fuel(file_fuel,lat,mui,mue,p_mol,vmr,file_opa)
  call read_pta(file_pta,p_mol,vmr,p,T,profil)
  
  call read_spe(file_spe,wave,spec_obs,error)
  call Interpolate(wave, fwhm)
  allocate(spec_syn(size(wave)))
  allocate(kk(size(wave),nlevel))
  call info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,kk,tau_aer)

!  ----- Initialisation for inversion
  i = 0
  ochi2 = sum(((spec_obs-spec_syn)/error)**2)
  print*, i, ochi2
  !  ------ Iteration for inversion
  do
     i = i + 1
     call coeur_ch4(dvmr,sigma,p,size(wave),kk,spec_obs,spec_syn, error,A, dr)
     profil(:,5) = exp(alog(profil(:,5))+ dvmr)
     call info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,kk,tau_aer)
     chi2 = sum(((spec_obs-spec_syn)/error)**2)
     print*,'Actual difference (%) -->', (abs(chi2-ochi2)/chi2),';  Degrees of freedom -->', dr
     write(*,*)i, chi2
     if (abs(chi2-ochi2)/chi2 < 0.05) exit
     ochi2 = chi2
  enddo
     
  call write_inv(file_inv,p,profil,sigma)
  call write_res(file_res,size(wave),wave,spec_obs,spec_syn,error)
  
end program inverse_ch4
