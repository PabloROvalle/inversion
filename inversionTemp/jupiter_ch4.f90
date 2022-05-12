program jupiter_ch4

  use declaration, only : nmol, nlevel
  use read_file
  use info_content
  
  implicit none

  character (len=7) :: gas
  character (len=50), parameter :: file_fuel='jupiter_ch4.fuel'
  character (len=50), parameter :: file_pta='jupiter_ch4.pta'
  character (len=50), parameter :: file_spe='jupiter_ch4.spe'
  character (len=50), parameter :: file_res='jupiter_ch4.res'
  character (len=68), dimension(nmol) :: file_opa
  integer :: p_mol,j
  real :: mui, mue, lat, dens
  real, dimension(:,:), allocatable :: kk
  real, dimension(nmol) :: vmr
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t, tau_haze, z, grav
  real, dimension(:), allocatable :: wave, spec_obs, spec_syn, error, fwhm
  integer :: i, a
  real, allocatable :: f(:), k(:)
  integer :: filesize, inputIndex

  !real, dimension(:), allocatable :: fwhm

!==============================================================================
  call read_fuel(file_fuel,lat,mui,mue,p_mol,vmr,file_opa)
  call read_pta(file_pta,p_mol,vmr,p,T,profil)

  call read_spe(file_spe,wave,spec_obs,error)
  allocate(spec_syn(size(wave)))
  
  call Interpolate(wave, fwhm)
  
  call info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn)

  call write_res(file_res,wave,spec_obs,spec_syn,error)
  
end program jupiter_ch4
