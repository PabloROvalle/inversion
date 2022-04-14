program jupiter_ch4

  use declaration, only : nmol, nlevel
  use read_file
  use info_content
  
  implicit none

  character (len=*), parameter :: file_fuel='jupiter_ch4.fuel'
  character (len=*), parameter :: file_pta='jupiter_ch4.pta'
  character (len=*), parameter :: file_spe='jupiter_ch4.spe'
  character (len=*), parameter :: file_res='jupiter_ch4.res'
  character (len=100), dimension(nmol) :: file_opa
  integer :: i, p_mol,j, n_inv, nummol, m, imol_inv, l, xm,xt
  integer :: n_invx
  integer, dimension(nmol) :: inv, aux
  real :: mui, mue, lat, ochi2, chi2
  real, dimension(:,:,:), allocatable :: kk
  real, dimension(nmol) :: vmr
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t,sigma, tau_haze, z, grav, dt
  real, dimension(:,:), allocatable :: dvmr, sigma_vmr
  real, dimension(:), allocatable :: wave, spec_obs, spec_syn, error, fwhm
  real, dimension(:,:,:), allocatable :: A
!==============================================================================
  call read_fuel(file_fuel,lat,mui,mue,p_mol,vmr,inv,file_opa, nummol)
  call read_pta(file_pta,p_mol,vmr,p,T,profil)

  
  call read_spe(file_spe,wave,spec_obs,error)
  call Interpolate(wave, fwhm)
  allocate(spec_syn(size(wave)))

!  print*,'nummol ', nummol
   n_inv= count(inv > 0)

!  print*,'n_inv ', n_inv
    
  if (n_inv == 0) then
    allocate(dvmr(nlevel,count(inv>1)))
    allocate(sigma_vmr(nlevel,count(inv>1)))
    allocate(A(nlevel,nlevel,count(inv>0)))
    call info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,n_inv,inv,imol_inv, xm,xt)
  else   
    allocate(kk(size(wave),nlevel,n_inv))
    allocate(dvmr(nlevel,count(inv>1)))
    allocate(sigma_vmr(nlevel,count(inv>1)))
    allocate(A(nlevel,nlevel,count(inv>0)))

    i = 0
    ochi2 = 1e38
    do
      call info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,n_inv,inv,imol_inv,xm,xt,kk)
      chi2 = sum(((spec_syn-spec_obs)/error)**2)
      print*, 'inverse --> ', i, ochi2, chi2
      print*,'chi2 --> ', abs(chi2-ochi2)/chi2
      if (abs(chi2-ochi2)/chi2 < 0.01) exit
      ochi2 = chi2
      call coeur_ch4(dvmr,dt,sigma,sigma_vmr,A,p,size(wave),p_mol,inv,xm,xt,kk,spec_obs,spec_syn,error)
      m = 1
      do l=1, p_mol
            if (inv(l)==1) then
	      	t = t + dt
            elseif (inv(l)==2) then
	      	profil(:,l) = exp(alog(profil(:,l))+dvmr(:,m))
                m = m+1
            endif
      enddo

      i = i + 1
    enddo  
  end if  

  call write_res('jupiter_ch4.res',size(wave),wave,spec_obs,spec_syn,error)
  call write_inv('jupiter_ch4.inv',p,profil,sigma_vmr, nummol, vmr, t, sigma, inv)
  if (n_inv /= 0) call write_kernel('jupiter_ch4.avg',A)

end program jupiter_ch4
