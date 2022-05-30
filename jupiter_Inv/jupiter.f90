program jupiter_ch4

  use declaration, only : nmol, nlevel
  use read_file
  use info_content
  
  implicit none

  character (len=*), parameter :: file_fuel='jupiter_ch4.fuel'
  character (len=*), parameter :: file_pta='jupiter_ch4.pta'
  character (len=*), parameter :: file_spe='jupiter_ch4.spe'
  character (len=*), parameter :: file_res='jupiter_ch4.res'
  character (len=*), parameter :: file_inv='jupiter_ch4.inv'
  character (len=*), parameter :: file_sum='jupiter_ch4.sum'
  character (len=*), parameter :: file_avg='jupiter_ch4.avg'
  character (len=100), dimension(nmol) :: file_opa
  integer :: i, p_mol,j, n_inv, l, ii, shift, numinv, niter, cloudstate
  integer, dimension(nmol) :: inv
  integer, dimension(nmol) :: arr_inv
  integer, dimension(:), allocatable :: inv_pos
  real :: mui, mue, lat, ochi2, chi2
  real, dimension(:,:,:), allocatable :: kk
  real, dimension(nmol) :: vmr
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t, taucloud, grav
  real, dimension(:), allocatable :: wave, spec_obs, spec_syn, error, fwhm, dr
  real, dimension(:,:), allocatable :: delta, sigma_delta
  real, dimension(:,:,:), allocatable :: A
!==============================================================================
  call read_fuel(file_fuel,lat,mui,mue,p_mol,vmr,inv,arr_inv,cloudstate,file_opa)
  call read_pta(file_pta,p_mol,vmr,p,T,profil,taucloud)

  call read_spe(file_spe,wave,spec_obs,error)
  call Interpolate(wave, fwhm)
  allocate(spec_syn(size(wave)))

  n_inv= count(inv > 0)

  allocate(delta(nlevel,n_inv))
  allocate(sigma_delta(nlevel,n_inv))
  allocate(inv_pos(n_inv))
  allocate(dr(n_inv))

  inv_pos = 0
  ii = 1
  do l=1, nmol
       if (inv(l)>1) then
           inv_pos(ii) = l
           ii = ii + 1
       end if
  end do    

  if (n_inv == 0) then
    call info_content_ch4(p,t,profil,p_mol,shift,taucloud,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,n_inv,inv)
  else   
    allocate(kk(size(wave),nlevel,n_inv))
    allocate(A(nlevel,nlevel,n_inv))

    i = 0
!    ochi2 = 1e38
    ochi2 = sum(((spec_obs-spec_syn)/error)**2)
    do
      call info_content_ch4(p,t,profil,p_mol,shift,taucloud,lat,mue,fwhm,file_opa,size(wave),wave,spec_syn,n_inv,inv,kk)
      chi2 = sum(((spec_syn-spec_obs)/error)**2)
      print*, 'inverse --> ', i, ochi2, chi2
      print*,'chi2 --> ', abs(chi2-ochi2)/chi2
      if (abs(chi2-ochi2)/chi2 < 0.01) exit
      ochi2 = chi2
      call coeur_ch4(delta,sigma_delta,A,dr,p,size(wave),n_inv,kk,spec_obs,spec_syn,error)
      
      if (shift==0) then
      	do ii = 1, n_inv
      	     profil(:,inv_pos(ii)) = exp(alog(profil(:,(inv_pos(ii))))+delta(:,ii))
      	enddo
      elseif (shift ==1) then
        t = t + delta(:,1)
        do ii = 1, n_inv-1
             profil(:,inv_pos(ii)) = exp(alog(profil(:,(inv_pos(ii))))+delta(:,ii+shift))
        enddo
      end if

    i = i + 1
    niter = i
    enddo  
  end if  

!======= Write the summary of the inversion: number iterations, chi2 and degrees of freedom
  if (n_inv /= 0) then
    open(129,file=file_sum,status='unknown',form='formatted')
    write(129,'(I4)')(niter)
    write(129,'(F12.3)')(chi2)
    write(129,'(10F6.3)')(dr(:))
    close(129)
  end if

  call write_res(file_res,size(wave),wave,spec_obs,spec_syn,error)
  call write_inv(file_inv,p,t,profil,sigma_delta, n_inv,inv_pos, vmr)
  if (n_inv /= 0) call write_kernel(file_avg,n_inv,A)

end program jupiter_ch4
