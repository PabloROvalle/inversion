module read_file

  implicit none
  
contains

  subroutine read_fuel(file_in,lat,mui,mue,p_mol,vmr,inv,arr_inv,cloudstate,file_opa)
    
    use declaration, only : nmol, file_root
    implicit None

    character (len=7) :: gas, cloud
    character (len=*), intent(in) :: file_in
    character (len=*), dimension(nmol), intent(out) :: file_opa
    integer :: l, eof, ninv2
    integer, intent(out) :: p_mol, cloudstate
    integer, dimension(nmol), intent(out) :: inv
    integer, dimension(nmol), intent(out) :: arr_inv
    real, intent(out) :: mui, mue, lat
    real, dimension(nmol), intent(out) :: vmr

!=============================================================================
    open(15,file=file_in,status='old',form='formatted')
    read(15,'(F7.3)')lat ; read(15,'(F7.3)')mui ; read(15,'(F7.3)')mue
    read(15,*)
    read(15,'(A7, I2)')cloud, cloudstate
    print*, 'cloudinv', cloudstate
    read(15,*)
    
    do l=1, nmol
       read(15,'(A7,E10.3,I2)',iostat=eof)gas, vmr(l), inv(l)
       if (eof /= 0) exit
       file_opa(l) = file_root // trim(adjustl(gas)) // '.opa'
       p_mol = l
    end do
    close(15)
  ninv2 = count(inv>0)

  do l=1, nmol
       if (inv(l)>0) arr_inv(l) = 1
  end do
  end subroutine read_fuel

  subroutine read_pta(file_in,p_mol,vmr,p,T,profil,cloudstate,cloudinv, taucloud)

    use declaration, only : nmol, nlevel

    implicit none
    character (len=*), intent(in) :: file_in
    integer, intent(in) :: p_mol, cloudstate
    integer :: k,l
    real, dimension(nmol), intent(in) :: vmr
    real, dimension(nlevel), intent(out) :: p, T, cloudinv,taucloud
    real, dimension(nlevel,nmol), intent(out) :: profil

    !====================================================
    open(9,file=file_in,status='old',form='formatted')
    do k=1, nlevel
       read(9,'(E9.2,F6.1,10E9.2)') p(k), T(k), (profil(k,l),l=1,p_mol), taucloud(k) 
    end do

    profil(:,1) = profil(:,1)*1.16667
    ! 1.33333, 1.333333

    do l=1, p_mol
       profil(:,l) = profil(:,l) * vmr(l)
    end do

    if (cloudstate==1) then
      print*, 'reading tropospheric kernel ...'
      open(686,file='nh3_clouds.txt',status='old',form='formatted')
      do k=1, nlevel
         read(686,'(E9.2,F2.1)') p(k), cloudinv(k)
      end do

    elseif (cloudstate==2) then
      print*, 'reading stratospheric kernel ...'
      open(689,file='ch4_hazes.txt',status='old',form='formatted')
      do k=1, nlevel
         read(689,'(E9.2,F2.1)') p(k), cloudinv(k)
      end do
    end if

  end subroutine read_pta

  subroutine read_spe(file_in,wave,spectre,error)

    implicit none
    character (len=*), intent(in) :: file_in
    integer :: nflux,j
    real, dimension(:), allocatable, intent(out) :: wave, spectre, error

!===============================================================    
    open(12,file=file_in,status='old',form='formatted')
    read(12,*) nflux
    allocate(wave(nflux),spectre(nflux),error(nflux))
    do j=1, nflux
       read(12,*) wave(j), spectre(j), error(j)
    enddo
    close(12)

  end subroutine read_spe


  subroutine write_res(file_out,nflux,wave,spectre,synthetic,error)

    implicit none

    integer :: j
    integer, intent(in) :: nflux
    character (len=*), intent(in) :: file_out
    real, dimension(nflux), intent(in) :: wave, spectre, synthetic, error

    open(12,file=file_out,status='unknown',form='formatted')
    do j=1, nflux
       write(12,*) wave(j), spectre(j), synthetic(j), error(j)
    enddo
    close(12)

  end subroutine write_res
  
  subroutine write_inv(file_out,p,t,profil,sigma,cloudstate,taucloud,n_inv,inv_pos,vmr, inv, p_mol)

    use declaration, only: nlevel, nmol
    implicit none

    integer :: j, flag
    character (len=*), intent(in) :: file_out
    integer, intent(in) :: p_mol
    integer, intent(in), dimension(p_mol) :: inv
    integer, intent(in) :: n_inv, cloudstate
    real, dimension(nlevel), intent(in) :: p, t, taucloud
    real, dimension(nlevel,n_inv), intent(in) :: sigma
    integer, dimension(n_inv), intent(in) :: inv_pos
    integer, dimension(n_inv-1) :: inv_pos_red
    real, dimension(nlevel, nmol), intent(in) :: profil
    real, dimension(nmol), intent(in) :: vmr

	flag = 0
    	do j=1, nmol
         if (inv(j)==3) then
	      flag = 1
	   end if
	end do
    
    if (flag == 0) then
    	if (cloudstate == 0) then
    		open(12,file=file_out,status='unknown',form='formatted')
    		do j=1, nlevel
       	   	write(12,'(E10.3,F6.1,16E10.3)') p(j),t(j),profil(j,inv_pos)/vmr(inv_pos), sigma(j,:)
    		enddo
    	elseif (cloudstate /= 0) then
        	open(12,file=file_out,status='unknown',form='formatted')
        	do j=1, nlevel
           		write(12,'(E10.2,F6.1,16E10.3)') p(j),t(j),profil(j,inv_pos)/vmr(inv_pos),taucloud(j), sigma(j,:)
        	enddo
    	end if
    elseif (flag == 1) then
    	if (cloudstate == 0) then
    		open(12,file=file_out,status='unknown',form='formatted')
    		do j=1, nlevel
       	   	write(12,'(E10.3,F6.1,16E10.3)') p(j),t(j),profil(j,inv_pos)/vmr(inv_pos),profil(j,2)/vmr(1), sigma(j,:)
    		enddo
    	elseif (cloudstate /= 0) then
        	open(12,file=file_out,status='unknown',form='formatted')
        	do j=1, nlevel
           		write(12,'(E10.2,F6.1,16E10.3)') p(j),t(j),profil(j,inv_pos)/vmr(inv_pos),taucloud(j),profil(j,2)/vmr(1), sigma(j,:)
        	enddo
    	end if
    end if
    close(12)

  end subroutine write_inv
  
  subroutine write_kernel(file_out,n_inv,A)

    use declaration, only: nlevel
    implicit none

    integer i, j, ii, l
    character (len=*), intent(in) :: file_out
    integer, intent(in) :: n_inv
    real, dimension(nlevel,nlevel, n_inv), intent(in) :: A
    integer, dimension(nlevel, n_inv) :: kernel
    real, dimension(nlevel) :: aux

    do ii = 1, n_inv
    	do l = 1, nlevel
    	   aux = A(:,l,ii)
    	   kernel(l,ii) = maxloc(aux, DIM = 1)
    	end do
    end do

    open(192,file=file_out,status='unknown',form='formatted')
    do j=1, nlevel
       write(192,*) j, kernel(j,:)
    enddo
    close(192)

  end subroutine write_kernel

  subroutine Interpolate(wave, fwhm)
	
	implicit none
	integer :: i, a
	real, allocatable :: f(:), k(:)
	real, allocatable :: xVal(:)
	integer :: filesize, inputIndex
	real, dimension(:), intent(in) :: wave
	real, dimension(:), allocatable, intent(out) :: fwhm 
       
	
	open (unit=33, file='resolution.txt', status='old', action='read')
	read(33,*) filesize
	allocate (f(filesize), k(filesize))
	do i=1,filesize
           read(33,*) f(i), k(i)
        enddo
	close(1)
  
 	 allocate(xVal(size(wave)))
 	 xVal = wave
 	 allocate(fwhm(size(wave)))
  
  	 do inputIndex = 1, size(xVal)
   	    do a = 1, size(f)-1
	       if (xVal(inputIndex)<=f(a) .AND. (xVal(inputIndex)>=f(a+1))) then
	          fwhm(inputIndex) = ((xVal(inputIndex) - f(a))/(f(a+1)-f(a)))*(k(a+1)-k(a))+k(a)
               end if
            enddo
         enddo
  end subroutine Interpolate 

end module read_file
