module read_file

  implicit none
  
contains

  subroutine read_fuel(file_in,lat,mui,mue,p_mol,vmr,inv,file_opa, nummol)
    
    use declaration, only : nmol, file_root
    implicit None

    character (len=7) :: gas
    character (len=*), intent(in) :: file_in
    character (len=*), dimension(nmol), intent(out) :: file_opa
    integer :: l, eof
    integer, intent(out) :: p_mol, nummol
    integer, dimension(nmol), intent(out) :: inv
    real, intent(out) :: mui, mue, lat
    real, dimension(nmol), intent(out) :: vmr
!=============================================================================
    open(15,file=file_in,status='old',form='formatted')
    read(15,'(F7.3)')lat ; read(15,'(F7.3)')mui ; read(15,'(F7.3)')mue
    read(15,*)
    read(15,*)
    read(15,*)
    
    do l=1, nmol
       read(15,'(A7,E10.3,I2)',iostat=eof)gas, vmr(l), inv(l)
       if (inv(l)==2) nummol = l
       if (eof /= 0) exit
       file_opa(l) = file_root // trim(adjustl(gas)) // '.opa'
       p_mol = l
    end do
    close(15)

  end subroutine read_fuel

  subroutine read_pta(file_in,p_mol,vmr,p,T,profil)

    use declaration, only : nmol, nlevel

    implicit none
    character (len=*), intent(in) :: file_in
    integer, intent(in) :: p_mol
    integer :: k,l
    real, dimension(nmol), intent(in) :: vmr
    real, dimension(nlevel), intent(out) :: p, T
    real, dimension(nlevel,nmol), intent(out) :: profil

    !====================================================
    open(9,file=file_in,status='old',form='formatted')
    do k=1, nlevel
       read(9,'(E9.2,F6.1,10E9.2)') p(k), T(k), (profil(k,l),l=1,p_mol)
    end do

    do l=1, p_mol
       profil(:,l) = profil(:,l) * vmr(l)
    end do

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
  
  subroutine write_inv(file_out,p,profil,sigma_vmr, nummol, vmr, t, sigma, inv)

    use declaration, only: nlevel, nmol
    implicit none

    integer :: j, molecule
    character (len=*), intent(in) :: file_out
    integer, intent(in) :: nummol
    integer, dimension(nmol), intent(in) :: inv
    real, dimension(nlevel), intent(in) :: p, sigma, t
    real, dimension(nlevel, nmol), intent(in) :: profil
    real, dimension(nmol), intent(in) :: vmr
    real, dimension(nlevel,count(inv>1)), intent(in) :: sigma_vmr
    
    molecule = 5
    
    open(12,file=file_out,status='unknown',form='formatted')
    do j=1, nlevel
       write(12,'(E10.2,F6.1,F6.1,5E10.2,5E10.2)') p(j), t(j),sigma(j), profil(j,molecule)/vmr(molecule)
    enddo
    close(12)

  end subroutine write_inv
  
  subroutine write_kernel(file_out,A)

    use declaration, only: nlevel
    implicit none

    integer i, j
    character (len=*), intent(in) :: file_out
    real, dimension(nlevel,nlevel), intent(in) :: A

    open(12,file=file_out,ACTION = 'write')
       do i=1, nlevel
         write(12, '(1000E13.5)')( real(A(i,j)) ,j=1,nlevel)
       end do
    close(12)

  end subroutine write_kernel

  subroutine Interpolate(wave, fwhm)
	
	implicit none
	integer :: i, a
	real, allocatable :: f(:), k(:)
	real, allocatable :: xVal(:)
	integer :: filesize, inputIndex
	real, dimension(:), intent(in) :: wave
	real, dimension(:), allocatable, intent(out) :: fwhm 
       
	
	open (unit=33, file='jwst_resol.txt', status='old', action='read')
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
