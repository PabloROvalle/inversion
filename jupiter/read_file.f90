module read_file

  implicit none
  
contains

  subroutine read_fuel(file_in,lat,mui,mue,p_mol,vmr,file_opa)
    
    use declaration, only : nmol, file_root
    implicit None

    character (len=7) :: gas
    character (len=50), intent(in) :: file_in
    character (len=68), dimension(nmol), intent(out) :: file_opa
    integer :: l, eof
    integer, intent(out) :: p_mol
    real, intent(out) :: mui, mue, lat
    real, dimension(nmol), intent(out) :: vmr
!=============================================================================
    open(15,file=file_in,status='old',form='formatted')
    read(15,'(F7.3)')lat ; read(15,'(F7.3)')mui ; read(15,'(F7.3)')mue
    read(15,*)
    read(15,*)
	read(15,*)
    
    do l=1, nmol
       read(15,'(A7,E10.3)',iostat=eof)gas, vmr(l)
       if (eof /= 0) exit
       file_opa(l) = file_root // trim(adjustl(gas)) // '.opa'
       p_mol = l
    end do

    close(15)

  end subroutine read_fuel

  subroutine read_pta(file_in,p_mol,vmr,p,T,profil)

    use declaration, only : nmol, nlevel

    implicit none
    character (len=50), intent(in) :: file_in
    integer, intent(in) :: p_mol
    integer :: k,l
    real, dimension(nmol), intent(in) :: vmr
    real, dimension(nlevel), intent(out) :: p, T
    real, dimension(nlevel,nmol), intent(out) :: profil

    !====================================================
    open(9,file=file_in,status='old',form='formatted')
    do k=1, nlevel
       read(9,'(E9.2,F6.1,6E9.2)') p(k), T(k), (profil(k,l),l=1,p_mol)
    end do

    do l=1, p_mol
       profil(:,l) = profil(:,l) * vmr(l)
    end do

  end subroutine read_pta

  subroutine read_spe(file_in,wave,spectre,error)

    implicit none
    character (len=50), intent(in) :: file_in
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


  subroutine write_res(file_out,wave,spectre,synthetic,error)

    implicit none

    integer :: j
    character (len=50), intent(in) :: file_out
    real, dimension(:), intent(in) :: wave, spectre, synthetic, error

    print*,size(wave)
    !open(12,file=file_out,status='unknown',form='formatted')
    open(12,file=file_out,status='unknown')
    do j=1, size(wave)
       write(12,'(F10.5,E19.10,E19.10,E19.10)') wave(j), spectre(j), synthetic(j), error(j)
    enddo
    close(12)

  end subroutine write_res


!PABLO!


  subroutine Interpolate(wave, fwhm)
!--------------------------------------------------------------!
!------------- Defining variables -----------------------------!
	
	implicit none
		integer :: i, a
		real, allocatable :: f(:), k(:)
		real, allocatable :: xVal(:)
		integer :: filesize, inputIndex
		real, dimension(:), intent(in) :: wave
		real, dimension(:), allocatable, intent(out) :: fwhm

!--------------------------------------------------------------!
!------------------- Reading file -----------------------------!
        
	
	open (unit=33, file='jwst_resol.txt', status='old', action='read') 	
	read(33,*) filesize
	allocate (f(filesize), k(filesize))
	do i=1,filesize
        read(33,*) f(i), k(i)
    enddo  
	close(1)

!--------------------------------------------------------------!
!------------- Interpolation routine --------------------------!
  
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


!--------------------------------------------------------------!
!------------- Writing interpolated data ----------------------!

open (unit=99, file='jwst_resol_1.txt', status='unknown', action='write') 	
	do i=1,size(xVal)
        write(99,*) xVal(i), fwhm(i)
    enddo  
close(99)



end subroutine Interpolate 
  
end module read_file
