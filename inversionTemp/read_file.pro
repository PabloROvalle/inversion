subroutine read_fuel(lat,mui,dz,rayleigh,apod_type,v_corr,pds_apriori,H2,He,&
     &tnmol,vmr,fichier)

  use declaration, only : nmol
  implicit None

  character (len=7) :: gas
  character (len=7), intent(out) :: apod_type
  character (len=68), dimension(nmol), intent(out) :: fichier
  character (len=68) :: dump
  integer, intent(out) :: tnmol
  integer :: l
  real, intent(out) :: mui, rayleigh, lat, dz, H2, He, v_corr, pds_apriori
  real, dimension(nmol), intent(out) :: vmr
!=============================================================================
  open(15,file='saturne.fuel',status='old',form='formatted')
  read(15,'(F7.3)')lat ; read(15,'(F7.3)')mui ; read(15,'(F7.3)')dz
  read(15,*) ; read(15,*)
  read(15,'(F7.3)')rayleigh ; read(15,'(A7)') apod_type ;  read(15,*) 
  read(15,'(F7.3)')v_corr ; read(15,'(F7.3)') pds_apriori ;  read(15,*) 
  read(15,*) H2 ; read(15,*) He
  read(15,'(I2)')tnmol
  write(*,'(F7.3)')lat ; write(*,'(F7.3)')mui ; write(*,'(F7.3)')dz
  write(*,*) ; write(*,*)
  write(*,'(F7.3)')rayleigh ; write(*,'(A7)') apod_type ;  write(*,*) 
  write(*,'(F7.3)')v_corr ; write(*,'(F7.3)') pds_apriori ;  write(*,*) 
  write(*,*) H2 ; write(*,*) He
  write(*,'(I2)')tnmol
  do l = 1, tnmol
     read(15,'(A7,E10.3)')gas,vmr(l)
     write(*,'(A7,1PE10.3)')gas,vmr(l)
     dump = '/data/opacite/saturn/saturne_'
     dump = '/home/tfouchet/data/opacite/cirs/saturne_'
     dump(len_trim(dump)+1:len_trim(dump)+7) = gas
!     dump(30:34) = gas
     fichier(l) = dump
  end do
  close(15)

end subroutine read_fuel
