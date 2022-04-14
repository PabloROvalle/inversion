program inverse

  use declaration
  implicit None
  integer :: i, j, k, l, tnmol, err
  integer :: nflux

  logical :: calc

  character (len=7) :: gas
  character (len=38), dimension(nmol) :: fichier
  character (len=38) :: dump
  real :: mui, mue, rayleigh, chi2, lat
  real, dimension(:,:), allocatable :: kk
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t, dt
  real, dimension(:), allocatable :: wavesup, waveinf, spectre, nesr, flux
  interface
     subroutine info_content(p,t,profil,tnmol,fichier,nflux,rayleigh,wavesup,&
          &waveinf,spectre,kk,calc)
       use declaration
       implicit none

       integer, intent(in) :: nflux, tnmol
       character (len=38), dimension(nmol), intent(in) :: fichier
       logical, intent(in) :: calc
       real, intent(in) :: rayleigh
       real, dimension(nlevel,nmol), intent(in) ::  profil
       real, dimension(nflux,nlevel), intent(inout) :: kk
       real, dimension(nlevel), intent(in) :: p, T
       real, dimension(nflux), intent(in) :: wavesup, waveinf
       real, dimension(nflux), intent(inout) :: spectre
     end subroutine info_content
     subroutine coeur(dt,p,kk,error,spectre,flux,nflux)
       use declaration
       implicit none
       integer, intent(in) :: nflux
       real, dimension(nflux), intent(in) :: error, flux, spectre
       real, dimension(nlevel), intent(inout) :: dt
       real, dimension(nlevel), intent(in) :: p
       real, dimension(nflux,nlevel), intent(in) :: kk
     end subroutine coeur
  end interface
!==============================================================================
  open(15,file='jupiter.fuel',status='old',form='formatted')
  read(15,'(F7.3)')mui ; read(15,'(F7.3)')mue
  read(15,*) ; read(15,'(F7.3)')rayleigh
  read(15,*)
  l = 1
  do
     read(15,'(A7)',iostat=err)gas
     if (err == -1) exit
     dump = repeat(' ', 38)
     dump = '/data/opacite/jupiter/jupiter_                     '
     dump(len_trim(dump)+1:len_trim(dump)+8) = gas
     dump(len_trim(dump)+1:len_trim(dump)+5) = '.opa'
     fichier(l) = dump
     l = l + 1
  end do
  tnmol = l-1 
  close(15)

!----Profil pta
  open(9,file='jupiter.pta',status='old',form='formatted')
  do k=1, nlevel
     read(9,'(E9.2,F6.1,6E9.2)') p(k), T(k), (profil(k,l), l=1,tnmol)
  enddo

!----Spectre a inverser
  open(12,file='jupiter.spe',status='old',form='formatted')
  read(12,*) nflux
  allocate(flux(nflux),spectre(nflux))
  allocate(wavesup(nflux),waveinf(nflux),nesr(nflux))
  do j=1, nflux
     read(12,*) waveinf(j), wavesup(j), flux(j), nesr(j) 
  enddo
  close(12)
!----Information content
  calc = .false.
  allocate(kk(nflux,nlevel))
  kk = 0.
  t(0:280) = t(0:280) - 10
  call info_content(p,t,profil,tnmol,fichier,nflux,rayleigh,wavesup,waveinf,&
       &spectre,kk,calc)
!  open(12,file='spectre.txt',form='formatted',status='old')
!  read(12,*)spectre
!  close(12)
!  write(*,*)ztan(4)
!  write(*,*)spectre(1+nflux4*3:nflux4*4)
!  stop
!!$  write(*,*)'1bis'
!!$  open(12,file='save.txt',form='formatted',status='new')
!!$  write(12,*)spectre
!!$  write(12,*)kk
!!$  close(12)
!!$  stop
!!$  do k=121,nlevel
!!$     kk(:,k+nlevel) = kk(:,k+nlevel)*(p(k)/1e-2)
!!$  end do
!!$  write(*,*)size(kk),size(spectre),size(flux),size(nesr),nalt3,nflux3,nalt4
!!$  do k=1,120
!!$     kk(:,k+nlevel) = kk(:,k+nlevel)*(1e-2/p(k))
!!$  end do
  write(*,*) 0, sum(((spectre-flux)/nesr)**2.)
!----Iteration
  do i=1,3
     call coeur(dt,p,kk,nesr,spectre,flux,nflux)
     t = t + dt
     if (i == 3) calc = .false.
     call info_content(p,t,profil,tnmol,fichier,nflux,rayleigh,wavesup,&
          &waveinf,spectre,kk,calc)
     chi2 = sum(((spectre-flux)/nesr)**2.)
     write(*,*)i, chi2
  enddo

!---- Grille d'altitude et de gravite
  open(12,file='jupiter.inv',form='formatted',status='unknown')
  do k=1,nlevel
     write(12,*)p(k), t(k)
  enddo
  close(12)

10 open(17,file='jupiter.res',form='formatted',status='unknown')
  do j=1, nflux
     write(17,*) waveinf(j), wavesup(j), flux(j), spectre(j)
  enddo
  close(17)

end program inverse