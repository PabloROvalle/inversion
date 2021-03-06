program abs_coef

  use declaration, only : planck, cvel, boltz, nmol, nlevel, Tref, nt, rg, T_tab, nrep, nfreq, gnu00, dgnu
  use read_file
  
  implicit none

  integer :: l, err, tnmol, k, nbande, i, r, nraie, ii, j, it, n, rec
  integer, dimension(1) :: minfnu, maxfnu

  character (len=*), parameter :: file_fuel='jupiter_ch4.fuel'
  character (len=*), parameter :: file_pta='jupiter_ch4.pta'
  

  character (len=100), dimension(nmol) :: file_opa, file_partf, file_spectro
  character (len=100) :: dump
  integer :: p_mol, pos
  real :: mui, mue, fwhm, lat, dens
  real :: hck, aux_ref, zmol, exp_rot, icoup, aux3, aux4!, zde
  real, dimension(nlevel,nmol) :: profil
  real, dimension(nlevel) :: p, t, tau_haze, z, grav
  real, dimension(nmol) :: vmr
  real (kind=8), dimension(nrep*nfreq) :: fnu
  real, dimension(nrep*nfreq) :: gnu1, total
  real, dimension(7) :: nu, d
  real, dimension(:), allocatable :: gnuc, dlor, E, expo
  real, dimension(:,:), allocatable :: zi, dnud, dnul
  interface
     function voigt(x,y)
       implicit none
       real, dimension(:), intent(in) :: x
       real, intent(in) :: y
       real, dimension(size(x)) :: voigt
    end function voigt
  end interface
!--------------------------------------------------------------------------
  hck = planck * cvel / boltz

  
  call read_fuel(file_fuel,lat,mui,mue,fwhm,p_mol,vmr,file_opa)
  call read_pta(file_pta,p_mol,vmr,p,T,profil)

  do l=1,p_mol
     pos=scan(file_opa(l), '.')
     dump = file_opa(l)
     file_partf(l) = dump(1:pos) // 'partf'
     file_spectro(l) = dump(1:pos) // 'geisa'

  end do

  
!----Calcul coefficients d'absorption
  do l = 1,p_mol
     open(7+l,file=file_opa(l),access='direct',recl=nfreq,status='new',&
          form='unformatted')
     write(*,*) file_opa(l)
  enddo

  aux_ref =  -hck / Tref
  do k=1, nlevel
     do l=1, p_mol
!----------Lecture des raies et fonction de partition-------------------------
        
        open(14, file=file_partf(l),status='old',form='formatted')
        read(14,*) zmol
        read(14,*) exp_rot
        read(14,*) icoup
        read(14,*) nbande
        if (nbande > 8) then
           write(*,*) 'Probleme avec Qv'
           stop
        endif
        do i=1,nbande
           read(14,*) nu(i), d(i)
        enddo
        close(14)
 

        open(10,file=file_spectro(l),status='old',form='formatted')
        r = 0
        do 
           Read(10,*,iostat=err)
           if (err == -1) exit
           r = r + 1
        enddo
        close(10)
        print*, r
        nraie = r
        print*, nraie
        allocate(zi(nraie,nt),dnud(nraie,nt),dnul(nraie,nt))
        allocate(gnuc(nraie),dlor(nraie),E(nraie),expo(nraie))!

        open(10,file=file_spectro(l),status='old',form='formatted')
        do r=1, nraie
           Read(10,'(F10.3,E10.3,F5.3,F10.3,36x,F4.2)') gnuc(r), zi(r,1), &
                dlor(r), E(r), expo(r)
        enddo
        close(10)
        do i=nt, 1, -1
           aux3 = -hck * (1./T_tab(i,k)-1./Tref)
           aux4 = -hck / T_tab(i,k)
           zi(:,i) = zi(:,1) * 2.6868e19 &
                * (Tref/T_tab(i,k))**exp_rot * exp(E*aux3) &! Partition rota.
                * (1.-exp(gnuc*aux4)) / (1.-exp(gnuc*aux_ref)) !Emis. Ind.
           do ii=1,nbande ! Partition vibra.
              zi(:,i) = zi(:,i) * ((1.-exp(nu(ii)*aux4)) &
                   / (1.-exp(nu(ii)*aux_ref)))**d(ii)
           enddo
           dnul(:,i) = dlor * ((Tref/T_tab(i,k))**expo) * (p(k)/1.01325)
           dnud(:,i) = gnuc * ((1.e+02/cvel) * SQRT(2.*rg*T_tab(i,k)/zmol))
        enddo

        fnu = gnu00 + dgnu*(/(j-1, j=1,nrep*nfreq)/)
        do it=1, nt
           total = 0.
           do r=1, nraie
              gnu1 = abs(fnu-gnuc(r))
              maxfnu = minloc(abs(fnu-gnuc(r)-icoup))
              minfnu = minloc(abs(fnu-gnuc(r)+icoup))
              !write(*,*)k,l,it,r,nraie, minfnu, maxfnu
              total(minfnu(1):maxfnu(1)) = total(minfnu(1):maxfnu(1)) &
                   + voigt(gnu1(minfnu(1):maxfnu(1))/dnud(r,it),&
                   dnul(r,it)/dnud(r,it)) * (zi(r,it)/dnud(r,it))
           enddo
           
           rec = (it + (k-1)*nt -1) * nrep
           do n = 1, nrep
              write(7+l,rec=rec+n) total(1+(n-1)*nfreq:n*nfreq)
           enddo
        enddo

        deallocate(zi,dnud,dnul)
        deallocate(gnuc,dlor,E,expo)

     enddo
  enddo
      
  do l = 1,tnmol
     close(7+l)
  enddo

end program abs_coef
