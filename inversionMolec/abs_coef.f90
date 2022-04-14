program abs_coef

  use declaration
  implicit none

  integer :: i, ii, j, k, l, n, r, nraie, nbande, rec, tnmol, err
  integer, dimension(1) :: minfnu, maxfnu

  character (len=14), dimension(nmol) :: spectro, partition
  character (len=38), dimension(nmol) :: fichier
  character (len=7) :: gas
  character (len=38) :: dump
  real :: icoup, zde, zmol, exp_rot, hck, aux_ref, aux3, aux4
  real (kind=8), dimension(nrep*nfreq) :: fnu
  real, dimension(2) :: nu, d
  real, dimension(nlevel) :: p, t
  real, dimension(:), allocatable :: dlor, expo, gnuc, E, dnud, dnul, zi
  real, dimension(nfreq*nrep) :: total, gnu1
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

  open(15,file='jupiter.fuel',status='old',form='formatted')
  read(15,'(F7.3)') ; read(15,'(F7.3)') ; read(15,'(F7.3)')
  read(15,*) ; read(15,*)
  l = 1
  do
     read(15,'(A7)',iostat=err)gas
     if (err == -1) exit
     dump = repeat(' ', 38)
     dump = '/data/opacite/jupiter/jupiter_                     '
     dump(len_trim(dump)+1:len_trim(dump)+8) = gas
     dump(len_trim(dump)+1:len_trim(dump)+5) = '.opa'
     fichier(l) = dump
     dump = gas
     dump(len_trim(dump)+1:len_trim(dump)+7) = '.geisa'
     spectro(l) = dump(1:14)
     dump = gas
     dump(len_trim(dump)+1:len_trim(dump)+7) = '.partf'
     partition(l) = dump(1:14)
     print*,fichier(l)
     print*,spectro(l)
     print*,partition(l)
     l = l + 1
  end do
  tnmol = l-1 
  close(15)

!----Profil pta
  open(9,file='jupiter.pta',status='old',form='formatted')
  do k=1, nlevel
     read(9,'(E9.2,F6.1)') p(k), t(k)
  enddo

!----Calcul coefficients d'absorption
  do l = 1,tnmol
     open(7+l,file=fichier(l),access='direct',recl=nfreq,status='new',&
          form='unformatted')
     print*, fichier(l)
  enddo

  aux_ref =  -hck / Tref
  do k=1, nlevel
     do l=1, tnmol
!----------Lecture des raies et fonction de partition-------------------------
        open(14, file=partition(l),status='old',form='formatted')
        read(14,*) zmol
        read(14,*) exp_rot
        read(14,*) icoup
        read(14,*) nbande
        if (nbande > 2) then
           print*, 'Probleme avec Qv'
           stop
        endif
        do i=1,nbande
           read(14,*) nu(i), d(i)
        enddo
        close(14)

        open(10,file=spectro(l),status='old',form='formatted')
        r = 0
        do 
           Read(10,*,iostat=err)
           if (err == -1) exit
           r = r + 1
        enddo
        close(10)
        nraie = r
        allocate(zi(nraie),gnuc(nraie),dlor(nraie),E(nraie),expo(nraie))
        allocate(dnud(nraie),dnul(nraie))

        open(10,file=spectro(l),status='old',form='formatted')
        do r=1, nraie
           Read(10,'(F10.3,E10.3,F5.3,F10.3,36x,F4.2)') gnuc(r), zi(r), &
                dlor(r), E(r), expo(r)
        enddo
        close(10)
        aux3 = -hck * (1./T(k)-1./Tref)
        aux4 = -hck / T(k)
        zi = zi * 2.6868e19 &
             * (Tref/T(k))**exp_rot * exp(E*aux3) &! Partition rota.
             * (1.-exp(gnuc*aux4)) / (1.-exp(gnuc*aux_ref)) !Emis. Ind.
        do ii=1,nbande ! Partition vibra.
           zi = zi * ((1.-exp(nu(ii)*aux4)) &
                / (1.-exp(nu(ii)*aux_ref)))**d(ii)
        enddo
        dnul = (dlor*((Tref/T(k))**expo)) * (p(k)/1.01325)
        dnud = gnuc * ((1.e+02/cvel) * SQRT(2.*rg*T(k)/zmol))

        print*, nraie, spectro(l), partition(l)

        fnu = gnu00 + dgnu*(/(j-1, j=1,nrep*nfreq)/)
        total = 0.
        do r=1, nraie
           gnu1 = abs(fnu-gnuc(r))
           maxfnu = minloc(abs(fnu-gnuc(r)-icoup))
           minfnu = minloc(abs(fnu-gnuc(r)+icoup))
           total(minfnu(1):maxfnu(1)) = total(minfnu(1):maxfnu(1)) &
                + voigt(gnu1(minfnu(1):maxfnu(1))/dnud(r),&
                dnul(r)/dnud(r)) * (zi(r)/dnud(r))
        enddo
        rec = (k-1) * nrep
        do n = 1, nrep
           write(7+l,rec=rec+n) total(1+(n-1)*nfreq:n*nfreq)
        enddo
        deallocate(zi,gnuc,dlor,E,expo)
        deallocate(dnud,dnul)

     enddo
  enddo
      
  do l = 1,tnmol
     close(7+l)
  enddo
  
end program abs_coef
