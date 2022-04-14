module info_content

  implicit none

contains
  
  subroutine info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,nflux,wave,synthetic)

    use declaration, only : nfreq, nrep, nlevel, nmol, nt, T_tab, gnu00,dgnu,dgnu_H2, dens, planck, rg, avo, cvel, boltz, H2, He
    use gravity
    use convolution

    implicit none

    integer i, j, j1, k, k1, l, it, rec, nin, n
    integer, intent(in) :: p_mol, nflux
    integer :: repmin, repmax

    character (len=68), dimension(nmol), intent(in) :: file_opa

    logical :: calc
    
    real, intent(in) :: lat, mue
    real :: hc2, hck, dt
    real, dimension(2000) :: ftab
    real, dimension(2000,nlevel) :: tabH2
    real, dimension(nlevel,nmol), intent(in) ::  profil
    real, dimension(nlevel,p_mol) ::  colonne
    !  real, dimension(:,:) :: kk
    real, dimension(nlevel), intent(in) :: p, T
    real, dimension(nlevel) :: aux_T, grav, z
    real, dimension(nflux), intent(in) :: wave
    real, dimension(nflux), intent(inout) :: synthetic
    real, dimension(:,:), allocatable ::  tau, radiance
    real, dimension(:), allocatable ::  fnu, aux
    real, dimension(nfreq,2) :: coef
	real, dimension(nflux), intent(in) :: fwhm

!------------------------------------------------------------------------------
    calc=.false.

!---- Domaine de longueur d'onde et nom des fichiers

!----Bornes des calculs
    do repmin=2,nrep
       if((gnu00+repmin*nfreq*dgnu) > (wave(1)-3*fwhm(1))) exit
    enddo
    repmin = repmin - 1
    do repmax=nrep-1, 1, -1
       if((gnu00+repmax*nfreq*dgnu) < (wave(nflux)+3*fwhm(nflux))) exit
    enddo
    repmax = repmax + 1
!---- Grille d'altitude et de gravite
    z(41) = 0.
    do k=41, nlevel-1
       call newgrav(lat,z(k)/1e3,grav(k))
       z(k+1) = z(k) + (rg*(T(k)+T(k+1))/2./dens/grav(k)) * alog(p(k)/p(k+1))
    enddo
    call newgrav(lat,z(nlevel)/1e3,grav(nlevel))
    do k=41, 2, -1
       call newgrav(lat,z(k)/1e3,grav(k))
       z(k-1) = z(k) + (rg*(T(k)+T(k-1))/2./dens/grav(k)) * alog(p(k)/p(k-1))
    enddo
    call newgrav(lat,z(1)/1e3,grav(1))
!----Colonne densite
    colonne = reshape(p,(/nlevel,p_mol/),pad=p)
    colonne = (colonne-eoshift(colonne,shift=1)) * profil
    colonne = (avo*1e5/2.6868e19/1e4/dens) * colonne &
         / reshape(grav,(/nlevel,p_mol/),pad=grav) 
    !----Epaissseur optique
    nin = ceiling(nfreq*(repmax-repmin+1)*dgnu / dgnu_H2) + 2
    print*, nin,'H2'
    if (nin > 2000) stop
    ftab(1:nin) = gnu00 + ((repmin-1)*nfreq)*dgnu + (/(i-2,i=1,nin)/)*dgnu_H2
    do k = 1,nlevel-1
       call cont_h2he(p,T,H2,He,dens,k,grav(k),ftab,nin,tabH2)
    enddo
    allocate(tau((repmin-1)*nfreq+1:repmax*nfreq,nlevel))
    allocate(fnu((repmin-1)*nfreq+1:repmax*nfreq))
    fnu = gnu00 + dgnu*(/(j, j=(repmin-1)*nfreq,repmax*nfreq-1)/)

    i = 1
    do j = (repmin-1)*nfreq+1, repmax*nfreq
       do k=i+1,i,-1
          if (ftab(k) < fnu(j)) exit
       enddo
       i = k
       tau(j,:) =  tabH2(i,:) + ((fnu(j)-ftab(i))/(ftab(i+1)-ftab(i))) &
            * (tabH2(i+1,:)-tabH2(i,:))
    enddo
!---- et les autres molecules
    do l = 1, p_mol
       write(*,*) file_opa(l)
       open(7+l,file=file_opa(l),access='direct',recl=nfreq,status='old',&
            form='unformatted')
    enddo
    do k = 1, nlevel
       do it = 1, nt
          if (T(k)<t_tab(it,k)) exit
       enddo
       it = it - 1
       if (it==0) then
          dt = 0.
          it = 1
          write(*,*)'!!!!Attention T < T_tab, niveau: ', k
       else if (it==nt) then
          dt = 1.
          it = nt - 1
          write(*,*)'!!!!Attention T > T_tab, niveau: ', k
       else
          dt = (T(k) - t_tab(it,k)) / (t_tab(it+1,k)-t_tab(it,k))
       end if
       rec = (it + (k-1)*nt - 1) * nrep
!         write(*,*)k,it,rec
       do l = 1,p_mol
          do n = repmin, repmax
             j1 = (n-1) * nfreq
             read(7+l, rec=rec+n) coef(:,1)
             read(7+l, rec=rec+n+nrep) coef(:,2)
             tau(j1+1:j1+nfreq,k) = tau(j1+1:j1+nfreq,k) + (coef(:,1)**(1-dt))&
                  &*(coef(:,2)**dt)*colonne(k,l)
          enddo
       enddo
    enddo
    do l = 1, p_mol
       close(7+l)
    enddo

    do k = nlevel-1, 1, -1
       tau(:,k) = tau(:,k+1) + tau(:,k)
    end do

    do k = 1, nlevel
       tau(:,k) = exp(-tau(:,k)/mue)
    end do
    
    hc2 = 2 * planck * cvel**2
    hck = planck * cvel / boltz
    aux_T = 0.5 * (T+eoshift(T,1,T(size(T)))) / hck
    dt = 1. / hck
    
    if (.not. calc) goto 20
!==============================================Derivee temperature
    deallocate(fnu)
    allocate(fnu((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    fnu = gnu00 + dgnu*10*(/(j, j=(repmin-1)*nfreq/10,repmax*nfreq/10-1)/)
    allocate(radiance((repmin-1)*nfreq/10+1:repmax*nfreq/10,nlevel))
    allocate(aux((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    radiance = 0. 
    do k = 1, nlevel-1
       aux = (hc2*fnu**3) / (exp(fnu/aux_T(k)) - 1.) * (tau(:,k+1) - tau(:,k))
       do k1 = 1, nlevel
          if (k1 == k) then
             radiance(:,k1) = radiance(:,k1) + (hc2*fnu**3) &
                  / (exp(fnu/(aux_T(k1)+dt)) - 1.) * (tau(:,k+1) - tau(:,k))
          else
             radiance(:,k1) = radiance(:,k1) + aux
          endif
       enddo
    enddo
!  call convol_jwst(size(radiance(:,nlevel)),radiance(:,nlevel),fwhm,fnu,nflux,wave,synthetic)
    do k1=1,nlevel-1
       !     call convol_jwst(radiance(:,k1),n_fnu,fnu,nflux,wave,&
       !          &kk(:,k1))
       !     kk(:,k1) = kk(:,k1) - synthetic
    enddo
    deallocate(radiance)
    deallocate(aux)

20  deallocate(fnu)
    allocate(fnu((repmin-1)*nfreq+1:repmax*nfreq))
    fnu = gnu00 + dgnu*(/(j, j=(repmin-1)*nfreq,repmax*nfreq-1)/)
    allocate(radiance((repmin-1)*nfreq+1:repmax*nfreq,1))
    allocate(aux((repmin-1)*nfreq+1:repmax*nfreq))
    radiance = 0. 
    do k = 1, nlevel-1
       aux = (hc2*fnu**3) / (exp(fnu/aux_T(k)) - 1.) * (tau(:,k+1) - tau(:,k))
       radiance(:,1) = radiance(:,1) + aux
    enddo
    call convol_jwst(size(radiance(:,1)),radiance(:,1),fwhm,fnu,nflux,wave,synthetic)
    deallocate(radiance)
    deallocate(aux)

  end subroutine info_content_ch4

end module info_content
