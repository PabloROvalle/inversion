module info_content

  implicit none

contains
  
  subroutine info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,nflux,wave,synthetic,kk, tau_aer)

    use declaration, only : nfreq, nrep, nlevel, nmol, nt, T_tab, gnu00,dgnu,dgnu_H2, dens, planck, rg, avo, cvel, boltz, H2, He
    use gravity
    use convolution

    implicit none

    integer i, j, j1, k, k1, l, it, rec, nin, n, etat
    integer, intent(in) :: p_mol, nflux
    integer :: repmin, repmax

    character (len=68), dimension(nmol), intent(in) :: file_opa

    real, intent(in) :: lat, mue
    real :: hc2, hck, dt
    real, dimension(2000) :: ftab
    real, dimension(2000,nlevel) :: tabH2
    real, dimension(nlevel,nmol), intent(in) ::  profil
    real, dimension(nlevel,p_mol) ::  colonne
    real, dimension(nflux,nlevel), optional :: kk
    real, dimension(nlevel), intent(in) :: p, T, tau_aer
    real, dimension(nlevel) :: aux_T, grav, z, nz, f
    real, dimension(nflux), intent(in) :: wave
    real, dimension(nflux), intent(inout) :: synthetic
    real, dimension(:,:), allocatable ::  tau, kappa
    real, dimension(:), allocatable ::  fnu, aux, radiance, auxR
    real, dimension(nfreq,2) :: coef
	real, dimension(nflux), intent(in) :: fwhm
!------------------------------------------------------------------------------

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
    do k = 1, nlevel-1
       call cont_h2he(p,T,H2,He,dens,k,grav(k),ftab,nin,tabH2)
    enddo
    allocate(tau((repmin-1)*nfreq+1:repmax*nfreq,nlevel))
    allocate(kappa((repmin-1)*nfreq+1:repmax*nfreq,nlevel))
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

!----------------------------------------------------------------------
             if (l == 5) then  ! ch4, ch3d, nh3, ph3, c2h2, c2h6
                kappa(j1+1:j1+nfreq,k) = (coef(:,1)**(1-dt))&
                  &*(coef(:,2)**dt)*colonne(k,l)
             end if
!----------------------------------------------------------------------

             tau(j1+1:j1+nfreq,k) = tau(j1+1:j1+nfreq,k) + (coef(:,1)**(1-dt))&
                  &*(coef(:,2)**dt)*colonne(k,l)
          enddo
       enddo
    enddo
    do l = 1, p_mol
       close(7+l)
    enddo

!!!!! Hazes !!!!!!

!    do k=1, nlevel-1
!       tau(:,k)=tau(:,k)+tau_aer(k)
!    end do
!    print*, tau_aer
!!!!!!!!!!!!!!!!!!

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

!**********************************************************************
    do i=1,nlevel
       nz(i) = i-1
    end do

    do i=1,nlevel
    	f(i) = (TANH((nz(i)-220)/(27*0.5))+1)/2
    end do
    f = abs(1-f)
!**********************************************************************

    if (.not. present(kk)) goto 20
!==============================================Derivee temperature
    deallocate(fnu)
    allocate(fnu((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    fnu = gnu00 + dgnu*(/(j, j=(repmin-1)*nfreq,repmax*nfreq-1,10)/)
    allocate(radiance((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    allocate(auxR((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    radiance = 0. 
    do k = 1, nlevel-1
!        auxR = (hc2*fnu**3)*(1./(exp(fnu/(aux_T(k)+dt))-1.)-1./(exp(fnu/(aux_T(k)))-1.))*tau(:,k)
!        auxR = (((hc2*fnu**3) / (exp(fnu*hck/T(k+1)) - 1.)) - ((hc2*fnu**3) / (exp(fnu*hck/T(k)) - 1.))) * tau(:,k) !for C2H2
        auxR = (hc2*fnu**3)*(1./(exp(fnu/(aux_T(k+1)))-1.)-1./(exp(fnu/(aux_T(k)))-1.)) * tau(:,k)
        radiance = radiance + auxR
!        call convol_jwst(size(radiance),radiance*kappa((repmin-1)*nfreq+1:repmax*nfreq:10,k),fwhm,fnu,nflux,wave,kk(:,k))
        call convol_jwst(size(radiance),radiance*kappa(:,k)/mue,fwhm,fnu,nflux,wave,kk(:,k))
    end do

!**************************
!    do k = 1, nlevel-1
!       kk(:,k) = f(k)*kk(:,k)
!    enddo
!**************************

    deallocate(radiance)

20  deallocate(fnu)
    allocate(fnu((repmin-1)*nfreq+1:repmax*nfreq),stat=etat)
    fnu = gnu00 + dgnu*(/(j, j=(repmin-1)*nfreq,repmax*nfreq-1)/)
    allocate(radiance((repmin-1)*nfreq+1:repmax*nfreq))
    radiance = 0.
    do k = 1, nlevel-1
       radiance = radiance + (hc2*fnu**3) / (exp(fnu/aux_T(k)) - 1.) * (tau(:,k+1) - tau(:,k))
    enddo
    call convol_jwst(size(radiance),radiance,fwhm,fnu,nflux,wave,synthetic)
    deallocate(radiance)
    deallocate(fnu)

  end subroutine info_content_ch4


  subroutine coeur_ch4(dvmr,sigma,p,nflux,kk,spec_obs,spec_syn,error,A,dr)

    use declaration
    implicit none

    integer :: i, j, k
    integer, intent(in) :: nflux

    Real, parameter :: c=1, factor=30e-1
    real :: alpha, beta, dr
    real :: trace_temp, trace_haze, trace_error
    real, dimension(nflux), intent(in) :: error, spec_obs, spec_syn
    real, dimension(nlevel), intent(out) :: dvmr, sigma
    real, dimension(nlevel), intent(in) :: p
    real, dimension(nflux,nlevel), intent(in) :: kk
    real, dimension(nlevel,nlevel) :: s
    real, dimension(nlevel,nflux) ::  w, aux3
    real, dimension(1,nflux) ::  v
    real, dimension(nflux,nlevel) :: aux
    real, dimension(nflux,nflux) :: aux1, aux2
    real, dimension(nlevel,nlevel), intent(out) :: A
!------------------------------------------------------------------------------
    do i=1,nlevel
       do j=1,nlevel
          s(i,j) = exp( -alog(p(i)/p(j))**2 / (2*c**2))
       enddo
    enddo

    trace_error = sum((/(error(i)**2,i=1,nflux)/) / sqrt(real(nflux,kind=4)))
      
    aux = matmul(kk,s) !--- Temperature
    aux1 = matmul(aux,transpose(kk))
    trace_temp = sum((/(aux1(i,i),i=1,nflux)/))
    alpha = factor * (trace_error/trace_temp)
    
    aux1 = alpha*aux1
    do i=1,nflux
       aux1(i,i) = aux1(i,i) + error(i)**2
    enddo

    call matinv(aux1,nflux,nflux,aux2)

    aux3 = matmul(transpose(kk),aux2)
    w = matmul(s,aux3)
    A = alpha * matmul(w,kk)
    dr = sum((/ (A(i,i), i=1,nlevel) /))
    write(*,*)'coeur', factor, dr
    dvmr = alpha * matmul(w,(spec_obs-spec_syn))

    !-------------------------------------------
!    open(unit=414, file='KERNEL.txt', ACTION="write", STATUS="replace")
!    do i=1, nlevel
!      write(414, '(1000E13.5)')( real(A(i,j)) ,j=1,nlevel)
!    end do
    !-------------------------------------------


    do i=1, nlevel
       sigma(i) = alpha * sqrt(dot_product(w(i,:)**2,error**2))
    end do

  end subroutine coeur_ch4

  
end module info_content
