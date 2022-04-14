module info_content

  implicit none

contains
  
  subroutine info_content_ch4(p,t,profil,p_mol,lat,mue,fwhm,file_opa,nflux,wave,synthetic,n_inv,inv,imol_inv,xm,xt,kk)

    use declaration, only : nfreq, nrep, nlevel, nmol, nt, T_tab, gnu00,dgnu,dgnu_H2, dens, planck, rg, avo, cvel, boltz, H2, He
    use gravity
    use convolution

    implicit none

    integer :: i, j, j1, k, l, it, nin, rec, n, etat, ii
    integer, intent(out) :: xm,xt
    integer, intent(in) :: p_mol, nflux, n_inv
    integer, intent(in), dimension(p_mol) :: inv
    integer :: repmin, repmax
    integer, intent(out) :: imol_inv

    character (len=*), dimension(nmol), intent(in) :: file_opa

    real, intent(in) :: lat, mue
    real :: hc2, hck, dt
    real, dimension(2000) :: ftab
    real, dimension(2000,nlevel) :: tabH2
    real, dimension(nlevel,nmol), intent(in) ::  profil
    real, dimension(nlevel,p_mol) ::  colonne
    real, dimension(nflux,nlevel,n_inv), optional :: kk
    real, dimension(nlevel), intent(in) :: p, T
    real, dimension(nlevel) :: aux_T, grav, z, nz, f
    real, dimension(nflux), intent(in) :: wave, fwhm
    real, dimension(nflux), intent(inout) :: synthetic
    real, dimension(:,:,:), allocatable ::  kappa
    real, dimension(:,:), allocatable ::  tau
    real, dimension(:), allocatable ::  fnu, aux, radiance
    real, dimension(nfreq,2) :: coef
!------------------------------------------------------------------------------
    xt = count(inv>0)
    xm = count(inv>1)
!---- Domaine de longueur d'onde et nom des fichiers

!----Bornes des calculs
    do repmin=2,nrep
       if((gnu00+repmin*nfreq*dgnu) > (minval(wave)-3*fwhm(1))) exit
    enddo
    repmin = repmin - 1
    do repmax=nrep-1, 1, -1
       if((gnu00+repmax*nfreq*dgnu) < (maxval(wave)+3*fwhm(nflux))) exit
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
    if (maxval(inv) == 2) allocate(kappa((repmin-1)*nfreq+1:repmax*nfreq,nlevel,xm))
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
       ii = 1
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

       do l = 1,p_mol
         do n = repmin, repmax
            j1 = (n-1) * nfreq
            read(7+l, rec=rec+n) coef(:,1)
            read(7+l, rec=rec+n+nrep) coef(:,2)
            if (inv(l) == 2) then
                kappa(j1+1:j1+nfreq,k,ii) = (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)
	    end if
            tau(j1+1:j1+nfreq,k) = tau(j1+1:j1+nfreq,k) + (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)
          enddo
       if (inv(l)==2) then
          ii = ii+1
       endif
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

!******************** CUTOFF UPPER ATMOSPHERE *************************
    do i=1,nlevel
       nz(i) = i-1
    end do

    do i=1,nlevel
    	f(i) = (TANH((nz(i)-220)/(27*0.5))+1)/2
    end do
    f = abs(1-f)
!**********************************************************************

!==============================================Calcul spectre
    allocate(radiance((repmin-1)*nfreq+1:repmax*nfreq))
    radiance = 0.
    do k = 1, nlevel-1
       radiance = radiance + (hc2*fnu**3) / (exp(fnu/aux_T(k)) - 1.) * (tau(:,k+1) - tau(:,k))
    enddo
    call convol_jwst(size(radiance),radiance,fwhm,fnu,nflux,wave,synthetic)
    deallocate(radiance)
    deallocate(fnu)
    ii = 1

    if (n_inv == 0) return
    
    allocate(fnu((repmin-1)*nfreq/10+1:repmax*nfreq/10))
    fnu = gnu00 + dgnu*(/(j, j=(repmin-1)*nfreq,repmax*nfreq-1,10)/)
    allocate(radiance((repmin-1)*nfreq/10+1:repmax*nfreq/10))

    do l=1, p_mol  
     if (inv(l) == 1) then
!==============================================Derivee temperature
         do k = 1, nlevel-1
            radiance = (hc2*fnu**3) * (1./(exp(fnu/(aux_T(k)+dt))-1.) - 1./(exp(fnu/(aux_T(k)))-1.)) * (tau((repmin-1)*nfreq+1:repmax*nfreq:10,k+1) - tau((repmin-1)*nfreq+1:repmax*nfreq:10,k))
            call convol_jwst(size(radiance),radiance,fwhm,fnu,nflux,wave,kk(:,k,ii))
         end do
         ii = ii+1
      elseif (inv(l) == 2) then
!==============================================Derivee abondance
         imol_inv = 1
         radiance = 0
         do k=1, nlevel-2
            radiance = radiance - (hc2*fnu**3) * tau((repmin-1)*nfreq+1:repmax*nfreq:10,k) *(1./(exp(fnu/aux_T(k))-1.) - 1./(exp(fnu/aux_T(k+1))-1.))
            call convol_jwst(size(radiance),radiance*kappa(::10,k,imol_inv)/mue,fwhm,fnu,nflux,wave,kk(:,k,ii))
         enddo
         imol_inv = imol_inv+1
         ii = ii+1
      endif
     enddo
     

     open(unit=202, file='ContFunc.txt', ACTION="write", STATUS="replace")
     do i=1, nflux
       write(202, '(1000E13.5)')( real(kk(i,j,1)) ,j=1,nlevel)
     end do

     open(unit=272, file='ContFunc2.txt', ACTION="write", STATUS="replace")
     do i=1, nflux
       write(272, '(1000E13.5)')( real(kk(i,j,2)) ,j=1,nlevel)
     end do
  
!**************************
!    do k = 1, nlevel-1
!       kk(:,k,:) = f(k)*kk(:,k,:)
!    enddo
!**************************

  end subroutine info_content_ch4


  subroutine coeur_ch4(dvmr,dt,sigma,sigma_vmr,A,p,nflux,p_mol,inv,xm,xt,kk,spec_obs,spec_syn,error)

    use declaration
    implicit none

    integer :: i, j, k, m
    integer, intent(in) :: p_mol
    integer, intent(in) :: nflux, xm, xt
    integer, intent(in), dimension(p_mol) :: inv

    Real, parameter :: c=0.75, factor=30e-1
    real :: beta
    real :: trace_haze, trace_error
    real, dimension(xt) :: trace_temp, alpha
    real, dimension(nflux), intent(in) :: error, spec_obs, spec_syn
    real, dimension(nlevel), intent(out) :: dt, sigma
    real, dimension(nlevel, xm), intent(out) :: dvmr, sigma_vmr
    real, dimension(nlevel), intent(in) :: p
    real, dimension(nflux,nlevel, xt), intent(in) :: kk
    real, dimension(nlevel,nlevel) :: s
    real, dimension(nlevel,nflux, xt) ::  w, aux3
    real, dimension(nlevel,nlevel, xt), intent(out) :: A
    real, dimension(1,nflux) ::  v
    real, dimension(nflux,nlevel, xt) :: aux
    real, dimension(nflux,nflux, xt) :: aux1, aux2
!------------------------------------------------------------------------------
    do i=1,nlevel
       do j=1,nlevel
          s(i,j) = exp( -alog(p(i)/p(j))**2 / (2*c**2))
       enddo
    enddo

    trace_error = sum((/(error(i)**2,i=1,nflux)/) / sqrt(real(nflux,kind=4)))
 
    do m = 1, xt 
    	aux(:,:,m) = matmul(kk(:,:,m),s(:,:))
    	aux1(:,:,m) = matmul(aux(:,:,m),transpose(kk(:,:,m)))
    	trace_temp(m) = sum((/(aux1(i,i,m),i=1,nflux)/))
    	alpha(m) = factor * (trace_error/trace_temp(m))  
     
    	aux1(:,:,m) = alpha(m)*aux1(:,:,m)
    	do i=1,nflux
       	    aux1(i,i,m) = aux1(i,i,m) + error(i)**2
    	enddo

    	call matinv(aux1(:,:,m),nflux,nflux,aux2(:,:,m))

    	aux3(:,:,m) = matmul(transpose(kk(:,:,m)),aux2(:,:,m))
    	w(:,:,m) = matmul(s,aux3(:,:,m))
        if (m==1) then
      	    dt = alpha(m) * matmul(w(:,:,m),(spec_obs-spec_syn))
        elseif (m>1) then
            dvmr(:,m) = alpha(m) * matmul(w(:,:,m),(spec_obs-spec_syn))
        endif
    	
    	A(:,:,m) = alpha(m) * matmul(w(:,:,m),kk(:,:,m))

    	do i=1, nlevel
            if (m==1) then
       		sigma(i) = alpha(m) * sqrt(dot_product(w(i,:,m)**2,error**2))
            elseif (m>1) then
            	sigma_vmr(i,m) = alpha(m) * sqrt(dot_product(w(i,:,m)**2,error**2))
            endif
    	end do
    end do
  end subroutine coeur_ch4  

end module info_content
