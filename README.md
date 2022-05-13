# Radiative transfer and inversion code
Here I explain the most basic things of the radiative transfer code, plus the maths behind the inversion algorithm. Finally I show you how to use it with a little example, and how the outputs means and look like.
## Radiative transfer
The radiative transfer code has its basis on GEISA and HITRAN databases. These databases provide high resolution spectroscopy data for several molecules. We use this
data to know the opacities of the molecules in our model. These files use de .opa denomination, and are not included due to its large size of tens of GB for each molecule (for more info about the format and info in this files ask to the mail in my welcome directory). In order to create a synthetic spectra, we need some a priori parameters, such as the temperature profile, and the abundances of the molecules included in the model. For my range of study, I use methane, ammonia, phosphine, ethane and acetylene, along with the hydrogen-helium continuum. The opacity is going to define the values of the optical thickness of every molecule (tau). We also calculate the optical thickness independent (kappa) for every pressure layer in our model (since our model is based on an stratified and plane paralell asumption of an atmosphere).

    if (inv(l) == 2) then
        kappa(j1+1:j1+nfreq,k,ii) = (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)
    end if

    tau(j1+1:j1+nfreq,k) = tau(j1+1:j1+nfreq,k) + (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)

The main calculation for the radiance is based on the integration of the blackbody radiation in our atmosphere, multiplied by the transmitance in that specific height, which in a numerical calculation means that we will sum the blackbody of every layer multiplied by the gradient of optical thickness in that layer (see the code to check the previous modifications in the optical thickness before arriving to this step).

    allocate(radiance((repmin-1)*nfreq+1:repmax*nfreq))
    radiance = 0.
    do k = 1, nlevel-1
       radiance = radiance + (hc2*fnu**3) / (exp(fnu/aux_T(k)) - 1.) * (tau(:,k+1) - tau(:,k))
    enddo
    call convol_jwst(size(radiance),radiance,fwhm,fnu,nflux,wave,synthetic)
    deallocate(radiance)
    deallocate(fnu)


As can be seen, the convolution.f90 function must be invoked, since the spectra is generated at a super high resolution due to the data from the databases, but since we want to compare the results with the spectra obtained by an instrument (for example the JWST), we must convolve the data. Due to the JWST has a variable resolution along the wavelength range, we need a file called resolution.txt where we specify the resolution of the instrument, for a wavenumber range BIGGER than the range studied, since the subroutine interpolation (in read_file.f90) needs this condition for the interpolation.
Here you see the convolution and the interpolation routines without the part of the code where we define the variables (thank you, Fortran :unamused: )

    do i=1, nflux
         where (abs(fnu - wave(i)) < 2.5*fwhm(i))
            aux = 2**(-((fnu-wave(i))/(fwhm(i)/2.))**2)
         elsewhere
            aux = 0.
         end where
         synthetic(i) = sum(radiance*aux) / sum(aux)
    enddo
!==========================================================

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

As an introduction to the inversion algorithm (we use the code implemented in Conrath et al. 1998) we have to calculate the jacobian matrixes of the vectors we want to retrieve (for example the temperature and the abundances). The way to calculate several vectors, called x and y, they are calculated by following the next equation:

![\\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}](https://latex.codecogs.com/svg.latex?\Large&space;I_{i}=\sum_{j=1}^{n}%20\frac{\partial%20I_{i}}{\partial%20x_{1,%20j}}%20\Delta%20x_{1,%20j}+\sum_{j=1}^{n}%20\frac{\partial%20I_{i}}{\partial%20x_{2,%20j}}%20\Delta%20x_{2,%20j})


If we develop the equation, we can calculate the jacobian and store it in the kk array, a 3D array in which every slice stores the jacobian of the vectors that want to be inverted. The calculation of the jacobians fot the temperature and the abundances is:

    do l=1, p_mol
     if (inv(l) == 1) then
!==============================================Derivee temperature
    
         do k = 1, nlevel-1
            radiance = (hc2*fnu**3) * (1./(exp(fnu/(aux_T(k)+dt))-1.) - 1./(exp(fnu/(aux_T(k)))-1.)) *(tau((repmin-1)*nfreq+1:repmax*nfreq:10,k+1) - tau((repmin-       1)*nfreq+1:repmax*nfreq:10,k))
            call convol_jwst(size(radiance),radiance,fwhm,fnu,nflux,wave,kk(:,k,1))
         end do

     elseif (inv(l) == 2) then
!==============================================Derivee abondance

         radiance = 0
         do k=1, nlevel-2
            radiance = radiance - (hc2*fnu**3) * tau((repmin-1)*nfreq+1:repmax*nfreq:10,k) *(1./(exp(fnu/aux_T(k))-1.) - 1./(exp(fnu/aux_T(k+1))-1.))
            call convol_jwst(size(radiance),radiance*kappa(::10,k,ii)/mue,fwhm,fnu,nflux,wave,kk(:,k,ii+shift))
         enddo
       ii = ii + 1
     endif
    enddo

## Inversion algorithm

From Conrath et al. 1998, we can finally meet to the final equation needed to retrieve the 'correct' profiles:


![U=\alpha_1SK_1^T(\alpha_1K_1SK_1^T+\alpha_2K_2SK_2^T+E^2)^{-1}](https://latex.codecogs.com/svg.latex?U=\alpha_1SK_1^T(\alpha_1K_1SK_1^T+\alpha_2K_2SK_2^T+E^2)^{-1})

![V=\alpha_2SK_2^T(\alpha_1K_1SK_1^T+\alpha_2K_2SK_2^T+E^2)^{-1}](https://latex.codecogs.com/svg.latex?U=\alpha_2SK_2^T(\alpha_1K_1SK_1^T+\alpha_2K_2SK_2^T+E^2)^{-1})
