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

Where alpha is a value that can change its definition, and will affect to the degrees of freedom of our inversion. In our case we define it as the trace of the error divided by the trace of the array K, all multiplied by a factor f that is fitted after trying various inversions with a value around 3. 
The S array will make that the free vectors will change with a maximum variation of 0.75 scale heights (we can change this by changing the 'c' parameter). The K array is the jacobian, with the sub index identifying it as the jacobian of an specific vector. There will be as many K's as free vectors in our inversion. Other important parameter is the degree of freedom, resulting of the trace of U times K. The amount delta that we have to change our profile can be calculated as U*(spec_obs-spec_syn). The part of the code that does this is the subroutine coeur, in info_content.

    do i=1,nlevel
       do j=1,nlevel
           s(i,j) = exp( -alog(p(i)/p(j))**2 / (2*c**2))
       enddo
    enddo
    trace_error = sum((/(error(i)**2,i=1,nflux)/) / sqrt(real(nflux,kind=4)))
    do ii=1, n_inv  
    	aux(:,:,ii) = matmul(kk(:,:,ii),s) 
    	aux1(:,:,ii) = matmul(aux(:,:,ii),transpose(kk(:,:,ii)))
    	trace_temp(ii) = sum((/(aux1(i,i,ii),i=1,nflux)/))
    	alpha(ii) = factor * (trace_error/trace_temp(ii))
        aux1(:,:,ii) = alpha(ii)*aux1(:,:,ii)
    end do
    aux1x = 0
    do ii =1, n_inv
    	aux1x = aux1x + aux1(:,:,ii)
    end do
    do i=1,nflux
       aux1x(i,i) = aux1x(i,i) + error(i)**2
    enddo
    call matinv(aux1x,nflux,nflux,aux2)
    do ii=1, n_inv
    	aux3(:,:,ii) = matmul(transpose(kk(:,:,ii)),aux2)
    	w(:,:,ii) = matmul(s,aux3(:,:,ii))
    	delta(:,ii) = alpha(ii) * matmul(w(:,:,ii),(spec_obs-spec_syn))
    	A(:,:,ii) = alpha(ii) * matmul(w(:,:,ii),kk(:,:,ii))
        dr(ii) = sum((/ (A(i,i,ii), i=1, size(A, 1)) /))
    end do    
    do ii = 1, n_inv
    	do i=1, nlevel
       	   sigma(i,ii) = alpha(ii) * sqrt(dot_product(w(i,:,ii)**2,error**2))
    	end do
    end do

With this, the only thing we have to do is to add delta to its corresponding pta profile, and then redo the algorithm as long as you want (in my case as long as the chi^2 between the last and the current iteration is lower than the 1%). There are several parts of the code that are not explained, but they are out of the scope of this explanation and you are free to investigate by yourself.

## Example

Lets say that we want to invert, at the same time, temperature and ammonia abundances. First we need the spectra, which is stored in the .spe file. There the first column is the wavenumber, the second the radiance in W cm^-2 sr^-1 / cm^-1 and the third one the error. The first value indicates the length of the file (number of rows). 

      871
      862.06  9.03333e-08  3.00000e-10
      862.34  8.79264e-08  3.00000e-10
      862.63  8.46567e-08  3.00000e-10
      862.92  8.21619e-08  3.00000e-10
      863.21  8.39743e-08  3.00000e-10
      863.50  8.65949e-08  3.00000e-10
      863.79  8.73493e-08  3.00000e-10
      864.07  8.58774e-08  3.00000e-10
                 ...
      
As mentioned, in the .pta we have data of the molecules and temperature. The order is Pressure (bar), Temperature (K), abond: Methane, Methane deuterated, Ammonia, Phosphine, Acetylene and Ethane.

	1.00E+01 339.1 2.10E-03 1.60E-07 2.36E-04 7.00E-07 4.60E-13 2.91E-09 0.00E+00
	9.44E+00 333.4 2.10E-03 1.60E-07 2.36E-04 7.00E-07 4.80E-13 3.03E-09 0.00E+00
	8.91E+00 327.7 2.10E-03 1.60E-07 2.36E-04 7.00E-07 5.01E-13 3.15E-09 0.00E+00
	8.41E+00 322.0 2.10E-03 1.60E-07 2.36E-04 7.00E-07 5.22E-13 3.28E-09 0.00E+00
	7.94E+00 316.6 2.10E-03 1.60E-07 2.36E-04 7.00E-07 5.45E-13 3.42E-09 0.00E+00
	7.50E+00 311.3 2.10E-03 1.60E-07 2.36E-04 7.00E-07 5.68E-13 3.56E-09 0.00E+00
	7.08E+00 305.9 2.10E-03 1.60E-07 2.36E-04 7.00E-07 5.93E-13 3.70E-09 0.00E+00
	6.68E+00 300.6 2.10E-03 1.60E-07 2.36E-04 7.00E-07 6.18E-13 3.85E-09 0.00E+00
	6.31E+00 295.5 2.10E-03 1.60E-07 2.36E-04 7.00E-07 6.45E-13 4.02E-09 0.00E+00
	5.96E+00 290.5 2.10E-03 1.60E-07 2.36E-04 7.00E-07 6.72E-13 4.18E-09 0.00E+00
                                             ...

Finally, we specify which molecules you want to invert. For that, if you want to invert temperature, you must put a 1 in the ch5 row, and then, if you want to invert a specific molecule, write a 2 in its row. You can only invert temperature by putting only a 1 in ch5, or invert only abondances by putting only 2 in t=all the molecules you want to invert. If you just want to create a synthetic spectra based on your .pta profiles, put all to zero. In this case we want to invert temperature (1 in ch5) and ammonia (2 in nh3):

	01.758 ! latitude (degree)
	1.000 ! cos(incidence angle)
	1.080 ! cos(emergence angle)
	0.500 ! fwhm (cm-1)
	ch5 1.080E+00 1
	ch3e 1.754E+03 0
	nh3 1.080E+00 2
	ph3 1.080E+00 0
	c2h3 1.080E+00 0
	c2h6 1.080E+00 0

Finally, after the inversion you will get new .res, .inv, .avg and .sum files. In .res you will get your real and synthetic spectra. In .avg you will get as many columns as free vectors, and they are the kernel of that inversion (you can use it to see where the information comes from). The .sum tells you the number of iterations, the final chi2 and the degrees of freedom of every free vector. Finally, in the .inv file, we will get a first pressure column, followed by a temperature column (these two columns will always appear, independently from inverting the temperature or not). The next columns will be the inversion of your molecules, then a 'NaN' column, and finally the sigma of every inversion. So in our case the first column is the pressure, the second the temperature (in our case inverted), the third the inverted ammonia profile, then the 'NaN' and then the sigma of T and the sigma of NH3.

 	0.10E+02 339.1 0.238E-03 NaN       0.651E-04 0.332E-04
  	0.94E+01 333.4 0.239E-03 NaN       0.831E-04 0.427E-04
  	0.89E+01 327.7 0.240E-03 NaN       0.105E-03 0.546E-04
  	0.84E+01 322.0 0.241E-03 NaN       0.133E-03 0.693E-04
  	0.79E+01 316.5 0.242E-03 NaN       0.167E-03 0.874E-04
  	0.75E+01 311.2 0.244E-03 NaN       0.207E-03 0.109E-03
  	0.71E+01 305.8 0.246E-03 NaN       0.256E-03 0.136E-03
  	0.67E+01 300.5 0.248E-03 NaN       0.316E-03 0.169E-03
  	0.63E+01 295.4 0.251E-03 NaN       0.386E-03 0.208E-03
                               ...


NOVEMBER 2022 UPDATE (jupiter_inv3)

This new version has spherical geometry, whic allows the retrieval of the spectra at high latitudes (but not very close to the limb). Other several changes were done, the model has two cloud layers that can be cativated in .fuel with cloud = 1 activating the retrieval of tropospheric opacity of the clouds and cloud = 2 activating the retrieval of stratospheric opacity of hazes.

The .avg file gives us the A matrix to know where the information comes from, the .sum file summarizes some information on the retrieval.

We can now invert also the temperature and methane profile (not recommended) by putting a 3 instead of a 1 or 2 in the .fuel ch4 row.

Some masking is applied to the jacobian to remove some non physical values close to the border of the height grid.

In info_content.f90, in coeur routine we have two functions (f and f2) to activate a threshold to let some profiles to change or not. For the retrieval of the hydrocarbons, f must be activated so 'delta(:,1) = delta(:,1)*f(:)' must be activated.


With this you know the basis of the code. 
If used, please cite some of the original authors of the code (Ex: Fouchet et al. 2000)

Last mod: 15th May 2022 Pablo Rodriguez Ovalle
