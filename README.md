# Radiative transfer and inversion code
## Radiative transfer
The radiative transfer code has its basis on GEISA and HITRAN databases. These databases provide high resolution spectroscopy data for several molecules. We use this
data to know the opacities of the molecules in our model. These files use de .opa denomination, and are not included due to its large size of tens of GB for each molecule (for more info about the format and info in this files ask to the mail in my welcome directory). In order to create a synthetic spectra, we need some a priori parameters, such as the temperature profile, and the abundances of the molecules included in the model. For my range of study, I use methane, ammonia, phosphine, ethane and acetylene, along with the hydrogen-helium continuum. The opacity is going to define the values of the optical thickness of every molecule. We also calculate the optical thickness independent for every pressure layer in our model (since our model is based on an stratified and plane paralell asumption of an atmosphere).


            if (inv(l) == 2) then
                 kappa(j1+1:j1+nfreq,k,ii) = (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)
            end if

            tau(j1+1:j1+nfreq,k) = tau(j1+1:j1+nfreq,k) + (coef(:,1)**(1-dt))*(coef(:,2)**dt)*colonne(k,l)
