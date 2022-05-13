# Radiative transfer and inversion code
## Radiative transfer
The radiative transfer code has its basis on GEISA and HITRAN databases. These databases provide high resolution spectroscopy data for several molecules. We use this
data to know the opacities of the molecules in our model. These files use de .opa denomination, and are not included due to its large size of tens of GB for each molecule (for more info about the format and info in this files ask to the mail in my welcome directory). In order to create a synthetic spectra, we need some a priori parameters, such as the temperature profile, and the abundances of the molecules included in the model. For my range of study, I use $CH_4$
![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)
