# jupiter_jwst

radiative transfer for JWST

Created inv_pos array. It tells you how many molecules are being inverted (length of the array) and which one are they (position)

Ex: knowing that CH4:1  ;  CH3D:2  ;  NH3:3  ;  PH3:4  ;  C2H2:5  ;  C2H6:6,
the inv_pos [2,4,5,6] means that we are inverting 4 molecules simultaneously: CH3D, PH3, C2H2 and C2H6

Developed the algorithm to retrieve simultaneously several molecular profiles. For that, the Jacobian array KK is now a 3D array 
with information of one specific molecule for every of its slices. Later, the coeur routine now separates and calculates
different alphas (alpha,beta,gamma...) and traces as a function of the number of molecules that you want to invert.

Write_kernel now prints the kernel array for every of the profiles inverted (temperature or abundances).

Currently the write_inv routine writes the pressure column along with the temperature (always shown) inverted profiles columns and
then the sigmas for the molecules inverted, separated by a 'NaN' column (can be removed by changing inv_pos for inv_pos_red in write_inv routine) .

HOW IT WORKS: The code now is much more simple to use. You will need a .pta file with your a priori profiles, and a .spe file with
your real spectra. In .fuel you will specify which molecules you want to invert marking them with a 2. CAUTION: If you are inverting
a molecule in a spectral region where there is NO information of it, you will get an error, so please try to add only the molecules
which will have an impact on the region you are studying.

last mod: 13th May 2022 by Pablo Rodriguez Ovalle
