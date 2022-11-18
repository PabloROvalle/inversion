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
then the sigmas for the molecules inverted, separated by a 'NaN' column.

HOW IT WORKS: The code now is much more simple to use. You will need a .pta file with your a priori profiles, and a .spe file with
your real spectra. In .fuel you will specify which molecules you want to invert marking them with a 2. CAUTION: If you are inverting
a molecule in a spectral region where there is NO information of it, you will get an error, so please try to add only the molecules
which will have an impact on the region you are studying.

*NOVEMBER 2022 UPDATE (jupiter_inv2)*

This new version has spherical geometry, whic allows the retrieval of the spectra at high latitudes (but not very close to the limb). Other several changes were done, the model has two cloud layers that can be cativated in .fuel with cloud = 1 activating the retrieval of tropospheric opacity of the clouds and cloud = 2 activating the retrieval of stratospheric opacity of hazes.

The .avg file gives us the A matrix to know where the information comes from, the .sum file summarizes some information on the retrieval.

We can now invert also the temperature and methane profile (not recommended) by putting a 3 instead of a 1 or 2 in the .fuel ch4 row.

Some masking is applied to the jacobian to remove some non physical values close to the border of the height grid.

In info_content.f90, in coeur routine we have two functions (f and f2) to activate a threshold to let some profiles to change or not. For the retrieval of the hydrocarbons, f must be activated so 'delta(:,1) = delta(:,1)*f(:)' must be activated.


With this you know the basis of the code. 
If used, please cite some of the original authors of the code (Ex: Fouchet et al. 2000)

Last mod: 22th November 2022 Pablo Rodriguez Ovalle

