      integer nlevel
      parameter(nlevel=361)

      integer nfreq, nrep
      parameter(nfreq=10000, nrep=17)

      real gnu00, dgnu, dgnu_H2
      parameter (gnu00=1, dgnu=1e-2, dgnu_H2=2)

      integer nmol
      parameter (nmol=4)

      integer nraie_max
      parameter (nraie_max=70000)

      integer nt
      parameter(nt=1)
      real t_tab(nt,nlevel)

      real dens, H2, He, CH4, Rsat, J2sat, J3sat, J4sat, omega, ellip
      real MassSat
      parameter (dens=2.35e-3, H2=0.86, He=0.1355, CH4=4.5e-3)
      parameter (Rsat=71492e3,J2sat=16.45,J3sat=0.,J4sat=-1000.0)
      parameter (omega=0.4136,ellip=0.0648744,MassSat=1898.2)

      real avo, boltz, cvel, planck, rg, Tref
      parameter (avo=6.02214e23,boltz=1.38066e-23,cvel=2.99792e+10)
      parameter (planck=6.62608e-34,rg=8.31451,Tref=296.)

