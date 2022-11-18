module gravity

  implicit none
  
contains


!*******************************************************************************
!*******************************************************************************
!------------------------------------------------------------------------------
!
!                           SUBROUTINE NEWGRAV
!
!      Sets the gravitational acceleration based on selected planet and 
!      latitude. Gravity is calculated normal to the surface, and corrected 
!      for rotational effects. The input latitude is assumed planetographic.
!
!      The radial vector component of the gravitational acceleration is
!      given by Lindal et al., 1985, Astr. J., 90 (6), 1136-1146, as the
!      expression
!                      ___inf
!          G M  (    \                    (R)**2n                  )
!      g = ---- ( 1 - \     (2n + 1) J_2n (-)     P_2n (sin(latc)) )
!          r**2 (     /__n=1              (r)                      )
!
!                       2
!                    -  - omega**2 r (1 - P_2(sin(latc)))
!                       3
!
!      where r = radial distance, G = gravitational constant, M = planet
!      mass, J_2n is the 2nth zonal harmonic coefficient, P_2n is the
!      Legendre polynomial of order 2n, and latc is the planetocentric
!      latitude. Note that the notation is using spherical coords based
!      on the planet centre. The last part of the expression is simply
!      the centrifugal component and is more simply written as
!
!                           omega**2 r cos(latc)**2.
!
!      The latitudinal vector component is given by the expression
!
!                 ___inf
!          G M  ( \            (R)**2n  d P_2n (sin(latc)) )
!      g = ---- (  \      J_2n (-)      ------------------ )
!          r**2 (  /__n=1      (r)            d latc       )
!
!                      1            d P_2 (sin(latc))
!                    + - omega**2 r -----------------
!                      3                 d latc
!
!      where the notation is as above. The derivatives of the Legendre
!      polynomials can be found through the recurrence relation
!
!                        d P_2n (z)
!             (z**2 - 1) ---------- = 2 n (z P_2n (z) - P_(2n-1) (z)).
!                           d z
!
!      Using this relation, it can be seen that the last expression of
!      the latitudinal component represents centrifugal acceleration and
!      is more clearly written as
!
!                           omega**2 r cos(latc) sin(latc).
!
!      These two gravitational vector components are normal to each other
!      and must be summed to give the overall gravity acting normal to the
!      planetary surface.
!
!      Original version    A.Weir
!      Adapted for general use 29/4/96  Pat Irwin
!-------------------------------------------------------------------------------

  subroutine newgrav (lat_in, h, g, radius) 

    use declaration, only : RJup, J2Jup, J4Jup, J6Jup, omega, ellip, MassJup,&
         Grav, pi
    implicit none

    ! Latitude (degrees), Height above reference surface (km)
    real, intent(in) :: lat_in, h
    !gravity acc. (m/s2)
    real, intent(out):: g, radius
    
    integer :: I
    real :: lat, latc, slatc, r, clatc, gtheta, Rr
    real :: xgm, xradius, xellip, xomega
    real, dimension(3) :: xcoeff
    real, dimension(6) :: pol

    !=============================================
    xgm = MassJup*Grav*1e24*1e6
    xcoeff(1) = J2Jup/1e3
    xcoeff(2) = J4Jup/1e6
    xcoeff(3) = J6Jup/1e6
    xradius = Rjup*1e2 ! en cm
    xellip = 1./(1.-ellip)
    xomega = 2*pi/(omega*24.0*3600.0)
    lat = 2*pi * lat_in/360.
    latc = thetagc(lat, xellip)
    slatc = sin(latc)
    clatc = cos(latc)

    ! Rr is the ratio of radius at equator to radius at current latitude
    Rr = sqrt(clatc**2 + (xellip**2 * slatc**2))
    r = (xradius+h*1e5)/Rr ! 1e5 pour km-> cm
    radius = (xradius/Rr)*1e-5

    do I = 1, 6
       pol(I) = legpol(I,slatc)
    enddo

!-------------------------------------------------------------------------------
!
!      Evaluate radial contribution from summation for first three terms,
!      then subtract centrifugal effect.
!
!-------------------------------------------------------------------------------

    g = 1.
    do I = 1, 3
       g = g - ((2*I+1) * Rr**(2 * I) * xcoeff(I) * pol(2*I))
    enddo
    g = (g * xgm/r**2) - (r * xomega**2 * clatc**2)

!-------------------------------------------------------------------------------
!
!      Evaluate latitudinal contribution for first three terms, then add
!      centrifugal effects.
!
!-------------------------------------------------------------------------------

    gtheta = 0.
    do I = 1, 3
       gtheta = gtheta - (4 * I**2 * Rr**(2 * I) * xcoeff(I) * (pol(2*I-1) - slatc * pol(2*I))/clatc)
    enddo
    gtheta = (gtheta * xgm/r**2) + (r * xomega**2 *  clatc * slatc)
    
!-------------------------------------------------------------------------------
!
!      Combine the two components and write the result.
!
!-------------------------------------------------------------------------------

    g = sqrt(g**2 + gtheta**2)*0.01 ! unité de cm/s2 -> m/s2
!    print*, lat_in, h, g 

  contains

!********************************************************************************
!********************************************************************************
!-------------------------------------------------------------------------------
!
!                           FUNCTION LEGPOL
!
!      Calculates zero order Legendre polynomials based on the recurrence
!      relation 
!               (n) P   (z) = (2n-1) z P    (z) - (n-1) P     (z).
!                    (n)                (n-1)            (n-2)
!
!-------------------------------------------------------------------------------

    function legpol (n, z)
      
      implicit none
      integer :: I
      integer, intent(in) :: n
      real, intent(in) :: z
      real :: legpol
      real, dimension(3) :: pol

      pol(1) = 1.
      pol(2) = z

      if ((n==0).or.(n==1)) then
         legpol = pol(n+1)
      else
         I = 1
10       I = I + 1
         
         if (I>50) then
            write (*,*) '* * * * * * * * * * * * *'
            write (*,*) '*                       *'
            write (*,*) '*   Looping in LEGPOL   *'
            write (*,*) '*   Program terminated  *'
            write (*,*) '*                       *'
            write (*,*) '* * * * * * * * * * * * *'
         endif

         pol(3) = (((2 * I - 1) * pol(2) * z)  - ((I - 1) * pol(1)))/float(I)
         if (I==n) then
            legpol = pol(3)
            return
         else
            pol(1) = pol(2)
            pol(2) = pol(3)
            goto 10
         endif
      endif

    end function legpol


!*******************************************************************************
!*******************************************************************************

!-------------------------------------------------------------------------------
!
!                           FUNCTION THETAGC
!
!      Converts planetographic latitude to planetocentric. Angles must be
!      supplied in radians.
!
!-------------------------------------------------------------------------------

    function thetagc (lat, e)

      implicit none
      real, intent(in) :: e, lat
      real :: thetagc
      
      thetagc = atan(tan(lat)/e**2)

    end function thetagc


!*******************************************************************************
!*******************************************************************************

  end subroutine newgrav

end module gravity
