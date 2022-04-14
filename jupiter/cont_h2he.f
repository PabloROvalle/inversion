      Subroutine CONT_H2He(p,T,H2,He,dens,k,grav,ftab,nin,tabH2)


c-----VERSION  21 - 08 - 1996----------------------------------------------

c-----In this subroutine the continuum of H2 and He is calculated at a given
c-----temperature and pressure, and for a given set of wavenumbers ,  ftab.
c-----The subroutines that follow are written by  Jacek Borysow and--------
c-----Lothar Frommhold. This part of the program is an interface between---
c-----the original program called  PROF_BOR.F  and the radiative transfer--
c-----code  NIMS_RT.F -----------------------------------------------------


c--------------------------------------------------------------------------
c---BEGIN--PARAMETER DECLARATION-------------------------------------------

      implicit none
      
      Integer  i ,              ! counter
     +         k ,              ! the level in the effective atmosphere
     +                          ! being considered
     +         nin

      Real  zn, aux, grav, H2, He, CH4, dens, avo, dgnu_H2

      real p(361), T(361)

      Real*8  abCH4( 801) ,        ! the opacity as a function of wavenumber
     +                            ! at level  k  due to CH4
     +        abHel( 801) ,        ! the opacity as a function of wavenumber
     +                            ! at level  k  due to He
     +        abHyd( 801)          ! the opacity as a function of wavenumber
     +                            ! at level  k  due to H2


      real ftab(2000), tabH2(2000,361)



c-----END--PARAMETER DECLARATION-------------------------------------------
c--------------------------------------------------------------------------

      avo=6.02214e23
      dgnu_H2=2
	CH4=1-H2-He
      if (k.eq.321) then
         zn = 273.16 * p(k) / T(k) / 1.01325
         aux = p(321)*1e5 / (grav*dens) * avo/2.6868e19 / 1e4
      else
         zn = 273.16 * 2* sqrt(p(k)*p(k+1)) / (T(k)+T(k+1)) / 1.01325
         aux = (p(k)-p(k+1))*1e5 / (grav*dens) * avo/2.6868e19 / 1e4
      endif
      Call H2CIA(abHyd,abHel,abCH4,ftab(1),ftab(nin),dgnu_H2,T(k))
      Do i = 1 , nin
         tabH2(i,k) = aux*H2 * zn *
     +      ( H2*abHyd(i) + He*abHel(i) + CH4*abCH4(i)) 
      Enddo
      
      Return
      End

c------------------------------------------------------------------------------
      Subroutine H2CIA(abhyd,abhel,abch4,gnu0,gnumax,dgnu,temp)

      Real*8 ABHYD( 801),ABHEL( 801),ABCH4( 801)

      Call addemh2(abhyd,gnu0,gnumax,dgnu,temp)
      Call addemhe(abhel,gnu0,gnumax,dgnu,temp)
      Call addemch4(abch4,gnu0,gnumax,dgnu,temp)

      Return
      End

c--------------------------------------------------------------------------

      Subroutine ADDEMh2(alfatot,gnu0,gnumax,dgnu,tp)

C     THIS PROGRAM GENERATES THE H2-H2 TRANSLATIONAL/ROTATIONL
C     CIA SPECTRA. IT IS BASED ON QUANTUM LINE SHAPE COMPUTATIONS AND
C     THE AB INITIO DIPOLE DATA BE W. MEYER. DIMER FINE STRUCTURES ARE
C     SUPPRESSED. THE PROGRAM WAS WRITTEN BY ALEKSANDRA BORYSOW AND
C     LOTHAR FROMMHOLD. THIS IS THE NOVEMBER 1985 VERSION

C     H2-H2 COMPUTATIONS REFERENCE: MEYER, FROMMHOLD AND BIRNBAUM,
C     TO BE PUBLISHED IN PHYS.REV.A IN 1985;
C     THE H2-H2 MODELING USED HERE IS BEING PUBLISHED: J.BORYSOW,
C     L.TRAFTON, L.FROMMHOLD, G.BIRNBAUM, AP.J. (1985)

C     TAPE3 IS OUTPUT: HEADER PLUS ABSORPTION COEFF. ALPHA(NU)

       Implicit double precision (a-h,o-z)
       Dimension alfatot( 801)
       Real gnu0,gnumax,dgnu,tp

      Common /RESULT1/ NF
      Common /result2/ FREQ( 801),ABSCOEF( 801)
      Common /H2PART1/ AUXIL(5)
      Common /h2part2/ IAUX,NORMAL

       Y(X,A,B,C)=A*DEXP((C*X+B)*X)

       normal=0
       temp=tp 

       fnumin=gnu0
       fnumax=gnumax
       dnu=dgnu

       NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
       IF (NF.GT. 801) NF= 801
       FNUMAX=FNUMIN+DFLOAT(NF-1)*DNU

       Call PARTSUM (TEMP)


C     THE H2-H2 SPECTRA
C     =================


       X=DLOG(TEMP)

       DO I = 1,NF

        FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU
        alfatot(i)=0.d0
        ABSCOEF(I)=0.

       Enddo


C     THE LAMBDA1,LAMBDA2,LAMBDA,L = 2023 AND 0223 COMPONENTS:
C     (QUADRUPOLE PLUS OVERLAP INDUCED COMPONENT)
C	sum of 2023+0223:

       S=Y(X,2.881d-61,-1.1005d0,0.1310d0)
       E=Y(X,7.3485d0,-1.3874d0,0.1660d0)
       T1=Y(X,7.883d-13,-.3652d0,-.0271d0)
       T2=Y(X,3.803d-13,-.4048d0,-.0091d0)
       T3=Y(X,1.0922d-12,-.4810d0,-.0127d0)
       T4=Y(X,5.174d-12,-.9841d0,0.0483d0)

       Call  ADDSPEC1 (S,E,T1,T2,T3,T4,TEMP,0,1,0,2,2,3,0,0,1.)

       Do i = 1,nf

	alfatot(i) = alfatot(i) + abscoef(i)

       Enddo 


C     PARAMETERS FOR 4045 AND 0445 (PURE HEXADECAPOLE) COMPONENTS

       S=Y(X,2.404d-64,-1.4161d0,0.1847d0)
       E=Y(X,-.8033d0,-.4474d0,-.0235d0)
       T1=Y(X,3.873d-13,-.4226d0,-.0183d0)
       T2=Y(X,2.743d-13,-.3566d0,-.0140d0)
       T3=Y(X,4.171d-13,-.5223d0,0.0097d0)
       T4=Y(X,2.2725d-12,-1.1056d0,0.0139d0)

       Call ADDSPEC1(S,E,T1,T2,T3,T4,TEMP,0,1,4,0,4,5,0,0,1.)

       Do i = 1,nf

        alfatot(i) = alfatot(i) + abscoef(i)

       Enddo

    
C     PARAMETERS FOR 0221 AND 2021 (PURE OVERLAP) COMPONENTS

       S=Y(X,6.393d-63,-1.5964d0,0.2359d0)
       E=Y(X,21.414d0,-1.2511d0,0.1178d0)
       T1=Y(X,1.876d-13,-.4615d0,-.0012d0)
       T2=Y(X,4.839d-13,-.5158d0,0.0075d0)
       T3=Y(X,4.550d-13,-.5507d0,0.0095d0)
       T4=Y(X,2.045d-12,-.5266d0,-.0240d0)

       Call  ADDSPEC1(S,E,T1,T2,T3,T4,TEMP,0,1,0,2,2,1,0,0,1.)

       Do i = 1,nf

        alfatot(i) = alfatot(i) + abscoef(i)

       Enddo


C     PARAMETERS FOR 2233 QUADRUPOLE INDUCED DOUBLE TRANSITIONS

       S=Y(X,5.965d-63,-1.0394d0,0.1184d0)
       E=Y(X,6.674d0    ,-.9459d0,0.1124d0)
       T1=Y(X,4.764d-13,-.1725d0,-.0450d0)
       T2=Y(X,4.016d-13,-.3802d0,-.0134d0)
       T3=Y(X,1.0752d-12,-.4617d0,-.0085d0)
       T4=Y(X,1.1405d-11,-1.2991d0,0.0729d0)

       Call ADDSPEC1(S,E,T1,T2,T3,T4,TEMP,0,1,2,2,3,3,0,0,1.)
       Do i = 1,nf

        alfatot(i) = alfatot(i) + abscoef(i)

       Enddo

c  160  Format( f7.1, e12.4)

c  140  Format (/,' Total ABSORPTION COEFFICIENT ALPHA(fnu)',/,
c     +  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
c     +  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)

c  150  Format(6E12.4)


c   30  Format( 34H1ABSORPTION SPECTRA OF HYDROGEN AT,F8.1,  2H K,/1X,43(
c     11H=),/, 11H MIN.FREQ.=,F8.1,  5H CM-1,10X, 10HMAX.FREQ.=,F8.1,  5H
c    2 CM-1,10X, 15HFREQ.INCREMENT=,F8.2,  5H CM-1,5X,  2HIN,I5,  6H STE
c    3PS,//)
c   40  Format(51H1ABSORPTION SPECTRUM OF HYDROGEN-HELIUM MIXTURES AT,F8
c     +.1,  2H K,/1X,59(1H=)/3F15.1,I10//)

       Return	
       End



c--------------------------------------------------------------------------

      Subroutine ADDSPEC1(G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,MP,LIKE,
     +                    lambda1,LAMBDA2,LAMBDA,LVALUE,NVIB1,
     +                    NVIB2,FACTOR)

C     THIS PROGRAM GENERATES A LISTING OF THE CIA TR ALFA(OMEGA)
C     IF BOTH LAMBDA1 AND LAMBDA2 ARE NEGATIVE: SINGLE TRANSITIONS;
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.
C     MP=1 GIVES LISTINGS OF INTERMEDIATE RESULTS.
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C

      Implicit double precision (a-h,o-z)

      Common /H2PART1/ Q,WH2(2),B0,D0
      Common /h2part2/ JRANGE1,NORMAL
      Common /RESULT1/ NF
      Common /result2/ FREQ( 801),ABSCOEF( 801)

      DATA CLOSCHM,BOLTZWN/2.68675484E19,.6950304/
      DATA HBAR,PI,CLIGHT/1.054588757d-27,3.1415926535898,2.9979250E10/

      EH2(N,I)=4395.34*(DFLOAT(N)+0.5)-117.905*(DFLOAT(N)+0.5)**2
     +  +(60.809-
     +  2.993*(DFLOAT(N)+0.5)+.025*(DFLOAT(N)+.5)**2)*DFLOAT(I)-
     + (.04648-.00134*(DFLOAT(N)+.5))*DFLOAT(I*I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*
     +  DEXP(-1.4387859/T*EH2(0,J*(J+1)))

      Do i = 1,nf

       abscoef(i) = 0.d0

      Enddo

      TWOPIC=2.*PI*CLIGHT

c      IF (MP.NE.1) then 	! suppressed here because mp is a constant
c      IF (LIKE .NE. 1) LIKE=0  ! suppressed here because like is a constant

      CALIB = TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB = CALIB/DFLOAT(1+LIKE)
      BETA = 1./(BOLTZWN*TEMP)
      LIST = NF



C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************


      IF ((LAMBDA1.LT.0).AND.(LAMBDA2.LT.0)) GO TO 60
      JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 50 I1=1,JRANGE1
         J1=I1-1
      DO 50 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR(J1,LAMBDA1,JP1)
         IF (CG1S) 50,50,10
   10    P1=PH2(J1,TEMP)/Q
         IF (P1.LT.0.001) GO TO 50
         OMEGA1=EH2(NVIB1,JP1*IP1)-EH2(0,J1*I1)
         DO 40 I2=1,JRANGE1
            J2=I2-1
         DO 40 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR(J2,LAMBDA2,JP2)
            IF (CG2S) 40,40,20
   20       P2=PH2(J2,TEMP)/Q
            IF (P2.LT.0.001) GO TO 40
            OMEGA2=EH2(NVIB2,JP2*IP2)-EH2(0,J2*I2)
            FAC=CALIB*P1*P2*CG1S*CG2S
            DO 30 I=1,LIST
               FRQ=FREQ(I)-OMEGA1-OMEGA2
               WKI=FREQ(I)*(1.-DEXP(-BETA*FREQ(I)))
               WKF=WKI*FAC
               XBG=G0*BGAMA(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
               ABSCOEF(I)=ABSCOEF(I)+XBG*WKF
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      GO TO 100
   60 JPLUSL=JRANGE1+LAMBDA
      DO 90 I=1,JRANGE1
         J=I-1
      DO 90 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR(J,LAMBDA,JP)
         IF (CGS) 90,90,70
   70    P=PH2(J,TEMP)/Q
         IF (P.LT.0.001) GO TO 90
         OMEGA1=EH2(NVIB1,JP*IP)-EH2(0,J*I)
         FAC=CALIB*P*CGS
         DO 80 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))
            WKF=WKI*FAC
            XBG=G0*BGAMA(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   80    CONTINUE
   90 CONTINUE
  100 CONTINUE
      RETURN
C
c  110 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 20H COM
c     1PONENT INCLUDED.,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X/ 
c     25HG(0)=,E12.3/)
c  140 FORMAT (/,' ABSORPTION COEFFICIENT ALPHA(fnu)',/,
c     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
c     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
c  150 FORMAT (  6E12.4)
C
      END


c--------------------------------------------------------------------------



      FUNCTION CLEBSQR (L,LAMBDA,LP)
C
C     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
C     FOR INTEGER ARGUMENTS ONLY
C     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
C
	implicit double precision (a-h,o-z)
      FC=DFLOAT(2*LP+1)
      GO TO 10
C

      ENTRY THREEJ2

C
C     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP
C                                                 0    0    0
C     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT
C     (LIMITATION TO INTEGER ARGUMENTS ONLY)
C
C     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE
C     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT
C     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS.
C
      FC=1.
   10 CLEBSQR=0.
      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
     1) RETURN
      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN
      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN
      F=1./DFLOAT(L+LP+1-LAMBDA)
      IF (LAMBDA.EQ.0) GO TO 30
      I1=(L+LP+LAMBDA)/2
      I0=(L+LP-LAMBDA)/2+1
      DO 20 I=I0,I1
   20 F=F*DFLOAT(I)/DFLOAT(2*(2*I+1))
   30 P=FC*F*FCTL(LAMBDA+L-LP)*FCTL(LAMBDA+LP-L)
      CLEBSQR=P/(FCTL((LAMBDA+L-LP)/2)*FCTL((LAMBDA+LP-L)/2))**2
      RETURN
C
      END


c--------------------------------------------------------------------------


      FUNCTION FCTL (N)
	implicit double precision (a-h,o-z)
      P(Z)=((((-2.294720936d-4)/Z-(2.681327160d-3))/Z+(3.472222222d-3))/
     1Z+(8.333333333d-2))/Z+1.
      FCTL=1.
      IF (N.LE.1) RETURN
      IF (N.GT.15) GO TO 20
      J=1
      DO 10 I=2,N
   10 J=J*I
      FCTL=DFLOAT(J)
      RETURN
   20 Z=DFLOAT(N+1)
      FCTL=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
      RETURN
      END

c--------------------------------------------------------------------------

      FUNCTION BGAMA(FNU,T1,T2,EPS,T3,T4,TEMP)
C
C     EQUATION 13, SO-CALLED EBC MODEL, OF BORYSOW,TRAFTON,FROMMHOLD,
C     AND BIRNBAUM, ASTROPHYS.J., TO BE PUBLISHED (1985)
C     NOTE THAT BGAMA REDUCES TO THE BC PROFILE FOR EPS=0.
C
	implicit double precision (a-h,o-z)
      REAL K0
      DATA PI,CLIGHT/3.1415926535898,2.99792458E10/
      DATA HBAR,BOLTZ/1.0545887d-27,1.380662d-16/
      P1(X)=((((((.0045813*X+.0360768)*X+.2659732)*X+1.2067492)*X+3.0899
     1424)*X+3.5156229)*X+1.)
      P2(X)=((((((.00000740*X+.00010750)*X+.00262698)*X+.03488590)*X+.23
     1069756)*X+.42278420)*X-.57721566)
      P3(X)=((((((.00032411*X+.00301532)*X+.02658733)*X+.15084934)*X+.51
     1498869)*X+.87890594)*X+.5)
      P4(X)=((((((-.00004686*X-.00110404)*X-.01919402)*X-.18156897)*X-.6
     17278579)*X+.15443144)*X+1.)
      P5(X)=((((((.00053208*X-.00251540)*X+.00587872)*X-.01062446)*X+.02
     1189568)*X-.07832358)*X+1.25331414)
      P6(X)=((((((-.00068245*X+.00325614)*X-.00780353)*X+.01504268)*X-.0
     13655620)*X+.23498619)*X+1.25331414)
C
      OMEGA=2.*PI*CLIGHT*FNU
      T0=HBAR/(2.*BOLTZ*TEMP)
      Z=SQRT((1.+(OMEGA*T1)**2)*(T2*T2+T0*T0))/T1
      IF (Z-2.) 10,10,20
   10 XK1=Z*Z*DLOG(Z/2.)*P3((Z/3.75)**2)+P4((Z/2.)**2)
      GO TO 30
   20 XK1=SQRT(Z)*DEXP(-Z)*P6(2./Z)
   30 IF (EPS.EQ.0.) GO TO 70
      ZP=SQRT((1.+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4
      IF (ZP-2.) 40,40,50
   40 K0=-DLOG(ZP/2.)*P1((ZP/3.75)**2)+P2((ZP/2.)**2)
      GO TO 60
   50 K0=DEXP(-ZP)*P5(2./ZP)/SQRT(ZP)
   60 BGAMA=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*
     2 XK1/(1.+(T1*OMEGA)**2)+EPS*(T3/
     1PI)*DEXP(T3/T4+T0*OMEGA)*K0)/(1.+EPS)
      RETURN
   70 BGAMA=(T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/(1.+(OMEGA*T1)**2)
      RETURN
C
      END


c--------------------------------------------------------------------------


      SUBROUTINE PARTSUM (TEMP)
C
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).
C
	implicit double precision (a-h,o-z)
      COMMON /H2PART1/ Q,WH2(2),B0,D0
      common /h2part2/ JRANGE1,NORMAL
C
C     NORMAL=0 INITIATES EQUILIBRIUM HYDROGEN;
C     NORMAL=1 ADJUSTS FOR ORTHO/PARA HYDROGEN RATIO OF 3:1
C
      DATA B0,D0,WH2(1),WH2(2)/59.3392,0.04599,1.,3./
      EH2(N,I)=4395.34*(DFLOAT(N)+.5)-117.905*(DFLOAT(N)+.5)**2
     2  +(60.809-2.993*(DFLOAT(N)+.5)+.025*(DFLOAT(N)+.5)**2)*
     2  DFLOAT(I)- (.04648-.00134*(DFLOAT(N)+.5))*DFLOAT(I*I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859*
     1  EH2(0,J*(J+1)
     1)/T)
C
C     Q,B0,D0,WH2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR H2
C
      Q=0.
      J=0
   10 DQ=PH2(J,TEMP)
      Q=Q+DQ
      J=J+1
      IF (DQ.GT.Q/900.) GO TO 10
      JRANGE1=J
      IF (NORMAL) 20,20,30
   20       RETURN
   30 J=-1
      S=0.
      SEV=0.
   40 J=J+1
      DS=PH2(J,TEMP)
      S=S+DS
      SEV=SEV+DS
      J=J+1
      S=S+PH2(J,TEMP)
      IF (J.GT.JRANGE1) GO TO 50
      GO TO 40
   50 SODD=S-SEV
      WH2(2)=WH2(2)*3.*SEV/SODD
      Q=4.*SEV
      RETURN
C
   60 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,10X,
     1 7HJ MAX =,I3,10X, 26HTHERMAL EQUILIBRIUM H2; T=,F10.1,  2H K,/)
   70 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,E12.4,4X,
     16HJ MAX=,I3,4X, 25HN O R M A L  HYDROGEN, T=,F6.1,  2H K,4X,  6HW1
     2,W2=,2F6.2/)
C
      END


c--------------------------------------------------------------------------


      SUBROUTINE ADDEMHE(alfatot,gnu0,gnumax,dgnu,tp)
C
C     THIS PROGRAM GENERATES THE H2-He ROTATIONAL/TRANSLATIONAL 
C     CIA SPECTRA. IT IS BASED ON QUANTUM LINE SHAPE COMPUTATIONS AND
C     THE AB INITIO DIPOLE DATA BY W. MEYER. 
C
C     File out.h2he IS OUTPUT: HEADER PLUS ABSORPTION COEFF. ALPHA(NU)
C
	implicit double precision (a-h,o-z)
	dimension alfatot( 801)
	real gnu0,gnumax,dgnu,tp
      COMMON /RESULT1/ NF
      common /result2/ FREQ( 801),ABSCOEF( 801)
      Y(X,A,B,C,D)=A*DEXP(B*X + C*X*X + D*X*X*X)
	fnumin=gnu0
	fnumax=gnumax
	dnu=dgnu
	temp=tp
c
C
      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
      IF (NF.GT. 801) NF= 801
      FNUMAX=FNUMIN+DFLOAT(NF-1)*DNU
	if(temp.le.40.d0 .or. temp.gt.3000.d0) stop 4444
      CALL PARTSUM (TEMP)
C
C     THE H2-He SPECTRA: 40K -- 3000K
C     =================
C
      X=DLOG(TEMP)
      DO 10 I=1,NF
      FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU
	alfatot(i)=0.d0
   10 ABSCOEF(I)=0.d0
C
C     THE LAMBDA,L = 23 COMPONENT:
C     (QUADRUPOLE PLUS OVERLAP INDUCED COMPONENT)
      XM0=Y(X,.262d-62, -0.71411d0, 0.14587d0, -0.00352d0 )
      XM1=Y(X, 0.269d-49, -0.98315d0, 0.21988d0, -0.00729d0)
      XM2=Y(X, 0.406d-35, -2.25664d0, 0.50098d0, -0.01925d0)
	call BCparam(temp, XM0,XM1,XM2, t1, t2)
	ifun = 1
      CALL  ADDSPEC2 (XM0,T1,T2,TEMP,0,0,-1,-1,2,3,0,0,ifun)
	do 111 i=1, nf
111	alfatot(i) = alfatot(i) + abscoef(i)

C     PARAMETERS FOR 21 (PURE OVERLAP) COMPONENT
      XM0=Y(X, 0.424d-64, -0.37075d0, 0.17473d0, -0.00489d0)
      XM1=Y(X, 0.174d-49, -1.89232d0, 0.44399d0, -0.02029d0)
      XM2=Y(X, 0.874d-35, -3.27717d0, 0.69166d0, -0.02865d0)
	call K0param(temp, XM0, XM1, XM2, t1, t2)
	ifun = 0
      CALL  ADDSPEC2(XM0,T1,T2,TEMP,0,0,-1,-1,2,1,0,0,ifun)
	do 311 i=1, nf
311	alfatot(i) = alfatot(i) + abscoef(i)
C
C     PARAMETERS FOR 01 (PURE OVERLAP) COMPONENT
      XM0=Y(X, 0.223d-61, -1.89198d0, 0.45505d0, -0.02238d0)
      XM1=Y(X, 0.625d-48, -1.96486d0, 0.47059d0, -0.02402d0)
      XM2=Y(X, 0.316d-33, -3.40400d0, 0.72793d0, -0.03277d0)
 	S=XM0
	call K0param(temp, XM0, XM1, XM2, t1, t2)
	ifun = 0
      CALL  ADDSPEC2(XM0, T1,T2,TEMP,0,0,-1,-1,0,1,0,0, ifun)
	do 319 i=1, nf
319	alfatot(i) = alfatot(i) + abscoef(i)

160	format( f7.1, e12.4)

  140 FORMAT (/,' Total ABSORPTION COEFFICIENT ALPHA(fnu)',/,
     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
  150 FORMAT (  6e12.4)

C
   30 FORMAT ( 34H1ABSORPTION SPECTRA OF HYDROGEN AT,F8.1,  2H K,/1X,43(
     11H=),/, 11H MIN.FREQ.=,F8.1,  5H CM-1,10X, 10HMAX.FREQ.=,F8.1,  5H
     2 CM-1,10X, 15HFREQ.INCREMENT=,F8.2,  5H CM-1,5X,  2HIN,I5,  6H STE
     3PS,//)
   40 FORMAT ( 51H1ABSORPTION SPECTRUM OF HYDROGEN-HELIUM MIXTURES AT,F8
     1.1,  2H K,/1X,59(1H=)/3F15.1,I10//)

	return
      END

c--------------------------------------------------------------------------


      SUBROUTINE ADDSPEC2(G0,TAU1,TAU2,TEMP,MP,LIKE,LAMBDA1
     1,LAMBDA2,LAMBDA,LVALUE,NVIB1,NVIB2,ifun)
C
C     THIS PROGRAM GENERATES A LISTING OF THE CIA TR ALFA(OMEGA)
C     IF BOTH LAMBDA1 AND LAMBDA2 ARE NEGATIVE: SINGLE TRANSITIONS;
C     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE.
C     MP=1 GIVES LISTINGS OF INTERMEDIATE RESULTS.
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
	implicit double precision (a-h,o-z)
      COMMON /H2PART1/ Q,WH2(2),B0,D0
	COMMON /H2PART2/ JRANGE1,NORMAL
      COMMON /RESULT1/ NF
      common /result2/ FREQ( 801),ABSCOEF( 801)
      DATA CLOSCHM,BOLTZWN/2.68675484E19,.6950304/
      DATA HBAR,PI,CLIGHT/1.054588757d-27,3.1415926535898,2.9979250E10/
      EH2(N,I)=4395.34*(DFLOAT(N)+0.5)-117.905*(DFLOAT(N)+0.5)**2
     2  +(60.809-
     1  2.993*(DFLOAT(N)+0.5)+.025*(DFLOAT(N)+.5)**2)*DFLOAT(I)-
     3 (.04648-.00134*(DFLOAT(N)+.5))*DFLOAT(I*I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*
     1  DEXP(-1.4387859/T*EH2(0,J*(J+1)))

	do 777 i=1, nf
	abscoef(i) = 0.d0
777	continue

      TWOPIC=2.*PI*CLIGHT
c      IF (MP.NE.1) MP=0
c      IF (LIKE.NE.1) LIKE=0
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB=CALIB/DFLOAT(1+LIKE)
      BETA=1./(BOLTZWN*TEMP)
      LIST=NF
C
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************
C
  110 FORMAT (/, ' LAMBDA,LVALUE=',2I3,'  COMPONENT INCLUDED.',/,
     1  ' LINE SHAPE PARAMETERS:',3E12.3 ,/ )
C
      IF ((LAMBDA1.LT.0).AND.(LAMBDA2.LT.0)) GO TO 60
      JPLUSL=JRANGE1+MAX0(LAMBDA1,LAMBDA2)
      DO 50 I1=1,JRANGE1
         J1=I1-1
      DO 50 IP1=1,JPLUSL
         JP1=IP1-1
         CG1S=CLEBSQR(J1,LAMBDA1,JP1)
         IF (CG1S) 50,50,10
   10    P1=PH2(J1,TEMP)/Q
         IF (P1.LT.0.001) GO TO 50
         OMEGA1=EH2(NVIB1,JP1*IP1)-EH2(0,J1*I1)
         DO 40 I2=1,JRANGE1
            J2=I2-1
         DO 40 IP2=1,JPLUSL
            JP2=IP2-1
            CG2S=CLEBSQR(J2,LAMBDA2,JP2)
            IF (CG2S) 40,40,20
   20       P2=PH2(J2,TEMP)/Q
            IF (P2.LT.0.001) GO TO 40
            OMEGA2=EH2(NVIB2,JP2*IP2)-EH2(0,J2*I2)
            FAC=CALIB*P1*P2*CG1S*CG2S
            DO 30 I=1,LIST
               FRQ=FREQ(I)-OMEGA1-OMEGA2
               WKI=FREQ(I)*(1.-DEXP(-BETA*FREQ(I)))
               WKF=WKI*FAC
               XBG=G0*BGAMA2(FRQ,TAU1,TAU2,TEMP,ifun)
               ABSCOEF(I)=ABSCOEF(I)+XBG*WKF
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      GO TO 100
c	Single transitions here:
   60 JPLUSL=JRANGE1+LAMBDA
      DO 90 I=1,JRANGE1
         J=I-1
      DO 90 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR(J,LAMBDA,JP)
         IF (CGS) 90,90,70
   70    P=PH2(J,TEMP)/Q
         IF (P.LT.0.001) GO TO 90
         OMEGA1=EH2(NVIB1,JP*IP)-EH2(0,J*I)
         FAC=CALIB*P*CGS
         DO 80 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))
            WKF=WKI*FAC
            XBG=G0*BGAMA2(FRQ,TAU1,TAU2,TEMP,ifun)
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   80    CONTINUE
   90 CONTINUE
  100 CONTINUE
      RETURN
C
  140 FORMAT (/,' ABSORPTION COEFFICIENT ALPHA(fnu)',/,
     3  ' FROM',F8.1,' CM-1 TO',F8.1,' CM-1, AT',F6.2,
     3  ' CM-1 INCREMENTS',' IN UNITS OF CM-1 AMAGAT-2',/)
  150 FORMAT (  6E12.4)
C
      END


c--------------------------------------------------------------------------

	subroutine BCparam(temp, G0, G1, G2, tau1, tau2)
	implicit double precision (a-h,o-z)
c	G0,G1,G2 are tree lowest translational moments	
      data HBAR/1.05458875D-27  /
      data BOLTZK/1.38054D-16/
	T = temp
      TAU0=HBAR/(2.*BOLTZK*T)
      TTA=G0*TAU0/G1
      TAU1=DSQRT((G2*TTA-G0*(1.+TAU0**2/TTA))/(G0*(TAU0/TTA)**2))
      TAU2=TTA/TAU1
c	write(3,10)  G0, G1, G2
10	format(' M0, M1, M2:', 3e12.4)
	return
	end


c--------------------------------------------------------------------------

	subroutine K0param(temp, G0,G1,G2, tau1, tau2)
	implicit double precision (a-h,o-z)
c	G0,G1,G2 are tree lowest translational moments	
      data HBAR/1.05458875D-27  /
      data BOLTZK/1.38054D-16/
	T = temp
      TAU0=HBAR/(2.*BOLTZK*T)
      DELT=(TAU0*G1)**2 -4.*(G1*G1/G0+G1/TAU0-G2)*TAU0*TAU0*G0
      if(delt.le.0.d0) go to 88
      TAU1=(-DSQRT(DELT)-TAU0*G1)/(2.*(G1*G1/G0+G1/TAU0-G2))
      if(tau1.le.0.d0) go to 88
      TAU1=DSQRT(TAU1)
      TAU2=TAU0*TAU1*G0/(G1*TAU1*TAU1-TAU0*G0)
c	write(3,10)  G0, G1, G2
10	format(' M0, M1, M2:', 3e12.4)
	return
88    write (6,177) delt, tau10
177   format(' Problem: one of the following is negative for K0, ',
     1  /,' delt, tau1:', 2e13.4)
      stop 9955   
	end

c--------------------------------------------------------------------------

c
      FUNCTION BGAMA2 (FNU,T1,T2,TEMP,ifun)
      implicit double precision (a-h,o-z)
	if(ifun.eq.0) bgama2=bgama0(fnu,t1,t2,temp)
	if(ifun.eq.1) bgama2=bgama1(fnu,t1,t2,temp)
      RETURN
      END

c--------------------------------------------------------------------------

c
      FUNCTION BGAMA0(FNU,TAU5,TAU6,TEMP)
C     K0 LINE SHAPE MODEL
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
      implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/

      TAU4=dSQRT(TAU5*TAU5+(HBOK/(TEMP*2.d0))**2)
      OMEGA=TWOPIC*FNU
      XNOM=1.d0/(TAU6*TAU6)+OMEGA*OMEGA
      X=TAU4*DSQRT(XNOM)
      TAU56=TAU5/TAU6
      TAU56=DMIN1(TAU56,430.d0)
      BGAMA0=(TAU5/PI)*DEXP(TAU56+FNU/(2.d0*BKW*TEMP))*
     1 XK0(X)
      RETURN
      END

c--------------------------------------------------------------------------

      FUNCTION BGAMA1(FNU,TAU1,TAU2,TEMP)
C     BIRNBAUM S CIA LINE SHAPE MODEL (K1)
C     NORMALIZATION SUCH THAT 0TH MOMENT EQUALS UNITY
C     FNU IS THE FREQUENCY IN CM-1; TEMP IN KELVIN.
       implicit double precision (a-h,o-z)
      DATA TWOPIC,BKW,HBOK,PI/1.883651568d11,.6950304256d0,
     1 7.638280918d-12, 3.141592654d0/
      TAU3=dSQRT(TAU2*TAU2+(HBOK/(TEMP*2.))**2)
      OMEGA=TWOPIC*FNU
      DENOM=1.d0+(OMEGA*TAU1)**2
      X=(TAU3/TAU1)*dSQRT(DENOM)
      AAA=TAU2/TAU1
      AAA=dMIN1(AAA,430.d0)
      BGAMA1=(TAU1/PI)*dEXP(AAA+FNU/(2.d0*BKW*TEMP))*
     1 XK1(X)/DENOM
      RETURN
      END
      FUNCTION XK0(X)
C     MODIFIED BESSEL FUNCTION K0(X)
C     ABRAMOWITZ AND STEGUN P.379
       implicit double precision (a-h,o-z)
      IF(X-2.d0) 10,10,20
   10 T=(X/3.75d0)**2
      FI0=(((((.0045813*T+.0360768)*T+.2659732)*T
     1 +1.2067492)*T+3.0899424)*T+3.5156229)*T+1.
      T=(X/2.)**2
      P=(((((.00000740*T+.00010750)*T+.00262698)*T
     1 +.03488590)*T+.23069756)*T+.42278420)*T+
     2 (-.57721566)
      X=DABS(X)
      XK0=-DLOG(X/2.)*FI0+P
      RETURN
   20 T=(2./X)
      P=(((((.00053208*T-.00251540)*T+.00587872)*T
     1 -.01062446)*T+.02189568)*T-.07832358)*T+
     2 1.25331414
      X=DMIN1(X,330.d0)
      XK0=DEXP(-X)*P/DSQRT(X)
      RETURN
      END

c--------------------------------------------------------------------------
 
      FUNCTION XK1(X)
C     MODIFIED BESSEL FUNCTION K1(X) TIMES X
C     PRECISION IS BETTER THAN 2.2e-7 EVERYWHERE.
C     ABRAMOWITZ AND S,TEGUN, P.379; TABLES P.417.
       implicit double precision (a-h,o-z)
      IF(X-2.) 10,10,20
   10 T=(X/3.75)**2
      FI1=X*((((((.00032411*T+.00301532)*T+.02658733)*T+.15084934)
     1 *T+.51498869)*T+.87890594)*T+.5)
      T=(X/2.)**2
      P=(((((-.00004686*T-.00110404)*T-.01919402)*T-.18156897)*T-
     1 .67278579)*T+.15443144)*T+1.
      XK1=X*dLOG(X/2)*FI1+P
      RETURN
   20  T=2./X
      P=(((((-.00068245*T+.00325614)*T-.00780353)*T+.01504268)*T-
     1 .03655620)*T+.23498619)*T+1.25331414
      X=dMIN1(X,330.d0)
      XK1=dSQRT(X)*dEXP(-X)*P
      RETURN
      END

c--------------------------------------------------------------------------
       SUBROUTINE ADDEMCH4(alfatot,gnu0,gnumax,dgnu,tp)
C
C       ****************************************************************
C       PROGRAM PREPARED BY ALEKSANDRA BORYSOW, UNIVERSITY OF TEXAS @AUS
C       AND JOINT INSTITUTE FOR LABORATORY ASTROPHYSICS, UNIV. COLORADO
C       LAST CORRECTION DATE: 14 NOVEMBER 1988
C       THIS PROGRAM IS COMPATIBLE WITH PAPER: A. BORYSOW & L. FROMMHOLD
C       ASTROPHYSICAL JOURNAL; VOL. 304, PP.849-865; 1986.
C       ****************************************************************
C      PROGRAM GENERATES H2-CH4 FREE-FREE, BOUND-FREE & BOUND-BOUND
C      COLLISION INDUCED ABSORPTION SPECTRA
C
C     OUTPUT: file output
C     File corrch4 contains STATISTICAL (NUCLEAR) CORRECTIONS
C     FOR METHANE INDUCTION
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Character*10 LGAS
      real gnu0,gnumax,dgnu,tp
      COMMON /BLOCKIN/ TEMP,FNUMIN,FNUMAX,DNU
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /RSILO/ RSILO(201)
      COMMON /BOU43/ INITB
      COMMON /BB/ OMEG,RSI,RSIGG,ALFA,SCAL,NSOL
      COMMON /BF/ G0BF,DELBF,OM0
      COMMON /STATI/ QW3(51,3),QW4(51,4)
      COMMON /LIKE/ LIKE,LGAS
      DIMENSION FREQ( 801), ABSCOEF( 801), ALFATOT( 801)
      DIMENSION RSI(201), RSIGG(201), TT(2), SS(1), OMEG(201), AIG(201)
c      DATA TEMP/70.D0/
c      DATA FNUMIN/10.D0/
c      DATA FNUMAX/400.d0/
      DATA SLIT/4.3d0/
c      DATA DNU/5.0D0/
      DATA LIKE/0/
      DATA LGAS/' H2 - CH4' /
      Y(X,A,B,C)=A*DEXP((C*X+B)*X)
C
	OPEN(UNIT=33, FILE='corrch4.txt', STATUS='OLD')
        
      temp=tp
c
      fnumin=gnu0
      fnumax=gnumax
      dnu=dgnu

      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1
      IF (NF.GT. 801) NF= 801
      FNUMAX=FNUMIN+DFLOAT(NF-1)*DNU
      CALL PARTSUMCH (TEMP)
C
C     THE H2 - CH4 SPECTRA   FOR 50-300K
C     =================
C
C     ONLY B-F + F-F TERMS INCLUDED IN MODELLING
C     TEMPERATURE INTERPOLATION GOOD FOR 45-300K
C
      X=DLOG(TEMP)
      DO 10 I=1,NF
         FREQ(I)=FNUMIN+DFLOAT(I-1)*DNU
         ALFATOT(I)=0.0
   10 ABSCOEF(I)=0.
C
C     READ FROM FILE corrch4 /HERE TAPE 33/ THE QUANTUM CORRECTIONS
C     (12 Y(J,JP,LAM)-1)
C
      REWIND 33
      DO 30 N=1,3
         DO 20 I=1,5
            K=10*(I-1)+1
            K1=10*I
            READ (33,160) (QW3(J,N),J=K,K1)
   20    CONTINUE
         READ (33,150) QW3(51,N)
   30 CONTINUE
      DO 50 N=1,4
         DO 40 I=1,5
            K=(I-1)*10+1
            K1=10*I
            READ (33,160) (QW4(J,N),J=K,K1)
   40    CONTINUE
         READ (33,160) QW4(51,N)
   50 CONTINUE
	CLOSE(33)
C
C     THESE STATISTICAL FACTORS ARE VALID ONLY FOR DELTA J>0
C     FOR DELTA J<0 THESE WILL BE EQUAL TO 1.
C     HYDROGEN'S QUADRUPOLE
C
      EPS=1.D-5
      TT(1)=10.D0
      CALL BOUND32 (TEMP,RSI,NSOL)
      DO 60 I=1,NSOL
         RSILO(I)=DLOG(RSI(I)*1.D80)
   60 OMEG(I)=DFLOAT(I-1)*DX
      CALL SPLINE (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
      SCAL=1.D0
      S=Y(X,.22916D-59,-1.27134D0, .12423D0)
      E=Y(X,-.2019D1,-.02259d0,-.05891D0)
      T1=Y(X,.48254D-13,.6421d0,-.10109D0)
      T2=Y(X,.97826D-12,-.48654d0,-.0146D0)
      T3=Y(X,.35247D-12,.10164d0,-.07879D0)
      T4=Y(X,.13961D-13,1.11146d0,-.09561D0)
C
C     THIS PART FOR MODELING A LOW FREQUENCY PART OF BOUND-FREE
C     TRANSLATIONAL SPECTRAL FUNCTION, BY A DESYMMETRIZED GAUSSIAN
C     PROFILE
C
      G0BF=Y(X,0.73536D-72,-.79815D0,-.0585D0)
      DELBF=Y(X,3.123D0,-.00178D0,0.00021D0)
      OM0=Y(X,8.6922D0,0.00785D0,-0.00054D0)
      CALL ADDSPEC (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,2,0,2,3)
      DO 70 I=1,NF
   70 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)
C
C     PARAMETERS FOR 4045  (HYDROGEN HEXADECAPOLE) COMPONENTS
C
      CALL BOUN54C (TEMP,RSI,NSOL)
      DO 80 I=1,NSOL
         RSILO(I)=DLOG(RSI(I)*1.D80)
   80 OMEG(I)=DFLOAT(I-1)*DX
      CALL SPLINE (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
      S=Y(X,(2.7321D-60/9.)*0.038,-2.32012D0,.23082D0)
      E=Y(X,-1.8198D0,-.00665D0,-.05626D0)
      T1=Y(X,1.3492D-13,.14472D0,-.06506D0)
      T2=Y(X,1.6596D-12,-.77679D0,.01401D0)
      T3=Y(X,5.9914D-13,-.16208D0,-.05895D0)
      T4=Y(X,1.9405D-14,.95965D0,-.10246D0)
      SCAL=0.038D0
      CALL ADDSPEC (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,4,0,4,5)
      DO 90 I=1,NF
   90 ALFATOT(I)=ALFATOT(I)+ABSCOEF(I)
C
C     METHANE INDUCED COMPONENTS
C     ===============================
C     OCTOPOLE- INDUCED TERM (43)
C
      EPS=1.D-5
      TT(1)=10.D0
      CALL BOUND43 (TEMP,RSI,NSOL)
      DO 100 I=1,NSOL
         RSILO(I)=DLOG(RSI(I)*1.D80)
  100 OMEG(I)=FLOAT(I+INITB-2)*DX
      CALL SPLINE (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
      SCAL=1.D0
      S=Y(X,.91196D-60/7.,-1.56529D0,.15284D0)
      E=Y(X,-.5D0,0.D0,0.D0)
      T1=Y(X,.51456D-13,.60523D0,-.10856D0)
      T2=Y(X,.62514D-12,-.51384D0,.00754D0)
      T3=Y(X,.55346D-12,-.40381D0,.00208D0)
      T4=Y(X,.13804D-13,1.9307D0,-.27921D0)
      CALL ADSPEC1 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,3,3,4)
      DO 110 I=1,NF
  110 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)
C
C     HEXADECAPOLE-INDUCED TERM
C     =============================
C
      EPS=1.D-5
      TT(1)=10.D0
      CALL BOUN54C (TEMP,RSI,NSOL)
      DO 120 I=1,NSOL
         RSILO(I)=DLOG(RSI(I)*1.D80)
  120 OMEG(I)=DFLOAT(I-1)*DX
      CALL SPLINE (NSOL,1,0,EPS,OMEG,RSILO,TT,SS,SI,NR,RSIGG)
      SCAL=1.D0
      S=Y(X,2.7321D-60/9.,-2.32012D0,.23082D0)
      E=Y(X,-1.8198D0,-.00665D0,-.05626D0)
      T1=Y(X,1.3492D-13,.14472D0,-.06506D0)
      T2=Y(X,1.6596D-12,-.77679D0,.01401D0)
      T3=Y(X,5.9914D-13,-.16208D0,-.05895D0)
      T4=Y(X,1.9405D-14,.95965D0,-.10246D0)
      CALL ADSPEC1 (S,E,T1,T2,T3,T4,TEMP,NF,FREQ,ABSCOEF,0,LIKE,0,4,4,5)
      DO 130 I=1,NF
  130 ALFATOT(I)=ABSCOEF(I)+ALFATOT(I)
C
  140 FORMAT (' ABSORPTION SPECTRA OF ', A10,' AT', F8.1,
     1  'K',/1X,
     1 ' MIN.FREQ.=', F8.1, ' CM-1', 10X, ' MAX.FREQ.=', F8.1,
     2 ' CM-1',/, ' FREQ.INCREMENT=', F8.2, ' CM-1', 5X, 'IN', I5,
     1  ' STEPS',/)
  150 FORMAT (F8.5)
  160 FORMAT (10F8.5)
  170 FORMAT (/, ' ABSORPTION COEFFICIENT ALPHA(FNU), FROM ',F5.1,
     1' CM-1 TO', F7.1,  9H CM-1, AT, F6.2, /,
     1  23H CM-1 INCREMENTS, AT T=,F7.2
     2, 29H K, IN UNITS OF CM-1 AMAGAT-2,/)
  180 FORMAT (  1H ,10E13.5)
  190 FORMAT (4F10.5,I5)
C
      RETURN
      END
      SUBROUTINE ADDSPEC (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,ABSCOEF
     1,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
C
C     THIS ENTRY FOR L I N E A R   M O L E C U L E!!!!!
C     SET LAMBDA1 NONZERO
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)
C     IF EITHER LAMBDA1 OR LAMBDA2 EQUAL TO ZERO - SINGLE TRANSITIONS;
C     LAMBDA1 CORRESPONDS TO H2, LAMBDA2 TO CH4
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BB/ OMEG(201),RSI(201),RSIGG(201),BETA,SCAL,NSOL
      COMMON /RSILO/ RSILO(201)
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /H2PART/ Q,WH2(2),B0,D0,JRANGE1
      COMMON /BF/ G0BF,DELBF,OM0
      COMMON /CHPART/ Q1,WCH(2),B01,D01,JRANGE2
      DIMENSION ABSCOEF(NF), FREQ(NF)
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/
      EH2(I)=(B0-DFLOAT(I)*D0)*DFLOAT(I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859/T*
     1 EH2(J*(J+1)))
      TWOPIC=2.*PI*CLIGHT
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB=CALIB/DFLOAT(1+LIKE)
      BETA=1./(BOLTZWN*TEMP)
      LIST=NF
      DO 10 I=1,LIST
   10 ABSCOEF(I)=0.0
C
C     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   *******************
C
      IF (LAMBDA1.EQ.0) STOP 6666
C
C     FOR THIS ENTRY LAMBDA1 HAS TO BE NONZERO
C     SINGLE TRANSITIONS ON HYDROGEN'S ROTATIONAL FREQUENCIES.
C     LAMBDA IS EQUAL TO LAMBDA1 FOR SINGLE TRANSITIONS
C
      JPLUSL=JRANGE1+LAMBDA
      DO 40 I=1,JRANGE1
         J=I-1
      DO 40 IP=1,JPLUSL
         JP=IP-1
         CGS=CLEBSQR(J,LAMBDA,JP)
         IF (CGS) 40,40,20
   20    P=PH2(J,TEMP)/Q
         OMEGA1=EH2(JP*IP)-EH2(J*I)
         FAC=CALIB*P*CGS

         DO 30 IQ=1,LIST
            FRQ=FREQ(IQ)-OMEGA1
            WKI=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))
            WKF=WKI*FAC
            XBG=G0*BGAMACH(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)


            IF ( DABS(FRQ).LE.WNRMAX3 )
     1       XBG=XBG+SCAL*SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)

            IF (LVALUE.EQ.3 .AND. G0BF.NE.0.D0)
     1      XBG = XBG + BGAUS(FRQ, G0BF, DELBF, OM0, TEMP)

C
C     THIS MODELES THE PART OF THE BOUND-FREE SPECTRAL FUNCTION BY
C     MEANS OF THE DESYMMETRIZED GAUSSIAN PROFILE
C
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   30    CONTINUE
   40 CONTINUE
      RETURN
C
   50 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 12H COM
     1PONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,  5HG(0)=,E12.3
     2/)
   60 FORMAT ((1X,10E12.4,/))
C
      END
      SUBROUTINE ADSPEC1 (G0,EP,TAU1,TAU2,TAU5,TAU6,TEMP,NF,FREQ,ABSCOEF
     1,MP,LIKE,LAMBDA1,LAMBDA2,LAMBDA,LVALUE)
C
C     FOR INDUCTION BY A  T E T R A H E D R A L   M O L E C U L E!!
C     THIS PROGRAM GENERATES LISTING OF CIA TR ALFA(OMEGA)
C     LAMBDA1 CORRESPONDS TO H2, LAMBDA2 TO CH4
C     FOR THIS ENTRY SET LAMBDA2 NONZERO
C     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE.
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BB/ OMEG(201),RSI(201),RSIGG(201),BETA,SCAL,NSOL
      COMMON /RSILO/ RSILO(201)
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /H2PART/ Q,WH2(2),B0,D0,JRANGE1
      COMMON /CHPART/ Q1,WCH(2),B01,D01,JRANGE2
      COMMON /STATI/ QW3(51,3),QW4(51,4)
      DIMENSION ABSCOEF(NF), FREQ(NF)
      DATA CLOSCHM,BOLTZWN/2.68675484D19,.6950304/
      DATA HBAR,PI,CLIGHT/1.054588757D-27,3.1415926535898,2.9979250D10/
      ECH4(I)=B01*DFLOAT(I)
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))
      TWOPIC=2.*PI*CLIGHT
      CALIB=TWOPIC*((4.*PI**2)/(3.*HBAR*CLIGHT))*CLOSCHM**2
      CALIB=CALIB/DFLOAT(1+LIKE)
      BETA=1./(BOLTZWN*TEMP)
      LIST=NF
      DO 10 I=1,LIST
   10 ABSCOEF(I)=0.0
C
C     ROTATIONAL SPECTRUM: DETAILED LISTING   *******************
C
C
      DO 60 I=1,JRANGE2
         J=I-1
         P=DFLOAT(2*J+1)*PCH4(J,TEMP)/Q1
         DO 40 IP=1,LAMBDA
C
C     POSITIVE DELTA J
C
            JP=J+IP
            OMEGA1=ECH4(JP*(JP+1))-ECH4(J*I)
            IF (LAMBDA.EQ.3) CC=(1.+QW3(I,IP)/4.)
            IF (LAMBDA.EQ.4) CC=(1.+QW4(I,IP)/4.)
            FAC=CALIB*P*CC
            DO 20 IQ=1,LIST
               FRQ=FREQ(IQ)-OMEGA1
               WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1         (2.*DFLOAT(JP)+1.)
               XBG=G0*BGAMACH(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)

               IF ( DABS(FRQ).LE.WNRMAX3 )
     1         XBG= XBG + SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
               ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   20       CONTINUE
C
C     NEGATIVE DELTA J
C
            JP=J-IP
            IF (JP.LT.0) GO TO 40
            OMEGA1=ECH4(JP*(JP+1))-ECH4(J*I)
            FAC=CALIB*P
            DO 30 IQ=1,LIST
               FRQ=FREQ(IQ)-OMEGA1
               WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1         (2.*DFLOAT(JP)+1.)
               XBG=G0*BGAMACH(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)

               IF ( DABS(FRQ).LE.WNRMAX3 )
     1     XBG=XBG+SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
               ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   30       CONTINUE
   40    CONTINUE
C
C     DELTA J=0
C
         JP=J
         FAC=CALIB*P
         DO 50 IQ=1,LIST
            FRQ=FREQ(IQ)
            WKF=FREQ(IQ)*(1.-DEXP(-BETA*FREQ(IQ)))*FAC*
     1      (2.*DFLOAT(JP)+1.)
            XBG=G0*BGAMACH(FRQ,TAU1,TAU2,EP,TAU5,TAU6,TEMP)

            IF (DABS(FRQ).LE.WNRMAX3)
     1   XBG=XBG+SPECFCT(FRQ,OMEG,RSILO,RSIGG,NSOL,BETA)
            ABSCOEF(IQ)=ABSCOEF(IQ)+XBG*WKF
   50    CONTINUE
   60 CONTINUE
      RETURN
C
   70 FORMAT (/, 32H LAMBDA1,LAMBDA2, LAMBDA,LVALUE=,2I3,2X,2I3, 12H COM
     1PONENT .,/15X, 22HLINE SHAPE PARAMETERS:,6E12.3,5X,  5HG(0)=,E12.3
     2/)
   80 FORMAT ((1X,10E12.4,/))
C
      END
c      FUNCTION CLEBSQR (L,LAMBDA,LP)
cC
cC     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0)
cC     FOR INTEGER ARGUMENTS ONLY
cC     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO.
cC
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      FC=DFLOAT(2*LP+1)
c      GO TO 10
cC
c      ENTRY THREEJ2
cC
cC     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP
cC                                                 0    0    0
cC     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT
cC     (LIMITATION TO INTEGER ARGUMENTS ONLY)
cC
cC     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE
cC     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT
cC     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS.
cC
c      FC=1.
c   10 CLEBSQR=0.
c      IF (((L+LAMBDA).LT.LP).OR.((LAMBDA+LP).LT.L).OR.((L+LP).LT.LAMBDA)
c     1) RETURN
c      IF (MOD(L+LP+LAMBDA,2).NE.0) RETURN
c      IF ((L.LT.0).OR.(LP.LT.0).OR.(LAMBDA.LT.0)) RETURN
c      F=1./DFLOAT(L+LP+1-LAMBDA)
c      IF (LAMBDA.EQ.0) GO TO 30
c      I1=(L+LP+LAMBDA)/2
c      I0=(L+LP-LAMBDA)/2+1
c      DO 20 I=I0,I1
c   20 F=F*DFLOAT(I)/DFLOAT(2*(2*I+1))
c   30 P=FC*F*FCTL(LAMBDA+L-LP)*FCTL(LAMBDA+LP-L)
c      CLEBSQR=P/(FCTL((LAMBDA+L-LP)/2)*FCTL((LAMBDA+LP-L)/2))**2
c      RETURN
cC
c      END
c      FUNCTION FCTL (N)
c       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      P(Z)=((((-2.294720936D-4)/Z-(2.681327160D-3))/Z+(3.472222222D-3))/
c     1Z+(8.333333333D-2))/Z+1.
c      FCTL=1.
c      IF (N.LE.1) RETURN
c      IF (N.GT.15) GO TO 20
c      J=1
c      DO 10 I=2,N
c   10 J=J*I
c      FCTL=DFLOAT(J)
c      RETURN
c   20 Z=DFLOAT(N+1)
c      FCTL=(DEXP(-Z)*(Z**(Z-0.5))*P(Z)*2.5066282746310)
c      RETURN
cC
c      END
      FUNCTION BGAMACH (FNU,T1,T2,EPS,T3,T4,TEMP)
C
C     EQUATION 13, SO-CALLED EBC MODEL, OF BORYSOW,TRAFTON,FROMMHOLD,
C     AND BIRNBAUM, ASTROPHYS.J., (1985)
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 K0
      DATA PI,CLIGHT/3.1415926535898,29979245800./
      DATA HBAR,BOLTZ/1.0545887D-27,1.380662D-16/
      P1(X)=((((((.0045813*X+.0360768)*X+.2659732)*X+1.2067492)*X+3.0899
     1424)*X+3.5156229)*X+1.)
      P2(X)=((((((.00000740*X+.00010750)*X+.00262698)*X+.03488590)*X+.23
     1069756)*X+.42278420)*X-.57721566)
      P3(X)=((((((.00032411*X+.00301532)*X+.02658733)*X+.15084934)*X+.51
     1498869)*X+.87890594)*X+.5)
      P4(X)=((((((-.00004686*X-.00110404)*X-.01919402)*X-.18156897)*X-.6
     17278579)*X+.15443144)*X+1.)
      P5(X)=((((((.00053208*X-.00251540)*X+.00587872)*X-.01062446)*X+.02
     1189568)*X-.07832358)*X+1.25331414)
      P6(X)=((((((-.00068245*X+.00325614)*X-.00780353)*X+.01504268)*X-.0
     13655620)*X+.23498619)*X+1.25331414)
C
      OMEGA=2.*PI*CLIGHT*FNU
      T0=HBAR/(2.*BOLTZ*TEMP)
      Z=DSQRT((1.+(OMEGA*T1)**2)*(T2*T2+T0*T0))/T1
      IF (Z-2.) 10,10,20
   10 XK1=Z*Z*DLOG(Z/2.)*P3((Z/3.75)**2)+P4((Z/2.)**2)
      GO TO 30
   20 XK1=DSQRT(Z)*DEXP(-Z)*P6(2./Z)
   30 ZP=DSQRT((1.+(OMEGA*T4)**2)*(T3*T3+T0*T0))/T4
      IF (ZP-2.) 40,40,50
   40 K0=-DLOG(ZP/2.)*P1((ZP/3.75)**2)+P2((ZP/2.)**2)
      GO TO 60
   50 K0=DEXP(-ZP)*P5(2./ZP)/DSQRT(ZP)
   60 BGAMACH=((T1/PI)*DEXP(T2/T1+T0*OMEGA)*XK1/
     1   (1.+(T1*OMEGA)**2)+EPS*(T3/PI)*DEXP(T3/T4+T0*OMEGA)*K0)/
     1   (1.+EPS)
      RETURN
      END

      SUBROUTINE PARTSUMCH (TEMP)
C
C     H2 ROTATIONAL PARTITION SUM Q = Q(T).
C
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /H2PART/ Q,WH2(2),B0,D0,JRANGE1
      COMMON /CHPART/ Q1,WCH(2),B01,D01,JRANGE2
      DATA B0,D0,WH2(1),WH2(2)/59.3392,0.04599,1.,3./
      DATA B01,D01,WCH(1),WCH(2)/5.24,0.,1.,1./
      EH2(I)=(B0-DFLOAT(I)*D0)*DFLOAT(I)
      PH2(J,T)=DFLOAT(2*J+1)*WH2(1+MOD(J,2))*DEXP(-1.4387859*
     1 EH2(J*(J+1))/T)
      ECH4(I)=B01*DFLOAT(I)
      PCH4(J,T)=DEXP(-1.4387859/T*ECH4(J*(J+1)))
C
C     Q,B01,D01,WCH - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR CH4
C     Q,B0,D0,WH2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR H2
C
      Q=0.
      J=0
   10 DQ=PH2(J,TEMP)
      Q=Q+DQ
      J=J+1
      IF (DQ.GT.Q/900.) GO TO 10
      JRANGE1=J
C
C     *** PARTITION FUNCTION FOR CH4 *******************
C
      Q1=0.
      J=0
   20 DQ=PCH4(J,TEMP)*DFLOAT((2*J+1)**2)
      Q1=Q1+DQ
      J=J+1
      IF (DQ.GT.Q1/1000.) GO TO 20
      JRANGE2=J+3
      RETURN
C
   30 FORMAT (/, 40H ROTATIONAL PARTITION FUNCTION OF H2: Q=,F8.2,10X,
     17HJ MAX =,I3/)
C
      END
      SUBROUTINE PROFILE (X,Y)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BL3/ RSI(401)
C
C     ATRIANGULAR SLIT FUNCTION IS USED.
C
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      IF (Y) 50,60,10
   10 CONTINUE
      X0=(NSRI+1.)+X/DX
      NC=X0
      N1=NC+1
      SLOPE=Y/SLIT
      NU=X0-NS
      IF (NU.LT.1) NU=1
      IF (NU.GT.NSRIUP) RETURN
      NO=X0+NS
      IF (NO.GT.NSRIUP) NO=NSRIUP
      IF (NO.LT.1) RETURN
      IF (NC.GT.NSRIUP) NC=NSRIUP
      IF (NC.LE.1) GO TO 30
      DO 20 I=NU,NC
         XI=(I-1.)*DX-WNRMAX3
         DR=SLOPE*(XI-(X-SLIT))
         IF (DR.LE.0.) GO TO 20
         RSI(I)=RSI(I)+DR
   20 CONTINUE
   30 IF (NC.GE.NSRIUP) RETURN
      IF (N1.LT.1) N1=1
      DO 40 I=N1,NO
         XI=(I-1.)*DX-WNRMAX3
         DR=Y-SLOPE*(XI-X)
         IF (DR.LE.0.) GO TO 40
         RSI(I)=RSI(I)+DR
   40 CONTINUE
      RETURN
   50 WRITE(1,70) SLIT
   60 CONTINUE
      RETURN
C
   70 FORMAT (/, 30H A TRIANGULAR SLIT FUNCTION OF,F6.3, 23H CM-1 HALFWI
     1DTH IS USED,/)
C
      END
      FUNCTION SPECFCT (FREQ,OMEGA,PHI,PHI2,N,RTEMP)
C
C     THIS INTERPOLATES THE SPECTRAL FUNCTION PHI(FREQ) DEFINED AT
C     OMEGA(N) AS PHI(N). PHI2 IS THE SECOND DERIVATIVE AT OMEGA
C     WHICH MUST BE OBTAINED FIRST (USE SPLINE FOR THAT PURPOSE).
C     RTEMP IS THE RECIPROCAL TEMPERATURE IN CM-1 UNITS.
C     NOTE THAT WE INTERPOLATE 1.E80 TIMES THE LOGARITHM OF PHI(OMEGA)
C
	IMPLICIT  DOUBLE  PRECISION (A-h,o-z)
      DIMENSION PHI(N), PHI2(N), OMEGA(N)
	dimension f(1), gp(1)
c	correction above  1-MAR-1989 18:10:02 in all subroutine
      TFAC=0.
      F(1)=FREQ
      IF (F(1) ) 10,20,20
   10 F(1)=dABS(F(1))
      TFAC=(-RTEMP*F(1))
   20 IF (F(1).LE.OMEGA(N)) GO TO 30
      SPECFCT=DEXP(-(PHI(N-1)-PHI(N))*(F(1)-OMEGA(N))/
     1  (OMEGA(N)-OMEGA(N-1))+
     1PHI(N)+TFAC)*(1.D-80)
      RETURN
c	exchanged ixpolat to spline   1-MAR-1989 18:08:27
   30 CALL spline  (N,1,0,1.D-6,OMEGA,PHI,F,GP,SI,NR,PHI2)
c	f(1), gp(1)
      SPECFCT=DEXP(TFAC+GP(1))*(1.D-80)
      RETURN
C
      END

      SUBROUTINE SPLINE (L,M,K,EPS,X,Y,T,SS,SI,NR,S2)
C
C     SPLINE INTERPOLATION AND QUADRATURE, THIRD ORDER AFTER GREVILLE.
C     INPUT ARGUMENTS L...Y, OUTPUT SS...NR.
C     L DATA POINTS X(1), Y(1) ... X(L),Y(L)
C     EPS=ERROR CRITERION, TYPICALLY EPS=1.E-5 FOR 5 DECI. PLACES ACCURA
C     M ARGUMENTS T(1)..T(M) FOR WHICH FUNCTION VALUES SS(1)..SS(M), FOR
C     K=0; OR FIRST OR SECOND DERIVATIVE FOR K=1 OR -1, RESPECTIVELY.
C     NOTE THAT M HAS TO BE AT LEAST EQUAL TO 1.
C     SI=INTEGRAL (OVER WHOLE INTERVAL) FOR K=2 ONLY.
C     FOR 'NATURAL' SPLINE FUNCTIONS, S2(1)=S2(L)=0. MUST BE INPUT*NOTE*
C     N0 INDICATES THE NUMBER OF OUT-OF-RANGE CALLS. X(1)@T(I)@X(L)
C     EXTRAPOLATE WITH CAUTION. (ASSUMPTION D2Y/DX2 = 0.)
C     S2(I) IS THE 2ND DERIVATIVE AT X=X(I) AND IS COMPUTED INTERNALLY.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(L), Y(L), T(M), SS(M), S2(L)
      DELY(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))
      B(I)=(X(I)-X(I-1))*0.5/(X(I+1)-X(I-1))
      C(I)=3.*(DELY(I)-DELY(I-1))/(X(I+1)-X(I-1))
      N=L
      N1=N-1
      NR=0
      DO 10 I=2,N1
   10 S2(I)=C(I)/1.5
      OMEGA=1.0717968
      IC=0
C
C     'NATURAL' SPLINE FUNCTIONS OF THIRD ORDER.
C
      S2(N)=0.
      S2(1)=S2(N)
   20 ETA=0.
      IC=IC+1
      SM=DABS(S2(1))
      DO 30 I=2,N
         IF (DABS(S2(I)).GT.SM) SM=DABS(S2(I))
   30 CONTINUE
      EPSI=EPS*SM
      DO 50 I=2,N1
         W=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*OMEGA
         IF (DABS(W)-ETA) 50,50,40
   40    ETA=DABS(W)
   50 S2(I)=S2(I)+W
      IF (ETA-EPSI) 60,20,20
C      ENTRY IXPOLAT
C
C     THIS ENTRY USEFUL WHEN ITERATION PREVIOUSLY COMPLETED
C
C      N=L
C      N1=N-1
C      NR=0
C      IC=-1
   60 IF (K.EQ.2) GO TO 260
      GO TO 70
   70 DO 250 J=1,M
         I=1
         IF (T(J)-X(1)) 110,210,80
   80    IF (T(J)-X(N)) 100,190,150
   90    IF (T(J)-X(I)) 200,210,100
  100    I=I+1
         GO TO 90
  110    NR=NR+1
         HT1=T(J)-X(1)
         HT2=T(J)-X(2)
         YP1=DELY(1)+(X(1)-X(2))*(2.*S2(1)+S2(2))/6.
         IF (K) 140,130,120
  120    SS(J)=YP1+HT1*S2(1)
         GO TO 250
  130    SS(J)=Y(1)+YP1*HT1+S2(1)*HT1*HT1/2.
         GO TO 250
  140    SS(J)=S2(I)
         GO TO 250
  150    HT2=T(J)-X(N)
         HT1=T(J)-X(N1)
         NR=NR+1
         YPN=DELY(N1)+(X(N)-X(N1))*(S2(N1)+2.*S2(N))/6.
         IF (K) 180,170,160
  160    SS(J)=YPN+HT2*S2(N)
         GO TO 250
  170    SS(J)=Y(N)+YPN*HT2+S2(N)*HT2*HT2/2.
         GO TO 250
  180    SS(J)=S2(N)
         GO TO 250
  190    I=N
  200    I=I-1
  210    HT1=T(J)-X(I)
         HT2=T(J)-X(I+1)
         PROD=HT1*HT2
         S3=(S2(I+1)-S2(I))/(X(I+1)-X(I))
         SS2=S2(I)+HT1*S3
         DELSQS=(S2(I)+S2(I+1)+SS2)/6.
         IF (K) 240,220,230
  220    SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS
         GO TO 250
  230    SS(J)=DELY(I)+(HT1+HT2)*DELSQS+PROD*S3/6.
         GO TO 250
  240    SS(J)=SS2
  250 CONTINUE
  260 SI=0.
      DO 270 I=1,N1
         H=X(I+1)-X(I)
  270 SI=SI+0.5*H*(Y(I)+Y(I+1))-H**3*(S2(I)+S2(I+1))/24.
      IF (K.EQ.2) NR=IC
      RETURN
      END
      SUBROUTINE BOUND32 (TEMP,RSI,NSOL)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)
C
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0
C
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /BL3/ RSIBB(401)
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,
     1  -14.669,-9.0915,-2.564,4.7286,12.4947/
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./
      DATA TWOPIC/1.88365183D11/
      NSRI=190
      WNRMAX3=25.D0
C
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH
C     B-B SPECTRUM, IN [CM-1].
C
      NSRIUP=2*NSRI+1
      DX=WNRMAX3/DFLOAT(NSRI)
      DO 10 I=1,401
   10 RSIBB(I)=0.0
      NS=INT(SLIT/DX)
      DO 20 I=1,9
      DO 20 J=1,9
   20 A(I,J)=0.0
C
C     A(I,J) FOR L=32 CONTRIBUTION H2-CH4 SYSTEM
C
      A(1,4)=.14127D-39
      A(2,3)=.18414D-39
      A(2,5)=.23414D-39
      A(3,4)=.18534D-39
      A(3,6)=.30894D-39
      A(4,5)=.21745D-39
      A(4,7)=.36564D-39
      A(5,6)=.24659D-39
      A(5,8)=.40032D-39
      A(6,7)=.26656D-39
      A(7,8)=.27376D-39
      A(8,9)=.26805D-39
      A(6,9)=.41233D-39
      ALFA=1./(0.69519*TEMP)
      RM=1.7768*1.672649D-24
      PI=3.141592654
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)
C
C     DELTA V = 0 BELOW
C
      DO 40 L=1,8
         STOKE1=EB(1,L+1)-EB(1,L)
C
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0
C
         IF (L.GT.6) GO TO 30
         STOKE3=EB(1,L+3)-EB(1,L)
C
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0
C
         STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE3,STOKI)
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF
         CALL PROFILE (-STOKE3,STOKIP)
   30    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE1,STOKI)
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF
         CALL PROFILE (-STOKE1,STOKIP)
   40 CONTINUE
C
C     DELTA V=1 BELOW
C
      AV1=.39933D-41
      STOKE=22.5287
      STOKI=AV1*DEXP(-ALFA*EB(1,3))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=0.3079D-41
      STOKE=18.6921
      STOKI=AV1*DEXP(-ALFA*EB(1,4))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=0.4281D-41
      STOKE=14.544
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      DO 50 N=1,NSRIUP
   50 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT
      NSOL=NSRI+1
      DO 60 I=1,NSOL
      K=(NSRI+1)+(I-1)
   60 RSI(I)=RSIBB(K)
C
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS
C
      RETURN
C
      END
      SUBROUTINE BOUND43 (TEMP,RSI,NSOL)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)
C
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0
C
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /BL3/ RSIBB(401)
      COMMON /BOU43/ INITB
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,-14.669,-
     19.0915,-2.564,4.7286,12.4947/
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./
      DATA TWOPIC/1.88365183D11/
      NSRI=190
      WNRMAX3=30.
      INITB=1
C
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH
C     B-B SPECTRUM, IN [CM-1].
C
      NSRIUP=2*NSRI+1
      DX=WNRMAX3/DFLOAT(NSRI)
      DO 10 I=1,401
   10 RSIBB(I)=0.0
      NS=INT(SLIT/DX)
      DO 20 I=1,9
      DO 20 J=1,9
   20 A(I,J)=0.0
C
C     A(I,J) FOR L,Lambda={3,4} CONTRIBUTION H2-CH4 SYSTEM
C     MATRIX ELEMENTS AS IF M(ll')/7
C
      A(1,5)=.33038D-41
      A(2,4)=.45246D-41
      A(2,6)=.52213D-41
      A(3,5)=.42125D-41
      A(3,7)=.65899D-41
      A(4,6)=.46727D-41
      A(4,8)=.74203D-41
      A(5,7)=.50413D-41
      A(5,9)=.77286D-41
      A(6,8)=.51657D-41
      A(7,9)=.50297D-41
      ALFA=1./(0.69519*TEMP)
      RM=1.7768*1.672649D-24
      PI=3.141592654
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)
C
C     DELTA V = 0 BELOW
C
      DO 40 L=1,7
         STOKE1=EB(1,L+2)-EB(1,L)
C
C     FREQUENCY SHIFT FOR DELTA J=+2,DELTA V=0
C
         IF (L.GT.5) GO TO 30
         STOKE3=EB(1,L+4)-EB(1,L)
C
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+4, DELTA V=0
C
         STOKI=A(L,L+4)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE3,STOKI)
         STOKIP=A(L,L+4)*DEXP(-ALFA*EB(1,L+4))/PF
         CALL PROFILE (-STOKE3,STOKIP)
   30    STOKI=A(L,L+2)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE1,STOKI)
         STOKIP=A(L,L+2)*DEXP(-ALFA*EB(1,L+2))/PF
         CALL PROFILE (-STOKE1,STOKIP)
   40 CONTINUE
C
C     DELTA V=1 BELOW
C
      AV1=.12153D-42
      STOKE=19.0835
      STOKI=AV1*DEXP(-ALFA*EB(1,4))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=0.89592D-43
      STOKE=14.1528
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=0.1182D-42
      STOKE=8.9669
      STOKI=AV1*DEXP(-ALFA*EB(1,6))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      DO 50 N=1,NSRIUP
   50 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT
      NSOL=NSRI+1
      DO 60 I=1,NSOL
      K=(NSRI+1)+(I-1)
   60 RSI(I)=RSIBB(K)
C
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS
C
      I=0
   70 I=I+1
      IF (RSI(I).EQ.0.) GO TO 70
      IF (I.EQ.1) GO TO 90
      NSOL=NSOL-I+1
      DO 80 J=1,NSOL
   80 RSI(J)=RSI(J+I-1)
      INITB=I
   90 CONTINUE
      RETURN
      END

      SUBROUTINE BOUN54C (TEMP,RSI,NSOL)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EB(2,9), A(9,9), RSI(201)
C
C     EB(I,K) - BOUND ENERGIES FOR VIB.NR.V=I-1, ROT.NR J=K-1
C     A(I,J)=( C(I,L,J)**2) * (MTX.EL. (I,BETA,J)**2 )
C     I,J-DENOTE ROTATIONAL STATES FOR VIB. LEVEL V=0
C
      COMMON /APP3/ SLIT,DX,WNRMAX3,NSRI,NS,NSRIUP,JSLIT
      COMMON /BL3/ RSIBB(401)
      DATA (EB(1,J),J=1,9)/-26.1308,-24.9681,-22.6532,-19.2080,-14.669,-
     19.0915,-2.564,4.7286,12.4947/
      DATA (EB(2,J),J=1,9)/-0.516,-0.1246,0.,0.,0.,0.,0.,0.,0./
      DATA TWOPIC/1.88365183D11/
      NSRI=190
      WNRMAX3=34.
C
C     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN
C     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION
C     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH
C     B-B SPECTRUM, IN [CM-1].
C
      NSRIUP=2*NSRI+1
      DX=WNRMAX3/DFLOAT(NSRI)
      DO 10 I=1,401
   10 RSIBB(I)=0.0
      NS=INT(SLIT/DX)
      DO 20 I=1,9
      DO 20 J=1,9
   20 A(I,J)=0.0
C
C     A(I,J) FOR L=54(CH4) CONTRIBUTION H2-CH4 SYSTEM
C     MATRIX ELEMENTS AS IF M(LL')/9
C
      A(1,6)=.99184D-42
      A(2,5)=.14168D-41
      A(2,7)=.14969D-41
      A(3,4)=.16093D-41
      A(3,6)=.12455D-41
      A(4,5)=.12587D-41
      A(4,7)=.13120D-41
      A(5,6)=.12714D-41
      A(5,8)=.13379D-41
      A(6,7)=.12818D-41
      A(7,8)=.12378D-41
      A(6,9)=.12963D-41
      A(8,9)=.11352D-41
      A(3,8)=.17927D-41
      A(4,9)=.19144D-41
      ALFA=1./(0.69519*TEMP)
      RM=1.7768*1.672649D-24
      PI=3.141592654
      PF=((RM*(1.380662D-16)*TEMP*2.*PI/(6.626176D-27)**2)**1.5)
C
C     DELTA V = 0 BELOW
C
      DO 50 L=1,8
         STOKE1=EB(1,L+1)-EB(1,L)
C
C     FREQUENCY SHIFT FOR DELTA J=+1,DELTA V=0
C
         IF (L.GT.6) GO TO 40
         STOKE3=EB(1,L+3)-EB(1,L)
C
C     STOKE3- FREQUENCY SHIFT FOR DELTA J=+3, DELTA V=0
C
         IF (L.GT.4) GO TO 30
         STOKE5=EB(1,L+5)-EB(1,L)
         STOKI=A(L,L+5)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE5,STOKI)
         STOKIP=A(L,L+5)*DEXP(-ALFA*EB(1,L+5))/PF
         CALL PROFILE (-STOKE5,STOKIP)
   30    STOKI=A(L,L+3)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE3,STOKI)
         STOKIP=A(L,L+3)*DEXP(-ALFA*EB(1,L+3))/PF
         CALL PROFILE (-STOKE3,STOKIP)
   40    STOKI=A(L,L+1)*DEXP(-ALFA*EB(1,L))/PF
         CALL PROFILE (STOKE1,STOKI)
         STOKIP=A(L,L+1)*DEXP(-ALFA*EB(1,L+1))/PF
         CALL PROFILE (-STOKE1,STOKIP)
   50 CONTINUE
C
C     DELTA V=1 BELOW
C
      AV1=0.45814D-43
      STOKE=14.5442
      STOKI=AV1*DEXP(-ALFA*EB(1,5))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=.32581D-43
      STOKE=8.5755
      STOKI=AV1*DEXP(-ALFA*EB(1,6))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,1))/PF
      CALL PROFILE (-STOKE,STOKIP)
      AV1=.41005D-43
      STOKE=2.4395
      STOKI=AV1*DEXP(-ALFA*EB(1,7))/PF
      CALL PROFILE (STOKE,STOKI)
      STOKIP=AV1*DEXP(-ALFA*EB(2,2))/PF
      CALL PROFILE (-STOKE,STOKIP)
      DO 60 N=1,NSRIUP
   60 RSIBB(N)=RSIBB(N)/TWOPIC/SLIT
      NSOL=NSRI+1
      DO 70 I=1,NSOL
      K=(NSRI+1)+(I-1)
   70 RSI(I)=RSIBB(K)
C
C     RSI - BOUND-BOUND CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS
C
      RETURN
      END

      FUNCTION BGAUS (FNU,G0,DELTA,OMEGA0,TEMP)
C
C     THIS IS DESYMMETRIZED GAUSSIAN PROFILE, G0 IS A ZEROTH MOMENT
C     FNU,DELTA AND OMEGA0 ARE IN CM-1
C     G0 IN CGS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA BKW,TWOPI/0.6950304256D0,6.2831853D0/
      D=G0/(DELTA*DSQRT(TWOPI))
      DESYM=2./(1. + DEXP(-FNU/(BKW*TEMP)))
      FEXP=(FNU-OMEGA0)**2/(2.*DELTA**2)
      FEXP=DMIN1(FEXP,300.D0)
      BGAUS=D*DESYM*DEXP(-FEXP)
      RETURN
      END

