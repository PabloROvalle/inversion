	SUBROUTINE MATINV(A,N,NDIM,AINV)
C       $Id: matinv.f,v 1.1.1.1 2000/08/17 09:27:01 irwin Exp $
C-----------------------------------------------------------------------
C_TITLE: MATINV: matrix inversion routine
C
C_ARGS:  A id NxN matrix stored in a NDIMxNDIM array. AINV id also a NDIMxNDIM
C        array.
C
C_KEYS:
C
C_DESCR: on exit AINV contains the NxN inverse of A. A is destroyed.
C        The routine used Numerical Recipies subroutines and is copied almost
C        exactly from the book.
C
C_FILES: none
C
C_CALLS: none
C
C_BUGS:
C
C_HIST:  21jan88 SBC ORIGINAL VERSION
C        15apr93 PGJI modified
C
C_END:
C-----------------------------------------------------------------------
        INTEGER IDIM
	PARAMETER(IDIM=1024)
C	uses an internally declared array to save the bother of remembering
C	how to use Numerical Recipies routines
	INTEGER I,J,N,NDIM
        REAL D
	DIMENSION A(NDIM,NDIM),AINV(NDIM,NDIM)
        DIMENSION INDX(IDIM)
	IF(NDIM.GT.IDIM)THEN
	  WRITE(6,10)
10	  FORMAT(' ERROR - internal array in MATINV is too small')
	  STOP
	  END IF
	DO 20 J=1,N
 	  DO 30 I=1,N
	   AINV(I,J)=0.
30	  CONTINUE
	 AINV(J,J)=1.
20	CONTINUE
C
	CALL LUDCMP(A,N,NDIM,INDX,D)
	DO 40 J=1,N
 	  CALL LUBKSB(A,N,NDIM,INDX,AINV(1,J))
40	CONTINUE
	RETURN
	END
C-----------------------------------------------------------------------
