*
* $Id: r2dp.F,v 1.1.1.1 1996/04/01 15:01:56 mclareni Exp $
*
* $Log: r2dp.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:56  mclareni
* Mathlib gen
*
*
#include "pilot.h"
#if defined(CERNLIB_DOUBLE)
      FUNCTION C309R2(X,ETA,ZL,P,EPS,LIMIT,KIND,ERR,NITS,
     1                FPMAX,ACC8,ACC16)
C
C *** evaluate the HYPERGEOMETRIC FUNCTION
C                                        i
C            F (AA;BB; Z) = SUM  (AA)   Z / ( (BB)  i! )
C           1 1              i       i            i
C
C     to accuracy EPS with at most LIMIT terms.
C  If KIND = 0 : using extended precision but real arithmetic only,
C            1 : using normal precision in complex arithmetic,
C   or       2 : using normal complex arithmetic, but with WDIGAM factor
C
C  where AA, BB, and Z are defined below
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMPLEX*16 X,ETA,ZL,P,AA,BB,Z,C309R2,WDIGAM
      COMPLEX*16 DD,G,F,AI,BI,T
#if (!defined(CERNLIB_UNIX))&&(!defined(CERNLIB_QMALPH))
      REAL*16 AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN
#endif
#if defined(CERNLIB_UNIX)||defined(CERNLIB_QMALPH)
      DOUBLE PRECISION AR,BR,GR,GI,DR,DI,TR,TI,UR,UI,FI,FI1,DEN
#endif

      PARAMETER(TBBB = 3D0/2D0)

#if defined(CERNLIB_QF2C)
#include "defdr.inc"
#endif

      ABSC(AA)=ABS(DREAL(AA))+ABS(DIMAG(AA))
      NINTC(AA)=NINT(DREAL(AA))

      C309R2 = 0.

      AA=ZL+1-ETA*P
      BB=2*ZL+2
      Z=2*P*X
      IF(DREAL(BB) .LE. 0 .AND. ABS(BB-NINTC(BB)) .LT.
     1 SQRT(SQRT(ACC8))**3 .AND. DREAL(BB)+LIMIT .GE. TBBB) THEN
       NITS=-1
       RETURN
      END IF
      IF(LIMIT .LE. 0) THEN
       C309R2=0
       ERR=0
       NITS=1
       RETURN
      END IF
      TA=1
      RK=1
      IF(KIND .LE. 0 .AND. ABSC(Z)*ABSC(AA) .GT. ABSC(BB)*1.0) THEN
       DR=1
       DI=0
       GR=1
       GI=0
       AR=DREAL(AA)
       BR=DREAL(BB)
       FI=0
       DO 20 I = 2,LIMIT
       FI1=FI+1
       TR=BR*FI1
       TI=DIMAG(BB)*FI1
       DEN=1/(TR*TR+TI*TI)
       UR=(AR*TR+DIMAG(AA)*TI)*DEN
       UI=(DIMAG(AA)*TR-AR*TI)*DEN
       TR=UR*GR-UI*GI
       TI=UR*GI+UI*GR
       GR=DREAL(Z)*TR-DIMAG(Z)*TI
       GI=DREAL(Z)*TI+DIMAG(Z)*TR
       DR=DR+GR
       DI=DI+GI
       ERR=ABS(GR)+ABS(GI)
       IF(ERR .GT. FPMAX) GO TO 60
       RK=ABS(DR)+ABS(DI)
       TA=MAX(TA,RK)
       IF(ERR .LT. RK*EPS .OR. I .GE. 4 .AND. ERR .LT. ACC16) GO TO 30
       FI=FI1
       AR=AR+1
   20  BR=BR+1
C
   30  C309R2=DR+(0,1)*DI
       ERR=ACC16*TA/RK
      ELSE
C
C*    If REAL*16 arithmetic is not available, (or already using it!),
C*    then use KIND > 0
C
       G=1
       F=1
       IF(KIND .GE. 2) F=WDIGAM(AA)-WDIGAM(BB)-WDIGAM(G)
       DD=F
       DO 40 I = 2,LIMIT
       AI=AA+(I-2)
       BI=BB+(I-2)
       R=I-1
       G=G*Z*AI/(BI*R)
C
C                       multiply by (psi(a+r)-psi(b+r)-psi(1+r))
C
       IF(KIND .EQ. 2) F=F+1/AI-1/BI-1/R
       T=G*F
       DD=DD+T
       ERR=ABSC(T)
       IF(ERR .GT. FPMAX) GO TO 60
       RK=ABSC(DD)
       TA=MAX(TA,RK)
       IF(ERR .LT. RK*EPS .OR. ERR .LT. ACC8 .AND. I .GE. 4) GO TO 50
   40  CONTINUE

   50  ERR=ACC8*TA/RK
       C309R2=DD
      END IF
   60 NITS=I
      RETURN
      END
#endif
