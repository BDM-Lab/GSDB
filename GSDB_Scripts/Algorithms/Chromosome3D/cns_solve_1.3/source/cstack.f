C STACK routines for CHARACTER*4 arrays
C======================================================================
      INTEGER FUNCTION CALLST(NUMWRD)
C
C stack allocation routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C begin
      INCLUDE 'cns.inc'
      INCLUDE 'cstack.inc'
      INTEGER NUMWRD
C begin
      IF (NUMWRD.LT.0) THEN
      CALL WRNDIE(-4,'<CALLST>','Allocation less than zero')
      END IF
      CALLST=CLSTUS+1
      CLSTUS=CLSTUS+NUMWRD
      CMAXUS=MAX(CMAXUS,CLSTUS)
      IF (CLSTUS.GT.CSTKSI) THEN
      CALL WRNDIE(-4,'<CALLST>','C-Stack overflow, increase CSTKSI')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE CFREST(NUMWRD)
C
C stack free-up routine
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cstack.inc'
      INTEGER NUMWRD
C begin
      IF (NUMWRD .LT. 0) THEN
      CALL WRNDIE(-4,'<CFREST>','Number of words less than zero')
      END IF
      CLSTUS=CLSTUS-NUMWRD
      IF(CLSTUS.LT.0) CALL WRNDIE(-4,'<CFREST>','Stack underflow')
      RETURN
      END
C======================================================================
      SUBROUTINE CINIST(QSET)
C
C stack initialization
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cstack.inc'
      LOGICAL QSET
C local
      INTEGER I
C begin
      CLSTUS=0
      CMAXUS=0
      IF (QSET) THEN
      DO I=1,CSTKSI
      CSTACK(I)='STCK'
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE PRISTK
C
C PRINTS INFO ABOUT STACK
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'timer.inc'
C begin
      IF ( WRNLEV .GE.15.OR.CLSTUS.NE.0 ) THEN
      WRITE(6,'(A,I9,A,I9,A,I9)')
     &' CSTACK: size=',CSTKSI,' used=',CMAXUS,' current=',CLSTUS
      END IF
C
      RETURN
      END
C
