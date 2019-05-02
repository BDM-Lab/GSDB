      SUBROUTINE TESTCNS
C
C CNS tests for debugging purposes
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
C begin
      CALL NEXTWD('TEST-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-test')
C
      ELSE IF (WD(1:4).EQ.'FIRS') THEN
      CALL TESTFD(NATOM,X,Y,Z)
      ELSE IF (WD(1:4).EQ.'HEAP') THEN
      CALL TSTHEP
      ELSE
      CALL DSPERR('TEST','unknown qualifier')
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE TESTFD(NATOM,X,Y,Z)
C
C Routine to test analytical first derivatives by
C finite differences along all cartesian dimensions of all
C selected atoms.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*)
C local
      INTEGER IDXFD, IDXAD, IXOLD, IYOLD, IZOLD, ISLCT
C begin
      IDXFD=ALLHP(IREAL8(3*NATOM))
      IDXAD=ALLHP(IREAL8(3*NATOM))
      IXOLD=ALLHP(IREAL8(NATOM))
      IYOLD=ALLHP(IREAL8(NATOM))
      IZOLD=ALLHP(IREAL8(NATOM))
      ISLCT=ALLHP(INTEG4(NATOM))
      CALL TESTF2(NATOM,X,Y,Z,HEAP(IXOLD),HEAP(IYOLD),
     &       HEAP(IZOLD),HEAP(IDXAD),HEAP(IDXFD),HEAP(ISLCT))
      CALL FREHP(ISLCT,INTEG4(NATOM))
      CALL FREHP(IZOLD,IREAL8(NATOM))
      CALL FREHP(IYOLD,IREAL8(NATOM))
      CALL FREHP(IXOLD,IREAL8(NATOM))
      CALL FREHP(IDXAD,IREAL8(3*NATOM))
      CALL FREHP(IDXFD,IREAL8(3*NATOM))
      RETURN
      END
C======================================================================
      SUBROUTINE TESTF2(NATOM,X,Y,Z,XOLD,YOLD,ZOLD,
     &     DXAD,DXFD,ISLCT)
C
C see TESTFD above
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'funct.inc'
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XOLD(*)
      DOUBLE PRECISION YOLD(*), ZOLD(*)
      DOUBLE PRECISION DXAD(*),DXFD(*)
      INTEGER ISLCT(*)
C local
      INTEGER NAT3, IAT, IAT3, I, NSELCT, NOK, II
      LOGICAL QMASS
      DOUBLE PRECISION STEP, TOL1
      DOUBLE PRECISION CR(3), DEL, VALMIN, DIF, VAL
      CHARACTER*4 SID, RID, RESN, AC
      CHARACTER*2 SNAME(3)
C parameter
      DOUBLE PRECISION ONE, ZERO, TWO, P5, P1
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, ONE=1.0D0, P5=0.005D0)
      PARAMETER (P1=0.0001D0)
C begin
      SNAME(1)=' X'
      SNAME(2)=' Y'
      SNAME(3)=' Z'
C defaults
      STEP=P5
      TOL1=P1
      QMASS=.FALSE.
      CALL FILL4(ISLCT,NATOM,1)
      NSELCT=NATOM
C
C parsing
      CALL PUSEND('TEST-FIRST>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('TEST-FIRST>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-test-first')
C
      ELSE IF (WD(1:4).EQ.'STEP') THEN
      CALL NEXTF('STEP=',STEP)
      ELSE IF (WD(1:3).EQ.'TOL') THEN
      CALL NEXTF('TOLerance=',TOL1)
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL NEXTLO('MASS=',QMASS)
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(ISLCT,NSELCT,X,Y,Z,.TRUE.)
      ELSE
      CALL CHKEND('TEST-FIRST>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NSELCT.GT.0) THEN
C
C save current coordinates
      CALL COPYCO(NATOM,X,Y,Z,XOLD,YOLD,ZOLD)
C
C call energy
      CALL ENERGY
C
      VALMIN=ONE
      DEL=STEP
      NAT3=NATOM*3
C
      IAT3=0
      DO IAT=1,NATOM
      IF (ISLCT(IAT).EQ.1) THEN
      DXAD(IAT3+1)=DX(IAT)
      DXAD(IAT3+2)=DY(IAT)
      DXAD(IAT3+3)=DZ(IAT)
      IAT3=IAT3+3
      END IF
      END DO
C
      IAT3=0
      DO IAT=1,NATOM
      IF (ISLCT(IAT).EQ.1) THEN
      CR(1)=XOLD(IAT)
      CR(2)=YOLD(IAT)
      CR(3)=ZOLD(IAT)
      DO I=1,3
      CR(1)=ZERO
      CR(2)=ZERO
      CR(3)=ZERO
      CR(I)=DEL
      IAT3=IAT3+1
C
C     make finite diff approx in (IAT,I) direction (I: X<>1, Y<>2, Z<>3)
C
      DO II=1,NATOM
      X(II)=XOLD(II)
      Y(II)=YOLD(II)
      Z(II)=ZOLD(II)
      END DO
      X(IAT)=X(IAT)+CR(1)
      Y(IAT)=Y(IAT)+CR(2)
      Z(IAT)=Z(IAT)+CR(3)
C
C call energy
      CALL ENERGY
      DXFD(IAT3)=RENR(SSENER)
C
      DO II=1,NATOM
      X(II)=XOLD(II)
      Y(II)=YOLD(II)
      Z(II)=ZOLD(II)
      END DO
      X(IAT)=X(IAT)-CR(1)
      Y(IAT)=Y(IAT)-CR(2)
      Z(IAT)=Z(IAT)-CR(3)
C
C call energy
      CALL ENERGY
C
      DXFD(IAT3)=(DXFD(IAT3)-RENR(SSENER))/(TWO*DEL)
C
      END DO
      END IF
      END DO
C
      WRITE(6,1000) STEP, QMASS, TOL1
1000  FORMAT(' TESTFD: Parameters: STEP=',F10.5,'  MASSweighting=',L1,
     & /,' TESTFD: The following first derivatives',
     & ' differ by more than TOL=',F12.6,//,
     & '  DIM.',9X,'ATOM',17X,'ANALYTIC',8X,'FINITE-DIFF',5X,
     & 'DEVIATION')
      NOK=0
      IAT3=0
      DO IAT=1,NATOM
      IF (ISLCT(IAT).EQ.1) THEN
      CALL ATOMID(IAT,SID,RID,RESN,AC)
      DO I=1,3
      IAT3=IAT3+1
      DIF=DXAD(IAT3)-DXFD(IAT3)
      VAL=ABS(DXAD(IAT3))
      IF(VAL.LT.VALMIN) VAL=VALMIN
      IF (ABS(DIF/VAL).GE.TOL1) THEN
      WRITE(6,4000) IAT,SNAME(I),SID,RID,RESN,AC,
     &                DXAD(IAT3),DXFD(IAT3),ABS(DIF)
4000  FORMAT(1X,I4,A,' (',4(1X,A),')',3F15.6)
      ELSE
      NOK=NOK+1
      END IF
      END DO
      END IF
      END DO
      WRITE(6,5000) NOK
5000  FORMAT(
     & ' TESTFD: A total of',I6,' elements were within the tolerance')
C
C retrieve original coordinates
      CALL COPYCO(NATOM,XOLD,YOLD,ZOLD,X,Y,Z)
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE TSTHEP
C
C test heap allocation and deallocation
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER LEN, BASE, OLEV
C begin
C
C turn on debugging message mode
      OLEV=WRNLEV
      WRNLEV=15
C command parsing
      CALL PUSEND('TESTHP>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('TESTHP>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-test-heap')
C
      ELSE IF (WD(1:4).EQ.'ALLO') THEN
      CALL NEXTI('ALLOcate=',LEN)
      BASE=ALLHP(LEN)
      ELSE IF (WD(1:4).EQ.'FREE') THEN
      CALL NEXTI('FREE (base)=',BASE)
      CALL NEXTI('FREE (length)=',LEN)
      CALL FREHP(BASE,LEN)
      ELSE
      CALL CHKEND('CONNECT>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C restore old warning level
      WRNLEV=OLEV
      RETURN
      END
