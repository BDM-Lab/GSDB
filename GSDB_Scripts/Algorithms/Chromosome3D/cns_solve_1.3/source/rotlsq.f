      SUBROUTINE ROTLSQ(X1,Y1,Z1,NATOM1,X2,Y2,Z2,NATOM2,
     &      ATOMPR,NPAIR,LMASS,AMASS1,AMASS2,QWGHT,KWMULT,LNOROT)
C
C     THE PROGRAM ROTATES COORDINATE SET 2 RESULTING IN A 2 SUCH THAT
C     THE SUM OF THE SQUARE OF THE DISTANCE BETWEEN EACH COORDINATE IN 1
C     AND 2 IS A MINIMUM.
C     THE LEAST SQUARE MINIMIZATION IS DONE ONLY WITH RESPECT TO THE
C     ATOMS REFERRED TO IN THE PAIR ARRAY. THE ROTATION MATRIX THAT IS
C     CALCULATED IS APPLIED TO THE ENTIRE SECOND COORDINATE SET.
C
C     THIS ROUTINE WRITTEN BY B. BROOKS, ADAPTED FROM ACTA CRYST
C     (1976) A32,922 W. KABSCH
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION X1(*),Y1(*),Z1(*)
      INTEGER NATOM1
      DOUBLE PRECISION X2(*),Y2(*),Z2(*)
      INTEGER NATOM2, ATOMPR(2,*), NPAIR
      LOGICAL LMASS
      DOUBLE PRECISION    AMASS1(*),AMASS2(*)
      LOGICAL QWGHT
      DOUBLE PRECISION    KWMULT(*)
      LOGICAL LNOROT
C local
      INTEGER MASS
C begin
      MASS=ALLHP(IREAL8(NPAIR))
      CALL PKMASS(AMASS1,AMASS2,HEAP(MASS),ATOMPR,NPAIR,LMASS,
     1            QWGHT,KWMULT)
      CALL ROTLS1(X1,Y1,Z1,NATOM1,X2,Y2,Z2,NATOM2,ATOMPR,NPAIR,
     2     HEAP(MASS),.TRUE.,LNOROT)
      CALL FREHP(MASS,IREAL8(NPAIR))
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ROTLS1(XA,YA,ZA,NATOMA,XB,YB,ZB,NATOMB,ATOMPR,NPAIR,
     &                  BMASS,LPRINT,LNOROT)
C
C     THIS ROUTINE DOES THE ACTUAL ROTATION OF B TO MATCH A. THE NEW
C     ARRAY B IS RETURNED IN X.
C     THIS ROUTINE WRITTEN BY B. BROOKS , ADAPTED FROM ACTA CRYST
C     (1976) A32,922 W. KABSCH
C
      IMPLICIT NONE
C input/output
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      DOUBLE PRECISION XA(*),YA(*),ZA(*)
      INTEGER NATOMA
      DOUBLE PRECISION XB(*),YB(*),ZB(*)
      INTEGER NATOMB
      INTEGER NPAIR, ATOMPR(2,*)
      DOUBLE PRECISION BMASS(*)
      LOGICAL LPRINT,LNOROT
C local
      DOUBLE PRECISION R(9), U(9), ROT(3,3)
      DOUBLE PRECISION CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION CMXC, CMYC, CMZC
      DOUBLE PRECISION XI, YI, ZI, XJ, YJ, ZJ
      DOUBLE PRECISION TMASS, RMST, RMSV
      DOUBLE COMPLEX DBCOMP
      INTEGER I, K, KA, KB
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      CMXA=ZERO
      CMYA=ZERO
      CMZA=ZERO
      CMXB=ZERO
      CMYB=ZERO
      CMZB=ZERO
      TMASS=ZERO
      DO K=1,NPAIR
      KA=ATOMPR(1,K)
      KB=ATOMPR(2,K)
      IF (INITIA(KA,XA,YA,ZA).AND.INITIA(KB,XB,YB,ZB)) THEN
      CMXA=CMXA+XA(KA)*BMASS(K)
      CMYA=CMYA+YA(KA)*BMASS(K)
      CMZA=CMZA+ZA(KA)*BMASS(K)
      CMXB=CMXB+XB(KB)*BMASS(K)
      CMYB=CMYB+YB(KB)*BMASS(K)
      CMZB=CMZB+ZB(KB)*BMASS(K)
      TMASS=TMASS+BMASS(K)
      END IF
      END DO
      CMXA=CMXA/TMASS
      CMYA=CMYA/TMASS
      CMZA=CMZA/TMASS
      CMXB=CMXB/TMASS
      CMYB=CMYB/TMASS
      CMZB=CMZB/TMASS
      DO K=1,NATOMB
      IF (INITIA(K,XB,YB,ZB)) THEN
      XB(K)=XB(K)-CMXB
      YB(K)=YB(K)-CMYB
      ZB(K)=ZB(K)-CMZB
      END IF
      END DO
C
      IF (LNOROT) THEN
C
C     USE A UNIT ROTATION MATRIX. NO ROTATION IS SPECIFIED
C
      DO K=1,NATOMB
      IF (INITIA(K,XB,YB,ZB)) THEN
      XB(K)=XB(K)+CMXA
      YB(K)=YB(K)+CMYA
      ZB(K)=ZB(K)+CMZA
      END IF
      END DO
C
      IF (LPRINT) THEN
C
C get translation vector t for (x'=R*x)
      CMXC=CMXA-CMXB
      CMYC=CMYA-CMYB
      CMZC=CMZA-CMZB
      WRITE(6,'(/2A)')
     & ' Fitted coordinate set r'' related to original set r by',
     & ' r''=r + T'
      WRITE(6,'(A)')
     & ' This transformation is applied to all coordinates.'
      WRITE(6,'(A,3F10.4,A/)')
     & ' Translation vector T = (',CMXC,CMYC,CMZC,')'
C
C put info about translation into symbol table
      CALL DECLAR( 'X', 'DP', ' ', DBCOMP, CMXC )
      CALL DECLAR( 'Y', 'DP', ' ', DBCOMP, CMYC )
      CALL DECLAR( 'Z', 'DP', ' ', DBCOMP, CMZC )
      END IF
C
      ELSE
C
C     COMPUTE ROTATION MATRIX FROM LAGRANGIAN
C
      DO I=1,9
      R(I)=ZERO
      END DO
      DO K=1,NPAIR
      KA=ATOMPR(1,K)
      KB=ATOMPR(2,K)
      IF (INITIA(KA,XA,YA,ZA).AND.INITIA(KB,XB,YB,ZB)) THEN
      XI=XB(KB)*BMASS(K)
      YI=YB(KB)*BMASS(K)
      ZI=ZB(KB)*BMASS(K)
      XJ=XA(KA)-CMXA
      YJ=YA(KA)-CMYA
      ZJ=ZA(KA)-CMZA
      R(1)=R(1)+XI*XJ
      R(2)=R(2)+XI*YJ
      R(3)=R(3)+XI*ZJ
      R(4)=R(4)+YI*XJ
      R(5)=R(5)+YI*YJ
      R(6)=R(6)+YI*ZJ
      R(7)=R(7)+ZI*XJ
      R(8)=R(8)+ZI*YJ
      R(9)=R(9)+ZI*ZJ
      END IF
      END DO
C
      CALL FROTU(R,U)
C
      DO K=1,NATOMB
      IF (INITIA(K,XB,YB,ZB)) THEN
      CMXC =U(1)*XB(K)+U(4)*YB(K)+U(7)*ZB(K)+CMXA
      CMYC =U(2)*XB(K)+U(5)*YB(K)+U(8)*ZB(K)+CMYA
      ZB(K)=U(3)*XB(K)+U(6)*YB(K)+U(9)*ZB(K)+CMZA
      XB(K)=CMXC
      YB(K)=CMYC
      END IF
      END DO
C
C print out info about rotation and translation
      IF (LPRINT) THEN
      ROT(1,1)=U(1)
      ROT(1,2)=U(4)
      ROT(1,3)=U(7)
      ROT(2,1)=U(2)
      ROT(2,2)=U(5)
      ROT(2,3)=U(8)
      ROT(3,1)=U(3)
      ROT(3,2)=U(6)
      ROT(3,3)=U(9)
C
C get translation vector t for (x'=R*x + t)
      CMXC=CMXA-(ROT(1,1)*CMXB+ROT(1,2)*CMYB+ROT(1,3)*CMZB)
      CMYC=CMYA-(ROT(2,1)*CMXB+ROT(2,2)*CMYB+ROT(2,3)*CMZB)
      CMZC=CMZA-(ROT(3,1)*CMXB+ROT(3,2)*CMYB+ROT(3,3)*CMZB)
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(/2A)')
     & ' Fitted coordinate set r'' related to original set r by',
     & ' r''=R*r  + T'
      WRITE(6,'(A)')
     & ' This transformation is applied to all coordinates.'
      WRITE(6,'(A,3F10.4,A)')
     & ' Translation vector T = (',CMXC,CMYC,CMZC,')'
C
      CALL MATPRI(ROT)
      END IF
C
C put info about translation into symbol table
      CALL DECLAR( 'X', 'DP', ' ', DBCOMP, CMXC )
      CALL DECLAR( 'Y', 'DP', ' ', DBCOMP, CMYC )
      CALL DECLAR( 'Z', 'DP', ' ', DBCOMP, CMZC )
C
C put info about rotation into symbol table
      CALL MATDCL(ROT)
      END IF
      END IF
C
      IF (LPRINT) THEN
      RMST=ZERO
      DO K=1,NPAIR
      KA=ATOMPR(1,K)
      KB=ATOMPR(2,K)
      IF (INITIA(KA,XA,YA,ZA) .AND. INITIA(KB,XB,YB,ZB)) THEN
      RMST=RMST+BMASS(K)*((XB(KB)-XA(KA))**2+(YB(KB)-YA(KA))**2+
     1          (ZB(KB)-ZA(KA))**2)
      END IF
      END DO
      RMSV=SQRT(RMST/TMASS)
      WRITE(6,'(2A,F10.4)')
     & ' R.m.s. diff. between fitted set and comp. set for ',
     & 'selected atoms =',RMSV
      END IF
      CALL DECLAR( 'RMSD', 'DP', ' ', DBCOMP, RMSV )
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, RMSV )
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FROTU(R,U)
C
C     THIS ROUTINE SOLVES THE CONSTRAINED MINIMIZATION EQUATION
C     USING LAGRANGE MULTIPLIERS.
C       BERNARD R. BROOKS
C
      IMPLICIT NONE
C input/output
      DOUBLE PRECISION R(3,3), U(3,3)
C local
      DOUBLE PRECISION W(3,3),A(3,3),B(3,3),SCR(24)
      DOUBLE PRECISION DET, TRACE, EVS
      INTEGER JP, JQ, KP, KQ
      INTEGER I, I1, I2, J, K, IPT
C parameter
      DOUBLE PRECISION ZERO, ONE, SMALL, THREE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, SMALL=1.0D-6, THREE=3.0D0)
C begin
      DET=ZERO
      DO I=1,3
      I1=I+1
      IF(I1.GT.3) I1=I1-3
      I2=I+2
      IF(I2.GT.3) I2=I2-3
      DET=DET+R(I,1)*(R(I1,2)*R(I2,3)-R(I2,2)*R(I1,3))
      END DO
C
      IPT=0
      DO I=1,3
      DO J=1,3
      W(I,J)=ZERO
      DO K=1,3
      W(I,J)=W(I,J)+R(J,K)*R(I,K)
      END DO
      END DO
      END DO
C
      TRACE=W(1,1)+W(2,2)+W(3,3)
      IF (TRACE.LT.THREE*SMALL) THEN
      DO I=1,3
      DO J=1,3
      U(I,J)=ZERO
      END DO
      U(I,I)=ONE
      END DO
      RETURN
      END IF
C
      CALL DIAGSQ(3,3,W,A,SCR)
C
CCC      CALL DIAGQ(3,3,W,A,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),
CCC     1  SCR(16),SCR(19),SCR(22),0)
C
C
      DO I=1,3
      SCR(I)=SQRT(ABS(SCR(I)))
      IF(SCR(I).LT.SMALL) SCR(I)=SMALL
      END DO
C
      IF(DET.LT.ZERO) SCR(1)=-SCR(1)
C
      DO J=1,3
      EVS=SCR(J)
      DO I=1,3
      B(I,J)=ZERO
      DO K=1,3
      B(I,J)=B(I,J)+R(K,I)*A(K,J)/EVS
      END DO
      END DO
      END DO
C
      DET=ZERO
      DO I=1,3
      I1=I+1
      IF(I1.GT.3) I1=I1-3
      I2=I+2
      IF(I2.GT.3) I2=I2-3
      DET=DET+A(I,1)*(A(I1,2)*A(I2,3)-A(I2,2)*A(I1,3))
      END DO
C
      DO J=1,3
      IF (ABS(SCR(J)).LE.SMALL) THEN
      JP=J+1
      JQ=J+2
      IF(JP.GT.3) JP=JP-3
      IF(JQ.GT.3) JQ=JQ-3
      DO K=1,3
      KP=K+1
      KQ=K+2
      IF(KP.GT.3) KP=KP-3
      IF(KQ.GT.3) KQ=KQ-3
      B(K,J)=B(KP,JP)*B(KQ,JQ)-B(KP,JQ)*B(KQ,JP)
      IF(DET.LT.ZERO) B(K,J)=-B(K,J)
      END DO
      END IF
      CALL NORMAL(B(1,J),3)
      END DO
C
      DO J=1,3
      DO I=1,3
      U(I,J)=ZERO
      DO K=1,3
      U(I,J)=U(I,J)+A(I,K)*B(J,K)
      END DO
      END DO
      END DO
C
      DO J=1,3
      CALL NORMAL(U(1,J),3)
      END DO
C
C     CHECK TO INSURE UNITY (AS OPPOSED TO ANTI-UNITARY)
      DET=ZERO
      DO I=1,3
      I1=I+1
      IF(I1.GT.3) I1=I1-3
      I2=I+2
      IF(I2.GT.3) I2=I2-3
      DET=DET+U(I,1)*(U(I1,2)*U(I2,3)-U(I2,2)*U(I1,3))
      END DO
C
      IF(ABS(DET-ONE).GT.SMALL) WRITE(6,105) DET
 105  FORMAT(/' ***** WARNING ***** FROM FROTU. ROTATION MATRIX IS ',
     1   'NOT UNITARY.'/,'  THE DETERMINANT IS',F14.8/)
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE PKMASS(AMASS1,AMASS2,BMASS,ATOMPR,NPAIR,LMASS,
     &                  QWGHT,KWMULT)
C
C     THIS ROUTINE PACKS THE BMASS ARRAY WITH MASSES FROM SELECTED
C     ATOMS. 1.0 IS USED IF LMASS IS FALSE
C     IF QWGHT IS TRUE THE PRODUCT OF BMASS AND KWMULT IS FORMED AND
C     STORED IN BMASS
C
      IMPLICIT NONE
C input/output
      DOUBLE PRECISION AMASS1(*),AMASS2(*)
      DOUBLE PRECISION BMASS(*)
      INTEGER ATOMPR(2,*), NPAIR
      LOGICAL LMASS, QWGHT
      DOUBLE PRECISION KWMULT(*)
C local
      INTEGER K, KX, KY, NWARN
C begin
C
      IF (LMASS) THEN
      NWARN=0
      DO K=1,NPAIR
      KX=ATOMPR(1,K)
      KY=ATOMPR(2,K)
      BMASS(K)=AMASS1(KX)
      IF (AMASS1(KX).NE.AMASS2(KY)) THEN
      NWARN=NWARN+1
      BMASS(K)=SQRT(AMASS1(KX)*AMASS2(KY))
      END IF
      END DO
      IF(NWARN.GT.0) WRITE(6,45) NWARN
  45  FORMAT(/' *** WARNING *** MASSES DONT MATCH FOR THIS',
     1  ' HOMOLOGY FOR',I5,' ATOMS. RESULTS WILL USE GEOMETRIC MEAN.')
      ELSE
      DO K=1,NPAIR
      BMASS(K)=1.0D0
      END DO
      END IF
      IF (QWGHT) THEN
      DO K=1,NPAIR
      BMASS(K)=BMASS(K)*KWMULT(ATOMPR(1,K))
      END DO
      END IF
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ROTATE(XA,YA,ZA,NPAIR,XB,YB,ZB,U,LPRINT,NGROUP,IOP)
C
C     THE PROGRAM ROTATES COORDINATE SET 2 RESULTING IN A 2 SUCH THAT
C     THE SUM OF THE SQUARE OF THE DISTANCE BETWEEN EACH COORDINATE IN 1
C     AND 2 IS A MINIMUM.
C     THE LEAST SQUARE MINIMIZATION IS DONE ONLY WITH RESPECT TO THE
C     ATOMS REFERRED TO IN THE PAIR ARRAY. THE ROTATION MATRIX THAT IS
C     CALCULATED IS APPLIED TO THE ENTIRE SECOND COORDINATE SET.
C          BERNARD R. BROOKS
C
C     THIS ROUTINE DOES THE ACTUAL ROTATION OF B TO MATCH A. THE NEW
C     ARRAY B IS RETURNED IN X.
C     THIS ROUTINE WRITTEN BY B. BROOKS , ADAPTED FROM ACTA CRYST
C     (1976) A32,922 W. KABSCH
C
C Modification of B. Brooks' ROTLSQ routines for non-crystallographic
C symmetry restraints.  Only significant change is removal of option
C for doing true center-of-mass: all atoms passed for least squares
C superposition have equal weight (i.e. get simple average).
C Checks for matching atoms have already been done by non-crystallographic
C symmetry parsing routine.
C -- Bill Weis
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      DOUBLE PRECISION XA(*),YA(*),ZA(*)
      INTEGER NPAIR
      DOUBLE PRECISION XB(*),YB(*),ZB(*)
      INTEGER NGROUP,IOP
C local
      DOUBLE PRECISION R(9), U(9)
      DOUBLE PRECISION CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION CMXC, CMYC, CMZC
      DOUBLE PRECISION ICMXC, ICMYC, ICMZC
      DOUBLE PRECISION XI, YI, ZI, XJ, YJ, ZJ
      DOUBLE PRECISION RMST, RMSV, ROT(3,3), IROT(3,3)
      INTEGER I, K
      LOGICAL LPRINT
      CHARACTER*20 OPS
      INTEGER LENOPS
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      CMXA=ZERO
      CMYA=ZERO
      CMZA=ZERO
      CMXB=ZERO
      CMYB=ZERO
      CMZB=ZERO
      DO K=1,NPAIR
      CMXA=CMXA+XA(K)
      CMYA=CMYA+YA(K)
      CMZA=CMZA+ZA(K)
      CMXB=CMXB+XB(K)
      CMYB=CMYB+YB(K)
      CMZB=CMZB+ZB(K)
      END DO
      CMXA=CMXA/NPAIR
      CMYA=CMYA/NPAIR
      CMZA=CMZA/NPAIR
      CMXB=CMXB/NPAIR
      CMYB=CMYB/NPAIR
      CMZB=CMZB/NPAIR
      DO K=1,NPAIR
      XB(K)=XB(K)-CMXB
      YB(K)=YB(K)-CMYB
      ZB(K)=ZB(K)-CMZB
      END DO
C
C     COMPUTE ROTATION MATRIX FROM LAGRANGIAN
C
      DO I=1,9
      R(I)=ZERO
      END DO
      DO K=1,NPAIR
      XI=XB(K)
      YI=YB(K)
      ZI=ZB(K)
      XJ=XA(K)-CMXA
      YJ=YA(K)-CMYA
      ZJ=ZA(K)-CMZA
      R(1)=R(1)+XI*XJ
      R(2)=R(2)+XI*YJ
      R(3)=R(3)+XI*ZJ
      R(4)=R(4)+YI*XJ
      R(5)=R(5)+YI*YJ
      R(6)=R(6)+YI*ZJ
      R(7)=R(7)+ZI*XJ
      R(8)=R(8)+ZI*YJ
      R(9)=R(9)+ZI*ZJ
      END DO
C
      CALL FROTU(R,U)
C
      DO K=1,NPAIR
      CMXC =U(1)*XB(K)+U(4)*YB(K)+U(7)*ZB(K)+CMXA
      CMYC =U(2)*XB(K)+U(5)*YB(K)+U(8)*ZB(K)+CMYA
      ZB(K)=U(3)*XB(K)+U(6)*YB(K)+U(9)*ZB(K)+CMZA
      XB(K)=CMXC
      YB(K)=CMYC
      END DO
C
C
C print out info about rotation and translation
      IF (LPRINT) THEN
      ROT(1,1)=U(1)
      ROT(1,2)=U(4)
      ROT(1,3)=U(7)
      ROT(2,1)=U(2)
      ROT(2,2)=U(5)
      ROT(2,3)=U(8)
      ROT(3,1)=U(3)
      ROT(3,2)=U(6)
      ROT(3,3)=U(9)
C
C get translation vector t for (x'=R*x + t)
      CMXC=CMXA-(ROT(1,1)*CMXB+ROT(1,2)*CMYB+ROT(1,3)*CMZB)
      CMYC=CMYA-(ROT(2,1)*CMXB+ROT(2,2)*CMYB+ROT(2,3)*CMZB)
      CMZC=CMZA-(ROT(3,1)*CMXB+ROT(3,2)*CMYB+ROT(3,3)*CMZB)
C
      IROT(1,1)=ROT(1,1)
      IROT(1,2)=ROT(2,1)
      IROT(1,3)=ROT(3,1)
      IROT(2,1)=ROT(1,2)
      IROT(2,2)=ROT(2,2)
      IROT(2,3)=ROT(3,2)
      IROT(3,1)=ROT(1,3)
      IROT(3,2)=ROT(2,3)
      IROT(3,3)=ROT(3,3)
      ICMXC=-(IROT(1,1)*CMXC)-(IROT(1,2)*CMYC)-(IROT(1,3)*CMZC)
      ICMYC=-(IROT(2,1)*CMXC)-(IROT(2,2)*CMYC)-(IROT(2,3)*CMZC)
      ICMZC=-(IROT(3,1)*CMXC)-(IROT(3,2)*CMYC)-(IROT(3,3)*CMZC)
C
C report transformation from reference to molecule
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(/2A)')
     & ' Reference set r is superimposed on',
     & ' equivalence set by: rnew=R*r + T'
      WRITE(6,'(A,3F10.4,A)')
     & ' Translation vector T = (',ICMXC,ICMYC,ICMZC,')'
      CALL MATPRI(IROT)
      END IF
C
C report transformation from molecule to reference
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(/2A)')
     & ' Equivalence set r'' is superimposed on',
     & ' reference set by: r''new=R''*r'' + T'''
      WRITE(6,'(A,3F10.4,A)')
     & ' Translation vector T = (',CMXC,CMYC,CMZC,')'
      CALL MATPRI(ROT)
      END IF
      END IF
C
      RMST=ZERO
      DO K=1,NPAIR
      RMST=RMST+((XB(K)-XA(K))**2+(YB(K)-YA(K))**2+(ZB(K)-ZA(K))**2)
      END DO
      RMSV=SQRT(RMST/NPAIR)
      IF (LPRINT) THEN
      WRITE(6,'(2A,F10.4)')
     & ' R.m.s. diff. between molecule and reference for ',
     & 'selected atoms =',RMSV
      END IF
C
C declare rotation matrix, translation vector and RMS difference
C 22-MAR-95 JSJ
C modified to add group number to symbol - PDA 26-2-97
      IF (LPRINT) THEN
C symbols for transformation from reference to molecule
      CALL ENCODI(NGROUP,WDT,WDTMAX,WDTLEN)
      OPS='ROT_'//WDT(1:WDTLEN)
      LENOPS=4+WDTLEN
      CALL OPSDCL(OPS,LENOPS,IROT,ICMXC,ICMYC,ICMZC,RMSV,IOP)
C symbols for transformation molecule to reference
      CALL ENCODI(NGROUP,WDT,WDTMAX,WDTLEN)
      OPS='ROTINV_'//WDT(1:WDTLEN)
      LENOPS=7+WDTLEN
      CALL OPSDCL(OPS,LENOPS,ROT,CMXC,CMYC,CMZC,RMSV,IOP)
      END IF
C
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ORINTC(NAT,X,Y,Z,XCOMP,YCOMP,ZCOMP,AMASS,LMASS,LRMS,
     &     ATOMIN,ISLCT,LNORO)
C
C  THIS ROUTINE IS CALLED BY THE ORIENT COMMAND AND ROTATES THE
C  MOLECULE SO THAT IT SITS ON THE ORIGIN WITH NO OFF-DIAGONAL
C  MOMENTS. LMASS CAUSES A MASS WEIGHTING TO BE DONE.
C  LRMS = .TRUE. WILL CAUSE THE ROTATION TO BE WRT THE OTHER SET.
C
C     Authors: Bernie Brooks
C              Axel Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER NAT
      DOUBLE PRECISION X(*),Y(*),Z(*),XCOMP(*),YCOMP(*),ZCOMP(*)
      DOUBLE PRECISION AMASS(*)
      LOGICAL LMASS, LRMS
      INTEGER ATOMIN(2,*)
      INTEGER ISLCT(*)
      LOGICAL LNORO
C local
      INTEGER N, NPR, I, J, IPT
      DOUBLE PRECISION ROT(3,3), CMXC, CMYC, CMZC
      DOUBLE PRECISION U(9),SCR(24),AMOM(3,3),XC,YC,ZC,AMASST,AM,DET
      DOUBLE PRECISION XX,XY,XZ,YY,YZ,ZZ,XN,YN,ZN
      DOUBLE COMPLEX DBCOMP
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
      N=NAT
C
      IF (LRMS) THEN
      NPR=0
      DO I=1,N
      IF (ISLCT(I).EQ.1) THEN
      NPR=NPR+1
      ATOMIN(1,NPR)=I
      ATOMIN(2,NPR)=I
      END IF
      END DO
C
      CALL ROTLSQ(XCOMP,YCOMP,ZCOMP,N,X,Y,Z,N,ATOMIN,NPR,
     1           LMASS,AMASS,AMASS,.FALSE.,ZERO,LNORO)
C
      ELSE
C
      XC=ZERO
      YC=ZERO
      ZC=ZERO
      AMASST=ZERO
      DO I=1,N
      IF (ISLCT(I).EQ.1) THEN
      IF (LMASS) THEN
      AM=AMASS(I)
      ELSE
      AM=ONE
      END IF
      XC=XC+X(I)*AM
      YC=YC+Y(I)*AM
      ZC=ZC+Z(I)*AM
      AMASST=AMASST+AM
      END IF
      END DO
      XC=XC/AMASST
      YC=YC/AMASST
      ZC=ZC/AMASST
      DO I=1,N
C
C do it for all atoms !!!
      IF (INITIA(I,X,Y,Z)) THEN
      X(I)=X(I)-XC
      Y(I)=Y(I)-YC
      Z(I)=Z(I)-ZC
      END IF
      END DO
C
      IF (.NOT.(LNORO) )THEN
      XX=ZERO
      XY=ZERO
      XZ=ZERO
      YY=ZERO
      YZ=ZERO
      ZZ=ZERO
      DO I=1,N
      IF (ISLCT(I).EQ.1) THEN
      IF (LMASS) THEN
      AM=AMASS(I)
      ELSE
      AM=ONE
      END IF
      XX=XX+X(I)*X(I)*AM
      XY=XY+X(I)*Y(I)*AM
      XZ=XZ+X(I)*Z(I)*AM
      YY=YY+Y(I)*Y(I)*AM
      YZ=YZ+Y(I)*Z(I)*AM
      ZZ=ZZ+Z(I)*Z(I)*AM
      END IF
      END DO
C
      AMOM(1,1)=ZZ
      AMOM(1,2)=YZ
      AMOM(2,1)=YZ
      AMOM(1,3)=XZ
      AMOM(3,1)=XZ
      AMOM(2,2)=YY
      AMOM(2,3)=XY
      AMOM(3,2)=XY
      AMOM(3,3)=XX
      WRITE(6,105) ZZ, YZ, XZ, YY, XY, XX
 105  FORMAT(' MOMENTS'/3F16.8,/16X,2F16.8,/32X,F16.8/)
C
      CALL DIAGSQ(3,3,AMOM,U,SCR)
CCC      CALL DIAGQ(3,3,AMOM,U,SCR(4),SCR(7),SCR(10),SCR(13),SCR(1),
CCC     1  SCR(16),SCR(19),SCR(22),0)
C
      DO I=1,3
      DET=U(I)
      U(I)=U(I+6)
      U(I+6)=DET
      END DO
      DO I=1,3
      IPT=(I-1)*3
      DET=U(IPT+1)
      U(IPT+1)=U(IPT+3)
      U(IPT+3)=DET
      IF (U(IPT+I).LT.ZERO) THEN
      DO J=1,3
      IPT=IPT+1
      U(IPT)=-U(IPT)
      END DO
      END IF
      END DO
      DET=U(1)*(U(5)*U(9)-U(6)*U(8))+U(2)*(U(6)*U(7)-U(4)*U(9))+
     1  U(3)*(U(4)*U(8)-U(5)*U(7))
      IF (DET.LT.ZERO) THEN
      U(7)=-U(7)
      U(8)=-U(8)
      U(9)=-U(9)
      END IF
      IF(ABS(DET-ONE).GT.RSMALL) WRITE(6,203) DET
 203  FORMAT(/' ***** WARNING ***** FROM ORINTC. ROTATION MATRIX IS ',
     1   'NOT UNITARY.'/,' DETERMINANT=',F14.8/)
C
C print out info about rotation and translation
      ROT(1,1)=U(1)
      ROT(1,2)=U(2)
      ROT(1,3)=U(3)
      ROT(2,1)=U(4)
      ROT(2,2)=U(5)
      ROT(2,3)=U(6)
      ROT(3,1)=U(7)
      ROT(3,2)=U(8)
      ROT(3,3)=U(9)
C
C get translation vector t for (x'=R*x + t)
      CMXC=-(ROT(1,1)*XC+ROT(1,2)*YC+ROT(1,3)*ZC)
      CMYC=-(ROT(2,1)*XC+ROT(2,2)*YC+ROT(2,3)*ZC)
      CMZC=-(ROT(3,1)*XC+ROT(3,2)*YC+ROT(3,3)*ZC)
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(/2A)')
     & ' Oriented coordinate set r'' related to original set r by',
     & ' r''=R*r  + T'
      WRITE(6,'(A,3F10.4,A)')
     & ' Translation vector T = (',CMXC,CMYC,CMZC,')'
      CALL MATPRI(ROT)
      END IF
C
C put info about translation into symbol table
      CALL DECLAR( 'X', 'DP', ' ', DBCOMP, CMXC )
      CALL DECLAR( 'Y', 'DP', ' ', DBCOMP, CMYC )
      CALL DECLAR( 'Z', 'DP', ' ', DBCOMP, CMZC )
C
C put info about rotation into symbol table
      CALL MATDCL(ROT)
C
      DO I=1,N
C
C do it for all atoms !!!
      IF (INITIA(I,X,Y,Z)) THEN
      XN=U(1)*X(I)+U(2)*Y(I)+U(3)*Z(I)
      YN=U(4)*X(I)+U(5)*Y(I)+U(6)*Z(I)
      ZN=U(7)*X(I)+U(8)*Y(I)+U(9)*Z(I)
      X(I)=XN
      Y(I)=YN
      Z(I)=ZN
      END IF
      END DO
      END IF
C
      END IF
C
      RETURN
      END
C================================================================
      SUBROUTINE NORMAL(V,N)
C  NORMALIZES VECTOR V OF LENGTH N
C   B. R. BROOKS
C
      IMPLICIT NONE
C input/output
      INTEGER N
      DOUBLE PRECISION V(N)
C local
      DOUBLE PRECISION C
      INTEGER I
C begin
      C=0.0D0
      DO 10 I=1,N
  10   C=C+V(I)*V(I)
      IF(C.GT.1.0D-12) GOTO 15
      WRITE(6,13) C
  13  FORMAT(' %NORMAL-ERR: TRYING TO NORMALIZE A ZERO VECTOR',
     1    ' NORM=',E12.4/' IT WILL BE ZEROED.')
      DO 14 I=1,N
  14  V(I)=0
      RETURN
C
  15   C=1.0D0/SQRT(C)
      DO 20 I=1,N
  20   V(I)=V(I)*C
      RETURN
      END
