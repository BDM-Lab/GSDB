      SUBROUTINE EBOND(EB,EBV,IB,JB,INTERE,NBOND,ACBC,ACBB,QACBC,QACBB,
     &                 ANALYS,WGHT,VWGHT,IAC,SEGID,RESID,TYPE)
C
C Calculates bond energy.
C Front-end for actual energy routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'update.inc'
      INCLUDE 'param.inc'
      DOUBLE PRECISION EB, EBV
      DOUBLE PRECISION WGHT, VWGHT
      INTEGER IB(*),JB(*), NBOND, INTERE(*)
      DOUBLE PRECISION ACBC(*), ACBB(*)
      LOGICAL QACBC(*), QACBB(*)
      CHARACTER*4 ANALYS, IAC(*), SEGID(*), RESID(*), TYPE(*)
C
C begin
C
C fill the unkown ACBC and ACBB values with type-based parameters
      IF (UPBOND) THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,IAC,
     &            SEGID,RESID,TYPE)
      UPBOND=.FALSE.
      END IF
C
      CALL EBOND2(EB,EBV,IB,JB,INTERE,NBOND,ACBC,ACBB,ANALYS,
     &                  WGHT,VWGHT)
      RETURN
      END
C
      SUBROUTINE EBOND2(EB,EBV,IB,JB,INTERE,NBOND,ACBC,ACBB,ANALYS,
     &                  WGHT,VWGHT)
C
C Calculates bond energy force field in a vectorizable way.
C Front-end for actual energy routine
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'pick.inc'
      DOUBLE PRECISION EB, EBV
      DOUBLE PRECISION WGHT, VWGHT
      INTEGER IB(*),JB(*), NBOND, INTERE(*)
      DOUBLE PRECISION ACBC(*), ACBB(*)
      CHARACTER*4 ANALYS
C local
      INTEGER IBOND, ILOC, IL, II
      DOUBLE PRECISION  DB, RX, RY, RZ, DF, S, TWO, WTWO
      PARAMETER (TWO=2.0D0)
C====>VECDIM determines the optimal length of the vector loops
C====>on the CRAY it is 64
      INTEGER VECDIM
      PARAMETER (VECDIM=64)
      INTEGER LIND(VECDIM)
      DOUBLE PRECISION XI(VECDIM), YI(VECDIM), ZI(VECDIM)
      DOUBLE PRECISION XJ(VECDIM), YJ(VECDIM), ZJ(VECDIM)
      DOUBLE PRECISION LCBB(VECDIM), LCBC(VECDIM)
      DOUBLE PRECISION DXI(VECDIM), DYI(VECDIM), DZI(VECDIM)
      LOGICAL QANAL
C begin
C
      QANAL=.FALSE.
      IF (ANALYS.EQ.'ANAL') QANAL=.TRUE.
C
C initialization
      EB=0.0D0
      WTWO=TWO*WGHT
C
C loop over all active bonds
C make vector loops of length VECDIM
      IBOND=0
      DO WHILE (IBOND.LT.NBOND)
      ILOC=0
C
C gather the data into temporary arrays for following vector loop
      DO WHILE (ILOC.LT.VECDIM.AND.IBOND.LT.NBOND)
      IBOND=IBOND+1
      IF (QANAL.OR.ABS(INTERE(IB(IBOND))+INTERE(JB(IBOND))).LE.+1) THEN
      ILOC=ILOC+1
      LIND(ILOC)=IBOND
      END IF
      END DO
      DO IL=1,ILOC
      XI(IL)=X(IB(LIND(IL)))
      YI(IL)=Y(IB(LIND(IL)))
      ZI(IL)=Z(IB(LIND(IL)))
      XJ(IL)=X(JB(LIND(IL)))
      YJ(IL)=Y(JB(LIND(IL)))
      ZJ(IL)=Z(JB(LIND(IL)))
      LCBB(IL)=ACBB(LIND(IL))
      LCBC(IL)=ACBC(LIND(IL))
      END DO
C
C==================
C-BEGIN-VECTOR-LOOP
C==================
      DO IL=1,ILOC
      RX=XI(IL)-XJ(IL)
      RY=YI(IL)-YJ(IL)
      RZ=ZI(IL)-ZJ(IL)
      S=SQRT(MAX(RSMALL,RX*RX+RY*RY+RZ*RZ))
      DB=S-LCBB(IL)
      DF=LCBC(IL)*DB
      EB=EB+DF*DB
      DF=DF*WTWO/S
      DXI(IL)=RX*DF
      DYI(IL)=RY*DF
      DZI(IL)=RZ*DF
      END DO
C================
C-END-VECTOR-LOOP
C================
C
C now scatter the forces (non-vectorizable)
      DO IL=1,ILOC
      II=LIND(IL)
      DX(IB(II))=DX(IB(II))+DXI(IL)
      DY(IB(II))=DY(IB(II))+DYI(IL)
      DZ(IB(II))=DZ(IB(II))+DZI(IL)
      DX(JB(II))=DX(JB(II))-DXI(IL)
      DY(JB(II))=DY(JB(II))-DYI(IL)
      DZ(JB(II))=DZ(JB(II))-DZI(IL)
      END DO
C
      END DO
C
      EBV=EB*VWGHT
      EB=EB*WGHT
C
C do some analysis if appropriate
      IF (ANALYS.EQ.'ANAL') THEN
      PCDATA(PCGEOM)=S
      PCDATA(PCENER)=EB
      PCDATA(PCEQUI)=ACBB(1)
      PCDATA(PCCONS)=ACBC(1)
      PCDATA(PCDEVI)=DB
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,IAC,
     &                  SEGID,RESID,TYPE)
C
C Gets type-based parameters for bonds.
C
C Author: Axel T. Brunger
C =======================
C
C modified for atom-based parameters.
C Thomas Simonson, June 91.
C
      IMPLICIT NONE
C input/output
      INCLUDE 'timer.inc'
      INCLUDE 'param.inc'
      INTEGER NBOND, IB(*), JB(*)
      DOUBLE PRECISION ACBC(*), ACBB(*)
      LOGICAL QACBC(*), QACBB(*)
      CHARACTER*4 IAC(*), SEGID(*), RESID(*), TYPE(*)
C local
      INTEGER IFLAG, I, ICB
      INTEGER MARK
      PARAMETER (MARK=0)
C begin
C
      IFLAG=0
C
      DO I=1,NBOND
C
C ignore atom-based parameter entries
      IF (QACBC(I).OR.QACBB(I)) THEN
C
      CALL BINBON(ICB,IAC(IB(I)),IAC(JB(I)),NCB,KCB1,KCB2,MARK)
      IF (ICB.EQ.MARK) THEN
      IFLAG=IFLAG+1
      WRITE(6,'(A)')
     &' %CODBON-ERR: missing bond parameters %%%%%%%%%%%%%%%%%%%%%%%%%%'
      IF (QACBC(I)) WRITE(6,'(A)') '  bond energy constant missing.'
      IF (QACBB(I)) WRITE(6,'(A)') '  target bond length missing.'
      WRITE(6,'(9A)')
     & '  ATOM1: SEGId="',SEGID(IB(I)),'",  RESId="',RESID(IB(I)),
     & '",  NAME="',TYPE(IB(I)),'",  CHEMical="',IAC(IB(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM2: SEGId="',SEGID(JB(I)),'",  RESId="',RESID(JB(I)),
     & '",  NAME="',TYPE(JB(I)),'",  CHEMical="',IAC(JB(I)),'"'
      WRITE(6,'(A)')
     &' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ELSE
      IF (QACBC(I)) ACBC(I)=CBC(ICB)
      IF (QACBB(I)) ACBB(I)=CBB(ICB)
      END IF
      END IF
      END DO
      IF (IFLAG.GT.0) THEN
      CALL WRNDIE(-1,'CODBON','program will be aborted.')
      ELSE IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' CODBON: bond type-based parameters retrieved'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE BINBON(BINPAR,K1,K2,NCB,KCB1,KCB2,MARK)
C
C makes binary searches to retrieve type-based bond
C parameters.
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INTEGER BINPAR
      CHARACTER*4 K1, K2
      INTEGER NCB
      CHARACTER*4 KCB1(*), KCB2(*)
      INTEGER MARK
C local
      LOGICAL LTSTEQ
      EXTERNAL LTSTEQ
      CHARACTER*16 A, B
      INTEGER   I, J, N
C begin
      BINPAR=MARK
      IF (NCB.GT.0) THEN
C
C bonds ( permuation 1,2 <-> 2,1 )
      IF (LTSTEQ(K1,4,K2,4,.TRUE.)) THEN
      A=K1//K2
      ELSE
      A=K2//K1
      END IF
      I=1
      J=NCB
      N=(I+J)/2
      B=KCB1(N)//KCB2(N)
      DO WHILE (B.NE.A.AND..NOT.I.GE.J)
      IF (.NOT.LTSTEQ(B,8,A,8,.TRUE.)) THEN
      J=N-1
      ELSE
      I=N+1
      END IF
      N=MAX(1,(I+J)/2)
      B=KCB1(N)//KCB2(N)
      END DO
      IF (A.EQ.B) THEN
      BINPAR=N
      END IF
      END IF
      RETURN
      END
C
