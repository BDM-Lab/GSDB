      SUBROUTINE EANGLE(ET,ETV,IT,JT,KT,INTERE,NTHETA,
     &           ACTC,ACTB,ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,
     &           ANALYS,WGHT,VWGHT,IAC,SEGID,RESID,TYPE)
C
C Calculates angle energy.
C Front-end for actual energy routine.
C
C Author: Axel T. Brunger
C ========================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'update.inc'
      INCLUDE 'param.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ET,ETV
      DOUBLE PRECISION WGHT,VWGHT
      INTEGER IT(*), JT(*), KT(*), NTHETA, INTERE(*)
      DOUBLE PRECISION ACTC(*), ACTB(*), ACTUC(*), ACTUB(*)
      LOGICAL QACTC(*), QACTB(*), QACTUC(*), QACTUB(*)
      CHARACTER*4 ANALYS, IAC(*), SEGID(*), RESID(*), TYPE(*)
C local
      INTEGER I
      LOGICAL QUB
      DOUBLE PRECISION ETU, ETVU
C begin
C
C fill the unkown ACTC and ACTB values with type-based parameters
      IF (UPANGL) THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,ACTB,ACTUC,ACTUB,
     &            QACTC,QACTB,QACTUC,QACTUB,IAC,SEGID,RESID,TYPE)
      UPANGL=.FALSE.
      END IF
C
      CALL EANGLE2(ET,ETV,IT,JT,KT,INTERE,NTHETA,ACTC,ACTB,
     &             ANALYS,WGHT,VWGHT)
C
C check if we have to compute the Urey-Bradley term
      QUB=.FALSE.
      DO I=1,NTHETA
      QUB=QUB.OR.(ABS(ACTUC(I)).GT.RSMALL.AND.ABS(ACTUB(I)).GT.RSMALL)
      END DO
      IF (QUB) THEN
C compute the Urey-Bradley term
      CALL EBOND2(ETU,ETVU,IT,KT,INTERE,NTHETA,ACTUC,ACTUB,ANALYS,
     &                  WGHT,VWGHT)
      ET=ET+ETU
      ETV=ETV+ETVU
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE EANGLE2(ET,ETV,IT,JT,KT,INTERE,NTHETA,ACTC,ACTB,
     &                   ANALYS,WGHT,VWGHT)
C
C Calculates angle energy
C
C Author: Axel T. Brunger
C ========================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'pick.inc'
      DOUBLE PRECISION ET,ETV
      DOUBLE PRECISION WGHT,VWGHT
      INTEGER IT(*), JT(*), KT(*), NTHETA, INTERE(*)
      DOUBLE PRECISION ACTC(*), ACTB(*)
      CHARACTER*4 ANALYS
C local
      DOUBLE PRECISION PHI, SP, DI, DF, E, RECSP, CP
      DOUBLE PRECISION RIJX, RIJY, RIJZ, RIJ, RKJX, RKJY, RKJZ, RKJ
      INTEGER ITHETA, ILOC, IL, II
CCC      DOUBLE PRECISION WARNG
      DOUBLE PRECISION WTWO
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, MCONST, RADI
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, RADI=180.0D0/PI)
C===>parameter MCONST truncates the reciprocals and makes sure
C===>we don't devide by zero
      PARAMETER (MCONST=0.0001D0)
C local vector arrays
C===>VECDIM determines the optimal length of the vector loops
C===>on the CRAY this is 64
      INTEGER VECDIM
      PARAMETER (VECDIM=64)
      INTEGER LIND(VECDIM)
      DOUBLE PRECISION XI(VECDIM), XJ(VECDIM), XK(VECDIM)
      DOUBLE PRECISION YI(VECDIM), YJ(VECDIM), YK(VECDIM)
      DOUBLE PRECISION ZI(VECDIM), ZJ(VECDIM), ZK(VECDIM)
      DOUBLE PRECISION LCTC(VECDIM), LCTB(VECDIM)
      DOUBLE PRECISION DPRIJX(VECDIM), DPRIJY(VECDIM), DPRIJZ(VECDIM)
      DOUBLE PRECISION DPRKJX(VECDIM), DPRKJY(VECDIM), DPRKJZ(VECDIM)
      LOGICAL QANAL
C begin
C
      QANAL=.FALSE.
      IF (ANALYS.EQ.'ANAL') QANAL=.TRUE.
C
C initialization
      ET=0.0D0
CCC      WARNG=0.0D0
      WTWO=TWO*WGHT
C
C loop over all active angles
C make vector loops of length VECDIM
      ITHETA=0
      DO WHILE (ITHETA.LT.NTHETA)
      ILOC=0
C
C gather the data into temporary arrays suitable for vectorization
      DO WHILE (ILOC.LT.VECDIM.AND.ITHETA.LT.NTHETA)
      ITHETA=ITHETA+1
      IF (QANAL.OR.ABS(INTERE(IT(ITHETA))+INTERE(JT(ITHETA))
     &                +INTERE(KT(ITHETA))).LE.+2) THEN
      ILOC=ILOC+1
      LIND(ILOC)=ITHETA
      END IF
      END DO
      DO IL=1,ILOC
      XI(IL)=X(IT(LIND(IL)))
      YI(IL)=Y(IT(LIND(IL)))
      ZI(IL)=Z(IT(LIND(IL)))
      XJ(IL)=X(JT(LIND(IL)))
      YJ(IL)=Y(JT(LIND(IL)))
      ZJ(IL)=Z(JT(LIND(IL)))
      XK(IL)=X(KT(LIND(IL)))
      YK(IL)=Y(KT(LIND(IL)))
      ZK(IL)=Z(KT(LIND(IL)))
      LCTB(IL)=ACTB(LIND(IL))
      LCTC(IL)=ACTC(LIND(IL))
      END DO
C
C==================
C-BEGIN-VECTOR-LOOP
C==================
      DO IL=1,ILOC
C
C first compute the differences RIJ=RI-RJ and RKJ=RJ-RK
      RIJX=XI(IL)-XJ(IL)
      RIJY=YI(IL)-YJ(IL)
      RIJZ=ZI(IL)-ZJ(IL)
      RKJX=XK(IL)-XJ(IL)
      RKJY=YK(IL)-YJ(IL)
      RKJZ=ZK(IL)-ZJ(IL)
C
C compute the norm of RIJ, RKJ and set to MCONST if it is too small.
      RIJ=ONE/SQRT(MAX(MCONST,RIJX**2+RIJY**2+RIJZ**2))
      RKJ=ONE/SQRT(MAX(MCONST,RKJX**2+RKJY**2+RKJZ**2))
C
C normalize RIJ, RKJ
      RIJX=RIJX*RIJ
      RIJY=RIJY*RIJ
      RIJZ=RIJZ*RIJ
      RKJX=RKJX*RKJ
      RKJY=RKJY*RKJ
      RKJZ=RKJZ*RKJ
C
C compute CP = cos( phi )
      CP=RIJX*RKJX+RIJY*RKJY+RIJZ*RKJZ
C
C compute SP = sin( phi )   ( make sure  CP within boundaries )
      SP=SQRT(MAX(ZERO,ONE-CP**2))
C
C compute PHI ( make sure  CP within boundaries )
      PHI=ACOS(MIN(ONE,MAX(-ONE,CP)))
C
C======================================================================
C==> compute function (phi) and derivative
C==> replace this part if you want to use another functional form
C
C harmonic function: const* ( phi - phi0 )**2
      DI=PHI-LCTB(IL)
      DF=LCTC(IL)*DI
      E=DF*DI
      DF=WTWO*DF
C======================================================================
C
C accumulate energy
      ET=ET+E
C
C compute reciprocal of SP multiplied by the derivative of the
C function DF
      RECSP=DF/MAX(MCONST,SP)
C
C compute:
C    d (phi)
C DF -------,  etc.
C    d rij
      DPRIJX(IL)=RECSP*RIJ*(CP*RIJX-RKJX)
      DPRIJY(IL)=RECSP*RIJ*(CP*RIJY-RKJY)
      DPRIJZ(IL)=RECSP*RIJ*(CP*RIJZ-RKJZ)
      DPRKJX(IL)=RECSP*RKJ*(CP*RKJX-RIJX)
      DPRKJY(IL)=RECSP*RKJ*(CP*RKJY-RIJY)
      DPRKJZ(IL)=RECSP*RKJ*(CP*RKJZ-RIJZ)
C
      END DO
C ===============
C-END-VECTOR-LOOP
C ===============
C
C now scatter forces ( non-vectorizable )
      DO IL=1,ILOC
      II=LIND(IL)
      DX(IT(II))=DX(IT(II))+DPRIJX(IL)
      DY(IT(II))=DY(IT(II))+DPRIJY(IL)
      DZ(IT(II))=DZ(IT(II))+DPRIJZ(IL)
      DX(JT(II))=DX(JT(II))-DPRKJX(IL)-DPRIJX(IL)
      DY(JT(II))=DY(JT(II))-DPRKJY(IL)-DPRIJY(IL)
      DZ(JT(II))=DZ(JT(II))-DPRKJZ(IL)-DPRIJZ(IL)
      DX(KT(II))=DX(KT(II))           +DPRKJX(IL)
      DY(KT(II))=DY(KT(II))           +DPRKJY(IL)
      DZ(KT(II))=DZ(KT(II))           +DPRKJZ(IL)
      END DO
C
      END DO
C
      ETV=ET*VWGHT
      ET=ET*WGHT
C
      IF (ANALYS.EQ.'ANAL') THEN
      PCDATA(PCGEOM)=PHI*RADI
      PCDATA(PCENER)=ET
      PCDATA(PCEQUI)=ACTB(1)*RADI
      PCDATA(PCCONS)=ACTC(1)
      PCDATA(PCDEVI)=DI*RADI
      END IF
C
CCC      IF (WRNLEV.GT.0.AND.WARNG.GT.RSMALL) THEN
CCC   WRITE(6,'(A)') '%EANGLE-ERR: there are some linear bond angles.'
CCC      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE CODANG(NTHETA,IT,JT,KT,ACTC,ACTB,
     &                  ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,IAC,
     &                  SEGID,RESID,TYPE)
C
C Gets type-based parameters for bonds.
C
C Authors: Axel T. Brunger and Thomas Simonson
C ============================================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'timer.inc'
      INCLUDE 'param.inc'
      INTEGER NTHETA, IT(*), JT(*), KT(*)
      DOUBLE PRECISION ACTC(*), ACTB(*), ACTUC(*), ACTUB(*)
      LOGICAL QACTC(*), QACTB(*), QACTUC(*), QACTUB(*)
      CHARACTER*4 IAC(*), SEGID(*), RESID(*), TYPE(*)
C local
      INTEGER IFLAG, I, ICT
      INTEGER MARK
      PARAMETER (MARK=0)
C begin
C
      IFLAG=0
C
      DO I=1,NTHETA
C ignore atom-based parameter entries
      IF (QACTC(I).OR.QACTB(I).OR.QACTUC(I).OR.QACTUB(I)) THEN
C
      CALL BINANG(ICT,IAC(IT(I)),IAC(JT(I)),IAC(KT(I)),
     &            NCT,KCT1,KCT2,KCT3,MARK)
      IF (ICT.EQ.MARK) THEN
      IFLAG=IFLAG+1
      WRITE(6,'(A)')
     &' %CODANG-ERR: missing angle parameters %%%%%%%%%%%%%%%%%%%%%%%%'
      IF (QACTC(I)) WRITE(6,'(A)') '  angle energy constant missing.'
      IF (QACTB(I)) WRITE(6,'(A)') '  target angle value missing.'
      WRITE(6,'(9A)')
     & '  ATOM1: SEGId="',SEGID(IT(I)),'",  RESId="',RESID(IT(I)),
     & '",  NAME="',TYPE(IT(I)),'",  CHEMical="',IAC(IT(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM2: SEGId="',SEGID(JT(I)),'",  RESId="',RESID(JT(I)),
     & '",  NAME="',TYPE(JT(I)),'",  CHEMical="',IAC(JT(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM3: SEGId="',SEGID(KT(I)),'",  RESId="',RESID(KT(I)),
     & '",  NAME="',TYPE(KT(I)),'",  CHEMical="',IAC(KT(I)),'"'
      WRITE(6,'(A)')
     &' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ELSE
      IF (QACTC(I)) ACTC(I)=CTC(ICT)
      IF (QACTB(I)) ACTB(I)=CTB(ICT)
      IF (QACTUC(I)) ACTUC(I)=CTUC(ICT)
      IF (QACTUB(I)) ACTUB(I)=CTUB(ICT)
      END IF
      END IF
      END DO
      IF (IFLAG.GT.0) THEN
      CALL WRNDIE(-1,'CODANG','program will be aborted.')
      ELSE IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' CODANG: angle type-based parameters retrieved'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE BINANG(BINPAR,K1,K2,K3,NCT,KCT1,KCT2,KCT3,MARK)
C
C makes binary searches to retrieve type-based bond
C parameters.
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INTEGER BINPAR
      CHARACTER*4 K1, K2, K3
      INTEGER NCT
      CHARACTER*4 KCT1(*), KCT2(*), KCT3(*)
      INTEGER MARK
C local
      LOGICAL LTSTEQ
      EXTERNAL LTSTEQ
      CHARACTER*16 A, B
      INTEGER   I, J, N
C begin
      BINPAR=MARK
      IF (NCT.GT.0) THEN
C angles ( permutation 1,2,3 <-> 3,2,1 )
      IF (LTSTEQ(K1,4,K3,4,.TRUE.)) THEN
      A=K1//K2//K3
      ELSE
      A=K3//K2//K1
      END IF
      I=1
      J=NCT
      N=(I+J)/2
      B=KCT1(N)//KCT2(N)//KCT3(N)
      DO WHILE (B.NE.A.AND..NOT.I.GE.J)
      IF (.NOT.LTSTEQ(B,12,A,12,.TRUE.)) THEN
      J=N-1
      ELSE
      I=N+1
      END IF
      N=MAX(1,(I+J)/2)
      B=KCT1(N)//KCT2(N)//KCT3(N)
      END DO
      IF (A.EQ.B) THEN
      BINPAR=N
      END IF
      END IF
      RETURN
      END
C
