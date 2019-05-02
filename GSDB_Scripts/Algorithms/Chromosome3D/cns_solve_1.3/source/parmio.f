      SUBROUTINE PARRDR
C
C Parameter reader, including atom-based parameters and parameter
C learning.
C
C Authors: Axel Brunger and Thomas Simonson
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER IFLAGS, JFLAGS, KFLAGS, LFLAGS, NMAX, IC, NC
C begin
      NMAX=MAX(MAXCB,MAXCT,MAXCP,MAXCI)
C
      IFLAGS=ALLHP(INTEG4(NATOM))
      JFLAGS=ALLHP(INTEG4(NATOM))
      KFLAGS=ALLHP(INTEG4(NATOM))
      LFLAGS=ALLHP(INTEG4(NATOM))
      NC=ALLHP(INTEG4(NMAX))
      IC=ALLHP(INTEG4(NMAX))
      CALL PARRD2(HEAP(IFLAGS),HEAP(JFLAGS),HEAP(KFLAGS),HEAP(LFLAGS),
     & HEAP(NC),HEAP(IC))
      CALL FREHP(IC,INTEG4(NMAX))
      CALL FREHP(NC,INTEG4(NMAX))
      CALL FREHP(LFLAGS,INTEG4(NATOM))
      CALL FREHP(KFLAGS,INTEG4(NATOM))
      CALL FREHP(JFLAGS,INTEG4(NATOM))
      CALL FREHP(IFLAGS,INTEG4(NATOM))
      RETURN
      END
C
      SUBROUTINE PARRD2(IFLAGS,JFLAGS,KFLAGS,LFLAGS,NC,IC)
C
C Parameter reader
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'update.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'nbonds.inc'
      INTEGER IFLAGS(*), JFLAGS(*), KFLAGS(*), LFLAGS(*), NC(*), IC(*)
C local
      INTEGER I, I1, J1, J, SLOT
      INTEGER ISELCT, JSELCT, KSELCT, LSELCT
      DOUBLE PRECISION EPSL, SIG, EPS14, SIG14, EPSI, EPSI14
      DOUBLE PRECISION SIGI, SIGI14
      DOUBLE PRECISION A, B, A14, B14
C
      INTEGER NMULT, MULT, N
      PARAMETER (NMULT=10)
      DOUBLE PRECISION LCC(NMULT), LCB(NMULT)
      INTEGER LCD(NMULT)
      CHARACTER*4 AI, AJ, AK, AL, ACTION
      LOGICAL OLDSLOT, QCC(NMULT), QCB(NMULT)
      LOGICAL QCD(NMULT)
C parameter
      INTEGER MARK
      DOUBLE PRECISION ONE, TWO, FIVE, SIX, TWO56
      DOUBLE PRECISION RAD, R4BIG2, HALF, ZERO, T298
      PARAMETER (MARK=-9999, ONE=1.D0, TWO=2.D0, FIVE=5.D0, SIX=6.D0)
      PARAMETER (HALF=0.5D0, T298=298.D0)
      PARAMETER (RAD = PI/180.0D0)
      PARAMETER (ZERO=0.D0)
C
      TWO56=ONE/TWO**(FIVE/SIX)
      R4BIG2=R4BIG*TWO
C begin
C
C default values are initialized in main program
C
C set type-based parameter update flags each time when invoking
C this routine
      UPBOND=.TRUE.
      UPANGL=.TRUE.
      UPDIHE=.TRUE.
      UPIMPR=.TRUE.
      UPNBLK=.TRUE.
C
C read parameters from stream
      CALL PUSEND('PARRDR>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PARRDR>')
      CALL MISCOM('PARRDR>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-parameter')
C
C======================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL NEXTWD('PARRDR>')
      IF (WD(1:4).EQ.'TYPE') THEN
      ACTION='TYPE'
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
      ACTION='ATOM'
      ELSE IF (WD(1:4).EQ.'ALL ') THEN
      ACTION='ALL'
      ELSE
      CALL SAVEWD
      ACTION='ALL'
      END IF
C
      IF (ACTION.EQ.'TYPE'.OR.ACTION.EQ.'ALL') THEN
C
C reset type data base
      NCB=0
      NCT=0
      NNCCP=0
      NCI=0
      NCN=0
      END IF
C
      IF (ACTION.EQ.'ATOM'.OR.ACTION.EQ.'ALL') THEN
C
C reset atom data base
      DO I=1,NBOND
      QACBB(I)=.TRUE.
      QACBC(I)=.TRUE.
      END DO
      DO I=1,NTHETA
      QACTB(I)=.TRUE.
      QACTC(I)=.TRUE.
      QACTUB(I)=.TRUE.
      QACTUC(I)=.TRUE.
      END DO
      DO I=1,NPHI
      QACPB(I)=.TRUE.
      QACPC(I)=.TRUE.
      QACPD(I)=.TRUE.
      END DO
      DO I=1,NIMPHI
      QACIB(I)=.TRUE.
      QACIC(I)=.TRUE.
      QACID(I)=.TRUE.
      END DO
      DO I=1,NATOM
      LOOKUP(I)=0
      QLOOKU(I)=.TRUE.
      END DO
C reset special slot count for nonbonded atom-parameters
      NLJAT=0
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL PARSORT
      CALL PARWTR('OUTPUT','NORMAL')
C======================================================================
      ELSE IF (WD(1:4).EQ.'VERB') THEN
      CALL PARSORT
      CALL PARWTR('OUTPUT','VERBOSE')
C======================================================================
      ELSE IF (WD(1:4).EQ.'REDU') THEN
      CALL PARREDU(IFLAGS,IC,NC)
C======================================================================
      ELSE IF (WD(1:4).EQ.'LEAR') THEN
      CALL PARLEAR
C======================================================================
      ELSE IF (WD(1:4).EQ.'BOND') THEN
C
      CALL NEXTWD('BOND>')
      CALL SAVEWD
      IF (WD(1:1).EQ.'(') THEN
C atom-based parameter
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
C
C get energy constant. check for TOKEN.
      CALL NEXTA4('energy-constant=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('energy-constant=',LCC(1))
      DO I=1,NBOND
      IF ((IFLAGS(IB(I)).GT.0.AND.JFLAGS(JB(I)).GT.0).OR.
     &    (IFLAGS(JB(I)).GT.0.AND.JFLAGS(IB(I)).GT.0)) THEN
      ACBC(I)=LCC(1)
      QACBC(I)=.FALSE.
      END IF
      END DO
      END IF
C
C get minimum distance. check for TOKEN
      CALL NEXTA4('minimum-distance=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('minimum-distance=',LCB(1))
      DO I=1,NBOND
      IF ((IFLAGS(IB(I)).GT.0.AND.JFLAGS(JB(I)).GT.0).OR.
     &    (IFLAGS(JB(I)).GT.0.AND.JFLAGS(IB(I)).GT.0)) THEN
      ACBB(I)=LCB(1)
      QACBB(I)=.FALSE.
      END IF
      END DO
      END IF
C
      ELSE
C type-based parameter
C note the permutation I <-> J
      IF (NCB.GE.MAXCB) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCB (PARAM) parameter --> recompile program')
      ELSE
      NCB=NCB+1
      END IF
      CALL NEXTA4('BOND-I=',AI)
      CALL NEXTA4('BOND-J=',AJ)
      IF (LTSTEQ(AI,4,AJ,4,.TRUE.)) THEN
      KCB1(NCB)=AI
      KCB2(NCB)=AJ
      ELSE
      KCB1(NCB)=AJ
      KCB2(NCB)=AI
      END IF
C get energy constant and minimum distance.
      CALL NEXTF('energy-constant=',CBC(NCB))
      CALL NEXTF('minimum-distance=',CBB(NCB))
C
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'ANGL') THEN
C
      CALL NEXTWD('ANGLE>')
      CALL SAVEWD
      IF (WD(1:1).EQ.'(') THEN
C atom-based parameter
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(JFLAGS,KFLAGS,NATOM)
      CALL SELCTA(KFLAGS,KSELCT,X,Y,Z,.TRUE.)
C get energy constant. check for TOKEN.
      CALL NEXTA4('energy-constant=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('energy-constant=',LCC(1))
      DO I=1,NTHETA
      IF ((IFLAGS(IT(I)).GT.0.AND.JFLAGS(JT(I)).GT.0.
     &  AND.KFLAGS(KT(I)).GT.0).OR.(IFLAGS(KT(I)).GT.0.
     &  AND.JFLAGS(JT(I)).GT.0.AND.KFLAGS(IT(I)).GT.0)) THEN
      ACTC(I)=LCC(1)
      QACTC(I)=.FALSE.
      END IF
      END DO
      END IF
C get minimum angle. check for TOKEN.
      CALL NEXTA4('minimum-angle=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('minimum-angle=',LCB(1))
      DO I=1,NTHETA
      IF ((IFLAGS(IT(I)).GT.0.AND.JFLAGS(JT(I)).GT.0.
     &  AND.KFLAGS(KT(I)).GT.0).OR.(IFLAGS(KT(I)).GT.0.
     &  AND.JFLAGS(JT(I)).GT.0.AND.KFLAGS(IT(I)).GT.0)) THEN
      ACTB(I)=LCB(1)*RAD
      QACTB(I)=.FALSE.
      END IF
      END DO
      END IF
C check for Urey-Bradley terms
      CALL NEXTWD('PARRDR>')
      IF (WD(1:2).EQ.'UB') THEN
C get energy constant. check for TOKEN.
      CALL NEXTA4('ub-energy-constant=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('ub-energy-constant=',LCC(1))
      DO I=1,NTHETA
      IF ((IFLAGS(IT(I)).GT.0.AND.JFLAGS(JT(I)).GT.0.
     &  AND.KFLAGS(KT(I)).GT.0).OR.(IFLAGS(KT(I)).GT.0.
     &  AND.JFLAGS(JT(I)).GT.0.AND.KFLAGS(IT(I)).GT.0)) THEN
      ACTUC(I)=LCC(1)
      QACTUC(I)=.FALSE.
      END IF
      END DO
      END IF
C get minimum angle. check for TOKEN.
      CALL NEXTA4('ub-minimum-distance=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('ub-minimum-distance=',LCB(1))
      DO I=1,NTHETA
      IF ((IFLAGS(IT(I)).GT.0.AND.JFLAGS(JT(I)).GT.0.
     &  AND.KFLAGS(KT(I)).GT.0).OR.(IFLAGS(KT(I)).GT.0.
     &  AND.JFLAGS(JT(I)).GT.0.AND.KFLAGS(IT(I)).GT.0)) THEN
      ACTUB(I)=LCB(1)
      QACTUB(I)=.FALSE.
      END IF
      END DO
      END IF
      ELSE
C
C Urey-Bradley term is not specified -> set to zero!
      CALL SAVEWD
      DO I=1,NTHETA
      IF ((IFLAGS(IT(I)).GT.0.AND.JFLAGS(JT(I)).GT.0.
     &  AND.KFLAGS(KT(I)).GT.0).OR.(IFLAGS(KT(I)).GT.0.
     &  AND.JFLAGS(JT(I)).GT.0.AND.KFLAGS(IT(I)).GT.0)) THEN
      QACTUB(I)=.FALSE.
      QACTUC(I)=.FALSE.
      ACTUB(I)=ZERO
      ACTUC(I)=ZERO
      END IF
      END DO
      END IF
C
      ELSE
C type-based parameter
C note the permutation I,J,K <-> K,J,I
      IF (NCT.GE.MAXCT) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCT (PARAM) parameter --> recompile program')
      ELSE
      NCT=NCT+1
      END IF
      CALL NEXTA4('ANGLE-I=',AI)
      CALL NEXTA4('ANGLE-J=',KCT2(NCT))
      CALL NEXTA4('ANGLE-K=',AK)
      IF (LTSTEQ(AI,4,AK,4,.TRUE.)) THEN
      KCT1(NCT)=AI
      KCT3(NCT)=AK
      ELSE
      KCT1(NCT)=AK
      KCT3(NCT)=AI
      END IF
C get energy constant and minimum angle.
      CALL NEXTF('energy-constant=',CTC(NCT))
      CALL NEXTF('minimum-angle=',CTB(NCT))
      CTB(NCT)=CTB(NCT)*RAD
C check for Urey-Bradley terms
      CALL NEXTWD('PARRDR>')
      IF (WD(1:2).EQ.'UB') THEN
C get ub energy constant and minimum distance.
      CALL NEXTF('energy-constant=',CTUC(NCT))
      CALL NEXTF('minimum-angle=',CTUB(NCT))
      ELSE
C not specified - set UB term to zero
      CALL SAVEWD
      CTUC(NCT)=ZERO
      CTUB(NCT)=ZERO
      END IF
C
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'DIHE') THEN
C
      CALL NEXTWD('DIHE>')
      CALL SAVEWD
      IF (WD(1:1).EQ.'(') THEN
C atom-based parameter
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(JFLAGS,KFLAGS,NATOM)
      CALL SELCTA(KFLAGS,KSELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(KFLAGS,LFLAGS,NATOM)
      CALL SELCTA(LFLAGS,LSELCT,X,Y,Z,.TRUE.)
C get multiplicity
      CALL NEXTWD('DIHE>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('multiplicity=',MULT)
      IF (MULT.GT.NMULT) THEN
      CALL WRNDIE(-5,'DIHE>','multiplicity too large')
      END IF
      ELSE
      MULT=1
      CALL SAVEWD
      ENDIF
C
C read small array of parameters, of size MULT
      DO I=1,MULT
C get energy constant. check for TOKEN
      CALL NEXTA4('energy-constant=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('energy-constant=',LCC(I))
      QCC(I)=.TRUE.
      ELSE
      QCC(I)=.FALSE.
      END IF
C get periodicity. If not an integer, do not use.
      CALL NEXTA4('periodicity=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTI('periodicity=',LCD(I))
      QCD(I)=.TRUE.
      ELSE
      QCD(I)=.FALSE.
      END IF
C get minimum angle. If not a real number, do not use.
      CALL NEXTA4('minimum-angle=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('minimum-angle=',LCB(I))
      QCB(I)=.TRUE.
      ELSE
      QCB(I)=.FALSE.
      END IF
      END DO
C
      I=1
      DO WHILE (I.LE.NPHI)
C
      IF ((IFLAGS(IP(I)).GT.0.AND.JFLAGS(JP(I)).GT.0.
     1  AND.KFLAGS(KP(I)).GT.0.AND.LFLAGS(LP(I)).GT.0).OR.
     2 (IFLAGS(LP(I)).GT.0.AND.JFLAGS(KP(I)).GT.0.
     3  AND.KFLAGS(JP(I)).GT.0.AND.LFLAGS(IP(I)).GT.0)) THEN
C
      DO J=1,MULT
      IF (QCC(J)) THEN
      ACPC(I)=LCC(J)
      QACPC(I)=.FALSE.
      END IF
      IF (QCB(J)) THEN
      ACPB(I)=LCB(J)*RAD
      QACPB(I)=.FALSE.
      END IF
      IF (QCD(J)) THEN
      ACPD(I)=LCD(J)
      QACPD(I)=.FALSE.
      END IF
      I=I+1
      END DO
C
      ELSE
      I=I+1
      ENDIF
C
      END DO
C
      ELSE
C type-based parameter
C note the permutation I,J,K,L <-> L,K,J,I.
      IF (NNCCP.GE.MAXCP) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCP (PARAM) parameter --> recompile program')
      ELSE
      NNCCP=NNCCP+1
      END IF
      CALL NEXTA4('DIHEDRAL-I=',AI)
      CALL NEXTA4('DIHEDRAL-J=',AJ)
      CALL NEXTA4('DIHEDRAL-K=',AK)
      CALL NEXTA4('DIHEDRAL-L=',AL)
      IF (LTSTEQ(AJ//AI,8,AK//AL,8,.TRUE.)) THEN
      KCP1(NNCCP)=AI
      KCP2(NNCCP)=AJ
      KCP3(NNCCP)=AK
      KCP4(NNCCP)=AL
      ELSE
      KCP1(NNCCP)=AL
      KCP2(NNCCP)=AK
      KCP3(NNCCP)=AJ
      KCP4(NNCCP)=AI
      END IF
C get multiplicity
      CALL NEXTWD('DIHE>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('multiplicity=',MULT)
      IF (MULT.GT.NMULT) THEN
      CALL WRNDIE(-5,'DIHE>','multiplicity too large')
      END IF
      ELSE
      MULT=1
      CALL SAVEWD
      ENDIF
C get energy constant, periodicity, minimum angle.
      CALL NEXTF('energy-constant=',CPC(NNCCP))
      CALL NEXTI('periodicity=',CPD(NNCCP))
      CALL NEXTF('minimum-angle=',CPB(NNCCP))
      CPB(NNCCP)=CPB(NNCCP)*RAD
C get more energy constants, ...
      DO N=2,MULT
      IF (NNCCP.GE.MAXCP) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCP (PARAM) parameter --> recompile program')
      ELSE
      NNCCP=NNCCP+1
      KCP1(NNCCP)=KCP1(NNCCP-1)
      KCP2(NNCCP)=KCP2(NNCCP-1)
      KCP3(NNCCP)=KCP3(NNCCP-1)
      KCP4(NNCCP)=KCP4(NNCCP-1)
      CALL NEXTF('energy-constant=',CPC(NNCCP))
      CALL NEXTI('periodicity=',CPD(NNCCP))
      CALL NEXTF('minimum-angle=',CPB(NNCCP))
      CPB(NNCCP)=CPB(NNCCP)*RAD
      END IF
      END DO
C
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'IMPR') THEN
C
      CALL NEXTWD('IMPR>')
      CALL SAVEWD
      IF (WD(1:1).EQ.'(') THEN
C atom-based parameter
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(JFLAGS,KFLAGS,NATOM)
      CALL SELCTA(KFLAGS,KSELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(KFLAGS,LFLAGS,NATOM)
      CALL SELCTA(LFLAGS,LSELCT,X,Y,Z,.TRUE.)
C get multiplicity
      CALL NEXTWD('IMPR>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('multiplicity=',MULT)
      ELSE
      MULT=1
      CALL SAVEWD
      ENDIF
C
C read small array of parameters, of size MULT
      DO I=1,MULT
C get energy constant. check for TOKEN
      CALL NEXTA4('energy-constant=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('energy-constant=',LCC(I))
      QCC(I)=.TRUE.
      ELSE
      QCC(I)=.FALSE.
      END IF
C get periodicity. If not an integer, do not use.
      CALL NEXTA4('periodicity=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTI('periodicity=',LCD(I))
      QCD(I)=.TRUE.
      ELSE
      QCD(I)=.FALSE.
      END IF
C get minimum angle. If not a real number, do not use.
      CALL NEXTA4('minimum-angle=',ACTION)
      IF (ACTION.NE.'TOKE') THEN
      CALL SAVEWD
      CALL NEXTF('minimum-angle=',LCB(I))
      QCB(I)=.TRUE.
      ELSE
      QCB(I)=.FALSE.
      END IF
      END DO
C
      I=1
      DO WHILE (I.LE.NIMPHI)
C
      IF ((IFLAGS(IM(I)).GT.0.AND.JFLAGS(JM(I)).GT.0.
     1  AND.KFLAGS(KM(I)).GT.0.AND.LFLAGS(LM(I)).GT.0).OR.
     2 (IFLAGS(LM(I)).GT.0.AND.JFLAGS(KM(I)).GT.0.
     3  AND.KFLAGS(JM(I)).GT.0.AND.LFLAGS(IM(I)).GT.0)) THEN
C
      DO J=1,MULT
      IF (QCC(J)) THEN
      ACIC(I)=LCC(J)
      QACIC(I)=.FALSE.
      END IF
      IF (QCB(J)) THEN
      ACIB(I)=LCB(J)*RAD
      QACIB(I)=.FALSE.
      END IF
      IF (QCD(J)) THEN
      ACID(I)=LCD(J)
      QACID(I)=.FALSE.
      END IF
      I=I+1
      END DO
C
      ELSE
      I=I+1
      END IF
C
      END DO
C
      ELSE
C type-based parameter
C note the permutation  I,J,K,L <-> L,K,J,I
      IF (NCI.GE.MAXCI) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCI (PARAM) parameter --> recompile program')
      ELSE
      NCI=NCI+1
      END IF
      CALL NEXTA4('IMPROPER-I=',AI)
      CALL NEXTA4('IMPROPER-J=',AJ)
      CALL NEXTA4('IMPROPER-K=',AK)
      CALL NEXTA4('IMPROPER-L=',AL)
      IF (LTSTEQ(AI//AJ,8,AL//AK,8,.TRUE.)) THEN
      KCI1(NCI)=AI
      KCI2(NCI)=AJ
      KCI3(NCI)=AK
      KCI4(NCI)=AL
      ELSE
      KCI1(NCI)=AL
      KCI2(NCI)=AK
      KCI3(NCI)=AJ
      KCI4(NCI)=AI
      END IF
      CALL NEXTWD('IMPR>')
      IF (WD(1:4).EQ.'MULT') THEN
      CALL NEXTI('multiplicity=',MULT)
      ELSE
      MULT=1
      CALL SAVEWD
      ENDIF
C get energy constant, periodicity, minimum angle.
      CALL NEXTF('energy-constant=',CIC(NCI))
      CALL NEXTI('periodicity=',CID(NCI))
      CALL NEXTF('minimum-angle=',CIB(NCI))
      CIB(NCI)=CIB(NCI)*RAD
C
C get more energy constants, ...
      DO N=2,MULT
      IF (NCI.GE.MAXCI) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCI (PARAM) parameter --> recompile program')
      ELSE
      NCI=NCI+1
      KCI1(NCI)=KCI1(NCI-1)
      KCI2(NCI)=KCI2(NCI-1)
      KCI3(NCI)=KCI3(NCI-1)
      KCI4(NCI)=KCI4(NCI-1)
      CALL NEXTF('energy-constant=',CIC(NCI))
      CALL NEXTI('periodicity=',CID(NCI))
      CALL NEXTF('minimum-angle=',CIB(NCI))
      CIB(NCI)=CIB(NCI)*RAD
      END IF
      END DO
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'NBON') THEN
      CALL NBDSET('PARSE')
C
      ELSE IF (WD(1:4).EQ.'NONB') THEN
C
C This section reads Lennard-Jones parameters.
C
C Lennard-Jones potential:
C
C                   sigma  12       sigma   6
C     E = 4 eps  ( ( ----- )    -  ( ----- )     )
C                      R               R
C
C                   6 ---
C     R     = sigma / 2              E    = - eps
C      min                            min
C               12                               6
C     A= 4 sigma   eps                B = 4 sigma  eps
C
C
C Combination rules:
C
C                  sigma(11) + sigma(22)               ----------------
C     sigma(12) = ---------------------   eps(12) =  \/ eps(11)*eps(22)
C
C                           2
C
C
      CALL NEXTWD('NONB>')
      CALL SAVEWD
      IF (WD(1:1).NE.'(') THEN
C
C type-based parameter
C
      IF (NCN+NLJAT.GE.MAXCN) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCN (PARAM) parameter --> recompile program')
      END IF
C
      NCN=NCN+1
      CALL NEXTA4('Lennard-Jones-chemical-type=',CNAC(NCN))
      CALL NEXTF('epsilon=',EPSL)
      CALL NEXTF('sigma=',SIG)
      CALL NEXTF('1:4-epsilon=',EPS14)
      CALL NEXTF('1:4-sigma=',SIG14)
C repel parameters
      CNBVR(NCN,NCN)=TWO*SIG*TWO56
      CBVR14(NCN,NCN)=TWO*SIG14*TWO56
C diagonal Lennard-Jones parameters
      CALL LJABES(CNBA(NCN,NCN),CNBB(NCN,NCN),EPSL,SIG)
      CALL LJABES(CNBA14(NCN,NCN),CNBB14(NCN,NCN),EPS14,SIG14)
C apply combination rules to obtain missing cross terms
      DO I=1,NCN-1
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJABES(CNBA(I,NCN),CNBB(I,NCN),SQRT(EPSI*EPSL),
     &     HALF*(SIGI+SIG))
      CNBB(NCN,I)=CNBB(I,NCN)
      CNBA(NCN,I)=CNBA(I,NCN)
      CNBVR(I,NCN)=(SIGI+SIG)*TWO56
      CNBVR(NCN,I)=CNBVR(I,NCN)
C 1-4 cross terms
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJABES(CNBA14(I,NCN),CNBB14(I,NCN),
     &           SQRT(EPSI14*EPS14),HALF*(SIGI14+SIG14))
      CNBB14(NCN,I)=CNBB14(I,NCN)
      CNBA14(NCN,I)=CNBA14(I,NCN)
      CBVR14(I,NCN)=(SIGI14+SIG14)*TWO56
      CBVR14(NCN,I)=CBVR14(I,NCN)
      END DO
C cross terms between type-based and atom-based parameters
      DO I=MAXCN-NLJAT+1,MAXCN
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJABES(CNBA(I,NCN),CNBB(I,NCN),SQRT(EPSI*EPSL),
     &     HALF*(SIGI+SIG))
      CNBB(NCN,I)=CNBB(I,NCN)
      CNBA(NCN,I)=CNBA(I,NCN)
C 1-4 cross terms
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJABES(CNBA14(I,NCN),CNBB14(I,NCN),
     &           SQRT(EPSI14*EPS14),HALF*(SIGI14+SIG14))
      CNBB14(NCN,I)=CNBB14(I,NCN)
      CNBA14(NCN,I)=CNBA14(I,NCN)
      END DO
C
      ELSE
C
C atom-based parameter
C
C atom-based parameters are stored at the end of the parameter tables
C (CNBA, CNBB) starting from the top (MAXCN) down.
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL NEXTF('epsilon=',EPSL)
      CALL NEXTF('sigma=',SIG)
      CALL NEXTF('1:4-epsilon=',EPS14)
      CALL NEXTF('1:4-sigma=',SIG14)
      IF (ISELCT.GT.0) THEN
C check to see if we actually need a new slot
      DO J=1,NATOM
      IF (IFLAGS(J).GT.0) SLOT=LOOKUP(J)
      END DO
      OLDSLOT=SLOT.NE.0.AND.SLOT.GT.NCN
      DO J=1,NATOM
      IF (IFLAGS(J).GT.0) THEN
      OLDSLOT=OLDSLOT.AND.(SLOT.EQ.LOOKUP(J))
      ELSE
      OLDSLOT=OLDSLOT.AND.(SLOT.NE.LOOKUP(J))
      END IF
      END DO
      IF (OLDSLOT) THEN
C we can use the old slot -- the new atom selection is
C a one-to-one replacement of a previous one
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(A)')
     & ' PARRD2-INFO: replacing previous nonbonded slot.'
      END IF
      ELSE
C
      IF (WRNLEV.GT.5) THEN
      WRITE(6,'(A)')
     & ' PARRD2-INFO: setting up new nonbonded slot.'
      END IF
C
C set up a new slot in parameter tables
      IF (NCN+NLJAT.GE.MAXCN) THEN
      CALL WRNDIE(-5,'PARRDR',
     & 'exceeded MAXCN (PARAM) parameter --> recompile program')
      ELSE
      NLJAT=NLJAT+1
      END IF
      SLOT=MAXCN-NLJAT+1
      END IF
C
C now have LOOKUP point into that slot
      DO J=1,NATOM
      IF (IFLAGS(J).GT.0) THEN
      LOOKUP(J)=SLOT
      QLOOKU(J)=.FALSE.
      END IF
      END DO
C repel parameters
      CNBVR(SLOT,SLOT)=TWO*SIG*TWO56
      CBVR14(SLOT,SLOT)=TWO*SIG14*TWO56
C diagonal Lennard-Jones terms
      CALL LJABES(CNBA(SLOT,SLOT),CNBB(SLOT,SLOT),EPSL,SIG)
      CALL LJABES(CNBA14(SLOT,SLOT),CNBB14(SLOT,SLOT),EPS14,SIG14)
C apply combination rules to missing cross terms
      DO I=1,NCN
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJABES(CNBA(I,SLOT),CNBB(I,SLOT),SQRT(EPSI*EPSL),
     &     HALF*(SIGI+SIG))
      CNBB(SLOT,I)=CNBB(I,SLOT)
      CNBA(SLOT,I)=CNBA(I,SLOT)
      CNBVR(I,SLOT)=(SIGI+SIG)*TWO56
      CNBVR(SLOT,I)=CNBVR(I,SLOT)
C 1-4 cross terms
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJABES(CNBA14(I,SLOT),CNBB14(I,SLOT),
     &           SQRT(EPSI14*EPS14),HALF*(SIGI14+SIG14))
      CNBB14(SLOT,I)=CNBB14(I,SLOT)
      CNBA14(SLOT,I)=CNBA14(I,SLOT)
      CBVR14(I,SLOT)=(SIGI14+SIG14)*TWO56
      CBVR14(SLOT,I)=CBVR14(I,SLOT)
      END DO
C cross terms with atom-based parameters
      DO I=MAXCN-NLJAT+1,MAXCN
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJABES(CNBA(I,SLOT),CNBB(I,SLOT),SQRT(EPSI*EPSL),
     &     HALF*(SIGI+SIG))
      CNBB(SLOT,I)=CNBB(I,SLOT)
      CNBA(SLOT,I)=CNBA(I,SLOT)
      CNBVR(I,SLOT)=(SIGI+SIG)*TWO56
      CNBVR(SLOT,I)=CNBVR(I,SLOT)
C 1-4 terms
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJABES(CNBA14(I,SLOT),CNBB14(I,SLOT),
     &           SQRT(EPSI14*EPS14),HALF*(SIGI14+SIG14))
      CNBB14(SLOT,I)=CNBB14(I,SLOT)
      CNBA14(SLOT,I)=CNBA14(I,SLOT)
      CBVR14(I,SLOT)=(SIGI14+SIG14)*TWO56
      CBVR14(SLOT,I)=CBVR14(I,SLOT)
      END DO
C
      END IF
C
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'NBFI') THEN
C
C This section reads nonbond fixes
      CALL NEXTWD('NBFIx=')
      CALL SAVEWD
      IF (WD(1:1).EQ.'(') THEN
C atom-based section
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
      CALL COPYI4(IFLAGS,JFLAGS,NATOM)
      CALL SELCTA(JFLAGS,JSELCT,X,Y,Z,.TRUE.)
      I1=MARK
      J1=MARK
C
      IF (ISELCT.GT.0) THEN
C get the old slot id for first atom
      DO J=1,NATOM
      IF (IFLAGS(J).GT.0) SLOT=LOOKUP(J)
      END DO
      OLDSLOT=SLOT.NE.0.AND.SLOT.GT.NCN
      DO J=1,NATOM
      IF (IFLAGS(J).GT.0) THEN
      OLDSLOT=OLDSLOT.AND.(SLOT.EQ.LOOKUP(J))
      ELSE
      OLDSLOT=OLDSLOT.AND.(SLOT.NE.LOOKUP(J))
      END IF
      END DO
      IF (OLDSLOT) I1=SLOT
      END IF
C
      IF (JSELCT.GT.0) THEN
C get the old slot id for second atom
      DO J=1,NATOM
      IF (JFLAGS(J).GT.0) SLOT=LOOKUP(J)
      END DO
      OLDSLOT=SLOT.NE.0.AND.SLOT.GT.NCN
      DO J=1,NATOM
      IF (JFLAGS(J).GT.0) THEN
      OLDSLOT=OLDSLOT.AND.(SLOT.EQ.LOOKUP(J))
      ELSE
      OLDSLOT=OLDSLOT.AND.(SLOT.NE.LOOKUP(J))
      END IF
      END DO
      IF (OLDSLOT) J1=SLOT
      END IF
C
      IF (I1.EQ.MARK) THEN
      WRITE(6,'(2A)')
     & ' %NBFIX-ERR: atom-based nonbonded entry for first selection ',
     & 'not found.'
      END IF
      IF (J1.EQ.MARK) THEN
      WRITE(6,'(2A)')
     & ' %NBFIX-ERR: atom-based nonbonded entry for second selection ',
     & 'not found.'
      END IF
C
      ELSE
C type-based section
      CALL NEXTA4('NBFIX-first-atom=',AI)
      CALL NEXTA4('NBFIX-second-atom=',AJ)
      I1=SRCHC4(AI,CNAC,1,NCN,MARK)
      J1=SRCHC4(AJ,CNAC,1,NCN,MARK)
      IF (I1.EQ.MARK.OR.J1.EQ.MARK) THEN
      WRITE(6,'(5A)')
     & ' %NBFIX-ERR: nonbonded entries for atom types ',
     &  AI,' ',AJ,' not found'
      END IF
      END IF
      CALL NEXTF('A=',A)
      CALL NEXTF('B=',B)
      CALL NEXTF('A14=',A14)
      CALL NEXTF('B14=',B14)
C
      IF (I1.NE.MARK.AND.J1.NE.MARK) THEN
      CNBA(I1,J1)=A
      CNBB(I1,J1)=B
      CNBA(J1,I1)=CNBA(I1,J1)
      CNBB(J1,I1)=CNBB(I1,J1)
C
      CALL LJESAB(EPSI,SIGI,A,B)
      CNBVR(I1,J1)=TWO*SIGI*TWO56
      CNBVR(J1,I1)=CNBVR(I1,J1)
C
      CNBA14(I1,J1)=A14
      CNBB14(I1,J1)=B14
      CNBA14(J1,I1)=CNBA14(I1,J1)
      CNBB14(J1,I1)=CNBB14(I1,J1)
C
      CALL LJESAB(EPSI14,SIGI14,A14,B14)
      CBVR14(I1,J1)=TWO*SIGI14*TWO56
      CBVR14(J1,I1)=CBVR14(I1,J1)
C
      END IF
C
      ELSE
      CALL CHKEND('PARRDR>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C
      CALL PARSORT
      RETURN
      END
C=====================================================================
      SUBROUTINE PARSORT
C
C parameter sorting is performed for bonds, angles, dihedrals,
C and impropers to enable a binary search.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'param.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cstack.inc'
      INCLUDE 'timer.inc'
C local
      INTEGER I, J
      EXTERNAL ORDPAR
      LOGICAL ORDPAR, COND
C pointer
      INTEGER IDX, CWORK, WORK
C begin
      I=MAX(NCB,NCT,NNCCP,NCI,NCN)
      IDX=ALLHP(INTEG4(I))
      CWORK=CALLST(ICHAR4(I))
      WORK=ALLHP(IREAL8(I))
C
C sort bond parameters.
      IF (NCB.GT.0) THEN
      CALL SORTP(NCB,HEAP(IDX),ORDPAR,1,0,0,0,0,0,0,0)
      CALL AINDC4(HEAP(IDX),KCB1,NCB,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCB2,NCB,CSTACK(CWORK))
      CALL AINDR8(HEAP(IDX),CBB,NCB,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CBC,NCB,HEAP(WORK))
      END IF
C
C sort angle parameters
      IF (NCT.GT.0) THEN
      CALL SORTP(NCT,HEAP(IDX),ORDPAR,2,0,0,0,0,0,0,0)
      CALL AINDC4(HEAP(IDX),KCT1,NCT,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCT2,NCT,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCT3,NCT,CSTACK(CWORK))
      CALL AINDR8(HEAP(IDX),CTB,NCT,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CTC,NCT,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CTUB,NCT,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CTUC,NCT,HEAP(WORK))
      END IF
C
C sort dihedral parameters
      IF (NNCCP.GT.0) THEN
      CALL SORTP(NNCCP,HEAP(IDX),ORDPAR,3,0,0,0,0,0,0,0)
      CALL AINDC4(HEAP(IDX),KCP1,NNCCP,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCP2,NNCCP,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCP3,NNCCP,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCP4,NNCCP,CSTACK(CWORK))
      CALL AINDR8(HEAP(IDX),CPB,NNCCP,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CPC,NNCCP,HEAP(WORK))
      CALL AINDX4(HEAP(IDX),CPD,NNCCP,HEAP(WORK))
      END IF
C
C sort improper parameters
      IF (NCI.GT.0) THEN
      CALL SORTP(NCI,HEAP(IDX),ORDPAR,4,0,0,0,0,0,0,0)
      CALL AINDC4(HEAP(IDX),KCI1,NCI,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCI2,NCI,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCI3,NCI,CSTACK(CWORK))
      CALL AINDC4(HEAP(IDX),KCI4,NCI,CSTACK(CWORK))
      CALL AINDR8(HEAP(IDX),CIB,NCI,HEAP(WORK))
      CALL AINDR8(HEAP(IDX),CIC,NCI,HEAP(WORK))
      CALL AINDX4(HEAP(IDX),CID,NCI,HEAP(WORK))
      END IF
C
      CALL FREHP(WORK,IREAL8(I))
      CALL FREHP(IDX,INTEG4(I))
      CALL CFREST(ICHAR4(I))
C
C now check for repetitions
      DO I=1,NCB-1
      IF (KCB1(I).EQ.KCB1(I+1).AND.KCB2(I).EQ.KCB2(I+1)) THEN
      WRITE(6,'(4A)')
     1' %PARRDR-info: duplication of bond ',KCB1(I),' ',KCB2(I)
      END IF
      END DO
C
      DO I=1,NCT-1
      COND=(KCT1(I).EQ.KCT1(I+1).AND.KCT2(I).EQ.KCT2(I+1).AND.
     1      KCT3(I).EQ.KCT3(I+1))
      IF (COND) THEN
      WRITE(6,'(6A)')
     &' %PARRDR-info: duplication of angle ',KCT1(I),' ',KCT2(I),
     &' ',KCT3(I)
      END IF
      END DO
C
      DO I=1,NNCCP-1
      COND=(KCP1(I).EQ.KCP1(I+1).AND.KCP2(I).EQ.KCP2(I+1).AND.
     1      KCP3(I).EQ.KCP3(I+1).AND.KCP4(I).EQ.KCP4(I+1))
      IF (COND.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(8A)')
     &' PARRDR-info: multiple dihedral entry for ',KCP1(I),
     & ' ',KCP2(I),' ',KCP3(I),' ',KCP4(I)
      END IF
      END DO
C
      DO I=1,NCI-1
      COND=(KCI1(I).EQ.KCI1(I+1).AND.KCI2(I).EQ.KCI2(I+1).AND.
     1      KCI3(I).EQ.KCI3(I+1).AND.KCI4(I).EQ.KCI4(I+1))
      IF (COND.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(8A)')
     &' PARRDR-info: multiple improper entry for ',KCI1(I),
     & ' ',KCI2(I),' ',KCI3(I),' ',KCI4(I)
      END IF
      END DO
C
      DO I=1,NCN
      DO J=1,I-1
      IF (CNAC(I).EQ.CNAC(J)) THEN
      WRITE(6,'(2A)')
     &' %PARRDR-info: duplication of nonbonded entry ',CNAC(I)
      END IF
      END DO
      END DO
C
      RETURN
      END
C=====================================================================
      LOGICAL FUNCTION ORDPAR(K,L,MODE,X2,X3,X4,X5,X6,X7,X8)
C
C defines canonical order bond parameters
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'param.inc'
      INTEGER K, L, MODE, X2, X3, X4, X5, X6, X7, X8
C local
      EXTERNAL LTSTEQ
      LOGICAL  LTSTEQ
      CHARACTER*(16) AK, AL
C begin
      IF (MODE.EQ.1) THEN
C bonds
      AK=KCB1(K)//KCB2(K)
      AL=KCB1(L)//KCB2(L)
      ORDPAR=LTSTEQ(AK,8,AL,8,.TRUE.)
      IF (AK.EQ.AL) ORDPAR=(K.LE.L)
      ELSE IF (MODE.EQ.2) THEN
C angles
      AK=KCT1(K)//KCT2(K)//KCT3(K)
      AL=KCT1(L)//KCT2(L)//KCT3(L)
      ORDPAR=LTSTEQ(AK,12,AL,12,.TRUE.)
      IF (AK.EQ.AL) ORDPAR=(K.LE.L)
      ELSE IF (MODE.EQ.3) THEN
C dihedrals
      AK=KCP1(K)//KCP2(K)//KCP3(K)//KCP4(K)
      AL=KCP1(L)//KCP2(L)//KCP3(L)//KCP4(L)
      ORDPAR=LTSTEQ(AK,16,AL,16,.TRUE.)
      IF (AK.EQ.AL) ORDPAR=(K.LE.L)
      ELSE IF (MODE.EQ.4) THEN
C impropers
      AK=KCI1(K)//KCI2(K)//KCI3(K)//KCI4(K)
      AL=KCI1(L)//KCI2(L)//KCI3(L)//KCI4(L)
      ORDPAR=LTSTEQ(AK,16,AL,16,.TRUE.)
      IF (AK.EQ.AL) ORDPAR=(K.LE.L)
      ELSE IF (MODE.EQ.6) THEN
C Lennard-Jones index
      ORDPAR=LTSTEQ(CNAC(K),4,CNAC(L),4,.TRUE.)
      IF (CNAC(K).EQ.CNAC(L)) ORDPAR=(K.LE.L)
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE LJESAB(EPS,SIG,A,B)
C
C calculate Lennard-Jones epsilon, sigma from A, B.
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      DOUBLE PRECISION EPS, SIG, A, B
C local
      DOUBLE PRECISION ONE, SIX, QUART, ZERO, TWELVE, TOL
      PARAMETER (ONE=1.0D0, SIX=6.0D0, QUART=0.25D0, ZERO=0.0D0)
      PARAMETER (TWELVE=12.0D0, TOL=0.01D0)
C begin
      IF (ABS(A**(ONE/TWELVE)).GT.TOL
     &    .AND.ABS(B**(ONE/SIX)).GT.TOL) THEN
      SIG=(A/B)**(ONE/SIX)
      EPS=QUART*B*B/A
      ELSE
      SIG=ZERO
      EPS=ZERO
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE LJABES(A,B,EPS,SIG)
C
C calculates A, B from Lennard-Jones epsilon and sigma
      IMPLICIT NONE
C input/output
      DOUBLE PRECISION A, B, EPS, SIG
C local
      DOUBLE PRECISION FOUR
      PARAMETER (FOUR=4.0D0)
C begin
      A=FOUR*EPS*SIG**12
      B=FOUR*EPS*SIG**6
      RETURN
      END
C=====================================================================
      SUBROUTINE PARWTR(FILE,MODE)
C
C parameter file write routine
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'update.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) FILE, MODE
C local
      DOUBLE PRECISION EPS, SIG, EPS14, SIG14, SD
      DOUBLE PRECISION SIGI, EPSI, SIGI14, EPSI14
      DOUBLE PRECISION SIGJ, EPSJ, SIGJ14, EPSJ14
      INTEGER I, J, UNIT
      LOGICAL ERROR, COND, COND14
C
      INTEGER MULT, II
      DOUBLE PRECISION KTBY2
C parameter
      DOUBLE PRECISION RAD, T298, BIG99, TWO
      PARAMETER (RAD=PI/180.0D0, T298=298.D0, BIG99=999999.D0)
      PARAMETER (TWO=2.0D0)
C begin
C
      KTBY2=KBOLTZ*T298/TWO
C
      CALL ASSFIL(FILE,UNIT,'WRITE','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
C
      IF (FILE.NE.'OUTPUT') CALL WRTITL(UNIT,0)
C
      IF (MODE.EQ.'NORMAL') THEN
      IF (NCB.GT.0) THEN
      WRITE(UNIT,'(A)')
      DO I=1,NCB
      IF (CBC(I).GT.RSMALL) THEN
      SD=SQRT(KTBY2/CBC(I))
      ELSE
      SD=BIG99
      END IF
      WRITE(UNIT,'(5A,F10.3,A,F10.3,A,F10.3)') ' BOND  ',KCB1(I),' ',
     &      KCB2(I),'  ',CBC(I),' {sd=',SD,'} ',CBB(I)
      END DO
      END IF
C
      IF (NCT.GT.0) THEN
C
      WRITE(UNIT,'(A)')
      DO I=1,NCT
      IF (CTC(I).GT.RSMALL) THEN
      SD=SQRT(KTBY2/CTC(I))/RAD
      ELSE
      SD=BIG99
      END IF
      IF (ABS(CTUC(I)).LT.RSMALL.OR.ABS(CTUB(I)).LT.RSMALL) THEN
C no Urey-Bradley term present
      WRITE(UNIT,'(7A,F9.2,A,F10.3,A,F11.4)') ' ANGLe  ',KCT1(I),' ',
     &      KCT2(I),' ',KCT3(I),'  ',CTC(I),' {sd=',SD,'} ',CTB(I)/RAD
      ELSE
C Urey-Bradley term is present
      WRITE(UNIT,'(7A,F9.2,A,F10.3,A,F11.4,A,F10.3,A,F10.3)')
     & ' ANGLe  ',KCT1(I),' ',
     &      KCT2(I),' ',KCT3(I),'  ',CTC(I),' {sd=',SD,'} ',CTB(I)/RAD,
     &  '    UB ',CTUC(I),'  ', CTUB(I)
      END IF
      END DO
      END IF
C
      IF (NNCCP.GT.0) THEN
      WRITE(UNIT,'(A)')
C
      I=1
      DO WHILE (I.LE.NNCCP)
C
C determine multiplicity
      MULT=1
      COND=(I.LT.NNCCP)
      IF (COND) THEN
      COND=(KCP1(I+1).EQ.KCP1(I).AND.KCP2(I+1).EQ.KCP2(I).AND.
     &      KCP3(I+1).EQ.KCP3(I).AND.KCP4(I+1).EQ.KCP4(I))
      END IF
      DO WHILE (COND)
      MULT=MULT+1
      I=I+1
      IF (COND) THEN
      COND=(KCP1(I+1).EQ.KCP1(I).AND.KCP2(I+1).EQ.KCP2(I).AND.
     &      KCP3(I+1).EQ.KCP3(I).AND.KCP4(I+1).EQ.KCP4(I))
      END IF
      END DO
C
      IF (MULT.EQ.1) THEN
      IF (CPC(I).GT.RSMALL) THEN
      SD=SQRT(KTBY2/CPC(I))/RAD
      ELSE
      SD=BIG99
      END IF
      WRITE(UNIT,'(9A,F10.2,A,F10.3,A,I4,A,F11.4)')
     & ' DIHEdral  ',KCP1(I),' ',
     &      KCP2(I),' ',KCP3(I),' ',KCP4(I),'  ',CPC(I),
     &      ' {sd=',SD,'} ',CPD(I),' ',
     &      CPB(I)/RAD
      ELSE
      II=I-MULT+1
      WRITE(UNIT,'(9A,I2,A,F10.2,A,I4,A,F11.4)')
     & ' DIHEdral  ',KCP1(II),' ',
     &   KCP2(II),' ',KCP3(II),' ',KCP4(II),'   MULTiple=',MULT,' ',
     &   CPC(II),'  ',CPD(II),' ',
     &   CPB(II)/RAD
      DO II=I-MULT+2,I
      WRITE(UNIT,'(45X,F10.2,A,I4,A,F11.4)')
     &   CPC(II),'  ',CPD(II),' ',CPB(II)/RAD
      END DO
      END IF
C
      I=I+1
      END DO
      END IF
C
      IF (NCI.GT.0) THEN
      WRITE(UNIT,'(A)')
      I=1
      DO WHILE (I.LE.NCI)
C
C determine multiplicity
      MULT=1
      COND=(I.LT.NCI)
      IF (COND) THEN
      COND=(KCI1(I+1).EQ.KCI1(I).AND.KCI2(I+1).EQ.KCI2(I).AND.
     &      KCI3(I+1).EQ.KCI3(I).AND.KCI4(I+1).EQ.KCI4(I))
      END IF
      DO WHILE (COND)
      MULT=MULT+1
      I=I+1
      IF (COND) THEN
      COND=(KCI1(I+1).EQ.KCI1(I).AND.KCI2(I+1).EQ.KCI2(I).AND.
     &      KCI3(I+1).EQ.KCI3(I).AND.KCI4(I+1).EQ.KCI4(I))
      END IF
      END DO
C
      IF (MULT.EQ.1) THEN
      IF (CIC(I).GT.RSMALL) THEN
      SD=SQRT(KTBY2/CIC(I))/RAD
      ELSE
      SD=BIG99
      END IF
      WRITE(UNIT,'(9A,F10.2,A,F10.3,A,I4,A,F11.4)')
     & ' IMPRoper  ',KCI1(I),' ',
     &      KCI2(I),' ',KCI3(I),' ',KCI4(I),'  ',CIC(I),
     &      ' {sd=',SD,'} ',CID(I),' ',
     &      CIB(I)/RAD
      ELSE
      II=I-MULT+1
      WRITE(UNIT,'(9A,I2,A,F10.2,A,I4,A,F11.4)')
     & ' IMPRoper  ',KCI1(II),' ',
     &   KCI2(II),' ',KCI3(II),' ',KCI4(II),'   MULTiple=',MULT,' ',
     &   CIC(II),'  ',CID(II),' ',
     &   CIB(II)/RAD
      DO II=I-MULT+2,I
      WRITE(UNIT,'(45X,F10.2,A,I4,A,F11.4)')
     &   CIC(II),'  ',CID(II),' ',CIB(II)/RAD
      END DO
      END IF
C
      I=I+1
      END DO
      END IF
C
      IF (NCN.GT.0) THEN
      WRITE(UNIT,'(A)')
      WRITE(UNIT,'(A)')
     1' !                  eps     sigma       eps(1:4) sigma(1:4)'
      DO I=1,NCN
      CALL LJESAB(EPS,SIG,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPS14,SIG14,CNBA14(I,I),CNBB14(I,I))
      WRITE(UNIT,'(3A,F8.4,A,F8.4,A,F8.4,A,F8.4)') ' NONBonded  ',
     &      CNAC(I),'  ',EPS,' ',SIG,'    ',EPS14,' ',SIG14
      END DO
C
      WRITE(UNIT,'(A)')
C
C only write cross-terms for NBFIxes
      DO I=1,NCN
      DO J=I+1,NCN
C
      CALL LJESAB(EPS,SIG,CNBA(I,J),CNBB(I,J))
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPSJ,SIGJ,CNBA(J,J),CNBB(J,J))
      COND=ABS(SIG-(SIGI+SIGJ)/TWO).LT.R4SMAL.AND.
     1     ABS(EPS-SQRT(EPSI*EPSJ)).LT.R4SMAL
C
      CALL LJESAB(EPS14,SIG14,CNBA14(I,J),CNBB14(I,J))
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJESAB(EPSJ14,SIGJ14,CNBA14(J,J),CNBB14(J,J))
      COND14=ABS(SIG14-(SIGI14+SIGJ14)/TWO).LT.R4SMAL.AND.
     1       ABS(EPS14-SQRT(EPSI14*EPSJ14)).LT.R4SMAL
C
      IF (.NOT.COND.OR..NOT.COND14) THEN
      WRITE(UNIT,'(5A,F12.3,A,F12.3,A,F12.3,A,F12.3)') ' NBFIx  ',
     &      CNAC(I),' ',CNAC(J),' ',CNBA(I,J),' ',CNBB(I,J),'    ',
     &      CNBA14(I,J),' ',CNBB14(I,J)
      END IF
      END DO
      END DO
C
      END IF
C
      ELSE
C verbose mode: all parameters
      WRITE(UNIT,'(A)')
     & ' Listing of all atom-based (F) and type-based (T) parameters.'
C
C
      IF (NBOND.GT.0.AND.QENER(SSBOND)) THEN
      IF (UPBOND) THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,IAC,
     &            SEGID,RESID,TYPE)
      END IF
      WRITE(UNIT,'(A)')
      DO I=1,NBOND
      WRITE(UNIT,'(15A,F10.3,L1,A,F10.3,L1)') ' BOND ',
     & '(ATOM ',SEGID(IB(I)),' ',RESID(IB(I)),' ',TYPE(IB(I)),') ',
     & '(ATOM ',SEGID(JB(I)),' ',RESID(JB(I)),' ',TYPE(JB(I)),') ',
     &  ACBC(I),QACBC(I),' ',ACBB(I),QACBB(I)
      END DO
      END IF
C
      IF (NTHETA.GT.0.AND.QENER(SSANGL)) THEN
      IF (UPANGL) THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,ACTB,ACTUC,ACTUB,
     &            QACTC,QACTB,QACTUC,QACTUB,IAC,SEGID,RESID,TYPE)
      END IF
      WRITE(UNIT,'(A)')
      DO I=1,NTHETA
      WRITE(UNIT,
     & '(15A,/28X,7A,F10.3,L1,A,F10.3,L1,A,F10.3,L1,A,F10.3,L1)')
     & ' ANGL ',
     &  '(ATOM ',SEGID(IT(I)),' ',RESID(IT(I)),' ',TYPE(IT(I)),') ',
     &  '(ATOM ',SEGID(JT(I)),' ',RESID(JT(I)),' ',TYPE(JT(I)),') ',
     &  '(ATOM ',SEGID(KT(I)),' ',RESID(KT(I)),' ',TYPE(KT(I)),') ',
     &  ACTC(I),QACTC(I),' ',ACTB(I)/RAD,QACTB(I),
     & '  UB ',ACTUC(I),QACTUC(I),'  ',ACTUB(I),QACTUB(I)
      END DO
      END IF
C
      IF (NPHI.GT.0.AND.QENER(SSDIHE)) THEN
      IF (UPDIHE) THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      END IF
      WRITE(UNIT,'(A)')
      I=1
      DO WHILE (I.LE.NPHI)
C
C determine multiplicity
      MULT=1
      COND=(I.LT.NPHI)
      IF (COND) THEN
      COND=(IP(I+1).EQ.IP(I).AND.JP(I+1).EQ.JP(I).AND.
     &      KP(I+1).EQ.KP(I).AND.LP(I+1).EQ.LP(I))
      END IF
      DO WHILE (COND)
      MULT=MULT+1
      I=I+1
      COND=(I.LT.NPHI)
      IF (COND) THEN
      COND=(IP(I+1).EQ.IP(I).AND.JP(I+1).EQ.JP(I).AND.
     &      KP(I+1).EQ.KP(I).AND.LP(I+1).EQ.LP(I))
      END IF
      END DO
      WRITE(UNIT,'(15A,/6X,15A,I4)') ' DIHE ',
     &  '(ATOM ',SEGID(IP(I)),' ',RESID(IP(I)),' ',TYPE(IP(I)),') ',
     &  '(ATOM ',SEGID(JP(I)),' ',RESID(JP(I)),' ',TYPE(JP(I)),') ',
     &  '(ATOM ',SEGID(KP(I)),' ',RESID(KP(I)),' ',TYPE(KP(I)),') ',
     &  '(ATOM ',SEGID(LP(I)),' ',RESID(LP(I)),' ',TYPE(LP(I)),') ',
     &  '      MULT',MULT
      DO II=I-MULT+1,I
      WRITE (UNIT,'(50X,F10.2,L1,I4,L1,F11.2,L1)')
     &    ACPC(II),QACPC(II),ACPD(II),QACPD(II),ACPB(II)/RAD,QACPB(II)
      END DO
      I=I+1
      END DO
      END IF
C
      IF (NIMPHI.GT.0.AND.QENER(SSIMPR)) THEN
      IF (UPIMPR) THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,QACIC,QACIB,QACID,
     &            IAC,SEGID,RESID,TYPE)
      END IF
      WRITE(UNIT,'(A)')
      I=1
      DO WHILE (I.LE.NIMPHI)
C
C determine multiplicity
      MULT=1
      COND=(I.LT.NIMPHI)
      IF (COND) THEN
      COND=(IM(I+1).EQ.IM(I).AND.JM(I+1).EQ.JM(I).AND.
     &      KM(I+1).EQ.KM(I).AND.LM(I+1).EQ.LM(I))
      END IF
      DO WHILE (COND)
      MULT=MULT+1
      I=I+1
      COND=(I.LT.NIMPHI)
      IF (COND) THEN
      COND=(IM(I+1).EQ.IM(I).AND.JM(I+1).EQ.JM(I).AND.
     &      KM(I+1).EQ.KM(I).AND.LM(I+1).EQ.LM(I))
      END IF
      END DO
      WRITE(UNIT,'(15A,/6X,15A,I4)') ' IMPR ',
     &  '(ATOM ',SEGID(IM(I)),' ',RESID(IM(I)),' ',TYPE(IM(I)),') ',
     &  '(ATOM ',SEGID(JM(I)),' ',RESID(JM(I)),' ',TYPE(JM(I)),') ',
     &  '(ATOM ',SEGID(KM(I)),' ',RESID(KM(I)),' ',TYPE(KM(I)),') ',
     &  '(ATOM ',SEGID(LM(I)),' ',RESID(LM(I)),' ',TYPE(LM(I)),') ',
     &  '      MULT',MULT
      DO II=I-MULT+1,I
      WRITE (UNIT,'(50X,F10.2,L1,I4,L1,F11.2,L1)')
     &    ACIC(II),QACIC(II),ACID(II),QACID(II),ACIB(II)/RAD,QACIB(II)
      END DO
      I=I+1
      END DO
      END IF
C
C Lennard-Jones parameters
      IF ((NCN.GT.0.OR.NLJAT.GT.0).AND.
     &              (QENER(SSVDW).OR.QENER(SSPVDW))) THEN
      CALL NBUPDA
C
      WRITE(UNIT,'(A)')
      WRITE(UNIT,'(A)')
     1' !         slot      eps     sigma       eps(1:4) sigma(1:4)'
      WRITE(UNIT,'(/A)') ' TYPE-based parameters'
      DO I=1,NCN
      CALL LJESAB(EPS,SIG,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPS14,SIG14,CNBA14(I,I),CNBB14(I,I))
      WRITE(UNIT,'(A,I5,A,F8.4,A,F8.4,A,F8.4,A,F8.4)') ' NONBonded  ',
     &      I,'  ',EPS,' ',SIG,'    ',EPS14,' ',SIG14
      END DO
C only write cross-terms for NBFIxes
      DO I=1,NCN
      DO J=I+1,NCN
C
      CALL LJESAB(EPS,SIG,CNBA(I,J),CNBB(I,J))
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPSJ,SIGJ,CNBA(J,J),CNBB(J,J))
      COND=ABS(SIG-(SIGI+SIGJ)/TWO).LT.R4SMAL.AND.
     1     ABS(EPS-SQRT(EPSI*EPSJ)).LT.R4SMAL
C
      CALL LJESAB(EPS14,SIG14,CNBA14(I,J),CNBB14(I,J))
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJESAB(EPSJ14,SIGJ14,CNBA14(J,J),CNBB14(J,J))
      COND14=ABS(SIG14-(SIGI14+SIGJ14)/TWO).LT.R4SMAL.AND.
     1       ABS(EPS14-SQRT(EPSI14*EPSJ14)).LT.R4SMAL
C
      IF (.NOT.COND.OR..NOT.COND14) THEN
      WRITE(UNIT,'(A,I5,A,I5,A,F12.3,A,F12.3,A,F12.3,A,F12.3)')
     &    ' NBFIx  ',I,' ',J,' ',CNBA(I,J),' ',CNBB(I,J),'    ',
     &      CNBA14(I,J),' ',CNBB14(I,J)
      END IF
      END DO
      END DO
C
      WRITE(UNIT,'(/A)') ' ATOM-based parameters'
      DO I=MAXCN,MAXCN-NLJAT+1,-1
      CALL LJESAB(EPS,SIG,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPS14,SIG14,CNBA14(I,I),CNBB14(I,I))
      WRITE(UNIT,'(A,I5,A,F8.4,A,F8.4,A,F8.4,A,F8.4)') ' NONBonded  ',
     &      I,'  ',EPS,' ',SIG,'    ',EPS14,' ',SIG14
      END DO
C only write cross-terms for NBFIxes
      DO I=MAXCN,MAXCN-NLJAT+1,-1
      DO J=MAXCN,I+1,-1
C
      CALL LJESAB(EPS,SIG,CNBA(I,J),CNBB(I,J))
      CALL LJESAB(EPSI,SIGI,CNBA(I,I),CNBB(I,I))
      CALL LJESAB(EPSJ,SIGJ,CNBA(J,J),CNBB(J,J))
      COND=ABS(SIG-(SIGI+SIGJ)/TWO).LT.R4SMAL.AND.
     1     ABS(EPS-SQRT(EPSI*EPSJ)).LT.R4SMAL
C
      CALL LJESAB(EPS14,SIG14,CNBA14(I,J),CNBB14(I,J))
      CALL LJESAB(EPSI14,SIGI14,CNBA14(I,I),CNBB14(I,I))
      CALL LJESAB(EPSJ14,SIGJ14,CNBA14(J,J),CNBB14(J,J))
      COND14=ABS(SIG14-(SIGI14+SIGJ14)/TWO).LT.R4SMAL.AND.
     1       ABS(EPS14-SQRT(EPSI14*EPSJ14)).LT.R4SMAL
C
      IF (.NOT.COND.OR..NOT.COND14) THEN
      WRITE(UNIT,'(A,I5,A,I5,A,F12.3,A,F12.3,A,F12.3,A,F12.3)')
     &    ' NBFIx  ',I,' ',J,' ',CNBA(I,J),' ',CNBB(I,J),'    ',
     &      CNBA14(I,J),' ',CNBB14(I,J)
      END IF
      END DO
      END DO
C
      WRITE(UNIT,'(A)')
      WRITE(UNIT,'(/A)') ' Atom-selections'
      DO I=1,NCN
      DO J=1,NATOM
      IF (LOOKUP(J).EQ.I) THEN
      WRITE(UNIT,'(7A,I5)') '   ',
     &      SEGID(J),' ',RESID(J),' ',TYPE(J),'     slot=',I
      END IF
      END DO
      END DO
      DO I=MAXCN,MAXCN-NLJAT+1,-1
      DO J=1,NATOM
      IF (LOOKUP(J).EQ.I) THEN
      WRITE(UNIT,'(7A,I5)') '   ',
     &      SEGID(J),' ',RESID(J),' ',TYPE(J),'     slot=',I
      END IF
      END DO
      END DO
C
      END IF
C
      END IF
C
      CALL VCLOSE(UNIT,'KEEP',ERROR)
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE PARINI
C
C routine initializes parameters
C
      IMPLICIT NONE
C
      INCLUDE 'param.inc'
C begin
C
C type-base
      NCB=0
      NCT=0
      NNCCP=0
      NCI=0
      NCN=0
C atom-base (nonbonded slots)
      NLJAT=0
      RETURN
      END
C======================================================================
      SUBROUTINE PARREDU(IFLAGS,IC,NC)
C
C Routine reduces atom-based parameters to type-based parameters
C
C Authors: Thomas Simonson and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'param.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'update.inc'
      INCLUDE 'funct.inc'
      INTEGER IFLAGS(*), IC(*), NC(*)
C local
      INTEGER I, IFLAG, J, ISELCT
      LOGICAL QOVER, QCONST, COND
      CHARACTER*4 SMODE, AI, AJ, AK, AL
      DOUBLE PRECISION KTBY2, OFFSET, DIFF, ACPBB, ACIBB
      CHARACTER*16 STAX, STA1
C parameters
      DOUBLE PRECISION T298, TWO, ONE, ZERO, BIG99
      INTEGER MARK
      PARAMETER (T298=298.D0, TWO=2.0D0, ZERO=0.0D0, BIG99=999999.D0)
      PARAMETER (ONE=1.0D0, MARK=-9999)
C
C defaults
      DO I=1,NATOM
      IFLAGS(I)=1
      END DO
      ISELCT=NATOM
      QOVER=.TRUE.
      QCONST=.TRUE.
C
      CALL PUSEND('REDUCE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('REDUCE>')
      CALL MISCOM('REDUCE>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-parameter-reduce')
C
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      IF (QCONST) THEN
      SMODE='STAT'
      ELSE
      SMODE='AVER'
      END IF
      WRITE(6,'(A,I6,A,L1,2A)')
     & ' Number of selected atoms=',ISELCT,
     & '   OVERwrite=',QOVER,'   MODE=',SMODE
C======================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      DO I=1,NATOM
      IFLAGS(I)=0
      END DO
      CALL SELCTA(IFLAGS,ISELCT,X,Y,Z,.TRUE.)
C======================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      IF (QCONST) THEN
      SMODE='STAT'
      ELSE
      SMODE='AVER'
      END IF
      CALL NEXTA4('MODE=',SMODE)
      IF (SMODE.EQ.'STAT') THEN
      QCONST=.TRUE.
      ELSE
      QCONST=.FALSE.
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'OVER') THEN
      CALL NEXTLO('overwrite=',QOVER)
C
      ELSE
      CALL CHKEND('REDUCE>',DONE)
      ENDIF
      ENDIF
C
      ENDDO
      DONE=.FALSE.
C
      KTBY2=KBOLTZ*T298/TWO
C
C sort the type-based parameters for the type-based
C lookup (calls to CODBON, CODANG, CODDIHE, CODIMPR)
      CALL PARSORT
C
C now do the actual reduction of atomic parameters to
C chemical parameters
C
C BONDS
C ------
      IF (QENER(SSBOND)) THEN
C update parameter lookup lists if required
      IF (UPBOND) THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,
     &            IAC,SEGID,RESID,TYPE)
      END IF
C
C initialize chemical arrays
      IFLAG=0
      DO I=1,MAXCB
      NC(I)=0
      END DO
C
C Lookup `chemical-based' parameters in data base and initialize.
C If a new chemical bond is needed, add to lists.
C If QOVER is not set, only initialize the new bond types.
      DO I=1,NBOND
      IF (IFLAGS(IB(I)).GT.0.AND.IFLAGS(JB(I)).GT.0) THEN
      J=1
C
      AI=IAC(IB(I))
      AJ=IAC(JB(I))
      IF (.NOT.LTSTEQ(AI,4,AJ,4,.TRUE.)) THEN
      AJ=IAC(IB(I))
      AI=IAC(JB(I))
C
      END IF
      DO WHILE (J.LE.NCB.AND..NOT.(KCB1(J).EQ.AI.AND.KCB2(J).EQ.AJ))
      J=J+1
      END DO
      IC(I)=J
      IF (IC(I).GT.NCB) THEN
      NCB=NCB+1
      IF (NCB.GT.MAXCB)
     &  CALL WRNDIE(-1,'REDUCE','MAXCB exceeded.')
      IC(I)=NCB
      KCB1(NCB)=AI
      KCB2(NCB)=AJ
      CBB(IC(I))=ZERO
      CBC(IC(I))=ZERO
      IFLAG=IFLAG+1
      ELSE IF (QOVER) THEN
      CBB(IC(I))=ZERO
      CBC(IC(I))=ZERO
      ENDIF
      ENDIF
      END DO
      IF (IFLAG.GT.0) WRITE(6,'(A,I4,A)')
     &  ' REDUCE: creating ',IFLAG,' new bond parameters'
C
C accumulate chemical arrays
      DO I=1,NBOND
      IF (IFLAGS(IB(I)).GT.0.AND.IFLAGS(JB(I)).GT.0) THEN
      IF (IC(I).GT.NCB-IFLAG.OR.QOVER) THEN
      NC(IC(I))=NC(IC(I))+1
      CBB(IC(I))=CBB(IC(I))+ACBB(I)
      IF (QCONST) THEN
C for the statistics option accumulate the square
      CBC(IC(I))=CBC(IC(I))+ACBB(I)**2
      ELSE
C for the average option accumulate the average of the energy constants
      CBC(IC(I))=CBC(IC(I))+ACBC(I)
      END IF
      ENDIF
      ENDIF
      END DO
C
C calculate average values
      DO I=1,NCB
      IF (NC(I).GT.0) THEN
      CBB(I)=CBB(I)/NC(I)
      CBC(I)=CBC(I)/NC(I)
C
      IF (QCONST) THEN
C for the statistics option compute the variance
      IF (ABS(CBC(I)-CBB(I)**2).LT.RSMALL) THEN
      CBC(I)=BIG99
      ELSE
      CBC(I)=KTBY2/(CBC(I)-CBB(I)**2)
      END IF
      END IF
C
      END IF
      END DO
C
      ENDIF
C
C ANGLES
C -------
      IF (QENER(SSANGL)) THEN
C update parameter lookup lists if required
      IF (UPANGL) THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,ACTB,ACTUC,ACTUB,
     &            QACTC,QACTB,QACTUC,QACTUB,IAC,SEGID,RESID,TYPE)
      END IF
C
      IFLAG=0
      DO I=1,MAXCT
      NC(I)=0
      END DO
C
C Lookup `chemical-based' parameters in data base and initialize.
C If a new chemical angle is needed, add to lists.
      DO I=1,NTHETA
      IF (IFLAGS(IT(I)).GT.0.AND.IFLAGS(JT(I)).GT.0.AND.
     &    IFLAGS(KT(I)).GT.0) THEN
      J=1
C
      AI=IAC(IT(I))
      AJ=IAC(JT(I))
      AK=IAC(KT(I))
      IF (.NOT.LTSTEQ(AI,4,AK,4,.TRUE.)) THEN
      AK=IAC(IT(I))
      AI=IAC(KT(I))
      END IF
C
      DO WHILE (J.LE.NCT.AND..NOT.
     &          (AI.EQ.KCT1(J).AND.AJ.EQ.KCT2(J).AND.AK.EQ.KCT3(J)))
      J=J+1
      END DO
      IC(I)=J
      IF (IC(I).GT.NCT) THEN
      NCT=NCT+1
      IF (NCT.GT.MAXCT)
     &  CALL WRNDIE(-1,'REDUCE','MAXCT exceeded.')
      IC(I)=NCT
      KCT1(NCT)=AI
      KCT2(NCT)=AJ
      KCT3(NCT)=AK
      CTB(IC(I))=ZERO
      CTC(IC(I))=ZERO
      CTUB(IC(I))=ZERO
      CTUC(IC(I))=ZERO
      IFLAG=IFLAG+1
      ELSE IF (QOVER) THEN
      CTB(IC(I))=ZERO
      CTC(IC(I))=ZERO
      CTUB(IC(I))=ZERO
      CTUC(IC(I))=ZERO
      ENDIF
      ENDIF
      END DO
      IF (IFLAG.GT.0) WRITE(6,'(A,I4,A)')
     &  ' REDUCE: creating ',IFLAG,' new angle parameters'
C
C accumulate chemical arrays.
      DO I=1,NTHETA
      IF (IFLAGS(IT(I)).GT.0.AND.IFLAGS(JT(I)).GT.0.AND.
     &    IFLAGS(KT(I)).GT.0) THEN
      IF (IC(I).GT.NCT-IFLAG.OR.QOVER) THEN
      NC(IC(I))=NC(IC(I))+1
      CTB(IC(I))=CTB(IC(I))+ACTB(I)
      CTUB(IC(I))=CTUB(IC(I))+ACTUB(I)
      IF (QCONST) THEN
C for the statistics option accumulate the square
      CTC(IC(I))=CTC(IC(I))+ACTB(I)**2
      CTUC(IC(I))=CTUC(IC(I))+ACTUB(I)**2
      ELSE
C for the average option accumulate the average of the energy constants
      CTC(IC(I))=CTC(IC(I))+ACTC(I)
      CTUC(IC(I))=CTUC(IC(I))+ACTUC(I)
      END IF
      ENDIF
      ENDIF
      END DO
C
C calculate average values
      DO I=1,NCT
      IF (NC(I).GT.0) THEN
      CTB(I)=CTB(I)/NC(I)
      CTC(I)=CTC(I)/NC(I)
      CTUB(I)=CTUB(I)/NC(I)
      CTUC(I)=CTUC(I)/NC(I)
C
      IF (QCONST) THEN
C for the statistics option compute the variance
      IF (ABS(CTC(I)-CTB(I)**2).LT.RSMALL) THEN
      CTC(I)=BIG99
      ELSE
      CTC(I)=KTBY2/(CTC(I)-CTB(I)**2)
      END IF
C only compute standard deviation for UB term if it is actually used
      IF (ABS(CTUB(I)).GT.RSMALL.AND.ABS(CTUC(I)).GT.RSMALL) THEN
      IF (ABS(CTUC(I)-CTUB(I)**2).LT.RSMALL) THEN
      CTUC(I)=BIG99
      ELSE
      CTUC(I)=KTBY2/(CTUC(I)-CTUB(I)**2)
      END IF
      END IF
C
      END IF
C
      END IF
      END DO
C
      ENDIF
C
C DIHEDRALS
C ----------
      IF (QENER(SSDIHE)) THEN
C
C update parameter lookup lists if required
      IF (UPDIHE) THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      END IF
C
      IFLAG=0
      DO I=1,MAXCP
      NC(I)=0
      END DO
C
C Lookup `chemical-based' parameters in data base and initialize.
C If a new chemical DIHEDRAL is needed, add to lists.
      DO I=1,NPHI
      IF (IFLAGS(IP(I)).GT.0.AND.IFLAGS(JP(I)).GT.0.AND.
     &    IFLAGS(KP(I)).GT.0.AND.IFLAGS(LP(I)).GT.0) THEN
      COND=(I.GT.1)
      IF (COND) THEN
      COND=(IP(I-1).EQ.IP(I).AND.JP(I-1).EQ.JP(I).AND.
     &      KP(I-1).EQ.KP(I).AND.LP(I-1).EQ.LP(I))
      END IF
      IF (COND) THEN
      OFFSET=OFFSET+1
      ELSE
      OFFSET=1
      ENDIF
      J=1
C
      AI=IAC(IP(I))
      AJ=IAC(JP(I))
      AK=IAC(KP(I))
      AL=IAC(LP(I))
      IF (.NOT.LTSTEQ(AJ//AI,8,AK//AL,8,.TRUE.)) THEN
      AL=IAC(IP(I))
      AK=IAC(JP(I))
      AJ=IAC(KP(I))
      AI=IAC(LP(I))
      END IF
C
      STA1=AI//AJ//AK//AL
      STAX=KCP1(J)//KCP2(J)//KCP3(J)//KCP4(J)
      DO WHILE (J.LE.NNCCP.AND.STAX.NE.STA1)
      J=J+1
      STAX=KCP1(J)//KCP2(J)//KCP3(J)//KCP4(J)
      END DO
      IC(I)=J+OFFSET-1
      IF (IC(I).LE.NNCCP.AND.QOVER) THEN
      CPB(IC(I))=ZERO
      CPC(IC(I))=ZERO
      CPD(IC(I))=MARK
      ENDIF
C
      IF (IC(I).GT.NNCCP) THEN
      NNCCP=NNCCP+1
      IF (NNCCP.GT.MAXCP)
     &  CALL WRNDIE(-1,'REDUCE','MAXCP exceeded.')
      IC(I)=NNCCP
      KCP1(NNCCP)=AI
      KCP2(NNCCP)=AJ
      KCP3(NNCCP)=AK
      KCP4(NNCCP)=AL
      CPB(IC(I))=ZERO
      CPC(IC(I))=ZERO
      CPD(IC(I))=MARK
      IFLAG=IFLAG+1
      ENDIF
      ENDIF
      END DO
      IF (IFLAG.GT.0) WRITE(6,'(A,I4,A)')
     &  ' REDUCE: creating ',IFLAG,' new dihedral parameters'
C
C accumulate chemical arrays
      DO I=1,NPHI
      IF (IFLAGS(IP(I)).GT.0.AND.IFLAGS(JP(I)).GT.0.AND.
     &    IFLAGS(KP(I)).GT.0.AND.IFLAGS(LP(I)).GT.0) THEN
      IF (IC(I).GT.NNCCP-IFLAG.OR.QOVER) THEN
      NC(IC(I))=NC(IC(I))+1
C
C note the periodicity of the angle
      IF (NC(IC(I)).LE.1) THEN
      DIFF=-ACPB(I)
      ELSE
      DIFF=CPB(IC(I))/(NC(IC(I))-1)-ACPB(I)
      END IF
      IF (ABS(DIFF).GT.PI) THEN
      ACPBB=ACPB(I)+TWO*PI*SIGN(ONE,DIFF)
      ELSE
      ACPBB=ACPB(I)
      END IF
C
C average the angle
      CPB(IC(I))=CPB(IC(I))+ACPBB
C
C make sure the peridicities are compatible
      IF (CPD(IC(I)).EQ.MARK) THEN
      CPD(IC(I))=ACPD(I)
      ELSE IF (CPD(IC(I)).NE.ACPD(I)) THEN
CCC modification ATB 4/27/08
      WRITE(6,'(8A)')
     & ' %REDUce-ERR: incompatible dihedral periodicities for types ',
     & IAC(IP(I)),' ',IAC(JP(I)),' ',IAC(KP(I)),' ',IAC(LP(I))
      END IF
C
C treatment of the energy constant
      IF (QCONST) THEN
C for the statistics option accumulate the square
      CPC(IC(I))=CPC(IC(I))+ACPBB**2
      ELSE
C for the average option accumulate the average of the energy constants
      CPC(IC(I))=CPC(IC(I))+ACPC(I)
      END IF
C
      ENDIF
      ENDIF
      END DO
C
C calculate average values
      DO I=1,NNCCP
      IF (NC(I).GT.0) THEN
      CPB(I)=CPB(I)/NC(I)
      CPC(I)=CPC(I)/NC(I)
C
      IF (QCONST) THEN
C for the statistics option compute the variance
      IF (ABS(CPC(I)-CPB(I)**2).LT.RSMALL) THEN
      CPC(I)=BIG99
      ELSE
      CPC(I)=KTBY2/(CPC(I)-CPB(I)**2)
      END IF
      END IF
C
      END IF
      END DO
C
      ENDIF
C
C IMPROPERS
C ----------
      IF (QENER(SSIMPR)) THEN
C
C update parameter lookup lists if required
      IF (UPIMPR) THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,QACIC,QACIB,QACID,
     &            IAC,SEGID,RESID,TYPE)
      END IF
C
      IFLAG=0
      DO I=1,MAXCI
      NC(I)=0
      END DO
C
C initialize chemical arrays. If a new chemical improper is needed, add to lists.
      DO I=1,NIMPHI
      IF (IFLAGS(IM(I)).GT.0.AND.IFLAGS(JM(I)).GT.0.AND.
     &    IFLAGS(KM(I)).GT.0.AND.IFLAGS(LM(I)).GT.0) THEN
C multiple improper condition
      COND=(I.GT.1)
      IF (COND) THEN
      COND=(IM(I-1).EQ.IM(I).AND.JM(I-1).EQ.JM(I).AND.
     &      KM(I-1).EQ.KM(I).AND.LM(I-1).EQ.LM(I))
      END IF
      IF (COND) THEN
      OFFSET=OFFSET+1
      ELSE
      OFFSET=1
      ENDIF
C
      J=1
C
      AI=IAC(IM(I))
      AJ=IAC(JM(I))
      AK=IAC(KM(I))
      AL=IAC(LM(I))
      IF (.NOT.LTSTEQ(AI//AJ,8,AL//AK,8,.TRUE.)) THEN
      AL=IAC(IM(I))
      AK=IAC(JM(I))
      AJ=IAC(KM(I))
      AI=IAC(LM(I))
      END IF
C
      STA1=AI//AJ//AK//AL
      STAX=KCI1(J)//KCI2(J)//KCI3(J)//KCI4(J)
      DO WHILE (J.LE.NCI.AND.STAX.NE.STA1)
      J=J+1
      STAX=KCI1(J)//KCI2(J)//KCI3(J)//KCI4(J)
      END DO
      IC(I)=J+OFFSET-1
      IF (IC(I).LE.NCI.AND.QOVER) THEN
      CIB(IC(I))=ZERO
      CIC(IC(I))=ZERO
      CID(IC(I))=MARK
      ENDIF
      IF (IC(I).GT.NCI) THEN
      NCI=NCI+1
      IF (NCI.GT.MAXCI)
     &  CALL WRNDIE(-1,'REDUCE','MAXCI exceeded.')
      IC(I)=NCI
      KCI1(NCI)=AI
      KCI2(NCI)=AJ
      KCI3(NCI)=AK
      KCI4(NCI)=AL
      CIB(IC(I))=ZERO
      CIC(IC(I))=ZERO
      CID(IC(I))=MARK
      IFLAG=IFLAG+1
      END IF
      END IF
      END DO
      IF (IFLAG.GT.0) WRITE(6,'(A,I4,A)')
     &  ' REDUCE: creating ',IFLAG,' new improper parameters'
C
C accumulate chemical arrays
      DO I=1,NIMPHI
      IF (IFLAGS(IM(I)).GT.0.AND.IFLAGS(JM(I)).GT.0.AND.
     &    IFLAGS(KM(I)).GT.0.AND.IFLAGS(LM(I)).GT.0) THEN
      IF (IC(I).GT.NCI-IFLAG.OR.QOVER) THEN
      NC(IC(I))=NC(IC(I))+1
C note the periodicity of the angle
      IF (NC(IC(I)).LE.1) THEN
      DIFF=-ACIB(I)
      ELSE
      DIFF=CIB(IC(I))/(NC(IC(I))-1)-ACIB(I)
      END IF
      IF (ABS(DIFF).GT.PI) THEN
      ACIBB=ACIB(I)+TWO*PI*SIGN(ONE,DIFF)
      ELSE
      ACIBB=ACIB(I)
      END IF
C
C average the angle
      CIB(IC(I))=CIB(IC(I))+ACIBB
C
C make sure the peridicities are compatible
      IF (CID(IC(I)).EQ.MARK) THEN
      CID(IC(I))=ACID(I)
      ELSE IF (CID(IC(I)).NE.ACID(I)) THEN
CCC modification ATB 4/27/08
      WRITE(6,'(8A)')
     & ' %REDUce-ERR: incompatible improper periodicities for types ',
     & IAC(IM(I)),' ',IAC(JM(I)),' ',IAC(KM(I)),' ',IAC(LM(I))
      END IF
C
C treatment of the energy constant
      IF (QCONST) THEN
C for the statistics option accumulate the square
      CIC(IC(I))=CIC(IC(I))+ACIBB**2
      ELSE
C for the average option accumulate the average of the energy constants
      CIC(IC(I))=CIC(IC(I))+ACIC(I)
      END IF
C
      ENDIF
      ENDIF
      END DO
C
C calculate average values
      DO I=1,NCI
      IF (NC(I).GT.0) THEN
      CIB(I)=CIB(I)/NC(I)
      CIC(I)=CIC(I)/NC(I)
C
      IF (QCONST) THEN
C for the statistics option compute the variance
      IF (ABS(CIC(I)-CIB(I)**2).LT.RSMALL) THEN
      CIC(I)=BIG99
      ELSE
      CIC(I)=KTBY2/(CIC(I)-CIB(I)**2)
      END IF
      END IF
C
      END IF
      END DO
C
      ENDIF
C
      RETURN
      END
C======================================================================
      SUBROUTINE PARLEAR
C
C learns parameters from cartesian coordinates
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/ouput
      INCLUDE 'cns.inc'
      INCLUDE 'learn.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C local
      CHARACTER*4 SMODE
C begin
C
      CALL NEXTWD('LEARn>')
      IF (WD(1:4).EQ.'HELP') THEN
      CALL CNSHELP('cns-parameter-learn')
      ELSE IF (WD(1:4).EQ.'INIT') THEN
      IF (PRINIT) THEN
C
C the initialization flag is set ===> set the defaults
      PRINIT=.FALSE.
C defaults
      IF (NATOM.GT.0) THEN
      LPRIN=NATOM
      PRIND=ALLHP(INTEG4(LPRIN))
      CALL FILL4(HEAP(PRIND),NATOM,1)
      PRSELCT=NATOM
      END IF
C
      PQCONST=.TRUE.
      ELSE
      CALL WRNDIE(-5,'PARLEAR','INITialization already performed')
      END IF
C parsing
C
      CALL PUSEND('LEARn-INITialize>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LEARn-INITialize>')
      CALL MISCOM('LEARn-INITialize>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'?   ') THEN
      IF (PQCONST) THEN
      SMODE='STAT'
      ELSE
      SMODE='NOST'
      END IF
      WRITE(6,'(A,I6,2A)')
     & ' Number of selected atoms=',PRSELCT,
     & '    MODE=',SMODE
C======================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(HEAP(PRIND),PRSELCT,X,Y,Z,.TRUE.)
C======================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      IF (PQCONST) THEN
      SMODE='STAT'
      ELSE
      SMODE='NOST'
      END IF
      CALL NEXTA4('MODE=',SMODE)
      IF (SMODE.EQ.'STAT') THEN
      PQCONST=.TRUE.
      ELSE
      PQCONST=.FALSE.
      END IF
C======================================================================
      ELSE
      CALL CHKEND('LEARn-INITialize>',DONE)
      ENDIF
      ENDIF
C
      ENDDO
      DONE=.FALSE.
      CALL PARLEA2(PRSELCT,HEAP(PRIND),PQCONST,PNFRAM,'INIT')
C======================================================================
      ELSE IF (WD(1:4).EQ.'ACCU') THEN
      CALL PUSEND('LEARn-ACCUmulate>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LEARn-ACCUmulate>')
      CALL MISCOM('LEARn-ACCUmulate>',USED)
      IF (.NOT.USED) THEN
      CALL CHKEND('LEARn-ACCUmulate>',DONE)
      ENDIF
C
      ENDDO
      DONE=.FALSE.
      IF (PRINIT) THEN
      CALL WRNDIE(-5,'PARLEAR','INITialization missing')
      ELSE
      CALL PARLEA2(PRSELCT,HEAP(PRIND),PQCONST,PNFRAM,'ACCU')
      END IF
C
C======================================================================
      ELSE IF (WD(1:4).EQ.'TERM') THEN
      CALL PUSEND('LEARn-TERMinate>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LEARn-TERMinate>')
      CALL MISCOM('LEARn-TERMinate>',USED)
      IF (.NOT.USED) THEN
      CALL CHKEND('LEARn-TERMinate>',DONE)
      ENDIF
C
      ENDDO
      DONE=.FALSE.
C
      IF (PRINIT) THEN
      CALL WRNDIE(-5,'PARLEAR','INITialization missing')
      ELSE
      CALL PARLEA2(PRSELCT,HEAP(PRIND),PQCONST,PNFRAM,'TERM')
      END IF
C
C free-up the heap space
      CALL PRTFRE
      CALL PRTINI
C
      ELSE
      CALL DSPERR('PARAmeter LEARn','unknown qualifier')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE PRTINI
C
C initializes the parameter learning facility
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'learn.inc'
C begin
C
C set the initialization flag
      PRINIT=.TRUE.
      PRIND=0
      RETURN
      END
C======================================================================
      SUBROUTINE PRTFRE
C
C frees space for the parameter learning facility if required
C and resets the facility.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'funct.inc'
      INCLUDE 'learn.inc'
C begin
      IF (PRIND.NE.0) CALL FREHP(PRIND,INTEG4(LPRIN))
      PRINIT=.TRUE.
      PRIND=0
      RETURN
      END
C======================================================================
      SUBROUTINE PARLEA2(ISELCT,IFLAGS,PQCONST,PNFRAM,ACTION)
C
C parameter learning facility
C this routine actually does the job
C
C Authors: Thomas Simonson and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INTEGER ISELCT, IFLAGS(*)
      LOGICAL PQCONST
      CHARACTER*4 ACTION
      INTEGER PNFRAM
C local
      DOUBLE PRECISION RIJX, RIJY, RIJZ, RIJ, RJKX, RJKY, RJKZ, RJK
      DOUBLE PRECISION RKLX, RKLY, RKLZ, RIJ2
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION RAR, RBR, RCR, CP, SP, PHI
      DOUBLE PRECISION KTBY2, DIFF
      INTEGER I
      LOGICAL ERR
C parameter
      DOUBLE PRECISION T298, TWO, ONE, BIG99, MCONST, ZERO
      PARAMETER (T298=298.D0, TWO=2.0D0, ONE=1.0D0, BIG99=999999.D0)
      PARAMETER (MCONST=0.0001D0, ZERO=0.0D0)
C begin
C
      KTBY2=KBOLTZ*T298/TWO
C
C now execute selected action
      IF (ACTION.EQ.'INIT') THEN
C
      PNFRAM=0
      IF (ISELCT.GT.0) THEN
C
      IF (QENER(SSBOND)) THEN
      DO I=1,NBOND
      IF (IFLAGS(IB(I)).GT.0.AND.IFLAGS(JB(I)).GT.0) THEN
      ACBB(I)=0.
      QACBB(I)=.FALSE.
      IF (PQCONST) THEN
      ACBC(I)=0.
      QACBC(I)=.FALSE.
      END IF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSANGL)) THEN
      DO I=1,NTHETA
      IF (IFLAGS(IT(I)).GT.0.AND.IFLAGS(JT(I)).GT.0.AND.
     &    IFLAGS(KT(I)).GT.0) THEN
      ACTB(I)=0.
      QACTB(I)=.FALSE.
C note: we do not learn UB terms but set them to zero!
      QACTUB(I)=.FALSE.
      QACTUC(I)=.FALSE.
      ACTUB(I)=ZERO
      ACTUC(I)=ZERO
      IF (PQCONST) THEN
      ACTC(I)=0.
      QACTC(I)=.FALSE.
      END IF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSDIHE)) THEN
      DO I=1,NPHI
      IF (IFLAGS(IP(I)).GT.0.AND.IFLAGS(JP(I)).GT.0.AND.
     &    IFLAGS(KP(I)).GT.0.AND.IFLAGS(LP(I)).GT.0) THEN
      ACPB(I)=0.
      QACPB(I)=.FALSE.
      IF (PQCONST) THEN
      ACPC(I)=0.
      QACPC(I)=.FALSE.
      END IF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSIMPR)) THEN
      DO I=1,NIMPHI
      IF (IFLAGS(IM(I)).GT.0.AND.IFLAGS(JM(I)).GT.0.AND.
     &    IFLAGS(KM(I)).GT.0.AND.IFLAGS(LM(I)).GT.0) THEN
      ACIB(I)=0.
      QACIB(I)=.FALSE.
      IF (PQCONST) THEN
      ACIC(I)=0.
      QACIC(I)=.FALSE.
      END IF
      ENDIF
      END DO
      ENDIF
C
      ENDIF
C
      ELSE IF (ACTION.EQ.'ACCU') THEN
      PNFRAM=PNFRAM+1
      IF (ISELCT.GT.0) THEN
C
      ERR=.FALSE.
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z).AND.IFLAGS(I).GT.0) THEN
      ERR=.TRUE.
      WRITE(6,'(9A)')
     &  ' %PARLEA2-ERR: unknown coordinates for atom "',
     &    SEGID(I),'-',RESID(I),'-',RES(I),'-',TYPE(I),'"'
      END IF
      END DO
      IF (ERR) THEN
      CALL WRNDIE(-5,'PARLEA2','Unknown coordinates.')
      ELSE
C
      IF (QENER(SSBOND)) THEN
      DO I=1,NBOND
      IF (IFLAGS(IB(I)).GT.0.AND.IFLAGS(JB(I)).GT.0) THEN
      RIJ2=(X(IB(I))-X(JB(I)))**2+(Y(IB(I))-Y(JB(I)))**2+
     &         (Z(IB(I))-Z(JB(I)))**2
      ACBB(I)=ACBB(I)+SQRT(RIJ2)
      IF (PQCONST) ACBC(I)=ACBC(I)+RIJ2
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSANGL)) THEN
      DO I=1,NTHETA
      IF (IFLAGS(IT(I)).GT.0.AND.IFLAGS(JT(I)).GT.0.AND.
     &    IFLAGS(KT(I)).GT.0) THEN
      RIJX=X(IT(I))-X(JT(I))
      RIJY=Y(IT(I))-Y(JT(I))
      RIJZ=Z(IT(I))-Z(JT(I))
      RJKX=X(KT(I))-X(JT(I))
      RJKY=Y(KT(I))-Y(JT(I))
      RJKZ=Z(KT(I))-Z(JT(I))
C compute the norm of RIJ, RKJ and set to MCONST if it is too small.
      RIJ=ONE/SQRT(MAX(MCONST,RIJX**2+RIJY**2+RIJZ**2))
      RJK=ONE/SQRT(MAX(MCONST,RJKX**2+RJKY**2+RJKZ**2))
C normalize RIJ, RKJ
      RIJX=RIJX*RIJ
      RIJY=RIJY*RIJ
      RIJZ=RIJZ*RIJ
      RJKX=RJKX*RJK
      RJKY=RJKY*RJK
      RJKZ=RJKZ*RJK
C compute truncated CP = cos( phi )
      CP=RIJX*RJKX+RIJY*RJKY+RIJZ*RJKZ
      CP=SIGN(1.D0,CP)*MIN(ABS(CP),0.99D0)
C compute PHI ( make sure  CP within boundaries )
      PHI=ACOS(MIN(ONE,MAX(-ONE,CP)))
      ACTB(I)=ACTB(I)+PHI
      IF (PQCONST) THEN
      ACTC(I)=ACTC(I)+PHI**2
      END IF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSDIHE)) THEN
      DO I=1,NPHI
      IF (IFLAGS(IP(I)).GT.0.AND.IFLAGS(JP(I)).GT.0.AND.
     &    IFLAGS(KP(I)).GT.0.AND.IFLAGS(LP(I)).GT.0) THEN
C first compute differences RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
      RIJX=X(IP(I))-X(JP(I))
      RIJY=Y(IP(I))-Y(JP(I))
      RIJZ=Z(IP(I))-Z(JP(I))
      RJKX=X(JP(I))-X(KP(I))
      RJKY=Y(JP(I))-Y(KP(I))
      RJKZ=Z(JP(I))-Z(KP(I))
      RKLX=X(KP(I))-X(LP(I))
      RKLY=Y(KP(I))-Y(LP(I))
      RKLZ=Z(KP(I))-Z(LP(I))
C now compute A= RIJ X RJK, B= RJK X RKL, C= RJK X ( RIJ X RJK )
      AX=RIJY*RJKZ-RIJZ*RJKY
      AY=RIJZ*RJKX-RIJX*RJKZ
      AZ=RIJX*RJKY-RIJY*RJKX
      BX=RJKY*RKLZ-RKLY*RJKZ
      BY=RJKZ*RKLX-RKLZ*RJKX
      BZ=RJKX*RKLY-RKLX*RJKY
      CX=RJKY*AZ-RJKZ*AY
      CY=RJKZ*AX-RJKX*AZ
      CZ=RJKX*AY-RJKY*AX
C compute the norm of A, B, C and set to MCONST if too small
      RAR=ONE/SQRT(MAX(MCONST,AX*AX+AY*AY+AZ*AZ))
      RBR=ONE/SQRT(MAX(MCONST,BX*BX+BY*BY+BZ*BZ))
      RCR=ONE/SQRT(MAX(MCONST,CX*CX+CY*CY+CZ*CZ))
C normalize A, B and C
      AX=AX*RAR
      AY=AY*RAR
      AZ=AZ*RAR
      BX=BX*RBR
      BY=BY*RBR
      BZ=BZ*RBR
      CX=CX*RCR
      CY=CY*RCR
      CZ=CZ*RCR
C compute CP= cos ( phi )  and SP= sin ( phi )
      CP=AX*BX+AY*BY+AZ*BZ
      SP=CX*BX+CY*BY+CZ*BZ
C compute PHI in degrees(make sure CP within boundaries and get sign from SP )
      PHI=-SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
C note the periodicity
      IF (PNFRAM.LE.1) THEN
      DIFF=-PHI
      ELSE
      DIFF=ACPB(I)/(PNFRAM-1)-PHI
      END IF
      IF (ABS(DIFF).GT.PI) THEN
      PHI=PHI+TWO*PI*SIGN(ONE,DIFF)
      END IF
C
      ACPB(I)=ACPB(I)+PHI
      IF (PQCONST) ACPC(I)=ACPC(I)+PHI**2
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSIMPR)) THEN
      DO I=1,NIMPHI
      IF (IFLAGS(IM(I)).GT.0.AND.IFLAGS(JM(I)).GT.0.AND.
     &    IFLAGS(KM(I)).GT.0.AND.IFLAGS(LM(I)).GT.0) THEN
C first compute differences RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
      RIJX=X(IM(I))-X(JM(I))
      RIJY=Y(IM(I))-Y(JM(I))
      RIJZ=Z(IM(I))-Z(JM(I))
      RJKX=X(JM(I))-X(KM(I))
      RJKY=Y(JM(I))-Y(KM(I))
      RJKZ=Z(JM(I))-Z(KM(I))
      RKLX=X(KM(I))-X(LM(I))
      RKLY=Y(KM(I))-Y(LM(I))
      RKLZ=Z(KM(I))-Z(LM(I))
C now compute A= RIJ X RJK, B= RJK X RKL, C= RJK X ( RIJ X RJK )
      AX=RIJY*RJKZ-RIJZ*RJKY
      AY=RIJZ*RJKX-RIJX*RJKZ
      AZ=RIJX*RJKY-RIJY*RJKX
      BX=RJKY*RKLZ-RKLY*RJKZ
      BY=RJKZ*RKLX-RKLZ*RJKX
      BZ=RJKX*RKLY-RKLX*RJKY
      CX=RJKY*AZ-RJKZ*AY
      CY=RJKZ*AX-RJKX*AZ
      CZ=RJKX*AY-RJKY*AX
C compute the norm of A, B, C and set to MCONST if too small
      RAR=ONE/SQRT(MAX(MCONST,AX*AX+AY*AY+AZ*AZ))
      RBR=ONE/SQRT(MAX(MCONST,BX*BX+BY*BY+BZ*BZ))
      RCR=ONE/SQRT(MAX(MCONST,CX*CX+CY*CY+CZ*CZ))
C normalize A, B and C
      AX=AX*RAR
      AY=AY*RAR
      AZ=AZ*RAR
      BX=BX*RBR
      BY=BY*RBR
      BZ=BZ*RBR
      CX=CX*RCR
      CY=CY*RCR
      CZ=CZ*RCR
C compute CP= cos ( phi )  and SP= sin ( phi )
      CP=AX*BX+AY*BY+AZ*BZ
      SP=CX*BX+CY*BY+CZ*BZ
C compute PHI ( make sure CP within boundaries and get sign from SP )
      PHI=-SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
C note the periodicity
      IF (PNFRAM.LE.1) THEN
      DIFF=-PHI
      ELSE
      DIFF=ACIB(I)/(PNFRAM-1)-PHI
      END IF
      IF (ABS(DIFF).GT.PI) THEN
      PHI=PHI+TWO*PI*SIGN(ONE,DIFF)
      END IF
C
      ACIB(I)=ACIB(I)+PHI
      IF (PQCONST) ACIC(I)=ACIC(I)+PHI**2
      ENDIF
      END DO
      ENDIF
C
      ENDIF
C
      ENDIF
      ELSE IF (ACTION.EQ.'TERM') THEN
      IF (ISELCT.GT.0) THEN
C
      IF (QENER(SSBOND)) THEN
      DO I=1,NBOND
      IF (IFLAGS(IB(I)).GT.0.AND.IFLAGS(JB(I)).GT.0) THEN
      ACBB(I)=ACBB(I)/PNFRAM
      IF (PQCONST) THEN
      ACBC(I)=ACBC(I)/PNFRAM
      IF (ABS(ACBC(I)-ACBB(I)**2).LT.RSMALL) THEN
      ACBC(I)=BIG99
      ELSE
      ACBC(I)=KTBY2/(ACBC(I)-ACBB(I)**2)
      ENDIF
      ENDIF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSANGL)) THEN
      DO I=1,NTHETA
      IF (IFLAGS(IT(I)).GT.0.AND.IFLAGS(JT(I)).GT.0.AND.
     &    IFLAGS(KT(I)).GT.0) THEN
      ACTB(I)=ACTB(I)/PNFRAM
      IF (PQCONST) THEN
      ACTC(I)=ACTC(I)/PNFRAM
      IF (ABS(ACTC(I)-ACTB(I)**2).LT.RSMALL) THEN
      ACTC(I)=BIG99
      ELSE
      ACTC(I)=KTBY2/(ACTC(I)-ACTB(I)**2)
      ENDIF
      ENDIF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSDIHE)) THEN
      DO I=1,NPHI
      IF (IFLAGS(IP(I)).GT.0.AND.IFLAGS(JP(I)).GT.0.AND.
     &    IFLAGS(KP(I)).GT.0.AND.IFLAGS(LP(I)).GT.0) THEN
      ACPB(I)=ACPB(I)/PNFRAM
C set periodicity to zero.
      ACPD(I)=0
      QACPD(I)=.FALSE.
      IF (PQCONST) THEN
      ACPC(I)=ACPC(I)/PNFRAM
      IF (ABS(ACPC(I)-ACPB(I)**2).LT.RSMALL) THEN
      ACPC(I)=BIG99
      ELSE
      ACPC(I)=KTBY2/(ACPC(I)-ACPB(I)**2)
      ENDIF
      ENDIF
      ENDIF
      END DO
      ENDIF
C
      IF (QENER(SSIMPR)) THEN
      DO I=1,NIMPHI
      IF (IFLAGS(IM(I)).GT.0.AND.IFLAGS(JM(I)).GT.0.AND.
     &    IFLAGS(KM(I)).GT.0.AND.IFLAGS(LM(I)).GT.0) THEN
      ACIB(I)=ACIB(I)/PNFRAM
C set periodicity to zero.
      ACID(I)=0
      QACID(I)=.FALSE.
      IF (PQCONST) THEN
      ACIC(I)=ACIC(I)/PNFRAM
      IF (ABS(ACIC(I)-ACIB(I)**2).LT.RSMALL) THEN
      ACIC(I)=BIG99
      ELSE
      ACIC(I)=KTBY2/(ACIC(I)-ACIB(I)**2)
      ENDIF
      ENDIF
      ENDIF
      END DO
      ENDIF
C
      ENDIF
      ENDIF
      RETURN
      END
C
