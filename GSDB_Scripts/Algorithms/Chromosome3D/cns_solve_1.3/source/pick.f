      SUBROUTINE PRPICK(SOBJ,THRESH)
C
C Print "quick" information about angles, bonds, dihedrals and impropers
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'param.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'update.inc'
      INCLUDE 'heap.inc'
      CHARACTER*4 SOBJ
      DOUBLE PRECISION THRESH
C local
      DOUBLE PRECISION RMS, RMST
      INTEGER NVIOL, NVIOLT, NRMS, NRMST, N
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
      LOGICAL ERR
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C loop through PIGs and check for unknown coordinates
      RMST=ZERO
      NRMST=0
      NVIOLT=0
      DO N=1,NPIG
      CALL ATMCHK(HEAP(IINTER(N)),ERR)
      IF (.NOT.ERR) THEN
      CALL PRPIC2(SOBJ,THRESH,RMS,NVIOL,NRMS,HEAP(IINTER(N)))
      NVIOLT=NVIOLT+NVIOL
      NRMST=NRMST+NRMS
      RMST=RMST+RMS
      END IF
      END DO
C
      IF (THRESH.GT.ZERO) THEN
      WRITE(6,'(A,F8.3,A,I5)')
     &     ' Number of violations greater ', THRESH, ': ', NVIOLT
      END IF
      DBPREC=NVIOLT
      CALL DECLAR( 'VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
C
      IF (NRMST.GT.0) THEN
      DBPREC = SQRT(RMST/NRMST)
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
      CALL DECLAR( 'RMS', 'DP', ' ', DBCOMP, DBPREC )
      WRITE(6,'(A,F8.3)') ' RMS deviation=',DBPREC
      ELSE
      DBPREC = ZERO
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
      CALL DECLAR( 'RMS', 'DP', ' ', DBCOMP, DBPREC )
      END IF
C
      RETURN
      END
C
      SUBROUTINE PRPIC2(SOBJ,THRESH,RMS,NVIOL,NRMS,INTERE)
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'update.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      CHARACTER*4 SOBJ
      DOUBLE PRECISION THRESH
      DOUBLE PRECISION RMS
      INTEGER NVIOL, INTERE(*)
C local
      DOUBLE PRECISION E, EDUMMY
      DOUBLE PRECISION CBBL, CBCL, CICL, CIBL, CPCL, CPBL, CTBL, CTCL
      INTEGER CIDL, CPDL
      INTEGER I, J, K, L, M, NRMS
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      NVIOL=0
      NRMS=0
      RMS=ZERO
      IF (SOBJ.EQ.'BOND') THEN
C update parameter lookup lists if required
      IF (UPBOND) THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,
     &            IAC,SEGID,RESID,TYPE)
      UPBOND=.FALSE.
      END IF
C
      WRITE(PUNIT,'(2A/)') ' (atom-i        |atom-j        )',
     &        '    dist.   equil.   delta    energy   const. '
      DO M=1,NBOND
      IF (ABS(INTERE(IB(M))+INTERE(JB(M))).LE.+1) THEN
      CBBL=ACBB(M)
      CBCL=ACBC(M)
      I=IB(M)
      J=JB(M)
      CALL EBOND2(E,EDUMMY,I,J,INTERE,1,CBCL,CBBL,'ANAL',ONE,ONE)
      IF (ABS(PCDATA(PCDEVI)) .GT. THRESH) THEN
      WRITE(PUNIT,'(13A,5(1X,F8.3))')
     &  ' (',SEGID(I),' ',RESID(I),' ',TYPE(I),'|',SEGID(J),' ',
     &  RESID(J),' ',TYPE(J),')',
     &  PCDATA(PCGEOM),PCDATA(PCEQUI),PCDATA(PCDEVI),PCDATA(PCENER),
     &  PCDATA(PCCONS)
      NVIOL=NVIOL+1
      END IF
      NRMS=NRMS+1
      RMS=RMS+PCDATA(PCDEVI)**2
      ENDIF
      END DO
C
      ELSE IF (SOBJ.EQ.'ANGL') THEN
C update parameter lookup lists if required
      IF (UPANGL) THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,
     &            ACTB,ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,IAC,
     &            SEGID,RESID,TYPE)
      UPANGL=.FALSE.
      END IF
      WRITE(PUNIT,'(2A/)') ' (atom-i        |atom-j        |atom-k ',
     &     '       )  angle    equil.     delta    energy  const. '
      DO M=1,NTHETA
      IF (ABS(INTERE(IT(M))+INTERE(JT(M))+INTERE(KT(M))).LE.+2) THEN
      CTBL=ACTB(M)
      CTCL=ACTC(M)
      I=IT(M)
      J=JT(M)
      K=KT(M)
      CALL EANGLE2(E,EDUMMY,I,J,K,INTERE,1,CTCL,CTBL,'ANAL',ONE,ONE)
      IF (ABS(PCDATA(PCDEVI)) .GT. THRESH) THEN
      WRITE(PUNIT,'(19A,5(1X,F8.3))')
     &  ' (',SEGID(I),' ',RESID(I),' ',TYPE(I),'|',SEGID(J),' ',
     &  RESID(J),' ',TYPE(J),'|',SEGID(K),' ',RESID(K),' ',TYPE(K),')',
     &  PCDATA(PCGEOM),PCDATA(PCEQUI),PCDATA(PCDEVI),
     &  PCDATA(PCENER),PCDATA(PCCONS)
C
C write the Urey-Bradley term if required
      CBBL=ACTUB(M)
      CBCL=ACTUC(M)
      IF (ABS(CBBL).GT.RSMALL.AND.ABS(CBCL).GT.RSMALL) THEN
      CALL EBOND2(E,EDUMMY,I,K,INTERE,1,CBCL,CBBL,'ANAL',ONE,ONE)
      WRITE(PUNIT,'(A,5(1X,F8.3))')
     &  '                        Urey-Bradley bond term:',
     &  PCDATA(PCGEOM),PCDATA(PCEQUI),PCDATA(PCDEVI),PCDATA(PCENER),
     &  PCDATA(PCCONS)
      END IF
C
      NVIOL=NVIOL+1
      END IF
      NRMS=NRMS+1
      RMS=RMS+PCDATA(PCDEVI)**2
      ENDIF
      END DO
C
      ELSE IF (SOBJ.EQ.'DIHE') THEN
C update parameter lookup lists if required
      IF (UPDIHE) THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      UPDIHE=.FALSE.
      END IF
      WRITE(PUNIT,'(3A/)') ' (atom-i        |atom-j        |',
     &   'atom-k        |atom-L        )',
     &   '    angle    equil.   delta    energy   const.   period'
      DO M=1,NPHI
      IF (ABS(INTERE(IP(M))+INTERE(JP(M))
     &    +INTERE(KP(M))+INTERE(LP(M))).LE.+3) THEN
      CPBL=ACPB(M)
      CPCL=ACPC(M)
      CPDL=ACPD(M)
      I=IP(M)
      J=JP(M)
      K=KP(M)
      L=LP(M)
      CALL ETOR(E,EDUMMY,I,J,K,L,INTERE,1,CPCL,CPDL,CPBL,'NORMAL',
     &          ZERO,0,'ANAL',ONE,ONE)
      IF (ABS(PCDATA(PCDEVI)) .GT. THRESH) THEN
      WRITE(PUNIT,'(25A,5(1X,F8.3),1X,I3)')
     &  ' (',SEGID(I),' ',RESID(I),' ',TYPE(I),'|',SEGID(J),' ',
     &  RESID(J),' ',TYPE(J),'|',SEGID(K),' ',RESID(K),' ',TYPE(K),'|',
     &  SEGID(L),' ',RESID(L),' ',TYPE(L),')',
     &  PCDATA(PCGEOM),PCDATA(PCEQUI),PCDATA(PCDEVI),
     &  PCDATA(PCENER),PCDATA(PCCONS),INT(PCDATA(PCPERI))
      NVIOL=NVIOL+1
      END IF
      NRMS=NRMS+1
      RMS=RMS+PCDATA(PCDEVI)**2
      ENDIF
      END DO
C
      ELSE IF (SOBJ.EQ.'IMPR') THEN
C update parameter lookup lists if required
      IF (UPIMPR) THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,QACIC,QACIB,QACID,
     &            IAC,SEGID,RESID,TYPE)
      UPIMPR=.FALSE.
      END IF
      WRITE(PUNIT,'(3A/)') ' (atom-i        |atom-j        |',
     &   'atom-k        |atom-L        )',
     &   '    angle    equil.   delta    energy   const.   period'
      DO M=1,NIMPHI
      IF (ABS(INTERE(IM(M))+INTERE(JM(M))
     &    +INTERE(KM(M))+INTERE(LM(M))).LE.+3) THEN
      CIBL=ACIB(M)
      CICL=ACIC(M)
      CIDL=ACID(M)
      I=IM(M)
      J=JM(M)
      K=KM(M)
      L=LM(M)
      CALL ETOR(E,EDUMMY,I,J,K,L,INTERE,1,CICL,CIDL,CIBL,'NORMAL',
     &          ZERO,0,'ANAL',ONE,ONE)
      IF (ABS(PCDATA(PCDEVI)) .GT. THRESH) THEN
      WRITE(PUNIT,'(25A,5(1X,F8.3),1X,I3)')
     &  ' (',SEGID(I),' ',RESID(I),' ',TYPE(I),'|',SEGID(J),' ',
     &  RESID(J),' ',TYPE(J),'|',SEGID(K),' ',RESID(K),' ',TYPE(K),'|',
     &  SEGID(L),' ',RESID(L),' ',TYPE(L),')',
     &  PCDATA(PCGEOM),PCDATA(PCEQUI),PCDATA(PCDEVI),
     &  PCDATA(PCENER),PCDATA(PCCONS),INT(PCDATA(PCPERI))
      NVIOL=NVIOL+1
      END IF
      NRMS=NRMS+1
      RMS=RMS+PCDATA(PCDEVI)**2
      ENDIF
      END DO
C
      END IF
C
      RETURN
      END
C
      SUBROUTINE PICK(QTRAJ)
C
C Routine allows to pick bonds, angles, dihedrals, impropers,...
C for arbitrary sets of atoms.  It can also produce trajectory
C series of the above properties.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'deriv.inc'
      LOGICAL QTRAJ
C local
      INTEGER ISLCT, JSLCT, KSLCT, LSLCT, MAXPUC, Q, P
      PARAMETER (MAXPUC=8)
C begin
      IF (NATOM+4.GT.MAXA) THEN
      CALL WRNDIE(-5,'PICK',
     & 'exceeded MAXA --> recompile program')
      ELSE
      CALL ATMINI(NATOM+1,NATOM+4)
      ISLCT=ALLHP(INTEG4(NATOM))
      JSLCT=ALLHP(INTEG4(NATOM))
      KSLCT=ALLHP(INTEG4(NATOM))
      LSLCT=ALLHP(INTEG4(NATOM))
      Q=ALLHP(IREAL8(MAXPUC))
      P=ALLHP(IREAL8(MAXPUC))
      CALL PICK2(QTRAJ,HEAP(ISLCT),HEAP(JSLCT),HEAP(KSLCT),
     &           HEAP(LSLCT),MAXPUC,HEAP(Q),HEAP(P))
      CALL FREHP(P,IREAL8(MAXPUC))
      CALL FREHP(Q,IREAL8(MAXPUC))
      CALL FREHP(LSLCT,INTEG4(NATOM))
      CALL FREHP(KSLCT,INTEG4(NATOM))
      CALL FREHP(JSLCT,INTEG4(NATOM))
      CALL FREHP(ISLCT,INTEG4(NATOM))
      END IF
      RETURN
      END
C
      SUBROUTINE PICK2(QTRAJ,ISLCT,JSLCT,KSLCT,LSLCT,MAXPUC,Q,P)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'pick.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'update.inc'
      LOGICAL QTRAJ
      INTEGER ISLCT(*), JSLCT(*), KSLCT(*), LSLCT(*), MAXPUC
      DOUBLE PRECISION Q(*), P(*)
C local
      INTEGER MXLIST
      PARAMETER (MXLIST=20)
      INTEGER NUNIT, BEGIN, SKIP, STOP, UNIT, ISTEP
      INTEGER FLIST(MXLIST), OUNIT, NP, NQ
      LOGICAL START, VDONE, QPARAM, QCMS, COND, QFORM, QMULT
      INTEGER NISLCT, NJSLCT, NKSLCT, NLSLCT, I, J, K, L, IC
      DOUBLE PRECISION DELTA, E, QT, EDUMMY
      CHARACTER*4 SPROP, SOBJ, HDR
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC, EMULT
      DOUBLE PRECISION CBBL, CBCL, CICL, CIBL, CPCL, CPBL, CTBL, CTCL
      DOUBLE PRECISION CTBUL, CTCUL
      INTEGER CIDL, CPDL, JC
      CHARACTER*32 AMP,PHI,ADST
      INTEGER STLEN,ADLEN
      LOGICAL CLOOP
C parameter
      INTEGER MARK
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, MARK=-9999)
C begin
C
C
C defaults
      NISLCT=1
      NJSLCT=1
      NKSLCT=1
      NLSLCT=1
      QPARAM=.FALSE.
      QCMS=.FALSE.
C parsing of PICK statement
      CALL NEXTA4('PICK-object=',SOBJ)
C
      IF (SOBJ.EQ.'HELP') THEN
C
      CALL CNSHELP('cns-pick')
C
      ELSE IF (SOBJ.EQ.'BOND') THEN
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL NEXTA4('property=',SPROP)
C update parameter lookup lists if required
      IF (UPBOND.AND.SPROP.NE.'GEOM') THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,
     &            IAC,SEGID,RESID,TYPE)
      UPBOND=.FALSE.
      END IF
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
      IF (NISLCT.GT.1.OR.NJSLCT.GT.1) THEN
      QCMS=.TRUE.
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' PICK: bond between CM of selected atoms'
      END IF
      I=NATOM+1
      J=NATOM+2
C use dummy parameters
      CBCL=ONE
      CBBL=ONE
      ELSE IF (NISLCT.EQ.1.AND.NJSLCT.EQ.1) THEN
      I=ISLCT(1)
      J=JSLCT(1)
      IC=FIND52(IB,JB,0,0,0,I,J,0,0,0,NBOND,2,MARK)
      IF (IC.EQ.MARK) THEN
      IC=FIND52(IB,JB,0,0,0,J,I,0,0,0,NBOND,2,MARK)
      END IF
      IF (IC.EQ.MARK) THEN
C use dummy parameters
      CBCL=ONE
      CBBL=ONE
      ELSE
      CBCL=ACBC(IC)
      CBBL=ACBB(IC)
      QPARAM=.TRUE.
      END IF
      END IF
C
      ELSE IF (SOBJ.EQ.'ANGL') THEN
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(KSLCT,NKSLCT,X,Y,Z,.TRUE.)
      CALL NEXTA4('property=',SPROP)
      IF (UPANGL.AND.SPROP.NE.'GEOM') THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,
     &            ACTB,ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,IAC,
     &            SEGID,RESID,TYPE)
      UPANGL=.FALSE.
      END IF
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
      CALL MAKIND(KSLCT,NATOM,NKSLCT)
      IF (NISLCT.GT.1.OR.NJSLCT.GT.1.OR.NKSLCT.GT.1) THEN
      QCMS=.TRUE.
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' PICK: angle between CM of selected atoms'
      END IF
      I=NATOM+1
      J=NATOM+2
      K=NATOM+3
C use dummy parameters
      CTCL=ONE
      CTBL=ONE
      CTCUL=ZERO
      CTBUL=ZERO
      ELSE IF (NISLCT.EQ.1.AND.NJSLCT.EQ.1.AND.NKSLCT.EQ.1) THEN
      I=ISLCT(1)
      J=JSLCT(1)
      K=KSLCT(1)
      IC=FIND52(IT,JT,KT,0,0,I,J,K,0,0,NTHETA,3,MARK)
      IF (IC.EQ.MARK) THEN
      IC=FIND52(IT,JT,KT,0,0,K,J,I,0,0,NTHETA,3,MARK)
      END IF
      IF (IC.EQ.MARK) THEN
C use dummy parameters
      CTCL=ONE
      CTBL=ONE
      CTCUL=ZERO
      CTBUL=ZERO
      ELSE
C use actual parameters
      CTCL=ACTC(IC)
      CTBL=ACTB(IC)
      CTCUL=ACTUC(IC)
      CTBUL=ACTUB(IC)
      QPARAM=.TRUE.
      END IF
      END IF
C
      ELSE IF (SOBJ.EQ.'DIHE') THEN
      QMULT=.FALSE.
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(KSLCT,NKSLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(LSLCT,NLSLCT,X,Y,Z,.TRUE.)
      CALL NEXTA4('property=',SPROP)
      IF (UPDIHE.AND.SPROP.NE.'GEOM') THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      UPDIHE=.FALSE.
      END IF
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
      CALL MAKIND(KSLCT,NATOM,NKSLCT)
      CALL MAKIND(LSLCT,NATOM,NLSLCT)
      IF (NISLCT.GT.1.OR.NJSLCT.GT.1.OR.NKSLCT.GT.1.OR.NLSLCT.GT.1) THEN
      QCMS=.TRUE.
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' PICK: dihedral between CM of selected atoms'
      END IF
      I=NATOM+1
      J=NATOM+2
      K=NATOM+3
      L=NATOM+4
C use dummy parameters
      CPCL=ONE
      CPDL=ONE
      CPBL=ONE
      ELSE IF (NISLCT.EQ.1.AND.NJSLCT.EQ.1.AND.NKSLCT.EQ.1
     &      .AND.NLSLCT.EQ.1) THEN
      I=ISLCT(1)
      J=JSLCT(1)
      K=KSLCT(1)
      L=LSLCT(1)
      IC=FIND52(IP,JP,KP,LP,0,I,J,K,L,0,NPHI,4,MARK)
      IF (IC.EQ.MARK) THEN
      IC=FIND52(IP,JP,KP,LP,0,L,K,J,I,0,NPHI,4,MARK)
      END IF
      IF (IC.EQ.MARK) THEN
C use dummy parameters
      CPCL=ONE
      CPDL=ONE
      CPBL=ONE
      ELSE
C use actual parameters
      CPCL=ACPC(IC)
      CPDL=ACPD(IC)
      CPBL=ACPB(IC)
      QPARAM=.TRUE.
C check if this is a multiple dihedral
      IF (IC.LT.NPHI) THEN
      IF (IP(IC+1).EQ.I.AND.JP(IC+1).EQ.J.AND.KP(IC+1).EQ.K.AND.
     &    LP(IC+1).EQ.L) QMULT=.TRUE.
      IF (IP(IC+1).EQ.L.AND.JP(IC+1).EQ.K.AND.KP(IC+1).EQ.J.AND.
     &    LP(IC+1).EQ.I) QMULT=.TRUE.
      END IF
      END IF
      END IF
C
      ELSE IF (SOBJ.EQ.'IMPR') THEN
      QMULT=.FALSE.
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(KSLCT,NKSLCT,X,Y,Z,.TRUE.)
      CALL SELCTA(LSLCT,NLSLCT,X,Y,Z,.TRUE.)
      CALL NEXTA4('property=',SPROP)
      IF (UPIMPR.AND.SPROP.NE.'GEOM') THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,QACIC,QACIB,QACID,
     &            IAC,SEGID,RESID,TYPE)
      UPIMPR=.FALSE.
      END IF
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
      CALL MAKIND(KSLCT,NATOM,NKSLCT)
      CALL MAKIND(LSLCT,NATOM,NLSLCT)
      IF(NISLCT.GT.1.OR.NJSLCT.GT.1.OR.NKSLCT.GT.1.OR.NLSLCT.GT.1) THEN
      QCMS=.TRUE.
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' PICK: improper between CM of selected atoms'
      END IF
      I=NATOM+1
      J=NATOM+2
      K=NATOM+3
      L=NATOM+4
C use dummy parameters
      CICL=ONE
      CIDL=ONE
      CIBL=ONE
      ELSE IF (NISLCT.EQ.1.AND.NJSLCT.EQ.1.AND.NKSLCT.EQ.1.AND.
     &         NLSLCT.EQ.1) THEN
      I=ISLCT(1)
      J=JSLCT(1)
      K=KSLCT(1)
      L=LSLCT(1)
      IC=FIND52(IM,JM,KM,LM,0,I,J,K,L,0,NIMPHI,4,MARK)
      IF (IC.EQ.MARK) THEN
      IC=FIND52(IM,JM,KM,LM,0,L,K,J,I,0,NIMPHI,4,MARK)
      END IF
      IF (IC.EQ.MARK) THEN
C use dummy parameters
      CICL=ONE
      CIDL=ONE
      CIBL=ONE
      ELSE
C use actual parameters
      CICL=ACIC(IC)
      CIDL=ACID(IC)
      CIBL=ACIB(IC)
      QPARAM=.TRUE.
C check if this is a multiple improper
      IF (IC.LT.NIMPHI) THEN
      IF (IM(IC+1).EQ.I.AND.JM(IC+1).EQ.J.AND.KM(IC+1).EQ.K.AND.
     &    LM(IC+1).EQ.L) QMULT=.TRUE.
      IF (IM(IC+1).EQ.L.AND.JM(IC+1).EQ.K.AND.KM(IC+1).EQ.J.AND.
     &    LM(IC+1).EQ.I) QMULT=.TRUE.
      END IF
      END IF
      END IF
C
      ELSE IF (SOBJ.EQ.'RING') THEN
      NISLCT=0
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      CALL SELCTA(JSLCT,NJSLCT,X,Y,Z,.TRUE.)
      CALL MAKIND(JSLCT,NATOM,NJSLCT)
      IF (NJSLCT.NE.1) THEN
      WRITE(6,'(A)') ' %PICK-ERR: selection has to refer to one atom'
      ELSE
      NISLCT=NISLCT+1
      ISLCT(NISLCT)=JSLCT(1)
      END IF
      CALL NEXTWD('RING>')
      CALL SAVEWD
      IF (WD(1:1).NE.'(') CLOOP=.FALSE.
      END DO
      CALL NEXTA4('property=',SPROP)
      IF (NISLCT.LE.3) THEN
      WRITE(6,'(A)') ' %PICK-ERR: "RING" requires at least 4 atoms'
      NISLCT=0
      END IF
      IF (NISLCT.GT.MAXPUC+3) THEN
      WRITE(6,'(A)')' %PICK-ERR: max. number of atoms in ring exceeded'
      NISLCT=0
      END IF
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I4,A)')
     &   ' PICK: a ',NISLCT,' membered ring has been selected'
      END IF
C
      END IF
C
C
C ============================================
C command parsing for trajectory series option
C ============================================
      IF (QTRAJ) THEN
      CALL DYNSET('INIT',USED,BEGIN,SKIP,STOP,NUNIT,MXLIST,FLIST,QFORM)
      OFILE='OUTPUT'
C
      CALL PUSEND('PICK>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PICK>')
      CALL DYNSET('PARSE',USED,BEGIN,SKIP,STOP,NUNIT,MXLIST,FLIST,
     &           QFORM)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPUT=',OFILE)
      ELSE
      CALL CHKEND('PICK>',DONE)
      END IF
      END IF
      END DO
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      DONE=.FALSE.
      END IF
C
      COND=NISLCT.GT.0.AND.NJSLCT.GT.0.AND.NKSLCT.GT.0.AND.NLSLCT.GT.0
      IF (COND) THEN
C
C =======================================================
C loop over trajectory or just take single coordinate set
C =======================================================
      START=.TRUE.
      CLOOP=.TRUE.
      DO WHILE (CLOOP)
      IF (QTRAJ) THEN
      CALL READDC(NATOM,X,Y,Z,START,VDONE,EOF,ERROR,
     &     NUNIT,FLIST,UNIT,BEGIN,STOP,SKIP,ISTEP,DELTA,HDR,QFORM)
      ELSE
      VDONE=.TRUE.
      END IF
C
C ================
C compute property
C ================
      IF (SOBJ.EQ.'BOND') THEN
      IF (QCMS) THEN
      CALL CMSIND(ISLCT,NISLCT,X(I),Y(I),Z(I))
      CALL CMSIND(JSLCT,NJSLCT,X(J),Y(J),Z(J))
      END IF
      CALL EBOND2(E,EDUMMY,I,J,ISLCT,1,CBCL,CBBL,'ANAL',ONE,ONE)
C
      ELSE IF (SOBJ.EQ.'ANGL') THEN
      IF (QCMS) THEN
      CALL CMSIND(ISLCT,NISLCT,X(I),Y(I),Z(I))
      CALL CMSIND(JSLCT,NJSLCT,X(J),Y(J),Z(J))
      CALL CMSIND(KSLCT,NKSLCT,X(K),Y(K),Z(K))
      END IF
C compute Urey-Bradley term if required
      IF (ABS(CTCUL).GT.RSMALL.AND.ABS(CTBUL).GT.RSMALL) THEN
      CALL EBOND2(E,EDUMMY,I,K,ISLCT,1,CTCUL,CTBUL,'ANAL',ONE,ONE)
      EMULT=PCDATA(PCENER)
      ELSE
      EMULT=ZERO
      END IF
      CALL EANGLE2(E,EDUMMY,I,J,K,ISLCT,1,CTCL,CTBL,'ANAL',ONE,ONE)
      PCDATA(PCENER)=PCDATA(PCENER)+EMULT
C
      ELSE IF (SOBJ.EQ.'DIHE') THEN
      IF (QCMS) THEN
      CALL CMSIND(ISLCT,NISLCT,X(I),Y(I),Z(I))
      CALL CMSIND(JSLCT,NJSLCT,X(J),Y(J),Z(J))
      CALL CMSIND(KSLCT,NKSLCT,X(K),Y(K),Z(K))
      CALL CMSIND(LSLCT,NLSLCT,X(L),Y(L),Z(L))
      END IF
      CALL ETOR(E,EDUMMY,I,J,K,L,ISLCT,1,CPCL,CPDL,CPBL,'NORMAL',
     &          ZERO,0,'ANAL',ONE,ONE)
C
C multiple dihedral: loop over all instances and accumulate energy
      IF (QMULT) THEN
      EMULT=PCDATA(PCENER)
      JC=IC
      DO WHILE (JC.LT.NPHI.AND.(IP(JC+1).EQ.I.AND.JP(JC+1).EQ.J.AND.
     &          KP(JC+1).EQ.K.AND.LP(JC+1).EQ.L).OR.
     &         (IP(JC+1).EQ.L.AND.JP(JC+1).EQ.K.AND.
     &          KP(JC+1).EQ.J.AND.LP(JC+1).EQ.I))
      JC=JC+1
      CPCL=ACPC(JC)
      CPDL=ACPD(JC)
      CPBL=ACPB(JC)
      CALL ETOR(E,EDUMMY,I,J,K,L,ISLCT,1,CPCL,CPDL,CPBL,'NORMAL',
     &          ZERO,0, 'ANAL',ONE,ONE)
      EMULT=EMULT+PCDATA(PCENER)
      ENDDO
      PCDATA(PCENER)=EMULT
      ENDIF
C
      ELSE IF (SOBJ.EQ.'IMPR') THEN
      IF (QCMS) THEN
      CALL CMSIND(ISLCT,NISLCT,X(I),Y(I),Z(I))
      CALL CMSIND(JSLCT,NJSLCT,X(J),Y(J),Z(J))
      CALL CMSIND(KSLCT,NKSLCT,X(K),Y(K),Z(K))
      CALL CMSIND(LSLCT,NLSLCT,X(L),Y(L),Z(L))
      END IF
      CALL ETOR(E,EDUMMY,I,J,K,L,ISLCT,1,CICL,CIDL,CIBL,'NORMAL',
     &           ZERO,0,'ANAL',ONE,ONE)
C
C multiple improper: loop over all instances and accumulate energy
      IF (QMULT) THEN
      EMULT=PCDATA(PCENER)
      JC=IC
      DO WHILE (JC.LT.NIMPHI.AND.(IM(JC+1).EQ.I.AND.JM(JC+1).EQ.J.AND.
     &          KM(JC+1).EQ.K.AND.LM(JC+1).EQ.L).OR.
     &         (IM(JC+1).EQ.L.AND.JM(JC+1).EQ.K.AND.
     &          KM(JC+1).EQ.J.AND.LM(JC+1).EQ.I))
      JC=JC+1
      CICL=ACIC(JC)
      CIDL=ACID(JC)
      CIBL=ACIB(JC)
      CALL ETOR(E,EDUMMY,I,J,K,L,ISLCT,1,CICL,CIDL,CIBL,'NORMAL',
     &          ZERO,0, 'ANAL',ONE,ONE)
      EMULT=EMULT+PCDATA(PCENER)
      ENDDO
      PCDATA(PCENER)=EMULT
      ENDIF
C
      ELSE IF (SOBJ.EQ.'RING') THEN
      CALL PUCKER(NISLCT,ISLCT,QT,NQ,Q,NP,P)
C
      END IF
C
C the output is a bit perverse... (have to change that in the future)
      DBPREC = ZERO
      IF (SOBJ.NE.'RING') THEN
      IF (SPROP.EQ.'GEOM') THEN
      DBPREC=PCDATA(PCGEOM)
      ELSE IF (SPROP.EQ.'ENER'.AND.QPARAM) THEN
      DBPREC=PCDATA(PCENER)
      ELSE IF (SPROP.EQ.'EQUI'.AND.QPARAM) THEN
      DBPREC=PCDATA(PCEQUI)
      ELSE IF (SPROP.EQ.'DELT'.AND.QPARAM) THEN
      DBPREC=PCDATA(PCDEVI)
      ELSE IF (SPROP.EQ.'CONS'.AND.QPARAM) THEN
      DBPREC=PCDATA(PCCONS)
      ELSE IF (SPROP.EQ.'PERI'.AND.QPARAM) THEN
      DBPREC=PCDATA(PCPERI)
      ELSE
      WRITE(6,'(A)')
     &  ' %PICK-ERR: undefined parameters for picked object'
      END IF
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
C
      IF (QTRAJ) THEN
      WRITE(OUNIT,'(1X,F14.6,1X,F14.6)') ISTEP*DELTA*TIMFAC,DBPREC
      ELSE
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A,F14.6)')  ' ',SPROP,'=',DBPREC
      END IF
      END IF
      ELSE
C
C pucker output
      DO I=1,NQ
      AMP='AMPLITUDE'
      STLEN=9
      CALL ENCODI(I,ADST,6,ADLEN)
      CALL ADDST(AMP,15,STLEN,ADST,ADLEN)
      DBPREC=Q(I)
      CALL DECLAR( AMP(1:STLEN), 'DP', ' ', DBCOMP, DBPREC )
      END DO
      DBPREC=NQ
      CALL DECLAR( 'N_AMPLITUDE','DP',' ',DBCOMP,DBPREC)
C
      DO I=1,NP
      PHI='PHASE'
      STLEN=5
      CALL ENCODI(I,ADST,6,ADLEN)
      CALL ADDST(PHI,15,STLEN,ADST,ADLEN)
      DBPREC=P(I)
      CALL DECLAR( PHI(1:STLEN), 'DP', ' ', DBCOMP, DBPREC )
      END DO
      DBPREC=NP
      CALL DECLAR( 'N_PHASE','DP',' ',DBCOMP,DBPREC)
C
      IF (QTRAJ) THEN
      WRITE(OUNIT,'(9(1X,F14.6))') ISTEP*DELTA*TIMFAC,
     &           (Q(I),I=1,NQ),(P(I),I=1,NP)
      ELSE
      WRITE(6,'(A,9(1X,F14.6))') ' PUCKer=',
     &           (Q(I),I=1,NQ),(P(I),I=1,NP)
      END IF
C
      END IF
C
      IF (VDONE .OR. EOF .OR. ERROR) CLOOP=.FALSE.
      END DO
C
      ELSE
      WRITE(6,'(A)') ' %PICK-ERR: zero selection specified'
      END IF
C
C
      RETURN
      END
C
      SUBROUTINE CMSIND(ISLCT,NISLCT,XC,YC,ZC)
C
C Routine computes center of mass of selected atoms and stores
C them into XC, YC, ZC.  "ISLCT" is assumed to be an index list.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INTEGER ISLCT(*), NISLCT
      DOUBLE PRECISION XC, YC, ZC
C local
      INTEGER II
      DOUBLE PRECISION TMASS
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      XC=ZERO
      YC=ZERO
      ZC=ZERO
      TMASS=ZERO
      DO II=1,NISLCT
      XC=XC+X(ISLCT(II))*AMASS(ISLCT(II))
      YC=YC+Y(ISLCT(II))*AMASS(ISLCT(II))
      ZC=ZC+Z(ISLCT(II))*AMASS(ISLCT(II))
      TMASS=TMASS+AMASS(ISLCT(II))
      END DO
      IF (TMASS.LT.RSMALL) THEN
      WRITE(6,'(A)') ' %CMSIND-ERR: total mass of selected atoms zero'
      ELSE
      XC=XC/TMASS
      YC=YC/TMASS
      ZC=ZC/TMASS
      END IF
      RETURN
      END
C
      SUBROUTINE PUCKER(N,INDEX,QT,NQ,Q,NPHI,PHI)
C
C Compute pucker coordinates (qt,q,phi) for the
C N-atom ring system for selected coordinates (INDEX).
C Dimensions: q(n/2-1), phi((n-1)/2-1), xyz(3,n)
C
C Algorithm from: Cremer&Pople, JACS 97:6, 1354-1358, 1975.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INTEGER N, INDEX(*)
      DOUBLE PRECISION QT
      INTEGER NQ
      DOUBLE PRECISION Q(*)
      INTEGER NPHI
      DOUBLE PRECISION PHI(*)
C local
      DOUBLE PRECISION RCX, RCY, RCZ, R1X, R1Y, R1Z, NORM, P
      DOUBLE PRECISION R2X, R2Y, R2Z, RJX, RJY, RJZ, RNX, RNY, RNZ
      DOUBLE PRECISION SUM, FACT, CSUM, SSUM, DOT, S, C
      INTEGER I, J, N2, NHALF, M
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, TWOPI, R90, R180, R270, R360
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (TWOPI=PI*TWO)
      PARAMETER (R90=90.0D0, R180=180.0D0, R360=360.0D0, R270=270.D0)
C begin
C
C compute geometrical center
      RCX=ZERO
      RCY=ZERO
      RCZ=ZERO
      DO I=1,N
      RCX=RCX+X(INDEX(I))
      RCY=RCY+Y(INDEX(I))
      RCZ=RCZ+Z(INDEX(I))
      END DO
      RCX=RCX/N
      RCY=RCY/N
      RCZ=RCZ/N
C
C normal vector
      R1X=ZERO
      R1Y=ZERO
      R1Z=ZERO
      R2X=ZERO
      R2Y=ZERO
      R2Z=ZERO
C
      DO J=1,N
      S=SIN(TWOPI*(J-1)/N)
      C=COS(TWOPI*(J-1)/N)
      R1X=R1X+(X(INDEX(J))-RCX)*S
      R1Y=R1Y+(Y(INDEX(J))-RCY)*S
      R1Z=R1Z+(Z(INDEX(J))-RCZ)*S
      R2X=R2X+(X(INDEX(J))-RCX)*C
      R2Y=R2Y+(Y(INDEX(J))-RCY)*C
      R2Z=R2Z+(Z(INDEX(J))-RCZ)*C
      END DO
C
C compute normalized cross-product
      RNX=R1Y*R2Z-R2Y*R1Z
      RNY=R1Z*R2X-R2Z*R1X
      RNZ=R1X*R2Y-R2X*R1Y
      NORM=SQRT(RNX**2+RNY**2+RNZ**2)
      RNX=RNX/NORM
      RNY=RNY/NORM
      RNZ=RNZ/NORM
C
C coordinate eqns
      N2=(N-1)/2
      DO M=2,N2
      CSUM=ZERO
      SSUM=ZERO
      DO J=1,N
      RJX=X(INDEX(J))-RCX
      RJY=Y(INDEX(J))-RCY
      RJZ=Z(INDEX(J))-RCZ
      DOT=(RJX*RNX+RJY*RNY+RJZ*RNZ)
      CSUM=CSUM+DOT*COS(TWOPI*M*(J-1)/N)
      SSUM=SSUM+DOT*SIN(TWOPI*M*(J-1)/N)
      END DO
C
      CSUM=CSUM*SQRT(TWO/N)
      SSUM=-SSUM*SQRT(TWO/N)
C
C Sums done, solve for Q & PHI
C Q*COS(PHI)=CSUM
C Q*SIN(PHI)=SSUM
C and get 0<PHI<360
      IF (ABS(CSUM).LT.RSMALL) THEN
      Q(M-1)=ABS(SSUM)
      IF (ABS(SSUM).LT.RSMALL) THEN
      PHI(M-1)=ZERO
      ELSE IF (SSUM.GT.0) THEN
      PHI(M-1)=R90
      ELSE
      PHI(M-1)=R270
      END IF
      ELSE IF (ABS(SSUM).LT.RSMALL) THEN
      Q(M-1)=ABS(CSUM)
      IF (CSUM.LT.0) THEN
      PHI(M-1)=R180
      ELSE
      PHI(M-1)=ZERO
      END IF
      ELSE
      P=ATAN2(SSUM,CSUM)
      IF (P.LT.0) P=TWOPI+P
      PHI(M-1)=P*R360/TWOPI
      Q(M-1)=CSUM/COS(P)
      END IF
      END DO
C
C One more q for even N
      IF (MOD(N,2).EQ.0) THEN
      SUM=ZERO
      FACT=ONE
      DO J=1,N
      RJX=X(INDEX(J))-RCX
      RJY=Y(INDEX(J))-RCY
      RJZ=Z(INDEX(J))-RCZ
      DOT=(RJX*RNX+RJY*RNY+RJZ*RNZ)
      SUM=SUM+FACT*DOT
      FACT=-FACT
      END DO
      NHALF=N/2
      Q(NHALF-1)=SUM*SQRT(ONE/N)
      N2=N2+1
      END IF
C
      N2=N2-1
      QT=ZERO
      DO I=1,N2
      QT=QT+Q(I)**2
      END DO
      QT=SQRT(QT)
C
      NQ=N2
      NPHI=(N-1)/2-1
      RETURN
      END
C
