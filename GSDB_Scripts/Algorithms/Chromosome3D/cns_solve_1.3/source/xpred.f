      SUBROUTINE XPRED(NA,NB,NC,NAP,NBP,
     &   NCP,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,XRTR,XRINTR,XRVOL,HPANOMFLAG,
     &   HPATOF,
     &   HPINDF,XRNATF,XRFP,XRFDP,XRSN,
     &   XRSM,XRSA,XRSB,XRSC,QLOOK,ELIM,
     &   XSFMX,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,XSFGNAM,XSFGTYP,XSFGORD,
     &   XRMREF,XRNREF,HPH,HPK,HPL,XRSYGP,XRSYIV,
     &   XRE,
     &   MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &   XRSCAL,XCVTEST,
     &   TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &   TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &   DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &   DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &   HPMULT,HPTYPE,
     &   QFFT,HPDX,HPDY,HPDZ,HPDT,HPDQ,IHPFLAG,XNAMEAS,XASSOC,
     &   QASELE)
C
C Computes electron density map from atomic model,
C computes structure factors from atomic model by direct summation,
C and the summation of atomic scattering factors.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'mtf.inc'
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      INTEGER ADEPTH, ARPNMX, ARPNN, ARPNX
      CHARACTER*(*) ARPN(4,*)
      INTEGER ARPNL(4,*)
      DOUBLE COMPLEX ARPNDB(4,*)
      INTEGER ARPNMLT(*)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
      LOGICAL XRMAP
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRVOL
      INTEGER HPANOMFLAG, HPATOF, HPINDF, XRNATF
      INTEGER XRSN, XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      LOGICAL QLOOK
      DOUBLE PRECISION ELIM
      INTEGER XSFMX, XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER XRMREF, XRNREF, HPH, HPK, HPL
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      DOUBLE PRECISION XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX), TRPNLEV(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
      INTEGER DRPNMX, DRPNX
      CHARACTER*(*) DRPN(4,DRPNMX,*)
      INTEGER DRPNL(4,DRPNMX,*), DRPNN(*)
      DOUBLE COMPLEX DRPNDB(4,DRPNMX,*)
      INTEGER DRPNMLT(DRPNMX,*), DRPNLEV(DRPNMX,*)
      CHARACTER*2 DRPNTYP(DRPNMX,*), DRPNDOM(DRPNMX,*)
      INTEGER DDEPTH(*)
      INTEGER HPMULT, HPTYPE
      LOGICAL QFFT
      INTEGER HPDX,HPDY,HPDZ,HPDT,HPDQ, IHPFLAG
      INTEGER XNAMEAS
      CHARACTER*(*) XASSOC(*)
      LOGICAL QASELE
C local
      LOGICAL COND, ERR, QDERIV
      INTEGER INSYM, I, LEN, NATSELE, IFCALC, II
      INTEGER ISTO
      INTEGER XRNATO, NANOM
      DOUBLE PRECISION SECS, SECS2
      CHARACTER*4 MODE, FORCE
      CHARACTER*(WDMAX) OBJECT
      CHARACTER*(WDMAX) SASSOC
      INTEGER FNDIOBJ
      EXTERNAL FNDIOBJ
      DOUBLE COMPLEX DBCOMP
C pointer
      INTEGER HPRMAP, HPIMAP, XL, YL, ZL, USE, AUSE, LTSEL
      INTEGER QSELE, ATSELE, PFCALC
      INTEGER HPATOM, HPINDX, HPFQS
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, SIX, R100
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, SIX=6.0D0)
      PARAMETER (R100=100.D0)
C begin
C
C allocate space for atom selection array
      ATSELE=ALLHP(INTEG4(NATOM))
C
C allocate space for reflection selection
      QSELE=ALLHP(ILOGIC(XRNREF))
C
C default values
      NATSELE=0
      QDERIV=.FALSE.
      ERR=.FALSE.
      MODE='REAL'
      OBJECT='MAP1'
      SASSOC='FCALC'
      LEN=4
C
      CALL PUSEND('PREDict>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PREDict>')
      CALL MISCOM('PREDict>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-predict')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODE=',MODE)
      IF (MODE.NE.'REAL'.AND.MODE.NE.'RECI'.AND.MODE.NE.'FF2 '
     &    .AND.MODE.NE.'DTAR')
     & THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: unknown mode'
      CALL WRNDIE(-5,'XPRED','mode undefined.')
      ERR=.TRUE.
      END IF
C
C parse dtarget object name
      IF (WD(1:4).EQ.'DTAR') THEN
      CALL NEXTWD('DTARget=')
      IF (WD(1:1).EQ.'(') THEN
      CALL NEXTWD('DTARget(reciprocal-space-object=')
      CALL COPYST(SASSOC,XNAMEAS,I,WD,WDLEN)
      CALL NEXTWD('DTARget(reciprocal-space-object=')
      IF (WD(1:1).NE.')') THEN
      CALL DSPERR('DTARget','")" expected.')
      ERR=.TRUE.
      END IF
C
C default is FCALC
      ELSE
      CALL SAVEWD
      SASSOC='FCALC'
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TO') THEN
      CALL NEXTWD('object=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('object=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') ' TO= ',OBJECT(1:LEN)
      ELSE
      CALL COPYST(OBJECT,WDMAX,LEN,WD,WDLEN)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,HEAP(QSELE),MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
      CALL SELCTA(HEAP(ATSELE),NATSELE,X,Y,Z,.TRUE.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     & 'precict-------------------------------------'
      WRITE(6,'(2A)') ' |   MODE = ',MODE
      WRITE(6,'(2A)') ' |   TO = ',OBJECT(1:LEN)
      WRITE(6,'(2A)') ' |--------------------------',
     & '-----------------------------------------'
      ELSE
      CALL CHKEND('XPRED>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C
      ERR=.FALSE.
C===
      IF (MODE.EQ.'REAL') THEN
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (OBJECT(1:LEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      ISTO=I
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: real space object ',
     & OBJECT(1:LEN),' not declared.'
      CALL WRNDIE(-5,'XPRED','object undeclared.')
      ERR=.TRUE.
      ELSEIF (NATSELE.EQ.0) THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: zero atoms selected. '
      CALL WRNDIE(-5,'XPRED','zero atoms selectected.')
      ERR=.TRUE.
      END IF
C===
      ELSE IF (MODE.EQ.'RECI'.OR.MODE.EQ.'FF2 '
     &     .OR.MODE.EQ.'DTAR') THEN
      COND=.FALSE.
      DO I=1,XSFNUM
      IF (OBJECT(1:LEN).EQ.XSFNAM(I)) THEN
      COND=.TRUE.
      ISTO=I
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: reciprocal space object ',
     & OBJECT(1:LEN),' not declared.'
      CALL WRNDIE(-5,'XPRED','object undeclared.')
      ERR=.TRUE.
      ELSEIF (NATSELE.EQ.0) THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: zero atoms selected. '
      CALL WRNDIE(-5,'XPRED','zero atoms selectected.')
      ERR=.TRUE.
      END IF
      IF (.NOT.ERR) THEN
      IF (MODE.EQ.'FF2 '.AND.XSFTYPE(ISTO).NE.'REAL') THEN
      WRITE(6,'(2A)') ' %XPRED-ERR: object ',
     & 'TYPE must be REAL for MODE=FF2.'
      CALL WRNDIE(-5,'XPRED','type mismatch.')
      ERR=.TRUE.
      END IF
      IF (MODE.EQ.'RECI'.AND.XSFTYPE(ISTO).NE.'COMP') THEN
      WRITE(6,'(2A)') ' %XPRED-ERR: object ',
     & 'TYPE must be COMPLEX for MODE=RECI.'
      CALL WRNDIE(-5,'XPRED','type mismatch.')
      ERR=.TRUE.
      END IF
      IF (MODE.EQ.'DTAR'.AND.XSFTYPE(ISTO).NE.'COMP') THEN
      WRITE(6,'(2A)') ' %XPRED-ERR: object ',
     & 'TYPE must be COMPLEX for MODE=DARGET.'
      CALL WRNDIE(-5,'XPRED','type mismatch.')
      ERR=.TRUE.
      END IF
      END IF
C===
      END IF
C
      IF  (.NOT.ERR) THEN
C
C get atom flag arrays, etc.
      IF (MODE.EQ.'RECI') THEN
      SASSOC=OBJECT
      END IF
      HPATOM=ALLHP(INTEG4(NATOM))
      HPINDX=ALLHP(INTEG4(NATOM))
      HPFQS=ALLHP(IREAL8(NATOM))
      CALL XRASSOC(HEAP(ATSELE),HEAP(HPANOMFLAG),
     &             HEAP(HPATOF),HEAP(HPINDF),
     &             HEAP(HPATOM),HEAP(HPINDX),HEAP(HPFQS),
     &             XRNATO,NANOM,
     &             XNAMEAS,SASSOC,QASELE,XRNATF,XRFDP,QHERM,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &             XRTR,XRINTR)
C
C ===============
C MODE=REAL
C ===============
C for generating electron density map from atomic model
      IF (MODE.EQ.'REAL') THEN
C
      IF (XRNATO.LE.0) THEN
      CALL WRNDIE(-5,'XPRED2',
     &       'Zero atoms selected or form factors missing.')
      ELSE
C
C allocate space for local atomic arrays
      XL=ALLHP(IREAL8(NATOM))
      YL=ALLHP(IREAL8(NATOM))
      ZL=ALLHP(IREAL8(NATOM))
C
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C
      IF (HPRRHO(ISTO).EQ.0) THEN
      CALL XMAPAL(HPRRHO(ISTO),HPIRHO(ISTO),QHERM,NRHO,IRHO)
      END IF
      HPRMAP=HPRRHO(ISTO)
      HPIMAP=HPIRHO(ISTO)
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C initialize map to zero
      CALL XRHOIN(MAASY,NAASY,MBASY,NBASY,MCASY,NCASY,
     &            QHERM,HEAP(HPRMAP),HEAP(HPIMAP))
C
C Loop over non-crystallographic symmetry operators to generate entire
C asymmetric unit of density from protomer coordinates.  The
C transformed coordinates are held in the arrays XL, YL, ZL.
      IF (XNNSYM.GT.1) THEN
      WRITE(6,'(A)') ' XPRED2: applying non-crystallographic symmetry.'
      END IF
      DO INSYM=1,XNNSYM
      CALL NCSXYZ(NATOM,X,Y,Z,
     &            INSYM,HEAP(XL),HEAP(YL),HEAP(ZL))
C
      CALL XPREDX(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   HEAP(HPRMAP),HEAP(HPIMAP),HEAP(HPRHOMA),
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,XRNATO,HEAP(XL),HEAP(YL),HEAP(ZL),
     &   XRSM,HEAP(HPINDX),XRSA,
     &   XRSB,XRSC,XRFP,HEAP(HPFQS),QLOOK,ELIM,HEAP(HPATOM))
C
      END DO
C
C free-up space
      CALL FREHP(ZL,IREAL8(NATOM))
      CALL FREHP(YL,IREAL8(NATOM))
      CALL FREHP(XL,IREAL8(NATOM))
C
      END IF
C
C ===============
C MODE=RECIprocal
C ===============
C Computes structure factors from selected coordinates
C FFT method or direct summation
      ELSE IF (MODE.EQ.'RECI') THEN
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      IF (XRNATO.LE.0) THEN
      CALL WRNDIE(-5,'XPRED2',
     &       'Zero atoms selected or form factors missing.')
      ELSE
C
C allocate space for TO reciprocal space object if necessary
      IF (HPSF(ISTO).EQ.0) THEN
      CALL XSFAL(HPSF(ISTO),XRMREF,XSFTYPE(ISTO))
      END IF
C
      LTSEL=ALLHP(INTEG4(XRNREF))
      CALL XRPRE8(.FALSE.,HEAP(HPSF(ISTO)),1,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,
     &   XRTR,XRINTR,XRVOL,HPATOF,
     &   HPINDF,XRNATF,XRFP,XRFDP,XRSN,
     &   XRSM,HPINDX,XRSA,XRSB,XRSC,HPFQS,QLOOK,ELIM,
     &   XRNATO,NANOM,HPATOM,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &   XRMREF,XRNREF,HPH,HPK,HPL,XRSYGP,XRSYIV,
     &   XRE,
     &   MBINS,XBINLOW,XBINHIGH,BINSHELL,MAPR,
     &   XRSCAL,XCVTEST,
     &   TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &   TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &   1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &   HPMULT,HPTYPE,
     &   QFFT,HPDX,HPDY,HPDZ,HPDT,HPDQ,HEAP(QSELE),
     &   LTSEL,HEAP(LTSEL))
      CALL FREHP(LTSEL,INTEG4(XRNREF))
C
      END IF
C ===============
C MODE=DTARGET
C ===============
C Computes structure factors from selected coordinates
C FFT method or direct summation
      ELSE IF (MODE.EQ.'DTAR') THEN
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      IF (XRNATO.LE.0) THEN
      CALL WRNDIE(-5,'XPRED2',
     &       'Zero atoms selected or form factors missing.')
      ELSE
C
C get FCALC index.  If non-existent, declare a structure factor object.
      IFCALC = FNDIOBJ(XSFNUM, XSFNAM,SASSOC)
      IF (IFCALC.LE.0) THEN
C declare object
      IF (XSFNUM.GE.XSFMX) THEN
      CALL WRNDIE(-5,'PREDict',
     & 'exceeded XSFMX parameter --> recompile program')
      ELSE
      XSFNUM=XSFNUM+1
      XSFNAM(XSFNUM)='FCALC'
      XSFTYPE(XSFNUM)='COMPLEX'
      XSFGNAM(XSFNUM)=0
      XSFGTYP(XSFNUM)=' '
      XSFGORD(XSFNUM)=0
      HPSF(XSFNUM)=0
      IFCALC=XSFNUM
      END IF
      END IF
      CALL CHKSFTYPE(IFCALC, XSFNAM, XSFTYPE, 'COMP', 'PREDict')
C
      IF (IFCALC.GT.0) THEN
C
C allocate space if necessary
      IF (HPSF(IFCALC).EQ.0) THEN
      CALL XSFAL(HPSF(IFCALC), XRMREF, XSFTYPE(IFCALC))
      END IF
C
C define PFCALC pointer
      PFCALC=HPSF(IFCALC)
C
C
C allocate space for TO reciprocal space object if necessary
      IF (HPSF(ISTO).EQ.0) THEN
      CALL XSFAL(HPSF(ISTO),XRMREF,XSFTYPE(ISTO))
      END IF
C
C look up which derivative (dtarget) definition corresponds
C to FCALC
C check if ASSOciate object is defined
      II=0
      DO I=1,IHPFLAG
      IF (XASSOC(I).EQ.SASSOC) THEN
      II=I
      END IF
      END DO
C
      IF (II.EQ.0) THEN
      WRITE(6,'(3A)') ' %XOFUNC-ERR: Object ',SASSOC,
     &                ' has no dtarget expression associated.'
      CALL WRNDIE(-5,'XRAY',
     & 'object not associated with dtarget expression.')
      ELSE
C
      LTSEL=ALLHP(INTEG4(XRNREF))
      CALL XRPRE8(.TRUE.,HEAP(PFCALC),HEAP(HPSF(ISTO)),
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,
     &   XRTR,XRINTR,XRVOL,HPATOF,
     &   HPINDF,XRNATF,XRFP,XRFDP,XRSN,
     &   XRSM,HPINDX,XRSA,XRSB,XRSC,HPFQS,QLOOK,ELIM,
     &   XRNATO,NANOM,HPATOM,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &   XRMREF,XRNREF,HPH,HPK,HPL,XRSYGP,XRSYIV,
     &   XRE,
     &   MBINS,XBINLOW,XBINHIGH,BINSHELL,MAPR,
     &   XRSCAL,XCVTEST,
     &   TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &   TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &   DRPNMX,DRPNX,DRPN(1,1,II),DRPNL(1,1,II),
     &   DRPNN(II),DRPNDB(1,1,II),
     &   DRPNMLT(1,II),DRPNLEV(1,II),DRPNTYP(1,II),
     &   DRPNDOM(1,II),DDEPTH(II),
     &   HPMULT,HPTYPE,
     &   QFFT,HPDX,HPDY,HPDZ,HPDT,HPDQ,HEAP(QSELE),
     &   LTSEL,HEAP(LTSEL))
      CALL FREHP(LTSEL,INTEG4(XRNREF))
C
      END IF
      END IF
      END IF
C
C ===============
C MODE=FF2
C ===============
C Computes the sum of the squared atomic form factors for
C selected atoms.
      ELSE IF (MODE.EQ.'FF2 ') THEN
C
      IF (XRNATO.LE.0) THEN
      CALL WRNDIE(-5,'XPRED2',
     &       'Zero atoms selected or form factors missing.')
      ELSE
      IF (HPSF(ISTO).EQ.0) THEN
      FORCE=XSFTYPE(ISTO)
      CALL XSFAL(HPSF(ISTO),XRMREF,XSFTYPE(ISTO))
      END IF
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      USE=ALLHP(IREAL8(XRSN))
      AUSE=ALLHP(IREAL8(XRSN))
      CALL XRFF2(HEAP(HPINDF),HEAP(HPATOF),XRNATF,
     &             XRSN,XRSM,XRFP,XRFDP,XRSA,XRSB,XRSC,
     &             XNNSYM,XRNSYM,XRNREF,
     &             HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPSF(ISTO)),
     &             NANOM,XRNATO,
     &             HEAP(HPINDX),HEAP(HPFQS),
     &             HEAP(USE),HEAP(AUSE),XRTR,HEAP(HPATOM))
      CALL FREHP(USE,IREAL8(XRSN))
      CALL FREHP(AUSE,IREAL8(XRSN))
      END IF
C
C=============
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' XPRED: CPU-time:  predict-calc=',SECS2-SECS
      END IF
C
      CALL FREHP(HPATOM,INTEG4(NATOM))
      CALL FREHP(HPINDX,INTEG4(NATOM))
      CALL FREHP(HPFQS,IREAL8(NATOM))
      END IF
C
C
      CALL FREHP(QSELE,ILOGIC(XRNREF))
C
C deallocate space for atom selection array
      CALL FREHP(ATSELE,INTEG4(NATOM))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XPREDX(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   RRHO,IRHO,RHOMASK,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,XRNATO,XL,YL,ZL,
     &   XRSM,XRINDX,XRSA,XRSB,XRSC,XRFP,XRFQS,QLOOK,ELIM,
     &   XRATOM)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'coord.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*), XRTR(3,3), XRINTR(3,3)
      LOGICAL QHERM
      INTEGER XRNATO
      DOUBLE PRECISION XL(*), YL(*), ZL(*)
      INTEGER XRSM, XRINDX(*)
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(*), XRFQS(*)
      LOGICAL QLOOK
      DOUBLE PRECISION ELIM
      INTEGER XRATOM(*)
C local
      INTEGER CUSHNA, CUSHNB, CUSHNC, CUSHMA, CUSHMB, CUSHMC
      INTEGER IAT, NN, I, J, SYM, AS, BS, CS, AA, BB, CC
      LOGICAL COND, ERROR
      DOUBLE PRECISION DX, DY, DZ
      INTEGER AX, BX, CX, NABC(3)
      DOUBLE PRECISION RCUT, RCUT2
      DOUBLE PRECISION ABC(3,3), ABCINV(3,3), FTAS, FTBS, FTCS
      DOUBLE PRECISION AEQ1, AEQ2, AEQ3, AEQ4, AEQ5
      DOUBLE PRECISION AE1, AE2, AE3, AE4, AE5
      DOUBLE PRECISION BE1, BE2, BE3, BE4, BE5
      DOUBLE PRECISION EXPLIM
C pointer
      INTEGER RHO, EXPTBL
C parameters
C  NTBL is the size of the EXP lookup tables for the FFT method
      INTEGER NTBL
      PARAMETER (NTBL=4000)
      DOUBLE PRECISION ZERO, FOUR
      PARAMETER (ZERO=0.0D0, FOUR=4.0D0)
C begin
C error check
      ERROR=.FALSE.
      DO IAT=1,XRNATO
      ERROR=ERROR.OR.WMAIN(XRATOM(IAT)).LT.RSMALL
      END DO
      IF (ERROR) THEN
      WRITE(6,'(A)')
     & ' %XPREDX-ERR: the B-factor is less or equal to zero for',
     & ' some atoms.  The electron density generation program ',
     & ' does not support this situation.'
      CALL WRNDIE(-5,'XPREDX','B <=0 for some atoms.')
      END IF
C
C
C compute sublattice transformation matrices
      NABC(1)=NA
      NABC(2)=NB
      NABC(3)=NC
      DO I=1,3
      DO J=1,3
      ABC(I,J)=NABC(I)*XRTR(I,J)
      END DO
      END DO
      DO I=1,3
      DO J=1,3
      ABCINV(I,J)=XRINTR(I,J)/NABC(J)
      END DO
      END DO
C
C compute norm of A*, B*, C* in FFT grid box
      FTAS=XRCELL(7)*NA
      FTBS=XRCELL(8)*NB
      FTCS=XRCELL(9)*NC
C
C initialize extent of temporary map (asu+cushion)
      CUSHNA=NAASY
      CUSHNB=NBASY
      CUSHNC=NCASY
      CUSHMA=MAASY
      CUSHMB=MBASY
      CUSHMC=MCASY
C
      DO IAT=1,XRNATO
C
C determine indices of lattice point closest to atomic position x,y,z
      DX= ABC(1,1)*(XL(XRATOM(IAT)))
     &   +ABC(1,2)*(YL(XRATOM(IAT)))
     &   +ABC(1,3)*(ZL(XRATOM(IAT)))
      AX=NINT(DX)
      DY= ABC(2,1)*(XL(XRATOM(IAT)))
     &   +ABC(2,2)*(YL(XRATOM(IAT)))
     &   +ABC(2,3)*(ZL(XRATOM(IAT)))
      BX=NINT(DY)
      DZ= ABC(3,1)*(XL(XRATOM(IAT)))
     &   +ABC(3,2)*(YL(XRATOM(IAT)))
     &   +ABC(3,3)*(ZL(XRATOM(IAT)))
      CX=NINT(DZ)
C
C map AX, BX, CX into asymmetric unit
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*AX+XRSYMM(SYM,1,2)*BX+XRSYMM(SYM,1,3)*CX
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*AX+XRSYMM(SYM,2,2)*BX+XRSYMM(SYM,2,3)*CX
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*AX+XRSYMM(SYM,3,2)*BX+XRSYMM(SYM,3,3)*CX
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
C
C map coordinates into asymmetric unit.  The result
C will overwrite the XL(IAT), YL(IAT), ZL(IAT) array!
      XL(XRATOM(IAT))=XRSYMM(SYM,1,1)*DX+XRSYMM(SYM,1,2)*DY+
     &         XRSYMM(SYM,1,3)*DZ
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH + AS-AA
      YL(XRATOM(IAT))=XRSYMM(SYM,2,1)*DX+XRSYMM(SYM,2,2)*DY+
     &         XRSYMM(SYM,2,3)*DZ
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH + BS-BB
      ZL(XRATOM(IAT))=XRSYMM(SYM,3,1)*DX+XRSYMM(SYM,3,2)*DY+
     &         XRSYMM(SYM,3,3)*DZ
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH + CS-CC
      DX=XL(XRATOM(IAT))
      DY=YL(XRATOM(IAT))
      DZ=ZL(XRATOM(IAT))
C
C copy index
      AX=NINT(DX)
      BX=NINT(DY)
      CX=NINT(DZ)
      ELSE
      CALL WRNDIE(-5,'XPREDX','Internal Error No. 3')
      END IF
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
C determine radius for atom IAT
      AEQ1=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),1)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT)))))**3
      AEQ2=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),2)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT)))))**3
      AEQ3=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),3)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT)))))**3
      AEQ4=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),4)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT)))))**3
      AEQ5=XRFQS(XRATOM(IAT))*
     &     (XRSC(XRINDX(IAT))+XRFP(XRINDX(IAT)))  *(SQRT(FOUR*PI/(
     &                               WMAIN(XRATOM(IAT)))))**3
      AE1=AEQ1*QMAIN(XRATOM(IAT))
      AE2=AEQ2*QMAIN(XRATOM(IAT))
      AE3=AEQ3*QMAIN(XRATOM(IAT))
      AE4=AEQ4*QMAIN(XRATOM(IAT))
      AE5=AEQ5*QMAIN(XRATOM(IAT))
      BE1=FOUR*PI**2/(XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT)))
      BE2=FOUR*PI**2/(XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT)))
      BE3=FOUR*PI**2/(XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT)))
      BE4=FOUR*PI**2/(XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT)))
      BE5=FOUR*PI**2/(WMAIN(XRATOM(IAT)))
C
C determine individual cutoff for electron density calculation
      RCUT2=MAX(ZERO,
     &           ( ELIM+LOG(MAX( ABS(AE1) , RSMALL )) )/BE1,
     &           ( ELIM+LOG(MAX( ABS(AE2) , RSMALL )) )/BE2,
     &           ( ELIM+LOG(MAX( ABS(AE3) , RSMALL )) )/BE3,
     &           ( ELIM+LOG(MAX( ABS(AE4) , RSMALL )) )/BE4,
     &           ( ELIM+LOG(MAX( ABS(AE5) , RSMALL )) )/BE5)
C
C compute maximum atomic radius
      RCUT=SQRT(RCUT2)
C
      CUSHNA=MAX(CUSHNA,INT(RCUT*FTAS-DX+RSMALL)+AX)
      CUSHNB=MAX(CUSHNB,INT(RCUT*FTBS-DY+RSMALL)+BX)
      CUSHNC=MAX(CUSHNC,INT(RCUT*FTCS-DZ+RSMALL)+CX)
      CUSHMA=MIN(CUSHMA,-INT(RCUT*FTAS+DX+RSMALL)+AX)
      CUSHMB=MIN(CUSHMB,-INT(RCUT*FTBS+DY+RSMALL)+BX)
      CUSHMC=MIN(CUSHMC,-INT(RCUT*FTCS+DZ+RSMALL)+CX)
      END DO
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(6(A,I6))')
     &  '  CUSHA=',CUSHMA,',...,',CUSHNA,
     &  '  CUSHB=',CUSHMB,',...,',CUSHNB,
     &  '  CUSHC=',CUSHMC,',...,',CUSHNC
      END IF
C
C allocate space for temporary map
      NN=(-CUSHMA+CUSHNA+1)*
     &   (-CUSHMB+CUSHNB+1)*
     &   (-CUSHMC+CUSHNC+1)
      RHO=ALLHP(IREAL8(NN))
C
C setup exponential lookup tables
      IF (QLOOK) THEN
      EXPTBL=ALLHP(IREAL8(NTBL+1))
      CALL XFLOOK(NTBL,HEAP(EXPTBL),EXPLIM)
      END IF
C
      CALL XPREDX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RHOMASK,HEAP(RHO),
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,ABC,
     &           ABCINV,XRNATO,FTAS,FTBS,FTCS,NTBL,EXPLIM,
     &           HEAP(EXPTBL),
     &           XRSM,XRINDX,XRSA,XRSB,XRSC,XRFP,
     &           XRFQS,QLOOK,ELIM,
     &           XRATOM,XL,YL,ZL)
C
      IF (QLOOK) THEN
      CALL FREHP(EXPTBL,IREAL8(NTBL+1))
      END IF
      CALL FREHP(RHO,IREAL8(NN))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XPREDX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RHOMASK,RHO,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,ABC,
     &           ABCINV,XRNATO,FTAS,FTBS,FTCS,NTBL,EXPLIM,EXPTBL,
     &           XRSM,XRINDX,XRSA,XRSB,XRSC,XRFP,XRFQS,QLOOK,ELIM,
     &           XRATOM,XL,YL,ZL)
C
C Computes electron density map of the atomic model
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER CUSHMC, CUSHNC, CUSHMB, CUSHNB, CUSHMA, CUSHNA
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION RHO(CUSHMA:CUSHNA,CUSHMB:CUSHNB,CUSHMC:CUSHNC)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION ABC(3,3)
      DOUBLE PRECISION ABCINV(3,3)
      INTEGER XRNATO
      DOUBLE PRECISION FTAS, FTBS, FTCS
      INTEGER NTBL
      DOUBLE PRECISION EXPLIM, EXPTBL(0:NTBL)
      INTEGER XRSM, XRINDX(*)
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFQS(*)
      LOGICAL QLOOK
      DOUBLE PRECISION ELIM
      INTEGER XRATOM(*)
      DOUBLE PRECISION XL(*), YL(*), ZL(*)
C local
      INTEGER A, B, C, AX, BX, CX, IAT, AA, BB, CC
      INTEGER AS, BS, CS
      DOUBLE PRECISION RCUT, RCUT2, DX, DY, DZ, DELX, DELY, DELZ, DEL2
      DOUBLE PRECISION DELXB, DELYB, DELZB, DELXC, DELYC, DELZC
      DOUBLE PRECISION AEQ1, AEQ2, AEQ3, AEQ4, AEQ5
      DOUBLE PRECISION AE1, AE2, AE3, AE4, AE5
      DOUBLE PRECISION BE1, BE2, BE3, BE4, BE5
      DOUBLE PRECISION BEE1, BEE2, BEE3, BEE4, BEE5
      INTEGER AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN, SYM
      LOGICAL COND
C parameters
      DOUBLE PRECISION ZERO, ONE, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, FOUR=4.0D0)
C begin
C
C
C initialize temporary density map
      DO C=CUSHMC,CUSHNC
      DO B=CUSHMB,CUSHNB
      DO A=CUSHMA,CUSHNA
      RHO(A,B,C)=ZERO
      END DO
      END DO
      END DO
C
C loop over all atoms.
      DO IAT=1,XRNATO
      AEQ1=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),1)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT)))))**3
      AEQ2=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),2)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT)))))**3
      AEQ3=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),3)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT)))))**3
      AEQ4=XRFQS(XRATOM(IAT))*XRSA(XRINDX(IAT),4)*(SQRT(FOUR*PI/(
     &           XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT)))))**3
      AEQ5=XRFQS(XRATOM(IAT))*
     &     (XRSC(XRINDX(IAT))+XRFP(XRINDX(IAT)))  *(SQRT(FOUR*PI/(
     &                               WMAIN(XRATOM(IAT)))))**3
      AE1=AEQ1*QMAIN(XRATOM(IAT))
      AE2=AEQ2*QMAIN(XRATOM(IAT))
      AE3=AEQ3*QMAIN(XRATOM(IAT))
      AE4=AEQ4*QMAIN(XRATOM(IAT))
      AE5=AEQ5*QMAIN(XRATOM(IAT))
      BE1=FOUR*PI**2/(XRSB(XRINDX(IAT),1)+WMAIN(XRATOM(IAT)))
      BE2=FOUR*PI**2/(XRSB(XRINDX(IAT),2)+WMAIN(XRATOM(IAT)))
      BE3=FOUR*PI**2/(XRSB(XRINDX(IAT),3)+WMAIN(XRATOM(IAT)))
      BE4=FOUR*PI**2/(XRSB(XRINDX(IAT),4)+WMAIN(XRATOM(IAT)))
      BE5=FOUR*PI**2/(WMAIN(XRATOM(IAT)))
C
C determine individual cutoff for electron density calculation
      RCUT2=MAX(ZERO,
     &           ( ELIM+LOG(MAX( ABS(AE1) , RSMALL )) )/BE1,
     &           ( ELIM+LOG(MAX( ABS(AE2) , RSMALL )) )/BE2,
     &           ( ELIM+LOG(MAX( ABS(AE3) , RSMALL )) )/BE3,
     &           ( ELIM+LOG(MAX( ABS(AE4) , RSMALL )) )/BE4,
     &           ( ELIM+LOG(MAX( ABS(AE5) , RSMALL )) )/BE5)
C
C compute maximum atomic radius
      RCUT=SQRT(RCUT2)
C
C precompute special coefficients for lookup-table mode
      IF (QLOOK) THEN
      BEE1=BE1*EXPLIM
      BEE2=BE2*EXPLIM
      BEE3=BE3*EXPLIM
      BEE4=BE4*EXPLIM
      BEE5=BE5*EXPLIM
      END IF
C
C determine indices of lattice point closest to atomic position x,y,z
      DX=XL(XRATOM(IAT))
      AX=NINT(DX)
      DY=YL(XRATOM(IAT))
      BX=NINT(DY)
      DZ=ZL(XRATOM(IAT))
      CX=NINT(DZ)
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
C convert this difference into real space
      DELX=ABCINV(1,3)*DZ+ABCINV(1,2)*DY+ABCINV(1,1)*DX
      DELY=ABCINV(2,3)*DZ+ABCINV(2,2)*DY+ABCINV(2,1)*DX
      DELZ=ABCINV(3,3)*DZ+ABCINV(3,2)*DY+ABCINV(3,1)*DX
C
C determine AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN of a
C box isomorphous to the unitcell geometry which surrounds
C the sphere with the atomic radius.  Compute maximum box.
      AMAXP=INT(RCUT*FTAS-DX+RSMALL)
      BMAXP=INT(RCUT*FTBS-DY+RSMALL)
      CMAXP=INT(RCUT*FTCS-DZ+RSMALL)
      AMAXN=INT(RCUT*FTAS+DX+RSMALL)
      BMAXN=INT(RCUT*FTBS+DY+RSMALL)
      CMAXN=INT(RCUT*FTCS+DZ+RSMALL)
C
      DO C=-CMAXN,CMAXP
      DELXC=DELX+ABCINV(1,3)*C
      DELYC=DELY+ABCINV(2,3)*C
      DELZC=DELZ+ABCINV(3,3)*C
C
      CC=C+CX
      DO B=-BMAXN,BMAXP
      DELXB=DELXC+ABCINV(1,2)*B
      DELYB=DELYC+ABCINV(2,2)*B
      DELZB=DELZC+ABCINV(3,2)*B
C
      BB=B+BX
      DO A=-AMAXN,AMAXP
      AA=A+AX
      DEL2=(ABCINV(1,1)*A+DELXB)**2
     &    +(ABCINV(2,1)*A+DELYB)**2
     &    +(ABCINV(3,1)*A+DELZB)**2
C
      IF (RCUT2.GT.DEL2) THEN
C
      IF (CC.LT.CUSHMC.OR.CC.GT.CUSHNC.OR.
     &    BB.LT.CUSHMB.OR.BB.GT.CUSHNB.OR.
     &    AA.LT.CUSHMA.OR.AA.GT.CUSHNA) THEN
      WRITE(6,'(A,I3)') ' aa,bb,cc=',AA,BB,CC
      WRITE(6,'(A,I3)') ' a,b,c=',A,B,C
      CALL WRNDIE(-5,'XPREDX','Internal Error No. 1')
      END IF
C
      IF (QLOOK) THEN
C use lookup tables
      RHO(AA,BB,CC)=RHO(AA,BB,CC)
     &    +AE1*EXPTBL(MIN(NTBL,INT(BEE1*DEL2)))
     &    +AE2*EXPTBL(MIN(NTBL,INT(BEE2*DEL2)))
     &    +AE3*EXPTBL(MIN(NTBL,INT(BEE3*DEL2)))
     &    +AE4*EXPTBL(MIN(NTBL,INT(BEE4*DEL2)))
     &    +AE5*EXPTBL(MIN(NTBL,INT(BEE5*DEL2)))
C
      ELSE
C don't use lookup tables
      RHO(AA,BB,CC)=RHO(AA,BB,CC)
     &    +AE1*EXP(-BE1*DEL2)
     &    +AE2*EXP(-BE2*DEL2)
     &    +AE3*EXP(-BE3*DEL2)
     &    +AE4*EXP(-BE4*DEL2)
     &    +AE5*EXP(-BE5*DEL2)
      END IF
C
      END IF
      END DO
      END DO
      END DO
      END DO
C
C map all grid points outside asymmetric unit into
C asymmetric unit
C =================================================
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
      IF (RHO(AA,BB,CC).GT.ZERO) THEN
C
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (COND) THEN
C
C just copy
      RRHO(AA,BB,CC)=RRHO(AA,BB,CC)+RHO(AA,BB,CC)
      ELSE
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      IF (COND) RRHO(AS,BS,CS)=RRHO(AS,BS,CS)+RHO(AA,BB,CC)
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'XPREDX','Internal error no. 2')
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XRHOIN(MAASY,NAASY,MBASY,NBASY,MCASY,NCASY,
     &                  QHERM,RRHO,IRHO)
C
C Routine initializes map to zero
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A, B, C
C parameter
      REAL ZERO
      PARAMETER (ZERO=0.0)
C begin
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      RRHO(A,B,C)=ZERO
      END DO
      END DO
      END DO
C
      IF (.NOT.QHERM) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IRHO(A,B,C)=ZERO
      END DO
      END DO
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRFF2(XRINDF,XRATOF,XRNATF,XRSN,XRSM,
     &                   XRFP,XRFDP,XRSA,XRSB,XRSC,XNNSYM,XRNSYM,
     &                   XRNREF,XRH,XRK,XRL,FFSQ,NANOM,XRNATO,
     &                   XRINDX,XRFQS,
     &                   USE,AUSE,XRTR,XRATOM)
C
C Computes the sum of squared atomic scattering factors
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INTEGER XRINDF(*), XRATOF(*), XRNATF
      INTEGER XRSN, XRSM
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      DOUBLE PRECISION XRSA(XRSM,*), XRSB(XRSM,*), XRSC(*)
      INTEGER XNNSYM, XRNSYM
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      DOUBLE PRECISION FFSQ(*)
      INTEGER NANOM, XRNATO
      INTEGER XRINDX(*)
      DOUBLE PRECISION XRFQS(*)
      DOUBLE PRECISION USE(*), AUSE(*), XRTR(3,3)
      INTEGER XRATOM(*)
C local
      INTEGER IAT , REFLCT, K
      DOUBLE PRECISION FFSUM, XRFFL
      DOUBLE PRECISION SQ
C parameter
      DOUBLE PRECISION ZERO, QUARTER
      PARAMETER (ZERO=0.0D0, QUARTER=0.25D0)
C begin
C
C make list of SCATTering factors actually in use
C and number of atoms for each scattering factor.
      DO K=1,XRSN
      USE(K)=ZERO
      AUSE(K)=ZERO
      END DO
C
      DO IAT=1,NANOM
      USE(XRINDX(IAT))=USE(XRINDX(IAT))
     &       +XRFQS(XRATOM(IAT))*QMAIN(XRATOM(IAT))
      END DO
C
      DO IAT=NANOM+1,XRNATO
      AUSE(XRINDX(IAT))=AUSE(XRINDX(IAT))
     &       +XRFQS(XRATOM(IAT))*QMAIN(XRATOM(IAT))
      END DO
C
C loop over all reflections
      DO REFLCT=1,XRNREF
      FFSUM=ZERO
C
C compute s**2/4 for this reflection
      SQ=  ((XRTR(1,1)*XRH(REFLCT)
     &     + XRTR(2,1)*XRK(REFLCT)
     &     + XRTR(3,1)*XRL(REFLCT))**2
     &     +(XRTR(1,2)*XRH(REFLCT)
     &     + XRTR(2,2)*XRK(REFLCT)
     &     + XRTR(3,2)*XRL(REFLCT))**2
     &     +(XRTR(1,3)*XRH(REFLCT)
     &     + XRTR(2,3)*XRK(REFLCT)
     &     + XRTR(3,3)*XRL(REFLCT))**2)*QUARTER
C
C compute atomic from factors
      DO IAT=1,XRSN
      IF (USE(IAT).GT.ZERO.OR.AUSE(IAT).GT.ZERO) THEN
C
      XRFFL=XRSC(IAT)+XRFP(IAT)
     &     +XRSA(IAT,1)*EXP(-XRSB(IAT,1)*SQ)
     &     +XRSA(IAT,2)*EXP(-XRSB(IAT,2)*SQ)
     &     +XRSA(IAT,3)*EXP(-XRSB(IAT,3)*SQ)
     &     +XRSA(IAT,4)*EXP(-XRSB(IAT,4)*SQ)
C
      IF (USE(IAT).GT.RSMALL) THEN
      FFSUM=FFSUM+XRFFL*XRFFL*USE(IAT)
      ELSEIF (AUSE(IAT).GT.RSMALL) THEN
      FFSUM=FFSUM+(XRFFL*XRFFL+XRFDP(IAT)*XRFDP(IAT))*AUSE(IAT)
      END IF
C
      END IF
      END DO
C
C multiply by number of NCS operators and number
C of crystallographic symmetry operators.
      FFSQ(REFLCT)=XNNSYM*XRNSYM*FFSUM
C
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE XRPRE8(QDERIV,FCALC,XDERIV,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,
     &   XRTR,XRINTR,XRVOL,HPATOF,
     &   HPINDF,XRNATF,XRFP,XRFDP,XRSN,
     &   XRSM,HPINDX,XRSA,XRSB,XRSC,HPFQS,QLOOK,ELIM,
     &   XRNATO,NANOM,HPATOM,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &   XRMREF,XRNREF,HPH,HPK,HPL,XRSYGP,XRSYIV,
     &   XRE,
     &   MBINS,XBINLOW,XBINHIGH,BINSHELL,MAPR,
     &   XRSCAL,XCVTEST,
     &   TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &   TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &   DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &   DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &   HPMULT,HPTYPE,
     &   QFFT,HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &   QSELE,LLTSEL,LTSEL)
C
C Computes fcalc and optionally dtarget/dfcalc
C from selected coordinates using fft method
C or direct summation method.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'mtf.inc'
      LOGICAL QDERIV
      DOUBLE COMPLEX FCALC(*), XDERIV(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRVOL
      INTEGER  HPATOF, HPINDF, XRNATF
      INTEGER XRSN, XRSM, HPINDX
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(*), XRFDP(*)
      INTEGER HPFQS
      LOGICAL QLOOK
      DOUBLE PRECISION ELIM
      INTEGER XRNATO, NANOM, HPATOM
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER XRMREF, XRNREF, HPH, HPK, HPL
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      DOUBLE PRECISION XRE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*), MAPR
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX), TRPNLEV(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
      INTEGER DRPNMX, DRPNX
      CHARACTER*(*) DRPN(4,DRPNMX)
      INTEGER DRPNL(4,DRPNMX), DRPNN
      DOUBLE COMPLEX DRPNDB(4,DRPNMX)
      INTEGER DRPNMLT(DRPNMX), DRPNLEV(DRPNMX)
      CHARACTER*2 DRPNTYP(DRPNMX), DRPNDOM(DRPNMX)
      INTEGER DDEPTH
      INTEGER HPMULT, HPTYPE
      LOGICAL QFFT
      INTEGER HPDX,HPDY,HPDZ,HPDT,HPDQ
      LOGICAL QSELE(*)
      INTEGER LLTSEL, LTSEL(*)
C
C local
      LOGICAL QDIRCT
      INTEGER REFLCT, COUNT
      DOUBLE COMPLEX DBCOMP
C pointers
      INTEGER XRNCSX, XRNCSY, XRNCSZ
      INTEGER XNCS, YNCS, ZNCS, DXNCS, DYNCS, DZNCS, EXPTBL
C parameters
C  NTBL is the size of the EXP lookup tables for the FFT method
      INTEGER NTBL
      PARAMETER (NTBL=4000)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      IF (XRNREF.GT.0.AND.XRNATO.GT.0) THEN
C
C decide whether direct summation method is required
      QDIRCT=.NOT.QFFT.OR.(XRNATO.GT.NANOM)
C
C allocate space for the FFT routines if required
      IF (QFFT) THEN
      XNCS=ALLHP(IREAL8(NATOM))
      YNCS=ALLHP(IREAL8(NATOM))
      ZNCS=ALLHP(IREAL8(NATOM))
      DXNCS=ALLHP(IREAL8(NATOM))
      DYNCS=ALLHP(IREAL8(NATOM))
      DZNCS=ALLHP(IREAL8(NATOM))
      EXPTBL=ALLHP(IREAL8(NTBL+1))
      END IF
C
C allocate space for the direct summation method if required
      IF (QDIRCT) THEN
      XRNCSX=ALLHP(IREAL8(NATOM))
      XRNCSY=ALLHP(IREAL8(NATOM))
      XRNCSZ=ALLHP(IREAL8(NATOM))
      END IF
C
C
C use reflection selection
      COUNT=0
      DO REFLCT=1,XRNREF
      IF (QSELE(REFLCT)) THEN
      COUNT=COUNT+1
      LTSEL(REFLCT)=1
      ELSE
      LTSEL(REFLCT)=0
      END IF
      END DO
C
      IF (COUNT.EQ.0) THEN
      WRITE(6,'(3A)') ' %XPRED-ERR: zero reflections selected.'
      CALL WRNDIE(-5,'XPRED','zero reflections selected.')
      ELSE
C
C initialize FCALC
      DO REFLCT=1,XRNREF
      FCALC(REFLCT)=DCMPLX(ZERO,ZERO)
      END DO
C
C compute FCALC
C
C Normal scatterers: call FFT method or direct summation
      IF (QFFT.AND.XRNATO.GT.0) THEN
      CALL XFFTUP(.FALSE.,.FALSE.,1,XRNATO,
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),FCALC,
     &     HEAP(HPATOM),HEAP(HPINDX),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XNCS),HEAP(YNCS),HEAP(ZNCS),
     &     HEAP(DXNCS),HEAP(DYNCS),HEAP(DZNCS),
     &     NTBL,HEAP(EXPTBL),QLOOK,XRSM,XRSA,XRSB,XRSC,XRFP,
     &     LTSEL,XRNREF,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &     QHERM,XRCELL,MAPR,XRTR,XRINTR,XRVOL)
C
C compute anomalous contribution: use direct summation
      IF (NANOM+1.LE.XRNATO) THEN
      CALL XRDRCT(LTSEL,.FALSE.,.FALSE.,NANOM+1,XRNATO,
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),FCALC,
     &     HEAP(HPATOM),HEAP(HPINDX),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.TRUE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      ELSEIF (.NOT.QFFT) THEN
      IF (NANOM.GT.0) THEN
      CALL XRDRCT(LTSEL,.FALSE.,.FALSE.,1,NANOM,
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),FCALC,
     &     HEAP(HPATOM),HEAP(HPINDX),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .FALSE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
C
C compute anomalous contribution: use direct summation
      IF (NANOM+1.LE.XRNATO) THEN
      CALL XRDRCT(LTSEL,.FALSE.,.FALSE.,NANOM+1,XRNATO,
     &     HEAP(HPH),HEAP(HPK),HEAP(HPL),FCALC,
     &     HEAP(HPATOM),HEAP(HPINDX),HEAP(HPFQS),
     &     HEAP(HPDX),HEAP(HPDY),HEAP(HPDZ),HEAP(HPDT),HEAP(HPDQ),
     &     HEAP(XRNCSX),HEAP(XRNCSY),HEAP(XRNCSZ),
     &     .TRUE.,.FALSE.,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRNREF)
      END IF
      END IF
C
C
C if QDERIV is true we have to compute the target function
C and its derivatives with respect to fcalc
      IF (QDERIV) THEN
      CALL XTARGETS(.TRUE.,.FALSE.,1,XRE,XDERIV,
     &           LLTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,1,
     &           XRCELL,XRVOL)
C
      END IF
C
C Free up heap space for FFT method
      IF (QFFT) THEN
      CALL FREHP(EXPTBL,IREAL8(NTBL+1))
      CALL FREHP(DZNCS,IREAL8(NATOM))
      CALL FREHP(DXNCS,IREAL8(NATOM))
      CALL FREHP(DYNCS,IREAL8(NATOM))
      CALL FREHP(ZNCS,IREAL8(NATOM))
      CALL FREHP(XNCS,IREAL8(NATOM))
      CALL FREHP(YNCS,IREAL8(NATOM))
      END IF
C
C free up heap space for direct summation method
      IF (QDIRCT) THEN
      CALL FREHP(XRNCSZ,IREAL8(NATOM))
      CALL FREHP(XRNCSY,IREAL8(NATOM))
      CALL FREHP(XRNCSX,IREAL8(NATOM))
      END IF
C
      END IF
      END IF
C
      RETURN
      END
