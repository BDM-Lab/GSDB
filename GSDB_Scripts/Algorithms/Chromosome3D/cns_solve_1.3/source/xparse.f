      SUBROUTINE XRAY
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
C begin
C
C check to see whether we have to allocate space for the
C XRAY atom list
      IF (XRUPAT) THEN
      XRUPAT=.FALSE.
      CALL XRAATM(NATOM)
C
C initialize the default associate selection
C ASSociate FCALC ( not hydrogen )
C
      IHPFLAG=1
      HPFLAG(1)=ALLHP(INTEG4(XRMATO))
      XASSOC(1)='FCALC'
      CALL XRINAT(HEAP(HPFLAG(1)),HEAP(HPANOMFLAG),NATOM)
      END IF
C
      CALL XREFI2
C
      RETURN
      END
C============================================================
      SUBROUTINE XREFI2
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'cnst.inc'
C local
      INTEGER NISLCT, NEWREF, IND
      INTEGER N, ITEMP, I, II
      INTEGER  SMAASY,SMBASY,SMCASY,SNAASY,SNBASY,SNCASY
      INTEGER  SHPRHOMA,SNRHO,SNMASK
      LOGICAL QANOM, OLHERM, ERR, QOK
      CHARACTER*4 SMETHO, SOBJ
      CHARACTER*(XNAMEAS) SASSOC
      DOUBLE PRECISION TEMP, TEMPS(2), ROLD, OMAPR
      DOUBLE PRECISION DBTEMP, XLOW, XHIGH
      DOUBLE COMPLEX DBCOMP
C pointer
      INTEGER XDERIV, ISLCT, PTATSELE
C parameters
      DOUBLE PRECISION ZERO, ONE, RAD, R001, R10000, M0001
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, RAD=PI/180.0D0)
      PARAMETER (R001=0.01D0, R10000=10000.D0, M0001=0.0001)
C begin
C
C
C allocate some default space
      IF (XRMREF.EQ.0) THEN
      CALL XRAREF(200,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
C
C all defaults are set in XRERES!
      CALL PUSEND('XRAY>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('XRAY>')
      CALL MISCOM('XRAY>',USED)
      XRQCHK=.TRUE.
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DECL') THEN
      CALL XDECLARE(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &   HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRNREF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RENA') THEN
      CALL XRENAME(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRMREF,XRNREF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UNDE') THEN
      CALL XUDECLA(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRMREF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'QUER') THEN
      CALL XQUERY(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   NRHO,IRHO,XRNREF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      CALL XDECGROUP(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UNGR') THEN
      CALL XDECUNGR(XNAMEMX,XSFMX,XRHOMX,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NREF') THEN
      NEWREF=XRMREF
      CALL NEXTI('NREFlections=',NEWREF)
      IF (NEWREF.GT.XRMREF) THEN
      WRITE(6,'(A,I7,A)')
     & ' XRAY: allocating space for ',NEWREF,' reflections.'
      CALL XRAREF(NEWREF,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'REFL') THEN
C
C check if a NREF command is specified.
      USED=.TRUE.
      DO WHILE (USED)
      CALL NEXTWD('REFLection>')
      CALL MISCOM('REFLection>',USED)
      END DO
      IF (WD(1:4).EQ.'NREF') THEN
      NEWREF=XRMREF
      CALL NEXTI('NREFlections=',NEWREF)
      IF (NEWREF.GT.XRMREF) THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I7,A)')
     & ' XRAY: increasing space allocation for up to ',
     &  NEWREF,' reflections.'
      END IF
C
      CALL XRAREF(NEWREF,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
      ELSE
      CALL SAVEWD
      END IF
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      CALL XRRRR(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &   HPSF,HPMULT,HPTYPE,XNAMEMX,XSFMX,
     &   QHERM,XRRED,XRREUP,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRINTR,
     &   XRSYGP,XRSYIV)
C
C set the target selection update flag
      XRREUP=.TRUE.
      XRRED=.TRUE.
      XQUICK=.FALSE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GENE') THEN
C
      DO I=1,2
      CALL NEXTWD('GENErate-res-limits=')
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTWD('GENErate-res-limits=')
      IF (WD(1:4).EQ.'INFI') THEN
      TEMPS(I)=-RSMALL
      ELSE
      CALL SAVEWD
      CALL NEXTF('GENErate-res-limits=',TEMPS(I))
      IF (TEMPS(I).LT.R001) THEN
      CALL DSPERR('XRAY',
     & 'resolution limit too small - set to 0.01. ')
      TEMPS(I)=R001
      END IF
      TEMPS(I)=ONE/TEMPS(I)
      END IF
      END DO
C
      IF (TEMPS(1).LT.TEMPS(2)) THEN
      XLOW=TEMPS(1)-RSMALL
      XHIGH=TEMPS(2)+RSMALL
      ELSE
      XLOW=TEMPS(2)-RSMALL
      XHIGH=TEMPS(1)+RSMALL
      END IF
      IF (XHIGH.EQ.ZERO) THEN
      CALL DSPERR('XRAY',
     & 'high resolution limit cant be infinity - set to 10000. ')
      XHIGH=R10000
      END IF
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C
      CALL XGENER(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRTR,XRINTR,XHIGH,XLOW)
C set the target selection update flag
      XRREUP=.TRUE.
      XRRED=.TRUE.
      XQUICK=.FALSE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DELE') THEN
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      ITEMP=XRNREF
      CALL XREFDEL(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,
     &           HPSF,HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
      IF (ITEMP.NE.XRNREF) THEN
      XRREUP=.TRUE.
      XRRED=.TRUE.
      XQUICK=.FALSE.
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'EXPA') THEN
C
C expand data
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C store current number of symmetry operators, number of reflections,
C and asymmetric parser info (no. of command lines).
C for quick reduction mode
      XRNREFRED=XRNREF
      XRNSYMRED=XRNSYM
      ARPNNRED=ARPNN
C
C only apply symmetry operators.
      CALL XEXPAN(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &      XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &      HPSF,HPMULT,HPTYPE,
     &      XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRTR,XRINTR)
C
C set number of symmetry operators to 1 (for P1)
      CALL XSYMRES(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &                   XRITSY,ARPNN)
C
C reduce reflections
      XRRED=.TRUE.
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C reset maps
      CALL XRMAPR(0)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' XRAY: data expanded and symmetry reset to P1.'
      END IF
C
C set the quick expansion/reduction info
      XQUICK=.TRUE.
      XQEXPAN=.TRUE.
      XRNREFEXP=XRNREF
C
C set the target selection update flag
      XRREUP=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UNEX') THEN
C
C quick unexpansion mode
      IF (XQUICK.AND.XQEXPAN) THEN
      XRNREF=XRNREFRED
      XRNSYM=XRNSYMRED
      XQEXPAN=.FALSE.
      DBTEMP=XRNSYM
      CALL DECLAR( 'SYMMETRY', 'DP', ' ', DBCOMP, DBTEMP )
C
C recover asymmetric parser info (no of command lines)
      ARPNN=ARPNNRED
C
C reset maps
      CALL XRMAPR(0)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' XRAY: performing quick unexpansion. '
      WRITE(6,'(A)') ' XRAY: data reduced and symmetry restored.'
      END IF
C
C reduce reflections
      XRRED=.TRUE.
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C set the target selection update flag
      XRREUP=.TRUE.
      ELSEIF (XQUICK.AND..NOT.XQEXPAN) THEN
      CALL WRNDIE(-5,'XRAY',
     & ' Diffraction data already reduced. ')
      ELSEIF (.NOT.XQUICK) THEN
      WRITE(6,'(A)')
     & ' %XRAY-err: Unexpansion impossible (due to intervening ',
     & '            SYMMmetry, ASYMmetry or REFL. statements).'
      CALL WRNDIE(-5,'XRAY',
     & ' Unexpansion impossible')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TARG') THEN
      CALL XEXPRDEF(TRPNMX,TRPNN,TRPNX,TRPN,TRPNL,TRPNDB,TRPNMLT,
     &             TRPNTYP,TRPNDOM,TRPNLEV,TDEPTH,XSFNUM,
     &             XSFNAM,XSFTYPE,'TARGet-expression>',QHERM,
     &             'EXPR')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DTAR') THEN
C
C check if associate object is specified
      ERR=.FALSE.
      CALL NEXTWD('DTARget=')
      IF (WD(1:1).EQ.'(') THEN
      CALL NEXTWD('DTARget(ASSOciate_object=')
      CALL COPYST(SASSOC,XNAMEAS,I,WD,WDLEN)
      CALL NEXTWD('DTARget(ASSOciate_object=')
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
C
      IF (.NOT.ERR) THEN
C
C check if ASSOciate object is defined
      II=0
      DO I=1,IHPFLAG
      IF (XASSOC(I).EQ.SASSOC) THEN
      II=I
      END IF
      END DO
C
      IF (II.EQ.0) THEN
      WRITE(6,'(3A)') ' %XRAY-ERR: Object ',SASSOC,
     &                ' is not ASSOciated.'
      CALL WRNDIE(-5,'XRAY',
     & 'object not associated with atom selection.')
      ELSE
      CALL XEXPRDEF(DRPNMX,DRPNN(II),DRPNX,DRPN(1,1,II),
     &             DRPNL(1,1,II),DRPNDB(1,1,II),DRPNMLT(1,II),
     &             DRPNTYP(1,II),DRPNDOM(1,II),DRPNLEV(1,II),DDEPTH(II),
     &             XSFNUM,
     &             XSFNAM,XSFTYPE,'DTARget-expression>',QHERM,
     &             'EXPR')
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MONI') THEN
      CALL XEXPRDEF(MRPNMX,MRPNN,MRPNX,MRPN,MRPNL,MRPNDB,MRPNMLT,
     &             MRPNTYP,MRPNDOM,MRPNLEV,MDEPTH,XSFNUM,
     &             XSFNAM,XSFTYPE,'MONItor-expression>',QHERM,
     &             'EXPR')
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TSEL') THEN
      CALL XEXPRDEF(SRPNMX,SRPNN,SRPNX,SRPN,SRPNL,SRPNDB,SRPNMLT,
     &             SRPNTYP,SRPNDOM,SRPNLEV,SDEPTH,XSFNUM,
     &             XSFNAM,XSFTYPE,'TSELection-expression>',QHERM,
     &             'SELE')
C set the target selection update flag
      XRREUP=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CVSE') THEN
      CALL XEXPRDEF(CRPNMX,CRPNN,CRPNX,CRPN,CRPNL,CRPNDB,CRPNMLT,
     &             CRPNTYP,CRPNDOM,CRPNLEV,CDEPTH,XSFNUM,
     &             XSFNAM,XSFTYPE,'CVSElection-expression>',QHERM,
     &             'SELE')
C set the target selection update flag
      XRREUP=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
      CALL XRQUER(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &            QHERM,XRCELL,XRTR)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OPTI') THEN
      CALL NEXTWD('optimize-qualifier=')
      IF (WD(1:4).EQ.'OVER') THEN
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C define target selection
      CALL XTSELSET(HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRREUP,HPTSEL,
     &           XRMREF,XRNREF,XCVTEST,
     &           XRTR,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C
      CALL XOVERB(HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTSEL))
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      CALL XGROUP
      ELSE IF (WD(1:4).EQ.'B-FA'.OR.WD(1:4).EQ.'BFAC') THEN
      CALL XROPTI
      ELSE
      CALL DSPERR('optimize-qualifier','unknown qualifier')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SCAL'.OR.WD(1:4).EQ.'MULT') THEN
C++++++++++++++++++++++++++++++++++++++++++
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
      CALL XSCALE(XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'LIDE') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
      PTATSELE = ALLHP(INTEG4(NATOM))
      CALL LIDENS(NA, NB, NC, NRHO,
     &            XRHONUM, XRHONAM, HPRRHO,
     &            XRNSYM, XRINTR, XRTR,
     &            NATOM, HEAP(PTATSELE), X, Y, Z)
      CALL FREHP(PTATSELE, INTEG4(NATOM))
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PSEA') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
      CALL PSEARCH(NA, NB, NC, NRHO,
     &             XRHONUM, XRHONAM, HPRRHO,
     &             XRNSYM, XRINTR, XRTR)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FMAP') THEN
      IF (.NOT. XQUICK) THEN
      CALL DSPERR('fmap','no previous expand')
      ELSE
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C get asymmetric unit for stored (inactive) symmetry
                 CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   SMAASY,SMBASY,SMCASY,SNAASY,SNBASY,SNCASY,
     &                   XRNSYMRED,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNNRED,
     &                   ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   SHPRHOMA,SNRHO,SNMASK,XRSYGP,XRSYIV)
      CALL FMAP(1, NA, NB, NC, NRHO, IRHO, QHERM,
     &          SMAASY, SMBASY, SMCASY,
     &          SNAASY, SNBASY, SNCASY,
     &          SHPRHOMA, SNMASK,
     &          XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &          XRNSYM, XRSYTH)
      IF (SHPRHOMA.NE.0) CALL FREHP(SHPRHOMA,INTEG4(SNRHO))
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MAPY') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
      CALL MAPYARD(1, NA, NB, NC, NRHO, IRHO, QHERM,
     &             XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &             XRNSYM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SMF ') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
      CALL SMF(NA, NB, NC, NRHO, IRHO, QHERM,
     &         XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &         XRNSYM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'IMF ') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
      PTATSELE = ALLHP(INTEG4(NATOM))
      CALL IMF(QHERM, NA, NB, NC, NRHO,
     &         XRHONUM, XRHONAM, HPRRHO,
     &         NATOM, XRNATF, HEAP(HPATOF), HEAP(HPINDF),
     &         X, Y, Z,
     &         XRSM, XRSA, XRSC, XRFP, XRFDP,
     &         XRNSYM, XRMSYM, XRSYTH, XRSYMM, XRITSY,
     &         XRTR, XRINTR,
     &         XNNSYM,
     &         HEAP(PTATSELE))
      CALL FREHP(PTATSELE, INTEG4(NATOM))
C=====================================================================
C #endif
C=====================================================================
      ELSE IF (WD(1:4).EQ.'HIST') THEN
      CALL XMHISTO(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,
     &  XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,
     &  IRHO,NMASK,XRMAP)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'AVER') THEN
      CALL XMAVER(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRTR,XRINTR,
     &     XRHONUM,XRHONAM,HPRRHO,HPIRHO,HPRHOMA,NRHO,IRHO,NMASK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SKEL') THEN
      CALL SKELETON(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,XRINTR,
     &  XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,IRHO,NMASK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FLIP') THEN
C flip all Friedel mates
      CALL XFLIP(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL))
C reduce reflections
      XRRED=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SEAR') THEN
      CALL NEXTWD('search-qualifier=')
      IF (WD(1:4).EQ.'ROTA') THEN
C
C real-space (map-based) rotation function
      CALL RSEARC
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (WD(1:4).EQ.'TSMA') THEN
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C define target selection
      CALL XTSELSET(HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRREUP,HPTSEL,
     &           XRMREF,XRNREF,XCVTEST,
     &           XRTR,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
      CALL TSMAP(HPH, HPK, HPL,
     &           XRNREF, HPTSEL, QHERM,
     &           NA, NB, NC, NRHO, IRHO,
     &           XRHONUM, XRHONAM, HPRRHO, HPIRHO,
     &           XSFNUM, XSFNAM, XSFTYPE, HPSF,
     &           XRNSYM, XRSYTH,
     &           XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           HPMULT,HPTYPE,XRMREF,XRCELL,XRVOL)
C
C rotation list clustering and sorting
      ELSE IF (WD(1:4).EQ.'CLUS') THEN
      CALL SCLUST
C=====================================================================
C #endif
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DIRE') THEN
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C reciprocal-space (direct) rotation function
      CALL DRSEAR
      ELSE
      CALL DSPERR('xray-search','unknown qualifier')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SPEC') THEN
      CALL XRSPECL(HPATOF,HPINDF,HPANOMFLAG,
     &     QASELE,
     &     XRNATF,XRFDP,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,XRTR,XRINTR)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'METH') THEN
      IF (QFFT) THEN
      SMETHO='FFT '
      ELSE
      SMETHO='DIRE'
      END IF
      CALL NEXTA4('METHod=',SMETHO)
      IF (SMETHO.EQ.'FFT ') THEN
      QFFT=.TRUE.
      ELSE IF (SMETHO.EQ.'DIRE') THEN
      QFFT=.FALSE.
      ELSE
      CALL DSPERR('XRAY>','only FFT or DIREct allowed')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'LOOK') THEN
      CALL NEXTLO('LOOKup=',QLOOK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FFT ') THEN
      CALL XFFT
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PROX') THEN
      CALL XPROX(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,XRHONUM,XRHONAM,HPRRHO,HPIRHO,HPRHOMA,NRHO,
     &  IRHO,NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &  HPH,HPK,HPL,
     &  XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &  HPMULT,HPTYPE,XRSYGP,XRSYIV)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MASK') THEN
      CALL XMASK(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,
     &  XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,
     &  IRHO,NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &  HPH,HPK,HPL,
     &  XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &  HPMULT,HPTYPE,XRSYGP,XRSYIV)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RESE') THEN
C
      WRITE(6,'(A)')
     & ' XRAY: whole xray database erased and reset'
      CALL XRAFRE
      CALL XRAATM(0)
      CALL XRMAPR(-1)
      CALL XRERES
      XRUPAT=.FALSE.
      CALL XRAATM(NATOM)
C
C initialize the default associate selection
C ASSociate FCALC ( not hydrogen )
C
      IHPFLAG=1
      HPFLAG(1)=ALLHP(INTEG4(XRMATO))
      XASSOC(1)='FCALC'
      CALL XRINAT(HEAP(HPFLAG(1)),HEAP(HPANOMFLAG),NATOM)
C
C
      CALL XRAREF(200,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'WA   ') THEN
      CALL NEXTF('WA=',XRSCAL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'A   ') THEN
      ROLD=XRCELL(1)
      CALL NEXTF('A=',XRCELL(1))
      IF (ROLD.NE.XRCELL(1)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C
      ELSE IF (WD(1:4).EQ.'B   ') THEN
      ROLD=XRCELL(2)
      CALL NEXTF('B=',XRCELL(2))
      IF (ROLD.NE.XRCELL(2)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C
      ELSE IF (WD(1:4).EQ.'C   ') THEN
      ROLD=XRCELL(3)
      CALL NEXTF('C=',XRCELL(3))
      IF (ROLD.NE.XRCELL(3)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C
      ELSE IF (WD(1:4).EQ.'ALPH') THEN
      ROLD=XRCELL(4)
      CALL NEXTF('ALPHa=',XRCELL(4))
      IF (ROLD.NE.XRCELL(4)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C
      ELSE IF (WD(1:4).EQ.'BETA') THEN
      ROLD=XRCELL(5)
      CALL NEXTF('BETA=',XRCELL(5))
      IF (ROLD.NE.XRCELL(5)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C
      ELSE IF (WD(1:4).EQ.'GAMM') THEN
      ROLD=XRCELL(6)
      CALL NEXTF('GAMMa=',XRCELL(6))
      IF (ROLD.NE.XRCELL(6)) THEN
      CALL XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      XRREUP=.TRUE.
      XRRED=.TRUE.
      CALL XRMAPR(0)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMM') THEN
C
      ITEMP=XRNSYM
      CALL XSYMPA(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &            XRITSY,ARPNN,XRSYGP,XRSYIV)
C
C
C need to update various flags if the symmetry operators have changed
      IF (ITEMP.NE.XRNSYM) THEN
      XRRED=.TRUE.
      XRREUP=.TRUE.
      XQUICK=.FALSE.
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ANOM'.OR.WD(1:4).EQ.'HERM') THEN
C
      OLHERM=QHERM
      QANOM=.NOT.QHERM
      IF (WD(1:4).EQ.'ANOM') THEN
      CALL NEXTLO('ANOMalous=',QANOM)
      QHERM=.NOT.QANOM
      ELSE
      CALL NEXTLO('HERMitian=',QHERM)
      QANOM=.NOT.QHERM
      END IF
C
C reset flags only if changed
      IF ((OLHERM.AND..NOT.QHERM).OR.(.NOT.OLHERM.AND.QHERM)) THEN
C
C===> flags have changed
C check if reflections are present
      IF (XRNREF.GT.0) THEN
C
      IF (QANOM) THEN
C
C===> means that we have turned-on the anomalous flag
C===> need to expand the existing data set
      WRITE(6,'(A)')
     & ' XPARSE: expanding data set (h,k,l)->{(h,k,l),(-h,-k,-l)}'
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,OLHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C expand data set using Friedel operator
      CALL XEXPFRIED(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYMM,XRTR,XRINTR)
C
      ELSE
C
C===> means that we have turned-off the anomalous flag
C===> need to average the Bijvoet mates and reduce the data
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(2A)')
     &  ' XPARSE: data will be reduced to hemisphere. ',
     &  '         Anomalous signal will be lost.'
      END IF
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,OLHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C average Bijvoet mates
      CALL XAVEFRIED(XRMREF,XRNREF,HPH,HPK,HPL,
     &           OLHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPTYPE)
C
      END IF
C
C reduce reflections
      XRRED=.TRUE.
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C reset map database
      CALL XRMAPR(0)
C set the target selection update flag
      XRREUP=.TRUE.
      XQUICK=.FALSE.
C
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FRIE') THEN
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
      CALL XAVEFRIED(XRMREF,XRNREF,HPH,HPK,HPL,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPTYPE)
C=====================================================================
      ELSE IF (WD(1:3).EQ.'DO') THEN
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      CALL XDO(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &         MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,
     &         ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &         XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &         HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &         XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &         XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &         HPMULT,
     &         HPTYPE,
     &         MBINS,XBINLOW,XBINHIGH,BINSHELL,XRSYGP,XRSYIV)
C set the target selection update flag
      XRREUP=.TRUE.
C
      ELSE IF (WD(1:4).EQ.'SHOW') THEN
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      CALL XSHOW(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRNSYM,XRMSYM,XRSYTH,
     &      XRSYMM,XRITSY,XRHONUM,XRHONAM,
     &      HPRRHO,HPIRHO,HPRHOMA,NRHO,NMASK,
     &      XRTR,XRNREF,
     &      HPH,HPK,HPL,
     &      XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &      HPMULT,HPTYPE,
     &      XRMREF,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C
      ELSE IF (WD(1:4).EQ.'ASYM') THEN
C
C disable quick expand mode
      XQUICK=.FALSE.
C
      CALL XASYMM(ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &            ARPNTYP,ARPNDOM,ARPNLEV,ADEPTH)
      CALL XRMAPR(0)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MAPR') THEN
      OMAPR=MAPR
      CALL NEXTF('MAPResolution=',MAPR)
      IF (MAPR.LT.R001) THEN
      CALL DSPERR('XRAY',
     & 'mapresolution limit too small - set to 0.01. ')
      MAPR=R001
      END IF
C
      IF (ABS(OMAPR-MAPR).GT.RSMALL) THEN
C set the map resolution
      CALL XRMAPR(0)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ASSO') THEN
      QOK=.FALSE.
      CALL NEXTWD('ASSOciate=')
      CALL COPYST(SASSOC,XNAMEAS,I,WD,WDLEN)
C
C reset: free up all associate atom selections
      IF (WD(1:4).EQ.'RESE') THEN
      DO I=1,IHPFLAG
      IF (HPFLAG(I).NE.0) CALL FREHP(HPFLAG(I),INTEG4(XRMATO))
      END DO
      IHPFLAG=0
      QOK=.TRUE.
      ELSEIF (WD(1:WDLEN).EQ.'?') THEN
      WRITE(6,'(A)')
     &  ' The following ASSOciate definitions are present:'
      WRITE(6,'(1X,8A)') (XASSOC(I),I=1,IHPFLAG)
      QOK=.TRUE.
      END IF
C
      IF (.NOT.QOK) THEN
C
C check if ASSOciate object is already defined:
      II=0
      DO I=1,IHPFLAG
      IF (XASSOC(I).EQ.SASSOC) THEN
      II=I
      END IF
      END DO
C
      IF (II.GT.0) THEN
C
C overwrite object
      CALL SELCTA(HEAP(HPFLAG(II)),NISLCT,X,Y,Z,.TRUE.)
C
      ELSE
      IF (IHPFLAG.LT.MAXHPFLAG) THEN
C
C create new object
      IHPFLAG=IHPFLAG+1
      HPFLAG(IHPFLAG)=ALLHP(INTEG4(XRMATO))
      CALL SELCTA(HEAP(HPFLAG(IHPFLAG)),NISLCT,X,Y,Z,.TRUE.)
      XASSOC(IHPFLAG)=SASSOC
      ELSE
      WRITE(6,'(A)') ' %XRAY-ERR: Too many ASSOciate statements.'
      CALL WRNDIE(-5,'XRAY',
     & 'exceeded MAXHPFLAG parameter --> recompile program')
      END IF
      END IF
      END IF
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ASEL') THEN
      CALL SELCTA(HEAP(HPANOMFLAG),NISLCT,X,Y,Z,.TRUE.)
      QASELE=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SCAT') THEN
      ISLCT=ALLHP(INTEG4(NATOM))
      CALL XRSCTP(HEAP(ISLCT),HEAP(HPINDF),HEAP(HPATOF),XRSN,
     &                  XRNATF,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
      CALL FREHP(ISLCT,INTEG4(NATOM))
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PRED') THEN
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      CALL XPRED(NA,NB,NC,NAP,NBP,
     &   NCP,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,XRTR,XRINTR,XRVOL,HPANOMFLAG,HPATOF,
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
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',XRLTOL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MBIN'.OR.WD(1:WDLEN).EQ.'BIN'.OR.
     &         WD(1:WDLEN).EQ.'BINS') THEN
      CALL NEXTI('BINs=',MBINS)
      IF (MBINS.LT.1) THEN
      MBINS=1
      CALL WRNDIE(-1,'XRAY',
     & ' number of bins must be larger than 0.  Reset to 1. ')
      ELSEIF (MBINS.GE.MAXMBIN) THEN
      CALL WRNDIE(-1,'XRAY',
     & ' number of bins exceeds maximum.  Reset to maximum.  ')
      MBINS=MAXMBIN-1
      END IF
      TEMP=MBINS
      CALL DECLAR( 'BIN_NUMBER', 'DP', ' ', DBCOMP, TEMP)
      IF (XBINLOW.NE.ZERO.AND.XBINHIGH.NE.ZERO) THEN
      CALL XMAKEBIN(MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE)
      END IF
C=====================================================================
      ELSE IF (WD(1:6).EQ.'BINRES') THEN
      CALL NEXTWD('BINRESolution-limits=')
      IF (WD(1:2).EQ.'? ') THEN
      IF (XBINLOW.EQ.ZERO.AND.XBINHIGH.EQ.0) THEN
      WRITE(6,'(A)') ' BINResolution is not set.'
      ELSEIF (XBINLOW.GT.ZERO) THEN
      WRITE(6,'(A,F7.2,A,F7.2)')
     & ' BINRESolution: ',ONE/XBINLOW,' to ',ONE/XBINHIGH
      WRITE(6,'(12F6.2)') (ONE/SQRT(BINSHELL(IND)),IND=MBINS+1,1,-1)
      ELSE
      WRITE(6,'(A,F7.2)')
     & ' BINRESolution range: INFInity to ',ONE/XBINHIGH
      WRITE(6,'(12F6.2)') (ONE/SQRT(BINSHELL(IND)),IND=MBINS,1,-1)
      END IF
      ELSE
C
      QBINSET=.TRUE.
      CALL SAVEWD
      DO I=1,2
      CALL NEXTWD('BINRESolution-limits=')
      IF (WD(1:WDLEN).EQ.'=') CALL NEXTWD('BINRESolution-limits=')
      IF (WD(1:4).EQ.'INFI') THEN
      TEMPS(I)=-RSMALL
      ELSE
      CALL SAVEWD
      CALL NEXTF('RESOlution-limits=',TEMPS(I))
      IF (TEMPS(I).LT.R001) THEN
      CALL DSPERR('XRAY',
     & 'resolution limit too small - set to 0.01. ')
      TEMPS(I)=R001
      END IF
      TEMPS(I)=ONE/TEMPS(I)
      END IF
      END DO
C
      IF (TEMPS(1).LT.TEMPS(2)) THEN
      XBINLOW=TEMPS(1)-RSMALL
      XBINHIGH=TEMPS(2)+RSMALL
      ELSE
      XBINLOW=TEMPS(2)-RSMALL
      XBINHIGH=TEMPS(1)+RSMALL
      END IF
      IF (XBINHIGH.EQ.ZERO) THEN
      CALL DSPERR('XRAY',
     & 'high resolution limit cant be infinity - set to 10000. ')
      XBINHIGH=R10000
      END IF
C
      IF (XBINLOW.GT.ZERO) THEN
      TEMP=ONE/XBINLOW
      CALL DECLAR('BIN_RESOLUTION_LOW', 'DP', ' ', DBCOMP, TEMP)
      ELSE
      CALL DECLAR('BIN_RESOLUTION_LOW', 'ST', 'INFINITY', DBCOMP, TEMP)
      END IF
      TEMP=ONE/XBINHIGH
      CALL DECLAR('BIN_RESOLUTION_HIGH', 'DP', ' ', DBCOMP, TEMP)
C
C generate array that contains the bin spacings
      CALL XMAKEBIN(MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE)
      END IF
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'STAT') THEN
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
      CALL XSTATS2(XRRED,XRREUP,
     &      XRSCAL,
     &      NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRNSYM,XRMSYM,XRSYTH,
     &      XRSYMM,XRITSY,XRTR,XRINTR,XRNREF,HPH,HPK,HPL,
     &      XSFNUM,XSFNAM,XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &      MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE,
     &      XRSYGP,XRSYIV,XRCELL,XRVOL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PRIN') THEN
      CALL NEXTA4('PRINt-object=',SOBJ)
C===
      IF (SOBJ.EQ.'TARG') THEN
C
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C define target selection
      CALL XTSELSET(HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRREUP,HPTSEL,
     &           XRMREF,XRNREF,XCVTEST,
     &           XRTR,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C use the reciprocal space targets
C
C for cross-validation testing print the target value for the test set
      IF (XCVTEST) THEN
C
C test set, no derivatives.  Using dummy parameters for dtarget expression.
      CALL XTARGETS(.FALSE.,.TRUE.,-1,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,0,
     &           XRCELL,XRVOL)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,F15.5)') ' TARGET value (test set)=',XRE
      END IF
      CALL DECLAR( 'TEST_TARGET', 'DP', ' ', DBCOMP, XRE )
      END IF
C
C working set, no derivatives.  Using dummy parameters for dtarget expression.
      CALL XTARGETS(.FALSE.,.TRUE.,1,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,0,
     &           XRCELL,XRVOL)
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,F15.5)') ' TARGET value (working set)=',XRE
      END IF
      CALL DECLAR( 'TARGET', 'DP', ' ', DBCOMP, XRE )
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, XRE )
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PEAK') THEN
      CALL XPEAKPIK(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,XRINTR,XRTR,
     &  XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,IRHO,NMASK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'READ') THEN
      CALL NEXTA4('read-object=',SOBJ)
      IF (SOBJ.EQ.'MAP') THEN
      CALL XMREAD(
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &           ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,XRSYGP,XRSYIV)
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSE IF (SOBJ.EQ.'MASK') THEN
      CALL XMSKRD(
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &           ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,XRSYGP,XRSYIV)
C=====================================================================
C #endif
C=====================================================================
      ELSE
      CALL DSPERR('read-object','unknown object')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'WRIT') THEN
      CALL NEXTA4('write-object=',SOBJ)
      IF (SOBJ.EQ.'REFL') THEN
C++++++++++++++++++++++++++++++++++++++++++
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
      CALL XWRIT(XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XSFMX,XRCELL,XRVOL)
C++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (SOBJ.EQ.'MAP') THEN
      CALL XMAPX('WRITE',
     &           NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &           XRTR,XRINTR,XRVOL,
     &           XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRRED,XRSYGP,XRSYIV)
C++++++++++++++++++++++++++++++++++++++++++
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSEIF (SOBJ.EQ.'MASK') THEN
      CALL XMSKWR(
     &            NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &            MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &            XRSYMM,XRITSY,QHERM,XRCELL,MAPR,
     &            XRHONUM,XRHONAM,HPRRHO,
     &            HPRHOMA,NRHO,NMASK,XRMAP,
     &            XRTR,XRINTR,XRVOL,
     &            XRMREF,XRNREF,HPH,HPK,HPL,
     &            XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &            HPMULT,
     &            HPTYPE,XRRED,XRSYGP,XRSYIV)
C=====================================================================
C #endif
C=====================================================================
      ELSE
      CALL DSPERR('write-object','unknown object')
      END IF
C=====================================================================
      ELSE
      CALL CHKEND('XRAY>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRSCTP(ISLCT,XRINDF,XRATOF,XRSN,
     &                  XRNATF,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP)
C
C Parsing of SCATter statement
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INTEGER ISLCT(*), XRINDF(*), XRATOF(*)
      INTEGER XRSN, XRNATF, XRSM
      DOUBLE PRECISION XRSA(XRSM,4), XRSB(XRSM,4), XRSC(XRSM)
      DOUBLE PRECISION XRFP(XRSM), XRFDP(XRSM)
C local
      INTEGER IAT, NISLCT, I, II
      LOGICAL QFIRST
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      CALL NEXTWD('SCATter=')
      IF (WD(1:4).EQ.'RESE') THEN
C
C reset atomic form factor data structure
      XRSN=0
      XRNATF=0
      ELSE IF (WD(1:1).EQ.'?') THEN
C
C print atomic form factors
      IF (XRSN.GT.0) WRITE(6,'(A)') ' Atomic form factors'
      DO I=1,XRSN
C
      NISLCT=0
      DO IAT=1,XRNATF
      IF (XRINDF(IAT).EQ.I) THEN
      NISLCT=NISLCT+1
      END IF
      END DO
C
      IF (NISLCT.GT.0) THEN
      II=0
      QFIRST=.TRUE.
      DO IAT=1,XRNATF
      IF (XRINDF(IAT).EQ.I) THEN
      II=II+1
      IF (QFIRST.AND.NISLCT.GT.1) THEN
      QFIRST=.FALSE.
      WRITE(PUNIT,'(7A)') ' SCATter ( (segid "',SEGID(XRATOF(IAT)),
     & '" and resid "',RESID(XRATOF(IAT)),
     & '" and name "',TYPE(XRATOF(IAT)),'" )'
      ELSEIF (QFIRST.AND.NISLCT.EQ.1) THEN
      WRITE(PUNIT,'(7A)') ' SCATter (segid "',SEGID(XRATOF(IAT)),
     & '" and resid "',RESID(XRATOF(IAT)),
     &  '" and name "',TYPE(XRATOF(IAT)),'" )'
      ELSEIF (.NOT.QFIRST.AND.II.LT.NISLCT) THEN
      WRITE(PUNIT,'(7A)') '        OR (segid "',SEGID(XRATOF(IAT)),
     & '" and resid "',RESID(XRATOF(IAT)),
     &  '" and name "',TYPE(XRATOF(IAT)),'" )'
      ELSEIF (.NOT.QFIRST.AND.II.EQ.NISLCT) THEN
      WRITE(PUNIT,'(7A)') '        OR (segid "',SEGID(XRATOF(IAT)),
     & '" and resid "',RESID(XRATOF(IAT)),
     &  '" and name "',TYPE(XRATOF(IAT)),'" ))'
      END IF
      END IF
      END DO
C
      WRITE(PUNIT,'(11(A,F7.3))')
     &  ' ',XRSA(I,1),' ',XRSB(I,1),
     &  ' ',XRSA(I,2),' ',XRSB(I,2),' ',XRSA(I,3),
     &  ' ',XRSB(I,3),' ',XRSA(I,4),' ',XRSB(I,4),
     &  ' ',XRSC(I),' FP=',XRFP(I),'  FDP=',XRFDP(I)
      END IF
      END DO
      ELSE
C
C do SCATter statement parsing
      CALL SAVEWD
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      IF (XRSN.GE.XRSM) THEN
      CALL WRNDIE(-5,'XRAY',
     & 'exceeded XRSM parameter --> recompile program')
      ELSE
      XRSN=XRSN+1
      END IF
      DO IAT=1,NATOM
      IF (ISLCT(IAT).EQ.1) THEN
      XRNATF=XRNATF+1
      XRINDF(XRNATF)=XRSN
      XRATOF(XRNATF)=IAT
      END IF
      END DO
      CALL NEXTF('SCATTER_A1=',XRSA(XRSN,1))
      CALL NEXTF('SCATTER_B1=',XRSB(XRSN,1))
      CALL NEXTF('SCATTER_A2=',XRSA(XRSN,2))
      CALL NEXTF('SCATTER_B2=',XRSB(XRSN,2))
      CALL NEXTF('SCATTER_A3=',XRSA(XRSN,3))
      CALL NEXTF('SCATTER_B3=',XRSB(XRSN,3))
      CALL NEXTF('SCATTER_A4=',XRSA(XRSN,4))
      CALL NEXTF('SCATTER_B4=',XRSB(XRSN,4))
      CALL NEXTF('SCATTER_C=',XRSC(XRSN))
      XRFP(XRSN)=ZERO
      XRFDP(XRSN)=ZERO
      CALL NEXTWD('XRAY>')
      IF (WD(1:4).EQ.'IMAG'.OR.WD(1:WDLEN).EQ.'FDP'
     &    .OR.WD(1:WDLEN).EQ.'SCATTER_FDP') THEN
      CALL NEXTF('FDP=',XRFDP(XRSN))
      CALL NEXTWD('XRAY>')
      IF (WD(1:4).EQ.'FP'
     &    .OR.WD(1:WDLEN).EQ.'SCATTER_FP') THEN
      CALL NEXTF('FP=',XRFP(XRSN))
      ELSE
      CALL SAVEWD
      END IF
      ELSEIF (WD(1:4).EQ.'FP'
     &    .OR.WD(1:WDLEN).EQ.'SCATTER_FP') THEN
      CALL NEXTF('FP=',XRFP(XRSN))
      CALL NEXTWD('XRAY>')
      IF (WD(1:4).EQ.'IMAG'.OR.WD(1:4).EQ.'FDP'
     &    .OR.WD(1:WDLEN).EQ.'SCATTER_FP') THEN
      CALL NEXTF('FDP=',XRFDP(XRSN))
      ELSE
      CALL SAVEWD
      END IF
      ELSE
      CALL SAVEWD
      END IF
C
C check for operlapping definitions
      DO IAT=1,NATOM
      ISLCT(IAT)=0
      END DO
      DO IAT=1,XRNATF
      ISLCT(XRATOF(IAT))=ISLCT(XRATOF(IAT))+1
      END DO
      DO IAT=1,NATOM
      IF (ISLCT(IAT).GT.1) THEN
      WRITE(6,'(9A)')
     & ' %XRAY-ERR: SCATter definition overlap for ( ',
     & SEGID(IAT),' ',RES(IAT),' ',RESID(IAT),' ',TYPE(IAT),' )'
      WRITE(6,'(A)')
     & ' The SCATter statements select at least one atom to be ',
     & ' in two or more classes.  Please check your SCATter ',
     & ' selections with respect to overlapping selections.'
      CALL WRNDIE(-1,'XRAY',
     & 'SCATter definition overlap --> change SCATter definition')
      END IF
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XFFT
C
C Parser routine for parameters of Fast Fourier computation for
C structure factors and derivatives.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
C parameters
      DOUBLE PRECISION ZERO, HALF
      PARAMETER (ZERO=0.0D0, HALF=0.5D0)
C begin
      CALL PUSEND('XFFT>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('XFFT>')
      CALL MISCOM('XFFT>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-fft')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GRID') THEN
      CALL NEXTF('GRIDsize=',GRID)
      IF (GRID.GT.HALF.OR.GRID.LT.ZERO) THEN
      WRITE(6,'(A)') ' %XFFT-ERR: grid should be <0.5 and >0.0'
      GRID=MIN(HALF,MAX(ZERO,GRID))
      END IF
      CALL XRMAPR(0)
      ELSE IF (WD(1:4).EQ.'XGRI') THEN
      CALL NEXTI('XGRIdfactor=',XGRIDF)
      ELSE IF (WD(1:4).EQ.'YGRI') THEN
      CALL NEXTI('YGRIdfactor=',YGRIDF)
      ELSE IF (WD(1:4).EQ.'ZGRI') THEN
      CALL NEXTI('ZGRIdfactor=',ZGRIDF)
      ELSE IF (WD(1:4).EQ.'BSCA') THEN
      CALL NEXTF('BSCAlefactor=',BSCAL)
      ELSE IF (WD(1:4).EQ.'ELIM') THEN
      CALL NEXTF('ELIMit=',ELIM)
      ELSE IF (WD(1:4).EQ.'MEMO') THEN
      CALL NEXTIT('MEMOry=',MEMORY)
      ELSE IF (WD(1:4).EQ.'AUTO') THEN
      CALL NEXTLO('AUTOmemory=',QMEMAUTO)
      ELSE IF (WD(1:4).EQ.'PRIM') THEN
      CALL NEXTI('PRIMe=',FTPRIM)
      ELSE IF (WD(1:4).EQ.'AVOI') THEN
      CALL NEXTI('AVOId=',FTAVOI)
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' ---------------------------',
     & '-----FFT-parameters--------------------------------'
      WRITE(6,'(3(A,F6.2),A,I9,A,L1)')
     & ' | GRID=',GRID,' BSCAle=',BSCAL,' ELIM=',ELIM,
     & ' MEMOry=',MEMORY,' AUTOmemory=',QMEMAUTO
      WRITE(6,'(2(A,I9))')
     & ' | PRIMe=',FTPRIM,' AVOId=',FTAVOI
      WRITE(6,'(A,I4)') ' | XGRIdfactor=', XGRIDF
      WRITE(6,'(A,I4)') ' | YGRIdfactor=', YGRIDF
      WRITE(6,'(A,I4)') ' | ZGRIdfactor=', ZGRIDF
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('XFFT>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
      RETURN
      END
C======================================================================
      SUBROUTINE XRQUER(XRNREF,XRH,XRK,XRL,QHERM,XRCELL,XRTR)
C
C General diffration data info.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(6), XRTR(3,3)
C local
      INTEGER I0, I1
      INTEGER REFLCT, HMAX, KMAX, LMAX, HMIN, KMIN, LMIN
      DOUBLE PRECISION D, RMIN, RMAX
C parameters
      DOUBLE PRECISION ZERO, ONE, THOUSD
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, THOUSD=1000.0D0)
C begin
C
      WRITE(6,'(2A)') ' ---------------------------',
     & 'diffraction-data-----------------------------------'
      WRITE(6,'(6(A,F6.2))')
     & ' | a=',XRCELL(1),', b=',XRCELL(2),', c=',XRCELL(3),
     & ', alpha=',XRCELL(4),', beta=',XRCELL(5),', gamma=',XRCELL(6)
C
      IF (QHERM) THEN
      WRITE(6,'(A)')
     & ' | ANOMalous=FALSe'
      ELSE
      WRITE(6,'(A)')
     & ' | ANOMalous=TRUE'
      END IF
C
C determine max and min resolution and hmax, kmax, lmax
      HMAX=0
      KMAX=0
      LMAX=0
      HMIN=0
      KMIN=0
      LMIN=0
      I0=0
      I1=0
      IF (XRNREF.EQ.0) THEN
      RMIN=ZERO
      RMAX=ZERO
      ELSE
      RMIN=THOUSD
      RMAX=-ONE
      END IF
      DO REFLCT=1,XRNREF
      HMAX=MAX(HMAX,XRH(REFLCT))
      KMAX=MAX(KMAX,XRK(REFLCT))
      LMAX=MAX(LMAX,XRL(REFLCT))
      HMIN=MIN(HMIN,XRH(REFLCT))
      KMIN=MIN(KMIN,XRK(REFLCT))
      LMIN=MIN(LMIN,XRL(REFLCT))
      CALL XRSSQ(XRH(REFLCT), XRK(REFLCT), XRL(REFLCT), D, XRTR)
      IF (D.GT.RSMALL) D=ONE/SQRT(D)
      RMIN=MIN(RMIN,D)
      RMAX=MAX(RMAX,D)
      END DO
C
      WRITE(6,'(A,I6)')
     & ' | Number or unique reflections read=',XRNREF
      WRITE(6,'(A,F10.2,A,F10.2)')
     & ' | resolution range covered by these reflections: ',
     &  RMAX,' to ',RMIN
      WRITE(6,'(6(A,I5))')
     & ' | Hmax=',HMAX,' Kmax=',KMAX,' Lmax=',LMAX,
     & ' Hmin=',HMIN,' Kmin=',KMIN,' Lmin=',LMIN
      WRITE(6,'(2A)') ' ---------------------------',
     & '---------------------------------------------------'
      RETURN
      END
C======================================================================
      SUBROUTINE XRINAT(FLAG,ANOMFLAG,NATOM)
C
C Routine initializes atom selection array to all non-H atoms.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'funct.inc'
      INTEGER FLAG(*), ANOMFLAG(*), NATOM
C local
      INTEGER IAT
C begin
      DO IAT=1,NATOM
      ANOMFLAG(IAT)=0
      IF (HYDROG(IAT)) THEN
      FLAG(IAT)=0
      ELSE
      FLAG(IAT)=1
      END IF
      END DO
      RETURN
      END
C========================================================================
      SUBROUTINE XMAKEBIN(MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE)
C
C Routine makes an array that defines resolution-dependent bins
C between the BINRESolution limits.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER BINMODE
C local
      DOUBLE PRECISION V1, V2, VSHELL
      INTEGER IND
C parameter
      DOUBLE PRECISION ZERO, ONE, THREE, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, THREE=3.0D0, FOUR=4.0D0)
C begin
C
      IF (XBINLOW.EQ.ZERO.AND.XBINHIGH.EQ.ZERO) THEN
      CALL WRNDIE(-1,'XMAKEBIN',
     & ' Fatal error -- BINResolution not specified. ')
      ELSEIF (BINMODE.EQ.1) THEN
C
C Divides the XBINLOW->XBINHIGH resolution range into MBINS shells
C of equal volume
C
C BINSHELL contains the bin borders stored as 1/d^2, where d is the
C Bragg spacing. BINSHELL entries stored as a decreasing-value array,
C first entry (BINSHELL(1)) is 1/dmin^2 and is the largest value.
C
      V1=FOUR*PI*XBINHIGH**3/THREE
      V2=FOUR*PI*XBINLOW**3/THREE
      VSHELL=(V2-V1)/MBINS
      DO IND=0,MBINS
      BINSHELL(IND+1)=
     & (((V1+IND*VSHELL)*THREE/(FOUR*PI))**(ONE/THREE))**2
      END DO
      BINSHELL(MBINS+1)=BINSHELL(MBINS+1)-RSMALL
      BINSHELL(1)=BINSHELL(1)+RSMALL
C
      IF (WRNLEV.GE.10) THEN
      IF (XBINLOW.GT.RSMALL) THEN
      WRITE(6,'(A,I8,A,/,F10.4,A,F10.4,A)')
     & ' XMAKEBIN: using ',MBINS,' equal-volume bins between',
     &  ONE/XBINHIGH,' and ',ONE/XBINLOW,' A resolution.'
      WRITE(6,'(12F6.2)') (ONE/SQRT(BINSHELL(IND)),IND=MBINS+1,1,-1)
      ELSE
      WRITE(6,'(A,I8,A,F10.4,A)')
     & ' XMAKEBIN: using ',MBINS,' equal-volume bins between',
     &  ONE/XBINHIGH,' and infinity A resolution.'
      WRITE(6,'(12F6.2)') (ONE/SQRT(BINSHELL(IND)),IND=MBINS,1,-1)
      END IF
      END IF
C
      ELSE
      CALL WRNDIE(-1,'XMAKEBIN',
     & ' Fatal error -- binmode unknown. ')
      END IF
C
      RETURN
      END
C========================================================================
      SUBROUTINE XEXPRDEF(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &                   RPNDOM,RPNLEV,DEPTH,XSFNUM,XSFNAM,XSFTYPE,
     &                   PROM,QHERM,MODE)
C
C Routine parses a syntactic expression statement for reciprocal space
C objects, such as TARGet=<expression> or TSELEection=<selection>.
C The expression or selection is converted into a RPN command stack.
C For mode='EXPR' a complex or real reciprocal space object expression
C is expected.  For mode='SELE', an reciprocal space object selection
C statement is expected.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C
      INTEGER RPNMX, RPNX
      CHARACTER*(*) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
      DOUBLE COMPLEX RPNDB(4,RPNMX)
      INTEGER RPNMLT(RPNMX), RPNLEV(RPNMX)
      CHARACTER*2 RPNTYP(RPNMX), RPNDOM(RPNMX)
      INTEGER DEPTH
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      CHARACTER*(*) PROM
      LOGICAL QHERM
      CHARACTER*(*) MODE
C local
      INTEGER I
      LOGICAL ERR
C
C definition of the functions
      EXTERNAL XDOFUNC
C
C definition of the structure factor operands
      EXTERNAL XDOOPER
C
C definition of function, operand, operation types and domains
      EXTERNAL XDOTYPE
C
      CHARACTER*2 TYPE, DOMAIN
C begin
C
      ERR=.FALSE.
C
C opening parenthesis
      CALL NEXTDO(PROM)
      IF (WD(1:1).EQ.'=') CALL NEXTDO(PROM)
      IF (WD(1:1).EQ.'?') THEN
      WRITE(6,'(A)') ' List of expression in RPN notation.'
      WRITE(6,'(2A)')
     & ' level | command     |   type |  #args | ',
     & '             constants'
      DO I=1,RPNN
      WRITE(6,'(A,I3,1X,13A,I2,4(A,F5.1,F5.1),A)') ' ',
     &         RPNLEV(I),'[',RPN(1,I)(1:4),';',
     &         RPN(2,I)(1:4),';',RPN(3,I)(1:4),';',RPN(4,I)(1:4),'] ',
     &         RPNTYP(I),' ',RPNDOM(I),' ',RPNMLT(I),
     &         ' [',DBLE(RPNDB(1,I)),DIMAG(RPNDB(1,I)),
     &         ' ;',DBLE(RPNDB(2,I)),DIMAG(RPNDB(2,I)),
     &         ' ;',DBLE(RPNDB(3,I)),DIMAG(RPNDB(3,I)),
     &         ' ;',DBLE(RPNDB(4,I)),DIMAG(RPNDB(4,I)),']'
      END DO
      ELSE
C
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR(PROM,'"(" expected')
      ERR=.TRUE.
      END IF
C
      CALL NEXTDO(PROM)
      IF (WD(1:1).EQ.')') THEN
C empty expression -> set to zero
      RPNN=0
      ELSE
      CALL SAVEWD
C
      IF (.NOT.ERR) THEN
C
C parse the selection (it is ok to use Hermitian symmetry here
C since we're not going to interpret structure factor or map
C operands)
      CALL EXRPN(PROM,
     &     RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &     RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,XDOFUNC,XDOOPER,XDOTYPE,
     &     QHERM,ERR,XSFNUM,XSFNAM,XSFTYPE,0,' ',' ')
      CALL NEXTDO(PROM)
      IF (WD(1:1).NE.')'.AND..NOT.ERR) THEN
      CALL DSPERR(PROM,
     &  'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
      IF (MODE.EQ.'EXPR'.AND.TYPE.NE.'DP'
     &   .AND.TYPE.NE.'DC'.AND..NOT.ERR) THEN
      CALL WRNDIE(-5,PROM,
     & 'Data type mismatch.  Expression must be real or complex.')
      ERR=.TRUE.
      ELSEIF
     & (MODE.EQ.'SELE'.AND.TYPE.NE.'LO'.AND..NOT.ERR)
     & THEN
      CALL WRNDIE(-5,PROM,
     & 'Data type mismatch.  Expression must be logical.')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (ERR) THEN
      RPNN=0
      DEPTH=1
      CALL WRNDIE(-5,'XEXPRDEF',
     & 'There were errors in expression.')
      END IF
C
      END IF
      END IF
C
      RETURN
      END
C====================================================
      SUBROUTINE XFLIP(XRNREF,HPH,HPK,HPL)
C
C Flips all Friedel mate indices.  The CNS reduction
C routine should be called afterwards. This operation
C will then appropriately map the phases and phase
C probability distributions.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER XRNREF, HPH(*), HPK(*), HPL(*)
C local
      INTEGER IREF
C begin
      DO IREF=1,XRNREF
      HPH(IREF)=-HPH(IREF)
      HPK(IREF)=-HPK(IREF)
      HPL(IREF)=-HPL(IREF)
      END DO
      RETURN
      END
