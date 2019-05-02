      SUBROUTINE XDO(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,
     &           ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &           XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRSYGP,XRSYIV)
C
C Structure factor and map manipulation routine. Front-end
C routine.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
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
      DOUBLE PRECISION XRTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C begin
C
C check for HELP
      CALL NEXTDO('XDO>')
      IF (WD(1:4).EQ.'HELP') THEN
      CALL CNSHELP('cns-xray-do')
      CALL CNSHELP('xray-object')
      CALL CNSHELP('xray-expression')
      CALL CNSHELP('xray-selection')
      ELSE
      CALL SAVEWD
C
      CALL XDO2(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,
     &           ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &           XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRSYGP,XRSYIV)
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XDO2(NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,MAASY,
     &           MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,
     &           ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,NRHO,IRHO,NMASK,XRMAP,
     &           XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRSYGP,XRSYIV)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'expression.inc'
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
      DOUBLE PRECISION XRTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C local
      LOGICAL ERR, OK, COND, QBYPASS, GOTIT
      INTEGER DEPTH, DEPTH2, RPNLH, I
      CHARACTER*2 TYPE, TYPE2, DOMAIN, DOMAIN2
C
C note: the RPN command stacks are stored in expression.inc
C
C left-hand-side operand
C ----------------------
      CHARACTER*(RPNX) LHS
      INTEGER LLHS
C
C type and name for LHS operand
      CHARACTER*2 ATYPE, ADOMAIN, FTYPE, FDOMAIN
      CHARACTER*(RPNX) F(4)
      INTEGER FL(4)
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
      LOGICAL QFT, QSPEC, QSPEC2
      INTEGER RPNFT, NN, FTSTART
C
      INTEGER HPSTORE, HPMAPR, HPMAPI
C
C parameter
      DOUBLE COMPLEX CZERO, ONE
      PARAMETER (ONE=1.0D0)
C
C begin
C
      CZERO=DCMPLX(0.0D0,0.0D0)
C
C set error flag
      ERR=.FALSE.
C
C opening parenthesis
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('DO>','"(" expected')
      ERR=.TRUE.
      END IF
C
C left-hand-side operand
      IF (.NOT.ERR) THEN
      CALL NEXTDO('DO>')
      CALL COPYST(LHS,RPNX,LLHS,WD,WDLEN)
C
C check if this is a valid operand
      CALL XDOOPER(OK,F,FL,.FALSE.,
     &             XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
C
      IF (.NOT.OK) THEN
      CALL DSPERR('DO>',
     &            'invalid operand for left-hand-side of equation.')
      ERR=.TRUE.
      ELSE
C
C get data and domain type
      CALL XDOTYPE(ERR,F(1),FL(1),1,0,
     &                   ATYPE,FTYPE,ADOMAIN,FDOMAIN,QHERM,
     &                   XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
      END IF
      END IF
C
C equal sign
      IF (.NOT.ERR) THEN
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.'=') THEN
      CALL DSPERR('DO>','"=" expected')
      ERR=.TRUE.
      END IF
      END IF
C
C convert right-hand-side expression into Reverse Polish Notation
      IF (.NOT.ERR) THEN
      CALL EXRPN('DO>',RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &     RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,XDOFUNC,XDOOPER,XDOTYPE,
     &     QHERM,ERR,XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
C
      END IF
C
C check if data type of lhs is compatible with expression
      IF (.NOT.ERR.AND..NOT.(
     &       FTYPE.EQ.'DC'.AND.TYPE.EQ.'DC'
     &   .OR.FTYPE.EQ.'DP'.AND.TYPE.EQ.'DP'
     &   .OR.FTYPE.EQ.'DC'.AND.TYPE.EQ.'DP')) THEN
CCC     &   .OR.FTYPE.EQ.'DP'.AND.TYPE.EQ.'DC')) THEN
      IF (FTYPE.EQ.'DP'.AND.TYPE.EQ.'DC') THEN
      WRITE(6,'(2A)')
     & ' Please use AMPLitude or REAL function to assign',
     & ' a complex number to a real number. '
      END IF
      CALL DSPERR('DO>',
     & 'Data type mismatch.  LHS incompatible with RHS of equation.')
      ERR=.TRUE.
      END IF
C
C check if domain type of lhs is compatible with expression
      IF (.NOT.ERR.AND..NOT.(
     &       FDOMAIN.EQ.'SF'.AND.DOMAIN.EQ.'SF'
     &   .OR.FDOMAIN.EQ.'SF'.AND.DOMAIN.EQ.'  '
     &   .OR.FDOMAIN.EQ.'MP'.AND.DOMAIN.EQ.'MP'
     &   .OR.FDOMAIN.EQ.'MP'.AND.DOMAIN.EQ.'  ')) THEN
      CALL DSPERR('DO>',
     & 'Domain type mismatch.  LHS incompatible with RHS of equation.')
      ERR=.TRUE.
      END IF
C
C parse closing parenthesis
      IF (.NOT.ERR) THEN
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.')') THEN
      CALL DSPERR('DO>',
     &  'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
      END IF
C
C check if a selection is specified
      IF (.NOT.ERR) THEN
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('DO>','selection expected.')
      ERR=.TRUE.
      ELSE
C
C parse the selection
      CALL EXRPN('DO>',RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,
     & RPNTYP2,RPNDOM2,RPNLEV2,TYPE2,DOMAIN2,DEPTH2,XDOFUNC,XDOOPER,
     & XDOTYPE,QHERM,ERR,XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
      DEPTH=MAX(DEPTH,DEPTH2)
C
C closing parenthesis
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.')') THEN
      CALL DSPERR('DO>',
     &  'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
C
C check data type
      IF (TYPE2.NE.'LO') THEN
      CALL DSPERR('DO>',
     & 'Data type mismatch.  Selection must be a logical expression.')
      ERR=.TRUE.
      END IF
C
C check that no FT is present in selection
      IF (.NOT.ERR) THEN
      DO NN=1,RPNN2
      ERR=ERR.OR.RPN2(1,NN)(1:RPNL2(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('DO',
     &  'No FT allowed in selections.')
      END IF
      END IF
C
C check if selection operation can be bypassed
      IF (RPNN2.EQ.1.AND.RPN2(1,1)(1:RPNL2(1,1)).EQ.'TRUE') THEN
      QBYPASS=.TRUE.
      ELSE
      QBYPASS=.FALSE.
      END IF
C
      END IF
      END IF
C
      IF (QBYPASS.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: bypassing selection.'
      END IF
C
C check if an FT operation is present
      IF (.NOT.ERR) THEN
      QFT=.FALSE.
      DO NN=1,RPNN
      IF (RPN(1,NN)(1:RPNL(1,NN)).EQ.'FT') THEN
      IF (QFT) ERR=.TRUE.
      QFT=.TRUE.
      RPNFT=NN
      END IF
      END DO
      IF (ERR) THEN
      CALL DSPERR('DO',
     &  'Only one FT operation allowed in expression.')
      QFT=.FALSE.
      ERR=.TRUE.
      END IF
      END IF
C
C check if a special structure factor
C operation is present in the selection
C or in the assignment statement.
C A special operation implies that all structure factor or map
C elements are required for the operation (e.g.,
C normalization of structure factors).
      IF (.NOT.ERR) THEN
C
      CALL XDOSPCL(RPNMX,RPNX,RPNN2,RPN2,RPNL2,QSPEC2)
      CALL XDOSPCL(RPNMX,RPNX,RPNN,RPN,RPNL,QSPEC)
      QSPEC=QSPEC.OR.QSPEC2
C
      END IF
C
C check if domain type of selection is compatible with expression
      IF (.NOT.ERR.AND..NOT.(.NOT.QFT.AND.(
     &       DOMAIN2.EQ.'SF'.AND.DOMAIN.EQ.'SF'
     &   .OR.DOMAIN2.EQ.'SF'.AND.DOMAIN.EQ.'  '
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'SF'
     &   .OR.DOMAIN2.EQ.'MP'.AND.DOMAIN.EQ.'MP'
     &   .OR.DOMAIN2.EQ.'MP'.AND.DOMAIN.EQ.'  '
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'MP'
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'  ') .OR.
     &         ( QFT.AND.(
     &       DOMAIN2.EQ.'SF'.AND.DOMAIN.EQ.'MP'
     &   .OR.DOMAIN2.EQ.'SF'.AND.DOMAIN.EQ.'  '
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'SF'
     &   .OR.DOMAIN2.EQ.'MP'.AND.DOMAIN.EQ.'SF'
     &   .OR.DOMAIN2.EQ.'MP'.AND.DOMAIN.EQ.'  '
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'MP'
     &   .OR.DOMAIN2.EQ.'  '.AND.DOMAIN.EQ.'  ')))) THEN
      CALL DSPERR('DO>',
     & 'Domain type mismatch. Expression incompatible with selection.')
      ERR=.TRUE.
      END IF
C
      IF (.NOT.ERR.AND.QFT) THEN
C
C create third command stack (operations left and right of FT())
C
C determine starting point
      NN=RPNFT
      FTSTART=1
      GOTIT=.FALSE.
      DO WHILE (NN.GT.1.AND..NOT.GOTIT)
      NN=NN-1
      GOTIT=RPNLEV(NN).LT.RPNLEV(RPNFT)
      END DO
      IF (GOTIT) FTSTART=NN+1
C
C copy all statements right of FT() or special operation
      RPNN3=0
      DO NN=1,FTSTART-1
      RPNN3=RPNN3+1
      RPN3(1,RPNN3)=RPN(1,NN)
      RPN3(2,RPNN3)=RPN(2,NN)
      RPN3(3,RPNN3)=RPN(3,NN)
      RPN3(4,RPNN3)=RPN(4,NN)
      RPNL3(1,RPNN3)=RPNL(1,NN)
      RPNL3(2,RPNN3)=RPNL(2,NN)
      RPNL3(3,RPNN3)=RPNL(3,NN)
      RPNL3(4,RPNN3)=RPNL(4,NN)
      RPNDB3(1,RPNN3)=RPNDB(1,NN)
      RPNDB3(2,RPNN3)=RPNDB(2,NN)
      RPNDB3(3,RPNN3)=RPNDB(3,NN)
      RPNDB3(4,RPNN3)=RPNDB(4,NN)
      RPNMLT3(RPNN3)=RPNMLT(NN)
      RPNTYP3(RPNN3)=RPNTYP(NN)
      RPNDOM3(RPNN3)=RPNDOM(NN)
      RPNLEV3(RPNN3)=RPNLEV(NN)
      END DO
      RPNN3=RPNN3+1
C
C put LHS
      RPNLH=RPNN3
      RPN3(1,RPNN3)=LHS
      RPNL3(1,RPNN3)=LLHS
      RPNMLT3(RPNN3)=0
      RPNTYP3(RPNN3)=FTYPE
      RPNDOM3(RPNN3)=FDOMAIN
      RPNLEV3(RPNN3)=1
C
C copy all statements left of FT() or special operation
      DO NN=RPNFT+1,RPNN
      RPNN3=RPNN3+1
      RPN3(1,RPNN3)=RPN(1,NN)
      RPN3(2,RPNN3)=RPN(2,NN)
      RPN3(3,RPNN3)=RPN(3,NN)
      RPN3(4,RPNN3)=RPN(4,NN)
      RPNL3(1,RPNN3)=RPNL(1,NN)
      RPNL3(2,RPNN3)=RPNL(2,NN)
      RPNL3(3,RPNN3)=RPNL(3,NN)
      RPNL3(4,RPNN3)=RPNL(4,NN)
      RPNDB3(1,RPNN3)=RPNDB(1,NN)
      RPNDB3(2,RPNN3)=RPNDB(2,NN)
      RPNDB3(3,RPNN3)=RPNDB(3,NN)
      RPNDB3(4,RPNN3)=RPNDB(4,NN)
      RPNMLT3(RPNN3)=RPNMLT(NN)
      RPNTYP3(RPNN3)=RPNTYP(NN)
      RPNDOM3(RPNN3)=RPNDOM(NN)
      RPNLEV3(RPNN3)=RPNLEV(NN)
      END DO
C
C copy all operations between FTSTART and RPNFT-1 into RPN
      RPNN=0
      DO NN=FTSTART,RPNFT-1
      RPNN=RPNN+1
      RPN(1,RPNN)=RPN(1,NN)
      RPN(2,RPNN)=RPN(2,NN)
      RPN(3,RPNN)=RPN(3,NN)
      RPN(4,RPNN)=RPN(4,NN)
      RPNL(1,RPNN)=RPNL(1,NN)
      RPNL(2,RPNN)=RPNL(2,NN)
      RPNL(3,RPNN)=RPNL(3,NN)
      RPNL(4,RPNN)=RPNL(4,NN)
      RPNDB(1,RPNN)=RPNDB(1,NN)
      RPNDB(2,RPNN)=RPNDB(2,NN)
      RPNDB(3,RPNN)=RPNDB(3,NN)
      RPNDB(4,RPNN)=RPNDB(4,NN)
      RPNMLT(RPNN)=RPNMLT(NN)
      RPNTYP(RPNN)=RPNTYP(NN)
      RPNDOM(RPNN)=RPNDOM(NN)
      RPNLEV(RPNN)=RPNLEV(NN)
      END DO
C
      END IF
C
C evaluate expression
      IF (.NOT.ERR) THEN
C
      IF (.NOT.QFT.AND.FDOMAIN.EQ.'SF') THEN
C
C----------------------------------------------------------------------
C
C no FT, operation is on structure factors
      CALL XDOSF(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &         RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &         RPNMLT2,RPNLEV2,
     &         RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,QBYPASS,
     &         XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &         XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &         HPMULT,HPTYPE,
     &         HPSTORE,.FALSE.,0,
     &         QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,QSPEC,
     &         XRCELL)
C
C----------------------------------------------------------------------
C
      ELSEIF (.NOT.QFT.AND.FDOMAIN.EQ.'MP') THEN
C
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C
C no FT, operation is on maps
      CALL XDOMP(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &           RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &           RPNMLT2,RPNLEV2,
     &           RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,
     &           NMASK,NA,NB,NC,MAASY,MBASY,MCASY,NAASY,
     &           NBASY,NCASY,QHERM,NRHO,IRHO,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,HPMAPR,HPMAPI,XRNSYM,QSPEC)
C
C----------------------------------------------------------------------
C
      ELSEIF (QFT.AND.FDOMAIN.EQ.'MP') THEN
C
C FT, structure factor to map
C
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: executing operations right of FT.'
      END IF
C
C allocate space for STORE and fill with zeros
      CALL XSFAL(HPSTORE,XRNREF,'SF')
C
C first operation is on structure factors. Store result in
C "STORE"
      CALL XDOSF(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &        RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &        RPNMLT2,RPNLEV2,
     &        RPNTYP2,RPNDOM2,DEPTH,'STORE',5,QBYPASS,
     &        XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &        XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &        HPMULT,HPTYPE,
     &        HPSTORE,.FALSE.,0,
     &        QHERM,XRNSYM,XRMSYM,XRSYTH,
     &        XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,QSPEC,
     &        XRCELL)
C
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,
     &                   NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C
C get the pointer for the LHS operand.
C check real space objects
      DO I=1,XRHONUM
      IF (LHS(1:LLHS).EQ.XRHONAM(I)) THEN
C
C allocate space for object if required
      IF (HPRRHO(I).EQ.0) THEN
      CALL XMAPAL(HPRRHO(I),HPIRHO(I),QHERM,NRHO,IRHO)
      END IF
      HPMAPR=HPRRHO(I)
      HPMAPI=HPIRHO(I)
      END IF
      END DO
C
C do FFT
      CALL XDOFT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     HEAP(HPRHOMA),HEAP(HPMAPR),HEAP(HPMAPI),
     &     XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPSTORE),
     &     NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
      CALL FREHP(HPSTORE,ICPLX8(XRNREF))
C
C bypass next operation if just operand
      IF (RPNN3.GT.1) THEN
C
C set selection to TRUE
      RPNN2=1
      RPN2(1,RPNN2)='TRUE'
      RPNL2(1,RPNN2)=4
      RPNMLT2(RPNN2)=0
      RPNTYP2(RPNN2)='LO'
      RPNDOM2(RPNN2)='  '
      RPNLEV2(RPNN2)=1
      CALL XDOMP(RPNMX,RPNX,RPN3,RPNL3,RPNN3,RPNDB3,RPNMLT3,RPNLEV3,
     &           RPNTYP3,RPNDOM3,RPN2,RPNL2,RPNN2,RPNDB2,
     &           RPNMLT2,RPNLEV2,
     &           RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,
     &           NMASK,NA,NB,NC,MAASY,MBASY,MCASY,NAASY,
     &           NBASY,NCASY,QHERM,NRHO,IRHO,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,HPMAPR,HPMAPI,XRNSYM,QSPEC)
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: executing operations left of FT.'
      END IF
C
      ELSE
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: bypassing operations left of FT.'
      END IF
      END IF
C
C----------------------------------------------------------------------
C
      ELSEIF (QFT.AND.FDOMAIN.EQ.'SF') THEN
C
C FT, map to structure factor
C
C bypass operation if it just consists of an operand
C get pointer for RHS operand
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (LHS(1:LLHS).EQ.XRHONAM(I)) THEN
C
C allocate space for object if required
      IF (HPRRHO(I).EQ.0) THEN
      WRITE(6,'(3A)') ' %XDO-ERR: real space object ',
     & LHS(1:LLHS),' undefined.'
      CALL WRNDIE(-5,'XDO','object undefined.')
      ELSE
      COND=.TRUE.
      HPMAPR=HPRRHO(I)
      HPMAPI=HPIRHO(I)
      END IF
      END IF
      END DO
C
      IF (.NOT.COND) THEN
C
C allocate space for the STORE map and zero fill the map.
      CALL XMAPAL(HPMAPR,HPMAPI,QHERM,NRHO,IRHO)
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: executing operations right of FT.'
      END IF
C
C first operation is on structure factors. Store result in
C "STORE"
      CALL XDOMP(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &           RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &           RPNMLT2,RPNLEV2,
     &           RPNTYP2,RPNDOM2,DEPTH,'STORE',5,
     &           NMASK,NA,NB,NC,MAASY,MBASY,MCASY,NAASY,
     &           NBASY,NCASY,QHERM,NRHO,IRHO,
     &           XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &           HPRHOMA,HPMAPR,HPMAPI,XRNSYM,QSPEC)
      END IF
C
      IF (COND.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: bypassing operations right of FT.'
      END IF
C
C get the pointer for the LHS operand.  Allocate space
C if not already done.
C
C use STORE
      RPN3(1,RPNLH)='STORE'
      RPNL3(1,RPNLH)=5
C
C allocate space for STORE and fill with zeros
      CALL XSFAL(HPSTORE,XRNREF,'SF')
C
C FT map to structure factor
C
C Fourier transformation of the map.  Result will be HPSTORE.
      CALL XDOIFT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     HEAP(HPRHOMA),HEAP(HPMAPR),HEAP(HPMAPI),
     &     XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPSTORE),
     &     NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
      IF (.NOT.COND) THEN
      CALL FREHP(HPMAPR,IREAL4(NRHO))
      CALL FREHP(HPMAPI,IREAL4(IRHO))
      END IF
C
C bypass next operation if just operand
      IF (RPNN3.GT.1.OR.RPN3(1,RPNLH).EQ.'STORE') THEN
C
C set selection to TRUE
      RPNN2=1
      RPN2(1,RPNN2)='TRUE'
      RPNL2(1,RPNN2)=4
      RPNMLT2(RPNN2)=0
      RPNTYP2(RPNN2)='LO'
      RPNDOM2(RPNN2)='  '
C
      CALL XDOSF(RPNMX,RPNX,RPN3,RPNL3,RPNN3,RPNDB3,RPNMLT3,RPNLEV3,
     &         RPNTYP3,RPNDOM3,RPN2,RPNL2,RPNN2,RPNDB2,
     &         RPNMLT2,RPNLEV2,
     &         RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,QBYPASS,
     &         XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &         XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &         HPMULT,HPTYPE,
     &         HPSTORE,.FALSE.,0,
     &         QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,QSPEC,
     &         XRCELL)
C
      IF (RPN3(1,RPNLH).EQ.'STORE') THEN
      CALL FREHP(HPSTORE,ICPLX8(XRNREF))
      END IF
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: executing operations left of FT.'
      END IF
C
      ELSE
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' DO: bypassing operations left of FT.'
      END IF
      END IF
C
C---------------------------------------------------------------------- C
      END IF
      END IF
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'SHOW',
     & 'There were errors in DO expression.')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSF(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &           RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &           RPNMLT2,RPNLEV2,RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,
     &           QBYPASS,
     &           XRTR,XRVOL,XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,
     &           HPSTORE,QSELE,
     &           SELE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           QSPEC,XRCELL)
C
C Routine evaluates a structure factor expression (stored
C in the RPN command stack) for selected elements (selection
C stored in the RPN2 command stack) and assignes it to
C specified left-hand-side operand (LHS, LLHS).
C
C If QSELE is set then store the selection in array SELE.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER RPNMX, RPNX
      CHARACTER*(*) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
      DOUBLE COMPLEX RPNDB(4,RPNMX)
      INTEGER RPNMLT(RPNMX), RPNLEV(RPNMX)
      CHARACTER*2 RPNTYP(RPNMX), RPNDOM(RPNMX)
      CHARACTER*(*) RPN2(4,RPNMX)
      INTEGER RPNL2(4,RPNMX), RPNN2
      DOUBLE COMPLEX RPNDB2(4,RPNMX)
      INTEGER RPNMLT2(RPNMX), RPNLEV2(RPNMX)
      CHARACTER*2 RPNTYP2(RPNMX), RPNDOM2(RPNMX)
      INTEGER DEPTH
      CHARACTER*(*) LHS
      INTEGER LLHS
      LOGICAL QBYPASS
      DOUBLE PRECISION XRTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, HPSTORE
      LOGICAL QSELE
      INTEGER SELE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL QSPEC
      DOUBLE PRECISION XRCELL(6)
C local
      INTEGER MSTACK, START, STOP, VLEVEL, NSELE, NINDEX, NNSELE, I
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C pointer
      INTEGER INDEX, VSTACK, LSTACK
C begin
C
      IF (QSELE) THEN
      DO I=1,XRNREF
      SELE(I)=0
      END DO
      END IF
C
      NNSELE=0
C
C allocate space for the variable stack.  We perform the operations
C in junks of MSTACK elements
C
      IF (QBYPASS) THEN
      NINDEX=XRNREF
      ELSE
      NINDEX=XRNREF
      END IF
C
C For special operations we have to allocate space for
C a stack that is big enough to hold all structure factor
C elements.  For all other operations we can do them
C in smaller junks.
C
      IF (QSPEC) THEN
      MSTACK=NINDEX
      ELSE
      MSTACK=MIN(NINDEX,10000)
      END IF
C
      INDEX=ALLHP(INTEG4(MSTACK))
      VSTACK=ALLHP(ICPLX8(MSTACK*DEPTH))
      LSTACK=ALLHP(ILOGIC(MSTACK*DEPTH))
C
C
Cbegin loop over junks  (size MSTACK, from START to STOP)
C---------------------
      START=1
      STOP=MIN(NINDEX,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP
      CALL XDOMAKD(HEAP(INDEX),START,STOP,NSELE)
C
C overwrite the selection if explicit selection is specified
      IF (.NOT.QBYPASS) THEN
C
C evaluate selection
      CALL XDOEVAL(RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,RPNLEV2,
     &     RPNTYP2,RPNDOM2,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,
     &     XRVOL)
C
C reduce the index according to the selection
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NSELE,HEAP(INDEX))
      END IF
C
C if QSELE is set then store the selection in SELE
      IF (QSELE) THEN
      CALL XDOMAKL(NSELE,SELE,HEAP(INDEX))
      END IF
C
      NNSELE=NNSELE+NSELE
C
C evaluate expression
      CALL XDOEVAL(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNLEV,
     &     RPNTYP,RPNDOM,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,
     &     XRVOL)
C
C assign to LHS
      CALL XDOASSN(LHS,LLHS,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF)
C
      START=START+MSTACK
      STOP=MIN(NINDEX,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I9,A)')
     & ' Total of ',NNSELE,' structure factor elements were selected.'
      END IF
      DBPREC=NNSELE
      CALL DECLAR( 'SELECT','DP',' ',DBCOMP,DBPREC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOMAKL(N,SELE,INDEX)
C
C fill selection array.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER N, SELE(*), INDEX(*)
C local
      INTEGER I
C begin
      DO I=1,N
      SELE(INDEX(I))=1
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOMP(RPNMX,RPNX,RPN,RPNL,RPNN,RPNDB,RPNMLT,RPNLEV,
     &                 RPNTYP,RPNDOM,RPN2,RPNL2,RPNN2,RPNDB2,
     &                 RPNMLT2,RPNLEV2,RPNTYP2,RPNDOM2,DEPTH,LHS,LLHS,
     &                 NMASK,NA,NB,NC,MAASY,MBASY,MCASY,NAASY,
     &                 NBASY,NCASY,QHERM,NRHO,IRHO,
     &                 XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &                 HPRHOMA,HPMAPR,HPMAPI,XRNSYM,QSPEC)
C
C Routine evaluates a map expression (stored
C in the RPN command stack) for selected elements (selection
C stored in the RPN2 command stack) and assignes it to
C specified left-hand-side operand (LHS, LLHS).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER RPNMX, RPNX
      CHARACTER*(*) RPN(4,RPNMX)
      INTEGER RPNL(4,RPNMX), RPNN
      DOUBLE COMPLEX RPNDB(4,RPNMX)
      INTEGER RPNMLT(RPNMX), RPNLEV(RPNMX)
      CHARACTER*2 RPNTYP(RPNMX), RPNDOM(RPNMX)
      CHARACTER*(*) RPN2(4,RPNMX)
      INTEGER RPNL2(4,RPNMX), RPNN2
      DOUBLE COMPLEX RPNDB2(4,RPNMX)
      INTEGER RPNMLT2(RPNMX), RPNLEV2(RPNMX)
      CHARACTER*2 RPNTYP2(RPNMX), RPNDOM2(RPNMX)
      INTEGER DEPTH
      CHARACTER*(*) LHS
      INTEGER LLHS
      INTEGER NMASK
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      INTEGER NRHO, IRHO
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, HPMAPR, HPMAPI, XRNSYM
      LOGICAL QSPEC
C local
      INTEGER MSTACK, A, B, C, START, STOP, NSELE, NN, VLEVEL, NNSELE
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C pointer
      INTEGER INDEXA, INDEXB, INDEXC, VSTACK, LSTACK
C begin
C
      NNSELE=0
C
C For special operations we have to allocate space for
C a stack that is big enough to hold all map
C elements.  For all other operations we can do them
C in smaller junks.
C
      IF (QSPEC) THEN
      MSTACK=NMASK
      ELSE
      MSTACK=1000
      END IF
C
C allocate space for the variable stack.  We perform the operations
C in junks of MSTACK elements
      INDEXA=ALLHP(INTEG4(MSTACK))
      INDEXB=ALLHP(INTEG4(MSTACK))
      INDEXC=ALLHP(INTEG4(MSTACK))
      VSTACK=ALLHP(ICPLX8(MSTACK*DEPTH))
      LSTACK=ALLHP(ILOGIC(MSTACK*DEPTH))
C
Cbegin loop over junks  (size MSTACK, from START to STOP)
C---------------------
      A=MAASY-1
      B=MBASY
      C=MCASY
      START=1
      STOP=MIN(NMASK,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP
      CALL XMDOMAK(A,B,C,HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),
     &           START,STOP,NSELE,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           HEAP(HPRHOMA))
C
C
C evaluate selection
      CALL XMDOEVA(RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,RPNLEV2,
     &             VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &             HEAP(LSTACK),NA,NB,NC,
     &             MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             QHERM,
     &             XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &             HPMAPR,HPMAPI,HPRHOMA,
     &             HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),XRNSYM)
C
C reduce the index according to the selection
      NN=NSELE
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NN,HEAP(INDEXA))
      NN=NSELE
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NN,HEAP(INDEXB))
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NSELE,HEAP(INDEXC))
C
      NNSELE=NNSELE+NSELE
C
C evaluate expression
      CALL XMDOEVA(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNLEV,
     &             VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &             HEAP(LSTACK),NA,NB,NC,
     &             MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             QHERM,
     &             XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &             HPMAPR,HPMAPI,HPRHOMA,
     &             HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),XRNSYM)
C
C assign to LHS
      CALL XMDOASS(LHS,LLHS,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     &             HEAP(LSTACK),
     &             MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             QHERM,NRHO,IRHO,
     &             XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &             HPMAPR,HPMAPI,
     &             HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC))
C
      START=START+MSTACK
      STOP=MIN(NMASK,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEXC,INTEG4(MSTACK))
      CALL FREHP(INDEXB,INTEG4(MSTACK))
      CALL FREHP(INDEXA,INTEG4(MSTACK))
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I9,A)')
     & ' Total of ',NNSELE,' map elements were selected.'
      END IF
      DBPREC=NNSELE
      CALL DECLAR( 'SELECT','DP',' ',DBCOMP,DBPREC)
C
      RETURN
      END
