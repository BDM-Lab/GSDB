      SUBROUTINE XSHOW(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRNSYM,XRMSYM,XRSYTH,
     &      XRSYMM,XRITSY,XRHONUM,XRHONAM,
     &      HPRRHO,HPIRHO,HPRHOMA,NRHO,NMASK,
     &      XRTR,XRNREF,
     &      HPH,HPK,HPL,
     &      XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &      HPMULT,HPTYPE,
     &      XRMREF,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C Parsing routine for SHOW facility for structure factors
C and maps.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'expression.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, NMASK, NRHO
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, XRMREF
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      LOGICAL ERR, QSPEC, QSPEC2
      INTEGER DEPTH, DEPTH2, NSELE, NN, START, STOP
      CHARACTER*2 TYPE, TYPE2, DOMAIN, DOMAIN2
      CHARACTER*4 SHOW
C
C note: command stacks are stored in expression.inc
C
C definition of the functions, operands, data, and domain types
      EXTERNAL XDOFUNC, XDOOPER, XDOTYPE
C
C variable stack pointers
      INTEGER VLEVEL, VSTACK, LSTACK
C index pointer
      INTEGER INDEX, INDEXA, INDEXB, INDEXC, A, B, C
C
C size of the variable stack
      INTEGER MSTACK
C
      INTEGER DUMMY
C
C show parameters
      DOUBLE COMPLEX SHOW1,SHOW2
      INTEGER NSHOW
C
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C set error flag
      ERR=.FALSE.
C
C check for HELP
      CALL NEXTDO('SHOW>')
      IF (WD(1:4).EQ.'HELP') THEN
      CALL CNSHELP('cns-xray-show')
      CALL CNSHELP('xray-object')
      CALL CNSHELP('xray-expression')
      CALL CNSHELP('xray-selection')
      ELSE
      CALL SAVEWD
C
C get mode for SHOW option
      CALL NEXTDO('SHOW>')
      IF (WD(1:3).EQ.'SUM') THEN
      SHOW='SUM'
      ELSEIF (WD(1:3).EQ.'MAX') THEN
      SHOW='MAX'
      ELSEIF (WD(1:4).EQ.'ELEM') THEN
      SHOW='ELEM'
      ELSEIF (WD(1:3).EQ.'MIN') THEN
      SHOW='MIN'
      ELSEIF (WD(1:3).EQ.'AVE'.OR.WD(1:4).EQ.'MEAN') THEN
      SHOW='AVE'
      ELSEIF (WD(1:4).EQ.'NORM') THEN
      SHOW='NORM'
      ELSEIF (WD(1:3).EQ.'RMS') THEN
      SHOW='RMS'
      ELSEIF (WD(1:4).EQ.'ELEM') THEN
      SHOW='ELEM'
      ELSEIF (WD(1:4).EQ.'(') THEN
      SHOW='ELEM'
      CALL SAVEWD
      ELSE
      CALL DSPERR('XSHOW','invalid show option.')
      ERR=.TRUE.
      END IF
C
C opening parenthesis
      IF (.NOT.ERR) THEN
      CALL NEXTDO('SHOW>')
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('SHOW>','"(" expected')
      ERR=.TRUE.
      END IF
      END IF
C
C
C convert right-hand-side expression into Reverse Polish Notation
      IF (.NOT.ERR) THEN
C
      CALL EXRPN('SHOW>',RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     & RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,XDOFUNC,XDOOPER,XDOTYPE,
     & QHERM,ERR,
     & XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
C
      END IF
C
C
C check if data type of expression is compatible with SHOW option
C
C Note: MIN, MAX, and RMS are performed on the real part of the complex number
C
      IF (.NOT.ERR.AND. .NOT.(
     &     SHOW.EQ.'MIN' .AND.TYPE.EQ.'DP'
     & .OR.SHOW.EQ.'MAX' .AND.TYPE.EQ.'DP'
     & .OR.SHOW.EQ.'SUM' .AND.(TYPE.EQ.'DP'.OR.TYPE.EQ.'DC')
     & .OR.SHOW.EQ.'AVE' .AND.(TYPE.EQ.'DP'.OR.TYPE.EQ.'DC')
     & .OR.SHOW.EQ.'NORM'.AND.(TYPE.EQ.'DP'.OR.TYPE.EQ.'DC')
     & .OR.SHOW.EQ.'RMS' .AND.TYPE.EQ.'DP'
     & .OR.SHOW.EQ.'ELEM'.AND.(TYPE.EQ.'DP'.OR.TYPE.EQ.'DC') ) ) THEN
      CALL DSPERR('SHOW','Expression incompatible with SHOW option.')
      ERR=.TRUE.
      END IF
C
C check that no FT operation is present
      IF (.NOT.ERR) THEN
      DO NN=1,RPNN
      ERR=ERR.OR.RPN(1,NN)(1:RPNL(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('SHOW','No FT operations allowed for SHOW.')
      END IF
      END IF
C
C parse closing parenthesis
      IF (.NOT.ERR) THEN
      CALL NEXTDO('SHOW')
      IF (WD(1:1).NE.')') THEN
      CALL DSPERR('SHOW',
     & 'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
      END IF
C
C opening parenthesis
      IF (.NOT.ERR) THEN
      CALL NEXTDO('SHOW>')
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('SHOW>','"(" expected')
      ERR=.TRUE.
      END IF
      END IF
C
      IF (.NOT.ERR) THEN
C
C parse the selection
      CALL EXRPN('SHOW',RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,
     & RPNTYP2,RPNDOM2,RPNLEV2,TYPE2,DOMAIN2,DEPTH2,XDOFUNC,
     & XDOOPER,XDOTYPE,QHERM,ERR,
     & XSFNUM,XSFNAM,XSFTYPE,XRHONUM,XRHONAM,' ')
      DEPTH=MAX(DEPTH,DEPTH2)
C
      CALL NEXTDO('SHOW>')
      IF (WD(1:1).NE.')'.AND..NOT.ERR) THEN
      CALL DSPERR('SHOW>',
     & 'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
C
C check data type compatibility
      IF (TYPE2.NE.'LO') THEN
      CALL DSPERR('SHOW',
     & 'Data type mismatch.  Selection must be a logical expression.')
      ERR=.TRUE.
      ELSEIF (DOMAIN.EQ.'SF'.AND.DOMAIN2.EQ.'MP') THEN
      CALL DSPERR('SHOW',
     &'Domain mismatch: structure factor expression vs map selection.')
      ERR=.TRUE.
      ELSEIF (DOMAIN.EQ.'MP'.AND.DOMAIN2.EQ.'SF') THEN
      CALL DSPERR('SHOW',
     & 'Domain mismatch: map selection vs struct. factor expression.')
      ERR=.TRUE.
      ELSEIF (DOMAIN.EQ.'  '.AND.DOMAIN2.EQ.'  ') THEN
      CALL DSPERR('SHOW',
     & 'Domain mismatch: neither maps or structure factors present.')
      ERR=.TRUE.
      END IF
      END IF
C
C check that no FT operation is present
      IF (.NOT.ERR) THEN
      DO NN=1,RPNN2
      ERR=ERR.OR.RPN2(1,NN)(1:RPNL2(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('SHOW','No FT operations allowed for SHOW.')
      END IF
      END IF
C
C evaluate expression
      IF (.NOT.ERR) THEN
C
C check if special structure factor or map operations are present
      CALL XDOSPCL(RPNMX,RPNX,RPNN2,RPN2,RPNL2,QSPEC2)
      CALL XDOSPCL(RPNMX,RPNX,RPNN,RPN,RPNL,QSPEC)
      QSPEC=QSPEC.OR.QSPEC2
C
C initialize parameters for SHOW option
      CALL XDOSHOI(SHOW,SHOW1,SHOW2,NSHOW)
C
      IF (DOMAIN.EQ.'MP'.OR.DOMAIN2.EQ.'MP') THEN
C
C evaluate show option in map domain
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
     & START,STOP,NSELE,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & HEAP(HPRHOMA))
C
C
C evaluate selection
      CALL XMDOEVA(RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,RPNLEV2,
     & VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     & HEAP(LSTACK),NA,NB,NC,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     & DUMMY,DUMMY,HPRHOMA,
     & HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),XRNSYM)
C
C reduce the index according to the selection
      NN=NSELE
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NN,HEAP(INDEXA))
      NN=NSELE
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NN,HEAP(INDEXB))
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NSELE,HEAP(INDEXC))
C
C evaluate expression
      CALL XMDOEVA(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNLEV,
     & VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     & HEAP(LSTACK),NA,NB,NC,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,DUMMY,DUMMY,HPRHOMA,
     & HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),XRNSYM)
C
C evaluate properties for SHOW option
      CALL XMDOSHO(SHOW,SHOW1,SHOW2,NSHOW,TYPE,VLEVEL,DEPTH,NSELE,
     & HEAP(VSTACK),HEAP(LSTACK),
     & HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & HEAP(HPRHOMA))
C
      START=START+MSTACK
      STOP=MIN(NMASK,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEXC,INTEG4(MSTACK))
      CALL FREHP(INDEXB,INTEG4(MSTACK))
      CALL FREHP(INDEXA,INTEG4(MSTACK))
C
      ELSE
C
C evaluate show option in structure factor domain
C
C allocate space for the variable stack.  We perform the operations
C in junks of MSTACK elements except when special
C structure factor operations are present.
      IF (QSPEC) THEN
      MSTACK=XRNREF
      ELSE
      MSTACK=MIN(XRNREF,10000)
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
      STOP=MIN(XRNREF,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP
      CALL XDOMAKD(HEAP(INDEX),START,STOP,NSELE)
C
C evaluate selection
      CALL XDOEVAL(RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,RPNLEV2,
     & RPNTYP2,RPNDOM2,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     & HEAP(LSTACK),HEAP(INDEX),XRTR,
     & HPH,HPK,HPL,
     & XSFNUM,XSFNAM,XSFTYPE,HPSF,
     & HPMULT,HPTYPE,
     & DUMMY,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C
C reduce the index according to the selection
      CALL XDOMAKV(VLEVEL,DEPTH,HEAP(LSTACK),NSELE,HEAP(INDEX))
C
C evaluate expression
      CALL XDOEVAL(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNLEV,
     & RPNTYP,RPNDOM,VLEVEL,DEPTH,NSELE,HEAP(VSTACK),
     & HEAP(LSTACK),HEAP(INDEX),XRTR,
     & HPH,HPK,HPL,
     & XSFNUM,XSFNAM,XSFTYPE,HPSF,
     & HPMULT,HPTYPE,
     & DUMMY,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C
C evaluate properties for SHOW option
      CALL XDOSHOW(SHOW,SHOW1,SHOW2,NSHOW,TYPE,VLEVEL,DEPTH,NSELE,
     & HEAP(VSTACK),HEAP(LSTACK),HEAP(INDEX),HEAP(HPH),
     & HEAP(HPK),HEAP(HPL))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
      END IF
C
C print/declare variables for SHOW option
      CALL XDOSHOF(SHOW,SHOW1,SHOW2,NSHOW,TYPE)
C
      END IF
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'SHOW',
     & 'There were errors in SHOW expression or selection.')
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSHOI(SHOW,SHOW1,SHOW2,NSHOW)
C
C initializes variables for SHOW
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      CHARACTER*(*) SHOW
      DOUBLE COMPLEX SHOW1, SHOW2
      INTEGER NSHOW
C local
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C always initalize SHOW2 to zero
      SHOW2 = DCMPLX(ZERO,ZERO)
C
      NSHOW=0
      IF (SHOW.EQ.'SUM') THEN
      SHOW1=DCMPLX(ZERO,ZERO)
      ELSEIF (SHOW.EQ.'AVE') THEN
      SHOW1=DCMPLX(ZERO,ZERO)
      ELSEIF (SHOW.EQ.'MIN') THEN
      SHOW1=R4BIG
      ELSEIF (SHOW.EQ.'MAX') THEN
      SHOW1=-R4BIG
      ELSEIF (SHOW.EQ.'NORM') THEN
      SHOW1=DCMPLX(ZERO,ZERO)
      ELSEIF (SHOW.EQ.'RMS') THEN
      SHOW1=DCMPLX(ZERO,ZERO)
      ELSEIF (SHOW.EQ.'ELEM') THEN
      CONTINUE
      ELSE
      CALL WRNDIE(-5,'SHOW','Internal error (routine XDOSHOI)')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSHOF(SHOW,SHOW1,SHOW2,NSHOW,TYPE)
C
C Computes final show properties, writes to PUNIT (depending on
C the message level), and declares symbols.
C
C TYPE is the data type of the expression.
C
C modification: WRITEs are done depending on WRNLEV, ATB, 12/04/08
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'timer.inc'
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) SHOW
      DOUBLE COMPLEX SHOW1, SHOW2
      INTEGER NSHOW
      CHARACTER*(*) TYPE
C local
      DOUBLE PRECISION AMP, PHAS, DBPREC
      DOUBLE COMPLEX DBCOMP
C parameters
      DOUBLE PRECISION S180, ZERO
      PARAMETER (S180=180.D0, ZERO=0.0D0)
C begin
      IF (SHOW.EQ.'SUM') THEN
      IF (TYPE.EQ.'DP') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Sum of ',NSHOW,' elements = ',DBLE(SHOW1)
      END IF
      CALL DECLAR( 'RESULT','DP',' ',DBCOMP,DBLE(SHOW1))
      ELSE
      CALL XPHASE(SHOW1,AMP,PHAS)
      PHAS=PHAS*S180/PI
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4,A,F22.4)')
     & ' Sum of ',NSHOW,' elements.  Amplitude=',AMP,
     & '  Phase=',PHAS
      END IF
      CALL DECLAR( 'RESULT','DC',' ',SHOW1,DBPREC)
      CALL DECLAR( 'AMPLITUDE','DP',' ',DBCOMP,AMP)
      CALL DECLAR( 'PHASE','DP',' ',DBCOMP,PHAS)
      END IF
C
      ELSEIF (SHOW.EQ.'AVE') THEN
      IF (NSHOW.EQ.0) THEN
      CALL WRNDIE(+5,'SHOW','Number of selected elements is zero.')
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, ZERO )
      ELSE
      IF (TYPE.EQ.'DP') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' SHOW: average of ',NSHOW,' elements= ',DBLE(SHOW1)/NSHOW
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBLE(SHOW1)/NSHOW)
      ELSE
      CALL XPHASE(SHOW1/NSHOW,AMP,PHAS)
      PHAS=PHAS*S180/PI
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4,A,F22.4)')
     & ' Average of ',NSHOW,' elements.  Amplitude=',AMP,
     & '  Phase=',PHAS
      END IF
      CALL DECLAR( 'RESULT', 'DC', ' ', SHOW1/NSHOW, DBPREC )
      CALL DECLAR( 'AMPLITUDE', 'DP', ' ', DBCOMP, AMP )
      CALL DECLAR( 'PHASE', 'DP', ' ', DBCOMP, PHAS )
      END IF
      END IF
C
      ELSEIF (SHOW.EQ.'MIN') THEN
      IF (NSHOW.NE.0) THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Minimum of ',NSHOW,' elements = ',DBLE(SHOW1)
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBLE(SHOW1))
      ELSE
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Minimum of ',0,' elements = ',ZERO
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,ZERO)
      END IF
C
      ELSEIF (SHOW.EQ.'MAX') THEN
      IF (NSHOW.NE.0) THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Maximum of ',NSHOW,' elements = ',DBLE(SHOW1)
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBLE(SHOW1))
      ELSE
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Maximum of ',0,' elements = ',ZERO
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,ZERO)
      END IF
C
      ELSEIF (SHOW.EQ.'NORM') THEN
      IF (NSHOW.EQ.0) THEN
      CALL WRNDIE(+5,'SHOW','Number of selected elements is zero.')
      CALL DECLAR('RESULT','DP',' ',DBCOMP,ZERO)
      ELSE
      DBPREC=SQRT(DBLE(SHOW1)/NSHOW)
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Norm of ',NSHOW,' elements = ',DBPREC
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBPREC)
      END IF
C
      ELSEIF (SHOW.EQ.'RMS') THEN
      IF (NSHOW.EQ.0) THEN
      CALL WRNDIE(+5,'SHOW','Number of selected elements is zero.')
      CALL DECLAR('RESULT','DP',' ',DBCOMP,ZERO)
      ELSE
      DBPREC=SQRT(DBLE(SHOW2)/NSHOW - (DBLE(SHOW1)/NSHOW)**2)
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A,F22.4)')
     & ' Rms of ',NSHOW,' elements = ',DBPREC
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBPREC)
      END IF
C
      ELSEIF (SHOW.EQ.'ELEM') THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I9,A)')
     & ' Total of ',NSHOW,' elements were selected.'
      END IF
C
      ELSE
      CALL WRNDIE(-5,'SHOW','Internal error (routine XDOSHOF)')
      END IF
      DBPREC=NSHOW
      CALL DECLAR('SELECT','DP',' ',DBCOMP,DBPREC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSHOW(SHOW,SHOW1,SHOW2,NSHOW,TYPE,
     &                   VLEVEL,VMAX,N,VSTACK,LSTACK,INDEX,
     &                   XRH,XRK,XRL)
C
C compute properties for SHOW option from stack.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      CHARACTER*(*) SHOW
      DOUBLE COMPLEX SHOW1, SHOW2
      INTEGER NSHOW
      CHARACTER*(*) TYPE
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,VMAX)
      LOGICAL LSTACK(N,VMAX)
      INTEGER INDEX(N), XRH(*), XRL(*), XRK(*)
C local
      DOUBLE COMPLEX DBCOMP, LSHOW1, LSHOW2
      DOUBLE PRECISION DBPREC, AMP, PHAS
      INTEGER I, LNSHOW
C parameters
      DOUBLE PRECISION S180
      PARAMETER (S180=180.D0)
C begin
      LSHOW1=SHOW1
      LSHOW2=SHOW2
      LNSHOW=NSHOW
      IF (SHOW.EQ.'SUM') THEN
      DO I=1,N
      LSHOW1=LSHOW1+VSTACK(I,VLEVEL)
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'AVE') THEN
      DO I=1,N
      LSHOW1=LSHOW1+VSTACK(I,VLEVEL)
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'MIN') THEN
      DO I=1,N
      IF (DBLE(LSHOW1).GT.DBLE(VSTACK(I,VLEVEL))) THEN
      LSHOW1=DBLE(VSTACK(I,VLEVEL))
      END IF
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'MAX') THEN
      DO I=1,N
      IF (DBLE(LSHOW1).LT.DBLE(VSTACK(I,VLEVEL))) THEN
      LSHOW1=DBLE(VSTACK(I,VLEVEL))
      END IF
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'NORM') THEN
      DO I=1,N
      LSHOW1=LSHOW1+ABS(VSTACK(I,VLEVEL))**2
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'RMS') THEN
      DO I=1,N
      LSHOW1=LSHOW1+DBLE(VSTACK(I,VLEVEL))
      LSHOW2=LSHOW2+(DBLE(VSTACK(I,VLEVEL)))**2
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'ELEM') THEN
      IF (TYPE.EQ.'DP') THEN
      DO I=1,N
C modification, ATB, 12/04/08, write depending on message level
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I5,A,I5,A,I5,A,F13.4)') ' [H=',XRH(INDEX(I)),
     &   ' K=',XRK(INDEX(I)),' L=',XRL(INDEX(I)),'] ',
     &    DBLE(VSTACK(I,VLEVEL))
      END IF
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBLE(VSTACK(I,VLEVEL)))
      LNSHOW=LNSHOW+1
      END DO
      ELSE
      DO I=1,N
      CALL XPHASE(VSTACK(I,VLEVEL),AMP,PHAS)
      PHAS=PHAS*S180/PI
C modification, ATB, 12/04/08, write depending on message level
      IF (WRNLEV.GE.5) THEN
      WRITE(PUNIT,'(A,I5,A,I5,A,I5,A,F13.4,A,F13.4)')
     &   ' [H=',XRH(INDEX(I)),
     &   ' K=',XRK(INDEX(I)),' L=',XRL(INDEX(I)),
     &   ']  Amplitude=',AMP,'  Phase=',PHAS
      END IF
      CALL DECLAR('RESULT','DC',' ',VSTACK(I,VLEVEL),DBPREC)
      CALL DECLAR('AMPLITUDE','DP',' ',DBCOMP,AMP)
      CALL DECLAR('PHASE','DP',' ',DBCOMP,PHAS)
      LNSHOW=LNSHOW+1
      END DO
      END IF
      ELSE
      CALL WRNDIE(-5,'SHOW','Internal error (routine XDOSHOW)')
      END IF
      SHOW1=LSHOW1
      SHOW2=LSHOW2
      NSHOW=LNSHOW
      RETURN
      END
C======================================================================
      SUBROUTINE XMDOSHO(SHOW,SHOW1,SHOW2,NSHOW,TYPE,
     &                   VLEVEL,VMAX,N,VSTACK,LSTACK,
     &                   INDEXA,INDEXB,INDEXC,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,RHOMASK)
C
C compute properties for SHOW option from stack.  Map version.
C
C Note: SUM, AVE, NORM, and RMS take special positions
C into account, e.g., the SUM over the ASU will be equal
C to the SUM over the whole unit cell divided by the
C number of symmetry operators.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      CHARACTER*(*) SHOW
      DOUBLE COMPLEX SHOW1, SHOW2
      INTEGER NSHOW
      CHARACTER*(*) TYPE
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,VMAX)
      LOGICAL LSTACK(N,VMAX)
      INTEGER INDEXA(N), INDEXB(*), INDEXC(*)
      INTEGER MAASY, MBASY, MCASY,NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      DOUBLE COMPLEX DBCOMP, LSHOW1, LSHOW2
      DOUBLE PRECISION DBPREC, AMP, PHAS
      INTEGER I, LNSHOW
C parameters
      DOUBLE PRECISION S180
      PARAMETER (S180=180.D0)
C begin
      LSHOW1=SHOW1
      LSHOW2=SHOW2
      LNSHOW=NSHOW
      IF (SHOW.EQ.'SUM') THEN
      DO I=1,N
      LSHOW1=LSHOW1+VSTACK(I,VLEVEL)
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'AVE') THEN
      DO I=1,N
      LSHOW1=LSHOW1+VSTACK(I,VLEVEL)
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'MIN') THEN
      DO I=1,N
      IF (DBLE(LSHOW1).GT.DBLE(VSTACK(I,VLEVEL))) THEN
      LSHOW1=DBLE(VSTACK(I,VLEVEL))
      END IF
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'MAX') THEN
      DO I=1,N
      IF (DBLE(LSHOW1).LT.DBLE(VSTACK(I,VLEVEL))) THEN
      LSHOW1=DBLE(VSTACK(I,VLEVEL))
      END IF
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'NORM') THEN
      DO I=1,N
      LSHOW1=LSHOW1+ABS(VSTACK(I,VLEVEL))**2
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'RMS') THEN
      DO I=1,N
      LSHOW1=LSHOW1+DBLE(VSTACK(I,VLEVEL))
      LSHOW2=LSHOW2+(DBLE(VSTACK(I,VLEVEL)))**2
      LNSHOW=LNSHOW+1
      END DO
      ELSEIF (SHOW.EQ.'ELEM') THEN
      IF (TYPE.EQ.'DP') THEN
      DO I=1,N
      WRITE(PUNIT,'(A,I5,A,I5,A,I5,A,F13.4)') ' [A=',INDEXA(I),
     &   ' B=',INDEXB(I),' C=',INDEXC(I),'] ',
     &    DBLE(VSTACK(I,VLEVEL))
      CALL DECLAR('RESULT','DP',' ',DBCOMP,DBLE(VSTACK(I,VLEVEL)))
      LNSHOW=LNSHOW+1
      END DO
      ELSE
      DO I=1,N
      CALL XPHASE(VSTACK(I,VLEVEL),AMP,PHAS)
      PHAS=PHAS*S180/PI
      WRITE(PUNIT,'(A,I5,A,I5,A,I5,A,F13.4,A,F13.4)')
     &   ' [A=',INDEXA(I),
     &   '  B=',INDEXB(I),' C=',INDEXC(I),
     &   ']  Amplitude=',AMP,'  Phase=',PHAS
      CALL DECLAR('RESULT','DC',' ',VSTACK(I,VLEVEL),DBPREC)
      CALL DECLAR('AMPLITUDE','DP',' ',DBCOMP,AMP)
      CALL DECLAR('PHASE','DP',' ',DBCOMP,PHAS)
      LNSHOW=LNSHOW+1
      END DO
      END IF
      ELSE
      CALL WRNDIE(-5,'SHOW','Internal error (routine XMDOSHO)')
      END IF
      SHOW1=LSHOW1
      SHOW2=LSHOW2
      NSHOW=LNSHOW
      RETURN
      END
