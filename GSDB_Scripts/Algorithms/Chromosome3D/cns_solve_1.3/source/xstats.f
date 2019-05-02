      SUBROUTINE XSTATS2(XRRED,XRREUP,
     &      XRSCAL,
     &      NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRNSYM,XRMSYM,XRSYTH,
     &      XRSYMM,XRITSY,XRTR,XRINTR,XRNREF,HPH,HPK,HPL,
     &      XSFNUM,XSFNAM,XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &      MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE,
     &      XRSYGP,XRSYIV,XRCELL,XRVOL)
C
C Parsing routine for STATistics facility for structure factor
C expressions.
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
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      LOGICAL XRRED, XRREUP
      DOUBLE PRECISION XRSCAL
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3),  XRINTR(3,3)
      INTEGER XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, XRMREF
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER BINMODE
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      LOGICAL ERR, QSPEC
C
      CHARACTER*2 TYPE, DOMAIN
      INTEGER DEPTH, FIRST, LAST, NSELE, NNSELE, OUNIT, NN, START, STOP
      INTEGER VLEVEL, MAXBIN, MDIM, HMAX, KMAX, LMAX, LBINS
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
C
C number of expressions
      INTEGER IEXPR, NEXPR
C
C maximum number of expressions
      INTEGER MEXPR
      PARAMETER (MEXPR=20)
C
C end point of expressions
      INTEGER EXPEND(MEXPR)
C
C depth of comand stacks
      INTEGER EXPDEP(MEXPR)
C
C selector function
      INTEGER XMODE(MEXPR)
C
C note: command stacks are stored in expression.inc
C
C definition of the functions, operands, data, and domain types
      EXTERNAL XDOFUNC, XDOOPER, XDOTYPE
C
C pointer
      INTEGER QSELE, SELE, ISHELL, NUMBIN, BINVAL
      INTEGER MSTACK, INDEX, VSTACK, LSTACK, DUMMY
      INTEGER MATRIX, COMPLT, NR
C
C parameter
      DOUBLE PRECISION M0001
      PARAMETER (M0001=0.0001)
C begin
C
C
C allocate space for reciprocal space selection arrays
      QSELE=ALLHP(ILOGIC(XRNREF))
      SELE=ALLHP(INTEG4(XRNREF))
C
C initialize number of expressions
      NEXPR=0
C
C initialize selection (default: ALL)
      CALL FILL4(HEAP(SELE),XRNREF,1)
      CALL MAKIND(HEAP(SELE),XRNREF,NSELE)
C
C initialize command stack 3
      RPNN3=0
C
C set error flag
      ERR=.FALSE.
C
C defaults
      OFILE=' '
C
C parsing
C
C capture the OVERall qualifier
      CALL NEXTWD('STATistics>')
      IF (WD(1:4).EQ.'OVER') THEN
      LBINS=1
C temporarily re-define the bins partition
      CALL XMAKEBIN(LBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE)
      ELSE
      LBINS=MBINS
      CALL SAVEWD
      END IF
C
      CALL PUSEND('STATistics>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('STATistics>')
      CALL MISCOM('STATistics>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-statistics')
C
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,HEAP(QSELE),LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
C
C convert logical selection array into list.
      CALL XSTATSE(NSELE,XRNREF,HEAP(QSELE),HEAP(SELE))
C----------------------------------------------------------------------
      ELSEIF (WD(1:1).EQ.'(') THEN
      CALL EXRPN('STATistics>',
     &     RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNTYP,
     &     RPNDOM,RPNLEV,TYPE,DOMAIN,DEPTH,XDOFUNC,XDOOPER,XDOTYPE,
     &     QHERM,ERR,XSFNUM,XSFNAM,XSFTYPE,0,' ',' ')
C
C check that no FT operation is present
      IF (.NOT.ERR) THEN
      DO NN=1,RPNN
      ERR=ERR.OR.RPN(1,NN)(1:RPNL(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('STATistics',
     &       'No FT operations allowed for STATistics.')
      END IF
      END IF
C
      CALL NEXTDO('DO>')
      IF (WD(1:1).NE.')') THEN
      CALL DSPERR('DO>',
     &  'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
C
      IF (TYPE.NE.'DP') THEN
      CALL DSPERR('DO>',
     &  'Only type "double precision" allowed.')
      ERR=.TRUE.
      END IF
C
      IF (NEXPR.GE.MEXPR) THEN
      CALL WRNDIE(-5,'XSTATS','Too many expressions specified.')
      ERR=.TRUE.
      ELSE
      NEXPR=NEXPR+1
      END IF
C
C Append all statements to command stack RPN3
C Note: this limits the length of the expressions because
C all expressions are stored in the third command stack,
C first expression, second expression, etc...
C
      IF (RPNN+RPNN3.GT.RPNMX) THEN
      CALL WRNDIE(-5,'XSTATS','Too many expressions specified.')
      ERR=.TRUE.
      ELSE
      DO NN=1,RPNN
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
C store last index of last command
      EXPEND(NEXPR)=RPNN3
C
C store depth of comand stack
      EXPDEP(NEXPR)=DEPTH
C
      XMODE(NEXPR)=0
C
      END IF
C
C----------------------------------------------------------------------
      ELSEIF (WD(1:4).EQ.'COMP') THEN
C
      IF (NEXPR.GE.MEXPR) THEN
      CALL WRNDIE(-5,'XSTATS','Too many expressions specified.')
      ERR=.TRUE.
      ELSE
      NEXPR=NEXPR+1
      END IF
C
      XMODE(NEXPR)=3
      IF (NEXPR.GT.1) THEN
      EXPEND(NEXPR)=EXPEND(NEXPR-1)
      ELSE
      EXPEND(NEXPR)=0
      END IF
C----------------------------------------------------------------------
      ELSE
      CALL CHKEND('STATistics>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      ELSE
      OUNIT=-1
      END IF
C
      IF (.NOT.ERR.AND..NOT.ERROR.AND.XRNREF.GT.0.AND.NSELE.GT.0) THEN
C
C allocate space for bins
      NUMBIN=ALLHP(INTEG4(LBINS+1))
      BINVAL=ALLHP(IREAL8((LBINS+1)*NEXPR))
      ISHELL=ALLHP(INTEG4(XRNREF))
C
C determine upper and lower bounds, fill the ISHELL array,
C compute the number of elements in each bin and initialize
C the BINVAL array.
      CALL XSTATBIN(NSELE,HEAP(SELE),XRNREF,HEAP(HPH),HEAP(HPK),
     &              HEAP(HPL),HEAP(ISHELL),XRTR,
     &              LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &              HEAP(NUMBIN),HEAP(BINVAL),NEXPR,
     &              MAXBIN)
C
C loop over all expressions
      FIRST=1
      LAST=EXPEND(1)
C
      DO IEXPR=1,NEXPR
      IF (XMODE(IEXPR).EQ.3) THEN
C
C completeness mode
C -----------------
C
C determine HMAX
      TEMPX=SQRT(XRINTR(1,1)**2+XRINTR(2,1)**2
     &      +XRINTR(3,1)**2)*XBINHIGH
      TEMPY=SQRT(XRINTR(1,2)**2+XRINTR(2,2)**2
     &      +XRINTR(3,2)**2)*XBINHIGH
      TEMPZ=SQRT(XRINTR(1,3)**2+XRINTR(2,3)**2
     &      +XRINTR(3,3)**2)*XBINHIGH
C
      IF (TEMPX.LT.M0001.OR.TEMPY.LT.M0001.OR.TEMPZ.LT.M0001) THEN
      CALL WRNDIE(-5,'PRRDAT',
     & 'high resolution limit too small or unit cell too large')
      END IF
      HMAX=INT(R4SMAL+TEMPX)+1
      KMAX=INT(R4SMAL+TEMPY)+1
      LMAX=INT(R4SMAL+TEMPZ)+1
C
C allocate space for the book-keeping matrix
      MDIM=(2*HMAX+1)*(2*KMAX+1)*(2*LMAX+1)
      MATRIX=ALLHP(INTEG4(MDIM))
C
C allocate heap space for binning
      COMPLT=ALLHP(INTEG4(LBINS))
      NR=ALLHP(INTEG4(LBINS))
C
      CALL STATCP(HEAP(QSELE),HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &            HEAP(COMPLT),HEAP(NR),HMAX,KMAX,LMAX,
     &            HEAP(MATRIX),
     &            XRNREF,LBINS,XBINHIGH,XBINLOW,BINSHELL,
     &            XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &            IEXPR,HEAP(NUMBIN),HEAP(BINVAL),XRTR)
C
C free up heap space for binning
      CALL FREHP(NR,INTEG4(LBINS))
      CALL FREHP(COMPLT,INTEG4(LBINS))
C free space for book-keeping matrix
      CALL FREHP(MATRIX,INTEG4(MDIM))
C------------------------------------------------------------------------
      ELSE
C
C normal statistics mode
C ----------------------
C
      IF (IEXPR.GT.1) THEN
      FIRST=EXPEND(IEXPR-1)+1
      END IF
      LAST=EXPEND(IEXPR)
C
      RPNN=0
      DO NN=FIRST,LAST
      RPNN=RPNN+1
      RPN(1,RPNN)=RPN3(1,NN)
      RPN(2,RPNN)=RPN3(2,NN)
      RPN(3,RPNN)=RPN3(3,NN)
      RPN(4,RPNN)=RPN3(4,NN)
      RPNL(1,RPNN)=RPNL3(1,NN)
      RPNL(2,RPNN)=RPNL3(2,NN)
      RPNL(3,RPNN)=RPNL3(3,NN)
      RPNL(4,RPNN)=RPNL3(4,NN)
      RPNDB(1,RPNN)=RPNDB3(1,NN)
      RPNDB(2,RPNN)=RPNDB3(2,NN)
      RPNDB(3,RPNN)=RPNDB3(3,NN)
      RPNDB(4,RPNN)=RPNDB3(4,NN)
      RPNMLT(RPNN)=RPNMLT3(NN)
      RPNTYP(RPNN)=RPNTYP3(NN)
      RPNDOM(RPNN)=RPNDOM3(NN)
      RPNLEV(RPNN)=RPNLEV3(NN)
      END DO
C
      DEPTH=EXPDEP(IEXPR)
C
C check if special structure factor operations are present
      CALL XDOSPCL(RPNMX,RPNX,RPNN,RPN,RPNL,QSPEC)
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
      STOP=MIN(NSELE,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP
      CALL XDOMAKS(HEAP(INDEX),START,STOP,HEAP(SELE),NNSELE)
C
C evaluate expression
      CALL XDOEVAL(RPNMX,RPNN,RPNX,RPN,RPNL,RPNDB,RPNMLT,RPNLEV,
     & RPNTYP,RPNDOM,VLEVEL,DEPTH,NNSELE,HEAP(VSTACK),
     & HEAP(LSTACK),HEAP(INDEX),XRTR,
     & HPH,HPK,HPL,
     & XSFNUM,XSFNAM,XSFTYPE,HPSF,
     & HPMULT,HPTYPE,
     & DUMMY,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,LBINS,XBINLOW,XBINHIGH,BINSHELL,
     & XRCELL,XRVOL)
C
C copy result in BINVAL array
      CALL XDOCOPS(VLEVEL,DEPTH,HEAP(VSTACK),HEAP(INDEX),NNSELE,
     &             HEAP(ISHELL),LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &             IEXPR,HEAP(BINVAL))
C
      START=START+MSTACK
      STOP=MIN(NSELE,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
C
      END IF
C--------------------------------------------------------------
C end loop over all expressions
      END DO
C
C print bin values if OUNIT is specified.
      IF (OUNIT.GE.1) THEN
      CALL XSTATPR(OUNIT,LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &             NEXPR,HEAP(NUMBIN),
     &             HEAP(BINVAL))
      END IF
C
C declare variables
      CALL XSTATVA(LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &             NEXPR,HEAP(NUMBIN),
     &             HEAP(BINVAL))
C
C de-allocate space for bins
      CALL FREHP(NUMBIN,INTEG4(LBINS+1))
      CALL FREHP(BINVAL,IREAL8((LBINS+1)*NEXPR))
      CALL FREHP(ISHELL,INTEG4(XRNREF))
C
      END IF
C
C de-allocate space for reciprocal space selection arrays
      CALL FREHP(QSELE,ILOGIC(XRNREF))
      CALL FREHP(SELE,INTEG4(XRNREF))
C
      IF (LBINS.NE.MBINS) THEN
C re-create the original bin partition
      CALL XMAKEBIN(MBINS,XBINLOW,XBINHIGH,BINSHELL,BINMODE)
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XDOMAKS(INDEX,START,STOP,SELE,NSELE)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER INDEX(*), START, STOP, SELE(*), NSELE
C local
      INTEGER I
C begin
      NSELE=0
      DO I=START,STOP
      NSELE=NSELE+1
      INDEX(NSELE)=SELE(I)
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOCOPS(VLEVEL,VMAX,VSTACK,INDEX,
     &             NNSELE,ISHELL,LBINS,XBINLOW,
     &             XBINHIGH,BINSHELL,IEXPR,BINVAL)
C
C routine copies stack element into appropriate bin
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER VLEVEL, VMAX, NNSELE
      DOUBLE COMPLEX VSTACK(NNSELE,VMAX)
      INTEGER INDEX(*), ISHELL(*)
      INTEGER LBINS, IEXPR
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION BINVAL(0:LBINS,*)
C local
      INTEGER I
C begin
C
      DO I=1,NNSELE
      BINVAL(ISHELL(INDEX(I)),IEXPR)=BINVAL(ISHELL(INDEX(I)),IEXPR)+
     &                               DBLE(VSTACK(I,VLEVEL))
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSTATBIN(NSELE,SELE,XRNREF,XRH,XRK,XRL,ISHELL,XRTR,
     &                    LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                    NUMBIN,BINVAL,NEXPR,MAXBIN)
C
C Routine determines the partitioning into equal-volume
C bins in reciprocal space.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER NSELE, SELE(*), XRNREF, XRH(*), XRK(*), XRL(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION XRTR(3,*)
      INTEGER LBINS, NUMBIN(0:LBINS)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION BINVAL(0:LBINS,*)
      INTEGER NEXPR, MAXBIN
C local
      LOGICAL COND, OUTSIDE
      INTEGER REFLCT, IEXPR, H, K, L, IND, BIN
      DOUBLE PRECISION SSQ
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0)
C begin
C
      IF (XBINLOW.EQ.ZERO.AND.XBINHIGH.EQ.ZERO) THEN
      CALL WRNDIE(-1,'XSTATBIN',
     & ' Fatal error -- BINResolution not specified. ')
      ELSE
      OUTSIDE=.FALSE.
C
      IF (NSELE.GT.0) THEN
C
C
C compute shell bin for each reflection
      DO REFLCT=1,XRNREF
      ISHELL(REFLCT)=0
      END DO
C
      DO REFLCT=1,NSELE
      H=XRH(SELE(REFLCT))
      K=XRK(SELE(REFLCT))
      L=XRL(SELE(REFLCT))
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      BIN=1
      COND=.FALSE.
      DO WHILE (.NOT.COND.AND.BIN.LE.LBINS)
      BIN=BIN+1
      COND= (SSQ.GE.BINSHELL(BIN).AND.SSQ.LE.BINSHELL(BIN-1))
      END DO
C
      IF (COND) THEN
      ISHELL(SELE(REFLCT))=BIN-1
      ELSE
      ISHELL(SELE(REFLCT))=0
      OUTSIDE=.TRUE.
      END IF
      END DO
C
      END IF
C
C initialize bins
      DO IND=0,LBINS
      NUMBIN(IND)=0
      END DO
      DO IEXPR=1,NEXPR
      DO IND=0,LBINS
      BINVAL(IND,IEXPR)=ZERO
      END DO
      END DO
C
      DO REFLCT=1,XRNREF
      NUMBIN(ISHELL(REFLCT))=NUMBIN(ISHELL(REFLCT))+1
      END DO
C
C
C determine maximum number of elements per bin
      MAXBIN=0
      DO IND=1,LBINS
      MAXBIN=MAX(MAXBIN,NUMBIN(IND))
      END DO
C
      IF (OUTSIDE) THEN
      WRITE(6,'(2A)' )
     & ' XSTATBIN-info: BINResolution range smaller than selection! ',
     & ' Elements outside range will be 0. '
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSTATSE(NSELE,XRNREF,QSELE,SELE)
C
C Routine converts logical selection array into list.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NSELE, XRNREF
      LOGICAL QSELE(*)
      INTEGER SELE(*)
C local
      INTEGER I
C begin
      NSELE=0
      DO I=1,XRNREF
      IF (QSELE(I)) THEN
      NSELE=NSELE+1
      SELE(NSELE)=I
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XSTATPR(OUNIT,LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                   NEXPR,NUMBIN,BINVAL)
C
C Routine prints the average values of all expressions
C to specified output file.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER OUNIT, LBINS, NEXPR
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER NUMBIN(0:*)
      DOUBLE PRECISION BINVAL(0:LBINS,*)
C local
      INTEGER I, IND
C parameter
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C begin
      WRITE(OUNIT,'(A)')
     &   ' #bin | resolution range | #refl | '
      DO IND=LBINS,1,-1
      IF (ABS(BINSHELL(IND+1)).GT.RSMALL) THEN
      WRITE(OUNIT,'(A,I3,A,F6.2,A,F6.2,A,3X,I6,A,20F10.4)')
     & '  ',LBINS-IND+1,' ',ONE/SQRT(BINSHELL(IND)),'  ',
     & ONE/SQRT(BINSHELL(IND+1)),
     & '  ',NUMBIN(IND),
     & '  ',(BINVAL(IND,I)/MAX(1,NUMBIN(IND)),I=1,NEXPR)
      ELSE
      WRITE(OUNIT,'(A,I3,A,F6.2,A,     A,3X,I6,A,20F10.4)')
     & '  ',LBINS-IND+1,' ',ONE/SQRT(BINSHELL(IND)),
     & ' infinity',
     & ' ',NUMBIN(IND),
     & '  ',(BINVAL(IND,I)/MAX(1,NUMBIN(IND)),I=1,NEXPR)
      END IF
C
      END DO
      RETURN
      END
C=====================================================================
      SUBROUTINE XSTATTAR(IND,ISHELL,TSEL,XRNREF,NN)
C
C
C reduces data to a resolution bin.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER IND, ISHELL(*)
      INTEGER TSEL(*), XRNREF, NN
C local
      INTEGER REFLCT
C begin
C
C select all reflections in that bin
      NN=0
      DO REFLCT=1,XRNREF
      IF (ISHELL(REFLCT).EQ.IND) THEN
      NN=NN+1
      TSEL(REFLCT)=1
      ELSE
      TSEL(REFLCT)=0
      END IF
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE XSTATVA(LBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                   NEXPR,NUMBIN,BINVAL)
C
C Declares variables $expression1, $expression2, ...
C for highest resolution bin.  $reflections contains the
C number of reflections in that bin.
C
C Also declaures the symbol structure:
C    $statistics.nbins
C    $statistics.<i>.nref
C    $statistics.<i>,high
C    $statistics.<i>.low
C    $statistics.<i>.value
C for the first "expression".
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
      INTEGER LBINS, NEXPR
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER NUMBIN(0:*)
      DOUBLE PRECISION BINVAL(0:LBINS,*)
C local
      INTEGER STMAX
      PARAMETER (STMAX=60)
      INTEGER I, IND
      CHARACTER*(STMAX) STRING,ADST
      INTEGER STLEN,ADLEN
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C begin
      DO I=1,NEXPR
      CALL ENCODI(I,ADST,STMAX,ADLEN)
      STRING='EXPRESSION'//ADST(1:ADLEN)
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      DBPREC=(BINVAL(1,I)/MAX(1,NUMBIN(1)))
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
      END DO
      DBPREC=NUMBIN(1)
      CALL DECLAR('REFLECTIONS','DP',' ',DBCOMP,DBPREC)
C
C modification 9/08/09 ATB
      STRING='STATISTICS.NBINS'
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      DBPREC=LBINS
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
C
      DO IND=LBINS,1,-1
C
      CALL ENCODI(IND,ADST,STMAX,ADLEN)
      STRING='STATISTICS.'//ADST(1:ADLEN)//'.HIGH'
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      DBPREC=ONE/SQRT(BINSHELL(IND))
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
C
      STRING='STATISTICS.'//ADST(1:ADLEN)//'.LOW'
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      IF (ABS(BINSHELL(IND+1)).GT.RSMALL) THEN      
      DBPREC=ONE/SQRT(BINSHELL(IND+1))
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
      ELSE
      DBPREC=9999.9D0
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
      END IF
C
      STRING='STATISTICS.'//ADST(1:ADLEN)//'.NREF'
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      DBPREC=NUMBIN(IND)
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
C  
      STRING='STATISTICS.'//ADST(1:ADLEN)//'.VALUE'
      STLEN=STMAX
      CALL TRIMM(STRING,STLEN)
      DBPREC=(BINVAL(IND,1)/MAX(1,NUMBIN(IND)))
      CALL DECLAR(STRING(1:STLEN),'DP',' ',DBCOMP,DBPREC)
C
      END DO
C end modification
C
      RETURN
      END
C===============================================================
      SUBROUTINE STATCP(QSELE,XRH,XRK,XRL,COMPLT,NR,HMAX,KMAX,LMAX,
     &                   MATRIX,
     &                   XRNREF,LBINS,XBINHIGH,XBINLOW,BINSHELL,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   IEXPR,NUMBIN,BINVAL,XRTR)
C
C Computes completeness of selected data in resolution bins.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'comand.inc'
      LOGICAL QSELE(*)
      INTEGER XRH(*), XRK(*), XRL(*), NR(*)
      INTEGER COMPLT(*)
      INTEGER HMAX, KMAX, LMAX
      INTEGER MATRIX(-HMAX:HMAX,-KMAX:KMAX,-LMAX:LMAX)
      INTEGER XRNREF
      INTEGER LBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER IEXPR
      INTEGER NUMBIN(0:*)
      DOUBLE PRECISION BINVAL(0:LBINS,*), XRTR(3,3)
C local
      INTEGER H, K, L, REFLCT, IISYM, IFRIED, HH, KK, LL
      INTEGER IND, I
      DOUBLE PRECISION SSQ, CC, S
      LOGICAL QERR, SYSAB
C parameters
      DOUBLE PRECISION ZERO, ONE, FOUR, THREE, S100
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, FOUR=4.0D0, THREE=3.0D0)
      PARAMETER (S100=100.D0)
C begin
C
C
C initialize the book-keeping matrix
      DO L=-LMAX,LMAX
      DO K=-KMAX,KMAX
      DO H=-HMAX,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
C
C initialize the bin counters
      DO IND=1,LBINS
      NR(IND)=0
      COMPLT(IND)=0
      END DO
C
C fill the present reflections
      QERR=.FALSE.
      DO REFLCT=1,XRNREF
      IF (QSELE(REFLCT)) THEN
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
C
C check if this is a systematic absence
      CALL XRSYSAB(H,K,L,SYSAB,
     &                  XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
      IF (.NOT.SYSAB) THEN
      MATRIX(H,K,L)=REFLCT
      CALL XRSSQ(XRH(REFLCT),XRK(REFLCT),XRL(REFLCT),SSQ,XRTR)
      QERR=QERR.OR.((H.LT.-HMAX.OR.H.GT.HMAX.OR.
     &              K.LT.-KMAX.OR.K.GT.KMAX.OR.
     &              L.LT.-LMAX.OR.L.GT.LMAX)).AND.
     &              SSQ.LE.BINSHELL(1).AND.SSQ.GE.BINSHELL(LBINS+1)
      IND=0
      DO I=1,LBINS
      IF (SSQ.LE.BINSHELL(I).AND.SSQ.GE.BINSHELL(I+1).AND.IND.EQ.0)
     &  IND=I
      END DO
      IF (IND.GT.0) THEN
      NR(IND)=NR(IND)+1
      END IF
      END IF
      END IF
      END DO
C
      IF (QERR) THEN
      CALL WRNDIE(-5,'STATCP','fatal coding error 1')
      END IF
C
C loop through the box determined by h,k,l max.
      DO H=-HMAX,HMAX
      DO K=-KMAX,KMAX
      DO L=-LMAX,LMAX
C
C compute s**2 for this reflection
      CALL XRSSQ(H,K,L,SSQ,XRTR)
C
C check resolution limits
      S=SQRT(SSQ)
      IF (S.LT.XBINHIGH.AND.S.GT.XBINLOW) THEN
C
C map this reflection into the asymmetric unit (returned in HH, KK, LL)
      CALL XRASYM(H,K,L,HH,KK,LL,IISYM,IFRIED,
     &                  XRNSYM,XRMSYM,XRITSY,QHERM)
C
C check if this is a systematic absence
      CALL XRSYSAB(HH,KK,LL,SYSAB,
     &                  XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C check whether this reflection not yet present
      IF (MATRIX(HH,KK,LL).EQ.0.AND..NOT.SYSAB) THEN
      MATRIX(HH,KK,LL)=1
C
C check resolution limits
      IND=0
      DO I=1,LBINS
      IF (SSQ.LE.BINSHELL(I).AND.SSQ.GE.BINSHELL(I+1)) IND=I
      END DO
      IF (IND.GT.0) THEN
      COMPLT(IND)=COMPLT(IND)+1
      END IF
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
C
      DO IND=LBINS,1,-1
      IF (COMPLT(IND)+NR(IND).GT.RSMALL) THEN
      CC=NR(IND)
      CC=CC/(COMPLT(IND)+NR(IND))
      ELSE
      CC=ZERO
      END IF
C multiply the completeness by number of elements in bin (NUMBIN)
C because the print routine will divide by this number.
      BINVAL(IND,IEXPR)=NUMBIN(IND)*CC
      END DO
C
      RETURN
      END
C======================================================================
