      SUBROUTINE XTARGETS(QDERIV,QPRINT,ITEST,XRE,XDERIV,
     &           HPTSEL,XRNREF,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &           DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &           MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &           MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &           HPH,HPK,HPL,XSFNUM,XSFNAM,
     &           XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,
     &           QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,IHPFLAG,
     &           XRCELL,XRVOL)
C
C Computes target functions for crystallographic
C refinement and crystallographic searches.
C
C If requested (QDERIV) routine computes the derivatives of the
C target with respect to FCALC.
C Routine returns the derivatives in FCALC.
C
C If requested, a monitor function is also computed.
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      LOGICAL QDERIV, QPRINT
      INTEGER ITEST
      DOUBLE PRECISION XRE
      INTEGER XRNREF
      DOUBLE COMPLEX XDERIV(XRNREF,*)
      INTEGER HPTSEL
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XRSCAL
      LOGICAL XCVTEST
C
      INTEGER TRPNMX, TRPNX
      CHARACTER*(*) TRPN(4,TRPNMX)
      INTEGER TRPNL(4,TRPNMX), TRPNN
      DOUBLE COMPLEX TRPNDB(4,TRPNMX)
      INTEGER TRPNMLT(TRPNMX), TRPNLEV(TRPNMX)
      CHARACTER*2 TRPNTYP(TRPNMX), TRPNDOM(TRPNMX)
      INTEGER TDEPTH
C
      INTEGER DRPNMX, DRPNX
      CHARACTER*(*) DRPN(4,DRPNMX,*)
      INTEGER DRPNL(4,DRPNMX,*), DRPNN(*)
      DOUBLE COMPLEX DRPNDB(4,DRPNMX,*)
      INTEGER DRPNMLT(DRPNMX,*), DRPNLEV(DRPNMX,*)
      CHARACTER*2 DRPNTYP(DRPNMX,*), DRPNDOM(DRPNMX,*)
      INTEGER DDEPTH(*)
C
      INTEGER MRPNMX, MRPNX
      CHARACTER*(*) MRPN(4,MRPNMX)
      INTEGER MRPNL(4,MRPNMX), MRPNN
      DOUBLE COMPLEX MRPNDB(4,MRPNMX)
      INTEGER MRPNMLT(MRPNMX), MRPNLEV(MRPNMX)
      CHARACTER*2 MRPNTYP(MRPNMX), MRPNDOM(MRPNMX)
      INTEGER MDEPTH
C
      INTEGER HPH, HPK, HPL, XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*), HPMULT, HPTYPE
      INTEGER XRMREF
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER IHPFLAG
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C local
      DOUBLE PRECISION MONITOR
      INTEGER COUNT
C pointer
      INTEGER XTARGT
C parameter
      DOUBLE PRECISION ZERO, ONE
      DOUBLE COMPLEX CPLXZERO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
      CPLXZERO=DCMPLX(ZERO,ZERO)
C
C
C compute general target function and derivatives (if QDERIV)
C for each ASSOciate structure factor object
C =============================================================
      XTARGT=ALLHP(ICPLX8(XRNREF))
C
      CALL XTAREXPR(QDERIV,ITEST,XRE,XDERIV,HEAP(XTARGT),
     &              TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &              TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &              DRPNMX,DRPNX,DRPN,DRPNL,
     &              DRPNN,DRPNDB,DRPNMLT,
     &              DRPNLEV,DRPNTYP,DRPNDOM,
     &              DDEPTH,
     &              XRTR,HPH,HPK,HPL,XSFNUM,XSFNAM,
     &              XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,HEAP(HPTSEL),
     &              XRNREF,QHERM,
     &              XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &              MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &              XRSCAL,IHPFLAG,XRCELL,XRVOL)
C
C compute the monitor function if required
      IF (MRPNN.GT.0.AND.QPRINT) THEN
      CALL XTAREXPR(.FALSE.,ITEST,MONITOR,XDERIV,HEAP(XTARGT),
     &              MRPNMX,MRPNX,MRPN,MRPNL,MRPNN,MRPNDB,
     &              MRPNMLT,MRPNLEV,MRPNTYP,MRPNDOM,MDEPTH,
     &              1,1,' ',1,0,CPLXZERO,1,1,' ',' ',1,
     &              XRTR,HPH,HPK,HPL,XSFNUM,XSFNAM,
     &              XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,HEAP(HPTSEL),
     &              XRNREF,QHERM,
     &              XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &              MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &              ONE,0,XRCELL,XRVOL)
      CALL XMONCOU(XRNREF,HEAP(HPTSEL),ITEST,COUNT)
      IF (COUNT.GT.0) THEN
CCC      IF (WRNLEV.GE.5) THEN
      IF (XCVTEST.AND.ITEST.EQ.1) THEN
      WRITE(6,'(A,F7.3)')
     &   ' XTAREXPR: ->[WORKING SET] monitor=',
     &    MONITOR/COUNT
      ELSEIF (XCVTEST.AND.ITEST.EQ.-1) THEN
      WRITE(6,'(A,F7.3)')
     &   ' XTAREXPR: ->[TEST SET]    monitor=',
     &    MONITOR/COUNT
      ELSEIF (.NOT.XCVTEST) THEN
      WRITE(6,'(A,F7.3)')
     &   ' XTAREXPR:                 monitor=',
     &    MONITOR/COUNT
      END IF
CCC      END IF
C
      IF (XCVTEST.AND.ITEST.EQ.1) THEN
      CALL DECLAR('MONITOR','DP',' ',CPLXZERO,MONITOR/COUNT)
      ELSEIF (XCVTEST.AND.ITEST.EQ.-1) THEN
      CALL DECLAR('TEST_MONITOR','DP',' ',CPLXZERO,MONITOR/COUNT)
      ELSEIF (.NOT.XCVTEST) THEN
      CALL DECLAR('MONITOR','DP',' ',CPLXZERO,MONITOR/COUNT)
      END IF
C
      END IF
C
      END IF
C
      CALL FREHP(XTARGT,ICPLX8(XRNREF))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XTAREXPR(QDERIV,ITEST,XRE,XDERIV,XTARGT,
     &                    TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &                    TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &                    DRPNMX,DRPNX,DRPN,DRPNL,DRPNN,DRPNDB,
     &                    DRPNMLT,DRPNLEV,DRPNTYP,DRPNDOM,DDEPTH,
     &                    XRTR,HPH,HPK,HPL,XSFNUM,XSFNAM,
     &                    XSFTYPE,HPSF,HPMULT,HPTYPE,XRMREF,TSEL,
     &                    XRNREF,QHERM,
     &                    XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                    MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                    XRSCAL,IHPFLAG,XRCELL,XRVOL)
C
C
C Target (target), derivative (dtarget), and monitor (monitor)
C evaluation routine.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      LOGICAL QDERIV
      INTEGER ITEST
      DOUBLE PRECISION XRE
      INTEGER XRNREF
      DOUBLE COMPLEX XDERIV(XRNREF,*)
      DOUBLE COMPLEX XTARGT(*)
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
      DOUBLE PRECISION XRTR(3,3)
      INTEGER HPH, HPK, HPL, XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*), HPMULT, HPTYPE
      INTEGER XRMREF, TSEL(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRSCAL
      INTEGER IHPFLAG
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C local
      INTEGER NSELE, NNSELE, MSTACK, START, STOP, DEPTH, VLEVEL
      INTEGER REFLCT, NN, II
      LOGICAL QSPEC, QSPEC2, ERR
C pointer
      INTEGER INDEX, VSTACK, LSTACK
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C initialize energy
      XRE=ZERO
C
C initialize derivatives
      IF (QDERIV) THEN
      DO II=1,IHPFLAG
      DO REFLCT=1,XRNREF
      XDERIV(REFLCT,II)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END IF
C
C initialize array that contains the result of the target expression
C calculation for each reflection
      DO REFLCT=1,XRNREF
      XTARGT(REFLCT)=DCMPLX(ZERO,ZERO)
      END DO
C
      NNSELE=0
C check if special operations are present
      CALL XDOSPCL(TRPNMX,TRPNX,TRPNN,TRPN,TRPNL,QSPEC)
      DEPTH=TDEPTH
      DO II=1,IHPFLAG
      CALL XDOSPCL(DRPNMX,DRPNX,DRPNN(II),DRPN(1,1,II),
     &             DRPNL(1,1,II),QSPEC2)
      DEPTH=MAX(DEPTH,DDEPTH(II))
      QSPEC=QSPEC.OR.QSPEC2
      END DO
C
C check that no FT operation is present
      ERR=.FALSE.
      DO NN=1,TRPNN
      ERR=ERR.OR.TRPN(1,NN)(1:TRPNL(1,NN)).EQ.'FT'
      END DO
      DO II=1,IHPFLAG
      DO NN=1,DRPNN(II)
      ERR=ERR.OR.DRPN(1,NN,II)(1:DRPNL(1,NN,II)).EQ.'FT'
      END DO
      END DO
      IF (ERR) THEN
      CALL DSPERR('XTAREXPR','No FT operations allowed for target.')
      ELSEIF (TRPNN.GT.0) THEN
C
C For special operations we have to allocate space for
C a stack that is big enough to hold all structure factor
C elements.  For all other operations we can do them
C in smaller junks.
C
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
C compute target function
C =======================
Cbegin loop over junks  (size MSTACK, from START to STOP)
C---------------------
      START=1
      STOP=MIN(XRNREF,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP in specified set
      CALL XDOMAKT(HEAP(INDEX),START,STOP,NSELE,
     &            ITEST,TSEL)
C
      NNSELE=NNSELE+NSELE
C
C evaluate expression
      CALL XDOEVAL(TRPNMX,TRPNN,TRPNX,TRPN,TRPNL,TRPNDB,TRPNMLT,
     &     TRPNLEV,TRPNTYP,TRPNDOM,
     &     VLEVEL,TDEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     0,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &     XRCELL,XRVOL)
C
C assign to target function array
      CALL XDOADC(VLEVEL,TDEPTH,HEAP(VSTACK),NSELE,XTARGT,HEAP(INDEX))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
C evaluate the sum of the target function array
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).EQ.ITEST) THEN
      XRE=XRE + XRSCAL*(DBLE(XTARGT(REFLCT))+DIMAG(XTARGT(REFLCT)))
      END IF
      END DO
C
      IF (QDERIV) THEN
C compute derivative of target function with respect to each ASSOciate
C object
C ===========================================================
C
Cbegin loop over associate objects
C---------------------------------
      DO II=1,IHPFLAG
      IF (DRPNN(II).GT.0) THEN
Cbegin loop over junks  (size MSTACK, from START to STOP)
C---------------------
      NNSELE=0
      START=1
      STOP=MIN(XRNREF,MSTACK)
      DO WHILE (START.LE.STOP)
C
C select all elements between START and STOP in specified set
      CALL XDOMAKT(HEAP(INDEX),START,STOP,NSELE,
     &             ITEST,TSEL)
C
      NNSELE=NNSELE+NSELE
C
C evaluate expression
      CALL XDOEVAL(DRPNMX,DRPNN(II),DRPNX,DRPN(1,1,II),
     &     DRPNL(1,1,II),DRPNDB(1,1,II),DRPNMLT(1,II),
     &     DRPNLEV(1,II),DRPNTYP(1,II),DRPNDOM(1,II),
     &     VLEVEL,DDEPTH(II),NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     0,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &     XRCELL,XRVOL)
C
C multiply result by XRSCAL
      CALL XDOMULL(VLEVEL,DDEPTH,HEAP(VSTACK),NSELE,XRSCAL,ZERO)
C copy into XDERIV(*,II) array
      CALL XDOADC(VLEVEL,DDEPTH,HEAP(VSTACK),NSELE,XDERIV(1,II),
     &            HEAP(INDEX))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
      ELSE
      WRITE(6,'(A)')
     & ' Some dtarget expressions for associate objects undefined.'
      CALL WRNDIE(-5,'XTAREXPR',
     &           'dtarget expression undefined')
C
      END IF
      END DO
Cend loop over associate objects
C-------------------------------
C
      END IF
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
C
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XDOMAKT(INDEX,START,STOP,NSELE,ITEST,TSEL)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER INDEX(*), START, STOP, NSELE, ITEST, TSEL(*)
C local
      INTEGER I
C begin
      NSELE=0
      DO I=START,STOP
      IF (TSEL(I).EQ.ITEST) THEN
      NSELE=NSELE+1
      INDEX(NSELE)=I
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOBINT(N,XRTR,XRH,XRK,XRL,MBINS,
     &                   XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C Creates array that assigns each reflection to a particular
C resolution-dependent bin.  Reflections that are outside
C the bin range will be given an index of 0.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER N
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER ISHELL(*)
C local
      INTEGER REFLCT
      LOGICAL COND, OUTSIDE
      INTEGER H, K, L, BIN
      DOUBLE PRECISION SSQ
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, THREE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
C begin
C
      IF (XBINLOW.EQ.ZERO.AND.XBINHIGH.EQ.ZERO) THEN
      CALL WRNDIE(-1,'XDOBINT',
     & ' Fatal error -- BINResolution not specified. ')
      ELSE
      OUTSIDE=.FALSE.
C
C compute shell bin for each reflection
      DO REFLCT=1,N
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      BIN=1
      COND=.FALSE.
      DO WHILE (.NOT.COND.AND.BIN.LE.MBINS)
      BIN=BIN+1
      COND= (SSQ.GE.BINSHELL(BIN).AND.SSQ.LE.BINSHELL(BIN-1))
      END DO
C
      IF (COND) THEN
      ISHELL(REFLCT)=BIN-1
      ELSE
      OUTSIDE=.TRUE.
      ISHELL(REFLCT)=0
      END IF
      END DO
C
      IF (OUTSIDE) THEN
      WRITE(6,'(2A)' )
     & ' XDOBINT-info: BINResolution range smaller than selection! ',
     & ' Elements outside range will be 0. '
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XTSELSET(HPH,HPK,HPL,
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
C Routine selects reflections based on target selection expression.
C
C This routine will also set the flag XCVTEST which indicates
C that cross-validation is active.
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      LOGICAL XRREUP
      INTEGER HPTSEL, XRMREF, XRNREF
      LOGICAL XCVTEST
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
C
      INTEGER SRPNMX, SRPNX
      CHARACTER*(*) SRPN(4,SRPNMX)
      INTEGER SRPNL(4,SRPNMX), SRPNN
      DOUBLE COMPLEX SRPNDB(4,SRPNMX)
      INTEGER SRPNMLT(SRPNMX), SRPNLEV(SRPNMX)
      CHARACTER*2 SRPNTYP(SRPNMX), SRPNDOM(SRPNMX)
      INTEGER SDEPTH
C
      INTEGER CRPNMX, CRPNX
      CHARACTER*(*) CRPN(4,CRPNMX)
      INTEGER CRPNL(4,CRPNMX), CRPNN
      DOUBLE COMPLEX CRPNDB(4,CRPNMX)
      INTEGER CRPNMLT(CRPNMX), CRPNLEV(CRPNMX)
      CHARACTER*2 CRPNTYP(CRPNMX), CRPNDOM(CRPNMX)
      INTEGER CDEPTH
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C
C local
      INTEGER XRIREF
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C begin
C
C only do the update if it is required
      IF (XRREUP) THEN
      XRREUP=.FALSE.
C
      CALL XRTES2(HPH,HPK,HPL,
     &            XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &            HPMULT,
     &            HPTYPE,
     &            XRMREF,XRNREF,XRIREF,XCVTEST,HEAP(HPTSEL),
     &            XRTR,
     &            SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &            SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &            CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &            CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &            QHERM,XRNSYM,XRMSYM,XRSYTH,
     &            XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &            XRCELL,XRVOL)
C
C
C print out
      IF (XRIREF.EQ.0) THEN
      WRITE(6,'(A)')
     & ' %XTSELSET-ERR: number of selected reflections is zero.'
      ELSE
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I7)')
     & ' XRTSELSET: number of selected reflections ',XRIREF
      END IF
      END IF
C
C define symbol
      DBPREC=XRIREF
      CALL DECLAR( 'TSELECT','DP',' ',DBCOMP,DBPREC)
C
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XRTES2(HPH,HPK,HPL,
     &                  XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &                  HPMULT,HPTYPE,XRMREF,XRNREF,XRIREF,XCVTEST,
     &                  TSEL,
     &                  XRTR,
     &                  SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &                  SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &                  CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &                  CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &                  QHERM,XRNSYM,XRMSYM,XRSYTH,
     &                  XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &                  XRCELL,XRVOL)
C
C Routine selects reflections based on resolution and F cutoff
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER XRMREF,XRNREF, XRIREF
      LOGICAL XCVTEST
      INTEGER TSEL(*)
      DOUBLE PRECISION XRTR(3,3)
C
      INTEGER SRPNMX, SRPNX
      CHARACTER*(*) SRPN(4,SRPNMX)
      INTEGER SRPNL(4,SRPNMX), SRPNN
      DOUBLE COMPLEX SRPNDB(4,SRPNMX)
      INTEGER SRPNMLT(SRPNMX), SRPNLEV(SRPNMX)
      CHARACTER*2 SRPNTYP(SRPNMX), SRPNDOM(SRPNMX)
      INTEGER SDEPTH
C
      INTEGER CRPNMX, CRPNX
      CHARACTER*(*) CRPN(4,CRPNMX)
      INTEGER CRPNL(4,CRPNMX), CRPNN
      DOUBLE COMPLEX CRPNDB(4,CRPNMX)
      INTEGER CRPNMLT(CRPNMX), CRPNLEV(CRPNMX)
      CHARACTER*2 CRPNTYP(CRPNMX), CRPNDOM(CRPNMX)
      INTEGER CDEPTH
C
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(3,3), XRVOL
C local
      INTEGER REFLCT, MSTACK
      LOGICAL ERR, QSPEC
C
      INTEGER VLEVEL, NN, START, STOP, NSELE
C pointer
      INTEGER INDEX, VSTACK, LSTACK, HPSTORE
C begin
C
      IF (SRPNN.GT.0) THEN
C
C a target selection expression is defined
C
C check that no FT operation is present
      ERR=.FALSE.
      DO NN=1,SRPNN
      ERR=ERR.OR.SRPN(1,NN)(1:SRPNL(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('TSELection',
     &  'No FT allowed in selections.')
      END IF
C
C check if special structure factor operations are present
      CALL XDOSPCL(SRPNMX,SRPNX,SRPNN,SRPN,SRPNL,QSPEC)
C
C
C initialize QSELE array
      DO REFLCT=1,XRNREF
      TSEL(REFLCT)=0
      END DO
C
C do it in MSTACK junks except when special structure factor
C operations are present.
      IF (QSPEC) THEN
      MSTACK=XRNREF
      ELSE
      MSTACK=MIN(XRNREF,10000)
      END IF
C
      INDEX=ALLHP(INTEG4(MSTACK))
      VSTACK=ALLHP(ICPLX8(MSTACK*SDEPTH))
      LSTACK=ALLHP(ILOGIC(MSTACK*SDEPTH))
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
      HPSTORE=0
      CALL XDOEVAL(SRPNMX,SRPNN,SRPNX,SRPN,SRPNL,SRPNDB,
     &     SRPNMLT,SRPNLEV,
     &     SRPNTYP,SRPNDOM,VLEVEL,SDEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &     XRCELL,XRVOL)
C
C copy selection into TSEL array
      CALL XDOIN(VLEVEL,SDEPTH,HEAP(LSTACK),NSELE,TSEL,HEAP(INDEX))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*SDEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*SDEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
C
      ELSE
C
C no TSELection specified: select all reflections
C
      DO REFLCT=1,XRNREF
      TSEL(REFLCT)=1
      END DO
C
      WRITE(6,'(A)')
     & ' **TSELSET-info**: no TSELection specified: using all refl.'
C
      END IF
C
C define cross-validation if defined
      IF (CRPNN.GT.0) THEN
C
C a cross-validation selection expression is defined
C
C check that no FT operation is present
      ERR=.FALSE.
      DO NN=1,CRPNN
      ERR=ERR.OR.CRPN(1,NN)(1:CRPNL(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('CVSElection',
     &  'No FT allowed in selections.')
      END IF
C
C check if special structure factor operations are present
      CALL XDOSPCL(CRPNMX,CRPNX,CRPNN,CRPN,CRPNL,QSPEC)
C
C do it in MSTACK junks except when special structure factor
C operations are present.
      IF (QSPEC) THEN
      MSTACK=XRNREF
      ELSE
      MSTACK=MIN(XRNREF,10000)
      END IF
C
      INDEX=ALLHP(INTEG4(MSTACK))
      VSTACK=ALLHP(ICPLX8(MSTACK*CDEPTH))
      LSTACK=ALLHP(ILOGIC(MSTACK*CDEPTH))
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
      HPSTORE=0
      CALL XDOEVAL(CRPNMX,CRPNN,CRPNX,CRPN,CRPNL,CRPNDB,
     &     CRPNMLT,CRPNLEV,
     &     CRPNTYP,CRPNDOM,VLEVEL,CDEPTH,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &     XRCELL,XRVOL)
C
C modify TSEL array according to selection
      CALL XDOINCV(VLEVEL,CDEPTH,HEAP(LSTACK),NSELE,TSEL,HEAP(INDEX))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*CDEPTH))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*CDEPTH))
      CALL FREHP(INDEX,INTEG4(MSTACK))
      END IF
C
C check if cross-validation is turned on (within selected set)
      XCVTEST=.FALSE.
      DO REFLCT=1,XRNREF
      IF (TSEL(REFLCT).EQ.-1) THEN
      XCVTEST=.TRUE.
      END IF
      END DO
C
C determine number of selected reflections
      XRIREF=0
      DO REFLCT=1,XRNREF
      IF (ABS(TSEL(REFLCT)).EQ.1) THEN
      XRIREF=XRIREF+1
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOIN(VLEVEL,VMAX,LSTACK,N,ARRAY,INDEX)
C
C assign integer array from stack.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER VLEVEL, VMAX, N
      LOGICAL LSTACK(N,*)
      INTEGER ARRAY(*)
      INTEGER INDEX(*)
C local
      INTEGER I
C begin
      DO I=1,N
      IF (LSTACK(I,VLEVEL)) THEN
      ARRAY(INDEX(I))=1
      ELSE
      ARRAY(INDEX(I))=0
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XMONCOU(XRNREF,TSEL,ITEST,COUNT)
C
C assign integer array from stack.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER XRNREF, TSEL(*), ITEST, COUNT
C local
      INTEGER I
C begin
      COUNT=0
      DO I=1,XRNREF
      IF (TSEL(I).EQ.ITEST) THEN
      COUNT=COUNT+1
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOINCV(VLEVEL,VMAX,LSTACK,N,ARRAY,INDEX)
C
C Set all selected elements to -1 if they were initially set to
C 1.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER VLEVEL, VMAX, N
      LOGICAL LSTACK(N,*)
      INTEGER ARRAY(*)
      INTEGER INDEX(*)
C local
      INTEGER I
C begin
      DO I=1,N
      IF (LSTACK(I,VLEVEL).AND.ARRAY(INDEX(I)).EQ.1) THEN
      ARRAY(INDEX(I))=-1
      END IF
      END DO
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDOE2E2(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                   HPTYPE,HPMULT,XRNSYM,QHERM,
     &                   XRTR,HPH,HPK,HPL,MBINS,
     &                   XBINLOW,XBINHIGH,BINSHELL)
C
C
C Routine to calculate the Xray target value from the
C correlation coefficient between Eobs^2 and Ecalc^2
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  e2e2(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL QHERM
C
C pointer
      INTEGER FOBAVG, FCAAVG, NUMAVG, ISHELL
C
C allocate space for temporary arrays
      FOBAVG=ALLHP(IREAL8(MBINS+1))
      FCAAVG=ALLHP(IREAL8(MBINS+1))
      NUMAVG=ALLHP(IREAL8(MBINS+1))
      ISHELL=ALLHP(INTEG4(N))
C
      CALL XDOE2E22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &              INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &              XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,HEAP(ISHELL),
     &              HEAP(FOBAVG),HEAP(FCAAVG),HEAP(NUMAVG))
C
C deallocate space
      CALL FREHP(ISHELL,INTEG4(N))
      CALL FREHP(NUMAVG,IREAL8(MBINS+1))
      CALL FREHP(FCAAVG,IREAL8(MBINS+1))
      CALL FREHP(FOBAVG,IREAL8(MBINS+1))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOE2E22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                    INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                    XRTR,XRH,XRK,XRL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL,ISHELL,
     &                    FOBAVG,FCAAVG,NUMAVG)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION FOBAVG(*), FCAAVG(*), NUMAVG(*)
C local
      INTEGER IND, I
      DOUBLE PRECISION WSUM
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, TEMP, EPSILON, WT, CSUM, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize averages
      DO IND=1,MBINS
      FOBAVG(IND)=ZERO
      FCAAVG(IND)=ZERO
      NUMAVG(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute the Fcalc^2 and Fobs^2 averages for each bin
         FOBAVG(ISHELL(I))=FOBAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FOBS**2
         FCAAVG(ISHELL(I))=FCAAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FCALC**2
         NUMAVG(ISHELL(I))=NUMAVG(ISHELL(I))+WT
      END DO
C
      DO IND=1,MBINS
      IF (NUMAVG(IND).GT.RSMALL.AND.
     &    ABS(FOBAVG(IND)).GT.RSMALL.AND.
     &    ABS(FCAAVG(IND)).GT.RSMALL) THEN
         FOBAVG(IND)=FOBAVG(IND)/NUMAVG(IND)
         FCAAVG(IND)=FCAAVG(IND)/NUMAVG(IND)
      ELSE
C
C set FOBAVG and FCAAVG to one to avoid division by zero
C (this is ok since the numerators will be zero)
         FOBAVG(IND)=ONE
         FCAAVG(IND)=ONE
      END IF
      END DO
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute E^2's
         IFOBS=FOBS**2/(EPSILON*FOBAVG(ISHELL(I)))
         IFCALC=FCALC**2/(EPSILON*FCAAVG(ISHELL(I)))
C
C accumulate information for correlation coefficients
         WSUM=WSUM+ONE
         CI=CI+IFOBS
         CJ=CJ+IFCALC
         CII=CII+IFOBS**2
         CJJ=CJJ+IFCALC**2
         CIJ=CIJ+IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CI).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOE2E2-error: sum over first argument is zero'
      ELSE IF (ABS(CJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOE2E2-error: sum over second argument is zero'
      END IF
C
C calculate denominator and numerator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate correlation - return 1-Correlation
C this is divided by the number of reflections because the
C target expression routine sums over all elements to give to total
C target value
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
         TEMP=(ONE-CORR)/N
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(TEMP,ZERO)
         END DO
      ELSE
         TEMP=ZERO
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODE2E2(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                    HPTYPE,HPMULT,XRNSYM,QHERM,
     &                    XRTR,HPH,HPK,HPL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL)
C
C Routine to calculate the derivative of the Xray target value from the
C correlation coefficient between Eobs^2 and Ecalc^2
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  de2e2(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL QHERM
C
C pointer
      INTEGER FOBAVG, FCAAVG, NUMAVG, DERAVG, ISHELL, DERIV
C
C allocate space for temporary arrays
      FOBAVG=ALLHP(IREAL8(MBINS+1))
      FCAAVG=ALLHP(IREAL8(MBINS+1))
      NUMAVG=ALLHP(IREAL8(MBINS+1))
      DERAVG=ALLHP(IREAL8(MBINS+1))
      ISHELL=ALLHP(INTEG4(N))
      DERIV=ALLHP(IREAL8(N))
C
C
      CALL XDODE2E22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &               INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &               XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &               XBINLOW,XBINHIGH,BINSHELL,HEAP(ISHELL),
     &               HEAP(FOBAVG),HEAP(FCAAVG),HEAP(NUMAVG),
     &               HEAP(DERAVG),HEAP(DERIV))
C
C deallocate space
      CALL FREHP(DERIV,IREAL8(N))
      CALL FREHP(ISHELL,INTEG4(N))
      CALL FREHP(DERAVG,IREAL8(MBINS+1))
      CALL FREHP(NUMAVG,IREAL8(MBINS+1))
      CALL FREHP(FCAAVG,IREAL8(MBINS+1))
      CALL FREHP(FOBAVG,IREAL8(MBINS+1))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODE2E22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                     INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                     XRTR,XRH,XRK,XRL,MBINS,
     &                     XBINLOW,XBINHIGH,BINSHELL,ISHELL,
     &                     FOBAVG,FCAAVG,NUMAVG,DERAVG,DERIV)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION FOBAVG(*), FCAAVG(*), NUMAVG(*), DERAVG(*)
      DOUBLE PRECISION DERIV(*)
C local
      INTEGER IND, I
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ, WSUM
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, CSUM, EPSILON, WT, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize averages
      DO IND=1,MBINS
      FOBAVG(IND)=ZERO
      FCAAVG(IND)=ZERO
      NUMAVG(IND)=ZERO
      DERAVG(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         DERIV(I)=ZERO
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute the Fcalc^2 and Fobs^2 averages for each bin
         FOBAVG(ISHELL(I))=FOBAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FOBS**2
         FCAAVG(ISHELL(I))=FCAAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FCALC**2
         NUMAVG(ISHELL(I))=NUMAVG(ISHELL(I))+WT
      END DO
C
      DO IND=1,MBINS
      IF (NUMAVG(IND).GT.RSMALL.AND.
     &    ABS(FOBAVG(IND)).GT.RSMALL.AND.
     &    ABS(FCAAVG(IND)).GT.RSMALL) THEN
         FOBAVG(IND)=FOBAVG(IND)/NUMAVG(IND)
         FCAAVG(IND)=FCAAVG(IND)/NUMAVG(IND)
      ELSE
C
C set FOBAVG and FCAAVG to one to avoid division by zero
C (this is ok since the numerators will be zero)
         FOBAVG(IND)=ONE
         FCAAVG(IND)=ONE
      END IF
      END DO
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
C calculate correlation sums
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute E^2's
         IFOBS=FOBS**2/(EPSILON*FOBAVG(ISHELL(I)))
         IFCALC=FCALC**2/(EPSILON*FCAAVG(ISHELL(I)))
C
C accumulate information for correlation coefficients
         WSUM=WSUM+ONE
         CI=CI+IFOBS
         CJ=CJ+IFCALC
         CII=CII+IFOBS**2
         CJJ=CJJ+IFCALC**2
         CIJ=CIJ+IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CI).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODE2E2-error: sum over first argument is zero'
      ELSE IF (ABS(CJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODE2E2-error: sum over second argument is zero'
      END IF
C
C compute numerator and denominator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate partial derivatives
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
C
         DO I=1,N
            FOBS=DBLE(VSTACK(I,VLEVEL))
            FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
            IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
               EPSILON=2*XRNSYM/MULT(INDEX(I))
               WT=TWO
            ELSE
C centric or not hermitian symmetry
               EPSILON=XRNSYM/MULT(INDEX(I))
               WT=ONE
            END IF
C
C compute E^2's
            IFOBS=FOBS**2/(EPSILON*FOBAVG(ISHELL(I)))
            IFCALC=FCALC**2/(EPSILON*FCAAVG(ISHELL(I)))
C
C derivative wrt. Ecalc^2
            DERIV(I)=-((IFOBS-CI/WSUM)/DSUM - (CORR/DSUM**2)*
     &                 (CII-CI**2/WSUM) * (IFCALC-CJ/WSUM))
C
C bin wise averages
            DERAVG(ISHELL(I))=DERAVG(ISHELL(I)) +
     &              DERIV(I) * FCALC**2/EPSILON
         END DO
C
C check that we can calculate the bin average
         DO IND=1,MBINS
            IF (NUMAVG(IND)*FCAAVG(IND).GT.RSMALL) THEN
               DERAVG(IND)=DERAVG(IND)/(NUMAVG(IND)*(FCAAVG(IND))**2)
            END IF
         END DO
C
C calculate derivatives wrt Fcalc
         DO I=1,N
C
C compute epsilon and weight
            IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
               EPSILON=2*XRNSYM/MULT(INDEX(I))
               WT=TWO
            ELSE
C centric or not hermitian symmetry
               EPSILON=XRNSYM/MULT(INDEX(I))
               WT=ONE
            END IF
C
C final derivative wrt Fcalc returned in VSTACK
            VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+1)*TWO*(DERIV(I)/
     &           (FCAAVG(ISHELL(I))*EPSILON) -
     &           DERAVG(ISHELL(I))*WT/EPSILON)
         END DO
      ELSE
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDOE1E1(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                   HPTYPE,HPMULT,XRNSYM,QHERM,
     &                   XRTR,HPH,HPK,HPL,MBINS,
     &                   XBINLOW,XBINHIGH,BINSHELL)
C
C
C Routine to calculate the Xray target value from the
C correlation coefficient between Eobs and Ecalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  e1e1(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL QHERM
C
C pointer
      INTEGER FOBAVG, FCAAVG, NUMAVG, ISHELL
C
C allocate space for temporary arrays
      FOBAVG=ALLHP(IREAL8(MBINS+1))
      FCAAVG=ALLHP(IREAL8(MBINS+1))
      NUMAVG=ALLHP(IREAL8(MBINS+1))
      ISHELL=ALLHP(INTEG4(N))
C
      CALL XDOE1E12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &              INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &              XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,HEAP(ISHELL),
     &              HEAP(FOBAVG),HEAP(FCAAVG),HEAP(NUMAVG))
C
C deallocate space
      CALL FREHP(ISHELL,INTEG4(N))
      CALL FREHP(NUMAVG,IREAL8(MBINS+1))
      CALL FREHP(FCAAVG,IREAL8(MBINS+1))
      CALL FREHP(FOBAVG,IREAL8(MBINS+1))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOE1E12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                    INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                    XRTR,XRH,XRK,XRL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL,ISHELL,
     &                    FOBAVG,FCAAVG,NUMAVG)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION FOBAVG(*), FCAAVG(*), NUMAVG(*)
C local
      INTEGER IND, I
      DOUBLE PRECISION WSUM
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, EPSILON, WT, TEMP, CSUM, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize averages
      DO IND=1,MBINS
      FOBAVG(IND)=ZERO
      FCAAVG(IND)=ZERO
      NUMAVG(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute the Fcalc^2 and Fobs^2 averages for each bin
         FOBAVG(ISHELL(I))=FOBAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FOBS**2
         FCAAVG(ISHELL(I))=FCAAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FCALC**2
         NUMAVG(ISHELL(I))=NUMAVG(ISHELL(I))+WT
      END DO
C
      DO IND=1,MBINS
      IF (NUMAVG(IND).GT.RSMALL.AND.
     &    ABS(FOBAVG(IND)).GT.RSMALL.AND.
     &    ABS(FCAAVG(IND)).GT.RSMALL) THEN
         FOBAVG(IND)=FOBAVG(IND)/NUMAVG(IND)
         FCAAVG(IND)=FCAAVG(IND)/NUMAVG(IND)
      ELSE
C
C set FOBAVG and FCAAVG to one to avoid division by zero
C (this is ok since the numerators will be zero)
         FOBAVG(IND)=ONE
         FCAAVG(IND)=ONE
      END IF
      END DO
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute E^2's
         IFOBS=SQRT(FOBS**2/(EPSILON*FOBAVG(ISHELL(I))))
         IFCALC=SQRT(FCALC**2/(EPSILON*FCAAVG(ISHELL(I))))
C
C accumulate information for correlation coefficients
         WSUM=WSUM+ONE
         CI=CI+IFOBS
         CJ=CJ+IFCALC
         CII=CII+IFOBS**2
         CJJ=CJJ+IFCALC**2
         CIJ=CIJ+IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CI).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOE1E1-error: sum over first argument is zero'
      ELSE IF (ABS(CJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOE1E1-error: sum over second argument is zero'
      END IF
C
C calculate denominator and numerator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate correlation - return 1-Correlation
C this is divided by the number of reflections because the
C target expression routine sums over all elements to give to total
C target value
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
         TEMP=(ONE-CORR)/N
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(TEMP,ZERO)
         END DO
      ELSE
         TEMP=ZERO
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODE1E1(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                    HPTYPE,HPMULT,XRNSYM,QHERM,
     &                    XRTR,HPH,HPK,HPL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL)
C
C Routine to calculate the derivative of the Xray target value from the
C correlation coefficient between Eobs and Ecalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  de1e1(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL QHERM
C
C pointer
      INTEGER FOBAVG, FCAAVG, NUMAVG, DERAVG, ISHELL, DERIV
C
C allocate space for temporary arrays
      FOBAVG=ALLHP(IREAL8(MBINS+1))
      FCAAVG=ALLHP(IREAL8(MBINS+1))
      NUMAVG=ALLHP(IREAL8(MBINS+1))
      DERAVG=ALLHP(IREAL8(MBINS+1))
      ISHELL=ALLHP(INTEG4(N))
      DERIV=ALLHP(IREAL8(N))
C
C
      CALL XDODE1E12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &               INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &               XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &               XBINLOW,XBINHIGH,BINSHELL,HEAP(ISHELL),
     &               HEAP(FOBAVG),HEAP(FCAAVG),HEAP(NUMAVG),
     &               HEAP(DERAVG),HEAP(DERIV))
C
C deallocate space
      CALL FREHP(DERIV,IREAL8(N))
      CALL FREHP(ISHELL,INTEG4(N))
      CALL FREHP(DERAVG,IREAL8(MBINS+1))
      CALL FREHP(NUMAVG,IREAL8(MBINS+1))
      CALL FREHP(FCAAVG,IREAL8(MBINS+1))
      CALL FREHP(FOBAVG,IREAL8(MBINS+1))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODE1E12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                     INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                     XRTR,XRH,XRK,XRL,MBINS,
     &                     XBINLOW,XBINHIGH,BINSHELL,ISHELL,
     &                     FOBAVG,FCAAVG,NUMAVG,DERAVG,DERIV)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION FOBAVG(*), FCAAVG(*), NUMAVG(*), DERAVG(*)
      DOUBLE PRECISION DERIV(*)
C local
      INTEGER IND, I
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ, WSUM
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, CSUM, EPSILON, WT, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize averages
      DO IND=1,MBINS
      FOBAVG(IND)=ZERO
      FCAAVG(IND)=ZERO
      NUMAVG(IND)=ZERO
      DERAVG(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         DERIV(I)=ZERO
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute the Fcalc^2 and Fobs^2 averages for each bin
         FOBAVG(ISHELL(I))=FOBAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FOBS**2
         FCAAVG(ISHELL(I))=FCAAVG(ISHELL(I)) +
     &                            (WT/EPSILON)*FCALC**2
         NUMAVG(ISHELL(I))=NUMAVG(ISHELL(I))+WT
      END DO
C
      DO IND=1,MBINS
      IF (NUMAVG(IND).GT.RSMALL.AND.
     &    ABS(FOBAVG(IND)).GT.RSMALL.AND.
     &    ABS(FCAAVG(IND)).GT.RSMALL) THEN
         FOBAVG(IND)=FOBAVG(IND)/NUMAVG(IND)
         FCAAVG(IND)=FCAAVG(IND)/NUMAVG(IND)
      ELSE
C
C set FOBAVG and FCAAVG to one to avoid division by zero
C (this is ok since the numerators will be zero)
         FOBAVG(IND)=ONE
         FCAAVG(IND)=ONE
      END IF
      END DO
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
C calculate correlation sums
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
         IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
            EPSILON=2*XRNSYM/MULT(INDEX(I))
            WT=TWO
         ELSE
C centric or not hermitian symmetry
            EPSILON=XRNSYM/MULT(INDEX(I))
            WT=ONE
         END IF
C
C compute E^2's
         IFOBS=SQRT(FOBS**2/(EPSILON*FOBAVG(ISHELL(I))))
         IFCALC=SQRT(FCALC**2/(EPSILON*FCAAVG(ISHELL(I))))
C
C accumulate information for correlation coefficients
         WSUM=WSUM+ONE
         CI=CI+IFOBS
         CJ=CJ+IFCALC
         CII=CII+IFOBS**2
         CJJ=CJJ+IFCALC**2
         CIJ=CIJ+IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CI).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODE1E1-error: sum over first argument is zero'
      ELSE IF (ABS(CJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODE1E1-error: sum over second argument is zero'
      END IF
C
C compute numerator and denominator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate partial derivatives
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
C
         DO I=1,N
            FOBS=DBLE(VSTACK(I,VLEVEL))
            FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute epsilon and weight
            IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
               EPSILON=2*XRNSYM/MULT(INDEX(I))
               WT=TWO
            ELSE
C centric or not hermitian symmetry
               EPSILON=XRNSYM/MULT(INDEX(I))
               WT=ONE
            END IF
C
C compute E^2's
            IFOBS=SQRT(FOBS**2/(EPSILON*FOBAVG(ISHELL(I))))
            IFCALC=SQRT(FCALC**2/(EPSILON*FCAAVG(ISHELL(I))))
C
C derivative of target wrt. E
            DERIV(I)=-((IFOBS-CI/WSUM)/DSUM - (CORR/DSUM**2)*
     &                 (CII-CI**2/WSUM) * (IFCALC-CJ/WSUM))
C
C derivative of E wrt E^2
            IF (ABS(IFCALC).GT.RSMALL) THEN
               DERIV(I)=DERIV(I)/(TWO*IFCALC)
            ELSE
               DERIV(I)=ZERO
            END IF
C
C bin wise averages
            DERAVG(ISHELL(I))=DERAVG(ISHELL(I)) +
     &              DERIV(I) * FCALC**2/EPSILON
         END DO
C
C check that we can calculate the bin average
         DO IND=1,MBINS
            IF (NUMAVG(IND)*FCAAVG(IND).GT.RSMALL) THEN
               DERAVG(IND)=DERAVG(IND)/(NUMAVG(IND)*(FCAAVG(IND))**2)
            END IF
         END DO
C
C calculate derivatives wrt Fcalc
         DO I=1,N
C
C compute epsilon and weight
            IF (TYPE(INDEX(I)).GE.1.AND.QHERM) THEN
C acentric
               EPSILON=2*XRNSYM/MULT(INDEX(I))
               WT=TWO
            ELSE
C centric or not hermitian symmetry
               EPSILON=XRNSYM/MULT(INDEX(I))
               WT=ONE
            END IF
C
C final derivative wrt Fcalc returned in VSTACK
            VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+1)*TWO*(DERIV(I)/
     &           (FCAAVG(ISHELL(I))*EPSILON) -
     &           DERAVG(ISHELL(I))*WT/EPSILON)
         END DO
      ELSE
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDOF2F2(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                   HPTYPE,HPMULT,XRNSYM,
     &                   XRTR,HPH,HPK,HPL,MBINS,
     &                   XBINLOW,XBINHIGH,BINSHELL)
C
C
C Routine to calculate the Xray target value from the
C correlation coefficient between Fobs^2 and Fcalc^2
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  f2f2(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL XQMULT
C
C default for multiplicity is false
      XQMULT=.FALSE.
C
C check for [MULT=TRUE|FALSE] flag
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XQMULT=.FALSE.
      ELSE
         XQMULT=.TRUE.
      END IF
C
C
      CALL XDOF2F22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &              INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,XQMULT,
     &              XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL)
C
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOF2F22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                    INDEX,TYPE,MULT,XRNSYM,XQMULT,
     &                    XRTR,XRH,XRK,XRL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL XQMULT
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
C local
      INTEGER I
      DOUBLE PRECISION WSUM, TEMP
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, WT, CSUM, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute weight
         IF (XQMULT) THEN
            WT=MULT(INDEX(I))
         ELSE
            WT=ONE
         END IF
C
C compute F^2's
         IFOBS=FOBS**2
         IFCALC=FCALC**2
C
C accumulate information for correlation coefficients
         WSUM=WSUM+WT
         CI=CI+WT*IFOBS
         CJ=CJ+WT*IFCALC
         CII=CII+WT*IFOBS**2
         CJJ=CJJ+WT*IFCALC**2
         CIJ=CIJ+WT*IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CII).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOF2F2-error: sum over first argument is zero'
      ELSE IF (ABS(CJJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOF2F2-error: sum over second argument is zero'
      END IF
C
C calculate denominator and numerator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate correlation - return 1-Correlation
C this is divided by the number of reflections because the
C target expression routine sums over all elements to give to total
C target value
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
         TEMP=(ONE-CORR)/N
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(TEMP,ZERO)
         END DO
      ELSE
         TEMP=ZERO
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODF2F2(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                    HPTYPE,HPMULT,XRNSYM,
     &                    XRTR,HPH,HPK,HPL,MBINS,
     &                    XBINLOW,XBINHIGH,BINSHELL)
C
C Routine to calculate the derivative of the Xray target value from the
C correlation coefficient between Fobs^2 and Fcalc^2
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  df2f2(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNSYM
      INTEGER HPTYPE, HPMULT, HPH, HPK, HPL
      INTEGER MBINS
      DOUBLE PRECISION XRTR(3,3)
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      LOGICAL XQMULT
C
C default for multiplicity is false
      XQMULT=.FALSE.
C
C check for [MULT=TRUE|FALSE] flag
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XQMULT=.FALSE.
      ELSE
         XQMULT=.TRUE.
      END IF
C
      CALL XDODF2F22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &               INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,XQMULT,
     &               XRTR,HEAP(HPH),HEAP(HPK),HEAP(HPL),MBINS,
     &               XBINLOW,XBINHIGH,BINSHELL)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODF2F22(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                     INDEX,TYPE,MULT,XRNSYM,XQMULT,
     &                     XRTR,XRH,XRK,XRL,MBINS,
     &                     XBINLOW,XBINHIGH,BINSHELL)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL XQMULT
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
C local
      INTEGER I
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ, WSUM
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, CSUM, WT, CORR, DERIV
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
      VLEVEL=VLEVEL-1
C
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute weight
         IF (XQMULT) THEN
            WT=MULT(INDEX(I))
         ELSE
            WT=ONE
         END IF
C
C compute F^2's
         IFOBS=FOBS**2
         IFCALC=FCALC**2
C
C accumulate information for correlation coefficients
         WSUM=WSUM+WT
         CI=CI+WT*IFOBS
         CJ=CJ+WT*IFCALC
         CII=CII+WT*IFOBS**2
         CJJ=CJJ+WT*IFCALC**2
         CIJ=CIJ+WT*IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CII).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODF2F2-error: sum over first argument is zero'
      ELSE IF (ABS(CJJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODF2F2-error: sum over second argument is zero'
      END IF
C
C compute numerator and denominator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate partial derivatives
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
C
         DO I=1,N
            FOBS=DBLE(VSTACK(I,VLEVEL))
            FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute weight
            IF (XQMULT) THEN
               WT=MULT(INDEX(I))
            ELSE
               WT=ONE
            END IF
C
C compute F^2's
            IFOBS=FOBS**2
            IFCALC=FCALC**2
C
C derivative wrt. Fcalc
            DERIV=-TWO*WT*((IFOBS-CI/WSUM)/DSUM - (CORR/DSUM**2)*
     &                 (CII-CI**2/WSUM) * (IFCALC-CJ/WSUM))
            VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+1) * DERIV
C
         END DO
      ELSE
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDOF1F1(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                   HPTYPE,HPMULT)
C
C
C Routine to calculate the Xray target value from the
C correlation coefficient between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  f1f1(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*)
      INTEGER HPTYPE, HPMULT
      LOGICAL XQMULT
C
C default for multiplicity is false
      XQMULT=.FALSE.
C
C check for [MULT=TRUE|FALSE] flag
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XQMULT=.FALSE.
      ELSE
         XQMULT=.TRUE.
      END IF
C
      CALL XDOF1F12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &              INDEX,HEAP(HPTYPE),HEAP(HPMULT),XQMULT)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOF1F12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                    INDEX,TYPE,MULT,XQMULT)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*)
      LOGICAL XQMULT
C local
      INTEGER I
      DOUBLE PRECISION WSUM
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, TEMP, WT, CSUM, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute weight
         IF (XQMULT) THEN
            WT=MULT(INDEX(I))
         ELSE
            WT=ONE
         END IF
C
C compute F's
         IFOBS=FOBS
         IFCALC=FCALC
C
C accumulate information for correlation coefficients
         WSUM=WSUM+WT
         CI=CI+WT*IFOBS
         CJ=CJ+WT*IFCALC
         CII=CII+WT*IFOBS**2
         CJJ=CJJ+WT*IFCALC**2
         CIJ=CIJ+WT*IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CII).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOF1F1-error: sum over first argument is zero'
      ELSE IF (ABS(CJJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDOF1F1-error: sum over second argument is zero'
      END IF
C
C calculate denominator and numerator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate correlation - return 1-Correlation
C this is divided by the number of reflections because the
C target expression routine sums over all elements to give to total
C target value
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
         TEMP=(ONE-CORR)/N
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(TEMP,ZERO)
         END DO
      ELSE
         TEMP=ZERO
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODF1F1(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                    HPTYPE,HPMULT)
C
C Routine to calculate the derivative of the Xray target value from the
C correlation coefficient between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  df1f1(Fobs, Fcalc) from script level
C
C Fobs: real
C Fcalc: complex
C
C arguments: Fobs, Fcalc
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*)
      INTEGER HPTYPE, HPMULT
      LOGICAL XQMULT
C
C pointer
      INTEGER DERIV
C
C default for multiplicity is false
      XQMULT=.FALSE.
C
C check for [MULT=TRUE|FALSE] flag
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XQMULT=.FALSE.
      ELSE
         XQMULT=.TRUE.
      END IF
C
C allocate space for temporary arrays
      DERIV=ALLHP(IREAL8(N))
C
C
      CALL XDODF1F12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &               INDEX,HEAP(HPTYPE),HEAP(HPMULT),XQMULT,
     &               HEAP(DERIV))
C
C deallocate space
      CALL FREHP(DERIV,IREAL8(N))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODF1F12(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                     INDEX,TYPE,MULT,XQMULT,
     &                     DERIV)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*)
      LOGICAL XQMULT
      DOUBLE PRECISION DERIV(*)
C local
      INTEGER I
      DOUBLE PRECISION CI, CJ, CII, CJJ, CIJ, WSUM
      DOUBLE PRECISION FOBS, FCALC, IFOBS, IFCALC
      DOUBLE PRECISION DSUM, CSUM, WT, CORR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
C initialize correlation coefficients
      WSUM=ZERO
      CI=ZERO
      CJ=ZERO
      CII=ZERO
      CJJ=ZERO
      CIJ=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-1
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         DERIV(I)=ZERO
C
C compute weight
         IF (XQMULT) THEN
            WT=MULT(INDEX(I))
         ELSE
            WT=ONE
         END IF
C
C compute F's
         IFOBS=FOBS
         IFCALC=FCALC
C
C accumulate information for correlation coefficients
         WSUM=WSUM+WT
         CI=CI+WT*IFOBS
         CJ=CJ+WT*IFCALC
         CII=CII+WT*IFOBS**2
         CJJ=CJJ+WT*IFCALC**2
         CIJ=CIJ+WT*IFOBS*IFCALC
      END DO
C
C do some checking
      IF (ABS(CII).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODF1F1-error: sum over first argument is zero'
      ELSE IF (ABS(CJJ).LT.RSMALL) THEN
         WRITE(6,'(A)')
     &        ' %XDODF1F1-error: sum over second argument is zero'
      END IF
C
C compute numerator and denominator
      DSUM=(CII-CI**2/WSUM)*(CJJ-CJ**2/WSUM)
      CSUM=CIJ-CI*CJ/WSUM
C
C calculate partial derivatives
      IF (DSUM.GT.RSMALL) THEN
         DSUM=SQRT(DSUM)
         CORR=CSUM/DSUM
C
         DO I=1,N
            FOBS=DBLE(VSTACK(I,VLEVEL))
            FCALC=ABS(VSTACK(I,VLEVEL+1))
C
C compute weight
            IF (XQMULT) THEN
               WT=MULT(INDEX(I))
            ELSE
               WT=ONE
            END IF
C
C compute F's
            IFOBS=FOBS
            IFCALC=FCALC
C
C derivative wrt. Fcalc
            IF (IFCALC.GT.RSMALL) THEN
               DERIV(I)=-WT*((IFOBS-CI/WSUM)/DSUM - (CORR/DSUM**2)*
     &              (CII-CI**2/WSUM) * (IFCALC-CJ/WSUM))/IFCALC
            ELSE
               DERIV(I)=ZERO
            END IF
C
C final derivative wrt Fcalc returned in VSTACK
            VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+1) * DERIV(I)
         END DO
      ELSE
         DO I=1,N
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END DO
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDORESI(VLEVEL,VMAX,VSTACK,N,ARGS)
C
C
C Routine to calculate the Xray target value from the
C standard least squares residual between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  resi(Fobs, Fcalc, Weight) from script level
C
C Fobs: real
C Fcalc: complex
C Weight: real
C
C arguments: Fobs, Fcalc, Weight
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
C local
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C default is floating scale factor
      XRFFK=ONE
      QKFIXED=.FALSE.
C
C check for [K=<real>]
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XRFFK=ONE
         QKFIXED=.FALSE.
      ELSE
         XRFFK=DBLE(ARGS(1))
         QKFIXED=.TRUE.
      END IF
C
      CALL XDORESI2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDORESI2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C local
      INTEGER I
      DOUBLE PRECISION FOBS, FCALC, WEIGHT, TARGET, SCALE
      DOUBLE PRECISION FOCSUM, FO2SUM, FC2SUM
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      FOCSUM=ZERO
      FO2SUM=ZERO
      FC2SUM=ZERO
      TARGET=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-2
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
C
C compute sums for scaling
         FOCSUM=FOCSUM+WEIGHT*FOBS*FCALC
         FO2SUM=FO2SUM+WEIGHT*FOBS*FOBS
         FC2SUM=FC2SUM+WEIGHT*FCALC*FCALC
      END DO
C
C if scale factor is not fixed then calculate - resets XRFFK
      IF (.NOT.QKFIXED) THEN
         IF (FC2SUM.GT.RSMALL) THEN
            XRFFK=FOCSUM/FC2SUM
         ELSE
            WRITE(6,'(A)')
     &           ' %XDORESI-error: sum over denominator is zero'
            CALL WRNDIE(-5,'XDORESI',
     &           'could not calculate scale factor XRFFK')
         END IF
      END IF
C
C calculate 1/Fo^2 normalization scale factor
      IF (FO2SUM.GT.RSMALL) THEN
         SCALE=ONE/FO2SUM
      ELSE
         SCALE=ONE
      END IF
C
C calculate residual target
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
         TARGET=TARGET+(WEIGHT*SCALE*(FOBS-XRFFK*FCALC)**2)
      END DO
C
C return target divided by the number of reflections in VSTACK
      DO I=1,N
         VSTACK(I,VLEVEL)=DCMPLX((TARGET/N),ZERO)
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODRESI(VLEVEL,VMAX,VSTACK,N,ARGS)
C
C
C Routine to calculate the derivative of the Xray target value from the
C standard least squares residual between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  dresi(Fobs, Fcalc, Weight) from script level
C
C Fobs: real
C Fcalc: complex
C Weight: real
C
C arguments: Fobs, Fcalc, Weight
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
C local
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C default is floating scale factor
      XRFFK=ONE
      QKFIXED=.FALSE.
C
C check for [K=<real>]
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XRFFK=ONE
         QKFIXED=.FALSE.
      ELSE
         XRFFK=DBLE(ARGS(1))
         QKFIXED=.TRUE.
      END IF
C
      CALL XDODRESI2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODRESI2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C local
      INTEGER I
      DOUBLE PRECISION FOBS, FCALC, WEIGHT, TEMP, SCALE
      DOUBLE PRECISION FOCSUM, FO2SUM, FC2SUM
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      FOCSUM=ZERO
      FO2SUM=ZERO
      FC2SUM=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-2
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
C
C compute sums for scaling
         FOCSUM=FOCSUM+WEIGHT*FOBS*FCALC
         FO2SUM=FO2SUM+WEIGHT*FOBS*FOBS
         FC2SUM=FC2SUM+WEIGHT*FCALC*FCALC
      END DO
C
C if scale factor is not fixed then calculate - resets XRFFK
      IF (.NOT.QKFIXED) THEN
         IF (FC2SUM.GT.RSMALL) THEN
            XRFFK=FOCSUM/FC2SUM
         ELSE
            WRITE(6,'(A)')
     &           ' %XDODRESI-error: sum over denominator is zero'
            CALL WRNDIE(-5,'XDODRESI',
     &           'could not calculate scale factor XRFFK')
         END IF
      END IF
C
C calculate 1/Fo^2 normalization scale factor
      IF (FO2SUM.GT.RSMALL) THEN
         SCALE=ONE/FO2SUM
      ELSE
         SCALE=ONE
      END IF
C
C calculate derivatives of the residual target
      DO I=1,N
         FOBS=DBLE(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
         TEMP=WEIGHT*SCALE*(FOBS-XRFFK*FCALC)
         IF (FCALC.GT.RSMALL) THEN
            VSTACK(I,VLEVEL)=-TWO*XRFFK*TEMP*
     &                        VSTACK(I,VLEVEL+1)/FCALC
         ELSE
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END IF
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDOVECT(VLEVEL,VMAX,VSTACK,N,ARGS)
C
C
C Routine to calculate the Xray target value from the
C least squares vector residual between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  vector(Fobs, Fcalc, Weight) from script level
C
C Fobs: complex
C Fcalc: complex
C Weight: real
C
C arguments: Fobs, Fcalc, Weight
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
C local
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C default is floating scale factor
      XRFFK=ONE
      QKFIXED=.FALSE.
C
C check for [K=<real>]
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XRFFK=ONE
         QKFIXED=.FALSE.
      ELSE
         XRFFK=DBLE(ARGS(1))
         QKFIXED=.TRUE.
      END IF
C
      CALL XDOVECT2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOVECT2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C local
      INTEGER I
      DOUBLE PRECISION FOBS, FCALC, WEIGHT, TARGET, SCALE
      DOUBLE PRECISION FOCSUM, FO2SUM, FC2SUM
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      FOCSUM=ZERO
      FO2SUM=ZERO
      FC2SUM=ZERO
      TARGET=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-2
      DO I=1,N
         FOBS=ABS(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
C
C compute sums for scaling
         FOCSUM=FOCSUM+WEIGHT*FOBS*FCALC
         FO2SUM=FO2SUM+WEIGHT*FOBS*FOBS
         FC2SUM=FC2SUM+WEIGHT*FCALC*FCALC
      END DO
C
C if scale factor is not fixed then calculate - resets XRFFK
      IF (.NOT.QKFIXED) THEN
         IF (FC2SUM.GT.RSMALL) THEN
            XRFFK=FOCSUM/FC2SUM
         ELSE
            WRITE(6,'(A)')
     &           ' %XDOVECT-error: sum over denominator is zero'
            CALL WRNDIE(-5,'XDOVECT',
     &           'could not calculate scale factor XRFFK')
         END IF
      END IF
C
C calculate 1/Fo^2 normalization scale factor
      IF (FO2SUM.GT.RSMALL) THEN
         SCALE=ONE/FO2SUM
      ELSE
         SCALE=ONE
      END IF
C
C calculate vector residual target
      DO I=1,N
         FOBS=ABS(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
         TARGET=TARGET+WEIGHT*SCALE*(ABS(VSTACK(I,VLEVEL)-
     &               XRFFK*VSTACK(I,VLEVEL+1)))**2
      END DO
C
      DO I=1,N
         VSTACK(I,VLEVEL)=DCMPLX((TARGET/N),ZERO)
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE XDODVECT(VLEVEL,VMAX,VSTACK,N,ARGS)
C
C
C Routine to calculate the derivative of the Xray target value from the
C standard least squares residual between Fobs and Fcalc
C
C Author: Axel Brunger
C
C Modified to create function by Paul Adams
C
C to be used as <sf-obj> =
C  dvector(Fobs, Fcalc, Weight) from script level
C
C Fobs: complex
C Fcalc: complex
C Weight: real
C
C arguments: Fobs, Fcalc, Weight
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
C local
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C default is floating scale factor
      XRFFK=ONE
      QKFIXED=.FALSE.
C
C check for [K=<real>]
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
         XRFFK=ONE
         QKFIXED=.FALSE.
      ELSE
         XRFFK=DBLE(ARGS(1))
         QKFIXED=.TRUE.
      END IF
C
      CALL XDODVECT2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDODVECT2(VLEVEL,VMAX,VSTACK,N,ARGS,QKFIXED,XRFFK)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION XRFFK
      LOGICAL QKFIXED
C local
      INTEGER I
      DOUBLE PRECISION FOBS, FCALC, WEIGHT, SCALE, FAB
      DOUBLE PRECISION FOCSUM, FO2SUM, FC2SUM
      DOUBLE PRECISION AOBS, BOBS, ACALC, BCALC
      DOUBLE COMPLEX CTEMP
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      FOCSUM=ZERO
      FO2SUM=ZERO
      FC2SUM=ZERO
      FAB=ZERO
C
C loop over all selected reflections
      VLEVEL=VLEVEL-2
      DO I=1,N
         FOBS=ABS(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
C
C compute sums for scaling
         FOCSUM=FOCSUM+WEIGHT*FOBS*FCALC
         FO2SUM=FO2SUM+WEIGHT*FOBS*FOBS
         FC2SUM=FC2SUM+WEIGHT*FCALC*FCALC
      END DO
C
C if scale factor is not fixed then calculate - resets XRFFK
      IF (.NOT.QKFIXED) THEN
         IF (FC2SUM.GT.RSMALL) THEN
            XRFFK=FOCSUM/FC2SUM
         ELSE
            WRITE(6,'(A)')
     &           ' %XDODVECT-error: sum over denominator is zero'
            CALL WRNDIE(-5,'XDODVECT',
     &           'could not calculate scale factor XRFFK')
         END IF
      END IF
C
C calculate 1/Fo^2 normalization scale factor
      IF (FO2SUM.GT.RSMALL) THEN
         SCALE=ONE/FO2SUM
      ELSE
         SCALE=ONE
      END IF
C
C compute term for XRFFK derivative
C
C do not calculate this if the scale factor is fixed
      IF (.NOT.QKFIXED) THEN
         DO I=1,N
            AOBS=DBLE(VSTACK(I,VLEVEL))
            BOBS=DIMAG(VSTACK(I,VLEVEL))
            ACALC=DBLE(VSTACK(I,VLEVEL+1))
            BCALC=DIMAG(VSTACK(I,VLEVEL+1))
            WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
            FAB=FAB+WEIGHT*((AOBS-XRFFK*ACALC)*ACALC +
     &                      (BOBS-XRFFK*BCALC)*BCALC)
         END DO
      ELSE
         FAB=ZERO
      END IF
C
C calculate derivative of vector residual target
      DO I=1,N
         FOBS=ABS(VSTACK(I,VLEVEL))
         FCALC=ABS(VSTACK(I,VLEVEL+1))
         WEIGHT=DBLE(VSTACK(I,VLEVEL+2))
         CTEMP=-TWO*WEIGHT*SCALE*XRFFK*
     &          (VSTACK(I,VLEVEL)-XRFFK*VSTACK(I,VLEVEL+1))
         IF (FCALC.GT.RSMALL) THEN
            IF (.NOT.QKFIXED) THEN
               VSTACK(I,VLEVEL)=CTEMP-TWO*WEIGHT*SCALE*FAB*
     &           (FOBS/FC2SUM-TWO*FOCSUM*FCALC/FC2SUM**2)*
     &            VSTACK(I,VLEVEL+1)/FCALC
            ELSE
C
C do not include the correction for the XRFFK derivative if fixed
               VSTACK(I,VLEVEL)=CTEMP
            END IF
         ELSE
            VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
         END IF
      END DO
C
      RETURN
      END
C
C=====================================================================
C
C
