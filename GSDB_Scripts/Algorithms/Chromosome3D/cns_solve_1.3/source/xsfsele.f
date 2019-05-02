      SUBROUTINE XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,QSELE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
C
C Parses selection expression and defines the logical array
C QSELE.  If QSELE(REFLCT) is true this reflection is selected.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'expression.inc'
      INCLUDE 'comand.inc'
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QSELE(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      LOGICAL ERR, QSPEC
      CHARACTER*2 TYPE2, DOMAIN2
      INTEGER DEPTH2, NN, REFLCT, NNSELE, MSTACK, START, STOP, VLEVEL
      INTEGER NSELE
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
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
C pointers
      INTEGER INDEX, VSTACK, LSTACK, HPSTORE
C
C begin
      ERR=.FALSE.
C
C check that selection is properly specified.
      CALL NEXTDO('SELEction>')
      IF (WD(1:1).EQ.'=') THEN
      CALL NEXTDO('SELEction>')
      END IF
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('DO>',
     &  ' "(" expected.')
      ERR=.TRUE.
      END IF
C
C parse the selection
      CALL EXRPN('SELEction>',
     & RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,
     & RPNTYP2,RPNDOM2,RPNLEV2,TYPE2,DOMAIN2,DEPTH2,XDOFUNC,XDOOPER,
     & XDOTYPE,QHERM,ERR,
     & XSFNUM,XSFNAM,XSFTYPE,0,' ',' ')
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
      CALL DSPERR('SELEction>',
     & 'Data type mismatch.  Selection must be a logical expression.')
      ERR=.TRUE.
      END IF
C
C check that no FT operation is present
      IF (.NOT.ERR) THEN
      DO NN=1,RPNN2
      ERR=ERR.OR.RPN2(1,NN)(1:RPNL2(1,NN)).EQ.'FT'
      END DO
      IF (ERR) THEN
      CALL DSPERR('SELEction',
     &  'No FT allowed in selections.')
      END IF
C
C check if special structure factor operations are present
      CALL XDOSPCL(RPNMX,RPNX,RPNN2,RPN2,RPNL2,QSPEC)
      END IF
C
C initialize QSELE array
      DO REFLCT=1,XRNREF
      QSELE(REFLCT)=.FALSE.
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
      VSTACK=ALLHP(ICPLX8(MSTACK*DEPTH2))
      LSTACK=ALLHP(ILOGIC(MSTACK*DEPTH2))
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
CC
C evaluate selection
      HPSTORE=0
      CALL XDOEVAL(RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,RPNLEV2,
     &     RPNTYP2,RPNDOM2,VLEVEL,DEPTH2,NSELE,HEAP(VSTACK),
     &     HEAP(LSTACK),HEAP(INDEX),XRTR,HPH,HPK,HPL,
     &     XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &     HPMULT,HPTYPE,
     &     HPSTORE,XRMREF,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &     XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,
     &     XRVOL)
C
C copy selection into QSELE array
      CALL XDOLO(VLEVEL,DEPTH2,HEAP(LSTACK),NSELE,QSELE,HEAP(INDEX))
C
      START=START+MSTACK
      STOP=MIN(XRNREF,STOP+MSTACK)
      END DO
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH2))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH2))
      CALL FREHP(INDEX,INTEG4(MSTACK))
C
      NNSELE=0
      DO REFLCT=1,XRNREF
      IF (QSELE(REFLCT)) NNSELE=NNSELE+1
      END DO
C
C print number of selected reflections and declare the symbol
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I9,A)')
     & ' Total of ',NNSELE,' structure factor elements were selected.'
      END IF
C
      DBPREC=NNSELE
      CALL DECLAR( 'SELECT','DP',' ',DBCOMP,DBPREC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOLO(VLEVEL,VMAX,LSTACK,N,ARRAY,INDEX)
C
C assign logical array from stack.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER VLEVEL, VMAX, N
      LOGICAL LSTACK(N,*), ARRAY(*)
      INTEGER INDEX(*)
C local
      INTEGER I
C begin
      DO I=1,N
      ARRAY(INDEX(I))=LSTACK(I,VLEVEL)
      END DO
      RETURN
      END
