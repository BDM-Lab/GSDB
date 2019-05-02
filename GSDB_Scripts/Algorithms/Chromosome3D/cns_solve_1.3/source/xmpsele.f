      SUBROUTINE XMPSELE(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &      HPRHOMA,NRHO,NMASK,
     &      XRNSYM,MPACK,NNSELE)
C
C Parse selection expression and pack the selected indices
C For map selection only
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C version 3.41,  12-DEC-93
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'expression.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, NMASK, NRHO
      INTEGER XRNSYM
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY), NNSELE
C local
      LOGICAL ERR, QSPEC, QSPEC2
      CHARACTER*2 TYPE2, DOMAIN2
      INTEGER DEPTH2, NSELE, NN, START, STOP, DUMMY
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DBCOMP
C
C note: command stacks are stored in expression.inc
C
C definition of the functions, operands, data, and domain types
      EXTERNAL XDOFUNC, XDOOPER, XDOTYPE
C
C variable stack pointers
      INTEGER VLEVEL, VSTACK, LSTACK
C index pointer
      INTEGER INDEXA, INDEXB, INDEXC, A, B, C
C
C size of the variable stack
      INTEGER MSTACK
C
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C begin
      ERR=.FALSE.
C
C the selection prompt is optional
      IF (WD(1:4).EQ.'SELE') THEN
      CALL NEXTDO('SELEction>')
      IF (WD(1:1).EQ.'=') THEN
      CALL NEXTDO('SELEction>')
      END IF
      END IF
C
C opening parenthesis
      IF (WD(1:1).NE.'(') THEN
      CALL DSPERR('SELE>',' "(" expected.')
      ERR=.TRUE.
      END IF
C
C
C parse the selection
      CALL EXRPN('SELEction>',
     & RPNMX,RPNN2,RPNX,RPN2,RPNL2,RPNDB2,RPNMLT2,
     & RPNTYP2,RPNDOM2,RPNLEV2,TYPE2,DOMAIN2,DEPTH2,XDOFUNC,
     & XDOOPER,XDOTYPE,QHERM,ERR,
     & 0,' ',' ',XRHONUM,XRHONAM,' ')
C
C closing parenthesis
      CALL NEXTDO('SELE>')
      IF (WD(1:1).NE.')'.AND..NOT.ERR) THEN
      CALL DSPERR('SELE>',
     & 'item cannot be interpreted or ")" expected.')
      ERR=.TRUE.
      END IF
C
C check data type
      IF (TYPE2.NE.'LO') THEN
      CALL DSPERR('SELE',
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
      CALL DSPERR('SELE','No FT operations allowed in selection.')
      END IF
      END IF
C
C check if special structure factor or map operations are present
      CALL XDOSPCL(RPNMX,RPNX,RPNN2,RPN2,RPNL2,QSPEC2)
      CALL XDOSPCL(RPNMX,RPNX,RPNN,RPN,RPNL,QSPEC)
      QSPEC=QSPEC.OR.QSPEC2
C
C
      IF (.NOT.ERR) THEN
C
C initialize the selection matrix
      NNSELE=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      MPACK(A,B,C)=0
      END DO
      END DO
      END DO
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
      VSTACK=ALLHP(ICPLX8(MSTACK*DEPTH2))
      LSTACK=ALLHP(ILOGIC(MSTACK*DEPTH2))
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
     & VLEVEL,DEPTH2,NSELE,HEAP(VSTACK),
     & HEAP(LSTACK),NA,NB,NC,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,DUMMY,DUMMY,HPRHOMA,
     & HEAP(INDEXA),HEAP(INDEXB),HEAP(INDEXC),XRNSYM)
C
C fill the selection matrix for this junk of data
      CALL XMPDEF(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &       MPACK,VLEVEL,DEPTH2,HEAP(LSTACK),NSELE,HEAP(INDEXA),
     &       HEAP(INDEXB),HEAP(INDEXC),NNSELE)
C
      START=START+MSTACK
      STOP=MIN(NMASK,STOP+MSTACK)
      END DO
Cend loop over junks
C-------------------
C
C deallocate space for the variable stack
      CALL FREHP(VSTACK,ICPLX8(MSTACK*DEPTH2))
      CALL FREHP(LSTACK,ILOGIC(MSTACK*DEPTH2))
      CALL FREHP(INDEXC,INTEG4(MSTACK))
      CALL FREHP(INDEXB,INTEG4(MSTACK))
      CALL FREHP(INDEXA,INTEG4(MSTACK))
C
C
C print number of selected map elements and declare the symbol
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I8,A)')
     & ' XMPSELE: total of ',NNSELE,' map elements were selected.'
      END IF
      DBPREC=NNSELE
      CALL DECLAR( 'SELECT','DP',' ',DBCOMP,DBPREC)
C
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMPPCK0(MPACK,NNSELE,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,RHOMASK)
C
C default selection:  all elements are set to one within the
C asymmetric unit, zero elsewhere
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER NNSELE
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      INTEGER A,B,C
C begin
      NNSELE=0
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      NNSELE=NNSELE+1
      MPACK(A,B,C)=1
      ELSE
      MPACK(A,B,C)=0
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMPDEF(MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &             MPACK,VLEVEL,VMAX,LSTACK,N,INDEXA,INDEXB,INDEXC,
     &             NNSELE)
C
C defines mask (as a result of a selection).
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INTEGER MAASY,MBASY,MCASY,NAASY,NBASY,NCASY
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER VLEVEL, VMAX, N
      LOGICAL LSTACK(N,*)
      INTEGER INDEXA(*), INDEXB(*), INDEXC(*), NNSELE
C local
      INTEGER I, NSELE
C begin
      NSELE=0
      DO I=1,N
      IF (LSTACK(I,VLEVEL)) THEN
      MPACK(INDEXA(I),INDEXB(I),INDEXC(I))=1
      NNSELE=NNSELE+1
      END IF
      END DO
C
      RETURN
      END
