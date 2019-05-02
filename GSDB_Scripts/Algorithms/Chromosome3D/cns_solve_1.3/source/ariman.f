C
C=====================================================================
C
      SUBROUTINE ARIMAN(NOEDIS,NOELOW,NOEHIG,NOEWGH,NOEVOL,NOECV)
C
C manipulate the NOE database
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vector.inc'
C
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
      DOUBLE PRECISION NOEVOL(*),NOEWGH(*)
      INTEGER NOECV(*)
C
C     local variable declarations
C
      INTEGER TPMAX, NNOE
      CHARACTER*4 MODE
      LOGICAL XRERR
      INTEGER IDUMMY
C
      XRERR = .FALSE.
      TPMAX = WDDMAX
C
      CALL NEXTQL('ARIA-DO>')
      IF ( WD(1:4) .EQ. 'HELP' ) THEN
C
         CALL CNSHELP('cns-aria-do')
C
      RETURN
      ELSE
      CALL SAVEWD
      ENDIF
C
C     get mode if not defaulted
C
      CALL NEXTQL('ARIA-DO>')
      IF ( WD(1:4) .EQ. '(   ' ) THEN
      CALL SAVEWD
      MODE = 'REAL'
      ELSE
      CALL DSPERR('ARIA-DO>','Illegal computational mode')
      RETURN
      ENDIF
C
C     parse the equation into the tables
C
      CALL HIER(MODE, XRERR)
      IF ( XRERR ) RETURN
      DO NNOE=1,NOENUM
C
C step a) store or obtain parsing tables
C
      CALL TABMGT( NNOE, XRERR )
      IF ( XRERR ) GOTO 999
C
C step b) fill values in for variables
C
      CALL NOTBVL( MODE, NNOE, XRERR, NOEDIS,NOELOW,NOEHIG,
     &  NOEWGH,NOEVOL,NOECV)
      IF ( XRERR ) GOTO 999
C
C         step c) evaluate the table
C
      CALL EVALTB( MODE, NNOE, IDUMMY, IDUMMY, XRERR )
      IF ( XRERR ) GOTO 999
C
C step d) assign the table result to the array on left hand side
C
      CALL NOASSN( NNOE, MODE, XRERR, NOEDIS,NOELOW,NOEHIG,
     &  NOEWGH,NOEVOL,NOECV)
      IF ( XRERR ) GOTO 999
      ENDDO
      RETURN
C
C     print out error message and return
C
999   CALL DSPERR('NOEMAN>','Assignment aborted')
      WRITE(6,'(A,I7)') ' Starting with reflection index: ', NNOE
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NOTBVL
     &  (MODE,NNOE,XRERR,NOEDIS,NOELOW,NOEHIG,NOEWGH,NOEVOL,NOECV)
C
      IMPLICIT NONE
C
C     if the variable is a double precision variable store the double
C     precision variable in the table's double precision variable.
C
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vector.inc'
C
      INTEGER NNOE
      CHARACTER*4 MODE
      LOGICAL XRERR
      DOUBLE PRECISION NOEDIS(*), NOELOW(*), NOEHIG(*)
      DOUBLE PRECISION NOEVOL(*), NOEWGH(*)
      INTEGER NOECV(*)
C
      INTEGER I
C
      CHARACTER*20 CNOEDIS,CNOEHIG,CNOELOW,CNOEVOL,CNOEWGH,CNOECV
      PARAMETER ( CNOEDIS = 'NOE_DISTANCE        ')
      PARAMETER ( CNOELOW = 'NOE_LOWER-ERROR     ')
      PARAMETER ( CNOEHIG = 'NOE_HIGHER-ERROR    ')
      PARAMETER ( CNOEVOL = 'NOE_VOLUME          ')
      PARAMETER ( CNOEWGH = 'NOE_WEIGHT          ')
      PARAMETER ( CNOECV  = 'NOE_TEST            ')
C
      DO I = 1,NUMOPS
      IF ( VAREQ(I) .EQ. CNOEDIS ) THEN
      RVAREQ(I) = NOEDIS(NNOE)
      ELSE IF (VAREQ(I).EQ.CNOEHIG) THEN
      RVAREQ(I)=NOEHIG(NNOE)
      ELSE IF (VAREQ(I).EQ.CNOELOW) THEN
      RVAREQ(I)=NOELOW(NNOE)
      ELSE IF (VAREQ(I).EQ.CNOEVOL) THEN
      RVAREQ(I)=NOEVOL(NNOE)
      ELSE IF (VAREQ(I).EQ.CNOEWGH) THEN
      RVAREQ(I)=NOEWGH(NNOE)
      ELSE IF (VAREQ(I).EQ.CNOECV) THEN
      RVAREQ(I)=NOECV(NNOE)
      ENDIF
      VARTYP(I) = 'DP'
      ENDDO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NOASSN(NNOE,MODE,XRERR,NOEDIS,NOELOW,NOEHIG,
     &                  NOEWGH,NOEVOL,NOECV)
C
      IMPLICIT NONE
C
C     the table has been evaluated leaving the final result in either
C     rvareq(1)
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'vector.inc'
C
C     argument declarations
C
      INTEGER NNOE
      CHARACTER*4 MODE
      LOGICAL XRERR
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*),
     &  NOEVOL(*),NOEWGH(*)
      INTEGER NOECV(*)
C
C     local variable declarations
C
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      CHARACTER*20 CNOEDIS,CNOEHIG,CNOELOW,CNOEVOL,CNOEWGH,CNOECV
      PARAMETER ( CNOEDIS = 'NOE_DISTANCE        ')
      PARAMETER ( CNOELOW = 'NOE_LOWER-ERROR     ')
      PARAMETER ( CNOEHIG = 'NOE_HIGHER-ERROR    ')
      PARAMETER ( CNOEVOL = 'NOE_VOLUME          ')
      PARAMETER ( CNOEWGH = 'NOE_WEIGHT          ')
      PARAMETER ( CNOECV  = 'NOE_TEST            ')
C
      IF ( MODE .EQ. 'REAL' ) THEN
C
      IF ( VAREQ(1) .EQ. CNOEDIS ) THEN
      NOEDIS(NNOE)=RVAREQ(1)
      ELSE IF (VAREQ(1).EQ.CNOEHIG) THEN
      NOEHIG(NNOE)=RVAREQ(1)
      ELSE IF (VAREQ(1).EQ.CNOELOW) THEN
      NOELOW(NNOE)=RVAREQ(1)
      ELSE IF (VAREQ(1).EQ.CNOEVOL) THEN
      NOEVOL(NNOE)=RVAREQ(1)
      ELSE IF (VAREQ(1).EQ.CNOEWGH) THEN
      NOEWGH(NNOE)=RVAREQ(1)
      ELSE IF (VAREQ(1).EQ.CNOECV) THEN
      NOECV(NNOE)=RVAREQ(1)
      ENDIF
C
      END IF
      RETURN
      END
C
C=====================================================================
C
