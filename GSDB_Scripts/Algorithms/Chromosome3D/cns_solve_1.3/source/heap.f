      SUBROUTINE QRYALLOC(VARID, BUF)
      IMPLICIT NONE
C I/O
      INTEGER    VARID
      CHARACTER  BUF*(*)
C
C heap statistics query routine
C
C Author: R.W. Grosse-Kunstleve
C
C local
      INTEGER  DIGITS(40), NDIGITS, DASH, I, J
C
C external
      INTEGER   CNSQALLOC
      EXTERNAL  CNSQALLOC
C
C begin
      BUF = ' '
C
          NDIGITS = CNSQALLOC(VARID, DIGITS, 40)
      IF (NDIGITS .EQ. 0) THEN
        CALL WRNDIE(-5, 'QRYALLOC', 'Fatal coding error.')
        CALL DIE
      END IF
C
      DASH = 0
      IF (NDIGITS .LT. 0) THEN
        DASH = 1
        NDIGITS = -NDIGITS
      END IF
C
      IF (DASH + NDIGITS .GT. LEN(BUF)) THEN
        DO I = 1, LEN(BUF)
          BUF(I:I) = '*'
        END DO
        RETURN
      END IF
C
      I = LEN(BUF)
      DO J = 1, NDIGITS
        BUF(I:I) = CHAR(ICHAR('0') + DIGITS(J))
        I = I - 1
      END DO
C
      IF (DASH .NE. 0) BUF(I:I) = '-'
C
      RETURN
      END
C
C=======================================================================
C
      LOGICAL FUNCTION IS0ALLOC()
      IMPLICIT NONE
C
C Check if heap is currently empty
C
C Author: R.W. Grosse-Kunstleve
C
C local
      CHARACTER*20  BUF
C
C begin
      IS0ALLOC = .FALSE.
C
      CALL QRYALLOC( 2, BUF)
      IF (BUF(LEN(BUF)-1:) .NE. ' 0') RETURN
C
      CALL QRYALLOC(-2, BUF)
      IF (BUF(LEN(BUF)-1:) .NE. ' 0') RETURN
C
      IS0ALLOC = .TRUE.
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE CHKALLOC(TYPE, STRLEN, LABEL, OLDELEM, NEWELEM, OFFS)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      CHARACTER*(*)  TYPE
      INTEGER        STRLEN
      CHARACTER*(*)  LABEL
      INTEGER        OLDELEM, NEWELEM, OFFS
C
C low-level heap allocation routine
C
C Author: R.W. Grosse-Kunstleve
C
C local
      INTEGER       STATUS, OLDOFFS
      CHARACTER*11  DBYTES, CURBYTES
      CHARACTER*11  DOVERH, CUROVERH
C
C external
      INTEGER   CNSIALLOC, CNSCALLOC
      EXTERNAL  CNSIALLOC, CNSCALLOC
C
C begin
! serialize heap manipulations
!$omp critical (heap_crit)
      OLDOFFS = OFFS
C
      IF      (TYPE(1:1) .EQ. 'I') THEN
        STATUS = CNSIALLOC(OLDELEM, NEWELEM, OFFS, HEAP)
      ELSE IF (TYPE(1:1) .EQ. 'C') THEN
        STATUS = CNSCALLOC(OLDELEM, NEWELEM, OFFS, STRLEN, C1HEAP)
      END IF
C
      IF (WRNLEV .GE. 15) THEN
        WRITE(6, '(1X, A, 1X, 3A, 1X, I6, 1X, 3A, 4(1X, I11))')
     &    '''@DM''', '''', TYPE, '''', STRLEN, '''', LABEL, '''',
     &    OLDELEM, NEWELEM, OLDOFFS, OFFS
      END IF
C
      IF      (STATUS .EQ. 1) THEN
        CALL QRYALLOC(1, DBYTES)
        WRITE(6, '(1X, 4A)')
     &    LABEL, ': request for ', DBYTES, ' bytes'
        WRITE(6,'(1X,A)')
     & '---------------------------------------------------------'
        WRITE(6,'(1X,A)')
     & 'There is not enough memory available to the program.'
        WRITE(6,'(1X,A)')
     & 'This may be because of too little physical memory (RAM)'
        WRITE(6,'(1X,A)')
     & 'or too little swap space on the machine. It could also be'
        WRITE(6,'(1X,A)')
     & 'the result of user or system limits. On most Unix systems'
        WRITE(6,'(1X,A)')
     & 'the "limit" command can be used to check the current user'
        WRITE(6,'(1X,A)')
     & 'limits. Please check that the datasize, memoryuse and'
        WRITE(6,'(1X,A)')
     & 'vmemoryuse limits are set at a large enough value.'
        WRITE(6,'(1X,A)')
     & '---------------------------------------------------------'
        CALL WRNDIE(-5, LABEL, 'not enough memory available')
        CALL DIE
      ELSE IF (STATUS .EQ. 2) THEN
        CALL WRNDIE(-5, LABEL, 'corrupt heap')
        CALL DIE
      ELSE IF (STATUS .EQ. 3) THEN
        CALL WRNDIE(-5, LABEL, 'out of INTEGER range')
        CALL DIE
      ELSE IF (STATUS .NE. 0) THEN
        WRITE(6, '(1X, A, I6)') 'dmemory error code = ', STATUS
        CALL WRNDIE(-5, LABEL, 'fatal coding error')
        CALL DIE
      END IF
C
      IF (TYPE(1:1) .EQ. 'I' .AND. NEWELEM .GE. 0) THEN
        IF (HEAP(OFFS + NEWELEM) .NE. NEWELEM) THEN
          CALL WRNDIE(-5, LABEL, 'corrupt offset')
          CALL DIE
        END IF
      END IF
C
      IF (WRNLEV .GE. 15) THEN
        CALL QRYALLOC( 1,   DBYTES)
        CALL QRYALLOC( 2, CURBYTES)
        CALL QRYALLOC(-1,   DOVERH)
        CALL QRYALLOC(-2, CUROVERH)
        WRITE(6, '(1X, 2A, I11, 8A)')
     &    LABEL, ': delem=',   MAX(0, NEWELEM)
     &                       - MAX(0, OLDELEM),
     &    ' dbytes=',   DBYTES,
     &    ' doverh=',   DOVERH,
     &    ' total=',  CURBYTES,
     &    ' overh=',  CUROVERH
      END IF
C
!$omp end critical (heap_crit)
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE CHKSHPA(ERROR)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      LOGICAL  ERROR
C
C check static heap arrays
C
C Author: Ralf W. Grosse-Kunstleve
C
C local
      INTEGER  I, J
C
C begin
      ERROR = .TRUE.
C
      IF (HEAP(0) .NE. 13) RETURN
      DO I = 1, 15
        IF (HEAP(I) .NE. MOD(HEAP(I-1) * 16807, 2147483647)) RETURN
      END DO
C
      IF (C1HEAP(0) .NE. CHAR(78)) RETURN
      DO I = 1, 15
        J = ICHAR(C1HEAP(I-1)) - 65
        IF (C1HEAP(I) .NE. CHAR(MOD(J * 5 + 3, 16) + 65)) RETURN
      END DO
C
      ERROR = .FALSE.
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE INITHP
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
C
C heap initialization routine
C
C Authors: Axel T. Brunger & Ralf W. Grosse-Kunstleve
C
C local
      INTEGER  I, J
      LOGICAL  ERROR
C
C begin
C initialize variables for statistics
      CALL CNS0ALLOC
C
C initialize static heap arrays with pseudo-random numbers
      HEAP(0) = 13
      DO I = 1, 15
        HEAP(I) = MOD(HEAP(I-1) * 16807, 2147483647)
      END DO
C
      C1HEAP(0) = CHAR(78)
      DO I = 1, 15
        J = ICHAR(C1HEAP(I-1)) - 65
        C1HEAP(I) = CHAR(MOD(J * 5 + 3, 16) + 65)
      END DO
C
      CALL CHKSHPA(ERROR)
      IF (ERROR) THEN
        WRITE(6, '(1X, A)') 'INITHP: Fatal coding error.'
        CALL DIE
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE PRINHP
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
C
C prints info about the HEAP and checks static heap arrays
C
C Author: Axel T. Brunger & Ralf W. Grosse-Kunstleve
C
C local
      LOGICAL       ERROR
      CHARACTER*11  CURBYTES, MAXBYTES
      CHARACTER*11  CUROVERH, MAXOVERH
C
C externals
      LOGICAL   IS0ALLOC
      EXTERNAL  IS0ALLOC
C
C begin
C
      IF (WRNLEV .GE. 15 .OR. .NOT. IS0ALLOC()) THEN
        CALL QRYALLOC( 2, CURBYTES)
        CALL QRYALLOC( 3, MAXBYTES)
        CALL QRYALLOC(-2, CUROVERH)
        CALL QRYALLOC(-3, MAXOVERH)
        WRITE(6, '(5A)')
     &  ' HEAP: maximum use      = ', MAXBYTES,
     &  ' current use      = ',       CURBYTES,
     &  ' bytes'
        WRITE(6, '(5A)')
     &  ' HEAP: maximum overhead = ', MAXOVERH,
     &  ' current overhead = ',       CUROVERH,
     &  ' bytes'
      END IF
C
      CALL CHKSHPA(ERROR)
      IF (ERROR) THEN
        CALL WRNDIE(-5, 'PRINHP', 'static heap array(s) corrupted.')
        CALL DIE
      END IF
C
      RETURN
      END
C
C======================================================================
C
      INTEGER FUNCTION ALLHP(SIZE)
      IMPLICIT NONE
C I/O
      INTEGER  SIZE
C
C heap allocation routine
C
C Author: Axel T. Brunger & R.W. Grosse-Kunstleve
C
C local
      INTEGER  OFFS
C
C begin
      IF (SIZE .LT. 0) THEN
        CALL WRNDIE(-5, 'ALLHP', 'requested length is less than zero')
        CALL DIE
      END IF
C
      OFFS = 0
      CALL CHKALLOC('I', 0, 'ALLHP', -1, SIZE, OFFS)
C
      ALLHP = OFFS
      RETURN
      END
C
C======================================================================
C
      INTEGER FUNCTION REAHP(OLDOFFS, OLDSIZE, NEWSIZE)
      IMPLICIT NONE
C I/O
      INTEGER  OLDOFFS, OLDSIZE, NEWSIZE
C
C heap re-allocation routine
C
C Author: R.W. Grosse-Kunstleve
C
C begin
      IF (NEWSIZE .LT. 0) THEN
        CALL WRNDIE(-5, 'REAHP', 'requested length is less than zero')
        CALL DIE
      END IF
C
      REAHP = OLDOFFS
      IF (OLDOFFS .EQ. 0 .AND. OLDSIZE .EQ. 0) THEN
        CALL CHKALLOC('I', 0, 'REAHP',      -1, NEWSIZE, REAHP)
      ELSE
        CALL CHKALLOC('I', 0, 'REAHP', OLDSIZE, NEWSIZE, REAHP)
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE FREHP(OFFS, SIZE)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER  OFFS, SIZE
C
C heap de-allocation routine
C
C Author: Axel T. Brunger & R.W. Grosse-Kunstleve
C
C begin
      CALL CHKALLOC('I', 0, 'FREHP', SIZE, -1, OFFS)
C
      RETURN
      END
C
C======================================================================
C
      INTEGER FUNCTION CALLHP(NSTR, STRLEN)
      IMPLICIT NONE
C I/O
      INTEGER  NSTR, STRLEN
C
C character heap allocation routine
C
C Author: R.W. Grosse-Kunstleve
C
C local
      INTEGER  OFFS
C
C begin
      IF (NSTR .LT. 0 .OR. STRLEN .LT. 0) THEN
        CALL WRNDIE(-5, 'CALLHP', 'requested length is less than zero')
        CALL DIE
      END IF
C
      OFFS = 0
      CALL CHKALLOC('C', STRLEN, 'CALLHP', -1, NSTR, OFFS)
C
      CALLHP = OFFS
      RETURN
      END
C
C======================================================================
C
      INTEGER FUNCTION CREAHP(OLDOFFS, OLDNSTR, NEWNSTR, STRLEN)
      IMPLICIT NONE
C I/O
      INTEGER  OLDOFFS, OLDNSTR, NEWNSTR, STRLEN
C
C character heap re-allocation routine
C
C Author: R.W. Grosse-Kunstleve
C
C begin
      IF (NEWNSTR .LT. 0 .OR. STRLEN .LT. 0) THEN
        CALL WRNDIE(-5, 'CREAHP', 'requested length is less than zero')
        CALL DIE
      END IF
C
      CREAHP = OLDOFFS
      IF (OLDOFFS .EQ. 0 .AND. OLDNSTR .EQ. 0) THEN
        CALL CHKALLOC('C', STRLEN, 'CREAHP',      -1, NEWNSTR, CREAHP)
      ELSE
        CALL CHKALLOC('C', STRLEN, 'CREAHP', OLDNSTR, NEWNSTR, CREAHP)
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE CFREHP(OFFS, NSTR, STRLEN)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER  OFFS, NSTR, STRLEN
C
C character heap de-allocation routine
C
C Author: R.W. Grosse-Kunstleve
C
C begin
      CALL CHKALLOC('C', STRLEN, 'CFREHP', NSTR, -1, OFFS)
C
      RETURN
      END
