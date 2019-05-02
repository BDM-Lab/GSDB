      SUBROUTINE HVFYSS(STR, LSTR, NSTR)
      IMPLICIT NONE
C I/O
      CHARACTER  STR(*)*(*)
      INTEGER    LSTR, NSTR
C
C Fill array of strings with characters A-Z.
C
C local
      INTEGER  I, J, K
C
C begin
      IF (LEN(STR(1)) .NE. LSTR) THEN
        CALL WRNDIE(-5, 'HVFYSS', 'LEN(STR(1)) .NE. LSTR')
        CALL DIE
      END IF
      K = 65
      DO I = 1, NSTR
        DO J = 1, LSTR
          STR(I)(J:J) = CHAR(K)
          K = K + 1
          IF (K .EQ. 91) K = 65
        END DO
      END DO
C
      RETURN
      END
C
C=======================================================================
C
      LOGICAL FUNCTION HVFYCS(STR, NSTR)
      IMPLICIT NONE
C I/O
      CHARACTER  STR(*)*(*)
      INTEGER    NSTR
C
C Check array of strings filled before with HVFYSS.
C
C local
      INTEGER  I, J, K
C
C begin
      HVFYCS = .TRUE.
      K = 65
      DO I = 1, NSTR
        DO J = 1, LEN(STR(I))
          IF (STR(I)(J:J) .NE. CHAR(K)) RETURN
          K = K + 1
          IF (K .EQ. 91) K = 65
        END DO
      END DO
      HVFYCS = .FALSE.
C
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE HEAPVFY(ERROR)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      LOGICAL  ERROR
C
C Exercise heap allocation routines.
C
C Author: R.W. Grosse-Kunstleve
C
C local
      INTEGER   NTESTS
      PARAMETER(NTESTS=32)
      INTEGER    N1, N2, N3, I
      INTEGER    OFFSI(NTESTS)
      INTEGER    OFFS0(NTESTS)
      CHARACTER  C5HEAP(0:1)*5
      EQUIVALENCE(C1HEAP(0), C5HEAP(0))
C
C externals
      LOGICAL   HVFYCS, IS0ALLOC
      EXTERNAL  HVFYCS, IS0ALLOC
C
C begin
      ERROR = .TRUE.
C
      N1 = 41
      N2 = 53
      N3 = 61
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 0'
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 1'
      OFFSI(1) = ALLHP(INTEG4(N1))
      CALL FILL4(HEAP(OFFSI(1)), N1, 0)
      CALL FREHP(OFFSI(1), INTEG4(N1))
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 2'
      OFFSI(1) = ALLHP(INTEG4(N1))
      CALL FILL4(HEAP(OFFSI(1)), N1, 0)
      OFFSI(1) = REAHP(OFFSI(1), INTEG4(N1), INTEG4(N2))
      CALL FILL4(HEAP(OFFSI(1)), N2, 0)
      CALL FREHP(OFFSI(1), INTEG4(N2))
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 3'
      OFFSI(1) = CALLHP(N1, 1)
      CALL HVFYSS(C1HEAP(OFFSI(1)), 1, N1)
      IF (HVFYCS(C1HEAP(OFFSI(1)), N1)) RETURN
      CALL CFREHP(OFFSI(1), N1, 1)
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 4'
      OFFSI(1) = CALLHP(N2, 4)
      CALL HVFYSS(C4HEAP(OFFSI(1)), 4, N2)
      IF (HVFYCS(C4HEAP(OFFSI(1)), N2)) RETURN
      CALL CFREHP(OFFSI(1), N2, 4)
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 5'
      OFFSI(1) = CALLHP(N3, 80)
      CALL HVFYSS(C80HEAP(OFFSI(1)), 80, N3)
      IF (HVFYCS(C80HEAP(OFFSI(1)), N3)) RETURN
      CALL CFREHP(OFFSI(1), N3, 80)
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 6'
      DO I = 1, NTESTS
        OFFSI(I) = CALLHP(N1 + I, 5)
        CALL HVFYSS(C5HEAP(OFFSI(I)), 5, N1 + I)
        IF (HVFYCS(C5HEAP(OFFSI(I)), N1 + I)) RETURN
        CALL CFREHP(OFFSI(I), N1 + I, 5)
        OFFS0(1) = ALLHP(0)
        CALL FREHP(OFFS0(1), 0)
        IF (.NOT. IS0ALLOC()) RETURN
      END DO
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 7'
      DO I = 1, NTESTS
        OFFSI(I) = CALLHP(N1 + I, 5)
        CALL HVFYSS(C5HEAP(OFFSI(I)), 5, N1 + I)
        IF (HVFYCS(C5HEAP(OFFSI(I)), N1 + I)) RETURN
        OFFS0(I) = ALLHP(0)
      END DO
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 8'
      DO I = 1, NTESTS
        OFFSI(I) = CREAHP(OFFSI(I), N1 + I, N3 + I, 5)
        IF (HVFYCS(C5HEAP(OFFSI(I)), MIN(N1, N3) + I)) RETURN
        CALL HVFYSS(C5HEAP(OFFSI(I)), 5, N3 + I)
        IF (HVFYCS(C5HEAP(OFFSI(I)), N3 + I)) RETURN
      END DO
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 9'
      DO I = NTESTS, 1, -1
        OFFSI(I) = CREAHP(OFFSI(I), N3 + I, N2 + I, 5)
        IF (HVFYCS(C5HEAP(OFFSI(I)), MIN(N3, N2) + I)) RETURN
        CALL HVFYSS(C5HEAP(OFFSI(I)), 5, N2 + I)
        IF (HVFYCS(C5HEAP(OFFSI(I)), N2 + I)) RETURN
      END DO
C
      IF (WRNLEV .GE. 10) WRITE(6, '(1X,A)') 'HEAP EXERCISE 10'
      DO I = 1, NTESTS
        CALL CFREHP(OFFSI(I), N2 + I, 5)
        CALL FREHP(OFFS0(I), 0)
      END DO
      IF (.NOT. IS0ALLOC()) RETURN
      OFFS0(1) = ALLHP(0)
      CALL FREHP(OFFS0(1), 0)
      IF (.NOT. IS0ALLOC()) RETURN
C
      ERROR = .FALSE.
C
      RETURN
      END
