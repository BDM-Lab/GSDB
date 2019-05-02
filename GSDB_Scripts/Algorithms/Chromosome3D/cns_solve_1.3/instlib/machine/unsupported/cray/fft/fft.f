C ======================================================================
C FFT routines ========================================================
C =====================================================================
      SUBROUTINE FFTPRP(BASE,PRIME,AVOID)
C
C Returns maximum factor PRIME allowed in integer
C dimensions of 3-d FFT.  Also returns BASE which specifies
C that BASE should be a factor, e.g. if BASE is two then
C the integer has to be even.
C
C ******************************************************
C UNICOS/COS version
C ******************************************************
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'xfft.inc'
      INTEGER BASE, PRIME, AVOID
C begin
      BASE=1
      PRIME=5
      AVOID=2
C
      FTCRBX=1
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE FFT3C(MX, MY, NX, NY, NZ, A, ERROR)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER         MX, MY, NX, NY, NZ
      DOUBLE COMPLEX  A(MX, MY, *)
      LOGICAL         ERROR
C
C 3D complex-to-complex transform
C
C local
      INTEGER  MTABLE
      INTEGER  MWORK
      INTEGER  ISYS(0:1)
C pointer
      INTEGER  TABLE
      INTEGER  WORK
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'FFT3C: Using Cray Scientific Library'
C
      ERROR = .FALSE.
C
      MTABLE = 100 + 2 * (NX + NY + NZ)
      TABLE = ALLHP(IREAL8(MTABLE))
C
      MWORK = 512 * MAX(NX, NY, NZ)
      WORK = ALLHP(IREAL8(MWORK))
C
      ISYS(0) = 0
C
      CALL CCFFT3D(0, NX, NY, NZ, 1.0, A, MX, MY, A, MX, MY,
     &             HEAP(TABLE), HEAP(WORK), ISYS)
      CALL CCFFT3D(1, NX, NY, NZ, 1.0, A, MX, MY, A, MX, MY,
     &             HEAP(TABLE), HEAP(WORK), ISYS)
C
      CALL FREHP(WORK, IREAL8(MWORK))
      CALL FREHP(TABLE, IREAL8(MTABLE))
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SFFT1C(N, A, ERROR, MWSP, WSP)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER  N
      COMPLEX           A(*)
      LOGICAL  ERROR
      INTEGER  MWSP
      INTEGER   WSP
C
C 1D complex-to-complex transform, single precision
C
C local
      INTEGER  MTABLE
      INTEGER  MWORK
      INTEGER  ISYS(0:1)
      LOGICAL  DOINIT
C pointer
      INTEGER  TABLE
      INTEGER  WORK
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1C: Using Cray Scientific Library'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        MTABLE = 100 + 8 * N
        MWORK = 8 * N
        ISYS(0) = 0
        DOINIT = .FALSE.
C
        IF (MWSP .EQ. 0) THEN
          MWSP = MTABLE + MWORK
           WSP = ALLHP(IREAL4(MWSP))
          DOINIT = .TRUE.
        END IF
C
        TABLE = WSP
        WORK  = WSP + MTABLE
C
        IF (DOINIT) THEN
          CALL CCFFT(0, N, 1.0, A, A,
     &               HEAP(TABLE), HEAP(WORK), ISYS)
        END IF
C
        CALL CCFFT(1, N, 1.0, A, A,
     &             HEAP(TABLE), HEAP(WORK), ISYS)
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, IREAL4(MWSP))
        MWSP = 0
         WSP = 0
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SFFT1CR(N, A, ERROR, MWSP, WSP)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER  N
      REAL              A(*)
      LOGICAL  ERROR
      INTEGER  MWSP
      INTEGER   WSP
C
C 1D complex-to-real transform, single precision
C
C local
      INTEGER  MTABLE
      INTEGER  MWORK
      INTEGER  ISYS(0:1)
      LOGICAL  DOINIT
C pointer
      INTEGER  TABLE
      INTEGER  WORK
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1CR: Using Cray Scientific Library'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        MTABLE = 100 + 4 * N
        MWORK = 4 + 4 * N
        ISYS(0) = 0
        DOINIT = .FALSE.
C
        IF (MWSP .EQ. 0) THEN
          MWSP = MTABLE + MWORK
           WSP = ALLHP(IREAL4(MWSP))
          DOINIT = .TRUE.
        END IF
C
        TABLE = WSP
        WORK  = WSP + MTABLE
C
        IF (DOINIT) THEN
          CALL CSFFT(0, N, 1.0, A, A,
     &               HEAP(TABLE), HEAP(WORK), ISYS)
        END IF
C
        CALL CSFFT(1, N, 1.0, A, A,
     &             HEAP(TABLE), HEAP(WORK), ISYS)
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, IREAL4(MWSP))
        MWSP = 0
         WSP = 0
      END IF
C
      RETURN
      END
