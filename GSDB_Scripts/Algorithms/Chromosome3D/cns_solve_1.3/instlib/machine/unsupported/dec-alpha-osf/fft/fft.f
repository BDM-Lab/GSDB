C======================================================================
C FFT routines ========================================================
C =====================================================================
      SUBROUTINE FFTPRP(BASE,PRIME,AVOID)
C
C Returns maximum factor PRIME allowed in integer
C dimensions of 3-d FFT.  Also returns BASE which specifies
C that BASE should be a factor, e.g. if BASE is two then
C the integer has to be even.
C
C ********************************************************
C Universal Version
C ********************************************************
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
      FTCRBX=2
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE FFT3C(MX,MY,NX,NY,NZ,A,ERROR)
C
C Computes 3-dimensional fast Fourier transform
C ************************************
C DXML version for complex arrays
C ************************************
C
C    X(x+1,y+1,z+1) = SUM {h=0...NX-1,k=0...NY-1,l=0...NZ-1}
C          A(h+1,k+1,l+1) *EXP(2*PI*i*(x*h/NX+y*k/NY+z*l/NZ))
C
C of a complex valued sequence A(x,y,z), where x=1,...,NX,
C y=1,...,NY, z=1,...,NZ/2.  The result X is returned in array
C A(x,y,z), where x=1,...,NX, y=1,...,NY, z=1,...,NZ.
C
C For conditions on NX, NY, see routine FFTPRP
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER MX, MY, NX, NY, NZ
      DOUBLE COMPLEX A(MX,MY,*)
      LOGICAL ERROR
C local
      INTEGER STATUS, I, J, K
C external
      INTEGER  ZFFT_3D
      EXTERNAL ZFFT_3D
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)')
     &    'FFT3C: Using Digital Extended Math Library (DXML)'
C
       ERROR = .FALSE.
C
       STATUS = ZFFT_3D ('C','C','B',A,A,NX,NY,NZ,MX,MY,1,1,1)
       IF (STATUS .NE. 0) THEN
         ERROR = .TRUE.
       ELSE
         DO K = 1,NZ
         DO J = 1,NY
         DO I = 1,NX
           A(I,J,K) = A(I,J,K) * DFLOAT(NX*NY*NZ)
         ENDDO
         ENDDO
         ENDDO
       END IF
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SFFT1C(N, A, ERROR, MWSP, WSP)
      IMPLICIT NONE
C I/O
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
      INTEGER STATUS, I
C
C external
      INTEGER  CFFT
      EXTERNAL CFFT
C
C begin
      ERROR = .FALSE.
C
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)')
     &    'SFFT1C: Using Digital Extended Math Library (DXML)'
C
      IF (N .GT. 0) THEN
        STATUS = CFFT('C', 'C', 'B', A, A, N, 1)
        IF (STATUS .NE. 0) THEN
          ERROR = .TRUE.
        ELSE
          DO I = 1, N
            A(I) = A(I) * FLOAT(N)
          END DO
        END IF
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
      INTEGER STATUS, I
C
C external
      INTEGER  SFFT
      EXTERNAL SFFT
C
C begin
      ERROR = .FALSE.
C
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)')
     &    'SFFT1CR: Using Digital Extended Math Library (DXML)'
C
      IF (N .GT. 0) THEN
        IF (MOD(N, 2) .NE. 0) THEN
          CALL WRNDIE(-5, 'SFFT1CR', 'Fatal Coding Error.')
          ERROR = .TRUE.
          RETURN
        END IF
C
        STATUS = SFFT('C', 'R', 'B', A, A, N, 1)
        IF (STATUS .NE. 0) THEN
          ERROR = .TRUE.
        ELSE
          DO I = 1, N
            A(I) = A(I) * FLOAT(N)
          END DO
        END IF
      END IF
C
      RETURN
      END
