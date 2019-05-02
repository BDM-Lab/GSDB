C**********************************************************************
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
      FTCRBX=1
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE FFT3C(MX,MY,NX,NY,NZ,A,ERROR)
C
C Computes 3-dimensional fast Fourier transform
C *******************
C Universal version for complex arrays
C *******************
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
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER MX, MY, NX, NY, NZ
      DOUBLE COMPLEX A(MX,MY,*)
      LOGICAL ERROR
C LOCAL
      INTEGER IX, IERR
C pointer
      INTEGER W1
C BEGIN
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'FFT3C: Using complib.sgimath'
C
      ERROR = .FALSE.
C
C WORK STORAGE ALLOCATION FOR FFT ROUTINE
      IX = NX+NY+NZ+3*15
      W1 = ALLHP(2*IREAL8(IX))
C
      CALL ZFFT3DI(NX,NY,NZ,HEAP(W1))
      CALL ZFFT3D(1,NX,NY,NZ,A,MX,MY,HEAP(W1))
      CALL ZFFT3DIFREE(NX,NY,NZ,HEAP(W1),IERR)
      IF (IERR .NE. 0) ERROR = .TRUE.
C
      CALL FREHP(W1,2*IREAL8(IX))
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
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1C: Using complib.sgimath'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        IF (MWSP .EQ. 0) THEN
          MWSP = N + 15
           WSP = ALLHP(ICPLX4(MWSP))
          CALL CFFT1DI(N, HEAP(WSP))
        END IF
        CALL CFFT1D(1, N, A, 1, HEAP(WSP))
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, ICPLX4(MWSP))
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
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1CR: Using complib.sgimath'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        IF (MWSP .EQ. 0) THEN
          MWSP = N + 15
          WSP = ALLHP(IREAL4(MWSP))
          CALL SFFT1DUI(N, HEAP(WSP))
        END IF
        CALL SFFT1DU(1, N, A, 1, HEAP(WSP))
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, IREAL4(MWSP))
        MWSP = 0
         WSP = 0
      END IF
C
      RETURN
      END
