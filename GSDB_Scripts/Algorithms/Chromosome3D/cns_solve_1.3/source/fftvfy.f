C{
      LOGICAL FUNCTION VFYGRF(GRIDF, MAXPRIME)
      IMPLICIT NONE
C I/O
      INTEGER  GRIDF, MAXPRIME
C
C Check if GridF contains a prime which is greater than MaxPrime
C
C begin
      INTEGER  N, P
C
      VFYGRF = .TRUE.
C
      N = GRIDF
C
      IF (N .LE. 1) RETURN
C
      IF (MAXPRIME .LT. 2) THEN
        VFYGRF = .FALSE.
        RETURN
      END IF
C
      DO WHILE (MOD(N, 2) .EQ. 0)
        N = N / 2
      END DO
      IF (N .EQ. 1) RETURN
C
      DO P = 3, MAXPRIME, 2
        DO WHILE (MOD(N, P) .EQ. 0)
          N = N / P
        END DO
        IF (N .EQ. 1) RETURN
      END DO
C
      VFYGRF = .FALSE.
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE VFYAGF(XGRIDF, YGRIDF, ZGRIDF,
     &                  XGF, YGF, ZGF,
     &                  FTPRIM)
      IMPLICIT NONE
C I/O
      INTEGER  XGRIDF, YGRIDF, ZGRIDF
      INTEGER  XGF, YGF, ZGF
      INTEGER  FTPRIM
C
C Check if [XYZ]GridF contains a prime which is greater than MaxPrime.
C Since WRNDIE does not always die, return "corrected" factors in [XYZ]GF.
C
C external
      LOGICAL   VFYGRF
      EXTERNAL  VFYGRF
C
C begin
      XGF = XGRIDF
      YGF = YGRIDF
      ZGF = ZGRIDF
C
      IF (.NOT. VFYGRF(XGF, FTPRIM)) THEN
        CALL WRNDIE(-5,'VFYAGF',
     &    'XGRIdfactor not consistent with PRIMe')
        XGF = 1
      END IF
C
      IF (.NOT. VFYGRF(YGF, FTPRIM)) THEN
        CALL WRNDIE(-5,'VFYAGF',
     &    'YGRIdfactor not consistent with PRIMe')
        YGF = 1
      END IF
C
      IF (.NOT. VFYGRF(ZGF, FTPRIM)) THEN
        CALL WRNDIE(-5,'VFYAGF',
     &    'ZGRIdfactor not consistent with PRIMe')
        ZGF = 1
      END IF
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE ICMAP(LX, LY, NX, NY, NZ, CMAP)
      IMPLICIT NONE
C I/O
      INTEGER         LX, LY, NX, NY, NZ
      DOUBLE COMPLEX  CMAP(LX, LY, *)
C
C Initialize double complex map for test of FFT library.
C
C local
      INTEGER  IX, IY, IZ
C
C begin
      DO IZ = 1, NZ
      DO IY = 1, NY
      DO IX = 1, NX
        CMAP(IX, IY, IZ) = CMPLX(DFLOAT(MOD(IX + IY + IZ, 11)),
     &                           DFLOAT(MOD(IX + IY + IZ, 13)))
      END DO
      END DO
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE CJGMAP(LX, LY, NX, NY, NZ, CMAP)
      IMPLICIT NONE
C I/O
      INTEGER         LX, LY, NX, NY, NZ
      DOUBLE COMPLEX  CMAP(LX, LY, *)
C
C Set map elements to complex-conjugate (b/o the way FFT3C works).
C
C local
      INTEGER  IX, IY, IZ
C
C begin
      DO IZ = 1, NZ
      DO IY = 1, NY
      DO IX = 1, NX
        CMAP(IX, IY, IZ) =   DCONJG(CMAP(IX, IY, IZ))
     &                     / DFLOAT(NX * NY * NZ)
      END DO
      END DO
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE VFYCMAP(LX, LY, NX, NY, NZ, CMAP, ERROR)
      IMPLICIT NONE
C I/O
      INTEGER         LX, LY, NX, NY, NZ
      DOUBLE COMPLEX  CMAP(LX, LY, *)
      LOGICAL         ERROR
C
C See if we got the right numbers back after two FFT's.
C
C local
      INTEGER  IX, IY, IZ
C
C begin
      ERROR = .TRUE.
C
      DO IZ = 1, NZ
      DO IY = 1, NY
      DO IX = 1, NX
        IF ( NINT(DBLE(CMAP(IX,IY,IZ))) .NE. MOD(IX+IY+IZ, 11)) RETURN
        IF (-NINT(DIMAG(CMAP(IX,IY,IZ))) .NE. MOD(IX+IY+IZ, 13)) RETURN
      END DO
      END DO
      END DO
C
      ERROR = .FALSE.
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE ISSEQ(N, SSEQ)
      IMPLICIT NONE
C I/O
      INTEGER  N
      COMPLEX  SSEQ(0:*)
C
C Initialize complex sequence for test of FFT library.
C
C local
      INTEGER  I
C
C begin
      DO I = 0, N - 1
        SSEQ(I) = CMPLX(DFLOAT(MOD(I + 4, 7)), DFLOAT(0))
      END DO
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE VFYSSEQ(N, SSEQ, ERROR)
      IMPLICIT NONE
C I/O
      INTEGER  N
      REAL     SSEQ(0:*)
      LOGICAL  ERROR
C
C See if we got the right numbers back after two FFT's.
C (SSEQ is transposed and multiplied by N.)
C
C local
      INTEGER  I
C
C begin
      ERROR = .TRUE.
C
      DO I = 0, N - 1
        IF (NINT(SSEQ(MOD(N - I, N)) / DFLOAT(N)) .NE. MOD(I + 4, 7))
     &    RETURN
      END DO
C
      ERROR = .FALSE.
C
      RETURN
      END
C}
C=======================================================================
C{
      SUBROUTINE VFYFFT
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'xfft.inc'
      INCLUDE 'timer.inc'
C
C Make sure FFT library actually works.
C
C local
      INTEGER  WLSAVE, I
      LOGICAL  BAILOUT, ERROR
      INTEGER  BASE, PRIME, AVOID
      INTEGER  NN(3), CM(3), PCM
      INTEGER         RM(1)
      INTEGER   MFFTW
      INTEGER  HPFFTW
C
C pointer
      INTEGER  HPCMAP, HPSSEQ
C
C external
      INTEGER   ILCM
      EXTERNAL  ILCM
      LOGICAL   VFYGRF
      EXTERNAL  VFYGRF
C
C begin
      CALL FFTPRP(BASE, PRIME, AVOID)
C
C Check if BASE and PRIME are consistent
      IF (     PRIME .LT. 2
     &    .OR. .NOT. VFYGRF(BASE,   PRIME)
     &    .OR. .NOT. VFYGRF(FTCRBx, PRIME)) THEN
        WRITE(6, '(1X, 2A)')
     &    'Fatal Coding Error: BASE & PRIME as returned by FFTPRP',
     &    ' are not consistent.'
        CALL DIE
      END IF
C
      BAILOUT = .FALSE.
C
C Set up double precision complex-to-complex transform
      NN(1) = 4*3
      NN(2) = 3*5
      NN(3) = 5*2
C
      DO I = 1, 3
        CALL XPRIME(NN(I), 1, BASE, PRIME)
        CM(I) = NN(I)
      END DO
C
      IF (AVOID .GT. 1) THEN
        IF (MOD(CM(1), AVOID) .EQ. 0) CM(1) = CM(1) + 1
        IF (MOD(CM(2), AVOID) .EQ. 0) CM(2) = CM(2) + 1
      END IF
C
      PCM = CM(1) * CM(2) * CM(3)
      HPCMAP = ALLHP(ICPLX8(PCM))
C
C FFT3C -> conjugate complex -> FFT3C -> VFYCMAP
      CALL  ICMAP(CM(1), CM(2), NN(1), NN(2), NN(3), HEAP(HPCMAP))
      WLSAVE = WRNLEV
      WRNLEV = MAX(11, WLSAVE)
      CALL  FFT3C(CM(1), CM(2), NN(1), NN(2), NN(3), HEAP(HPCMAP),
     &            ERROR)
      WRNLEV = WLSAVE
      IF (ERROR) THEN
        WRITE(6, '(1X, A)') 'Fatal Error: Call FFT3C 1 failed.'
        CALL FAQDIE
      END IF
      CALL CJGMAP(CM(1), CM(2), NN(1), NN(2), NN(3), HEAP(HPCMAP))
      CALL  FFT3C(CM(1), CM(2), NN(1), NN(2), NN(3), HEAP(HPCMAP),
     &            ERROR)
      IF (ERROR) THEN
        WRITE(6, '(1X, A)') 'Fatal Error: Call FFT3C 2 failed.'
        CALL FAQDIE
      END IF
C
C See if we got the right numbers back
      CALL VFYCMAP(CM(1), CM(2), NN(1), NN(2), NN(3), HEAP(HPCMAP),
     &             ERROR)
      IF (ERROR) THEN
        WRITE(6, '(1X, A)') 'Fatal Error: Test of FFT3C failed.'
        BAILOUT = .TRUE.
      END IF
C
      CALL FREHP(HPCMAP, ICPLX8(PCM))
C
C Test single precision complex-to-complex-to-real transform
C Try all N with 2 <= N <= 256 which are compatible with BASE and PRIME
C
      IF (FTCRBX .GT. 1) BASE = ILCM(BASE, FTCRBX)
C
      NN(1) = 1
      DO WHILE (NN(1) .LE. 256)
        NN(1) = NN(1) + 1
        CALL XPRIME(NN(1), 1, BASE, PRIME)
        CALL MFFTCR(1, NN, FTAvoi, RM)
        HPSSEQ = ALLHP(ICPLX4(NN(1)))
C
C SFFT1C -> SFFT1CR (with unique half) -> VFYSSEQ
        CALL  ISSEQ(NN(1), HEAP(HPSSEQ))
         MFFTW = 0
        HPFFTW = 0
        CALL SFFT1C(NN(1), HEAP(HPSSEQ), ERROR, MFFTW, HPFFTW)
        IF (ERROR) THEN
          WRITE(6, '(1X, A)') 'Fatal Error: CALL SFFT1C 1 failed.'
          CALL FAQDIE
        END IF
        CALL SFFT1C(    0,            0, ERROR, MFFTW, HPFFTW)
        IF (ERROR) THEN
          WRITE(6, '(1X, A)') 'Fatal Error: CALL SFFT1C 2 failed.'
          CALL FAQDIE
        END IF
C
        HPSSEQ = REAHP(HPSSEQ, ICPLX4(NN(1)), IREAL4(RM(1)))
C
        CALL SFFT1CR(NN(1), HEAP(HPSSEQ), ERROR, MFFTW, HPFFTW)
        IF (ERROR) THEN
          WRITE(6, '(1X, A)') 'Fatal Error: CALL SFFT1CR 1 failed.'
          CALL FAQDIE
        END IF
        CALL SFFT1CR(    0,            0, ERROR, MFFTW, HPFFTW)
        IF (ERROR) THEN
          WRITE(6, '(1X, A)') 'Fatal Error: CALL SFFT1CR 2 failed.'
          CALL FAQDIE
        END IF
        CALL VFYSSEQ(NN(1), HEAP(HPSSEQ), ERROR)
        IF (ERROR) THEN
          WRITE(6, '(1X, 2A, I3)')
     &      'Fatal Error: Test of SFFT1C & SFFT1CR failed for',
     &      ' N = ', NN(1)
          BAILOUT = .TRUE.
        END IF
C
        CALL FREHP(HPSSEQ, IREAL4(RM(1)))
      END DO
C
      IF (BAILOUT) THEN
         CALL FAQDIE
      END IF
C
      RETURN
      END
C}
C{
      SUBROUTINE FAQDIE
C
      IMPLICIT NONE
C
      WRITE(6, '(2A)')
     & 'Check the FAQ in the online documentation',
     & ' to see how to fix this problem'
      CALL DIE
C
      RETURN
      END
C}
