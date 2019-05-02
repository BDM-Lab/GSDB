C data structure XFFT
C ===================
C
C parameters for Fast Fourier computation of structure factors
C and derivatives.
C
C grid size relative to high-resolution
      DOUBLE PRECISION GRID
C required grid factors
      INTEGER XGRIDF, YGRIDF, ZGRIDF
C maximum memory allocation for electron density matrix
      INTEGER MEMORY
C automatically allocate memory for electron density matrix
      LOGICAL QMEMAUTO
C default memory allocation for AUTOmatic mode
      INTEGER DEFMEMA
C temperature factor scale to minimize aliasing
      DOUBLE PRECISION BSCAL
C ratio of form atomic form factor at maximum and at RCUT
      DOUBLE PRECISION ELIM
C volume of grid box, a*, b*, c*
      DOUBLE PRECISION FTVOL, FTAS, FTBS, FTCS
C
C FFT restrictions
      INTEGER FTBASE, FTCRBX, FTPRIM, FTAVOI
C
C parameter
      PARAMETER (DEFMEMA=2000000)
C
C     double precision block
C
      COMMON /COXFFT/ GRID,  BSCAL, ELIM,
     &                FTVOL, FTAS, FTBS, FTCS
C
C     integer block
C
      COMMON /IOXFFT/ XGRIDF, YGRIDF, ZGRIDF,
     &                MEMORY, FTBASE, FTCRBX, FTPRIM, FTAVOI
C
C     logical block
C
      COMMON /LOXFFT/ QMEMAUTO
C
      SAVE /COXFFT/
      SAVE /IOXFFT/
      SAVE /LOXFFT/
