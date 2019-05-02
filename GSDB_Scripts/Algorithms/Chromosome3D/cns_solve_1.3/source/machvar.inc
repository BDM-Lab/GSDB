C
C machvar.inc
C ===========
C
C Machine dependent variables, obtained empirically on startup.
C
C temporary variable
      DOUBLE PRECISION  DPTMPV
C
C machine epsilon
      DOUBLE PRECISION  FPEPS
C
C limit of epsilon (eps = 2^MXFPEPS2)
      INTEGER MXFPEPS2
      PARAMETER (MXFPEPS2=1024)
C
      COMMON /MACHDP/   FPEPS, DPTMPV
      SAVE   /MACHDP/
