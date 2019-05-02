C data structure COORDC
C =====================
C
C comparison cartesian coordinates, temperature factor and occupancies.
C
C MAXAT4 must match MAXA
C
      INTEGER MAXAT4
      PARAMETER (MAXAT4=MAX_ATOM)
      DOUBLE PRECISION XCOMP(MAXAT4), YCOMP(MAXAT4), ZCOMP(MAXAT4)
      DOUBLE PRECISION WCOMP(MAXAT4), QCOMP(MAXAT4)
      COMMON/COORDC/ XCOMP, YCOMP, ZCOMP, WCOMP, QCOMP
C
      SAVE /COORDC/
