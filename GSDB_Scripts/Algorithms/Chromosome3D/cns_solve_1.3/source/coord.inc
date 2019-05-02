C data structure COORD
C ====================
C
C main cartesian coordinates, temperature factor and occupancies.
C
C MAXAT3 must match MAXA
C
      INTEGER MAXAT3
      PARAMETER (MAXAT3=MAX_ATOM)
      DOUBLE PRECISION X(MAXAT3), Y(MAXAT3), Z(MAXAT3)
      DOUBLE PRECISION WMAIN(MAXAT3), QMAIN(MAXAT3)
      COMMON/COORD/ X, Y, Z, WMAIN, QMAIN
C
      SAVE /COORD/
