C Maximum number of maps in mapyard
      INTEGER   MMAPYARD
      PARAMETER(MMAPYARD=10)
C Labels of maps
      CHARACTER  LBMAPYARD(MMAPYARD)*16
C Dimensions of maps
      INTEGER    RNMAPYARD(3, MMAPYARD)
C Pointers for maps
      INTEGER    PTMAPYARD(MMAPYARD)
C
      COMMON /MYCHR/ LBMAPYARD
      COMMON /MYINT/ RNMAPYARD
      COMMON /MYPTR/ PTMAPYARD
C
      SAVE /MYCHR/
      SAVE /MYINT/
      SAVE /MYPTR/
