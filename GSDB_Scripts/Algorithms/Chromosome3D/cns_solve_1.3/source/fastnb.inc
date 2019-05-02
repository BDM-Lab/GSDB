C fastnb.inc
C ==========
C Fast nonbonded interactions
C
C Paul Adams 11/4/94
C
C grid based pairlist variables
      LOGICAL QFSTNB
      INTEGER MXBOX, MNBOX
      PARAMETER (MXBOX=2000,MNBOX=27)
      DOUBLE PRECISION RATMIN
      PARAMETER (RATMIN=2.5D0)
C
C PMTA variables
      CHARACTER*4 PMTAM
      DOUBLE PRECISION PTHETA
      INTEGER PLEVEL,PTERMS
      LOGICAL QPFFT
C
C ====================================================================
C
C     logical block
C
      COMMON /LFNB/ QFSTNB,QPFFT
C
C     double precision block
C
      COMMON /RFNB/ PTHETA
C
C     integer block
C
      COMMON /IFNB/ PLEVEL,PTERMS
C
C     character block
C
      COMMON /CFNB/ PMTAM
C
      SAVE /LFNB/
      SAVE /RFNB/
      SAVE /IFNB/
      SAVE /CFNB/
C
