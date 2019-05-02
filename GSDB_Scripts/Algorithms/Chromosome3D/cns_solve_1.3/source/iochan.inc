C iochan.inc
C ==========
C
C holds machine dependent information about FORTRAN units 1...99
C and their associated files
C
C
C list for free FORTRAN units
C see routine FREEUN
      INTEGER MAXUN
      PARAMETER (MAXUN=99)
      INTEGER IFREEU(MAXUN)
C
C for CRAY version:
C common
      CHARACTER*60 COSPDN(MAXUN)
      CHARACTER*8  COSDN(MAXUN)
C
C-COMMON-declarations
      COMMON /IOCHAN/ IFREEU
      COMMON /CIOCHA/ COSPDN, COSDN
      SAVE /IOCHAN/
      SAVE /CIOCHA/
