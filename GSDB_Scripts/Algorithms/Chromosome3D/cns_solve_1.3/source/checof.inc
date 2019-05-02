C checof.inc
C ==========
C
C Chebychev coefficients for maximum likelihood estimation
C
C R.J. Read and N.S. Pannu
C
C Modified Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      INTEGER MXCHVL
      PARAMETER (MXCHVL=50)
C
      INTEGER NCHC(8)
C
      DOUBLE PRECISION WCHC(8)
C
      DOUBLE PRECISION CHCO(MXCHVL,8)
C
      COMMON /ICHE/ NCHC
      COMMON /RCHE/ WCHC, CHCO
C
      SAVE /ICHE/
      SAVE /RCHE/
C
