C lbfgs.inc
C =========
C
C common parameters for Limited Memory BFGS minimizer
C
C Paul Adams and Axel Brunger
C
      DOUBLE PRECISION GTOL, STPMIN, STPMAX
C
C variables needed for reverse calling
      INTEGER INFOC
C
      LOGICAL BRACKT,STAGE1
C
      DOUBLE PRECISION DG, DGINIT, DGTEST, DGX, DGY
      DOUBLE PRECISION FINIT, FTEST1, FX, FY, STX, STY
      DOUBLE PRECISION STMIN, STMAX, WIDTH, WIDTH1
C
      COMMON /LBFGSL/ BRACKT, STAGE1
C
      COMMON /LBFGSI/ INFOC
C
      COMMON /LBFGSD/ GTOL, STPMIN, STPMAX,
     &                DG, DGINIT, DGTEST, DGX, DGY,
     &                FINIT, FTEST1, FX, FY, STX, STY,
     &                STMIN, STMAX, WIDTH, WIDTH1
C
      SAVE /LBFGSL/
      SAVE /LBFGSI/
      SAVE /LBFGSD/
C
