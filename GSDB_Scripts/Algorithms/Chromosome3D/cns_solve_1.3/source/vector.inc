C
C vector.inc
C
C common block for VECTOR command parsing and execution.
C
C
C most vflc's or operators allowed in a vector expression; where
C vflc's are: variables, functions, or literal constants.
C parameter MAXOPR: maximum number of operators in expression
      INTEGER MAXOPR
      PARAMETER(MAXOPR=40)
      INTEGER VNAMLN
      PARAMETER(VNAMLN=WORD_SIZE)
      INTEGER SVARMX
      PARAMETER(SVARMX=WORD_SIZE)
      CHARACTER*1 TEO(MAXOPR)
      CHARACTER*(VNAMLN) VAREQ(MAXOPR)
      CHARACTER*2  VARTYP(MAXOPR)
      INTEGER VARNUM(MAXOPR)
      DOUBLE PRECISION RVAREQ(MAXOPR)
      INTEGER VOP(MAXOPR)
      INTEGER NUMOPS
      DOUBLE COMPLEX CVAREQ(MAXOPR)
      CHARACTER*(SVARMX) SVAREQ(MAXOPR)
      INTEGER SVARLN(MAXOPR)
      INTEGER EPAR
      INTEGER XSTMAX
C
      CHARACTER*(VNAMLN) VNAMES(80)
C
      CHARACTER*2 VTYP(80)
      INTEGER VNUM(79), VUNIQ(79), SCNUMV, XRNUMV
C     to perform the normalization function we need to store
C     the sum of the squares
      DOUBLE PRECISION NSQSUM
C
C stores: pointer to heap (PTRSTO) and length of store (LENSTO)
C =============================================================
      INTEGER MAXSTO
      PARAMETER (MAXSTO=9)
      INTEGER PTRSTO, LENSTO
C
C COMMON BLOCKS
C =============
C
C     double precision block
C
      COMMON /SCBLCK/  RVAREQ, NSQSUM
C
C     double complex block
C
      COMMON /SCBLCK/  CVAREQ
C
C     integer block
C
      COMMON /SCBLCK/ SCNUMV, VARNUM, VOP, NUMOPS, VNUM,
     &                VUNIQ,  XRNUMV, EPAR, XSTMAX, SVARLN,
     &                PTRSTO(MAXSTO), LENSTO(MAXSTO)
C
C character common block
C
      COMMON /SCBLC2/ SVAREQ, VNAMES, VTYP,  TEO, VAREQ, VARTYP
C
