C learn.inc
C
C
C data structure for parameter learning facility
C
      INTEGER PRIND, LPRIN, PNFRAM, PRSELCT
      LOGICAL PRINIT, PQCONST
C
C     logical block
C
      COMMON/QLEARN/ PRINIT, PQCONST
C
C     integer block
C
      COMMON/ILEARN/ PRIND, LPRIN, PNFRAM, PRSELCT
C
      SAVE /QLEARN/
      SAVE /ILEARN/
