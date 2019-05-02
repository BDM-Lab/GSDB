C funct.inc
C
C     function declarations:
C
C==============================
C     declared in genic.s
C
      INTEGER  PATOM
      EXTERNAL PATOM
C
C==============================
C     declared in array.s
C
      INTEGER  FIND52,SRCHWD
      EXTERNAL FIND52,SRCHWD
      INTEGER  SRCHC4
      EXTERNAL SRCHC4
C
C==============================
C     declared in util.s
C
      INTEGER  GETATN
      EXTERNAL GETATN
      LOGICAL  HYDROG,INITIA
      EXTERNAL HYDROG,INITIA
C
C==============================
C     declared in sort.s
C
      LOGICAL  ORDER5,ORDER,ORDERR
      EXTERNAL ORDER5,ORDER,ORDERR
C
C==============================
C     declared in string.s
C
      LOGICAL  LTSTEQ
      EXTERNAL LTSTEQ
      INTEGER  DECODI
      EXTERNAL DECODI
      DOUBLE PRECISION DECODF
      EXTERNAL DECODF
      INTEGER RINDEX
      EXTERNAL RINDEX
C
C==============================
C     declared in cstack.s
C
      INTEGER  CALLST
      EXTERNAL CALLST
C
C==============================
C     declared in machine_f.s
C
      INTEGER  IREAL8,IREAL4,INTEG4,ICPLX4,ICPLX8
      EXTERNAL IREAL8,IREAL4,INTEG4,ICPLX4,ICPLX8
      INTEGER  ICHAR4,ILOGIC
      EXTERNAL ICHAR4,ILOGIC
C
C==============================
C     declared in heap.s
C
      INTEGER  ALLHP, REAHP, CALLHP, CREAHP
      EXTERNAL ALLHP, REAHP, CALLHP, CREAHP
C
