C===============
C
C Anisotropy.inc
C
C by John Kuszewski July 1993
C===============
      INTEGER MAXSANICLASSES
      PARAMETER (MAXSANICLASSES = 10)
C
C arrays that hold anisotropy info
C anisassndx tells ending index of the anis_x arrays (below)
C for each class.
C anisforces holds k1 and anisswitches holds k2 for each class
C anisclasstypes holds the type of each class
C
C these are small enough not to bother with the heap.
C
      INTEGER SANIASSNDX (MAXSANICLASSES), SANISWITCHES(MAXSANICLASSES)
      REAL SANIFORCES (MAXSANICLASSES), SANICOEF0 (MAXSANICLASSES)
      REAL SANICOEF1 (MAXSANICLASSES), SANICOEF2 (MAXSANICLASSES)
      CHARACTER*8 SANICLASSNAMES (MAXSANICLASSES)
      LOGICAL PRINTCLASS(MAXSANICLASSES)
C
C maxanis = number of slots set aside for anisotropy assignments
C nanis = total number of anisotropy entered
C
      INTEGER MAXSANI, NSANI, NCLASSES, CURCLASS
C
C pointers to arrays to hold atom numbers, Anis-obs, and
C error for each aniostropy assignment
C
      INTEGER SANIIPTR, SANIJPTR, SANIKPTR, SANILPTR, SANIMPTR,
     &        SANINPTR,
     &        SANIOBSPTR,
     &        SANIERRPTR, CALCSANIPTR
C
C these names are for convenience
C
      INTEGER RELAX, MISC
      PARAMETER (RELAX = 1)
      PARAMETER (MISC = 2)
C
C input modes
C
      INTEGER MODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
       INTEGER POTENTIAL, SQUARE, HARMONIC
       PARAMETER (SQUARE = 1)
       PARAMETER (HARMONIC = 2)
C
C parameters as set up in ETOR
C
C===>parameter MCONST truncates the reciprocals and makes sure
C===>we don't devide by zero
      DOUBLE PRECISION MCONST
      PARAMETER (MCONST=0.0001D0)
C===>parameter EPS determines the absolute value of SIN(PHI) at
C===>which we want to switch to the alternative computation for the
C===>derivative
      DOUBLE PRECISION EPS
      PARAMETER (EPS=0.1D0)
C
C
C common blocks
C
      COMMON /CSANI/ SANICLASSNAMES
      COMMON /ISANI1/ SANIASSNDX, POTENTIAL,
     &       MAXSANI, NSANI, SANISWITCHES, CURCLASS, NCLASSES
      COMMON /ISANI2/ SANIIPTR, SANIJPTR, SANIKPTR, SANILPTR,
     &       SANIMPTR, SANINPTR,
     &       SANIOBSPTR, SANIERRPTR, CALCSANIPTR, MODE
      COMMON /RSANI/ SANIFORCES, SANICOEF0, SANICOEF1, SANICOEF2
      COMMON /LSANI/ PRINTCLASS
      SAVE /CSANI/
      SAVE /ISANI1/
      SAVE /ISANI2/
      SAVE /RSANI/
      SAVE /LSANI/
