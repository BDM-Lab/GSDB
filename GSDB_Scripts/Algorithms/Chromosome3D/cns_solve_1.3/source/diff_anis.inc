C===============
C
C Diffusion Anisotropy.inc
C
C by Nico Tjandra June 1997
C===============
      INTEGER MAXDANICLASSES
      PARAMETER (MAXDANICLASSES = 10)
C
C arrays that hold anisotropy info
C anisassndx tells ending index of the anis_x arrays (below)
C for each class.
C anisforces holds k1 and anisswitches holds k2 for each class
C anisclasstypes holds the type of each class
C
C these are small enough not to bother with the heap.
C
      INTEGER DANIASSNDX (MAXDANICLASSES), DANISWITCHES(MAXDANICLASSES)
      REAL DANIFORCES (MAXDANICLASSES), DANICOEF0 (MAXDANICLASSES)
      REAL DANICOEF1 (MAXDANICLASSES), DANICOEF2 (MAXDANICLASSES)
      REAL DANICOEF3 (MAXDANICLASSES), DANICOEF4 (MAXDANICLASSES)
      CHARACTER*8 DANICLASSNAMES (MAXDANICLASSES)
      LOGICAL PRINTCLASS(MAXDANICLASSES)
C
C maxanis = number of slots set aside for anisotropy assignments
C ndani = total number of anisotropy entered
C
      INTEGER MAXDANI, NDANI, NCLASSES, CURCLASS
C
C pointers to arrays to hold atom numbers, Anis-obs, and
C error for each aniostropy assignment
C
      INTEGER DANIIPTR, DANIJPTR, DANIKPTR, DANILPTR, DANIMPTR,
     &        DANINPTR,
     &        DANIOBSPTR,
     &        DANIERRPTR, CALCDANIPTR
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
      COMMON /CDANI/ DANICLASSNAMES
      COMMON /IDANI1/ DANIASSNDX, POTENTIAL,
     &       MAXDANI, NDANI, DANISWITCHES, CURCLASS, NCLASSES
      COMMON /IDANI2/ DANIIPTR, DANIJPTR, DANIKPTR, DANILPTR,
     &       DANIMPTR, DANINPTR,
     &       DANIOBSPTR, DANIERRPTR, CALCDANIPTR, MODE
      COMMON /RDANI/ DANIFORCES, DANICOEF0, DANICOEF1,
     &       DANICOEF2, DANICOEF3, DANICOEF4
      COMMON /LDANI/ PRINTCLASS
      SAVE /CDANI/
      SAVE /IDANI1/
      SAVE /IDANI2/
      SAVE /RDANI/
      SAVE /LDANI/
