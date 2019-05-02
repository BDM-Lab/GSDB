C===============
C
C Rama.inc
C
C helpful hint:  MAKE SURE YOU CHANGE THE NAMES OF THE COMMON BLOCKS
C                IF YOU USE THIS FOR SOMETHING ELSE (JK Sep 94)
C
C by John Kuszewski May 1995
C===============
      INTEGER MAXRAMACLASSES
      PARAMETER (MAXRAMACLASSES = 1000)
C
C arrays that hold coupling info
C ramaassndx tells ending index of the rama_x arrays (below)
C for each class.
C ramaforces holds k1 for each class
C
C these are small enough not to bother with the heap.
C
      INTEGER RAMACLASSTYPE(MAXRAMACLASSES)
      INTEGER EXPRAMAPTRS(MAXRAMACLASSES), CALCRAMAPTR
      DOUBLE PRECISION RAMAFORCES(MAXRAMACLASSES)
      DOUBLE PRECISION RAMAERRORS(MAXRAMACLASSES)
      CHARACTER*50 RAMACLASSNAMES(MAXRAMACLASSES)
C
C (potentially) 12 phase corrections for each class
C
      DOUBLE PRECISION RAMAPHASEC(12, MAXRAMACLASSES)
C
C maxramas = number of slots set aside for chemical rama
C             assignments
C nramas = total number of coupling constants entered
C
      INTEGER MAXRAMAS, NRAMAS, NRAMACLASSES, CURRAMACLASS
C
C pointer to array to hold atom numbers and the class of each
C entry
C
      INTEGER RAMAATOMPTR
C
C overall scaling factor for force constants
C
      DOUBLE PRECISION RAMASCALE
C
C number of entries along the phi and psi dimensions
C of the expectation arrays
C
      INTEGER PHISTEPS(MAXRAMACLASSES), PSISTEPS(MAXRAMACLASSES)
      INTEGER CHISTEPS(MAXRAMACLASSES), THTSTEPS(MAXRAMACLASSES)
C
C input modes
C
      INTEGER RAMAMODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
      INTEGER RAMAPOTENTIAL, SQUARE, HARMONIC
      PARAMETER (SQUARE = 1)
      PARAMETER (HARMONIC = 2)
C
C use this to identify no-expectation slots
C This should be the same as in cshifts.inc
C
      DOUBLE PRECISION NOEXPECTATION
      PARAMETER (NOEXPECTATION = 1.0E20)
      DOUBLE PRECISION NOEXPTEST
      PARAMETER (NOEXPTEST = 1.0E17)
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
      LOGICAL ALLOCRAMA
C
C common blocks
C
      COMMON /CRAMA/ RAMACLASSNAMES
      COMMON /IRAMA/ RAMACLASSTYPE, EXPRAMAPTRS, RAMAPOTENTIAL,
     &               CALCRAMAPTR, MAXRAMAS, NRAMAS, CURRAMACLASS,
     &               NRAMACLASSES, PHISTEPS, PSISTEPS, CHISTEPS,
     &               THTSTEPS, RAMAMODE, RAMAATOMPTR
      COMMON /DRAMA/ RAMAFORCES, RAMAERRORS, RAMAPHASEC, RAMASCALE
      COMMON /LRAMA/ ALLOCRAMA
C
      SAVE /CRAMA/
      SAVE /IRAMA/
      SAVE /DRAMA/
      SAVE /LRAMA/
