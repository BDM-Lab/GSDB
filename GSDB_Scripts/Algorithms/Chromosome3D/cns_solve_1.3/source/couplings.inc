C===============
C
C Couplings.inc
C
C by John Kuszewski July 1993
C===============
      INTEGER MAXCOUPCLASSES
      PARAMETER (MAXCOUPCLASSES = 10)
C
C arrays that hold coupling info
C coupassndx tells ending index of the coup_x arrays (below)
C for each class.
C coupforces holds k1 and coupswitches holds k2 for each class
C coupclasstypes holds the type of each class
C
C these are small enough not to bother with the heap.
C
      INTEGER COUPASSNDX (MAXCOUPCLASSES)
      INTEGER COUPSWITCHES(MAXCOUPCLASSES)
      DOUBLE PRECISION COUPFORCES (2,MAXCOUPCLASSES)
      CHARACTER*8 COUPCLASSNAMES (MAXCOUPCLASSES)
      LOGICAL PRINTCLASS(MAXCOUPCLASSES)
      INTEGER COUPPOTENTIAL(MAXCOUPCLASSES)
C
C holders for the Karplus curve coefficients
C for each class
C
      DOUBLE PRECISION COUPAS(MAXCOUPCLASSES),
     &     COUPBS(MAXCOUPCLASSES),
     &     COUPCS(MAXCOUPCLASSES), COUPPHASES(MAXCOUPCLASSES)
C
C maxcoups = number of slots set aside for coupling constant
C            assignments
C ncoups = total number of coupling constants entered
C
      INTEGER MAXCOUPS, NCOUPS, NCLASSES, CURCLASS
C
C pointers to arrays to hold atom numbers, J-obs, and
C error for each coupling assignment
C
      INTEGER COUPIPTR, COUPJPTR, COUPKPTR, COUPLPTR,
     &     COUPMPTR, COUPNPTR, COUPPPTR, COUPQPTR,
     &     COUPJOBSPTR, COUPJERRPTR, CALCCOUPPTR
C
C  cross-validation
C  Alexandre Bonvin Dec 95
C
C cross-validation array
      INTEGER COUPCV
C
C cross-validation test number
      INTEGER JICV
C
C input modes
C
      INTEGER MODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
       INTEGER SQUARE, HARMONIC, MULTIPLE
       PARAMETER (SQUARE = 1)
       PARAMETER (HARMONIC = 2)
       PARAMETER (MULTIPLE = 3)
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
C common blocks
C
      COMMON /CCOUP/ COUPCLASSNAMES
      COMMON /ICOUP/ COUPASSNDX, COUPPOTENTIAL,
     &     MAXCOUPS, NCOUPS, COUPSWITCHES, CURCLASS, NCLASSES,
     &     COUPIPTR, COUPJPTR, COUPKPTR, COUPLPTR,
     &     COUPMPTR, COUPNPTR, COUPPPTR, COUPQPTR,
     &     COUPJOBSPTR, COUPJERRPTR, CALCCOUPPTR, MODE,
     &     JICV, COUPCV
      COMMON /RCOUP/ COUPFORCES,
     &     COUPAS, COUPBS, COUPCS, COUPPHASES
      COMMON /LCOUP/ PRINTCLASS
      SAVE /CCOUP/
      SAVE /ICOUP/
      SAVE /RCOUP/
      SAVE /LCOUP/
