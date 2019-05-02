C===============
C
C Angledb.inc
C
C helpful hint:  MAKE SURE YOU CHANGE THE NAMES OF THE COMMON BLOCKS
C                IF YOU USE THIS FOR SOMETHING ELSE (JK Sep 94)
C
C by John Kuszewski Apr 1996
C===============
      LOGICAL ANGLEDBFLAG
C
      INTEGER MAXANGLEDBCLASSES
      PARAMETER (MAXANGLEDBCLASSES = 200)
C
C arrays that hold angle database setup stuff
C angledbassndx tells ending index of the angledb_x arrays (below)
C for each class.
C angledbforces holds force const for each class
C
C these are small enough not to bother with the heap.
C
      INTEGER ANGLEDBASSNDX (MAXANGLEDBCLASSES),
     &        ANGLEDBCLASSTYPE(MAXANGLEDBCLASSES),
     &        EXPANGLEDBPTRS(MAXANGLEDBCLASSES),
     &        CALCANGLEDBPTR
      DOUBLE PRECISION ANGLEDBFORCES (MAXANGLEDBCLASSES),
     &                 ANGLEDBERRORS(MAXANGLEDBCLASSES)
      CHARACTER*20 ANGDBCLASSNAMES (MAXANGLEDBCLASSES)
      LOGICAL PRINTTHISCLASS(MAXANGLEDBCLASSES),
     &        ANGLEDBDERIV(MAXANGLEDBCLASSES)
C
C whether or not we need to write the remarks to the PDB file
C
C
C maxangledbs = number of slots set aside for angle database entries
C nangledbs = total number of angle database entries
C
      INTEGER MAXANGLEDBS, NANGLEDBS, NANGLEDBCLASSES, CURANGLEDBCLASS
C
C exponent to apply to energy surface
C
      DOUBLE PRECISION ANGLEDBEXPONENT
C
C pointers to arrays to hold atom numbers
C
      INTEGER ANGLEDBIPTR, ANGLEDBJPTR, ANGLEDBKPTR, ANGLEDBLPTR,
     &        ANGLEDBMPTR, ANGLEDBNPTR, ANGLEDBPPTR, ANGLEDBQPTR,
     &        ANGLEDBRPTR, ANGLEDBSPTR, ANGLEDBTPTR, ANGLEDBUPTR
C
C number of entries along the phi and psi dimensions
C of the expectation arrays
C
      INTEGER PHISTEPS(MAXANGLEDBCLASSES),
     &        PSISTEPS(MAXANGLEDBCLASSES)
C
C nine phase corrections for each class
C
      DOUBLE PRECISION ANGLEDBPHASEC(9, MAXANGLEDBCLASSES)
C
C input modes
C
      INTEGER ANGLEDBMODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
C
C potential types
C
      INTEGER ANGLEDBPOTENTIAL, SQUARE, HARMONIC
      PARAMETER (SQUARE = 1)
      PARAMETER (HARMONIC = 2)
C
C database types
C
      INTEGER ANGLE, TORSION
      PARAMETER (ANGLE = 1)
      PARAMETER (TORSION = 2)
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
C common blocks
C
      COMMON /CANGLEDB/ ANGDBCLASSNAMES
      COMMON /IANGLEDB1/ ANGLEDBASSNDX, ANGLEDBCLASSTYPE,
     &                   EXPANGLEDBPTRS, ANGLEDBPOTENTIAL,
     &                   CALCANGLEDBPTR, MAXANGLEDBS, NANGLEDBS,
     &                   CURANGLEDBCLASS, NANGLEDBCLASSES,
     &                   PHISTEPS, PSISTEPS, ANGLEDBMODE
      COMMON /IANGLEDB2/ ANGLEDBIPTR, ANGLEDBJPTR, ANGLEDBKPTR,
     &                   ANGLEDBLPTR, ANGLEDBMPTR, ANGLEDBNPTR,
     &                   ANGLEDBPPTR, ANGLEDBQPTR, ANGLEDBRPTR,
     &                   ANGLEDBSPTR, ANGLEDBTPTR, ANGLEDBUPTR
      COMMON /DANGLEDB/ ANGLEDBFORCES, ANGLEDBERRORS, ANGLEDBPHASEC,
     &                  ANGLEDBEXPONENT
      COMMON /LANGLEDB/ PRINTTHISCLASS, ANGLEDBDERIV, ANGLEDBFLAG
C
      SAVE /CANGLEDB/
      SAVE /IANGLEDB1/
      SAVE /IANGLEDB2/
      SAVE /DANGLEDB/
      SAVE /LANGLEDB/
