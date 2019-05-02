C ===============
C
C Hshifts.inc
C
C Stuff for the 1H chemical shift refinement package
C
C Note that all constants come from
C Williamson & Asakura, JMR B 101, 63 (1993)
C
C by John Kuszewski Nov 1994
C ===============
C
      INTEGER MAXHSHIFTCLASSES
      PARAMETER (MAXHSHIFTCLASSES = 10)
C
C arrays that hold proton shift info
C shiftassndx tells ending index of the shift_x arrays (below)
C for each class.
C shiftforces holds k1 and kconstantine for each class
C
C these are small enough not to bother with the heap.
C
      INTEGER HCLASSNDX (MAXHSHIFTCLASSES),
     &     HSHIFTDEGENERACY(MAXHSHIFTCLASSES),
     &     HSHIFTPOTENTIAL(MAXHSHIFTCLASSES)
      DOUBLE PRECISION HSHIFTFORCES (2, MAXHSHIFTCLASSES),
     &     HSHIFTERRORS(MAXHSHIFTCLASSES)
      CHARACTER*8 HSHIFTCLASSNAMES (MAXHSHIFTCLASSES)
      LOGICAL PRINTCLASS(MAXHSHIFTCLASSES)
      LOGICAL HCLASSCLOSED(MAXHSHIFTCLASSES)
C
C number of classes that actually exist, and which one
C we're working with currently
C
      INTEGER NHSHIFTCLASSES, CURHSHIFTCLASS
C
C maximum number of protons selected in all assignments
C
      INTEGER MAXPROTONS
      PARAMETER (MAXPROTONS = 12000)
C
      INTEGER PROTNATOM
C
C input modes
C
      INTEGER HSHIFTMODE, NEW, UPDATE
      PARAMETER (NEW = 1)
      PARAMETER (UPDATE = 2)
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
C potential types
C
       INTEGER SQUARE, HARMONIC
       PARAMETER (SQUARE = 1)
       PARAMETER (HARMONIC = 2)
C
C flag to say whether the RMSD array should be overwritten
C during a prot print... statement
C
        LOGICAL PRINTRMSD
C
C hold cutoff values for each proton shift term
C
      DOUBLE PRECISION ELECCUTOFF, ANISCUTOFF, RINGCUTOFF
      PARAMETER (ELECCUTOFF = 6.0)
      PARAMETER (ANISCUTOFF = 12.0)
      PARAMETER (RINGCUTOFF = 15.0)
      DOUBLE PRECISION ELECLISTCUTOFF, ANISLISTCUTOFF, RINGLISTCUTOFF
C
C hold the number of proton index entries
C
      INTEGER NPROTASSIGNS
      INTEGER NPROTONS
C
C hold the number of anisotropic bonds
C
      INTEGER NANIS
C
C hold the number of heavy atoms for elec calc
C
      INTEGER NHEAVIES
C
C hold the number of rings for the ring current calc
C
      INTEGER NRINGS
C
C hold the number of ring current exceptions (which
C are alphas and amides)
C
      INTEGER NRINGEXCEPS
C
C hold the number of amide protons (which are the
C elec field exceptions)
C
      INTEGER NAMIDES
C
C hold the number of heavy atoms for elec field stuff
C
C hold the constants for anisotropy & electric field stuff
C
      DOUBLE PRECISION SIGMAE, CODELTAX1, CODELTAX2
      DOUBLE PRECISION CNDELTAX1, CNDELTAX2, QHN, QC, QN, QO
      DOUBLE PRECISION COCENTER, CNCENTER
      PARAMETER (SIGMAE = 0.60)
      PARAMETER (CODELTAX1 = -13.0)
      PARAMETER (CODELTAX2 = -4.0)
      PARAMETER (CNDELTAX1 = -11.0)
      PARAMETER (CNDELTAX2 = 1.40)
      PARAMETER (QHN = 0.70)
      PARAMETER (QC = 1.60)
      PARAMETER (QN = -1.70)
      PARAMETER (QO = -2.30)
      PARAMETER (COCENTER = 0.905)
      PARAMETER (CNCENTER = 0.85)
C
C ring current constants
C
      DOUBLE PRECISION RINGCONSTB, PHERINGI, TYRRINGI
      DOUBLE PRECISION TRP6RINGI, TRP5RINGI, HISRINGI
      DOUBLE PRECISION ADE6RINGI, ADE5RINGI
      DOUBLE PRECISION GUA6RINGI, GUA5RINGI
      DOUBLE PRECISION CYTRINGI, THYRINGI, URARINGI
      PARAMETER (RINGCONSTB = 5.4548)
      PARAMETER (PHERINGI = 1.05)
      PARAMETER (TYRRINGI = 0.92)
      PARAMETER (TRP6RINGI = 1.04)
      PARAMETER (TRP5RINGI = 0.90)
      PARAMETER (HISRINGI = 0.43)
      PARAMETER (ADE6RINGI = 0.90)
      PARAMETER (ADE5RINGI = 0.66)
      PARAMETER (GUA6RINGI = 0.30)
      PARAMETER (GUA5RINGI = 0.64)
      PARAMETER (CYTRINGI = 0.28)
      PARAMETER (THYRINGI = 0.11)
      PARAMETER (URARINGI = 0.11)
C
C these are for keeping track of the number of
C worthwhile anis, elec, and ring entries for
C each atom
C
      INTEGER NANISENTRIES, NELECENTRIES, NRINGENTRIES
      INTEGER ALPTR, ALNDXPTR, APPTR, ELPTR, ELNDXPTR, RLPTR, RLNDXPTR
      INTEGER OLDXPTR, OLDYPTR, OLDZPTR
      DOUBLE PRECISION HSHIFTMOVECUTOFF
C
C heap pointers to various things
C
      INTEGER
     &     DHCDXPTR, DHCDYPTR, DHCDZPTR,
     &     DHC1DXPTR, DHC1DYPTR, DHC1DZPTR,
     &     DHC2DXPTR, DHC2DYPTR, DHC2DZPTR,
     &     DHOBS1, DHOBS2, DHST1, DHST2, DHFI1, DHFI2,
     &     AAPTR, ACPTR, ACENPTR, AVFPTR,
     &     EEPTR, EHPTR, ECPTR,
     &     RAPTR, REPTR,
     &     HRCPTR, PROTPTR, RCPTR
C
C flag to say if memory has already been allocated
C
      LOGICAL HSHIFTALLOCATED
C
C hold the printout cutoff value
C
      DOUBLE PRECISION PROTPRINTCUTOFF
C
C common blocks
C
      COMMON /DHSHIFT/ PROTPRINTCUTOFF, HSHIFTFORCES, HSHIFTERRORS,
     &     HSHIFTMOVECUTOFF,
     &     ELECLISTCUTOFF, ANISLISTCUTOFF, RINGLISTCUTOFF
      COMMON /IHSHIFT/ NPROTONS, NPROTASSIGNS,
     &     NANIS, NHEAVIES, NRINGS,
     &     DHCDXPTR, DHCDYPTR, DHCDZPTR,
     &     DHC1DXPTR, DHC1DYPTR, DHC1DZPTR,
     &     DHC2DXPTR, DHC2DYPTR, DHC2DZPTR,
     &     DHOBS1, DHOBS2, DHST1, DHST2, DHFI1, DHFI2,
     &     AAPTR, ACPTR, ACENPTR, AVFPTR, APPTR,
     &     EEPTR, EHPTR, ECPTR, RAPTR, REPTR, RCPTR,
     &     HRCPTR, PROTPTR, NAMIDES,
     &     NRINGEXCEPS, NANISENTRIES, NELECENTRIES, NRINGENTRIES,
     &     ALPTR, ALNDXPTR, ELPTR, ELNDXPTR, RLPTR, RLNDXPTR,
     &     OLDXPTR, OLDYPTR, OLDZPTR, PROTNATOM
      COMMON /IHSHIF1/ HCLASSNDX, HSHIFTDEGENERACY, HSHIFTPOTENTIAL,
     &       CURHSHIFTCLASS, NHSHIFTCLASSES, HSHIFTMODE
      COMMON /CHSHIF/ HSHIFTCLASSNAMES
      COMMON /LHSHIF/ HSHIFTALLOCATED, PRINTCLASS, PRINTRMSD,
     &     HCLASSCLOSED
      SAVE /DHSHIFT/
      SAVE /IHSHIFT/
      SAVE /CHSHIF/
      SAVE /IHSHIF1/
      SAVE /LHSHIF/
