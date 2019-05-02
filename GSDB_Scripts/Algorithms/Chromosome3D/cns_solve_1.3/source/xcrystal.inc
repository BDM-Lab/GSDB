C
C xcrystal.inc
C
C Data structure for diffraction data and objects belonging
C to a particular crystal.
C
C Author: Axel T. Brunger
C
C
C method flag for computation of atomic structure factors (FFT or direct)
      LOGICAL QFFT
C flag indicating that lookup tables are used in the FFT calculation
      LOGICAL QLOOK
C
C
C x-ray diffraction indices:
C XRMREF: max. allocation for reflections
C XRNREF: current number of reflections
C HPH, HPK, HPL: pointers for reflection indices
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
C
C unit cell
      DOUBLE PRECISION XRCELL(9), XRTR(3,3), XRINTR(3,3), XRVOL
C
C space group symmetry operators
C xrnsym: number of space group operators
C xrmsym: maximum number of operators
C xrsymm, xritsy, xrnsym: symmetry operators
C xrsygp: group matrix (id of symm_op1, id of symm_op2)
C         -> (id of symm_op1 * symm_op2)
C xrsyiv: inversion operation (id of symm_op1) -> (id of inverse of symm_op1)
      INTEGER XRMSYM, XRSYTH
      PARAMETER (XRMSYM=200, XRSYTH=24)
      INTEGER XRNSYM, XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER XRSYGP(XRMSYM,XRMSYM)
      INTEGER XRSYIV(XRMSYM)
      LOGICAL QHERM
C
C resolution-dependent binning
C number of bins for all resolution-dependent functions
      INTEGER MBINS
C bin resolution range
      DOUBLE PRECISION XBINLOW, XBINHIGH
C bin mode
      INTEGER BINMODE
C bin set flag
      LOGICAL QBINSET
C bin shells
      INTEGER MAXMBIN
      PARAMETER (MAXMBIN=1000)
      DOUBLE PRECISION BINSHELL(MAXMBIN)
C
C
C flag indicating that fast data expansion/reduction is possible
      LOGICAL XQUICK
C
C flag indicating if data is expanded/reduced (for XQUICK mode)
      LOGICAL XQEXPAN
C
C number of reflections in reduced mode
      INTEGER XRNREFRED
C
C number of reflections in expanded mode
      INTEGER XRNREFEXP
C
C number of symmetry operators in reduced mode
      INTEGER XRNSYMRED
C
C asymmetric unit stack (no of command lines) in reduced mode
      INTEGER ARPNNRED
C
C resolution of maps for this crystal
      DOUBLE PRECISION MAPR
C
C map dimension in P1 ( a, b, c ) for this crystal
      INTEGER NA, NB, NC
C map sublattice
      INTEGER NAP, NBP, NCP, NAPP, NBPP, NCPP
C brick that encloses the asymmetric unit (MAASY,...,NAASY;
C MBASY,...,NBASY; MCASY,...,NCASY).
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
C total number of grid points in the brick (size of
C RHOMASK matrix)
      INTEGER NRHO
C
C IRHO is the number of points in the imaginary part of the maps
      INTEGER IRHO
C pointer for integer mask defining the ASU and special positions
C (RHOMASK)
      INTEGER HPRHOMA
C number of non-zero elements in RHOMASK
      INTEGER NMASK
C
C user-defined real space objects
C ===============================
      INTEGER XRHOMX, XNAMEMX
C XNAMEMX must match RPNX in expression.inc!!
      PARAMETER (XRHOMX=100, XNAMEMX=20)
C number of defined objects
      INTEGER XRHONUM
C names of real space objects
      CHARACTER*(XNAMEMX) XRHONAM(XRHOMX)
C heap pointers for real and imaginary parts of real space objects
      INTEGER HPRRHO(XRHOMX), HPIRHO(XRHOMX)
C ===============================
C
C
C XRMAP flag: when XRMAP gets set, the asymmetric unit for maps
C needs to be defined.
      LOGICAL XRMAP
C
C XRRED flag: when XRRED gets set, the reflection data
C based needs to be reduced to an asymmetric unit
      LOGICAL XRRED
C
C
C user-defined reciprocal space objects
C =====================================
C
C maximum number of objects
      INTEGER XSFMX
      PARAMETER (XSFMX=300)
C number of defined objects
      INTEGER XSFNUM
C type of reciprocal space objects
      CHARACTER*4 XSFTYPE(XSFMX)
C XNAMEMX must match RPNX in expression.inc!!
C names of reciprocal space objects
      CHARACTER*(XNAMEMX) XSFNAM(XSFMX)
C type of reciprocal space object groups
      CHARACTER*4 XSFGTYP(XSFMX)
C id of reciprocal space object groups (allocated internally)
      INTEGER XSFGNAM(XSFMX)
C order of objects in reciprocal object space groups
      INTEGER XSFGORD(XSFMX)
C heap pointers for reciprocal space objects
      INTEGER HPSF(XSFMX)
C =====================================
C
C multiplicity for reach reflection
      INTEGER HPMULT
C
C type (-1: centric, 1: acentric, 0: systematic absence)
      INTEGER HPTYPE
C
C command stack that defines the ASU expression
C ---------------------------------------------
      INTEGER ARPNMX, ARPNX
      PARAMETER (ARPNMX=200, ARPNX=XRNAME_SIZE)
      CHARACTER*(ARPNX) ARPN(4,ARPNMX)
      INTEGER ARPNL(4,ARPNMX), ARPNN
C arpndb constains double precision/double complex constants
      DOUBLE COMPLEX ARPNDB(4,ARPNMX)
C arpnmlt is the number of arguments for functions
      INTEGER ARPNMLT(ARPNMX)
C arpntyp is the data type of the command, arpndom is the domain  type
      CHARACTER*2 ARPNTYP(ARPNMX), ARPNDOM(ARPNMX)
C arpnlev is the level of the operation
      INTEGER ARPNLEV(ARPNMX)
C adepth is the depth of the expression
      INTEGER ADEPTH
C
C     double precision common block
C
      COMMON /DXCRYST/ MAPR, ARPNDB,
     &                 XRCELL, XRTR, XRINTR,
     &                 XRVOL, XBINLOW,
     &                 XBINHIGH, BINSHELL
C
C     integer common block
C
      COMMON /IXCRYST/ NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP,
     &                 MAASY, MBASY, MCASY, NAASY, NBASY, NCASY,
     &                 NRHO, HPRHOMA, NMASK,
     &                 ARPNL, ARPNN, ARPNMLT, ARPNLEV, ADEPTH,
     &                 IRHO, HPMULT, HPTYPE, HPRRHO, HPIRHO, HPSF,
     &                 XSFNUM, XRHONUM, XSFGNAM,
     &                 XRMREF, XRNREF, HPH, HPK, HPL,
     &                 XRNSYM, XRSYMM, MBINS, XRITSY,
     &                 XRSYGP, XRSYIV, BINMODE,
     &                 XRNREFRED, XRNREFEXP, XRNSYMRED, ARPNNRED
C
C     logical common block
C
      COMMON /LXCRYST/ XRMAP, XRRED,
     &                 QFFT, QLOOK, QHERM, QBINSET,
     &                 XQUICK, XQEXPAN
C
C     character common block
C
      COMMON /CXCRYST/ ARPN, ARPNTYP, ARPNDOM, XRHONAM, XSFNAM, XSFTYPE,
     &                 XSFGTYP, XSFGORD
C
      SAVE /DXCRYST/
      SAVE /IXCRYST/
      SAVE /LXCRYST/
      SAVE /CXCRYST/
