C data structure NCS
C ====================
C
C Parameters for non-crystallographic symmetry re/constraints
C
C max. no. groups
      INTEGER MAXNGR
      PARAMETER (MAXNGR=50)
C=====================================================================
C The following are for REstraints option
C flag LNCSRE indicates that we're using NCS restraints
      LOGICAL LNCSRE
C heap pointer and length arrays for equivalent atoms
      INTEGER HPSPNT(MAXNGR), LPSPNT(MAXNGR)
C No. of equivalences w/in a group, no. of unique atoms w/in a group, no. groups
      INTEGER NUMEQV(MAXNGR),NATNCS(MAXNGR),NGROUP
C
C Positional and B factor weights.  WTNCS is effective force constant for
C positional restraints; SIGBNC is target deviation from average B.
      DOUBLE PRECISION WTNCS(MAXNGR),SIGBNC(MAXNGR)
C======================================================================
C These are for strict constraints option
C LNCSST indicates if we are using strict constraints.
      LOGICAL LNCSST
C Number of NCS operations for structure factor calculation and non-bonded
C interactions
      INTEGER XNNSYM,NNSYM
C NCS rotation matrices and translation vectors for structure factors and
C non-bonded interactions, MNCSYM is maximum number.
      INTEGER MNCSYM
      PARAMETER (MNCSYM=50)
      DOUBLE PRECISION NCSOP(MNCSYM,3,4)
C======================================================================
C
C     double precision block
C
      COMMON/NCS1/ WTNCS, SIGBNC, NCSOP
C
C     integer block
C
      COMMON/NCS2/ HPSPNT, NUMEQV, NATNCS, NGROUP, XNNSYM,
     &             NNSYM, LPSPNT
C
C     logical block
C
      COMMON /NCS3/ LNCSRE, LNCSST
C
      SAVE /NCS1/
      SAVE /NCS2/
      SAVE /NCS3/
