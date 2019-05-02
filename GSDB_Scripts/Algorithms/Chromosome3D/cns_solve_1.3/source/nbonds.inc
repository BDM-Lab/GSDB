C nbonds.inc
C =========
C Nonbonded interactions.
C
C Author: Axel T. Brunger
C modified by JJK 8/10/95 to support the attractive-repulsive term
C
C number of Pairs of Interacting Groups NPIGMA4.
C Note: NPIGMA4 must match NPIGMAX in cnst.inc and NPIGMA3 in ener.inc
      INTEGER NPIGMA4
      PARAMETER (NPIGMA4=MAX_PIGS)
C
C -> update flags for interacting pair N.
C UPNBEX(N) ... enforce update of exclusion list, 1-4 interaction lists
C            and initialize the other nonbonded interaction lists
      LOGICAL UPNBEX(NPIGMA4)
C UPNBLS(N) ... enforce update of atom pair nonbonded interaction lists
C            for intra-molecular interactions
      LOGICAL UPNBLS(NPIGMA4)
C UPNBSS(N) ... enforce update of atom pair nonbonded interaction lists
C            for symmetry interactions
      LOGICAL UPNBSS(NPIGMA4)
C UPNBLK ... enforce update of nonbonded parameter lookup table
C            (one update flag for all Pairs of Interacting Groups)
      LOGICAL UPNBLK
C
C -> automatic update.
C NBXTOL ... how far *any* atom is allowed to move before
C            an update is performed
      DOUBLE PRECISION NBXTOL
C
C NBX, NBY, and NBZ are the "previous positions" of the atoms
C since the last update. Note: there is only one previous
C position array for all PIGs.
C INBX ... heap pointer for the x coordinates (NBX)
C INBY ... heap pointer for the y coordinates (NBY)
C INBZ ... heap pointer for the z coordinates (NBZ)
C XNBX ... allocated lengths of the NBX, NBY, NBZ arrays
      INTEGER XNBX, INBX, INBY, INBZ
C
C -> cutoff options.
C CUTNB ... cutoff
      DOUBLE PRECISION CUTNB
C LGROUP ... FLAG: cutoff by groups/ by atoms
      LOGICAL LGROUP
C WMIN ... close contact warning distance
      DOUBLE PRECISION WMIN
C
C exclusion list options
C NBXMOD ... the nonbond exclusion mode, see MAKINB
      INTEGER NBXMOD
C
C NBQSPC ... value indicating the distance at which
C         to exclude zero contacts that arise due to
C         atoms at special positions in cryst. spacegroups
      DOUBLE PRECISION NBQSPC
C
C electrostatic options
C EPS ... dielectric constant (for direct interactions)
      DOUBLE PRECISION EPS
C E14FAC ... factor for 1-4 electrostatics
      DOUBLE PRECISION E14FAC
C LCONS ... FLAG: constant dielectric/"r-dep. screening-function"
      LOGICAL LCONS
C LSHFT ... FLAG: use a shifted/switched potential for electr.
      LOGICAL LSHFT
C LSHFO ... FLAG: force/potential shifting for electr.
      LOGICAL LSHFO
C BSHFO ... define force shifting function (1 or 4)
      DOUBLE PRECISION BSHFO
C
C van der Waals options
C LVSHFT ... FLAG: use a shifted/switched potential for VDW's
      LOGICAL LVSHFT
C FREPEL ... if > 0, multiply vdw radius and turn repel pot. on.
C        ... for harmonic, repulsive terms -- no electrostatics
      DOUBLE PRECISION FREPEL
C CREP ... force constant on repel potential
      DOUBLE PRECISION CREP
C EXPREP ... exponent for repel potential
      INTEGER EXPREP
C EXIREP ... exponent for repel potential
      INTEGER EXIREP
C
C attractive-repulsive term options
C added by JJK
C
      DOUBLE PRECISION RATT
      DOUBLE PRECISION DOW
C
C options for switching/shifting (both electrostatic/van der Waals)
C CTONNB      switching function turn on dist.
      DOUBLE PRECISION CTONNB
C CTOFNB      switching function and shifting turn off
      DOUBLE PRECISION CTOFNB
C QTRUNC      truncation flag
      LOGICAL QTRUNC
C
C inhibit option
      DOUBLE PRECISION INHIBT
C
C
C =================================================================
C nonbonding lists
C ================
C
C -> the intra-molecular nonbonded interaction list for double
C -> selection N.
C NNNB(N) ... total number of intra-molecular nonbonded interactions.
C IJNB(N) ... heap pointer for an array (JNB) that specifies the
C             number of nonbonded interactions with atom I.
C XJNB(N) ... length of the JNB heap allocation.
C LIST(4,N) ... linked list containing the second atom for all
C               selected nonbonded interactions.
C    LIST(1,N): heap pointer to first array
C    LIST(2,N): heap pointer to current array
C    LIST(3,N): position in current array
C    LIST(4,N): number of elements in each array
C
      INTEGER NNNB(NPIGMA4)
      INTEGER IJNB(NPIGMA4), XJNB(NPIGMA4), LIST(4,NPIGMA4)
C
C -> the symmetry-related nonbonded interaction list for interacting
C -> pair of groups N.
C NNSB(N) ... total number of symmetry-related nonbonded interactions
C IJSB(N) ... heap pointer for an array (JSB) that specifies the number
C             of nonbonded symmetry interactions for each atom I.
C XJSB(N) ... length of the JSB heap allocation
C SLIST(4,N) ... linked list containing the second atom for all
C                selected symmetry-related nonbonded interactions
C    SLIST(1,N): heap pointer to first array
C    SLIST(2,N): heap pointer to current array
C    SLIST(3,N): position in current array
C    SLIST(4,N): number of elements in each array
C
      INTEGER NNSB(NPIGMA4)
      INTEGER IJSB(NPIGMA4), XJSB(NPIGMA4), SLIST(4,NPIGMA4)
C
C -> the special 1-4 interaction list for interacting pair N.
C NNB14(N) ... total number of 1-4 interactions.
C
C IINB14(N) ... heap pointer for an array (INB14) that specifies
C               the second atom of all 1-4 interaction pairs of atoms.
C XINB14(N) ... allocated length of the array INB14.
C LIST14(4,N) ... linked list containing the second atom for all
C                 selected 1-4 nonbonded interactions
C    LIST14(1,N): heap pointer to first array
C    LIST14(2,N): heap pointer to current array
C    LIST14(3,N): position in current array
C    LIST14(4,N): number of elements in each array
      INTEGER NNB14(NPIGMA4)
      INTEGER XINB14(NPIGMA4), IINB14(NPIGMA4), LIST14(4,NPIGMA4)
C
C -> the 1-2, 1-3, 1-4 nonbonded "exclusion" list.
C NNBEX(N) ... total number of exclusions.
C IIBLOEX(N) ... heap pointer for an array (IBLOEX) that specifies
C                a pointer for each atom I.  The pointers point
C                into the IBLOEX array.  For an atom I the numbers
C                of all other atoms that form an exclusion with I are
C                listed between the pointer of atom I-1 (exclusive)
C                and the pointer of atom I (inclusive).
C XIBLOEX(N) ... allocated heap length of the array IBLOEX.
C IINBEX(N) ... heap pointer for an array (INBEX) that specifies
C               the second atom of all excluded interaction pairs
C               of atoms.
C XINBEX(N) ... allocated length of the array INBEX.
      INTEGER NNBEX(NPIGMA4)
      INTEGER XIBLOE(NPIGMA4), IIBLOE(NPIGMA4)
      INTEGER XINBEX(NPIGMA4), IINBEX(NPIGMA4)
C
C
C====================================================================
C
C operator for symmetry-related derivatives
      DOUBLE PRECISION XPOPER(3,4)
C
C====================================================================
C
C     double precision  block
C
      COMMON /NBND/ NBXTOL,  CUTNB,   WMIN,    EPS,    E14FAC,
     &              FREPEL,  CREP,    CTONNB,  CTOFNB, XPOPER,
     &              INHIBT,  BSHFO,   RATT, DOW, NBQSPC
C
C     integer block
C
      COMMON /IBND/  NBXMOD,  EXPREP, EXIREP,
     &               NNNB,   IJNB,  XJNB, LIST,
     &               NNSB,   IJSB,  XJSB, SLIST,
     &               NNB14,  XINB14, IINB14, LIST14,
     &               XIBLOE, IIBLOE,  NNBEX,  XINBEX, IINBEX,
     &               INBX,   INBY,    INBZ,   XNBX
C
C     logical block
C
      COMMON /LBND/  UPNBLS, UPNBSS, UPNBEX, UPNBLK, LSHFO,
     &               LGROUP, LCONS,  LSHFT,  LVSHFT, QTRUNC
C
      SAVE /IBND/
      SAVE /NBND/
      SAVE /LBND/
