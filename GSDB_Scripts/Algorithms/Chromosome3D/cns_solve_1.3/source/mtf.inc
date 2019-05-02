C mtf.inc
C
C
C molecular topology data structures
C
C     Purpose:
C
C     To store the lists of bonds, angles, torsions, and improper
C     torsions for the entire structure, as well as information needed
C     to generate the non-bonded and hydrogen bond lists.
C
C atoms
C     NATOM  Number of atoms
C     MAXA   Maximum number of atoms
C     AMASS  Mass of each atom
C     CG     Charge of each atom
C     IAC    chemical type
C     LOOKUP non-bonding parameter lookup index
C     QLOOKU logical flag indicating if atom-based parameter
C     IMOVE  Flag indicating whether this atom can move.
C            0 if mobile, 1 if fixed
C     TYPE   IUPAC Name
C     RES    Residue Name
C     RESID  Residue Identifier
C     SEGID  Segment Identifier
C     ALTID  Alternate Atom Identifier
C
      INTEGER MAXA
      PARAMETER (MAXA=MAX_ATOM)
      INTEGER NATOM
      DOUBLE PRECISION  CG(MAXA), AMASS(MAXA)
      INTEGER LOOKUP(MAXA), IMOVE(MAXA)
      LOGICAL QLOOKU(MAXA)
      CHARACTER*4 TYPE(MAXA), RESID(MAXA), RES(MAXA), SEGID(MAXA)
      CHARACTER*1 ALTID(MAXA)
      CHARACTER*4 IAC(MAXA)
C
C explicit nonbonded exclusions
C     NNB    Number of non-bonded exclusions
C     MAXNB  Maximum number of non-bonded RTF specific exclusions
C     IBLO   Points to INB giving last non-bonded exclusion
C            for this atom
C     INB    Non-bonded exclusions, indexed by IBLO
C
      INTEGER MAXNB
      PARAMETER (MAXNB=MAX_ATOM)
      INTEGER NNB, IBLO(MAXA), INB(MAXNB)
C
C bonds
C     MAXB   Maximum number of bonds
C     NBOND  Number of bonds
C     IB     First atom of a bond
C     JB     Second atom of bond
C     ACBB   Bond distance
C     ACBC   Force constant
C     QACBB  logical flag indicating if atom-based parameter is set
C     QACBC  logical flag indicating if atom-based parameter is set
C
      INTEGER MAXB
      PARAMETER (MAXB=MAX_ATOM)
      INTEGER  NBOND, IB(MAXB),JB(MAXB)
      DOUBLE PRECISION ACBB(MAXB), ACBC(MAXB)
      LOGICAL QACBB(MAXB), QACBC(MAXB)
C
C angles
C     MAXT    Maximum number of bond angles
C     NTHETA  Number of angles
C     IT      First atom of bond angle
C     JT      Second atom of bond angle
C     KT      Third atom of bond angle
C     ACTB    Angle
C     ACTC    Force constant
C     ACTUB   Urey-Bradley distance
C     ACTUC   Urey-Bradley force constant
C     QACTB  logical flag indicating if atom-based parameter is set
C     QACTC  logical flag indicating if atom-based parameter is set
C     QACTUB  logical flag indicating if atom-based parameter is set
C     QACTUC  logical flag indicating if atom-based parameter is set
C
      INTEGER MAXT
      PARAMETER (MAXT=2*MAX_ATOM)
      INTEGER NTHETA, IT(MAXT), JT(MAXT), KT(MAXT)
      DOUBLE PRECISION ACTB(MAXT), ACTC(MAXT), ACTUB(MAXT), ACTUC(MAXT)
      LOGICAL QACTB(MAXT), QACTC(MAXT), QACTUB(MAXT), QACTUC(MAXT)
C
C dihedrals
C     MAXP  Maximum number of dihedrals
C     NPHI  Number of dihedrals
C     IP    First atom of dihedral
C     JP    Second atom of dihedral
C     KP    Third atom of dihedral
C     JP    Fourth atom of dihedral
C     ACPB  Phase shift for dihedral
C     ACPC  Force constant
C     ACPD  Periodicity
C     QACPB  logical flag indicating if atom-based parameter is set
C     QACPC  logical flag indicating if atom-based parameter is set
C     QACPD  logical flag indicating if atom-based parameter is set
C
      INTEGER MAXP
      PARAMETER (MAXP=2*MAX_ATOM)
      INTEGER NPHI, IP(MAXP), JP(MAXP), KP(MAXP), LP(MAXP), ACPD(MAXP)
      DOUBLE PRECISION ACPB(MAXP), ACPC(MAXP)
      LOGICAL QACPB(MAXP), QACPC(MAXP), QACPD(MAXP)
C
C impropers
C     MAXIMP  Maximum number of improper torsions
C     NIMPHI  Number of improper torsions
C     IM      First atom of improper torsion
C     JM      Second atom of improper torsion
C     KM      Third atom of improper torsion
C     LM      Fourth atom of improper torsion
C     ACIB    Improper angle
C     ACIC    Force constant
C     ACID    Periodicity of the improper energy
C     QACIB  logical flag indicating if atom-based parameter is set
C     QACIC  logical flag indicating if atom-based parameter is set
C     QACID  logical flag indicating if atom-based parameter is set
C
      INTEGER MAXIMP
      PARAMETER (MAXIMP=MAX_ATOM)
      INTEGER NIMPHI, IM(MAXIMP), JM(MAXIMP), KM(MAXIMP), LM(MAXIMP)
      INTEGER ACID(MAXIMP)
      DOUBLE PRECISION ACIB(MAXIMP), ACIC(MAXIMP)
      LOGICAL QACIB(MAXIMP), QACIC(MAXIMP), QACID(MAXIMP)
C
C group linked list
C     MAXGRP   maximum number of groups
C     NGRP     number of groups
C     IGPBS    base pointer to first atom in each group.
C
      INTEGER MAXGRP
      PARAMETER (MAXGRP=MAX_ATOM)
      INTEGER NGRP
      INTEGER IGPBS(MAXGRP+1)
C
C common blocks
C =============
C
C     double precision block
C
      COMMON/MTF/ CG, AMASS, ACBB, ACBC, ACTB, ACTC, ACTUB, ACTUC,
     &    ACPB, ACPC, ACIB, ACIC
C
C     integer block
C
      COMMON/IMTF/
     &    NATOM, LOOKUP, IMOVE,
     &    NNB, INB, IBLO,
     &    NBOND, IB, JB,
     &    NTHETA, IT, JT, KT,
     &    NPHI, IP, JP, KP, LP, ACPD,
     &    NIMPHI, IM, JM, KM, LM, ACID,
     &    NGRP, IGPBS
C
C     character string block
C
      COMMON/CMTF/ TYPE, SEGID, ALTID, RES, RESID, IAC
C
C
C     logical block
C
      COMMON/QMTF/ QLOOKU, QACBB, QACBC, QACTB, QACTC, QACTUB, QACTUC,
     &             QACPB, QACPC, QACPD, QACIB, QACIC, QACID
C
      SAVE /MTF/
      SAVE /IMTF/
      SAVE /CMTF/
      SAVE /QMTF/
