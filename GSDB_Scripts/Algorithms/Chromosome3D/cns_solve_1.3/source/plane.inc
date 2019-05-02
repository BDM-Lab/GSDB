C data structure PLANE
C ====================
C
C Parameters for planarity restraints
C
C max. no. groups
      INTEGER MAXPGR
      PARAMETER (MAXPGR=500)
C=====================================================================
C The following are for REstraints option
C flag LPLNRE indicates that we're using PLN restraints
      LOGICAL LPLNRE
C heap pointer and length arrays for equivalent atoms
      INTEGER HPSPNP(MAXPGR), LPSPNP(MAXPGR)
C No. of atoms w/in a group, no. groups
      INTEGER NATPLN(MAXPGR),NGRUPP
C
C WTPLN is effective force constant for planarity restraints
      DOUBLE PRECISION WTPLN(MAXPGR)
C======================================================================
C
C     double precision block
C
      COMMON/PLN1/ WTPLN
C
C     integer block
C
      COMMON/PLN2/ HPSPNP, NATPLN, NGRUPP, LPSPNP
C
C     logical block
C
      COMMON /PLN3/ LPLNRE
C
      SAVE /PLN1/
      SAVE /PLN2/
      SAVE /PLN3/
