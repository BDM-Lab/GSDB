      SUBROUTINE PLNPAR
C
C Parse planarity restraints
C Author: Peter Schultze
C                           based on the NCS subroutines
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'plane.inc'
      INCLUDE 'mtf.inc'
C local
      INTEGER PLNFLG, POINTR
C
C begin
C allocate space on heap for temporary pointer arrays
      PLNFLG=ALLHP(INTEG4(NATOM))
      POINTR=ALLHP(INTEG4(NATOM))
      CALL RESPLN(HEAP(PLNFLG),HEAP(POINTR))
C free space on heap
      CALL FREHP(POINTR,INTEG4(NATOM))
      CALL FREHP(PLNFLG,INTEG4(NATOM))
C
      RETURN
      END
C
C----------------------------------------------------------------------
      SUBROUTINE PLNINR
C
C Initialization of planarity restraints
C
C Author: Peter Schultze
C
      IMPLICIT NONE
C I/O
      INCLUDE 'plane.inc'
C local
      INTEGER K
C begin
C
      LPLNRE=.FALSE.
C Restraints
      NGRUPP=0
C initialize heap pointers
      DO K=1,MAXPGR
      HPSPNP(K)=0
      LPSPNP(K)=0
      END DO
      RETURN
      END
C
C---------------------------------------------------------------------
      SUBROUTINE RESPLN(PLNFLG,POINTR)
C
C Parsing routine to set up pointer array POINTR to reference
C groups of atoms to be restrained to planarity.
C
C Author: Peter Schultze
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'plane.inc'
      INCLUDE 'mtf.inc'
      INTEGER PLNFLG(*), POINTR(*)
C local
      INTEGER I,NSELEC
C
C begin parsing
      CALL PUSEND('PLANe>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PLANe>')
      CALL MISCOM('PLANe>',USED)
      IF (.NOT.USED) THEN
C
C====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints-plane')
C
      ELSE IF (WD(1:4).EQ.'INIT') THEN
C
C free-up HEAP space if used before
      CALL PLNFIN
      CALL PLNINR
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
C
C turn on PLN restraints flag
      LPLNRE=.TRUE.
C
C make new group
      IF (NGRUPP.GT.MAXPGR) THEN
      CALL WRNDIE(-5,'PLANe',
     & 'MAXPGR (max. no. of PLANe groups) exceeded')
      ELSE
      NGRUPP=NGRUPP+1
      END IF
C
C fill default values
      NATPLN(NGRUPP)=0
      WTPLN(NGRUPP)=300.0D0
C
      CALL PUSEND('PLANe-GROUp>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PLANe-GROUp>')
      CALL MISCOM('PLANe-GROUp>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-restraints-plane-group')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(PLNFLG,NSELEC,X,Y,Z,.TRUE.)
      CALL MAKIND(PLNFLG,NATOM,NSELEC)
      NATPLN(NGRUPP)=NSELEC
C
C fill temporary pointer array
      DO I=1,NATPLN(NGRUPP)
      POINTR(I)=PLNFLG(I)
      END DO
C
      ELSE IF (WD(1:4).EQ.'WEIG') THEN
      CALL NEXTF('WEIGht=',WTPLN(NGRUPP))
C
      ELSE
      CALL CHKEND('PLANe-GROUp>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C end group parsing: check whether we've done anything sensible
      IF (NATPLN(NGRUPP).LT.4) THEN
      CALL WRNDIE(-1,'PLANe',
     & 'Fewer than 4 atoms selected!')
      NGRUPP=NGRUPP-1
      ELSE
C allocate permanent space on HEAP for list of atom pointers
      LPSPNP(NGRUPP)=NATPLN(NGRUPP)
      HPSPNP(NGRUPP)=ALLHP(INTEG4(LPSPNP(NGRUPP)))
C
C then copy temporary heap pointer into permanent list
      CALL COPYI4(POINTR,HEAP(HPSPNP(NGRUPP)),
     &            NATPLN(NGRUPP))
C
      END IF
C
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      IF (NGRUPP.EQ.0) THEN
      WRITE(6,'(A)') ' No PLANe groups are defined!!!'
      ELSE
C
      DO I=1,NGRUPP
      WRITE(6,'(/,A,I4,A,/,A,I6,/,A,F10.2,A/)')
     &' PLANe GROUP NO. ',I,':',
     &' Number of atoms = ', NATPLN(I),
     &' Effective force constant for PLANe positional restraints =',
     &  WTPLN(I),' Kcal/mol-A**2'
      END DO
C
      END IF
C=====================================================================
      ELSE
      CALL CHKEND('PLANe>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE PLNFIN
C
C Routine allocates space for restraints group lists on heap.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'plane.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER K
C begin
      DO K=1,MAXPGR
      IF (HPSPNP(K).NE.0) THEN
      CALL FREHP(HPSPNP(K),INTEG4(LPSPNP(K)))
      LPSPNP(K)=0
      HPSPNP(K)=0
      END IF
      END DO
      RETURN
      END
C=====================================================================
      SUBROUTINE ENEPLN(EPLN)
C
C Compute planarity restraint energy
C Loop over all groups to allocate space on heap for the
C temporary arrays needed.
C Author: Peter Schultze
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'plane.inc'
      DOUBLE PRECISION EPLN
C local
      INTEGER K
      DOUBLE PRECISION DBPREC
      DOUBLE COMPLEX DCVAL
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C begin
      EPLN=ZERO
C
C loop over all groups
      DO K=1,NGRUPP
C
      CALL ENEPL2(EPLN,HEAP(HPSPNP(K)),NATPLN(K),WTPLN(K))
C
      END DO
C make symbols for number of restraints (groups) NPLAN and
      DBPREC = NGRUPP
      CALL DECLAR( 'NUMBER', 'DP', '  ', DCVAL, DBPREC)
C
      RETURN
      END
C=====================================================================
      SUBROUTINE ENEPL2(EPLN,POINTR,NATPLN,WTPLN)
C
C Effective energy and derivatives for planarity restraints
C Author: Peter Schultze
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      DOUBLE PRECISION EPLN
      INTEGER POINTR(*),  NATPLN
      DOUBLE PRECISION WTPLN
C local
      DOUBLE PRECISION XAVR, YAVR, ZAVR, XDIF, YDIF, ZDIF
      DOUBLE PRECISION CONST, DIST
      DOUBLE PRECISION T(3,3), EIGVAL(3), EIGVEC(3,3)
      INTEGER M,MINEIG
      DOUBLE COMPLEX DCVAL
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C Calculate geometric centre of mass
      XAVR=ZERO
      YAVR=ZERO
      ZAVR=ZERO
C
      DO M=1,NATPLN
      XAVR=XAVR+X(POINTR(M))
      YAVR=YAVR+Y(POINTR(M))
      ZAVR=ZAVR+Z(POINTR(M))
      END DO
C
      XAVR=XAVR/NATPLN
      YAVR=YAVR/NATPLN
      ZAVR=ZAVR/NATPLN
C
C set up the moment of inertia tensor
      T(1,1)=ZERO
      T(1,2)=ZERO
      T(1,3)=ZERO
      T(2,1)=ZERO
      T(2,2)=ZERO
      T(2,3)=ZERO
      T(3,1)=ZERO
      T(3,2)=ZERO
      T(3,3)=ZERO
C
      DO M=1,NATPLN
      XDIF=X(POINTR(M))-XAVR
      YDIF=Y(POINTR(M))-YAVR
      ZDIF=Z(POINTR(M))-ZAVR
C
      T(1,1)=T(1,1)+X(POINTR(M))*XDIF
      T(1,2)=T(1,2)+X(POINTR(M))*YDIF
      T(1,3)=T(1,3)+X(POINTR(M))*ZDIF
      T(2,1)=T(2,1)+Y(POINTR(M))*XDIF
      T(2,2)=T(2,2)+Y(POINTR(M))*YDIF
      T(2,3)=T(2,3)+Y(POINTR(M))*ZDIF
      T(3,1)=T(3,1)+Z(POINTR(M))*XDIF
      T(3,2)=T(3,2)+Z(POINTR(M))*YDIF
      T(3,3)=T(3,3)+Z(POINTR(M))*ZDIF
      END DO
C
C the lowest Eigenvalue of T equals the sum of the squares of the orthogonal
C distances between atoms and plane. The corresponding Eigenvector
C is the normvector of the plane. Do the diagonalisation:
      CALL DIAGSQ(3,3,T,EIGVEC,EIGVAL)
C
C find pointer to the lowest Eigenvalue
      MINEIG=1
      IF( EIGVAL(2) .LT. EIGVAL(MINEIG) ) THEN
      MINEIG=2
      END IF
      IF( EIGVAL(3) .LT. EIGVAL(MINEIG) ) THEN
      MINEIG=3
      END IF
C
C determine energy
      EPLN=EPLN + WTPLN* EIGVAL(MINEIG)
C
C make symbols of the components of the normvector:
      CALL DECLAR( 'PLANX', 'DP', '  ', DCVAL, EIGVEC(1,MINEIG) )
      CALL DECLAR( 'PLANY', 'DP', '  ', DCVAL, EIGVEC(2,MINEIG) )
      CALL DECLAR( 'PLANZ', 'DP', '  ', DCVAL, EIGVEC(3,MINEIG) )
C
C find the constant for the equation of the plane by plugging in
C the centre of mass:
      CONST=EIGVEC(1,MINEIG)*XAVR +
     &      EIGVEC(2,MINEIG)*YAVR +
     &      EIGVEC(3,MINEIG)*ZAVR
C
C find the derivatives
C
      DO M=1,NATPLN
C
C find signed distance of atom from the plane and multiply with constant
C factors, ready for use in the derivatives:
      DIST=EIGVEC(1,MINEIG)* X(POINTR(M)) +
     &     EIGVEC(2,MINEIG)* Y(POINTR(M)) +
     &     EIGVEC(3,MINEIG)* Z(POINTR(M))
     &     -CONST
      DIST=DIST *2.0 * WTPLN
C
C sum up gradient components
      DX(POINTR(M))= DX(POINTR(M)) + EIGVEC(1,MINEIG) * DIST
      DY(POINTR(M))= DY(POINTR(M)) + EIGVEC(2,MINEIG) * DIST
      DZ(POINTR(M))= DZ(POINTR(M)) + EIGVEC(3,MINEIG) * DIST
      END DO
C
      RETURN
      END
