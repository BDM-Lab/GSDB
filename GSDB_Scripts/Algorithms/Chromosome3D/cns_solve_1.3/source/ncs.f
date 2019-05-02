C
      SUBROUTINE NCSPAR
C
C Parse non-crystallographic symmetry (NCS) qualifiers (re/constraints).
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'mtf.inc'
C local
      INTEGER NCSFLG, POINTR
C
C begin
      CALL NEXTWD('NCS-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
      CALL CNSHELP('cns-ncs')
      ELSE IF (WD(1:4).EQ.'REST') THEN
C allocate space on heap for temporary pointer arrays
      NCSFLG=ALLHP(INTEG4(NATOM))
      POINTR=ALLHP(INTEG4(NATOM))
      CALL RESNCS(HEAP(NCSFLG),HEAP(POINTR))
C free space on heap
      CALL FREHP(POINTR,INTEG4(NATOM))
      CALL FREHP(NCSFLG,INTEG4(NATOM))
C=====================================================================
      ELSE IF (WD(1:4).EQ.'STRI') THEN
      CALL STRNCS
C=====================================================================
      ELSE
      CALL DSPERR ('NCS','unknown qualifier')
      END IF
C
      RETURN
      END
C
C-------------------------------------------------------------------------
      SUBROUTINE NCSINR
C
C Initialization of non-crystallographic symmetry options for
C restraints.
C
C Author: Bill Weis
C
      IMPLICIT NONE
C I/O
      INCLUDE 'ncs.inc'
C local
      INTEGER K
C begin
C
      LNCSRE=.FALSE.
C Restraints
      NGROUP=0
C initialize heap pointers
      DO K=1,MAXNGR
      HPSPNT(K)=0
      LPSPNT(K)=0
      END DO
      RETURN
      END
C
      SUBROUTINE NCSINS
C
C Initialization of non-crystallographic symmetry options for
C strict symmetry constraints.
C
C Author: Bill Weis
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
C local
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DOUBLE PRECISION DBTEMP
      DOUBLE COMPLEX DBCOMP
C begin
      LNCSST=.FALSE.
C
C Constraints
C Number of operators for generating crystallographic asymmetric unit from
C protomer coordinates for structure factor calculations.
      XNNSYM=1
      DBTEMP=XNNSYM
      CALL DECLAR( 'NCS', 'DP', ' ', DBCOMP, DBTEMP )
C
C Number of operators for NCS non-bonded interactions
      NNSYM=1
      NCSOP(1,1,1)=ONE
      NCSOP(1,1,2)=ZERO
      NCSOP(1,1,3)=ZERO
      NCSOP(1,1,4)=ZERO
      NCSOP(1,2,1)=ZERO
      NCSOP(1,2,2)=ONE
      NCSOP(1,2,3)=ZERO
      NCSOP(1,2,4)=ZERO
      NCSOP(1,3,1)=ZERO
      NCSOP(1,3,2)=ZERO
      NCSOP(1,3,3)=ONE
      NCSOP(1,3,4)=ZERO
C
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE STRNCS
C
C Parsing for strict non-crystallographic symmetry constraints
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
C local
      INTEGER I,J,K,L,INSYM, N
      INTEGER LDUM(3),MDUM(3)
      LOGICAL QXNCS, LSKEW
      DOUBLE PRECISION SKWMAT(3,3),SKWINV(3,3),DET,SKWVEC(3)
      DOUBLE PRECISION NCSROT(3,3),NCSTRN(3), DBTEMP
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION PUP(3), PV(3)
C parameter
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C default (no skewing--identity matrix)
      LSKEW=.FALSE.
      DO I=1,3
      DO J=1,3
      SKWMAT(I,J)=ZERO
      SKWINV(I,J)=ZERO
      END DO
      SKWVEC(I)=ZERO
      END DO
C
      DO I=1,3
      SKWMAT(I,I)=ONE
      SKWINV(I,I)=ONE
      END DO
C
C begin parsing
      CALL PUSEND('NCS-strict>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('NCS-strict>')
      CALL MISCOM('NCS-strict>',USED)
      IF (.NOT.USED) THEN
C
C====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-ncs-strict')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'INIT') THEN
C
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
C
      CALL NCSINS
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SKEW') THEN
C
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
C
      LSKEW=.TRUE.
C
C begin parsing
      CALL PUSEND('SKEW>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SKEW>')
      CALL MISCOM('SKEW>',USED)
      IF (.NOT.USED) THEN
C
C====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-ncs-strict-skew')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAN') THEN
      CALL NEXTVF('TRANslation= ',SKWVEC)
C==================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &         WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &         WD(1:4).EQ.'AXIS') THEN
      CALL MATPAR(WD(1:4),SKWMAT)
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL MATPRI(SKWMAT)
      WRITE(6,'(A,3F10.4,A/)')
     & ' Translation vector = (',(SKWVEC(I),I=1,3),')'
C======================================================================
      ELSE
      CALL CHKEND('SKEW>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      DO I=1,3
      WRITE(6,'(10X,3F12.6)') (SKWMAT(I,J),J=1,3)
      END DO
C
      DO I=1,3
      DO J=1,3
      SKWINV(I,J)=SKWMAT(I,J)
      END DO
      END DO
      CALL MINVV(SKWINV,3,DET,LDUM,MDUM)
      WRITE(6,'(A,F12.6,/)')
     . ' Determinant of skew matrix = ', DET
      IF (DET.EQ.0.0D0) THEN
      CALL WRNDIE(-1,'NCS-strict',' Skew matrix is singular!!')
      ELSE
      WRITE(6,'(A)') ' INVERSE OF SKEW MATRIX:'
      DO I=1,3
      WRITE(6,'(10X,3F12.6)') (SKWINV(I,J),J=1,3)
      END DO
      END IF
C
      WRITE(6,'(A)') ' SKEW VECTOR:'
      WRITE(6,'(10X,3F12.6)') SKWVEC
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'XNCS' .OR. WD(1:4).EQ.'NCSR') THEN
C
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
C
C
C turn on NCS strict flag
      LNCSST=.TRUE.
      QXNCS=WD(1:4).EQ.'XNCS'
C
      IF (.NOT.LSKEW) THEN
      WRITE(6,'(A)')
     &' SKEW TRANSFORMATION HAS NOT BEEN GIVEN, ASSUMING IDENTITY'
      END IF
C
C initial
      DO I=1,3
      DO J=1,3
      NCSROT(I,J)=ZERO
      END DO
      NCSROT(I,I)=ONE
      NCSTRN(I)=ZERO
      END DO
C begin parsing
      CALL PUSEND('XNCS>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('XNCS>')
      CALL MISCOM('XNCS>',USED)
      IF (.NOT.USED) THEN
C
C====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-ncs-strict-ncsrelation')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAN') THEN
      CALL NEXTVF('TRANslation= ',NCSTRN)
C==================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &         WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &         WD(1:4).EQ.'AXIS') THEN
      CALL MATPAR(WD(1:4),NCSROT)
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL MATPRI(NCSROT)
      WRITE(6,'(A,3F10.4,A/)')
     & ' Translation vector = (',(NCSTRN(I),I=1,3),')'
C======================================================================
      ELSE
      CALL CHKEND('XNCS>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C If this is the identity transformation, ignore it, as we've already
C set it up for NNSYM=1.
      IF (.NOT. (NCSROT(1,1).EQ.ONE .AND. NCSROT(1,2).EQ.ZERO .AND.
     &    NCSROT(1,3).EQ.ZERO .AND. NCSROT(2,1).EQ.ZERO .AND.
     &    NCSROT(2,2).EQ.ONE .AND. NCSROT(2,3).EQ.ZERO .AND.
     &    NCSROT(3,1).EQ.ZERO .AND. NCSROT(3,2).EQ.ZERO .AND.
     &    NCSROT(3,3).EQ.ONE .AND. NCSTRN(1).EQ.ZERO .AND.
     &    NCSTRN(2).EQ.ZERO .AND. NCSTRN(3).EQ. ZERO) ) THEN
C
      IF (NNSYM.GE.MNCSYM) THEN
      CALL WRNDIE(-5,'NCS',
     & 'exceeded MNCSYM parameter --> recompile program')
      ELSE
      NNSYM=NNSYM+1
      IF (QXNCS) XNNSYM=XNNSYM + 1
      DBTEMP=XNNSYM
      CALL DECLAR( 'NCS', 'DP', ' ', DBCOMP, DBTEMP )
C
      END IF
C
C Set up operator for this relation.
C This is P U P**-1 (= R)
      DO I=1,3
      DO J=1,3
      NCSOP(NNSYM,I,J)=ZERO
      DO K=1,3
      DO L=1,3
      NCSOP(NNSYM,I,J)=NCSOP(NNSYM,I,J) +
     &       SKWMAT(I,K)*NCSROT(K,L)*SKWINV(L,J)
      END DO
      END DO
      END DO
      END DO
C
C This is  P U P**-1 SKWVEC
      DO I=1,3
      PUP(I)=ZERO
      DO J=1,3
      PUP(I)= PUP(I)+NCSOP(NNSYM,I,J)*SKWVEC(J)
      END DO
      END DO
C
C T = -P U P**-1 SKWVEC + P V + SKWVEC
      DO I=1,3
      PV(I)=ZERO
      DO J=1,3
      PV(I)=PV(I) + SKWMAT(I,J)*NCSTRN(J)
      END DO
      END DO
      DO I=1,3
      NCSOP(NNSYM,I,4)= -PUP(I) + PV(I) + SKWVEC(I)
      END DO
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      IF (NNSYM.EQ.1) THEN
      WRITE(6,'(A)') ' No NCS strict symmetry operators are defined.'
      ELSE
      WRITE(6,'(/,A,I2,A,/,A,I2,2A,/,A)')
     &' THERE ARE ',NNSYM,' NCS RELATIONS DEFINED, INCLUDING IDENTITY',
     &' The first ',XNNSYM,' will be used for structure factor',
     &' calculations;',
     &'     all will be used for non-bonded interactions.'
      WRITE(6,'(A)')
     &' The operators give the NCS related atoms in the',
     &' same frame as the input coordinates:'
      DO INSYM=1,NNSYM
      WRITE(6,'(//,A,I2,A)') ' RELATION NO. ',INSYM,':'
      DO I=1,3
      DO J=1,3
      NCSROT(I,J)=NCSOP(INSYM,I,J)
      END DO
      END DO
      CALL MATPRI(NCSROT)
      WRITE(6,'(A,3F10.4,A/)')
     & ' Translation vector = (',(NCSOP(INSYM,I,4),I=1,3),')'
      END DO
      END IF
C======================================================================
      ELSE
      CALL CHKEND('NCS-strict>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C always declare symbols even if ? is not given - PDA 7/98
      IF ( NNSYM .GT. 0 ) THEN
      DO INSYM=1,NNSYM
      DO I=1,3
      DO J=1,3
      NCSROT(I,J)=NCSOP(INSYM,I,J)
      END DO
      END DO
      CALL OPSDCL('NCSOP',5,NCSROT,NCSOP(INSYM,1,4),NCSOP(INSYM,2,4),
     &             NCSOP(INSYM,3,4),ZERO,INSYM)
      END DO
      END IF
C
      RETURN
      END
C
C----------------------------------------------------------------------
      SUBROUTINE OPSDCL(OPNAME,LENNAME,ROT,TX,TY,TZ,RMSD,IOP)
C
C Declare symbols for the operating matrix, translation vector,
C and RMS difference.
C
C The named operators are in the following format, e.g.
C "$OPNAME_$k_$i_$j" corresponding to OPNAME(IOP, row, column).
C where "$k" is the number of the operator; "$i" runs from 1 to 3;
C "$j" runs from 1 to 4.  When "$j=4" represents the translation
C vector.  "$OPNAME_$k_RMSD" represents the r.m.s difference
C of "$k" mate to the first mate.
C
C When IOP=0 the symbols just become "$OPNAME_$i_$j"
C
C Authors: J.-S. Jiang and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*20 OPNAME
      INTEGER LENNAME
      DOUBLE PRECISION ROT(3,3), TX, TY, TZ
      DOUBLE PRECISION RMSD
      INTEGER IOP
C local
      DOUBLE COMPLEX DBCOMP
      CHARACTER*20 OP, OPS
      INTEGER LENOP, LENOPS
C
C begin
      CALL ENCODI(IOP,WDT,WDTMAX,WDTLEN)
      IF (IOP.GT.0) THEN
      OPS=OPNAME(1:LENNAME)//'_'//WDT(1:WDTLEN)
      LENOPS=LENNAME+1+WDTLEN
      ELSE
      OPS=OPNAME(1:LENNAME)
      LENOPS=LENNAME
      END IF
      LENOP=LENOPS+4
C
C rotation matrix
      OP=OPS(1:LENOPS)//'_1_1'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(1,1) )
      OP=OPS(1:LENOPS)//'_1_2'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(1,2) )
      OP=OPS(1:LENOPS)//'_1_3'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(1,3) )
      OP=OPS(1:LENOPS)//'_2_1'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(2,1) )
      OP=OPS(1:LENOPS)//'_2_2'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(2,2) )
      OP=OPS(1:LENOPS)//'_2_3'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(2,3) )
      OP=OPS(1:LENOPS)//'_3_1'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(3,1) )
      OP=OPS(1:LENOPS)//'_3_2'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(3,2) )
      OP=OPS(1:LENOPS)//'_3_3'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, ROT(3,3) )
C
C translation vector
      OP=OPS(1:LENOPS)//'_1_4'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, TX )
      OP=OPS(1:LENOPS)//'_2_4'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, TY )
      OP=OPS(1:LENOPS)//'_3_4'
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, TZ )
C
C r.m.s. difference
      OP=OPS(1:LENOPS)//'_RMSD'
      LENOP=LENOPS+5
      CALL DECLAR( OP(1:LENOP), 'DP', ' ', DBCOMP, RMSD )
C
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE NCSXYZ(NAT,X,Y,Z,INSYM,XNCS,YNCS,ZNCS)
C
C Calculate non-crystallographic symmetry related coordinates XNCS, YNCS,
C ZNCS from input X,Y,Z in orthogonal Angstroms.   Called by FFT and non-
C bonding routines.
C
C Modification: operation only applied to known coordinates
C ATB 11/03/09
C
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'ncs.inc'
      INCLUDE 'funct.inc'
      INTEGER NAT,INSYM
      DOUBLE PRECISION X(*),Y(*),Z(*),XNCS(*),YNCS(*),ZNCS(*)
C
C local
      INTEGER IAT
      DOUBLE PRECISION ANUM
      PARAMETER (ANUM=9999.0D0)
C
C begin
      DO IAT=1,NAT
      IF (INITIA(IAT,X,Y,Z)) THEN
      XNCS(IAT)=NCSOP(INSYM,1,1)*X(IAT) + NCSOP(INSYM,1,2)*Y(IAT)
     .          + NCSOP(INSYM,1,3)*Z(IAT) + NCSOP(INSYM,1,4)
      YNCS(IAT)=NCSOP(INSYM,2,1)*X(IAT) + NCSOP(INSYM,2,2)*Y(IAT)
     .          + NCSOP(INSYM,2,3)*Z(IAT) + NCSOP(INSYM,2,4)
      ZNCS(IAT)=NCSOP(INSYM,3,1)*X(IAT) + NCSOP(INSYM,3,2)*Y(IAT)
     .          + NCSOP(INSYM,3,3)*Z(IAT) + NCSOP(INSYM,3,4)
      ELSE
      XNCS(IAT)=ANUM
      YNCS(IAT)=ANUM
      ZNCS(IAT)=ANUM
      END IF
      END DO
      RETURN
      END
C---------------------------------------------------------------------
      SUBROUTINE RESNCS(NCSFLG,POINTR)
C
C Parsing routine to set up pointer array POINTR to reference non-
C crystallographic symmetry (NCS) related atoms in the structure.  "GROUPS"
C refer to different sets of non-crystallographically related molecules (e.g.
C one set of atoms related by an NCS axis would constitute one group, while
C another set of atoms related by a different axis would constitute a second
C group). Note that in the general case, the same atom can be involved in more
C than one group.  "EQUIVALENCES" refer to NCS equivalent atoms within a group
C (i.e. there are N equivalent copies of an atom around an N-fold NCS axis).
C
C Author: Bill Weis
C =================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'mtf.inc'
      INTEGER NCSFLG(*), POINTR(*)
C local
      INTEGER I,M,NSELEC
      INTEGER XREF,YREF,ZREF,XROT,YROT,ZROT,NCOUNT
C
C begin parsing
      CALL PUSEND('NCS-restraints>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('NCS-restraints>')
      CALL MISCOM('NCS-restraints>',USED)
      IF (.NOT.USED) THEN
C
C====================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-ncs-restraints')
C
      ELSE IF (WD(1:4).EQ.'INIT') THEN
C
C free-up HEAP space if used before
      CALL NCSFIN
      CALL NCSINR
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
C
C turn on NCS restraints flag
      LNCSRE=.TRUE.
C
C make new group
CCC modification ATB 4/27/08
      IF (NGROUP.GE.MAXNGR) THEN
      CALL WRNDIE(-5,'NCS-restraints',
     & 'MAXNGR (max. no. of NCS groups) exceeded')
      ELSE
      NGROUP=NGROUP+1
      END IF
C
C fill default values
      NUMEQV(NGROUP)=0
      NATNCS(NGROUP)=0
      WTNCS(NGROUP)=300.0D0
      SIGBNC(NGROUP)=2.0D0
      NCOUNT=0
C
      CALL PUSEND('NCS-restraints-group>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('NCS-restraints-group>')
      CALL MISCOM('NCS-restraints-group>',USED)
      IF (.NOT.USED) THEN
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-ncs-restraints-group')
C
      ELSE IF (WD(1:4).EQ.'EQUI') THEN
      NUMEQV(NGROUP)=NUMEQV(NGROUP)+1
      CALL SELCTA(NCSFLG,NSELEC,X,Y,Z,.TRUE.)
C
C BEGIN MODIFICATION, make sure unknown atoms are not selected, ATB, 11/09/09
      M=0
      DO I=1,NATOM
      IF ((.NOT.INITIA(I,X,Y,Z)).AND.NCSFLG(I).EQ.1) THEN
      M=M+1
      NCSFLG(I)=0
      NSELEC=NSELEC-1
      END IF
      END DO
      IF (M.GT.0) THEN
      WRITE(6,'(A,/,A,/,A,I7,A,I7)') 
     & ' NCS-WRN: Unknown atoms present in equivalence selection.',
     & ' NCS-WRN: These atoms have been de-seleced.',
     & ' NCS-WRN',NSELEC,' atoms have now been selected out of',NATOM
      END IF
C END MODIFICATION
C
      CALL MAKIND(NCSFLG,NATOM,NSELEC)
C
C for first equivalence store the number of unique atoms
      IF (NUMEQV(NGROUP).EQ.1) THEN
      NATNCS(NGROUP)=NSELEC
      ERROR=.FALSE.
      ELSE
C
C otherwise make several consistency checks
      IF (NSELEC.EQ.0) THEN
      WRITE(6,'(A,/,A,I6,A,I6)')
     &' %NCS-WRN: Zero atoms in equivalence selection: ',
     &' empty equivalence set ',NUMEQV(NGROUP),' in group ',NGROUP
      END IF
C
      IF (NATNCS(NGROUP).NE.NSELEC) THEN
      WRITE(6,'(A,I2,A,/,A,I2,A,I6,A,I6,A)')
     &' %NCS-ERR: For NCS group ',NGROUP,
     &' number of atoms does not match:',
     &' # atoms in equivalence set ',NUMEQV(NGROUP),' = ',NSELEC,
     &' # atoms in equivalence set 1 = ',NATNCS(NGROUP),' !!!'
      ERROR=.TRUE.
      ELSE IF (NATNCS(NGROUP)*NUMEQV(NGROUP).GT.NATOM) THEN
      WRITE(6,'(A,/,A,I6,A,I6,A)')
     &' %NCS-ERR: Too many equivalence specified: ',
     &' number of equiv. (',NUMEQV(NGROUP),
     &') x  number of unique atoms (',
     & NATNCS(NGROUP),') > number of atoms in molecule'
      ERROR=.TRUE.
      END IF
C
C check internal consistency of equivalence specification
C modification in following statement ATB 5/28/08
      DO M=1,MIN(NATNCS(NGROUP),NSELEC)
      IF (RES(POINTR(M)).NE.RES(NCSFLG(M)) ) THEN
      WRITE(6,'(A,I2,A,I2,A,I5,A,/,A,A4,A,A4)')
     &' %NCS-ERR: For NCS group ',NGROUP,', equivalence set ',
     &  NUMEQV(NGROUP),
     &', atom ',M,' :',' RES = ',RES(POINTR(M)),
     &', while RES from set 1 = ',RES(NCSFLG(M))
      ERROR=.TRUE.
C
      ELSE IF (TYPE(POINTR(M)).NE.TYPE(NCSFLG(M)) ) THEN
      WRITE(6,'(A,I2,A,I2,A,I5,A,/,A,A4,A,A4)')
     &' %NCS-ERR: For NCS group ',NGROUP,', equivalence set ',
     & NUMEQV(NGROUP),
     &', atom ',M,' :',' TYPE = ',TYPE(POINTR(M)),
     &', while TYPE from set 1 = ',TYPE(NCSFLG(M))
C
      ERROR=.TRUE.
      END IF
      END DO
      END IF
C
      IF (ERROR) THEN
      CALL WRNDIE(-1,'NCS-restraints',
     & 'Improperly defined non-crystallographic symmetry')
      END IF
C
C fill temporary pointer array
      DO I=1,NATNCS(NGROUP)
      NCOUNT=NCOUNT+1
      POINTR(NCOUNT)=NCSFLG(I)
      END DO
C
      ELSE IF (WD(1:4).EQ.'WEIG') THEN
      CALL NEXTF('WEIGht-ncs=',WTNCS(NGROUP))
      ELSE IF (WD(1:4).EQ.'SIGB') THEN
      CALL NEXTF('SIGB-ncs=',SIGBNC(NGROUP))
C
      ELSE
      CALL CHKEND('NCS-restraints-group>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C end group parsing: check whether we've done anything at all
      IF (NCOUNT.EQ.0) THEN
      NGROUP=NGROUP-1
      ELSE
C allocate permanent space on HEAP for list of atom pointers
      LPSPNT(NGROUP)=NUMEQV(NGROUP)*NATNCS(NGROUP)
      HPSPNT(NGROUP)=ALLHP(INTEG4(LPSPNT(NGROUP)))
C
C then copy temporary heap pointer into permanent list
      CALL COPYI4(POINTR,HEAP(HPSPNT(NGROUP)),
     &            NUMEQV(NGROUP)*NATNCS(NGROUP))
C
      END IF
C
C======================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      IF (NGROUP.EQ.0) THEN
      WRITE(6,'(A)') ' No NCS equivalences are defined!!!'
      ELSE
C
      WRITE(6,'(2A)')
     & ' NOTE: The first EQUIvalence set defined in a GROUp is taken',
     & ' as the reference!'
C allocate temporary space to print matrices and vectors.
      DO I=1,NGROUP
      XREF=ALLHP(IREAL8(NATNCS(I)))
      YREF=ALLHP(IREAL8(NATNCS(I)))
      ZREF=ALLHP(IREAL8(NATNCS(I)))
      XROT=ALLHP(IREAL8(NATNCS(I)))
      YROT=ALLHP(IREAL8(NATNCS(I)))
      ZROT=ALLHP(IREAL8(NATNCS(I)))
      CALL PRINCS(NUMEQV(I),NATNCS(I),HEAP(HPSPNT(I)),
     &       HEAP(XREF),HEAP(YREF),HEAP(ZREF),
     &       HEAP(XROT),HEAP(YROT),HEAP(ZROT),I)
      CALL FREHP(ZROT,IREAL8(NATNCS(I)))
      CALL FREHP(YROT,IREAL8(NATNCS(I)))
      CALL FREHP(XROT,IREAL8(NATNCS(I)))
      CALL FREHP(ZREF,IREAL8(NATNCS(I)))
      CALL FREHP(YREF,IREAL8(NATNCS(I)))
      CALL FREHP(XREF,IREAL8(NATNCS(I)))
      END DO
C
      DO I=1,NGROUP
      WRITE(6,'(/,A,I2,A,/,A,F10.2,A,/,A,F8.3,A,/)')
     &' NCS GROUP NO. ',I,':',
     &' Effective force constant for NCS positional restraints =',
     &  WTNCS(I),' Kcal/mol-A**2',
     &' Target deviation of NCS related B factors from average = ',
     &  SIGBNC(I),' A**2'
      END DO
C
      END IF
C=====================================================================
      ELSE
      CALL CHKEND('NCS-restraints>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      RETURN
      END
C----------------------------------------------------------------
      SUBROUTINE NCSFIN
C
C Routine allocates space for restraints group lists on heap.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'ncs.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER K
C begin
      DO K=1,MAXNGR
      IF (HPSPNT(K).NE.0) THEN
      CALL FREHP(HPSPNT(K),INTEG4(LPSPNT(K)))
      LPSPNT(K)=0
      HPSPNT(K)=0
      END IF
      END DO
      RETURN
      END
C=====================================================================
      SUBROUTINE PRINCS(NUMEQV,NATNCS,POINTR,
     &                 XREF,YREF,ZREF,XROT,YROT,ZROT,NGROUP)
C
C Print current values of all non-crystallographic symmetry rotation and
C translation relationships for RESTraints option.
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'comand.inc'
      INTEGER NUMEQV,NATNCS,POINTR(*)
      DOUBLE PRECISION XREF(*),YREF(*),ZREF(*),XROT(*),YROT(*),ZROT(*)
      INTEGER NGROUP
C local
      INTEGER J,M, EQPOS
      DOUBLE PRECISION UU(9), ROT(3,3)
      CHARACTER*20 OPS
      INTEGER LENOPS
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      WRITE(6,'(A,I2,A,I2,A)') ' NCS group ',NGROUP,' has ',
     & NUMEQV,' sets of equivalent atoms:'
C
C The first subgroup is the reference.
      DO M=1,NATNCS
      XREF(M)=X(POINTR(M))
      YREF(M)=Y(POINTR(M))
      ZREF(M)=Z(POINTR(M))
      END DO
C
C declare symbols for identity operator - for completeness
      ROT(1,1)=ONE
      ROT(1,2)=ZERO
      ROT(1,3)=ZERO
      ROT(2,1)=ZERO
      ROT(2,2)=ONE
      ROT(2,3)=ZERO
      ROT(3,1)=ZERO
      ROT(3,2)=ZERO
      ROT(3,3)=ONE
      CALL ENCODI(NGROUP,WDT,WDTMAX,WDTLEN)
      OPS='ROT_'//WDT(1:WDTLEN)
      LENOPS=4+WDTLEN
      CALL OPSDCL(OPS,LENOPS,ROT,ZERO,ZERO,ZERO,ZERO,1)
      CALL ENCODI(NGROUP,WDT,WDTMAX,WDTLEN)
      OPS='ROTINV_'//WDT(1:WDTLEN)
      LENOPS=7+WDTLEN
      CALL OPSDCL(OPS,LENOPS,ROT,ZERO,ZERO,ZERO,ZERO,1)
C
C loop over all equivalences
      DO J=2,NUMEQV
      EQPOS=NATNCS*(J-1)
      DO M=1,NATNCS
      XROT(M)=X(POINTR(M+EQPOS))
      YROT(M)=Y(POINTR(M+EQPOS))
      ZROT(M)=Z(POINTR(M+EQPOS))
      END DO
C
      WRITE(6,'(//,A,I2,A,/)') '  Equivalence set ',J,' :'
C
      CALL ROTATE(XREF,YREF,ZREF,NATNCS,XROT,YROT,ZROT,UU,.TRUE.,
     &            NGROUP,J)
C
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------
      SUBROUTINE ENENCS(ENCS)
C
C Compute non-crystallographic symmetry energy (actually done in s/r ENENCS).
C Loop over all groups to allocate space on heap for the
C temporary arrays needed.
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'ncs.inc'
      DOUBLE PRECISION ENCS
C local
      INTEGER K
      INTEGER XREF,YREF,ZREF
      INTEGER XBAR,YBAR,ZBAR
      INTEGER XROT,YROT,ZROT
      INTEGER RX,RY,RZ
      INTEGER U
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
C begin
      ENCS=ZERO
C
C loop over all groups
      DO K=1,NGROUP
      XREF=ALLHP(IREAL8(NATNCS(K)))
      YREF=ALLHP(IREAL8(NATNCS(K)))
      ZREF=ALLHP(IREAL8(NATNCS(K)))
      XBAR=ALLHP(IREAL8(NATNCS(K)))
      YBAR=ALLHP(IREAL8(NATNCS(K)))
      ZBAR=ALLHP(IREAL8(NATNCS(K)))
      XROT=ALLHP(IREAL8(NATNCS(K)))
      YROT=ALLHP(IREAL8(NATNCS(K)))
      ZROT=ALLHP(IREAL8(NATNCS(K)))
      RX=ALLHP(IREAL8(NUMEQV(K)*NATNCS(K)))
      RY=ALLHP(IREAL8(NUMEQV(K)*NATNCS(K)))
      RZ=ALLHP(IREAL8(NUMEQV(K)*NATNCS(K)))
      U=ALLHP(IREAL8(9*NUMEQV(K)*NATNCS(K)))
C
      CALL ENENC2(ENCS,HEAP(HPSPNT(K)),NUMEQV(K),NATNCS(K),
     & HEAP(XREF),HEAP(YREF),HEAP(ZREF),HEAP(XBAR),
     & HEAP(YBAR),HEAP(ZBAR),HEAP(XROT),HEAP(YROT),HEAP(ZROT),
     & HEAP(RX),HEAP(RY),HEAP(RZ),HEAP(U),WTNCS(K))
C
      CALL FREHP(U,IREAL8(9*NUMEQV(K)*NATNCS(K)))
      CALL FREHP(RZ,IREAL8(NUMEQV(K)*NATNCS(K)))
      CALL FREHP(RY,IREAL8(NUMEQV(K)*NATNCS(K)))
      CALL FREHP(RX,IREAL8(NUMEQV(K)*NATNCS(K)))
      CALL FREHP(ZROT,IREAL8(NATNCS(K)))
      CALL FREHP(YROT,IREAL8(NATNCS(K)))
      CALL FREHP(XROT,IREAL8(NATNCS(K)))
      CALL FREHP(ZBAR,IREAL8(NATNCS(K)))
      CALL FREHP(YBAR,IREAL8(NATNCS(K)))
      CALL FREHP(XBAR,IREAL8(NATNCS(K)))
      CALL FREHP(ZREF,IREAL8(NATNCS(K)))
      CALL FREHP(YREF,IREAL8(NATNCS(K)))
      CALL FREHP(XREF,IREAL8(NATNCS(K)))
      END DO
C
      RETURN
      END
      SUBROUTINE ENENC2(ENCS,POINTR,NUMEQV,NATNCS,
     & XREF,YREF,ZREF,XBAR,YBAR,ZBAR,XROT,YROT,ZROT,
     & RX,RY,RZ,U,WTNCS)
C
C Effective energy and derivatives for Non-Crystallographic
C Symmetry restraints
C Author: Bill Weis
C =================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      DOUBLE PRECISION ENCS
      INTEGER POINTR(*), NUMEQV, NATNCS
      DOUBLE PRECISION XREF(*),YREF(*),ZREF(*)
      DOUBLE PRECISION XBAR(*),YBAR(*),ZBAR(*)
      DOUBLE PRECISION XROT(*),YROT(*),ZROT(*)
      DOUBLE PRECISION RX(NATNCS,*),RY(NATNCS,*),RZ(NATNCS,*)
      DOUBLE PRECISION U(9,*),WTNCS
C local
      DOUBLE PRECISION UU(9)
      DOUBLE PRECISION DELX,DELY,DELZ,UDX,UDY,UDZ,S2,DFAC
      INTEGER M,J,L,EQPOS,IPNT
C parameter
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
      U(1,1)=ONE
      U(2,1)=ZERO
      U(3,1)=ZERO
      U(4,1)=ZERO
      U(5,1)=ONE
      U(6,1)=ZERO
      U(7,1)=ZERO
      U(8,1)=ZERO
      U(9,1)=ONE
C
C Compute 'average' coordinates: superpose
C NUMEQV-1 NCS equivalent sets onto the reference set.  Keep running
C computation of the average superposed coordinates.  Load transformed
C coordinate sets into 2-D arrays RX, RY, and RZ.  The fast
C index in these arrays is the counter for atoms within the equivalence
C group, slow is the equivalence set number.
C
C The first subgroup is the reference.
      DO M=1,NATNCS
      XREF(M)=X(POINTR(M))
      YREF(M)=Y(POINTR(M))
      ZREF(M)=Z(POINTR(M))
C
C Initialize averages
      XBAR(M)=XREF(M)
      YBAR(M)=YREF(M)
      ZBAR(M)=ZREF(M)
C
      RX(M,1)=XREF(M)
      RY(M,1)=YREF(M)
      RZ(M,1)=ZREF(M)
C
      END DO
C
      DO J=2,NUMEQV
      EQPOS=NATNCS*(J-1)
      DO M=1,NATNCS
      XROT(M)=X(POINTR(M+EQPOS))
      YROT(M)=Y(POINTR(M+EQPOS))
      ZROT(M)=Z(POINTR(M+EQPOS))
      END DO
C
      CALL ROTATE(XREF,YREF,ZREF,NATNCS,XROT,YROT,ZROT,UU,
     &            .FALSE.,1,J)
C
      DO M=1,NATNCS
      XBAR(M)=XBAR(M)+XROT(M)
      YBAR(M)=YBAR(M)+YROT(M)
      ZBAR(M)=ZBAR(M)+ZROT(M)
C
      RX(M,J)=XROT(M)
      RY(M,J)=YROT(M)
      RZ(M,J)=ZROT(M)
      END DO
C
C Store orthogonal rotation matrix U, indexed by equivalence group.
      DO L=1,9
      U(L,J)=UU(L)
      END DO
      END DO
C            equivalences
C
      DO M=1,NATNCS
      XBAR(M)=XBAR(M)/NUMEQV
      YBAR(M)=YBAR(M)/NUMEQV
      ZBAR(M)=ZBAR(M)/NUMEQV
      END DO
C
C Compute energy and derivatives
      DFAC=TWO*WTNCS
      DO J=1,NUMEQV
      EQPOS=NATNCS*(J-1)
C
C we're able to ignore any vector dependencies here
C$DIR NO_RECURRENCE
CDIR$ IVDEP
      DO M=1,NATNCS
      DELX=RX(M,J)-XBAR(M)
      DELY=RY(M,J)-YBAR(M)
      DELZ=RZ(M,J)-ZBAR(M)
      S2=DELX*DELX + DELY*DELY + DELZ*DELZ
      ENCS=ENCS + WTNCS*S2
C
C Transform direction of forces.
      UDX=U(1,J)*DELX + U(2,J)*DELY + U(3,J)*DELZ
      UDY=U(4,J)*DELX + U(5,J)*DELY + U(6,J)*DELZ
      UDZ=U(7,J)*DELX + U(8,J)*DELY + U(9,J)*DELZ
C
      IPNT=POINTR(M+EQPOS)
      DX(IPNT)=DX(IPNT) + DFAC*UDX
      DY(IPNT)=DY(IPNT) + DFAC*UDY
      DZ(IPNT)=DZ(IPNT) + DFAC*UDZ
      END DO
      END DO
C            equivalences
C
      RETURN
      END
