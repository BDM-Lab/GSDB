      SUBROUTINE XSYMPA(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &                  XRITSY,ARPNN,XRSYGP,XRSYIV)
C
C Routine parses XRAY SYMMetry statement.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      LOGICAL QHERM
      INTEGER XRSYTH, XRMSYM
      INTEGER XRNSYM, XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER ARPNN, XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C local
      INTEGER N
      LOGICAL ERR
      DOUBLE PRECISION DBTEMP
      DOUBLE COMPLEX DBCOMP
      INTEGER MAXWORD, ADDLEN, WORDLEN
      PARAMETER (MAXWORD=80)
      CHARACTER*(MAXWORD) WORD, ADDWORD
C begin
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SYMMetry=')
      CALL MISCOM('XRAY>',USED)
      IF (.NOT.USED) THEN
C======================================================================
      IF (WD(1:4).EQ.'RESE') THEN
      CALL XSYMRES(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &                   XRITSY,ARPNN)
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      CALL XRMAPR(0)
      QHERM=.TRUE.
C======================================================================
      ELSE IF (WD(1:1).EQ.'?') THEN
      WRITE(6,'(A)') ' | Symmetry operators:'
      WRITE(6,'(A)') ' | -------------------'
      CALL XSYPRI(XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
C======================================================================
      ELSE
      CALL SAVEWD
C
      CALL NEXTEX('SYMMetry=',WDD,WDMAX,WDDLEN)
C
C
C define symbol
      WORD='SYMMETRY_OP_'
      WORDLEN=12
      CALL ENCODI(XRNSYM+1,ADDWORD,MAXWORD,ADDLEN)
      CALL ADDST(WORD,MAXWORD,WORDLEN,ADDWORD,ADDLEN)
      CALL DECLAR( WORD(1:WORDLEN),'ST',WDD(1:WDDLEN), DBCOMP, DBTEMP )
C
C interpret symmetry operator and store
      CALL XRSYPA(WDD,WDDLEN,XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
C
      DO N=1,NPIG
      UPNBSS(N)=.TRUE.
      ENDDO
      CALL XRMAPR(0)
      END IF
C
      CALL NEXTWD('XRAY>')
C
      IF (WD(1:4).EQ.'SYMM') THEN
      DONE=.FALSE.
      ELSE
      DONE=.TRUE.
      CALL SAVEWD
      END IF
C
      END IF
      END DO
      DONE=.FALSE.
C
C check consistency of symmetry operators
      CALL XSYMCHK(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,ERR,
     &                   XRSYGP,XRSYIV)
      IF (ERR) THEN
      WRITE(6,'(A)') ' %ERR: Errors in symmetry operators.'
      END IF
C
C declare symbols $SYMMETRY and $SYMMETRY2
      DBTEMP=XRNSYM
      CALL DECLAR( 'SYMMETRY_OP_1', 'ST', '(X,Y,Z)', DBCOMP, DBTEMP )
      CALL DECLAR( 'SYMMETRY', 'DP', ' ', DBCOMP, DBTEMP )
      RETURN
      END
C======================================================================
      SUBROUTINE XRSYPA(ST,STLEN,XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
C
C Parses and interprets a symmetry operator, such as (-x,y+1/2,-z)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      CHARACTER*(*) ST
      INTEGER STLEN
      INTEGER XRMSYM
      INTEGER XRNSYM, XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3), XRSYTH
C local
      INTEGER TPLEN, FIELD, FIRST(3), LAST(3), IL, IFF, XYZ, IC1, IC2
      INTEGER TPMAX, I, J, II, III, DET
      PARAMETER (TPMAX=40)
      LOGICAL OKI, OKII, QID
      CHARACTER*(TPMAX) TP
      CHARACTER*1 XYZOP(3)
C begin
      IF (XRNSYM.GE.XRMSYM) THEN
      CALL WRNDIE(-5,'XRSYPA',
     & 'exceeded XRMSYM parameter --> recompile program')
      ELSE
      XRNSYM=XRNSYM+1
C
C now parse the symmetry operator string
      ERROR=.FALSE.
C
C make sure that the string contains no blanks AND EVERYTHING UPPERCASE
      TPLEN=0
      DO I=1,MIN(TPMAX,STLEN)
      IF (ST(I:I).NE.' ') THEN
      TPLEN=TPLEN+1
      TP(TPLEN:TPLEN)=ST(I:I)
      END IF
      END DO
C
C get comma separators
      IC1=INDEX(TP(1:TPLEN),',')
      IC2=INDEX(TP(IC1+1:TPLEN),',')
      IF (IC1.EQ.0.OR.IC2.EQ.0) THEN
      WRITE(6,'(A)') ' %XSYMM-ERR: missing comma(s)'
      ERROR=.TRUE.
      ELSE
      TP(IC1:IC1)=' '
      TP(IC1+IC2:IC1+IC2)=' '
      END IF
      IC2=IC2+IC1
C
C perpare string boundary list for outer loop
C first field
      FIRST(1)=2
      LAST(1)=IC1-1
C second field
      FIRST(2)=IC1+1
      LAST(2)=IC2-1
C third field
      FIRST(3)=IC2+1
      LAST(3)=TPLEN-1
C
C prepare x-y-z operator list inner loop
      XYZOP(1)='X'
      XYZOP(2)='Y'
      XYZOP(3)='Z'
C
C outer loop over string fields
      DO FIELD=1,3
C
      IFF=FIRST(FIELD)
      IL=LAST(FIELD)
C
C inner loop over basic operators X, Y, Z
      DO XYZ=1,3
      I=INDEX(TP(IFF:IL),XYZOP(XYZ))
      IF (I.GT.0) THEN
      TP(I+IFF-1:I+IFF-1)=' '
      IF (TP(I+IFF-2:I+IFF-2).EQ.'-') THEN
      TP(I+IFF-2:I+IFF-2)=' '
      XRSYMM(XRNSYM,FIELD,XYZ)=-1
      ELSEIF (TP(I+IFF-2:I+IFF-2).EQ.'+') THEN
      TP(I+IFF-2:I+IFF-2)=' '
      XRSYMM(XRNSYM,FIELD,XYZ)=1
      ELSE
      XRSYMM(XRNSYM,FIELD,XYZ)=+1
      END IF
      ELSE
      XRSYMM(XRNSYM,FIELD,XYZ)=0
      END IF
      END DO
C
C check translation operator
      I=INDEX(TP(IFF:IL),'/')
      IF (I.GT.1) THEN
      TP(IFF+I-1:IFF+I-1)=' '
      II=DECODI(TP(I+IFF-2:I+IFF-2),1,OKI)
      III=DECODI(TP(I+IFF:I+IFF),1,OKII)
      IF (.NOT.OKI.OR..NOT.OKII.OR.III.EQ.0) THEN
      WRITE(6,'(A)') ' SYMM-ERR: error in interpreting fraction'
      ERROR=.TRUE.
      ELSE
      TP(I+IFF-2:I+IFF-2)=' '
      TP(I+IFF:I+IFF)=' '
      IF (TP(I+IFF-3:I+IFF-3).EQ.'-') THEN
      TP(I+IFF-3:I+IFF-3)=' '
      XRSYMM(XRNSYM,FIELD,4)=-(II*XRSYTH)/III
      ELSEIF (TP(I+IFF-3:I+IFF-3).EQ.'+') THEN
      TP(I+IFF-3:I+IFF-3)=' '
      XRSYMM(XRNSYM,FIELD,4)=(II*XRSYTH)/III
      ELSE
      XRSYMM(XRNSYM,FIELD,4)=(II*XRSYTH)/III
      END IF
      END IF
      ELSE
      XRSYMM(XRNSYM,FIELD,4)=0
      END IF
C
      END DO
C
C check if any uninterpreted characters remain:
      DO I=2,TPLEN-1
      IF (TP(I:I).NE.' ') THEN
      WRITE(6,'(3A)')
     & ' XRSYPA: could not interpret character "',TP(I:I),
     & '" in symmetry expression.'
      ERROR=.TRUE.
      END IF
      END DO
      IF (ERROR) THEN
      CALL WRNDIE(-5,'XRSYPA',
     & 'incorrect symmetry expression.')
      END IF
C
C check determinant
      DET=(XRSYMM(XRNSYM,1,1)*XRSYMM(XRNSYM,2,2)
     &     -XRSYMM(XRNSYM,1,2)*XRSYMM(XRNSYM,2,1))*XRSYMM(XRNSYM,3,3)
     &   +(XRSYMM(XRNSYM,2,1)*XRSYMM(XRNSYM,3,2)
     &     -XRSYMM(XRNSYM,2,2)*XRSYMM(XRNSYM,3,1))*XRSYMM(XRNSYM,1,3)
     &   +(XRSYMM(XRNSYM,3,1)*XRSYMM(XRNSYM,1,2)
     &     -XRSYMM(XRNSYM,3,2)*XRSYMM(XRNSYM,1,1))*XRSYMM(XRNSYM,2,3)
      IF (ABS(DET).NE.1) THEN
      WRITE(6,'(A,I3)') ' %XRSYPA-ERR: invalid determinant =',DET
      ERROR=.TRUE.
      END IF
C
C check whether symmetry operator is ID
      QID=     XRSYMM(XRNSYM,1,1).EQ.1.AND.XRSYMM(XRNSYM,1,2).EQ.0
     &    .AND.XRSYMM(XRNSYM,1,3).EQ.0.AND.XRSYMM(XRNSYM,1,4).EQ.0
     &    .AND.XRSYMM(XRNSYM,2,1).EQ.0.AND.XRSYMM(XRNSYM,2,2).EQ.1
     &    .AND.XRSYMM(XRNSYM,2,3).EQ.0.AND.XRSYMM(XRNSYM,2,4).EQ.0
     &    .AND.XRSYMM(XRNSYM,3,1).EQ.0.AND.XRSYMM(XRNSYM,3,2).EQ.0
     &    .AND.XRSYMM(XRNSYM,3,3).EQ.1.AND.XRSYMM(XRNSYM,3,4).EQ.0
      IF (ERROR.OR.QID) THEN
      XRNSYM=XRNSYM-1
      ELSE
C
C compute transpose of the inverse of the symmetry operator
C (this is the symmetry operation for (h,k,l))
      XRITSY(XRNSYM,1,1)=XRSYMM(XRNSYM,2,2)*XRSYMM(XRNSYM,3,3)
     &                  -XRSYMM(XRNSYM,2,3)*XRSYMM(XRNSYM,3,2)
      XRITSY(XRNSYM,2,1)=XRSYMM(XRNSYM,3,2)*XRSYMM(XRNSYM,1,3)
     &                  -XRSYMM(XRNSYM,3,3)*XRSYMM(XRNSYM,1,2)
      XRITSY(XRNSYM,3,1)=XRSYMM(XRNSYM,1,2)*XRSYMM(XRNSYM,2,3)
     &                  -XRSYMM(XRNSYM,1,3)*XRSYMM(XRNSYM,2,2)
      XRITSY(XRNSYM,1,2)=XRSYMM(XRNSYM,2,3)*XRSYMM(XRNSYM,3,1)
     &                  -XRSYMM(XRNSYM,2,1)*XRSYMM(XRNSYM,3,3)
      XRITSY(XRNSYM,2,2)=XRSYMM(XRNSYM,3,3)*XRSYMM(XRNSYM,1,1)
     &                  -XRSYMM(XRNSYM,3,1)*XRSYMM(XRNSYM,1,3)
      XRITSY(XRNSYM,3,2)=XRSYMM(XRNSYM,1,3)*XRSYMM(XRNSYM,2,1)
     &                  -XRSYMM(XRNSYM,1,1)*XRSYMM(XRNSYM,2,3)
      XRITSY(XRNSYM,1,3)=XRSYMM(XRNSYM,2,1)*XRSYMM(XRNSYM,3,2)
     &                  -XRSYMM(XRNSYM,2,2)*XRSYMM(XRNSYM,3,1)
      XRITSY(XRNSYM,2,3)=XRSYMM(XRNSYM,3,1)*XRSYMM(XRNSYM,1,2)
     &                  -XRSYMM(XRNSYM,3,2)*XRSYMM(XRNSYM,1,1)
      XRITSY(XRNSYM,3,3)=XRSYMM(XRNSYM,1,1)*XRSYMM(XRNSYM,2,2)
     &                  -XRSYMM(XRNSYM,1,2)*XRSYMM(XRNSYM,2,1)
      DO I=1,3
      DO J=1,3
      XRITSY(XRNSYM,I,J)=XRITSY(XRNSYM,I,J)/DET
      END DO
      END DO
      END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XSYMCHK(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,ERR,
     &                   XRSYGP,XRSYIV)
C
C Does a variety of checks to make sure that the crystallographic
C symmetry operators are well-defined.
C
C Author: Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL ERR
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      INTEGER LSY(3,4), SYM, SYM2, SYM3, I, J, K, LSYM3, DET
      LOGICAL LQID, QID
      DOUBLE PRECISION DBREAL
      DOUBLE COMPLEX DBCOMP
C begin
C
      CALL DECLAR('CENTRIC_SPACEGROUP', 'LO', 'FALSE', DBCOMP, DBREAL)
C make sure that the symmetry operators form a group
      ERR=.FALSE.
      DO SYM=1,XRNSYM
C
C apply MOD operation to translational part to make
C sure that the translation is within -1,..,+1
C (this is the case when entering the operators from
C the International Tables).
      DO I=1,3
      XRSYMM(SYM,1,4)=MOD(XRSYMM(SYM,1,4)+10000*XRSYTH,XRSYTH)
      XRSYMM(SYM,2,4)=MOD(XRSYMM(SYM,2,4)+10000*XRSYTH,XRSYTH)
      XRSYMM(SYM,3,4)=MOD(XRSYMM(SYM,3,4)+10000*XRSYTH,XRSYTH)
      END DO
      END DO
C
C
      DO SYM=1,XRNSYM
C check determinant
      DET=(XRSYMM(SYM,1,1)*XRSYMM(SYM,2,2)
     &     -XRSYMM(SYM,1,2)*XRSYMM(SYM,2,1))*XRSYMM(SYM,3,3)
     &   +(XRSYMM(SYM,2,1)*XRSYMM(SYM,3,2)
     &     -XRSYMM(SYM,2,2)*XRSYMM(SYM,3,1))*XRSYMM(SYM,1,3)
     &   +(XRSYMM(SYM,3,1)*XRSYMM(SYM,1,2)
     &     -XRSYMM(SYM,3,2)*XRSYMM(SYM,1,1))*XRSYMM(SYM,2,3)
      IF (ABS(DET).NE.1) THEN
      WRITE(6,'(A,I3)') ' %XSYMCHK-ERR: invalid determinant =',DET
      ERR=.TRUE.
      END IF
      IF (DET.EQ.-1) THEN
      CALL DECLAR('CENTRIC_SPACEGROUP', 'LO', 'TRUE', DBCOMP, DBREAL)
      END IF
      END DO
C
C make sure that the first symmetry operator is the identity
C operator
      IF (.NOT.( XRSYMM(1,1,1).EQ.1.AND.XRSYMM(1,1,2).EQ.0
     &      .AND.XRSYMM(1,1,3).EQ.0.AND.XRSYMM(1,1,4).EQ.0
     &      .AND.XRSYMM(1,2,1).EQ.0.AND.XRSYMM(1,2,2).EQ.1
     &      .AND.XRSYMM(1,2,3).EQ.0.AND.XRSYMM(1,2,4).EQ.0
     &      .AND.XRSYMM(1,3,1).EQ.0.AND.XRSYMM(1,3,2).EQ.0
     &      .AND.XRSYMM(1,3,3).EQ.1.AND.XRSYMM(1,3,4).EQ.0 ) ) THEN
      ERR=.TRUE.
      WRITE(6,'(A)')
     & ' Error: first symmetry operator must be the identity operator.'
      END IF
C
      DO SYM=1,XRNSYM
C
      DO SYM2=1,XRNSYM
C
C check for duplications
      IF (SYM.NE.SYM2
     & .AND.XRSYMM(SYM,1,1).EQ.XRSYMM(SYM2,1,1)
     & .AND.XRSYMM(SYM,1,2).EQ.XRSYMM(SYM2,1,2)
     & .AND.XRSYMM(SYM,1,3).EQ.XRSYMM(SYM2,1,3)
     & .AND.MOD(XRSYMM(SYM,1,4)-XRSYMM(SYM2,1,4)
     &      +10000*XRSYTH,XRSYTH).EQ.0
     & .AND.XRSYMM(SYM,2,1).EQ.XRSYMM(SYM2,2,1)
     & .AND.XRSYMM(SYM,2,2).EQ.XRSYMM(SYM2,2,2)
     & .AND.XRSYMM(SYM,2,3).EQ.XRSYMM(SYM2,2,3)
     & .AND.MOD(XRSYMM(SYM,2,4)-XRSYMM(SYM2,2,4)
     &     +10000*XRSYTH,XRSYTH).EQ.0
     & .AND.XRSYMM(SYM,3,1).EQ.XRSYMM(SYM2,3,1)
     & .AND.XRSYMM(SYM,3,2).EQ.XRSYMM(SYM2,3,2)
     & .AND.XRSYMM(SYM,3,3).EQ.XRSYMM(SYM2,3,3)
     & .AND.MOD(XRSYMM(SYM,3,4)-XRSYMM(SYM2,3,4)
     &     +10000*XRSYTH,XRSYTH).EQ.0)
     & THEN
      ERR=.TRUE.
      WRITE(6,'(A,2I4)')
     & ' Error: Duplication of symmetry operators. Check #',
     &  SYM, SYM2
      END IF
C
C compute product:  SYM * SYM2
      DO I=1,3
      DO J=1,4
      LSY(I,J)=0
      END DO
      END DO
C
      DO I=1,3
      DO J=1,3
      DO K=1,3
      LSY(I,K)=LSY(I,K)+XRSYMM(SYM,I,J)*XRSYMM(SYM2,J,K)
      END DO
      END DO
      END DO
C
      DO I=1,3
      LSY(I,4)=XRSYMM(SYM,I,4)
      DO J=1,3
      LSY(I,4)=LSY(I,4)+XRSYMM(SYM,I,J)*XRSYMM(SYM2,J,4)
      END DO
      END DO
C
      DO I=1,3
      LSY(I,4)=MOD(LSY(I,4)+10000*XRSYTH,XRSYTH)
      END DO
C
      QID=.FALSE.
      DO SYM3=1,XRNSYM
      LQID=
     &      XRSYMM(SYM3,1,1).EQ.LSY(1,1)
     & .AND.XRSYMM(SYM3,1,2).EQ.LSY(1,2)
     & .AND.XRSYMM(SYM3,1,3).EQ.LSY(1,3)
     & .AND.XRSYMM(SYM3,1,4).EQ.LSY(1,4)
     & .AND.XRSYMM(SYM3,2,1).EQ.LSY(2,1)
     & .AND.XRSYMM(SYM3,2,2).EQ.LSY(2,2)
     & .AND.XRSYMM(SYM3,2,3).EQ.LSY(2,3)
     & .AND.XRSYMM(SYM3,2,4).EQ.LSY(2,4)
     & .AND.XRSYMM(SYM3,3,1).EQ.LSY(3,1)
     & .AND.XRSYMM(SYM3,3,2).EQ.LSY(3,2)
     & .AND.XRSYMM(SYM3,3,3).EQ.LSY(3,3)
     & .AND.XRSYMM(SYM3,3,4).EQ.LSY(3,4)
      QID=QID.OR.LQID
      IF (LQID) LSYM3=SYM3
      END DO
C
      IF (.NOT.QID) THEN
      ERR=.TRUE.
      WRITE(6,'(A,I4,A,I4)')
     & ' Error: symmetry operators do not form a group. Check #',
     &  SYM,' *  #',SYM2
      ELSE
C
C fill symmetry operator matrix
      XRSYGP(SYM,SYM2)=LSYM3
C
C fill inverse operation matrix
      IF (LSYM3.EQ.1) THEN
      XRSYIV(SYM)=SYM2
      END IF
      END IF
C
      END DO
      END DO
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'XMDOAS3','Error involving symmetry operators.')
      END IF
C
      IF (WRNLEV.GE.15) THEN
      WRITE(6,'(A)') ' '
      WRITE(6,'(2A)')
     & ' Symmetry operator matrix ',
     & '(id of op1, id of op2)-> (id of op1 * op2 )'
      DO SYM=1,XRNSYM
      DO SYM2=1,XRNSYM
      WRITE(6,'(3I6)') SYM, SYM2, XRSYGP(SYM,SYM2)
      END DO
      END DO
C
      WRITE(6,'(A)') ' '
      WRITE(6,'(A)')
     & ' Inverse symmetry operator (id of op) -> (id of inverse op)'
      DO SYM=1,XRNSYM
      WRITE(6,'(2I6)') SYM, XRSYIV(SYM)
      END DO
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSYMRES(QHERM,XRSYTH,XRMSYM,XRNSYM,XRSYMM,
     &                   XRITSY,ARPNN)
C
C Routine resets symmetry operator database
C
C Axel T. Brunger
      IMPLICIT NONE
C I/O
      LOGICAL QHERM
      INTEGER XRSYTH, XRMSYM
      INTEGER XRNSYM, XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER ARPNN
C local
      DOUBLE PRECISION DBTEMP
      DOUBLE COMPLEX DBCOMP
      INTEGER I, J
C begin
      XRNSYM=1
C
C declare symbols $SYMMETRY* and $SYMMETRY2*
      DBTEMP=XRNSYM
      CALL DECLAR( 'SYMMETRY', 'DP', ' ', DBCOMP, DBTEMP)
C
C initialize the first symmetry operator to identity
      DO I=1,3
      DO J=1,4
      XRSYMM(1,I,J)=0
      END DO
      XRSYMM(1,I,I)=1
      DO J=1,3
      XRITSY(1,I,J)=0
      END DO
      XRITSY(1,I,I)=1
      END DO
C
C reset asymmetric unit for maps
      ARPNN=0
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSYPRI(XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
C
C Routine prints the current space group symmetry operators
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER XRMSYM
      INTEGER XRNSYM, XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3), XRSYTH
C local
      INTEGER TPMAX
      PARAMETER (TPMAX=40)
      CHARACTER*(TPMAX) TP
      CHARACTER*1 XYZOP(3)
      INTEGER FIELD, XYZ, DENOM, FRACT, ISYM, TPLEN, I
      LOGICAL QFIRST
      CHARACTER*4 TEMP
C begin
C
C make list of symmetry operators
      XYZOP(1)='X'
      XYZOP(2)='Y'
      XYZOP(3)='Z'
      DO ISYM=1,XRNSYM
      CALL COPYST(TP,TPMAX,TPLEN,'(',1)
      DO FIELD=1,3
C
C check X-Y-Z operators
      QFIRST=.TRUE.
      DO XYZ=1,3
      IF (ABS(XRSYMM(ISYM,FIELD,XYZ)).EQ.1) THEN
      IF (XRSYMM(ISYM,FIELD,XYZ).EQ.+1) THEN
      IF (.NOT.QFIRST) CALL ADDST(TP,TPMAX,TPLEN,'+',1)
      CALL ADDST(TP,TPMAX,TPLEN,XYZOP(XYZ),1)
      ELSE
      CALL ADDST(TP,TPMAX,TPLEN,'-',1)
      CALL ADDST(TP,TPMAX,TPLEN,XYZOP(XYZ),1)
      END IF
      IF (QFIRST) QFIRST=.FALSE.
      END IF
      END DO
C
C check translation
      IF (XRSYMM(ISYM,FIELD,4).NE.0) THEN
      FRACT=ABS(XRSYMM(ISYM,FIELD,4))
      DENOM=XRSYTH
      DO I=8,2,-1
      IF (MOD(FRACT,I).EQ.0.AND.MOD(DENOM,I).EQ.0) THEN
      FRACT=FRACT/I
      DENOM=DENOM/I
      END IF
      END DO
      IF (XRSYMM(ISYM,FIELD,4).GT.0) THEN
      CALL ADDST(TP,TPMAX,TPLEN,'+',1)
      ELSE
      CALL ADDST(TP,TPMAX,TPLEN,'-',1)
      END IF
      CALL ENCODI(FRACT,TEMP,4,I)
      CALL ADDST(TP,TPMAX,TPLEN,TEMP,I)
      CALL ADDST(TP,TPMAX,TPLEN,'/',1)
      CALL ENCODI(DENOM,TEMP,4,I)
      CALL ADDST(TP,TPMAX,TPLEN,TEMP,I)
      END IF
      IF (FIELD.NE.3) CALL ADDST(TP,TPMAX,TPLEN,',',1)
      END DO
      CALL ADDST(TP,TPMAX,TPLEN,')',1)
      WRITE(6,'(2A)') ' | SYMMetry=',TP(1:TPLEN)
      END DO
      RETURN
      END
C
C
      SUBROUTINE XSYMCEL(XRINTR, XRNSYM, XRMSYM, XRSYMM)
C
C checks consistency of cell dimensions with crystallographic
C symmetry operators
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION  XRINTR(3,3)
      INTEGER           XRNSYM, XRMSYM
      INTEGER           XRSYMM(XRMSYM,3,4)
C local
      LOGICAL           ERR
      INTEGER           ISYM
      DOUBLE PRECISION  LENBASE, SLENBASE
      DOUBLE PRECISION   BASEF(3),  BASEC(3)
      DOUBLE PRECISION  SBASEF(3), SBASEC(3)
      DOUBLE PRECISION  MAXDELTA, DELTA
C parameter
      DOUBLE PRECISION  RELERR
      PARAMETER (RELERR = 1.0D-4)
C begin
      ERR = .FALSE.
C
      BASEF(1) = DFLOAT(13) / DFLOAT(100)
      BASEF(2) = DFLOAT(17) / DFLOAT(100)
      BASEF(3) = DFLOAT(19) / DFLOAT(100)
      CALL TRVEC3(XRINTR, BASEF, BASEC)
      LENBASE =   BASEC(1)**2
     &          + BASEC(2)**2
     &          + BASEC(3)**2
      MAXDELTA = LENBASE * RELERR
C
      DO ISYM = 1, XRNSYM
        SBASEF(1) =   XRSYMM(ISYM,1,1) * BASEF(1)
     &              + XRSYMM(ISYM,1,2) * BASEF(2)
     &              + XRSYMM(ISYM,1,3) * BASEF(3)
        SBASEF(2) =   XRSYMM(ISYM,2,1) * BASEF(1)
     &              + XRSYMM(ISYM,2,2) * BASEF(2)
     &              + XRSYMM(ISYM,2,3) * BASEF(3)
        SBASEF(3) =   XRSYMM(ISYM,3,1) * BASEF(1)
     &              + XRSYMM(ISYM,3,2) * BASEF(2)
     &              + XRSYMM(ISYM,3,3) * BASEF(3)
        CALL TRVEC3(XRINTR, SBASEF, SBASEC)
        SLENBASE =   SBASEC(1)**2
     &             + SBASEC(2)**2
     &             + SBASEC(3)**2
        DELTA = ABS(SLENBASE - LENBASE)
        IF (DELTA .GT. MAXDELTA) ERR = .TRUE.
      END DO
C
      IF (ERR) THEN
        WRITE(6,'(2A)') ' %XSYMCEL-ERR: check spacegroup and ',
     &      ' cell dimensions.'
        CALL WRNDIE(-5,'XSYMCEL',
     &   'Spacegroup is incompatible with cell dimensions.')
      END IF
C
      RETURN
      END
