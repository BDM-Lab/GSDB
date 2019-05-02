      SUBROUTINE NBDSET(OPTION)
C
C Nonbonded interaction (electrostatic / van der Waals ) parsing
C routine
C
C For syntax see main parsing loop "HELP"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'nbonds.inc'
      CHARACTER*(*) OPTION
C local
      INTEGER N
C parameter
      DOUBLE PRECISION ZERO,ONE,TWO,FOUR
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FOUR=4.0D0)
C begin
      IF (OPTION.EQ.'INIT') THEN
      CALL NBDINI
      ELSE IF (OPTION.EQ.'PARSE') THEN
C
C parsing
      CALL PUSEND('NBDSET>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('NBDSET>')
      CALL MISCOM('NBDSET>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-parameter-nbonds')
C
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',NBXTOL)
      ELSE IF (WD(1:4).EQ.'ATOM') THEN
      LGROUP=.FALSE.
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      LGROUP=.TRUE.
      ELSE IF (WD(1:4).EQ.'CDIE') THEN
      LCONS=.TRUE.
      ELSE IF (WD(1:4).EQ.'RDIE') THEN
      LCONS=.FALSE.
      ELSE IF (WD(1:4).EQ.'SWIT') THEN
      LSHFT=.FALSE.
      QTRUNC=.FALSE.
      ELSE IF (WD(1:4).EQ.'SHIF') THEN
      LSHFT=.TRUE.
      LSHFO=.FALSE.
      QTRUNC=.FALSE.
      ELSE IF (WD(1:4).EQ.'SHFO') THEN
      LSHFO=.TRUE.
      LSHFT=.FALSE.
      WRITE(6,'(A)')
     &' Electrostatic force shifting selected.'
      ELSE IF (WD(1:4).EQ.'BSHF') THEN
      CALL NEXTF('BSHFo=',BSHFO)
      ELSE IF (WD(1:4).EQ.'VSWI') THEN
      LVSHFT=.FALSE.
      QTRUNC=.FALSE.
      ELSE IF (WD(1:3).EQ.'EPS') THEN
      CALL NEXTF('EPS=',EPS)
      ELSE IF (WD(1:4).EQ.'E14F') THEN
      CALL NEXTF('E14Fac=',E14FAC)
      ELSE IF (WD(1:5).EQ.'CUTNB') THEN
      CALL NEXTF('CUTNB=',CUTNB)
      ELSE IF (WD(1:6).EQ.'CTOFNB') THEN
      CALL NEXTF('CTOFNB=',CTOFNB)
      ELSE IF (WD(1:6).EQ.'CTONNB') THEN
      CALL NEXTF('CTONNB=',CTONNB)
      ELSE IF (WD(1:4).EQ.'WMIN') THEN
      CALL NEXTF('WMIN=',WMIN)
      ELSE IF (WD(1:4).EQ.'NBXM') THEN
      CALL NEXTI('NBXMode=',NBXMOD)
      ELSE IF (WD(1:4).EQ.'TRUN') THEN
      QTRUNC=.TRUE.
      ELSE IF (WD(1:4).EQ.'REPE') THEN
      CALL NEXTF('REPEl-factor=',FREPEL)
      ELSE IF (WD(1:4).EQ.'INHI') THEN
      CALL NEXTF('INHIbit-factor=',INHIBT)
      ELSE IF (WD(1:4).EQ.'REXP') THEN
      CALL NEXTI('REXPonent=',EXPREP)
      ELSE IF (WD(1:4).EQ.'IREX') THEN
      CALL NEXTI('IREXponent=',EXIREP)
      ELSE IF (WD(1:4).EQ.'RCON') THEN
      CALL NEXTF('RCONstant=',CREP)
C
C added by JJK to support the attractive-repulsive term
C
      ELSE IF (WD(1:4).EQ.'RATT') THEN
           CALL NEXTF('RATTractive=', RATT)
      ELSE IF (WD(1:5).EQ.'DEPTH') THEN
           CALL NEXTF('DEPTHofwell=', DOW)
      ELSE IF (WD(1:4).EQ.'SPEC') THEN
           CALL NEXTF('SPECial=', NBQSPC)
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      CALL NBDPRI
      ELSE
      CALL CHKEND('NBDSET>',DONE)
      END IF
C
      DO N=1,NPIGMAX
C update exclusion list, 1-4 interaction list, and initialize others
      UPNBEX(N)=.TRUE.
C update intra-molecular interaction list
      UPNBLS(N)=.TRUE.
C update symmetry-related interaction list
      UPNBSS(N)=.TRUE.
      ENDDO
C
      END IF
      END DO
      DONE=.FALSE.
C
C check consistency of cutoffs
      IF (CTOFNB-R4SMAL.LT.CTONNB) THEN
      WRITE(6,'(A)')
     &' %NBDSET-ERR: inconsistent CTONNB, CTOFNB given. Reset.'
      CTONNB=CTOFNB-1.0D0
      END IF
C
      IF (CUTNB.LT.CTOFNB+TWO*NBXTOL-RSMALL.AND..NOT.QTRUNC.AND..NOT.
     &     FREPEL.GT.0.0D0) THEN
      WRITE(6,'(2A)')
     &' %NBDSET-ERR: inconsistent CUTNB, TOLErance, CTOFNB given. ',
     &' Reset.'
      CUTNB=CTOFNB+TWO*NBXTOL
      END IF
C
      IF (LSHFO.AND.(BSHFO.EQ.ZERO)) THEN
      WRITE(6,'(2A)')
     &' %NBDSET-ERR: inconsistent force shifting constant, BSHFO ',
     &' Reset to ONE.'
      BSHFO=ONE
      END IF
C
      ELSE IF (OPTION.EQ.'PRINT') THEN
      CALL NBDPRI
      END IF
      RETURN
      END
C===================================================================
      SUBROUTINE NBDINI
C
C routine initializes information about nonbonded interactions
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'cnst.inc'
C local
      INTEGER N
C begin
C
C=============================================================
C enforce initialization and update of all lists.
      DO N=1,NPIGMAX
      UPNBEX(N)=.TRUE.
      UPNBLS(N)=.TRUE.
      UPNBSS(N)=.TRUE.
      ENDDO
C enforce update of atom van der Waals parameter lookup table
      UPNBLK=.TRUE.
C=============================================================
C
C the non-bonded data structure is initialized -- set all heap
C pointers to zero.
      DO N=1,NPIGMAX
      IJNB(N)=0
      IJSB(N)=0
      IINB14(N)=0
      IIBLOE(N)=0
      IINBEX(N)=0
      LIST(1,N)=0
      SLIST(1,N)=0
      LIST14(1,N)=0
      END DO
C
C set defaults for parameters and options
C nonbonding lists
      DO N=1,NPIGMAX
      NNNB(N)=0
      NNBEX(N)=0
      NNB14(N)=0
      NNSB(N)=0
      ENDDO
C exclusion options
      NBXMOD=5
C update options
      CUTNB=8.5D0
      NBXTOL=0.5D0
      LGROUP=.FALSE.
      WMIN=1.5D0
C electrostatic options
      EPS=1.0D0
      E14FAC=1.0D0
      LCONS=.TRUE.
      LSHFT=.TRUE.
      LSHFO=.FALSE.
      BSHFO=1.0D0
C v d Waals options
      FREPEL=0.0D0
      CREP=100.0D0
      EXPREP=2
      EXIREP=2
      LVSHFT=.FALSE.
C
C attractive/repulsive term options
C added by JJK
C
      RATT = 0.0D0
      DOW = 0.0D0
C shifting/switching options
      CTONNB=6.5D0
      CTOFNB=7.5D0
C inhibition option
      INHIBT=0.25D0
C trunction option
      QTRUNC=.FALSE.
C
C special positions: exclude zero distance interactions due to spec.
C positions
      NBQSPC=0.1D0
      RETURN
      END
C===================================================================
      SUBROUTINE NBDPRI
C
C routine prints info about current nonbonded options
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
C local
      CHARACTER*5 SGROUP
      CHARACTER*11 SCONS
      CHARACTER*6 SSHFT
      CHARACTER*5 SSHFO
      CHARACTER*7 SVSHFT
C begin
      IF (LGROUP) THEN
      SGROUP='GROUp'
      ELSE
      SGROUP='ATOM '
      END IF
      IF (LCONS) THEN
      SCONS='CDIElectric'
      ELSE
      SCONS='RDIElectric'
      END IF
      IF (LSHFT.OR.LSHFO) THEN
      SSHFT='SHIFt'
      ELSE
      SSHFT='SWITch'
      END IF
      IF (LSHFO) THEN
      SSHFO='FORCE'
      ELSE
      SSHFO='POTEN'
      END IF
      IF (LVSHFT) THEN
      SVSHFT='VSHIft'
      ELSE
      SVSHFT='VSWItch'
      END IF
      IF (QTRUNC) THEN
      SSHFT='TRUNc'
      SVSHFT='TRUNc'
      END IF
C
      WRITE(6,'(2A)') ' -----nonbonded-list-options-----------',
     & '--------------------'
      WRITE(6,'(A,F8.3,A,F8.3,A,F8.3,3A)')
     & ' | CUTNb=',CUTNB,' TOLErance=',NBXTOL,' WMIN=',WMIN,' ',SGROUP,
     & '  |'
      WRITE(6,'(A,F8.3,A)')
     & ' | INHIbit=',INHIBT,
     & '                                       |'
      IF (FREPEL.GT.0.0D0) THEN
      WRITE(6,'(2A)') ' -----no electrostatic-----------------',
     & '--------------------'
      WRITE(6,'(2A)') ' ------------- repulsive potential ----',
     & '--------------------'
      WRITE(6,'(A,F8.3,A,F8.3,A,I4,A)')
     & ' | REPEl',FREPEL,' RCONst=',CREP,' REXPonent=',EXPREP,
     & '           |'
      WRITE(6,'(A,I4,A)')
     & ' | IREXponent=',EXIREP,
     & '                                        |'
      ELSE
      WRITE(6,'(2A)') ' -----electrostatic options------------',
     & '--------------------'
      WRITE(6,'(A,F8.3,A,F8.3,7A)')
     & ' | EPS=',EPS,' E14Fac=',E14FAC,' ',SCONS,' ',SSHFO,' ',
     & SSHFT,'  |'
      IF (LSHFO) THEN
      WRITE(6,'(A,F8.3,2A)')
     & ' | BSHForce=',BSHFO,'                                ',
     & '      |'
      END IF
      WRITE(6,'(2A)') ' -----van der Waals options------------',
     & '--------------------'
      WRITE(6,'(3A)') ' | ',SVSHFT,
     & '                                                |'
      IF (CTONNB.LT.99.0D0) THEN
      WRITE(6,'(2A)') ' -----switching /shifting parameters---',
     & '--------------------'
      WRITE(6,'(A,F8.3,A,F8.3,A)')
     & ' | CTONNB=',CTONNB,' CTOFNB=',CTOFNB,
     & '                        |'
      END IF
      END IF
      WRITE(6,'(2A)') ' -----exclusion list options-----------',
     & '--------------------'
      WRITE(6,'(A,I4,A)')
     & ' | NBXMOD=',NBXMOD,
     & '                                            |'
      WRITE(6,'(2A)') ' --------------------------------------',
     & '--------------------'
C
      RETURN
      END
C=====================================================================
      SUBROUTINE RENBND(L1,L2,L3,L4,LXYZ,L5,L6,L7,L8,N)
C
C resizes nonbonded lists:
C  positive entry: free list and allocate new space
C  zero entry: just free list
C  negative entry: don't touch it
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'funct.inc'
      INTEGER L1, L2, L3, L4, LXYZ, L5, L6, L7, L8
      INTEGER N
C begin
C intra-molecular interaction list
      IF (L1.GE.0) THEN
      IF (IJNB(N).NE.0) THEN
      CALL FREHP(IJNB(N),INTEG4(XJNB(N)))
      IJNB(N)=0
      XJNB(N)=0
      END IF
      END IF
C
      IF (L2.GE.0) THEN
      IF (LIST(1,N).NE.0) THEN
      CALL FREENB(LIST(1,N))
      END IF
      END IF
C
C symmetry-image interaction list
      IF (L3.GE.0) THEN
      IF (IJSB(N).NE.0) THEN
      CALL FREHP(IJSB(N),INTEG4(XJSB(N)))
      IJSB(N)=0
      XJSB(N)=0
      END IF
      END IF
C
      IF (L4.GE.0) THEN
      IF (SLIST(1,N).NE.0) THEN
      CALL FREENB(SLIST(1,N))
      END IF
      END IF
C
C list of last atomic positions
      IF (INBX.NE.0.AND.LXYZ.GE.0) THEN
      CALL FREHP(INBZ,IREAL8(XNBX))
      CALL FREHP(INBY,IREAL8(XNBX))
      CALL FREHP(INBX,IREAL8(XNBX))
      INBX=0
      INBY=0
      INBZ=0
      XNBX=0
      END IF
C
C 1-4 interaction list
      IF (IINB14(N).NE.0.AND.L5.GE.0) THEN
      CALL FREHP(IINB14(N),INTEG4(XINB14(N)))
      IINB14(N)=0
      XINB14(N)=0
      END IF
C
      IF (L6.GE.0) THEN
      IF (LIST14(1,N).NE.0) THEN
      CALL FREENB(LIST14(1,N))
      END IF
      END IF
C
C exclusion list
      IF (IIBLOE(N).NE.0.AND.L7.GE.0) THEN
      CALL FREHP(IIBLOE(N),INTEG4(XIBLOE(N)))
      IIBLOE(N)=0
      XIBLOE(N)=0
      END IF
      IF (IINBEX(N).NE.0.AND.L8.GE.0) THEN
      CALL FREHP(IINBEX(N),INTEG4(XINBEX(N)))
      IINBEX(N)=0
      XINBEX(N)=0
      END IF
C
      IF (L8.GT.0) THEN
      XINBEX(N)=L8
      IINBEX(N)=ALLHP(INTEG4(L8))
      END IF
      IF (L7.GT.0) THEN
      XIBLOE(N)=L7
      IIBLOE(N)=ALLHP(INTEG4(L7))
      END IF
      IF (L6.GT.0) THEN
      CALL ININB(LIST14(1,N),L6)
      END IF
      IF (L5.GT.0) THEN
      XINB14(N)=L5
      IINB14(N)=ALLHP(INTEG4(L5))
      END IF
      IF (LXYZ.GT.0) THEN
      XNBX=LXYZ
      INBX=ALLHP(IREAL8(LXYZ))
      INBY=ALLHP(IREAL8(LXYZ))
      INBZ=ALLHP(IREAL8(LXYZ))
      END IF
      IF (L4.GT.0) THEN
      CALL ININB(SLIST(1,N),L4)
      END IF
      IF (L3.GT.0) THEN
      XJSB(N)=L3
      IJSB(N)=ALLHP(INTEG4(L3))
      END IF
      IF (L2.GT.0) THEN
      CALL ININB(LIST(1,N),L2)
      END IF
      IF (L1.GT.0) THEN
      XJNB(N)=L1
      IJNB(N)=ALLHP(INTEG4(L1))
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE PUTNB(FIRST,LAST,PIECE,LIST)
C
C Routine copies the FIRST...LAST elements of the array PIECE
C into the linked LIST starting at the current position of the LIST.
C The current position of the LIST is then reset to be after the N
C copied elements. The current position is reset by calling subroutine
C RESNB.  Prior to the first time the LIST is used it has to
C be initialized by subroutine ININB.
C
C The contents of LIST actually resides on the HEAP as a linked
C list of arrays of length LIST(4).  The last element of each array
C points to the next array.  Termination of the list is indicated by
C a zero pointer.
C
C LIST(1): heap pointer to first array
C LIST(2): heap pointer to current array
C LIST(3): position in current array
C LIST(4): fragmentation, that is, number of elements in each array
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER FIRST, LAST, PIECE(*), LIST(4)
C local
      INTEGER ILAST
      INTEGER LFIRST, LLAST
C pointer
      INTEGER NEWPTR
C begin
      LFIRST=FIRST
      LLAST=LAST-1
      DO WHILE (LLAST.NE.LAST)
      ILAST=MIN(LIST(4),LIST(3)+LAST-LFIRST)
      LLAST=LFIRST+(ILAST-LIST(3))
      CALL COPYIS(PIECE,LFIRST,LLAST,HEAP(LIST(2)),LIST(3))
      LFIRST=LLAST+1
      LIST(3)=ILAST+1
      IF (LLAST.NE.LAST) THEN
C need more space
C
C first see whether there is already something linked to the
C current array (if so, the heap pointer for the next array
C is stored in element LIST(4)+1 of the current array).
      CALL COPYIS(HEAP(LIST(2)),LIST(4)+1,LIST(4)+1,NEWPTR,1)
      IF (NEWPTR.EQ.0) THEN
C have to allocate another array
      NEWPTR=ALLHP(INTEG4(LIST(4)+1))
C link it to the current array
      CALL COPYIS(NEWPTR,1,1,HEAP(LIST(2)),LIST(4)+1)
C terminate the new array
      CALL COPYIS(0,1,1,HEAP(NEWPTR),LIST(4)+1)
      END IF
C set the current heap pointer to the new array
      LIST(2)=NEWPTR
C reset the position variable
      LIST(3)=1
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE GETNB(N,PIECE,LIST)
C
C Routine copies the next N elements from the linked LIST into
C the array PIECE starting at the current position of the LIST.
C The current position of the LIST is then reset to be after the N
C copied elements. The current position is reset by calling subroutine
C RESNB.
C
C The contents of LIST actually resides on the HEAP as a linked
C list of arrays of length LIST(4).  The last element of each array
C points to the next array.  Termination of the list is indicated by
C a zero pointer.
C
C LIST(1): heap pointer to first array
C LIST(2): heap pointer to current array
C LIST(3): position in current array
C LIST(4): fragmentation, that is, number of elements in each array
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INTEGER N, PIECE(*), LIST(4)
C local
      INTEGER ILAST
      INTEGER LFIRST, LLAST
C pointer
      INTEGER NEWPTR
C begin
      LFIRST=1
      LLAST=N-1
      DO WHILE (LLAST.NE.N)
      ILAST=MIN(LIST(4),LIST(3)+N-LFIRST)
      LLAST=LFIRST+(ILAST-LIST(3))
      CALL COPYIS(HEAP(LIST(2)),LIST(3),ILAST,PIECE,LFIRST)
      LFIRST=LLAST+1
      LIST(3)=ILAST+1
      IF (LLAST.NE.N) THEN
C get next fragment from linked list
      CALL COPYIS(HEAP(LIST(2)),LIST(4)+1,LIST(4)+1,NEWPTR,1)
      IF (NEWPTR.EQ.0) THEN
C
      WRITE(6,'(A)')
     & ' Fatal error: reading beyond end of nonbonded LIST'
      CALL WRNDIE(-5,'NBGET','fatal coding error')
      END IF
C set the current heap pointer to the new array
      LIST(2)=NEWPTR
C reset the position variable
      LIST(3)=1
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE RESNB(LIST)
C
C Routine resets the current heap pointer and the current position.
C
C LIST(1): heap pointer to first array
C LIST(2): heap pointer to current array
C LIST(3): position in current array
C LIST(4): fragmentation, that is, number of elements in each array
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
C
C I/O
      INTEGER LIST(4)
C begin
C
C set current heap pointer to LIST(1)
      LIST(2)=LIST(1)
C
C reset the position variable
      LIST(3)=1
      RETURN
      END
C======================================================================
      SUBROUTINE ININB(LIST,NFRAG)
C
C Routine allocates the first array for the linked LIST.
C
C LIST(1): heap pointer to first array
C LIST(2): heap pointer to current array
C LIST(3): position in current array
C LIST(4): fragmentation, that is, number of elements in each array
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C I/O
      INTEGER LIST(4), NFRAG
C begin
      LIST(4)=NFRAG
C allocate an array of size NFRAG+1
      LIST(1)=ALLHP(INTEG4(NFRAG+1))
C the (NFRAG+1)th element in this list is the heap pointer to the next
C fragment of the linked list.  Since we're only allocating
C one fragment here, we set this pointer to zero.
      CALL COPYIS(0,1,1,HEAP(LIST(1)),NFRAG+1)
C reset the current heap counter and position
      CALL RESNB(LIST)
      RETURN
      END
C======================================================================
      SUBROUTINE FREENB(LIST)
C
C Routine recursively frees up all heap space associated with the
C linked LIST.
C
C LIST(1): heap pointer to first array
C LIST(2): heap pointer to current array
C LIST(3): position in current array
C LIST(4): fragmentation, that is, number of elements in each array
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER LIST(4)
C local
      INTEGER NEWPTR
C begin
      LIST(2)=LIST(1)
      DO WHILE (LIST(2).NE.0)
C get the heap pointer to the next fragment of the linked list
      CALL COPYIS(HEAP(LIST(2)),LIST(4)+1,LIST(4)+1,NEWPTR,1)
C free the current fragment
      CALL FREHP(LIST(2),INTEG4(LIST(4)+1))
      LIST(2)=NEWPTR
      END DO
C
C set the pointer to zero.
      LIST(1)=0
      RETURN
      END
C======================================================================
      SUBROUTINE COPYIS(A,ASTART,AEND,COPY,CSTART)
C
C  Copies array integer A[ASTART,AEND] into COPY[CSTART:]
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER A(*), ASTART, AEND, COPY(*), CSTART
C local
      INTEGER I, II
C begin
      IF (AEND.GE.ASTART) THEN
      II=CSTART
      DO I=ASTART,AEND
      COPY(II)=A(I)
      II=II+1
      END DO
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XPSEAR(N,NNPIG)
C
C Routine searches for nonbonded interactions.
C
C QSYMM vs QNCS mode (exclusive options):
C  QSYMM mode uses crystallographic symmetry operators to generate
C        images of the primary molecule stored in the structure file.
C        All symmetry related molecules are generated and the
C        distance vectors to the primary molecule are computed.
C        Then the minimum image convention is applied to all distance
C        vectors in the following way: The distance vector is
C        converted to fractional coordinates, then the standard
C        minimum image convention is applied, and then the distance
C        vector is converted back to orthogonal A coordinates.
C        Then the distance cutoff criterion is applied to the
C        distance vector.  The information about symmetry
C        operators and the unit-cell is communicated via xcrystal.inc.
C
C  QNCS mode uses non-crystallographic symmetry operators to generate
C       images of the primary molecule stored in the structure file.
C        All symmetry related molecules are generated and the
C        distance vectors to the primary molecule are computed.
C        Then the distance cutoff criterion is applied to the
C        distance vector.  The information about symmetry
C        operators is communicated via ncs.inc.
C
C The routine splits the interactions into intra-molecular interactions
C (JNB, LIST) and symmetry interactions (JSB, SLIST).  These lists
C correspond to the energy terms ELEC,VDW and PELE,PVDW, respectively.
C Depending on which energy terms are turned on, the appropriate lists
C are generated.  Only for intra-molecular interactions, the
C non-bonded exclusion list is checked.
C
C Intra-molecular interactions between primary molecule atoms
C are contained as "trivial" cases in both QSYMM and QNCS mode:
C the first symmetry operator is always the identity.
C The case that only intra-molecular interactions are required
C (PELE, PVDW false) is treated as a special QNCS case.
C
C GROUP vs. ATOM mode (exclusive options)
C  GROUP mode excludes or includes whole groups as specified in the
C  structure file.  Distances are checked between group centroids (after
C  applying the minimum image convention in QSYMM mode).  The standard
C  atom-atom interaction lists are generated for allowed group-group
C  interactions.  Then, the minimum image convention is applied to
C  individual atom-atom interactions, and they are put either in the
C  list of intra-molecular interactions (JNB, INBLO) or the list of symmetry
C  related interactions (JSB, ISBLO).  The minimum image convention
C  is used again during the actual energy calculation on an atom-by
C  atom basis.  This can cause problems for situations where one of the
C  unit-cell lenghts is less than twice the non-bonded cutoff.  There
C  is no work-around for this problem yet.
C
C ATOM mode excludes or includes individual atom-atom interactions.
C  This mode uses the groups to make the calculation more efficient
C  (rectangles around group-centroid criterion).   Note, that the
C  results are independent of the size of the groups!
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INTEGER N, NNPIG
C local
      LOGICAL QSYMM, QINTRA, QNCS
      INTEGER NLEN, NSLEN
C pointers
      INTEGER QMIM
      INTEGER QTEST, QSTEST, XTEMP, YTEMP, ZTEMP, XF, YF, ZF
      INTEGER XJF, YJF, ZJF, INTERJ
      INTEGER XCENT, YCENT, ZCENT, XCENTF, YCENTF, ZCENTF
      INTEGER XSIZ, YSIZ, ZSIZ, SPOINT, INTERG, INTERC, ILOCAL
C begin
      IF (NATOM.GT.0) THEN
C
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: at entry')
C
C
C Determine whether we have to compute the intra or
C symmetry interactions (or
C non-crystallographic symmetry interactions).  If we are using strict
C non-crystallographic symmetry, and the packing terms are on, then
C crystallographic symmetry is automatically shut off, and
C the packing energies used for NCS related non-bonded interactions.
      IF (.NOT. LNCSST) THEN
      QSYMM=QENER(SSPVDW).OR.QENER(SSPELE)
      QNCS=.FALSE.
      ELSE
      QSYMM=.FALSE.
      QNCS=QENER(SSPVDW).OR.QENER(SSPELE)
      END IF
      QINTRA=QENER(SSVDW).OR.QENER(SSELEC)
C
C Determine lengths of some of the arrays
      IF (QINTRA) THEN
      NLEN=NATOM
      ELSE
      NLEN=0
      END IF
      IF (QSYMM .OR. QNCS) THEN
C Note that xtal packing and ncs images are mutually exclusive
      IF (QSYMM) THEN
      NSLEN=NATOM*XRNSYM
      ELSE
      NSLEN=NATOM*NNSYM
      END IF
      ELSE
      NSLEN=0
      END IF
C
C Initialization.
C ==============
      IF (UPNBEX(N)) THEN
      CALL RENBND(0,0,0,0,0,0,0,0,0,N)
C
      IF (QINTRA) THEN
C update exclusion list.
      WRITE(6,'(A,I2)')
     &' NBONDS: generating intra-molecular exclusion list with mode=',
     & NBXMOD
      CALL UPINB(N)
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: after excl. ')
      END IF
C allocate space for various lists (the fragmentation of the
C linked lists is NATOM).  Note that
C for the last Pair of Interacting Groups (PIG), allocate the
C 'previous position' arrays.
      IF (N.EQ.NNPIG) THEN
      CALL RENBND(NLEN,NATOM,NSLEN,NATOM,NATOM,-1,-1,-1,-1,N)
      ELSE
      CALL RENBND(NLEN,NATOM,NSLEN,NATOM,-1,-1,-1,-1,-1,N)
      END IF
      END IF
C
C allocate temporary space for search routine
      QTEST=ALLHP(ILOGIC(NATOM))
      QSTEST=ALLHP(ILOGIC(NATOM))
      XTEMP=ALLHP(IREAL8(NATOM))
      YTEMP=ALLHP(IREAL8(NATOM))
      ZTEMP=ALLHP(IREAL8(NATOM))
      QMIM=ALLHP(INTEG4(NATOM))
      XF=ALLHP(IREAL8(NATOM))
      YF=ALLHP(IREAL8(NATOM))
      ZF=ALLHP(IREAL8(NATOM))
      XJF=ALLHP(IREAL8(NATOM))
      YJF=ALLHP(IREAL8(NATOM))
      ZJF=ALLHP(IREAL8(NATOM))
      INTERJ=ALLHP(INTEG4(NATOM))
      XCENT=ALLHP(IREAL8(NGRP))
      YCENT=ALLHP(IREAL8(NGRP))
      ZCENT=ALLHP(IREAL8(NGRP))
      XCENTF=ALLHP(IREAL8(NGRP))
      YCENTF=ALLHP(IREAL8(NGRP))
      ZCENTF=ALLHP(IREAL8(NGRP))
      XSIZ=ALLHP(IREAL8(NGRP))
      YSIZ=ALLHP(IREAL8(NGRP))
      ZSIZ=ALLHP(IREAL8(NGRP))
      SPOINT=ALLHP(INTEG4(NATOM))
      INTERG=ALLHP(INTEG4(NGRP))
      INTERC=ALLHP(INTEG4(NGRP))
      ILOCAL=ALLHP(INTEG4(NATOM))
C
      CALL XPSEA2(QSYMM,QINTRA,QNCS,HEAP(IJNB(N)),
     &            LIST(1,N),HEAP(IJSB(N)),NATOM,
     &            SLIST(1,N),HEAP(IIBLOE(N)),
     &            HEAP(IINBEX(N)),HEAP(QTEST),HEAP(QSTEST),
     &            HEAP(XTEMP),HEAP(YTEMP),HEAP(ZTEMP),HEAP(QMIM),
     &            HEAP(XF),HEAP(YF),HEAP(ZF),HEAP(XJF),HEAP(YJF),
     &            HEAP(ZJF),HEAP(INTERJ),
     &            HEAP(XCENT),HEAP(YCENT),HEAP(ZCENT),
     &            HEAP(XCENTF),HEAP(YCENTF),HEAP(ZCENTF),
     &            HEAP(XSIZ),HEAP(YSIZ),HEAP(ZSIZ),HEAP(SPOINT),
     &            HEAP(INTERG),HEAP(INTERC),
     &            HEAP(IINTER(N)),NNNB(N),NNSB(N),
     &            CUTNB,WMIN,LGROUP,HEAP(ILOCAL),NBQSPC)
C
C free temporary space for search routine
      CALL FREHP(ILOCAL,INTEG4(NATOM))
      CALL FREHP(INTERC,INTEG4(NGRP))
      CALL FREHP(INTERG,INTEG4(NGRP))
      CALL FREHP(SPOINT,INTEG4(NATOM))
      CALL FREHP(ZSIZ,IREAL8(NGRP))
      CALL FREHP(YSIZ,IREAL8(NGRP))
      CALL FREHP(XSIZ,IREAL8(NGRP))
      CALL FREHP(ZCENTF,IREAL8(NGRP))
      CALL FREHP(YCENTF,IREAL8(NGRP))
      CALL FREHP(XCENTF,IREAL8(NGRP))
      CALL FREHP(ZCENT,IREAL8(NGRP))
      CALL FREHP(YCENT,IREAL8(NGRP))
      CALL FREHP(XCENT,IREAL8(NGRP))
      CALL FREHP(INTERJ,INTEG4(NATOM))
      CALL FREHP(ZJF,IREAL8(NATOM))
      CALL FREHP(YJF,IREAL8(NATOM))
      CALL FREHP(XJF,IREAL8(NATOM))
      CALL FREHP(ZF,IREAL8(NATOM))
      CALL FREHP(YF,IREAL8(NATOM))
      CALL FREHP(XF,IREAL8(NATOM))
      CALL FREHP(QMIM,INTEG4(NATOM))
      CALL FREHP(ZTEMP,IREAL8(NATOM))
      CALL FREHP(YTEMP,IREAL8(NATOM))
      CALL FREHP(XTEMP,IREAL8(NATOM))
      CALL FREHP(QSTEST,ILOGIC(NATOM))
      CALL FREHP(QTEST,ILOGIC(NATOM))
C
      IF (TIMER.GT.1) CALL DSPCPU(' NBONDS: at exit')
C
C now update INBX, INBY, INBZ for next tolerance check
      IF (N.EQ.NNPIG)
     &  CALL COPYCO(NATOM,X,Y,Z,HEAP(INBX),HEAP(INBY),HEAP(INBZ))
C
      END IF
      RETURN
      END
C
      SUBROUTINE XPSEA2(QSYMM,QINTRA,QNCS,JNB,LIST,
     &            JSB,LEN,SLIST,IBLOEX,INBEX,
     &            QTEST,QSTEST,XTEMP,YTEMP,ZTEMP,QMIM,
     &            XF,YF,ZF,XJF,YJF,ZJF,
     &            INTERJ,XCENT,YCENT,
     &            ZCENT,XCENTF,YCENTF,ZCENTF,XSIZ,YSIZ,ZSIZ,SPOINT,
     &            INTERG,INTERC,
     &            INTERE,NNNB,NNSB,CUTNB,WMIN,LGROUP,
     &            ILOCAL,NBQSPC)
C
C see routine above
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      LOGICAL QSYMM, QINTRA, QNCS
      INTEGER JNB(*), LIST(4)
      INTEGER LEN, JSB(LEN, *), SLIST(4)
      INTEGER INBEX(*), IBLOEX(*)
      LOGICAL QTEST(*), QSTEST(*)
      DOUBLE PRECISION XTEMP(*), YTEMP(*), ZTEMP(*)
      INTEGER QMIM(*)
      DOUBLE PRECISION XF(*), YF(*), ZF(*), XJF(*), YJF(*), ZJF(*)
      INTEGER INTERJ(*)
      DOUBLE PRECISION XCENT(*), YCENT(*), ZCENT(*)
      DOUBLE PRECISION XCENTF(*), YCENTF(*), ZCENTF(*)
      DOUBLE PRECISION XSIZ(*), YSIZ(*), ZSIZ(*)
      INTEGER SPOINT(*), INTERG(*), INTERC(*)
      INTEGER INTERE(*),NNNB,NNSB
      DOUBLE PRECISION CUTNB,WMIN
      LOGICAL LGROUP
      INTEGER ILOCAL(*)
      DOUBLE PRECISION NBQSPC
C local
      DOUBLE PRECISION CUTNB2, WMIN2, DBPREC
      DOUBLE COMPLEX DBCOMP
      LOGICAL CHCK
      INTEGER ISYM, I, IAT, J, JJ, NBSYM, NXI, NXIMAX
      INTEGER JSYM, NNIAT, III, NVIOL
      DOUBLE PRECISION NBQSPC2
C parameter
      DOUBLE PRECISION ZERO, EEPS, PADDNG, ONE, TWO
      PARAMETER (ZERO=0.0D0, EEPS=0.1D0, PADDNG=3)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
C
      NVIOL=0
      WMIN2=WMIN**2
      CUTNB2=CUTNB**2
      NBQSPC2=NBQSPC**2
C
C set initial number of interactions to zero
      NNSB=0
      NNNB=0
C
C reset the nonbonded linked list pointers
      CALL RESNB(LIST)
      CALL RESNB(SLIST)
C
C compute group centers, size and interaction criterium
C for all groups of asymmetric unit
      CALL XPGRUP(X,Y,Z,XCENT,YCENT,ZCENT,XSIZ,YSIZ,ZSIZ,INTERG,INTERC,
     &            INTERE)
C
C
C get number of symmetry operators which we want to use.
      IF (.NOT.QSYMM.AND..NOT.QNCS) THEN
C
C no symmetry interactions required: just use the first symmetry operator
C which is the identity.
      JSYM=1
      ELSE IF (QNCS) THEN
C
C NCS symmetry
      JSYM=NNSYM
      ELSE
C
C crystallographic symmetry
      JSYM=XRNSYM
      END IF
C
      DO ISYM=1,JSYM
C
      IF (QSYMM) THEN
C
C crystallographic-symmetry-mode:
C we have to apply
C the symmetry operators to the j-set atoms and groups.  --> return in
C XF, YF, ZF, XCENTF, YCENTF, ZCENTF (fractional coordinates!).
      CALL XPSYMO(NATOM,X,Y,Z,ISYM,XF,YF,ZF)
      CALL XPSYMO(NGRP,XCENT,YCENT,ZCENT,ISYM,XCENTF,YCENTF,ZCENTF)
C
      ELSE
C
C NCS-symmetry mode (includes pure-intramolecular mode):
C we have to apply the NCS symmetry operator
C to the j-set atoms and groups.  --> return in XF, YF, ZF, XCENTF,
C YCENTF, ZCENTF (orthogonal coordinates!).
      CALL NCSXYZ(NATOM,X,Y,Z,ISYM,XF,YF,ZF)
      CALL NCSXYZ(NGRP,XCENT,YCENT,ZCENT,ISYM,XCENTF,YCENTF,ZCENTF)
      END IF
C
C
C loop over i groups (we only have to search groups j=i,...,ngrp
C because of symmetry, i.e., if G[i] and S G{j} are interacting
C then by application of another symmetry operator we must have
C S' G[i] and G[j] interacting.
      DO I=1,NGRP
C
C check whether group contains atoms that interact at all
      IF (ABS(INTERG(I)).LT.+9999) THEN
C
      IF (QSYMM) THEN
C
C crystallographic symmetry mode:
C compute minimum image distance between groups i,j
C
C
      CALL XPIMAG(XCENT(I),YCENT(I),ZCENT(I),I,NGRP,
     &            XCENTF,YCENTF,ZCENTF,XTEMP,YTEMP,ZTEMP,.FALSE.,QMIM)
C
      ELSE
C
C NCS symmetry mode (includes pure intra-molecular mode):
C compute direct distance between groups i, j.
      DO J=I,NGRP
      XTEMP(J)=XCENT(I)-XCENTF(J)
      YTEMP(J)=YCENT(I)-YCENTF(J)
      ZTEMP(J)=ZCENT(I)-ZCENTF(J)
      END DO
      END IF
C
      IF (LGROUP) THEN
C
C in group-by-group-cutoff mode compute distances between group centers
      DO J=I,NGRP
      XJF(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
      ELSE
C
C in atom-by-atom-cutoff mode
C compute minimum distances between solid rectangles
      DO J=I,NGRP
      XJF(J)= MAX(ZERO,ABS(XTEMP(J))-XSIZ(I)-XSIZ(J))**2
     &     +  MAX(ZERO,ABS(YTEMP(J))-YSIZ(I)-YSIZ(J))**2
     &     +  MAX(ZERO,ABS(ZTEMP(J))-ZSIZ(I)-ZSIZ(J))**2
      END DO
      END IF
C
C check cutoff criterion and interaction criterion between groups i,j
      DO J=I,NGRP
      QSTEST(J)=XJF(J).LT.CUTNB2 .AND.
     &         ABS(INTERG(I)+INTERG(J)) .LT. INTERC(I)+INTERC(J)
      END DO
C
C set QSTEST(I) always to true (this is needed below)
      QSTEST(I)=.TRUE.
C
C now that we know all groups j which are interacting with group i
C we generate the index list SPOINT of all atoms of groups j that
C are (possibly) interacting with group i.
      NBSYM=0
      DO J=I,NGRP
      IF (QSTEST(J)) THEN
      DO JJ=IGPBS(J)+1,IGPBS(J+1)
      NBSYM=NBSYM+1
      SPOINT(NBSYM)=JJ
      END DO
      END IF
      END DO
C
C to make the necessary atom-atom search vectorizable we
C gather the "j" interacting atom coordinates and
C interaction array
      DO J=1,NBSYM
      XJF(J)=XF(SPOINT(J))
      YJF(J)=YF(SPOINT(J))
      ZJF(J)=ZF(SPOINT(J))
      INTERJ(J)=INTERE(SPOINT(J))
      END DO
C
C now we have to search the individual atom i to atom j interactions
C
C loop over atoms of group i
      III=0
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
C
C keep track of the number of atoms in this group
C we use this to exclude self-interactions within the
C same group (NOTE: QSTEST(I) is always true)
      III=III+1
C
C If ISYM is 1, the symmetry operation is the identity.  Therefore
C interactions between atoms i, j have to be check whether they are
C subject to the exclusions.
C We update the exclusion list pointers for atom iat.
CCC      IF (ISYM.EQ.1) THEN  ! ATB 2/18/93
      IF (QINTRA.AND.ISYM.EQ.1) THEN
      IF (IAT.GT.1) THEN
      NXI=IBLOEX(IAT-1)+1
      ELSE
      NXI=1
      END IF
      NXIMAX=IBLOEX(IAT)
      END IF
C
      IF (QSYMM) THEN
C
C crystallographic-symmetry mode:
C determine minimum image distance vectors to any "j" atom
      CHCK=ISYM.EQ.1
      CALL XPIMAG(X(IAT),Y(IAT),Z(IAT),III,NBSYM,
     &            XJF,YJF,ZJF,XTEMP,YTEMP,ZTEMP,CHCK,QMIM)
      ELSE
C
C NCS-symmetry mode (includes pure intramolecular mode):
C direct distance vectors to any "j" atom
      DO J=III,NBSYM
      XTEMP(J)=X(IAT)-XJF(J)
      YTEMP(J)=Y(IAT)-YJF(J)
      ZTEMP(J)=Z(IAT)-ZJF(J)
      END DO
C
      END IF
C
C compute distance vector lenghts
      DO J=III,NBSYM
      XTEMP(J)=XTEMP(J)**2+YTEMP(J)**2+ZTEMP(J)**2
      END DO
C
      IF (LGROUP) THEN
C
C group-by-group-cutoff mode:
C   just check interaction criterion
      DO J=III,NBSYM
      QSTEST(J)=ABS(INTERE(IAT)+INTERJ(J)).LE.+1
      END DO
      ELSE
C
C atom-by-atom-cutoff mode:
C    check cutoff and interaction criterion
      DO J=III,NBSYM
      QSTEST(J)=XTEMP(J).LT.CUTNB2
     &           .AND. ABS(INTERE(IAT)+INTERJ(J)).LE.+1
      END DO
      END IF
C
C for the first symmetry operator we have to separate
C intra-molecular and inter-molecular interactions
      IF (ISYM.EQ.1) THEN
C
C crystallographic-symmetry mode:
      IF (QSYMM) THEN
      DO J=III,NBSYM
      QTEST(J)=QSTEST(J).AND.QMIM(J).EQ.0
      QSTEST(J)=QSTEST(J).AND..NOT.QTEST(J)
      END DO
      ELSE
C
C NCS-symmetry mode (includes pure intramolecular mode):
      DO J=III,NBSYM
      QTEST(J)=QSTEST(J)
      QSTEST(J)=.FALSE.
      END DO
      END IF
C
C if required, compute intra-molecular interactions
      IF (QINTRA) THEN
C
C set number of interactions with atom IAT to zero
      NNIAT=0
C
C exclude self-interactions  (III+1!!)
      DO J=III+1,NBSYM
      IF (QTEST(J)) THEN
C
C Check for exclusions
35    IF (NXI.GT.NXIMAX) GOTO 40
      IF (SPOINT(J).EQ.INBEX(NXI)) GOTO 55
      IF (SPOINT(J).LT.INBEX(NXI)) GOTO 40
      NXI=NXI+1
      GOTO 35
C
40    CONTINUE
C found an atom-atom interaction.
      NNNB=NNNB+1
      NNIAT=NNIAT+1
C
C store this interaction in the temporary array
      ILOCAL(NNIAT)=SPOINT(J)
C
C give warning message if appropriate
      IF (XTEMP(J).LT.WMIN2) THEN
      NVIOL=NVIOL+1
      WRITE(PUNIT,'(17A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),'-',RES(IAT),
     &  '-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),'-',
     & RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '" only ',SQRT(XTEMP(J)),' A apart'
      END IF
C     Skip interaction label
55    CONTINUE
      END IF
      END DO
C
C copy number of interactions with atom IAT into JNB
      JNB(IAT)=NNIAT
C copy list of interactions with atom IAT into LIST
      IF (NNIAT.GT.0) THEN
      CALL PUTNB(1,NNIAT,ILOCAL,LIST)
      END IF
      END IF
      END IF
C
C if required, compute symmetry related interactions
      IF (QSYMM.OR.QNCS) THEN
C
C test for special crystallographic positions
      IF (NBQSPC.GT.ZERO) THEN
      IF (ISYM.NE.1) THEN
      IF (ABS(XTEMP(III)).LT.NBQSPC2) THEN
      QSTEST(III)=.FALSE.
      END IF
      END IF
      END IF
C
C set number of interactions with atom IAT to zero
      NNIAT=0
C
C include self-interactions (III !!)
      DO J=III,NBSYM
      IF (QSTEST(J)) THEN
C
C found a symmetry-related atom-atom interaction
      NNSB=NNSB+1
      NNIAT=NNIAT+1
C
C store this interaction in the temporary array
      ILOCAL(NNIAT)=SPOINT(J)
C
C give warning message if appropriate
      IF (XTEMP(J).LT.WMIN2) THEN
      NVIOL=NVIOL+1
C
      IF (QSYMM) WRITE(PUNIT,'(17A,I3,A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),
     &  '-',RES(IAT),'-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),
     & '-',RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '"(XSYM#',ISYM,') only ',SQRT(XTEMP(J)),' A apart'
C
      IF (QNCS) WRITE(PUNIT,'(17A,I3,A,F5.2,A)')
     &  ' %atoms "', SEGID(IAT),'-',RESID(IAT),
     &  '-',RES(IAT),'-',TYPE(IAT),
     & '" and "',SEGID(SPOINT(J)),'-',RESID(SPOINT(J)),
     & '-',RES(SPOINT(J)),'-',TYPE(SPOINT(J)),
     & '"(NCS#',ISYM,') only ',SQRT(XTEMP(J)),' A apart'
      END IF
      END IF
      END DO
C
C copy number of interactions with atom IAT into JSB
      JSB(IAT,ISYM)=NNIAT
C copy list of interactions with atom IAT into SLIST
      IF (NNIAT.GT.0) THEN
      CALL PUTNB(1,NNIAT,ILOCAL,SLIST)
      END IF
      END IF
      END DO
C
      ELSE
C
C this group has no interactions whatsoever, but we still have
C to fill the JNB, JSB arrays for the atoms of this group
      IF (QSYMM.OR.QNCS) THEN
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
      JSB(IAT,ISYM)=0
      END DO
      END IF
      IF (QINTRA.AND.ISYM.EQ.1) THEN
      DO IAT=IGPBS(I)+1,IGPBS(I+1)
      JNB(IAT)=0
      END DO
      END IF
      END IF
C
      END DO
      END DO
C
      IF (QINTRA) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNNB,' intra-atom interactions'
      END IF
      IF (QSYMM) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNSB,' symmetry-atom interactions'
      END IF
      IF (QNCS) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NNSB,' NCS-atom interactions'
      END IF
C
C declare number of violations
      IF (QINTRA.OR.QSYMM.OR.QNCS) THEN
      IF (NVIOL.GT.0) THEN
      WRITE(PUNIT,'(A,I8,A,I8,A)')
     & ' NBONDS: found ',NVIOL,' nonbonded violations'
      END IF
      DBPREC=NVIOL
      CALL DECLAR( 'VIOLATIONS', 'DP', ' ', DBCOMP, DBPREC )
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE XPIMAG(XFROM,YFROM,ZFROM,N,NTO,XTO,YTO,ZTO,DX,DY,DZ,
     &                  CHCK,QMIM)
C
C routine computes minimum image distance vectors between
C reference point i and points j.
C
C
C Input:
C
C   XFROM, YFROM, ZFROM: point i (orthogonal coordinates)
C
C   XTO(*), YTO(*), ZTO(*): array of points j (fractional coordinates)
C   N, NTO: indices specifying the start and the end of the XTO, YTO, ZTO
C   arrays.
C
C
C
C   XRTR, XRINTR: orthogonal-to-fractional and fractional-to-orthogonal
C
C                 matrices.  (passed in common block xcrystal.inc)
C
C
C Output:
C   DX(*), DY(*), DZ(*): is minimum image vector between i and j
C                        (orthogonal)
C   QMIM(*): logical array indicating whether interaction is direct or
C
C            whether it involves a unit-cell translation.  (only when
C            CHCK true -- otherwise no output.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      DOUBLE PRECISION XFROM, YFROM, ZFROM
      INTEGER N, NTO
      DOUBLE PRECISION XTO(*), YTO(*), ZTO(*)
      DOUBLE PRECISION DX(*), DY(*), DZ(*)
      LOGICAL CHCK
      INTEGER QMIM(*)
C local
      DOUBLE PRECISION XFFROM, YFFROM, ZFFROM
      DOUBLE PRECISION XTEMP, YTEMP, ZTEMP
      INTEGER IAT
      LOGICAL QTEST
      DOUBLE PRECISION ZERO, ONE, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, HALF=0.5D0)
C begin
C
C fractionalize XFROM, YFROM, XFROM coordinates (primary cell)
      XFFROM=XRTR(1,1)*XFROM +XRTR(1,2)*YFROM +XRTR(1,3)*ZFROM
      YFFROM=XRTR(2,1)*XFROM +XRTR(2,2)*YFROM +XRTR(2,3)*ZFROM
      ZFFROM=XRTR(3,1)*XFROM +XRTR(3,2)*YFROM +XRTR(3,3)*ZFROM
C
C the reason for duplicating the code is performance on some machines
      IF (CHCK) THEN
C
C loop over all XTO, YTO, ZTO coordinates
      DO IAT=N,NTO
C
C compute difference vector
      XTEMP=XFFROM-XTO(IAT)
      YTEMP=YFFROM-YTO(IAT)
      ZTEMP=ZFFROM-ZTO(IAT)
C
C set the test flag
      QTEST=((ABS(XTEMP).GE.HALF).OR.(ABS(YTEMP).GE.HALF).OR.
     &       (ABS(ZTEMP).GE.HALF))
C
C compute minimum image difference
C
      QMIM(IAT)=0
      IF (QTEST) THEN
C
      XTEMP=INT(ABS(XTEMP)+HALF)*SIGN(ONE,-XTEMP) +XTEMP
      YTEMP=INT(ABS(YTEMP)+HALF)*SIGN(ONE,-YTEMP) +YTEMP
      ZTEMP=INT(ABS(ZTEMP)+HALF)*SIGN(ONE,-ZTEMP) +ZTEMP
C if any of XM, YM, ZM are not equal zero, a unit translation along
C x, y or z is applied to the difference vector.  Set the QMIM flag.
      QMIM(IAT)=1
      END IF
C
C orthogonalize distance vector and store in DX, DY, DZ
      DX(IAT)=XRINTR(1,1)*XTEMP +XRINTR(1,2)*YTEMP +XRINTR(1,3)*ZTEMP
      DY(IAT)=XRINTR(2,1)*XTEMP +XRINTR(2,2)*YTEMP +XRINTR(2,3)*ZTEMP
      DZ(IAT)=XRINTR(3,1)*XTEMP +XRINTR(3,2)*YTEMP +XRINTR(3,3)*ZTEMP
      END DO
      ELSE
C
C loop over all XTO, YTO, ZTO coordinates
      DO IAT=N,NTO
C
C compute difference vector
      XTEMP=XFFROM-XTO(IAT)
      YTEMP=YFFROM-YTO(IAT)
      ZTEMP=ZFFROM-ZTO(IAT)
C
C set the test flag
      QTEST=((ABS(XTEMP).GE.HALF).OR.(ABS(YTEMP).GE.HALF).OR.
     &       (ABS(ZTEMP).GE.HALF))
C
C compute minimum image difference
C
      IF (QTEST) THEN
      XTEMP=INT(ABS(XTEMP)+HALF)*SIGN(ONE,-XTEMP) +XTEMP
      YTEMP=INT(ABS(YTEMP)+HALF)*SIGN(ONE,-YTEMP) +YTEMP
      ZTEMP=INT(ABS(ZTEMP)+HALF)*SIGN(ONE,-ZTEMP) +ZTEMP
      END IF
C
C orthogonalize distance vector and store in DX, DY, DZ
      DX(IAT)=XRINTR(1,1)*XTEMP +XRINTR(1,2)*YTEMP +XRINTR(1,3)*ZTEMP
      DY(IAT)=XRINTR(2,1)*XTEMP +XRINTR(2,2)*YTEMP +XRINTR(2,3)*ZTEMP
      DZ(IAT)=XRINTR(3,1)*XTEMP +XRINTR(3,2)*YTEMP +XRINTR(3,3)*ZTEMP
      END DO
      END IF
C
      RETURN
      END
C=======================================================================
      SUBROUTINE XPGRUP(X,Y,Z,XCENT,YCENT,ZCENT,XSIZ,YSIZ,ZSIZ,
     &                  INTERG,INTERC,INTERE)
C routine computes group centers, size and interaction criterium
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION X(*), Y(*), Z(*), XCENT(*), YCENT(*), ZCENT(*)
      DOUBLE PRECISION XSIZ(*), YSIZ(*), ZSIZ(*)
      INTEGER INTERG(*), INTERC(*), INTERE(*)
C local
      INTEGER I, IS, IL, IAT, NAT
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, ZERO
      DOUBLE PRECISION XC, YC, ZC
      PARAMETER (ZERO=0.0D0)
C begin
      DO I=1,NGRP
      IS=IGPBS(I)+1
      IL=IGPBS(I+1)
      XMIN=X(IS)
      XMAX=XMIN
      YMIN=Y(IS)
      YMAX=YMIN
      ZMIN=Z(IS)
      ZMAX=ZMIN
      XC=ZERO
      YC=ZERO
      ZC=ZERO
      NAT=0
      INTERG(I)=0
      INTERC(I)=0
      DO IAT=IS,IL
      IF (INITIA(IAT,X,Y,Z)) THEN
      XC=XC+X(IAT)
      YC=YC+Y(IAT)
      ZC=ZC+Z(IAT)
      NAT=NAT+1
      XMIN=MIN(XMIN,X(IAT))
      YMIN=MIN(YMIN,Y(IAT))
      ZMIN=MIN(ZMIN,Z(IAT))
      XMAX=MAX(XMAX,X(IAT))
      YMAX=MAX(YMAX,Y(IAT))
      ZMAX=MAX(ZMAX,Z(IAT))
C
C only count interacting atoms in group
      IF (INTERE(IAT).LT.+9999) THEN
      INTERG(I)=INTERG(I)+INTERE(IAT)
      INTERC(I)=INTERC(I)+1
      END IF
      END IF
      END DO
C
C center of geometry of group
      IF (NAT.GT.0) THEN
      XC=XC/NAT
      YC=YC/NAT
      ZC=ZC/NAT
      END IF
      XCENT(I)=XC
      YCENT(I)=YC
      ZCENT(I)=ZC
C
C size of rectangular box surrounding group
      XSIZ(I)=MAX(XMAX-XC,XC-XMIN)
      YSIZ(I)=MAX(YMAX-YC,YC-YMIN)
      ZSIZ(I)=MAX(ZMAX-ZC,ZC-ZMIN)
C
C set INTERG to +9999 if it has no interacting atoms at all
      IF (INTERC(I).EQ.0) INTERG(I)=+9999
      END DO
C
      RETURN
      END
C
      SUBROUTINE XPSYMO(NATOM,X,Y,Z,ISYM,XF,YF,ZF)
C
C routine computes symmetry related molecule, symmetry op. number ISYM
C input is X,Y,Z in orthogonal coordinates, output is XF,YF,ZF
C in fractional (!!) coordinates.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INTEGER NATOM
      DOUBLE PRECISION X(NATOM), Y(NATOM), Z(NATOM)
      INTEGER ISYM
      DOUBLE PRECISION XF(NATOM), YF(NATOM), ZF(NATOM)
C local
      DOUBLE PRECISION  SYFR(3,4), RTH, ZERO, ONE, TEN
      INTEGER I, J, K
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TEN=10.0D0)
C begin
C
C concatenate symmetry operator and orthogonal-to-fractional operator
      DO I=1,3
      DO J=1,3
      SYFR(I,J)=ZERO
      DO K=1,3
      SYFR(I,J)=SYFR(I,J)+XRSYMM(ISYM,I,K)*XRTR(K,J)
      END DO
      END DO
      END DO
      RTH=ONE/XRSYTH
      DO I=1,3
      SYFR(I,4)=XRSYMM(ISYM,I,4)*RTH
      END DO
C
C apply SYFR operator to coords and make sure they are in primary cell
      DO I=1,NATOM
      XF(I)=SYFR(1,1)*X(I) +SYFR(1,2)*Y(I)
     &         +SYFR(1,3)*Z(I) +SYFR(1,4)
      YF(I)=SYFR(2,1)*X(I) +SYFR(2,2)*Y(I)
     &         +SYFR(2,3)*Z(I) +SYFR(2,4)
      ZF(I)=SYFR(3,1)*X(I) +SYFR(3,2)*Y(I)
     &         +SYFR(3,3)*Z(I) +SYFR(3,4)
      END DO
      RETURN
      END
C==================================================================
      SUBROUTINE UPINB(N)
C
C Do an INB update as specified. Front end for MKINB.
C
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      INTEGER  N
C local
      INTEGER NATBON, IATMXB, I, IATBON
      INTEGER IPK, JPK, IPK14, JPK14
      INTEGER NPAIR, NPAI14
      LOGICAL CMPLTD
C begin
C
      IATMXB=12
C
C first, generate an exclusion pair list and a 1-4 interaction
C pair list. Make estimates for size of pair lists until
C completed.
      I=3*NBOND
      IF(IABS(NBXMOD).GT.3) I=I*2
      IF(NBXMOD.GT.0) I=MAX(NNB+1,I)
      CMPLTD = .FALSE.
      DO WHILE (.NOT. CMPLTD)
      IPK=ALLHP(INTEG4(I))
      JPK=ALLHP(INTEG4(I))
      IPK14=ALLHP(INTEG4(I/2))
      JPK14=ALLHP(INTEG4(I/2))
C
C allocate space for the temporary bond lists
      NATBON=ALLHP(INTEG4(NATOM))
      IATBON=ALLHP(INTEG4(IATMXB*NATOM))
      CALL MAKINB(NATOM,HEAP(IINTER(N)),IB,JB,NBOND,INB,IBLO,NNB,
     &     NBXMOD,HEAP(NATBON),HEAP(IATBON),IATMXB,CMPLTD,
     &     I,NPAIR,HEAP(IPK),HEAP(JPK),
     &     I/2,NPAI14,HEAP(IPK14),HEAP(JPK14))
      CALL FREHP(IATBON,INTEG4(IATMXB*NATOM))
      CALL FREHP(NATBON,INTEG4(NATOM))
      IF (.NOT.(CMPLTD) )THEN
      WRITE(6,'(A)') ' MAKINB: Ran out of space. RESIZING'
      CALL FREHP(JPK14,INTEG4(I/2))
      CALL FREHP(IPK14,INTEG4(I/2))
      CALL FREHP(JPK,INTEG4(I))
      CALL FREHP(IPK,INTEG4(I))
      I=I*1.5D0+10
      END IF
      END DO
C
C now generate the lists IBLOEX,INBEX and LIST14,INB14
C first, resize the non-bonding lists
      NNNB(N)=0
      NNSB(N)=0
      NNB14(N)=NPAI14
      NNBEX(N)=NPAIR
      CALL RENBND(0,0,0,0,0,NATOM,NNB14(N),NATOM,NNBEX(N),N)
      CALL RESNB(LIST14(1,N))
      CALL GENINB(NBXMOD,NATOM,
     &            NNB14(N),LIST14(1,N),HEAP(IINB14(N)),
     &            NNBEX(N),HEAP(IIBLOE(N)),HEAP(IINBEX(N)),
     &            NPAIR,HEAP(IPK),HEAP(JPK),
     &            NPAI14,HEAP(IPK14),HEAP(JPK14))
C
      CALL FREHP(JPK14,INTEG4(I/2))
      CALL FREHP(IPK14,INTEG4(I/2))
      CALL FREHP(JPK,INTEG4(I))
      CALL FREHP(IPK,INTEG4(I))
C
      RETURN
      END
C================================================================
      SUBROUTINE MAKINB(NATOM,INTERE,IB,JB,NBOND,INB,IBLO,NNB,
     &     MODE,NATBON,IATBON,IATBMX,CMPLTD,
     &     MAXWRK,NPAIR,IPK,JPK,MXWRK4,NPAIR4,IPK14,JPK14)
C
C Routine generates a exclusion pair list as well as the special
C 1-4 interaction list (if requested). Generation is just
C based on molecular topology bond list (connectivity).
C
C definition of MODE:
C     MODE = +- 1  exclude nothing
C     MODE = +- 2  exclude only 1-2 (bond) interactions
C     MODE = +- 3  exclude 1-2 and 1-3 (bond and angle)
C     MODE = +- 4  include 1-2 1-3 and 1-4'S
C     MODE = +- 5  include 1-2 1-3 and 1-4'S put on as 1-4 interactions
C  A positive mode value causes the explicit exclusion array INB
C  to be added.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'funct.inc'
      INTEGER NATOM, INTERE(*), IB(*), JB(*), NBOND, INB(*), IBLO(*)
      INTEGER NNB, MODE, NATBON(*), IATBMX, IATBON(IATBMX,*)
      LOGICAL CMPLTD
      INTEGER MAXWRK,NPAIR,IPK(*),JPK(*)
      INTEGER MXWRK4,NPAIR4,IPK14(*),JPK14(*)
C local
      EXTERNAL EXCH5
      INTEGER MODEX, ILAST, I, J, I2, J2, IBT, JBT, IJ, IJA, IK, K, IKA
      INTEGER II, IOLD, JOLD, I14, II14
C begin
C
      CMPLTD=.FALSE.
      MODEX=ABS(MODE)
C
      NPAIR=0
      NPAIR4=0
C
C For mode greater than zero include the explicit INB/IBLO lists,
      IF (MODE.GT.0) THEN
      ILAST=0
      DO I=1,NATOM
      IF (IBLO(I).GT.ILAST) THEN
      DO J=ILAST+1,IBLO(I)
      I2=MIN(I,INB(J))
      J2=MAX(I,INB(J))
      IF (I2.GT.0.AND.I2.NE.J2) THEN
      NPAIR=NPAIR+1
      IF(NPAIR.GT.MAXWRK) RETURN
      IPK(NPAIR)=I2
      JPK(NPAIR)=J2
      END IF
      END DO
      ILAST=IBLO(I)
      END IF
      END DO
      END IF
C
C include 1-2 exclusions
      IF (MODEX.GT.1) THEN
      DO I=1,NBOND
      I2=MIN(IB(I),JB(I))
      J2=MAX(IB(I),JB(I))
      IF (I2.GT.0.AND.I2.NE.J2) THEN
      NPAIR=NPAIR+1
      IF(NPAIR.GT.MAXWRK) RETURN
      IPK(NPAIR)=I2
      JPK(NPAIR)=J2
      END IF
      END DO
      END IF
C
C Next make a list of all the possible 1-3 and 1-4 interactions.
C This is based on the bond list.  First make a list of all bonds
C for every atom.
      IF (MODEX.GT.2) THEN
      DO I=1,NATOM
      NATBON(I)=0
      END DO
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      IF (IBT.GT.0 .AND. JBT.GT.0) THEN
      NATBON(IBT)=NATBON(IBT)+1
      IF (NATBON(IBT).GT.IATBMX) THEN
      WRITE(6,335) IBT
335   FORMAT(' %MAKINB-ERR: Too many bonds for atom',I5,' Check code')
      CALL DIE
      END IF
      IATBON(NATBON(IBT),IBT)=I
      NATBON(JBT)=NATBON(JBT)+1
      IF (NATBON(JBT).GT.IATBMX) THEN
      WRITE(6,335) JBT
      CALL DIE
      END IF
      IATBON(NATBON(JBT),JBT)=-I
      END IF
      END DO
C
C Now make the unsorted list of 1-3 interactions by taking bonds
C and extending in a direction one bond.
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      DO J=1,NATBON(IBT)
      IJ=IATBON(J,IBT)
      IF (ABS(IJ) .GT. I) THEN
      IF (IJ.GT.0) THEN
      IJA=JB(IJ)
      ELSE
      IJA=IB(ABS(IJ))
      END IF
      I2=MIN(IJA,JBT)
      J2=MAX(IJA,JBT)
      IF (I2.GT.0) THEN
      NPAIR=NPAIR+1
      IF(NPAIR.GT.MAXWRK) RETURN
      IPK(NPAIR)=I2
      JPK(NPAIR)=J2
      END IF
      END IF
      END DO
      DO J=1,NATBON(JBT)
      IJ=IATBON(J,JBT)
      IF (ABS(IJ) .GT. I) THEN
      IF (IJ.GT.0) THEN
      IJA=JB(IJ)
      ELSE
      IJA=IB(ABS(IJ))
      END IF
      I2=MIN(IJA,IBT)
      J2=MAX(IJA,IBT)
      IF (I2.GT.0.AND.I2.NE.J2) THEN
      NPAIR=NPAIR+1
      IF(NPAIR.GT.MAXWRK) RETURN
      IPK(NPAIR)=I2
      JPK(NPAIR)=J2
      END IF
      END IF
      END DO
      END DO
      END IF
C
C Now make the unsorted list of 1-4 interactions by taking bonds
C and extending in each direction one bond.
      IF (MODEX.GT.3) THEN
      DO I=1,NBOND
      IBT=IB(I)
      JBT=JB(I)
      DO J=1,NATBON(IBT)
      IJ=IATBON(J,IBT)
      IF (ABS(IJ) .NE. I) THEN
      IF (IJ.GT.0) THEN
      IJA=JB(IJ)
      ELSE
      IJA=IB(ABS(IJ))
      END IF
      DO K=1,NATBON(JBT)
      IK=IATBON(K,JBT)
      IF (ABS(IK) .NE. I) THEN
      IF (IK.GT.0) THEN
      IKA=JB(IK)
      ELSE
      IKA=IB(ABS(IK))
      END IF
      I2=MIN(IJA,IKA)
      J2=MAX(IJA,IKA)
      IF (I2.GT.0.AND.I2.NE.J2) THEN
      NPAIR=NPAIR+1
      IF(NPAIR.GT.MAXWRK) RETURN
      IPK(NPAIR)=I2
      JPK(NPAIR)=J2
C
C ... store 1-4 interaction in special pair list if requested
      IF (MODEX.EQ.5) THEN
C
C check for interaction atom pairs
      IF (ABS(INTERE(I2)+INTERE(J2)).LE.+1) THEN
      NPAIR4=NPAIR4+1
      IF(NPAIR4.GT.MXWRK4) RETURN
      IPK14(NPAIR4)=I2
      JPK14(NPAIR4)=J2
      END IF
      END IF
      END IF
      END IF
      END DO
      END IF
      END DO
      END DO
      END IF
C
C include a mark at the end of 1-4 list
      NPAIR4=NPAIR4+1
      IF(NPAIR4.GT.MXWRK4) RETURN
      IPK14(NPAIR4)=NATOM
      JPK14(NPAIR4)=NATOM
C
C Next sort the list of all the possible 1-4 interactions.
      CALL SORT(NPAIR4,EXCH5,ORDER5,IPK14,JPK14,0,0,0,0,0,2)
C
C remove all fixed atom exclusions from pair list
      II=0
      DO I=1,NPAIR
      I2=IPK(I)
      J2=JPK(I)
      IF (ABS(INTERE(I2)+INTERE(J2)).LE.+1) THEN
      II=II+1
      IPK(II)=I2
      JPK(II)=J2
      END IF
      END DO
      NPAIR=II
C
C Sort the pair list.
      CALL SORT(NPAIR,EXCH5,ORDER5,IPK,JPK,0,0,0,0,0,2)
C
C remove duplications in IPK,JPK as well as in IPK14, JPK14.
C If there is a duplication in IPK, JPK remove the corresponding
C 1-4 interaction if necessary. We assume that the 1-4 list
C is a subset of the IPK, JPK list.
      II=0
      IOLD=0
      JOLD=0
      I14=1
      II14=0
      DO I=1,NPAIR
C
C only do this if IPK14(I14), JPK14(I14) are well-defined, ATB 2/22/10
      IF (I14.GT.0.AND.I14.LE.NPAIR4) THEN
      IF (IPK(I).EQ.IPK14(I14).AND.JPK(I).EQ.JPK14(I14)) THEN
      II14=II14+1
      IPK14(II14)=IPK14(I14)
      JPK14(II14)=JPK14(I14)
      I14=I14+1
      END IF
      END IF
C
      IF (IPK(I).NE.IOLD.OR.JPK(I).NE.JOLD) THEN
      II=II+1
      IPK(II)=IPK(I)
      JPK(II)=JPK(I)
      IOLD=IPK(I)
      JOLD=JPK(I)
      ELSE
C we have a duplication
C
C the following lines were taken out.  Otherwise, 1-4
C interactions in 6-membered rings were incorrectly removed.
CCCCCC      DO WHILE (IPK14(II14).EQ.IOLD.AND.JPK14(II14).EQ.JOLD)
CCCCCC      II14=II14-1
CCCCCC      END DO
CCC modification ATB 4/17/92
C
C only do this if IPK14(II14), JPK14(II14) are well-defined, ATB 2/22/10
      IF (II14.GT.0.AND.II14.LE.NPAIR4) THEN
      IF (IPK14(II14).EQ.IOLD.AND.JPK14(II14).EQ.JOLD) THEN
      II14=II14-1
      END IF
      END IF
C
      END IF
      END DO
C
      NPAIR=II
      NPAIR4=II14
C
      WRITE(6,432) MODE,NPAIR-NPAIR4,NPAIR4
 432  FORMAT(' MAKINB: mode',I4,' found',I7,' exclusions and',I7,
     1     ' interactions(1-4)')
C
      CMPLTD=.TRUE.
      RETURN
      END
C====================================================================
      SUBROUTINE GENINB(MODE,NATOM,NNB14,LIST14,INB14,
     &            NNBEX,IBLOEX,INBEX,
     &            NPAIR,IPK,JPK,NPAI14,IPK14,JPK14)
C
C Converts pair list IPK,JPK to IBLOEX,INBEX and
C IPK14,JPK14 to LIST14,INB14.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER MODE, NATOM, NNB14, LIST14(4), INB14(*), NNBEX
      INTEGER IBLOEX(*), INBEX(*), NPAIR, IPK(*), JPK(*), NPAI14
      INTEGER IPK14(*), JPK14(*)
C local
      INTEGER IATOM, I, II, J
C begin
      NNBEX=NPAIR
      NNB14=NPAI14
C
      IATOM=1
      DO I=1,NPAIR
      IF (IPK(I).GT.IATOM) THEN
      DO J=IATOM,IPK(I)-1
      IBLOEX(J)=I-1
      END DO
      IATOM=IPK(I)
      END IF
      INBEX(I)=JPK(I)
      END DO
      DO J=IATOM,NATOM
      IBLOEX(J)=NPAIR
      END DO
C
      DO I=1,NATOM
      INB14(I)=0
      END DO
C
      IF (NPAI14.GT.0) THEN
      II=1
      IATOM=IPK14(1)
      DO I=2,NPAI14
      IF (IATOM.EQ.IPK14(I)) THEN
      II=II+1
      ELSE
      INB14(IATOM)=II
      CALL PUTNB(I-II,I-1,JPK14,LIST14)
      II=1
      IATOM=IPK14(I)
      END IF
      END DO
      INB14(IATOM)=II
      CALL PUTNB(NPAI14-II+1,NPAI14,JPK14,LIST14)
      END IF
C
      RETURN
      END
