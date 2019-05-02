      SUBROUTINE XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
C
C Routine reduces all reflections.
C   0. checks self-consistency of symmetry operators
C   1. maps reflections into asymmetric unit as defined by
C      routine XRASYM
C   2. sorts the data
C   3. discards symmetry-related reflections,
C      Friedel mates if hermitian symmetry is turned on.
C   4. removes systematic absences
C   5. computes multiplicity and type for each reflection
C
C definition of type:
C    positive for acentric reflections
C    negative for centric reflections
C    for reflections with a Friedel mate: the absolute value
C    of type contains the reflection index and the symmetry
C    operator index.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      LOGICAL XRRED
      INTEGER XRMREF, XRNREF
      INTEGER HPTSEL, HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
C local
      INTEGER LDIM
      LOGICAL ERR
C pointer
      INTEGER MAP, MAP2, MARK, HHH, KKK, LLL
      INTEGER GRPPTR, GRPOBJ
C begin
C
      IF (XRRED) THEN
C
C reset the XRRED flag
      XRRED=.FALSE.
C
C
      ERR=.FALSE.
C
C check self-consistency of symmetry operators
      CALL XSYMCHK(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,ERR,
     &             XRSYGP,XRSYIV)
C check consistency of symmetry operators and cell dimensions
      CALL XSYMCEL(XRINTR,XRNSYM,XRMSYM,XRSYMM)
C
      IF (XRNREF.GT.0.AND..NOT.ERR) THEN
C
C
      LDIM=XRNREF
C
      MAP=ALLHP(INTEG4(LDIM))
      MAP2=ALLHP(INTEG4(LDIM))
      MARK=ALLHP(INTEG4(LDIM))
      HHH=ALLHP(INTEG4(LDIM))
      KKK=ALLHP(INTEG4(LDIM))
      LLL=ALLHP(INTEG4(LDIM))
      GRPPTR=ALLHP(INTEG4(XSFNUM+1))
      GRPOBJ=ALLHP(INTEG4(XSFNUM))
      CALL XREDU2(XRMREF,XRNREF,HEAP(HPTSEL),
     &           HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HEAP(HPMULT),HEAP(HPTYPE),
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           HEAP(MAP),HEAP(MAP2),
     &           HEAP(MARK),HEAP(HHH),HEAP(KKK),HEAP(LLL),XRSYGP,XRSYIV,
     &           HEAP(GRPPTR), HEAP(GRPOBJ))
      CALL FREHP(GRPOBJ,INTEG4(XSFNUM))
      CALL FREHP(GRPPTR,INTEG4(XSFNUM+1))
      CALL FREHP(KKK,INTEG4(LDIM))
      CALL FREHP(LLL,INTEG4(LDIM))
      CALL FREHP(HHH,INTEG4(LDIM))
      CALL FREHP(MARK,INTEG4(LDIM))
      CALL FREHP(MAP2,INTEG4(LDIM))
      CALL FREHP(MAP,INTEG4(LDIM))
      END IF
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE XREDU2(XRMREF,XRNREF,TSEL,XRH,XRK,XRL,
     &           MULT,TYPE,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           MAP,MAP2,MARK,HHH,KKK,LLL,XRSYGP,XRSYIV,
     &           GRPPTR,GRPOBJ)
C
C see routine XREDUC above
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, XRNREF, TSEL(*), XRH(*), XRK(*), XRL(*)
      INTEGER MULT(*), TYPE(*)
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER MAP(*), MAP2(*), MARK(*), HHH(*), KKK(*), LLL(*)
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      INTEGER GRPPTR(*), GRPOBJ(*)
C local
      INTEGER REFLCT, HH, KK, LL, N, IFRIED, IISYM, MUL, I
      INTEGER DELREF, IREF, IREF2, ITEMP
      DOUBLE PRECISION RTH, PHAS, TEMP
      DOUBLE COMPLEX SHIFT, CTEMP
      LOGICAL QMARK, SYSAB, QFIRS
      INTEGER HMAX, KMAX, LMAX, HMIN, KMIN, LMIN, MDIM
      INTEGER PA, PB, PC, PD, NGROUP, IGROUP, COUNTER
      LOGICAL FOUND
      DOUBLE PRECISION PHASINV, PHASSHIFT
C pointer
      INTEGER IPTR
      INTEGER MATRIX
C parameter
      DOUBLE PRECISION ONE, TWO, ZERO, FOUR, S180
      PARAMETER (ONE=1.0D0, TWO=2.0D0, ZERO=0.0D0, FOUR=4.0D0)
      PARAMETER (S180=180.D0)
C begin
C
      RTH=XRSYTH
C
      QFIRS=.TRUE.
C
C get object group list
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
C
C initialize the removal flags for all reflections and copy
C indices into HHH, KKK, LLL set
      DO REFLCT=1,XRNREF
      MARK(REFLCT)=0
      HHH(REFLCT)=XRH(REFLCT)
      KKK(REFLCT)=XRK(REFLCT)
      LLL(REFLCT)=XRL(REFLCT)
      END DO
C
C map all reflections into the asymmetric unit
      DO REFLCT=1,XRNREF
      CALL XRASYM(XRH(REFLCT),XRK(REFLCT),XRL(REFLCT),HH,KK,LL,
     &    IISYM,IFRIED,XRNSYM,XRMSYM,XRITSY,QHERM)
      XRH(REFLCT)=HH
      XRK(REFLCT)=KK
      XRL(REFLCT)=LL
C
C apply phase shift if the indeces are different (symmetry operator
C and/or Friedel operator was required to map the reflection into
C the asymmetric unit)
C
      IF (IISYM.NE.1.OR.IFRIED.EQ.-1) THEN
C
      IF (QFIRS) THEN
      QFIRS=.FALSE.
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(2A)')
     & ' XREDUC: some reflection(s) converted',
     & ' to CNS standard asymm. unit.'
      END IF
      END IF
C
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*HH
     &              +XRSYMM(IISYM,2,4)*KK
     &              +XRSYMM(IISYM,3,4)*LL ) / RTH )
      PHASSHIFT=PHAS
C
      IF (IFRIED.EQ.-1) THEN
      PHAS=-PHAS
      PHASINV=-ONE
      ELSE
      PHASINV=ONE
      END IF
C
      SHIFT=DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
C
C map reciprocal space objects except grouped arrays
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFGNAM(I).EQ.0) THEN
      IF (XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(I)),REFLCT)
      CTEMP=CTEMP*SHIFT
      IF (IFRIED.EQ.-1) CTEMP=DCONJG(CTEMP)
      CALL XCOPY(HEAP(HPSF(I)),REFLCT,CTEMP,1)
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
C do nothing
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
C do nothing
      END IF
      END IF
      END DO
C
C loop over all groups
      DO IGROUP=1,NGROUP
C
C check if all elements are well-defined
      FOUND=.TRUE.
      DO COUNTER=GRPPTR(IGROUP)+1,GRPPTR(IGROUP+1)
      FOUND=FOUND.AND.HPSF(GRPOBJ(COUNTER)).NE.0
      END DO
C
      IF (FOUND) THEN
C
C map this set of HL coefficients
      PA=GRPOBJ(GRPPTR(IGROUP)+1)
      PB=GRPOBJ(GRPPTR(IGROUP)+2)
      PC=GRPOBJ(GRPPTR(IGROUP)+3)
      PD=GRPOBJ(GRPPTR(IGROUP)+4)
      CALL XMAPHL(REFLCT,HEAP(HPSF(PA)),HEAP(HPSF(PB)),
     &           HEAP(HPSF(PC)),HEAP(HPSF(PD)),PHAS,
     &           ONE,PHASINV)
      END IF
C
      END DO
C
      END IF
      END DO
C
C
C find the permutation "map" that sorts all reflections.  Note
C that SORTP conserves the order of reflections if all indices
C are equal.  Thus, among all symmetry-related and redundant
C reflections, the one that was read the first time is always
C on top of the list.
      CALL SORTP(XRNREF,MAP,ORDER5,XRL,XRK,XRH,LLL,KKK,HHH,0,6)
C
C check for multiple entries, removing systematic absences,
C and re-define the "map"
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(2A)')
     &  ' XREDU2: checking for multiple entries. ',
     &  'Removing systematic absences.'
      END IF
C
      QMARK=.FALSE.
      DELREF=0
      N=0
      MUL=0
      DO REFLCT=1,XRNREF
C
C check if this is a systematic absence
      CALL XRSYSAB(XRH(MAP(REFLCT)),XRK(MAP(REFLCT)),
     &             XRL(MAP(REFLCT)),SYSAB,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
      IF (REFLCT.EQ.1) THEN
      IF (.NOT.SYSAB) THEN
C
C first entry -- no systematic absence
C ------------------------------------
      N=N+1
      MAP2(1)=MAP(1)
      MARK(MAP(1))=0
      ELSE
C first entry -- systematic absence
C ---------------------------------
      MARK(MAP(REFLCT))=0
      QMARK=.TRUE.
      DELREF=DELREF+1
      END IF
C
      ELSE
      IF (SYSAB) THEN
C
C systematic absence
C ------------------
      MARK(MAP(REFLCT))=0
      QMARK=.TRUE.
      DELREF=DELREF+1
C
      ELSEIF (.NOT.SYSAB.AND.
     &        (XRH(MAP(REFLCT-1)).EQ.XRH(MAP(REFLCT))
     &    .AND.XRK(MAP(REFLCT-1)).EQ.XRK(MAP(REFLCT))
     &    .AND.XRL(MAP(REFLCT-1)).EQ.XRL(MAP(REFLCT)))) THEN
C
C this is a multiple entry
C ------------------------
      QMARK=.TRUE.
      DELREF=DELREF+1
      IF (MARK(MAP(REFLCT-1)).NE.0) THEN
      MARK(MAP(REFLCT))=MARK(MAP(REFLCT-1))
      ELSE
      MUL=MUL+1
      MARK(MAP(REFLCT))=MUL
      MARK(MAP(REFLCT-1))=MUL
      END IF
      ELSE
C
C no systematic absence, no multiple entry
C ----------------------------------------
      N=N+1
      MAP2(N)=MAP(REFLCT)
      MARK(MAP(REFLCT))=0
      END IF
      END IF
C
      END DO
C
C
      IF (QMARK) THEN
      WRITE(6,'(A)')
     & ' XREDU2: multiple instances found for some reflections.'
      WRITE(6,'(A)')
     & ' XREDU2: deleting redundant reflections.'
C
C multiple reflections found.  Need to remove the
C redundancies.  We also sort the list of reflections.
C
C print reflections that will be discarded
      IF (WRNLEV.GE.10) THEN
      MUL=0
      DO REFLCT=1,XRNREF
      IF (MARK(MAP(REFLCT)).NE.0) THEN
      IF (MUL.NE.MARK(MAP(REFLCT))) THEN
      MUL=MARK(MAP(REFLCT))
      WRITE(6,'(A)') ' ========================================'
      WRITE(6,'(A,3I5)')
     &  ' XREDU2: multiple occurrence of : ',
     &    HHH(MAP(REFLCT)),KKK(MAP(REFLCT)),LLL(MAP(REFLCT))
      ELSE
      WRITE(6,'(A,3I5)')
     &  ' XREDU2: multiple entry being merged: ',
     &    HHH(MAP(REFLCT)),KKK(MAP(REFLCT)),LLL(MAP(REFLCT))
      END IF
      END IF
      END DO
      END IF
C
C merge multiple entries
      MUL=0
      DO REFLCT=1,XRNREF
      IF (MARK(MAP(REFLCT)).NE.0) THEN
      IF (MUL.NE.MARK(MAP(REFLCT))) THEN
      MUL=MARK(MAP(REFLCT))
      IREF=MAP(REFLCT)
      ELSE
C
      IREF2=MAP(REFLCT)
C
C merge everything into reflection marked IREF
C reciprocal space objects
C
C Any non-zero entries will be copied into the reflection
C that was first read among the set of symmetry-related
C and redundant reflections.  This is the reflection
C that will be kept and all other reflections related by
C symmetry or redundancy will be removed after this operation.
C
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0) THEN
C
      IF (XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(I)),IREF)
      IF (CTEMP.EQ.DCMPLX(ZERO,ZERO)) THEN
      CALL XCOPY(HEAP(HPSF(I)),IREF,HEAP(HPSF(I)),IREF2)
      END IF
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      CALL XCOPYR(TEMP,1,HEAP(HPSF(I)),IREF)
      IF (TEMP.EQ.ZERO) THEN
      CALL XCOPYR(HEAP(HPSF(I)),IREF,HEAP(HPSF(I)),IREF2)
      END IF
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      CALL XCOPYI(ITEMP,1,HEAP(HPSF(I)),IREF)
      IF (ITEMP.EQ.0) THEN
      CALL XCOPYI(HEAP(HPSF(I)),IREF,HEAP(HPSF(I)),IREF2)
      END IF
      END IF
      END IF
      END DO
C
C
      END IF
      END IF
      END DO
C
C
C now we throw out multiple entries by mapping all data
      IPTR=ALLHP(INTEG4(N))
      CALL AINDX4(MAP2,XRH,N,HEAP(IPTR))
      CALL AINDX4(MAP2,XRL,N,HEAP(IPTR))
      CALL AINDX4(MAP2,XRK,N,HEAP(IPTR))
      CALL AINDX4(MAP2,MULT,N,HEAP(IPTR))
      CALL AINDX4(MAP2,TYPE,N,HEAP(IPTR))
      CALL AINDX4(MAP2,TSEL,N,HEAP(IPTR))
      CALL FREHP(IPTR,INTEG4(N))
C
C reciprocal space objects
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0) THEN
C
      IF (XSFTYPE(I).EQ.'COMP') THEN
      IPTR=ALLHP(ICPLX8(N))
      CALL AINDC8(MAP2,HEAP(HPSF(I)),N,HEAP(IPTR))
      CALL FREHP(IPTR,ICPLX8(N))
C
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      IPTR=ALLHP(IREAL8(N))
      CALL AINDR8(MAP2,HEAP(HPSF(I)),N,HEAP(IPTR))
      CALL FREHP(IPTR,IREAL8(N))
C
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      IPTR=ALLHP(INTEG4(N))
      CALL AINDX4(MAP2,HEAP(HPSF(I)),N,HEAP(IPTR))
      CALL FREHP(IPTR,INTEG4(N))
      END IF
C
      END IF
      END DO
C
C
C update number of reflections
      XRNREF=N
      WRITE(6,'(A,I8,A)') ' XREDU2: ',DELREF,
     &  ' reflections have been deleted.'
      WRITE(6,'(A)') ' XREDU2: sorting all reflections.'
      END IF
C
C determine multiplicity and type for each reflection
C +++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C allocate space for a bookkeeping matrix
      CALL XRRR3(XRNREF,XRH,XRK,XRL,HMAX,KMAX,LMAX,
     &                 HMIN,KMIN,LMIN)
C
      IF (HMAX.EQ.0) HMAX=1
      IF (KMAX.EQ.0) KMAX=1
      IF (LMAX.EQ.0) LMAX=1
      MDIM=(HMAX-HMIN+1)*(KMAX-KMIN+1)*(LMAX-LMIN+1)
      MATRIX=ALLHP(INTEG4(MDIM))
      CALL XREDU4(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                  QHERM,XRNREF,XRH,XRK,XRL,TYPE,MULT,
     &                  HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,HEAP(MATRIX))
      CALL FREHP(MATRIX,INTEG4(MDIM))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XREDU4(XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &                  QHERM,XRNREF,XRH,XRK,XRL,
     &                  TYPE,MULT,HMAX,KMAX,LMAX,
     &                  HMIN,KMIN,LMIN,MATRIX)
C
C Routine computes multiplicity, type (centric or acentric),
C pointer and symmetry operator to produce the Friedel mate if
C it is in the data base.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*), MULT(*)
      INTEGER HMAX, KMAX, LMAX, HMIN, KMIN, LMIN
      INTEGER MATRIX(HMIN:HMAX,KMIN:KMAX,LMIN:LMAX)
C local
      INTEGER K, L, H, REFLCT, MUL, MUL2, ISYM, HH, KK, LL
      DOUBLE PRECISION RTH
      INTEGER ICEN
C parameter
C
C begin
C
C determine number of centering operations for this spacegroup
      ICEN=0
      DO ISYM=1,XRNSYM
      IF (XRSYMM(1,1,1).EQ.XRSYMM(ISYM,1,1).AND.
     &    XRSYMM(1,1,2).EQ.XRSYMM(ISYM,1,2).AND.
     &    XRSYMM(1,1,3).EQ.XRSYMM(ISYM,1,3).AND.
     &    XRSYMM(1,2,1).EQ.XRSYMM(ISYM,2,1).AND.
     &    XRSYMM(1,2,2).EQ.XRSYMM(ISYM,2,2).AND.
     &    XRSYMM(1,2,3).EQ.XRSYMM(ISYM,2,3).AND.
     &    XRSYMM(1,3,1).EQ.XRSYMM(ISYM,3,1).AND.
     &    XRSYMM(1,3,2).EQ.XRSYMM(ISYM,3,2).AND.
     &    XRSYMM(1,3,3).EQ.XRSYMM(ISYM,3,3)) THEN
      ICEN=ICEN+1
      END IF
      END DO
C
C
      RTH=XRSYTH
C
C setup bookkeeping matrix
      DO L=LMIN,LMAX
      DO K=KMIN,KMAX
      DO H=HMIN,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
      DO REFLCT=1,XRNREF
      MATRIX(XRH(REFLCT),XRK(REFLCT),XRL(REFLCT))=REFLCT
      END DO
C
      DO REFLCT=1,XRNREF
      MUL=0
C
C initialize the reflection type to acentric, no Friedel mate in
C date base.
      TYPE(REFLCT)=1
C go through all symmetry operators
      DO ISYM=1,XRNSYM
C generate symmetry related reflections
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      HH=XRITSY(ISYM,1,1)*H + XRITSY(ISYM,1,2)*K + XRITSY(ISYM,1,3)*L
      KK=XRITSY(ISYM,2,1)*H + XRITSY(ISYM,2,2)*K + XRITSY(ISYM,2,3)*L
      LL=XRITSY(ISYM,3,1)*H + XRITSY(ISYM,3,2)*K + XRITSY(ISYM,3,3)*L
C
C symmetry operator produces same reflection:
      IF (H.EQ.HH.AND.K.EQ.KK.AND.L.EQ.LL) MUL=MUL+1
C
C symmetry operator produces Friedel mate:
      IF (H.EQ.-HH.AND.K.EQ.-KK.AND.L.EQ.-LL) THEN
C
C we've got a centric reflection !
      TYPE(REFLCT)=-ABS(TYPE(REFLCT))
      IF (QHERM) MUL=MUL+1
      END IF
C
      IF (HMIN.LE.-HH.AND.-HH.LE.HMAX.AND.
     &    KMIN.LE.-KK.AND.-KK.LE.KMAX.AND.
     &    LMIN.LE.-LL.AND.-LL.LE.LMAX.AND.ABS(TYPE(REFLCT)).EQ.1) THEN
      IF (MATRIX(-HH,-KK,-LL).NE.0) THEN
C
C found the Friedel mate for this reflection
C
C encode the reflection index and the symmetry operator
C id into a single integer.  The sign of this
C integer determines if the reflection is acentric (+)
C or centric (-).  If the integer is +-1 this means that
C there is no Friedel mate for this reflection in the data
C base.
      TYPE(REFLCT)=TYPE(REFLCT)*(MATRIX(-HH,-KK,-LL)*XRNSYM+ISYM)
      END IF
      END IF
C
      END DO
C
C
C We are storing the quantity ICEN * multiplicity in the
C array MULT.  This simplifies the calculation of epsilon
C throughout the code.  However, in order to get the actual
C multiplicity for a given reflection one needs to divide
C the MULT array by ICEN.
C
      IF (QHERM) THEN
      MULT(REFLCT)=ICEN*2*XRNSYM/MUL
      ELSE
      MULT(REFLCT)=ICEN*XRNSYM/MUL
      END IF
C
C
C compute epsilon (do a self-consistency test)
      MUL=0
      DO ISYM=1,XRNSYM
C generate symmetry related reflections
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      HH=XRITSY(ISYM,1,1)*H + XRITSY(ISYM,1,2)*K + XRITSY(ISYM,1,3)*L
      KK=XRITSY(ISYM,2,1)*H + XRITSY(ISYM,2,2)*K + XRITSY(ISYM,2,3)*L
      LL=XRITSY(ISYM,3,1)*H + XRITSY(ISYM,3,2)*K + XRITSY(ISYM,3,3)*L
      IF (H.EQ.HH.AND.K.EQ.KK.AND.L.EQ.LL) THEN
      MUL=MUL+1
      END IF
      END DO
C
      MUL=MUL/ICEN
C MUL is now epsilon
      IF (TYPE(REFLCT).GE.1.AND.QHERM) THEN
C acentric
      MUL2=2*XRNSYM/MULT(REFLCT)
      ELSE
C centric or not hermitian symmetry
      MUL2=XRNSYM/MULT(REFLCT)
      END IF
C
      IF (MUL.NE.MUL2) THEN
      WRITE(6,'(A,I6,I6)') ' XREDU4: internal error no. 6.',MUL,MUL2
      WRITE(6,'(3I6)') H,K,L
      END IF
C
      END DO
C
      RETURN
      END
C===================================================================
      SUBROUTINE XAVEFRIED(XRMREF,XRNREF,HPH,HPK,HPL,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPTYPE)
C
C Routine averages Bijvoet mates.  It operates on
C real, complex, integer, and HL reciprocal space objects
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XRMREF, XRNREF, HPH, HPK, HPL
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*), HPTYPE
C pointer
      INTEGER GRPPTR, GRPOBJ, MARK
C begin
C
      GRPPTR=ALLHP(INTEG4(XSFNUM+1))
      GRPOBJ=ALLHP(INTEG4(XSFNUM))
      MARK=ALLHP(INTEG4(XRNREF))
C
      CALL XAVEFRIE2(XRMREF,XRNREF,
     &           HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           HEAP(HPTYPE),
     &           HEAP(MARK),HEAP(GRPPTR),HEAP(GRPOBJ))
C
      CALL FREHP(MARK,INTEG4(XRNREF))
      CALL FREHP(GRPOBJ,INTEG4(XSFNUM))
      CALL FREHP(GRPPTR,INTEG4(XSFNUM+1))
C
      RETURN
      END
C===================================================================
      SUBROUTINE XAVEFRIE2(XRMREF,XRNREF,XRH,XRK,XRL,
     &           QHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           TYPE,
     &           MARK,GRPPTR,GRPOBJ)
C
C Routine averages Bijvoet mates.  It operates on
C real, complex, integer, and HL reciprocal space objects
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, XRNREF
      INTEGER XRH(*), XRK(*), XRL(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*), TYPE(*), MARK(*)
      INTEGER GRPPTR(*), GRPOBJ(*)
C local
      INTEGER REFLCT, TEMP, IISYM, R2, H, K, L, I
      LOGICAL FOUND
      DOUBLE PRECISION PHAS, RTH
      DOUBLE COMPLEX FDSHIFT, CTEMP, CAVE, CMATE
      DOUBLE PRECISION RTEMP, RAVE, RMATE
      INTEGER ITEMP, IAVE, IMATE
      INTEGER IGROUP, NGROUP, COUNTER, PA, PB, PC, PD
      DOUBLE PRECISION RPA, RPB, RPC, RPD
      DOUBLE PRECISION PAL, PBL, PCL, PDL
      DOUBLE PRECISION RMATEPA, RMATEPB, RMATEPC, RMATEPD
      DOUBLE PRECISION RAVEPA, RAVEPB, RAVEPC, RAVEPD
C parameter
      DOUBLE PRECISION ONE, TWO
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
C
      RTH=XRSYTH
C
C only required for anomalous data
      IF (.NOT.QHERM.AND.XRNREF.GT.0) THEN
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(2A)')
     &  ' XAVEFRIED: Friedel mates of all recipr. space obj.',
     &  ' will be averaged.'
      WRITE(6,'(A)')
     &  '    (also operates on phases and HL coefficients)'
      END IF
C
C get object HL group list
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
C
C initialize the mark array
      DO REFLCT=1,XRNREF
      MARK(REFLCT)=0
      END DO
C
      DO REFLCT=1,XRNREF
C
      IF (MARK(REFLCT).EQ.1) THEN
C
C reflection has already been averaged
      CONTINUE
      ELSEIF (TYPE(REFLCT).EQ.0) THEN
C
C keep centric reflections
      CONTINUE
      ELSEIF (ABS(TYPE(REFLCT)).EQ.1) THEN
C
C reflection without Bijvoet mate
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,3I6,A)')
     &  ' XAVEFRIED-warning: Bijvoet mate for reflection ',
     &  XRH(REFLCT),XRK(REFLCT),XRL(REFLCT),' not found.'
      END IF
C
      ELSEIF (TYPE(REFLCT).GT.1) THEN
C
C this is an acentric reflection with existing Bijvoet mate
C
C get index of Bijvoet mate
      FOUND=ABS(TYPE(REFLCT)).GT.1
      IF (FOUND) THEN
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(REFLCT))-1
      IISYM=MOD(TEMP,XRNSYM)+1
C
C R2 is the index of the Bijvoet mate of reflection I
      R2=TEMP/XRNSYM
C
C compute phase shift
      H=XRH(R2)
      K=XRK(R2)
      L=XRL(R2)
      PHAS=-TWO*PI*(( XRSYMM(IISYM,1,4)*H
     &               +XRSYMM(IISYM,2,4)*K
     &               +XRSYMM(IISYM,3,4)*L ) / RTH )
C
C this is the phase shift
      FDSHIFT=DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
C mark this reflections and its Bijvoet mate as done
      MARK(REFLCT)=1
      MARK(R2)=1
C
C
C loop over all reciprocal space objects except those
C that are in groups
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFGNAM(I).EQ.0) THEN
      IF (XSFTYPE(I).EQ.'COMP') THEN
C
C complex data type:
C   F+ -> [ F+  +   (F-)* ]/2
C   F- -> [ F-  +   (F+)* ]/2
C
      CALL XCOPY(CTEMP,1,HEAP(HPSF(I)),REFLCT)
      CALL XCOPY(CMATE,1,HEAP(HPSF(I)),R2)
      CAVE=(CTEMP+DCONJG(CMATE*FDSHIFT))/TWO
      CALL XCOPY(HEAP(HPSF(I)),REFLCT,CAVE,1)
C
CC CAVE=DCONJG(CAVE)/FDSHIFT
      CAVE=(CMATE+DCONJG(CTEMP)*DCONJG(FDSHIFT))/TWO
      CALL XCOPY(HEAP(HPSF(I)),R2,CAVE,1)
C
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
C
C real data type
C   F+ -> [ F+ + F- ]/2
C   F- -> [ F+ + F- ]/2
      CALL XCOPYR(RTEMP,1,HEAP(HPSF(I)),REFLCT)
      CALL XCOPYR(RMATE,1,HEAP(HPSF(I)),R2)
      RAVE=(RTEMP+RMATE)/TWO
      CALL XCOPYR(HEAP(HPSF(I)),REFLCT,RAVE,1)
      CALL XCOPYR(HEAP(HPSF(I)),R2,RAVE,1)
C
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
C
C integer data type
C   F+ -> [ F+ + F- ]/2
C   F- -> [ F+ + F- ]/2
      CALL XCOPYI(ITEMP,1,HEAP(HPSF(I)),REFLCT)
      CALL XCOPYI(IMATE,1,HEAP(HPSF(I)),R2)
      IAVE=(ITEMP+IMATE)/2
      CALL XCOPYI(HEAP(HPSF(I)),REFLCT,IAVE,1)
      CALL XCOPYI(HEAP(HPSF(I)),R2,IAVE,1)
      END IF
      END IF
      END DO
C
C loop over all HL groups
      DO IGROUP=1,NGROUP
C
C check if all elements are well-defined
      FOUND=.TRUE.
      DO COUNTER=GRPPTR(IGROUP)+1,GRPPTR(IGROUP+1)
      FOUND=FOUND.AND.HPSF(GRPOBJ(COUNTER)).NE.0
      END DO
C
      IF (FOUND) THEN
C
C map this set of HL coefficients
      PA=GRPOBJ(GRPPTR(IGROUP)+1)
      PB=GRPOBJ(GRPPTR(IGROUP)+2)
      PC=GRPOBJ(GRPPTR(IGROUP)+3)
      PD=GRPOBJ(GRPPTR(IGROUP)+4)
      CALL XCOPYR(RPA,1,HEAP(HPSF(PA)),REFLCT)
      CALL XCOPYR(RPB,1,HEAP(HPSF(PB)),REFLCT)
      CALL XCOPYR(RPC,1,HEAP(HPSF(PC)),REFLCT)
      CALL XCOPYR(RPD,1,HEAP(HPSF(PD)),REFLCT)
      CALL XCOPYR(RMATEPA,1,HEAP(HPSF(PA)),R2)
      CALL XCOPYR(RMATEPB,1,HEAP(HPSF(PB)),R2)
      CALL XCOPYR(RMATEPC,1,HEAP(HPSF(PC)),R2)
      CALL XCOPYR(RMATEPD,1,HEAP(HPSF(PD)),R2)
C
C apply phase shift to Bijvoet mate HL and complex conjugation
      PAL=RMATEPA
      PBL=RMATEPB
      PCL=RMATEPC
      PDL=RMATEPD
      CALL XMAPHL(1,PAL,PBL,PCL,PDL,PHAS,ONE,-ONE)
C
C average HL coefficients
      RAVEPA=(RPA+PAL)/TWO
      RAVEPB=(RPB+PBL)/TWO
      RAVEPC=(RPC+PCL)/TWO
      RAVEPD=(RPD+PDL)/TWO
C
C copy to primary reflection
      CALL XCOPYR(HEAP(HPSF(PA)),REFLCT,RAVEPA,1)
      CALL XCOPYR(HEAP(HPSF(PB)),REFLCT,RAVEPB,1)
      CALL XCOPYR(HEAP(HPSF(PC)),REFLCT,RAVEPC,1)
      CALL XCOPYR(HEAP(HPSF(PD)),REFLCT,RAVEPD,1)
C
C now do the same thing for the Bijvoet mate
C apply negative phase shift to primary reflection + compl. conjg
      PAL=RPA
      PBL=RPB
      PCL=RPC
      PDL=RPD
      CALL XMAPHL(1,PAL,PBL,PCL,PDL,-PHAS,-ONE,ONE)
C
C average HL coefficients
      RAVEPA=(RMATEPA+PAL)/TWO
      RAVEPB=(RMATEPB+PBL)/TWO
      RAVEPC=(RMATEPC+PCL)/TWO
      RAVEPD=(RMATEPD+PDL)/TWO
C
C copy to primary reflection
      CALL XCOPYR(HEAP(HPSF(PA)),R2,RAVEPA,1)
      CALL XCOPYR(HEAP(HPSF(PB)),R2,RAVEPB,1)
      CALL XCOPYR(HEAP(HPSF(PC)),R2,RAVEPC,1)
      CALL XCOPYR(HEAP(HPSF(PD)),R2,RAVEPD,1)
C
C
      END IF
C
      END DO
C
C
      END IF
      END IF
C
      END DO
C
      END IF
C
      RETURN
      END
C===================================================================
      SUBROUTINE XMAPHL(REFLCT,PA,PB,PC,PD,SHIFT,INVPRE,INVPOST)
C
C Routine applies a phase shift to Hendrickson-Lattman
C coefficients
C
      IMPLICIT NONE
C I/O
      INTEGER REFLCT
      DOUBLE PRECISION PA(*), PB(*), PC(*), PD(*), SHIFT
      DOUBLE PRECISION INVPRE, INVPOST
C local
      DOUBLE PRECISION PAL, PBL, PCL, PDL
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C begin
      PAL=PA(REFLCT)
      PBL=INVPRE*PB(REFLCT)
      PCL=PC(REFLCT)
      PDL=INVPRE*PD(REFLCT)
      PA(REFLCT)=PAL*COS(SHIFT)-PBL*SIN(SHIFT)
      PB(REFLCT)=INVPOST*(PAL*SIN(SHIFT)+PBL*COS(SHIFT))
      PC(REFLCT)=PCL*COS(TWO*SHIFT)-PDL*SIN(TWO*SHIFT)
      PD(REFLCT)=INVPOST*(PCL*SIN(TWO*SHIFT)+PDL*COS(TWO*SHIFT))
C
      RETURN
      END
