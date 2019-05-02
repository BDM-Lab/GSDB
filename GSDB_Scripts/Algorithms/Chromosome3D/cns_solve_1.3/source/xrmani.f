      SUBROUTINE XRRRR(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XNAMEMX,XSFMX,
     &           QHERM,XRRED,XRREUP,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRINTR,
     &           XRSYGP,XRSYIV)
C
C Routine reads and/or merges reflections
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, HPTSEL, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER XNAMEMX, XSFMX
      LOGICAL QHERM, XRRED, XRREUP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      INTEGER HMAX, KMAX, LMAX, HMIN, KMIN, LMIN
      INTEGER MDIM, MATRIX
C begin
C
C compute the maximum absolute values of the H,K,L indices
      CALL XRRR3(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HMAX,KMAX,LMAX,
     &                 HMIN,KMIN,LMIN)
C
      HMAX=MAX(1,ABS(HMAX),ABS(HMIN))
      KMAX=MAX(1,ABS(KMAX),ABS(KMIN))
      LMAX=MAX(1,ABS(LMAX),ABS(LMIN))
C
      MDIM=(2*HMAX+1)*(2*KMAX+1)*(2*LMAX+1)
      MATRIX=ALLHP(INTEG4(MDIM))
      CALL XRRR2(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           HMAX,KMAX,LMAX,
     &           HEAP(MATRIX),XNAMEMX,XSFMX,
     &           QHERM,XRRED,XRREUP,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRINTR,
     &           XRSYGP,XRSYIV)
      CALL FREHP(MATRIX,INTEG4(MDIM))
      RETURN
      END
C======================================================================
      SUBROUTINE XRRR2(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           HMAX,KMAX,LMAX,
     &           MATRIX,XNAMEMX,XSFMX,
     &           QHERM,XRRED,XRREUP,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRTR,XRINTR,
     &           XRSYGP,XRSYIV)
C
C actually parses information to read and/or merge reflections
C
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
      INTEGER XRMREF, XRNREF, HPTSEL
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT, HPTYPE
      INTEGER HMAX, KMAX, LMAX
      INTEGER MATRIX(-HMAX:HMAX,-KMAX:KMAX,-LMAX:LMAX)
      INTEGER XNAMEMX, XSFMX
      LOGICAL QHERM, XRRED, XRREUP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      INTEGER REFLCT, H, K, L, NEWREF, I, ITEMP, LPROMPT, LPROMPT2
      INTEGER NEWALL
      DOUBLE PRECISION AMPLTD, PHASE, AFOBS, PFOBS, TEMP
      DOUBLE COMPLEX CTEMP
      LOGICAL OK, ECHOLD, COND, OKL, QANOM, OLHERM
CCC modification ATB 4/27/08
      CHARACTER*100 PROMPT, PROMPT2
C parameter
      DOUBLE PRECISION RAD, ZERO, ONE
      INTEGER MARK
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, ONE=1.0D0)
      PARAMETER (MARK=-99999)
C
C setup merge matrix
      CALL XRFMX(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,-HMAX,-KMAX,-LMAX,MATRIX)
C
      IF (XRNREF.GT.0) THEN
      WRITE(6,'(A)')
     & ' REFLection: data will be merged with existing data.'
      END IF
C
      NEWREF=0
C
      ECHOLD=QECHO
C
      CALL PUSEND('REFLection>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('REFLection>')
      CALL MISCOM('REFLection>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-reflection')
C
C------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'RESE') THEN
C
      WRITE(6,'(A)')
     & ' XRRR: all existing reciprocal data will be scratched.'
C
      CALL XRAFRE
      CALL XRAREF(200,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
C
C
C reset merge matrix
      DO L=-LMAX,LMAX
      DO K=-KMAX,KMAX
      DO H=-HMAX,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
C------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      CALL XDECGROUP(XNAMEMX,XSFMX,0,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF)
C------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'DECL') THEN
      CALL XDECLARE(XNAMEMX,XSFMX,0,
     &   XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,HPSF,
     &   0,' ',0,0,
     &   0,0,XRNREF)
C------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'ANOM'.OR.WD(1:4).EQ.'HERM') THEN
C
      OLHERM=QHERM
      QANOM=.NOT.QHERM
      IF (WD(1:4).EQ.'ANOM') THEN
      CALL NEXTLO('ANOMalous=',QANOM)
      QHERM=.NOT.QANOM
      ELSE
      CALL NEXTLO('HERMitian=',QHERM)
      QANOM=.NOT.QHERM
      END IF
C
C reset flags only if changed
      IF ((OLHERM.AND..NOT.QHERM).OR.(.NOT.OLHERM.AND.QHERM)) THEN
C
C===> flags have changed
C check if reflections are present
      IF (XRNREF.GT.0) THEN
      IF (QANOM) THEN
C
C===> means that we have turned-on the anomalous flag
C===> need to expand the existing data set
      WRITE(6,'(A)')
     & ' XPARSE: expanding data set (h,k,l)->{(h,k,l),(-h,-k,-l)}'
C
C expand data set using Friedel operator
      CALL XEXPFRIED(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYMM,XRTR,XRINTR)
      ELSE
C
C===> means that we have turned-off the anomalous flag
C===> need to average the Bijvoet mates and reduce the data
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(2A)')
     &  ' XRRRR: data will be reduced to hemisphere. ',
     &  '         Anomalous signal will be lost.'
      END IF
C average Bijvoet mates
      CALL XAVEFRIED(XRMREF,XRNREF,HPH,HPK,HPL,
     &           OLHERM,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPTYPE)
C reduce reflections
      CALL XREDUC(XRRED,
     &           XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRSYGP,XRSYIV,XRTR,XRINTR)
      END IF
      END IF
C
C
C reset map database
      CALL XRMAPR(0)
C
C initialize reduction operation
      XRRED=.TRUE.
      XRREUP=.TRUE.
C
C setup merge matrix again
      CALL XRFMX(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,-HMAX,-KMAX,-LMAX,MATRIX)
C
      END IF
C------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'INDE') THEN
C
C turn off echo
      QECHO=.FALSE.
C
C parse miller indices
      CALL NEXTI('H=',H)
      CALL NEXTI('K=',K)
      CALL NEXTI('L=',L)
C
C check whether reflection contains already information that
C was previously read
      IF (H.GT.HMAX.OR.K.GT.KMAX.OR.L.GT.LMAX .OR.
     &    H.LT.-HMAX.OR.K.LT.-KMAX.OR.L.LT.-LMAX) THEN
      REFLCT=0
      ELSE
      REFLCT=MATRIX(H,K,L)
      END IF
C
      IF (REFLCT.EQ.0) THEN
C
C allocate new list index for that reflection
      IF (XRNREF.GE.XRMREF) THEN
      NEWALL=XRMREF+10000
      CALL XRAREF(NEWALL,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
C
      XRNREF=XRNREF+1
      NEWREF=NEWREF+1
      REFLCT=XRNREF
C
      CALL XCOPYI(HEAP(HPH),REFLCT,H,1)
      CALL XCOPYI(HEAP(HPK),REFLCT,K,1)
      CALL XCOPYI(HEAP(HPL),REFLCT,L,1)
C
C set target selection array
      CALL XCOPYI(HEAP(HPTSEL),REFLCT,1,1)
C
C set defaults for all reciprocal space objects
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(HEAP(HPSF(I)),REFLCT,DCMPLX(ZERO,ZERO),1)
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'REAL') THEN
      CALL XCOPYR(HEAP(HPSF(I)),REFLCT,ZERO,1)
      ELSEIF (HPSF(I).NE.0.AND.XSFTYPE(I).EQ.'INTE') THEN
      CALL XCOPYI(HEAP(HPSF(I)),REFLCT,0,1)
      END IF
      END DO
      END IF
C
C ok, now that we know the list index of the reflection we
C can parse the remaining information if any.
C
C
C check for information
      OK=.FALSE.
      DO WHILE (.NOT.OK)
      CALL NEXTWD('REFLection>')
C
C check reciprocal space objects
      COND=.FALSE.
      DO I=1,XSFNUM
      IF (.NOT.COND.AND.WD(1:WDLEN).EQ.XSFNAM(I)) THEN
      COND=.TRUE.
      IF (HPSF(I).EQ.0) CALL XSFAL(HPSF(I),XRMREF,XSFTYPE(I))
C
C special treatment for complex FOBS array
      IF (WD(1:WDLEN).EQ.'FOBS'.AND.XSFTYPE(I).EQ.'COMP') THEN
      PFOBS=ZERO
      CALL NEXTF('FOBS_amplitude=',AFOBS)
C is the next word a number?  If so, interpret it as the
C phase of FOBS.
      CALL NEXTWD('FOBS_phase (optional)=')
      CALL CHKNUMN(WD,WDLEN,OKL)
      IF (OKL) THEN
      TEMP=DECODF(WD,WDLEN,OKL)
      END IF
      IF (OKL) THEN
      PFOBS=TEMP
      ELSE
      CALL SAVEWD
      END IF
      CALL XAB(CTEMP,AFOBS,PFOBS*RAD)
      CALL XCOPY(HEAP(HPSF(I)),REFLCT,CTEMP,1)
C
      ELSEIF (XSFTYPE(I).EQ.'COMP') THEN
      PROMPT=WD(1:WDLEN)//'_amplitude='
      LPROMPT=WDLEN+11
      PROMPT2=WD(1:WDLEN)//'_phase='
      LPROMPT2=WDLEN+7
      CALL NEXTF(PROMPT(1:LPROMPT),AMPLTD)
      PROMPT=WD(1:WDLEN)//'_phase='
      CALL NEXTF(PROMPT2(1:LPROMPT2),PHASE)
      CALL XAB(CTEMP,AMPLTD,PHASE*RAD)
      CALL XCOPY(HEAP(HPSF(I)),REFLCT,CTEMP,1)
C
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      PROMPT=WD(1:WDLEN)//'_real='
      CALL NEXTF(PROMPT(1:WDLEN+6),AMPLTD)
      CALL XCOPYR(HEAP(HPSF(I)),REFLCT,AMPLTD,1)
C
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      PROMPT=WD(1:WDLEN)//'_integer='
      CALL NEXTI(PROMPT(1:WDLEN+9),ITEMP)
      CALL XCOPYI(HEAP(HPSF(I)),REFLCT,ITEMP,1)
      END IF
      END IF
      END DO
C
      IF (.NOT.COND.AND.WD(1:4).NE.'INDE'
     &             .AND.WD(1:3).NE.'END') THEN
      WRITE(6,'(2A)') ' %XRRR-err: unknown reciprocal space object: ',
     & WD(1:WDLEN)
      CALL WRNDIE(-5,'XRRR','unknown reciprocal space object.')
      ELSEIF (.NOT.COND) THEN
      CALL SAVEWD
      OK=.TRUE.
      END IF
C
      END DO
C
C============================================================
      ELSE
      CALL CHKEND('REFLection>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C re-establish old echo
      QECHO=ECHOLD
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I8,A)') ' XRRR2: ',NEWREF,
     &  ' new h,k,l indices have been added.'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRRR3(XRNREF,XRH,XRK,XRL,HMAX,KMAX,LMAX,
     &                 HMIN,KMIN,LMIN)
C
C Routine determines HMAX, KMAX, LMAX, HMIN, KMIN, LMIN
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*), HMAX, KMAX, LMAX
      INTEGER HMIN, KMIN, LMIN
C local
      INTEGER REFLCT
C begin
      IF (XRNREF.GT.0) THEN
      HMAX=XRH(1)
      KMAX=XRK(1)
      LMAX=XRL(1)
      HMIN=XRH(1)
      KMIN=XRK(1)
      LMIN=XRL(1)
      ELSE
      HMAX=0
      KMAX=0
      LMAX=0
      HMIN=0
      KMIN=0
      LMIN=0
      END IF
      DO REFLCT=2,XRNREF
      HMAX=MAX(HMAX,XRH(REFLCT))
      KMAX=MAX(KMAX,XRK(REFLCT))
      LMAX=MAX(LMAX,XRL(REFLCT))
      HMIN=MIN(HMIN,XRH(REFLCT))
      KMIN=MIN(KMIN,XRK(REFLCT))
      LMIN=MIN(LMIN,XRL(REFLCT))
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XCOPY(ARRAY,REFLCT,CTEMP,REFLCT2)
C
C Copies element into array (double complex version).
C
C Author Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE COMPLEX ARRAY(*), CTEMP(*)
      INTEGER REFLCT, REFLCT2
C local
C begin
      ARRAY(REFLCT)=CTEMP(REFLCT2)
      RETURN
      END
C======================================================================
      SUBROUTINE XCOPYR(ARRAY,REFLCT,TEMP,REFLCT2)
C
C Copies element into array (double precision version).
C
C Author Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION ARRAY(*), TEMP(*)
      INTEGER REFLCT, REFLCT2
C local
C begin
      ARRAY(REFLCT)=TEMP(REFLCT2)
      RETURN
      END
C======================================================================
      SUBROUTINE XCOPYI(ARRAY,REFLCT,TEMP,REFLCT2)
C
C Copies element into array (integer version).
C
C Author Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER ARRAY(*), TEMP(*)
      INTEGER REFLCT, REFLCT2
C local
C begin
      ARRAY(REFLCT)=TEMP(REFLCT2)
      RETURN
      END
C======================================================================
      SUBROUTINE XCOPYCH(ARRAY1, I1, ARRAY2, I2)
      IMPLICIT NONE
C
C Copies element into array (character version).
C
C I/O
      CHARACTER  ARRAY1(*)*(*), ARRAY2(*)*(*)
      INTEGER    I1, I2
C
C begin
      ARRAY1(I1) = ' '
      ARRAY1(I1) = ARRAY2(I2)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRFMX(XRNREF,XRH,XRK,XRL,
     &           HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,MATRIX)
C
C Fill book-keeping matrix with all reflections.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*), HMAX, KMAX, LMAX
      INTEGER HMIN, KMIN, LMIN
      INTEGER MATRIX(HMIN:HMAX,KMIN:KMAX,LMIN:LMAX)
C local
      INTEGER H, K, L, REFLCT
C begin
C
C setup merge matrix
      DO L=LMIN,LMAX
      DO K=KMIN,KMAX
      DO H=HMIN,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
C
      DO REFLCT=1,XRNREF
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      IF (H.GE.HMIN.AND.H.LE.HMAX.AND.
     &    K.GE.KMIN.AND.K.LE.KMAX.AND.
     &    L.GE.LMIN.AND.L.LE.LMAX) THEN
      MATRIX(H,K,L)=REFLCT
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRFMXF(XRNREF,XRH,XRK,XRL,
     &           HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,MATRIX,QHERM)
C
C Fill book-keeping matrix with all reflections and
C their Friedel-mates if QHERM=TRUE.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*), HMAX, KMAX, LMAX
      INTEGER HMIN, KMIN, LMIN
      INTEGER MATRIX(HMIN:HMAX,KMIN:KMAX,LMIN:LMAX)
      LOGICAL QHERM
C local
      INTEGER H, K, L, REFLCT
C begin
C
C setup merge matrix
      DO L=LMIN,LMAX
      DO K=KMIN,KMAX
      DO H=HMIN,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
C
      DO REFLCT=1,XRNREF
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      IF (H.GE.HMIN.AND.H.LE.HMAX.AND.
     &    K.GE.KMIN.AND.K.LE.KMAX.AND.
     &    L.GE.LMIN.AND.L.LE.LMAX) THEN
      MATRIX(H,K,L)=REFLCT
      IF (QHERM) THEN
      MATRIX(-H,-K,-L)=REFLCT
      END IF
      ELSE
      CALL WRNDIE(-5,'XRFMX','fatal coding error')
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRFMXS(XRNREF,XRH,XRK,XRL,
     &           HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,MATRIX,
     &           XRNSYM,XRMSYM,XRSYMM)
C
C Fill book-keeping matrix with all reflections
C and all symmetry mates.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*), HMAX, KMAX, LMAX
      INTEGER HMIN, KMIN, LMIN
      INTEGER MATRIX(HMIN:HMAX,KMIN:KMAX,LMIN:LMAX)
      INTEGER XRNSYM, XRMSYM, XRSYMM(XRMSYM,3,4)
C local
      INTEGER H, K, L, HH, KK, LL, REFLCT, ISYM
C begin
C
C setup merge matrix
      DO L=LMIN,LMAX
      DO K=KMIN,KMAX
      DO H=HMIN,HMAX
      MATRIX(H,K,L)=0
      END DO
      END DO
      END DO
C
      DO ISYM=1,XRNSYM
C
      DO REFLCT=1,XRNREF
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
C apply transpose of symmetry operator to H, K, L
      HH=XRSYMM(ISYM,1,1)*H+XRSYMM(ISYM,2,1)*K+XRSYMM(ISYM,3,1)*L
      KK=XRSYMM(ISYM,1,2)*H+XRSYMM(ISYM,2,2)*K+XRSYMM(ISYM,3,2)*L
      LL=XRSYMM(ISYM,1,3)*H+XRSYMM(ISYM,2,3)*K+XRSYMM(ISYM,3,3)*L
      IF (HH.GE.HMIN.AND.HH.LE.HMAX.AND.
     &    KK.GE.KMIN.AND.KK.LE.KMAX.AND.
     &    LL.GE.LMIN.AND.LL.LE.LMAX) THEN
      MATRIX(HH,KK,LL)=REFLCT
      ELSE
      CALL WRNDIE(-5,'XRFMXS','fatal coding error')
      END IF
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XGENER(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRTR,XRINTR,XRHIGH,XRLOW)
C
C This routine complements all reflections to
C a full asymmetric unit for the selected resolution range.
C The asymmetric unit is defined in subroutine XRASYM.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER XRMREF, HPTSEL, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRHIGH, XRLOW
C local
      INTEGER HMAX, KMAX, LMAX, MDIM
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
C pointers
      INTEGER MATRIX
C parameters
      DOUBLE PRECISION M0001
      PARAMETER (M0001=0.0001D0)
C begin
C
      IF (XRNREF.GT.0) THEN
C first we reduce all reflections to an asymmetric unit
      WRITE(6,'(A)')
     & ' XGENER: existing reflections will be conserved',
     & '         appending new reflections to produce',
     & '         a full set for the specified resolution range'
      ELSE
      WRITE(6,'(A)')
     & ' XGENER: generating reflections to produce a full set',
     & '         for the specified resolution range.'
      END IF
C
C determine HMAX, ...
C the following condition can be derived as follows:
C
C we have
C    s = [F]t h
C where [F]t is the transpose of the orthogonal
C to fractional coordinate transformation matrix.
C Thus,
C    h = [F]t-1 s
C where [F]t-1 is the inverse of the transpose of F.
C The high resolution limit condition is
C    s < sh
C With sh = 1/dhigh.  This implies
C    s[1] < sh, s[2] < sh, s[3] < sh
C for the components s[1], s[2], s[3] of the vector s.
C We also have
C    ([F]t-1)[ij] >= 0
C for all components of the matrix [F]t-1.
C Thus,
C   h < [F]t-1 * [sh,sh,sh]t.
C
C
      TEMPX=SQRT(XRINTR(1,1)**2+XRINTR(2,1)**2+XRINTR(3,1)**2)*XRHIGH
      TEMPY=SQRT(XRINTR(1,2)**2+XRINTR(2,2)**2+XRINTR(3,2)**2)*XRHIGH
      TEMPZ=SQRT(XRINTR(1,3)**2+XRINTR(2,3)**2+XRINTR(3,3)**2)*XRHIGH
C
      IF (TEMPX.LT.M0001.OR.TEMPY.LT.M0001.OR.TEMPZ.LT.M0001) THEN
      CALL WRNDIE(-5,'XGENER',
     & 'high resolution limit too small or unit cell too large')
      END IF
      HMAX=INT(R4SMAL+TEMPX)+1
      KMAX=INT(R4SMAL+TEMPY)+1
      LMAX=INT(R4SMAL+TEMPZ)+1
C
C
C allocate space for the book-keeping matrix
      MDIM=(2*HMAX+1)*(2*KMAX+1)*(2*LMAX+1)
      MATRIX=ALLHP(INTEG4(MDIM))
C
C call routine that actually does the work
      CALL XGENE2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFTYPE,HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRTR,XRHIGH,XRLOW,
     &           HMAX,KMAX,LMAX,HEAP(MATRIX),
     &           HPMULT,HPTYPE,HPH,HPK,HPL)
C
C free space for book-keeping matrix
      CALL FREHP(MATRIX,INTEG4(MDIM))
      RETURN
      END
C======================================================================
      SUBROUTINE XGENE2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFTYPE,HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           XRTR,XRHIGH,XRLOW,
     &           HMAX,KMAX,LMAX,MATRIX,
     &           HPMULT,HPTYPE,HPH,HPK,HPL)
C
C See routine XGENER above.
C Author:  Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, XRNREF, HPTSEL
      INTEGER XSFNUM
      CHARACTER*(*) XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3), XRHIGH, XRLOW
      INTEGER HMAX, KMAX, LMAX
      INTEGER MATRIX(-HMAX:HMAX,-KMAX:KMAX,-LMAX:LMAX)
      INTEGER HPMULT, HPTYPE, HPH, HPK, HPL
C local
      INTEGER HH, KK, LL, H, K, L, IISYM, IFRIED, I
      INTEGER NEWREF, NEWALL
      LOGICAL SYSAB
      DOUBLE PRECISION  SSQ, D
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C
C setup merge matrix
      CALL XRFMX(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,-HMAX,-KMAX,-LMAX,MATRIX)
C
C
      NEWREF=0
C
C loop through the box determined by h,k,l max.
      DO H=-HMAX,HMAX
      DO K=-KMAX,KMAX
      DO L=-LMAX,LMAX
C
C compute s**2 for this reflection
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
C
C check resolution limits
      D=SQRT(SSQ)
      IF (D.LT.XRHIGH.AND.D.GT.XRLOW) THEN
C
C map this reflection into the asymmetric unit (returned in HH, KK, LL)
      CALL XRASYM(H,K,L,HH,KK,LL,IISYM,IFRIED,
     &                  XRNSYM,XRMSYM,XRITSY,QHERM)
C
C check if this is a systematic absence
      CALL XRSYSAB(HH,KK,LL,SYSAB,
     &                  XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C check consistency
      IF (HH.LT.-HMAX.OR.HH.GT.HMAX.OR.
     &    KK.LT.-KMAX.OR.KK.GT.KMAX.OR.
     &    LL.LT.-LMAX.OR.LL.GT.LMAX) THEN
      CALL WRNDIE(-5,'XGENE2','fatal coding error 2')
      END IF
C
C check whether this reflection is already present
C omit systematic absences
      IF (MATRIX(HH,KK,LL).EQ.0.AND..NOT.SYSAB) THEN
C
C add this reflection.
C
C allocate new list index for that reflection
      IF (XRNREF.GE.XRMREF) THEN
      NEWALL=XRMREF+10000
      CALL XRAREF(NEWALL,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
      NEWREF=NEWREF+1
      XRNREF=XRNREF+1
C
      CALL XCOPYI(HEAP(HPH),XRNREF,HH,1)
      CALL XCOPYI(HEAP(HPK),XRNREF,KK,1)
      CALL XCOPYI(HEAP(HPL),XRNREF,LL,1)
C
C set target selection array
      CALL XCOPYI(HEAP(HPTSEL),XRNREF,1,1)
C
C set reciprocal space objects to zero
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0) THEN
      IF (XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(HEAP(HPSF(I)),XRNREF,DCMPLX(ZERO,ZERO),1)
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      CALL XCOPYR(HEAP(HPSF(I)),XRNREF,ZERO,1)
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      CALL XCOPYI(HEAP(HPSF(I)),XRNREF,0,1)
      END IF
      END IF
      END DO
C
      CALL XCOPYI(HEAP(HPMULT),XRNREF,1,1)
      CALL XCOPYI(HEAP(HPTYPE),XRNREF,1,1)
C
C
C
C update book-keeping matrix
      MATRIX(HH,KK,LL)=XRNREF
      END IF
      END IF
      END DO
      END DO
      END DO
C
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I8,A)') ' XGENE2: ',NEWREF,
     &  ' new reflections have been generated.'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XRASYM(H,K,L,HH,KK,LL,IISYM,IFRIED,
     &                  XRNSYM,XRMSYM,XRITSY,QHERM)
C
C Maps reflection H,K,L into the asymmetric unit
C (returned in HH, KK, LL).
C
C An order is defined for two symmetry related reflections (r1, r2).
C First, both reflections are mapped into a hemisphere if hermitian
C symmetry is turned on (see FUNCTION RECHEM).  Then the order is
C decided by FUNCTION CMPHKL.
C
C All symmetry related reflections of H,K,L are generated.  The
C reflection which is largest according to the above definition
C is returned in HH, KK, LL.  The largest reflection is determined
C by a binary tree search.
C
C IISYM returns the number of the symmetry operator that transforms
C H,K,L into HH,KK,LL and IRIEDL=-1 indicates whether HH,KK,LL is
C obtained by application of the hermitian symmetry.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER H, K, L, HH, KK, LL, IISYM, IFRIED
      INTEGER XRNSYM, XRMSYM
      INTEGER XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
C
      INTEGER   RECHEM, CMPHKL
      EXTERNAL  RECHEM, CMPHKL
C local
      INTEGER ISYM, HHH, KKK, LLL, IIFRID
C parameter
C
C go through all symmetry related reflections and take the one
C that is in the asymmetric unit.
      IISYM=1
      IFRIED=1
      DO ISYM=1,XRNSYM
C
C apply transpose of the inverse of the symmetry opertor
      HHH=XRITSY(ISYM,1,1)*H + XRITSY(ISYM,1,2)*K + XRITSY(ISYM,1,3)*L
      KKK=XRITSY(ISYM,2,1)*H + XRITSY(ISYM,2,2)*K + XRITSY(ISYM,2,3)*L
      LLL=XRITSY(ISYM,3,1)*H + XRITSY(ISYM,3,2)*K + XRITSY(ISYM,3,3)*L
C
C if hermitian symmetry is present then map the reflection into the
C l>0 hemisphere
      IF (QHERM .AND. RECHEM(HHH, KKK, LLL) .LT. 0) THEN
      HHH=-HHH
      KKK=-KKK
      LLL=-LLL
      IIFRID=-1
      ELSE
      IIFRID=+1
      END IF
C
      IF (ISYM.EQ.1) THEN
      HH=HHH
      KK=KKK
      LL=LLL
      IFRIED=IIFRID
      ELSE
      IF (CMPHKL(HH, KK, LL, HHH, KKK, LLL) .GT. 0) THEN
      HH=HHH
      KK=KKK
      LL=LLL
      IISYM=ISYM
      IFRIED=IIFRID
      END IF
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XRSYSAB(H,K,L,SYSAB,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C
C Checks if reflection H,K,L represents a systematic
C absence.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER H, K, L
      LOGICAL SYSAB
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER ISYM
      INTEGER HH, KK, LL
      DOUBLE PRECISION RTH, SHIFT, SINSFT, COSSFT
C parameter
      DOUBLE PRECISION ONE, TWO
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C
      SYSAB=.FALSE.
C loop over all symmetry operators
      RTH=XRSYTH
      DO ISYM=1,XRNSYM
C
C apply transpose of inverse symmetry operator to H, K, L
      HH=XRITSY(ISYM,1,1)*H+XRITSY(ISYM,1,2)*K+XRITSY(ISYM,1,3)*L
      KK=XRITSY(ISYM,2,1)*H+XRITSY(ISYM,2,2)*K+XRITSY(ISYM,2,3)*L
      LL=XRITSY(ISYM,3,1)*H+XRITSY(ISYM,3,2)*K+XRITSY(ISYM,3,3)*L
C
C compute phase shift of translational part of symmetry operator,
C use TRANSFORMED indices
      SHIFT=TWO*PI*(( XRSYMM(ISYM,1,4)*HH
     &              +XRSYMM(ISYM,2,4)*KK
     &              +XRSYMM(ISYM,3,4)*LL ) / RTH )
      COSSFT=DCOS(SHIFT)
      SINSFT=DSIN(SHIFT)
C
C test for systematic absence
      IF (HH.EQ.H.AND.
     &    KK.EQ.K.AND.
     &    LL.EQ.L.AND.
     &    (ABS(COSSFT-ONE).GT.R4SMAL.OR.ABS(SINSFT).GT.R4SMAL)) THEN
      SYSAB=.TRUE.
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XEXPAN(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRTR,XRINTR)
C
C This routine applies the current symmetry operators to all
C reflections and appends the symmetry-related reflections
C to the reflection list (except when they're related by
C Friedel symmetry to existing reflections in QHERM=TRUE)
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER XRMREF, HPTSEL, XRNREF
      INTEGER HPH, HPK, HPL
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
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
C local
      INTEGER HMAX, KMAX, LMAX, MDIM
      DOUBLE PRECISION RMIN, TEMPX, TEMPY, TEMPZ
C parameters
      DOUBLE PRECISION M0001
      PARAMETER (M0001=0.0001D0)
C pointers
      INTEGER MATRIX
      INTEGER GRPPTR, GRPOBJ
C begin
C
      IF (XRNREF.GT.0) THEN
C
      CALL XEXPA3(RMIN,XRNREF,HEAP(HPH),HEAP(HPL),HEAP(HPK),XRTR)
C
      TEMPX=SQRT(XRINTR(1,1)**2+XRINTR(2,1)**2+XRINTR(3,1)**2)*RMIN
      TEMPY=SQRT(XRINTR(1,2)**2+XRINTR(2,2)**2+XRINTR(3,2)**2)*RMIN
      TEMPZ=SQRT(XRINTR(1,3)**2+XRINTR(2,3)**2+XRINTR(3,3)**2)*RMIN
C
      IF (TEMPX.LT.M0001.OR.TEMPY.LT.M0001.OR.TEMPZ.LT.M0001) THEN
      CALL WRNDIE(-5,'XEXPAN',
     & 'high resolution limit too small or unit cell too large')
      END IF
      HMAX=INT(R4SMAL+TEMPX)+1
      KMAX=INT(R4SMAL+TEMPY)+1
      LMAX=INT(R4SMAL+TEMPZ)+1
C
C
C allocate space for the book-keeping matrix
      MDIM=(2*HMAX+1)*(2*KMAX+1)*(2*LMAX+1)
      MATRIX=ALLHP(INTEG4(MDIM))
C
      GRPPTR=ALLHP(INTEG4(XSFNUM+1))
      GRPOBJ=ALLHP(INTEG4(XSFNUM))
C
C call the routine that actually does the work
      CALL XEXPA2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           HMAX,KMAX,LMAX,HEAP(MATRIX),
     &           HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HEAP(GRPPTR), HEAP(GRPOBJ))
C
      CALL FREHP(GRPOBJ,INTEG4(XSFNUM))
      CALL FREHP(GRPPTR,INTEG4(XSFNUM+1))
C free space for book-keeping matrix
      CALL FREHP(MATRIX,INTEG4(MDIM))
      END IF
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XEXPA2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &           HMAX,KMAX,LMAX,MATRIX,
     &           HPH,HPK,HPL,HPMULT,HPTYPE,
     &           GRPPTR,GRPOBJ)
C
C =====================================================================
C For each symmetry operator in real space
C     rho(S*r+d) = rho(r)
C the corresponding symmetry operator in Fourier space is
C given by
C     F( T*H )= F(H) exp(i T*H*d)
C where T is the transpose of the inverse of the symmetry operator S.
C
C The hermitian symmetry operator is simply given by
C     F(h,k,l)=F*(-h,-k,-l)
C =====================================================================
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, XRNREF, HPTSEL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER HMAX, KMAX, LMAX
      INTEGER MATRIX(-HMAX:HMAX,-KMAX:KMAX,-LMAX:LMAX)
      INTEGER HPH, HPK, HPL, HPTYPE, HPMULT
      INTEGER GRPPTR(*), GRPOBJ(*)
C local
      INTEGER ISYM, NEWREF, I, TT
      INTEGER REFLCT, H, K, L, HH, KK, LL, XROREF
      DOUBLE PRECISION RTH, PHAS
      DOUBLE COMPLEX SHIFT, CTEMP
      INTEGER NEWALL
      INTEGER PA, PB, PC, PD, NGROUP, IGROUP, COUNTER
      LOGICAL FOUND
C parameter
      DOUBLE PRECISION ONE, TWO, FOUR, S180
      PARAMETER (ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0, S180=180.D0)
C begin
C
C setup book-keeping matrix for all reflections and all Friedel
C mates if QHERM=TRUE
      CALL XRFMXF(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,-HMAX,-KMAX,-LMAX,MATRIX,QHERM)
C
C
C get object group list
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
      NEWREF=0
      XROREF=XRNREF
C
C loop over all symmetry operators
      RTH=XRSYTH
      DO ISYM=1,XRNSYM
C
C loop over all reflections
      DO REFLCT=1,XROREF
C
C apply transpose of inverse symmetry operator to H, K, L
C
      CALL XCOPYI(HH,1,HEAP(HPH),REFLCT)
      CALL XCOPYI(KK,1,HEAP(HPK),REFLCT)
      CALL XCOPYI(LL,1,HEAP(HPL),REFLCT)
C
      H=XRITSY(ISYM,1,1)*HH+XRITSY(ISYM,1,2)*KK+XRITSY(ISYM,1,3)*LL
      K=XRITSY(ISYM,2,1)*HH+XRITSY(ISYM,2,2)*KK+XRITSY(ISYM,2,3)*LL
      L=XRITSY(ISYM,3,1)*HH+XRITSY(ISYM,3,2)*KK+XRITSY(ISYM,3,3)*LL
C
C compute phase shift of translational part of symmetry operator,
C use TRANSFORMED indices
      PHAS=TWO*PI*(( XRSYMM(ISYM,1,4)*H
     &              +XRSYMM(ISYM,2,4)*K
     &              +XRSYMM(ISYM,3,4)*L ) / RTH )
C
C check whether this reflection (or its Friedel mate if QHERM=TRUE)
C is already present
      IF (MATRIX(H,K,L).EQ.0) THEN
C
C add this reflection.
      IF (XRNREF.GE.XRMREF) THEN
      NEWALL=XRMREF+10000
      CALL XRAREF(NEWALL,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
C
      NEWREF=NEWREF+1
      XRNREF=XRNREF+1
C
C complex-multiply everything and store factor in temporary vector FLOC
      SHIFT=DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
      CALL XCOPYI(HEAP(HPH),XRNREF,H,1)
      CALL XCOPYI(HEAP(HPK),XRNREF,K,1)
      CALL XCOPYI(HEAP(HPL),XRNREF,L,1)
C
      CALL XCOPYI(TT,1,HEAP(HPTSEL),REFLCT)
      CALL XCOPYI(HEAP(HPTSEL),XRNREF,TT,1)
C
      CALL XCOPYI(TT,1,HEAP(HPMULT),REFLCT)
      CALL XCOPYI(HEAP(HPMULT),XRNREF,TT,1)
      CALL XCOPYI(TT,1,HEAP(HPTYPE),REFLCT)
      CALL XCOPYI(HEAP(HPTYPE),XRNREF,TT,1)
C
C map reciprocal space objects except grouped arrays
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFGNAM(I).EQ.0) THEN
      IF (XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(I)),REFLCT)
      CTEMP=CTEMP*SHIFT
      CALL XCOPY(HEAP(HPSF(I)),XRNREF,CTEMP,1)
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      CALL XCOPYR(HEAP(HPSF(I)),XRNREF,HEAP(HPSF(I)),REFLCT)
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      CALL XCOPYI(HEAP(HPSF(I)),XRNREF,HEAP(HPSF(I)),REFLCT)
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
      CALL XCOPYR(HEAP(HPSF(PA)),XRNREF,HEAP(HPSF(PA)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PB)),XRNREF,HEAP(HPSF(PB)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PC)),XRNREF,HEAP(HPSF(PC)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PD)),XRNREF,HEAP(HPSF(PD)),REFLCT)
      CALL XMAPHL(XRNREF,HEAP(HPSF(PA)),HEAP(HPSF(PB)),
     &           HEAP(HPSF(PC)),HEAP(HPSF(PD)),PHAS,
     &           ONE,ONE)
      END IF
C
      END DO
C
C update the book-keeping matrix
      IF (H.LT.-HMAX.OR.H.GT.HMAX.OR.
     &    K.LT.-KMAX.OR.K.GT.KMAX.OR.
     &    L.LT.-LMAX.OR.L.GT.LMAX) THEN
      CALL WRNDIE(-5,'XEXPA2','fatal coding error')
      ELSE
      MATRIX(H,K,L)=XRNREF
      IF (QHERM) MATRIX(-H,-K,-L)=XRNREF
      END IF
C
      END IF
      END DO
      END DO
C
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I8,A)') ' XEXPA2: ',NEWREF,
     &  ' new reflections have been generated.'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XEXPFRIED(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRNSYM,XRMSYM,XRSYMM,XRTR,XRINTR)
C
C This routine applies the hermitian operator to all
C reflections and appends the Friedel-related reflections
C to the reflection list.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER XRMREF, HPTSEL, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE, XRNSYM, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4)
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3)
C local
      INTEGER HMAX, KMAX, LMAX, MDIM
      DOUBLE PRECISION RMIN, TEMPX, TEMPY, TEMPZ
C parameters
      DOUBLE PRECISION M0001
      PARAMETER (M0001=0.0001D0)
C pointers
      INTEGER MATRIX
      INTEGER GRPPTR, GRPOBJ
C begin
C
      IF (XRNREF.GT.0) THEN
C
      CALL XEXPA3(RMIN,XRNREF,HEAP(HPH),HEAP(HPL),HEAP(HPK),XRTR)
C
      TEMPX=SQRT(XRINTR(1,1)**2+XRINTR(2,1)**2+XRINTR(3,1)**2)*RMIN
      TEMPY=SQRT(XRINTR(1,2)**2+XRINTR(2,2)**2+XRINTR(3,2)**2)*RMIN
      TEMPZ=SQRT(XRINTR(1,3)**2+XRINTR(2,3)**2+XRINTR(3,3)**2)*RMIN
C
      IF (TEMPX.LT.M0001.OR.TEMPY.LT.M0001.OR.TEMPZ.LT.M0001) THEN
      CALL WRNDIE(-5,'XEXPFRIED',
     & 'high resolution limit too small or unit cell too large')
      END IF
      HMAX=INT(R4SMAL+TEMPX)+1
      KMAX=INT(R4SMAL+TEMPY)+1
      LMAX=INT(R4SMAL+TEMPZ)+1
C
C
C allocate space for the book-keeping matrix
      MDIM=(2*HMAX+1)*(2*KMAX+1)*(2*LMAX+1)
      MATRIX=ALLHP(INTEG4(MDIM))
C
      GRPPTR=ALLHP(INTEG4(XSFNUM+1))
      GRPOBJ=ALLHP(INTEG4(XSFNUM))
C
C call the routine that actually does the work
      CALL XEXPFR2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYMM,
     &           HMAX,KMAX,LMAX,HEAP(MATRIX),
     &           HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HEAP(GRPPTR), HEAP(GRPOBJ))
C
      CALL FREHP(GRPOBJ,INTEG4(XSFNUM))
      CALL FREHP(GRPPTR,INTEG4(XSFNUM+1))
C free space for book-keeping matrix
      CALL FREHP(MATRIX,INTEG4(MDIM))
      END IF
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XEXPFR2(XRMREF,XRNREF,HPTSEL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,
     &           XRNSYM,XRMSYM,XRSYMM,
     &           HMAX,KMAX,LMAX,MATRIX,
     &           HPH,HPK,HPL,HPMULT,HPTYPE,
     &           GRPPTR,GRPOBJ)
C
C =====================================================================
C For each symmetry operator in real space
C     rho(S*r+d) = rho(r)
C the corresponding symmetry operator in Fourier space is
C given by
C     F( T*H )= F(H) exp(i T*H*d)
C where T is the transpose of the inverse of the symmetry operator S.
C
C The hermitian symmetry operator is simply given by
C     F(h,k,l)=F*(-h,-k,-l)
C =====================================================================
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER XRMREF, XRNREF, HPTSEL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER XRNSYM, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4)
      INTEGER HMAX, KMAX, LMAX
      INTEGER MATRIX(-HMAX:HMAX,-KMAX:KMAX,-LMAX:LMAX)
      INTEGER HPH, HPK, HPL, HPTYPE, HPMULT
      INTEGER GRPPTR(*), GRPOBJ(*)
C local
      INTEGER NEWREF, I, TT, ISYM, HH, KK, LL
      INTEGER REFLCT, H, K, L, XROREF
      DOUBLE COMPLEX CTEMP
      INTEGER NEWALL
      INTEGER PA, PB, PC, PD, NGROUP, IGROUP, COUNTER
      LOGICAL FOUND
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
C begin
C
C setup book-keeping matrix for all reflections and all symmetry mates
      CALL XRFMXS(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,-HMAX,-KMAX,-LMAX,MATRIX,
     &           XRNSYM,XRMSYM,XRSYMM)
C
C
C get object group list
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
      NEWREF=0
C
      XROREF=XRNREF
C
C loop over all reflections
      DO REFLCT=1,XROREF
C
      CALL XCOPYI(H,1,HEAP(HPH),REFLCT)
      CALL XCOPYI(K,1,HEAP(HPK),REFLCT)
      CALL XCOPYI(L,1,HEAP(HPL),REFLCT)
C
C generate the Friedel mate
      H=-H
      K=-K
      L=-L
C
C check whether this reflection (or any of its symmetry mates)
C is already present
      IF (MATRIX(H,K,L).EQ.0) THEN
C
C add this reflection.
      IF (XRNREF.GE.XRMREF) THEN
      NEWALL=XRMREF+10000
      CALL XRAREF(NEWALL,XRMREF,XRNREF,HPTSEL,
     &                  HPH,HPK,HPL,HPMULT,HPTYPE,
     &                  XSFNUM,HPSF,XSFTYPE)
      END IF
C
      NEWREF=NEWREF+1
      XRNREF=XRNREF+1
C
      CALL XCOPYI(HEAP(HPH),XRNREF,H,1)
      CALL XCOPYI(HEAP(HPK),XRNREF,K,1)
      CALL XCOPYI(HEAP(HPL),XRNREF,L,1)
C
      CALL XCOPYI(TT,1,HEAP(HPTSEL),REFLCT)
      CALL XCOPYI(HEAP(HPTSEL),XRNREF,TT,1)
C
      CALL XCOPYI(TT,1,HEAP(HPMULT),REFLCT)
      CALL XCOPYI(HEAP(HPMULT),XRNREF,TT,1)
      CALL XCOPYI(TT,1,HEAP(HPTYPE),REFLCT)
      CALL XCOPYI(HEAP(HPTYPE),XRNREF,TT,1)
C
C map reciprocal space objects except grouped arrays
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0.AND.XSFGNAM(I).EQ.0) THEN
      IF (XSFTYPE(I).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(I)),REFLCT)
      CTEMP=DCONJG(CTEMP)
      CALL XCOPY(HEAP(HPSF(I)),XRNREF,CTEMP,1)
      ELSEIF (XSFTYPE(I).EQ.'REAL') THEN
      CALL XCOPYR(HEAP(HPSF(I)),XRNREF,HEAP(HPSF(I)),REFLCT)
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      CALL XCOPYI(HEAP(HPSF(I)),XRNREF,HEAP(HPSF(I)),REFLCT)
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
      CALL XCOPYR(HEAP(HPSF(PA)),XRNREF,HEAP(HPSF(PA)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PB)),XRNREF,HEAP(HPSF(PB)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PC)),XRNREF,HEAP(HPSF(PC)),REFLCT)
      CALL XCOPYR(HEAP(HPSF(PD)),XRNREF,HEAP(HPSF(PD)),REFLCT)
      CALL XMAPHL(XRNREF,HEAP(HPSF(PA)),HEAP(HPSF(PB)),
     &           HEAP(HPSF(PC)),HEAP(HPSF(PD)),ZERO,
     &           ONE,-ONE)
      END IF
C
      END DO
C
C update the book-keeping matrix
      DO ISYM=1,XRNSYM
C apply transpose of symmetry operator to H, K, L
      HH=XRSYMM(ISYM,1,1)*H+XRSYMM(ISYM,2,1)*K+XRSYMM(ISYM,3,1)*L
      KK=XRSYMM(ISYM,1,2)*H+XRSYMM(ISYM,2,2)*K+XRSYMM(ISYM,3,2)*L
      LL=XRSYMM(ISYM,1,3)*H+XRSYMM(ISYM,2,3)*K+XRSYMM(ISYM,3,3)*L
      IF (HH.LT.-HMAX.OR.HH.GT.HMAX.OR.
     &    KK.LT.-KMAX.OR.KK.GT.KMAX.OR.
     &    LL.LT.-LMAX.OR.LL.GT.LMAX) THEN
      CALL WRNDIE(-5,'XEXPFR2','fatal coding error')
      ELSE
      MATRIX(HH,KK,LL)=XRNREF
      END IF
      END DO
      END IF
      END DO
C
C modification, ATB, 12/04/08
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,I8,A)') ' XEXPFR2: ',NEWREF,
     &  ' new reflections have been generated.'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XEXPA3(RMIN,XRNREF,XRH,XRL,XRK,XRTR)
C
C Routine determines rmin=1/dmin
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION RMIN
      INTEGER XRNREF, XRH(*), XRL(*), XRK(*)
      DOUBLE PRECISION XRTR(3,3)
C local
      INTEGER REFLCT, H, K, L
      DOUBLE PRECISION D, SSQ
C parameter
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D0, ZERO=0.0D0)
C begin
C
C determine minimum resolution
      DO REFLCT=1,XRNREF
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      IF (SSQ.GT.RSMALL) THEN
      D=SQRT(SSQ)
      ELSE
      D=ZERO
      END IF
      IF (REFLCT.EQ.1) THEN
      RMIN=D
      ELSE
      RMIN=MAX(RMIN,D)
      END IF
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XWRIT(XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XSFMX,XRCELL,XRVOL)
C
C Routine writes reflection files
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER XSFMX
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
C pointer
      INTEGER QSELE, IOBJ
      INTEGER GRPPTR, GRPOBJ
C begin
C
      QSELE=ALLHP(ILOGIC(XRNREF))
      IOBJ=ALLHP(INTEG4(XSFNUM))
      GRPPTR=ALLHP(INTEG4(XSFMX+1))
      GRPOBJ=ALLHP(INTEG4(XSFMX))
      CALL XWRIT2(XRMREF,XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,HEAP(QSELE),
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,HEAP(IOBJ),
     &           HEAP(GRPPTR), HEAP(GRPOBJ),XRCELL,XRVOL)
      CALL FREHP(GRPOBJ,INTEG4(XSFMX))
      CALL FREHP(GRPPTR,INTEG4(XSFMX+1))
      CALL FREHP(IOBJ,INTEG4(XSFNUM))
      CALL FREHP(QSELE,ILOGIC(XRNREF))
      RETURN
      END
C======================================================================
      SUBROUTINE XWRIT2(XRMREF,XRNREF,XRH,XRK,XRL,
     &           HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,XSFGTYP,XSFGORD,
     &           HPSF,HPMULT,HPTYPE,
     &           XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QSELE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,IOBJ,
     &           GRPPTR,GRPOBJ,XRCELL,XRVOL)
C
C Routine writes reflection files
C
C See routine XWRIT above
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XRMREF, XRNREF
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER XSFGNAM(*)
      CHARACTER*(*) XSFGTYP(*)
      INTEGER XSFGORD(*)
      INTEGER HPSF(*)
      INTEGER HPMULT, HPTYPE
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QSELE(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      INTEGER IOBJ(*)
      INTEGER GRPPTR(*), GRPOBJ(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      INTEGER LENG, I, II, REFLCT
      LOGICAL OK, QQSELE
      DOUBLE PRECISION AFCALC, PFCALC, TEMP
      DOUBLE COMPLEX CTEMP
      INTEGER OUNIT, IS, IL, ITEMP
      LOGICAL FOUND
      INTEGER NOBJ, COUNTER, IGROUP, NGROUP
      INTEGER  MSTAR, NSTAR
      INTEGER  STARLBL, STAROBJ, STARVAL
      INTEGER  STARTRN, STARTRPI, STARTRPS
      INTEGER  STARBUF, LSTARBUF
C parameters
      INTEGER MAXLEN
      PARAMETER (MAXLEN=80)
      DOUBLE PRECISION RAD, ONE
      PARAMETER (RAD=PI/180.0D0, ONE=1.0D0)
C
      CHARACTER*(MAXLEN) LINE
C
C begin
      NSTAR = 0
      MSTAR = 0
      STARLBL  = 0
      STAROBJ  = 0
      STARVAL  = 0
      STARTRN  = 0
      STARTRPI = 0
      STARTRPS = 0
C
C defaults
      OFILE='OUTPUT'
      NOBJ=0
      QQSELE=.FALSE.
C
C parsing
      CALL PUSEND('WRITE>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITE>')
      CALL MISCOM('WRITE>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-write-reflection')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,QSELE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
      QQSELE=.TRUE.
C
      ELSE IF (WD(1:4).EQ.'STAR') THEN
        IF (NSTAR .EQ. MSTAR) THEN
          MSTAR = MSTAR + MAX(10, XSFNUM)
          STARLBL  = CREAHP(STARLBL, NSTAR, MSTAR, WORD_SIZE)
          STAROBJ  =  REAHP(STAROBJ,  INTEG4(NSTAR), INTEG4(MSTAR))
          STARVAL  = CREAHP(STARVAL, NSTAR, MSTAR, 4)
          STARTRN  =  REAHP(STARTRN,  INTEG4(NSTAR), INTEG4(MSTAR))
          STARTRPI =  REAHP(STARTRPI, INTEG4(NSTAR), INTEG4(MSTAR))
          STARTRPS =  REAHP(STARTRPS, INTEG4(NSTAR), INTEG4(MSTAR))
        END IF
        CALL XWRSTARP(CWSHEAP(STARLBL), HEAP(STAROBJ), C4HEAP(STARVAL),
     &                HEAP(STARTRN), HEAP(STARTRPI), HEAP(STARTRPS),
     &                NSTAR,
     &                XSFNUM, XSFNAM, XSFTYPE, HPSF)
C
      ELSE
C
C check reciprocal space objects
      OK=.FALSE.
      DO I=1,XSFNUM
      IF (WD(1:WDLEN).EQ.XSFNAM(I)) THEN
      OK=.TRUE.
      II=I
      END IF
      END DO
C
      IF (OK) THEN
      IF (HPSF(II).EQ.0) THEN
      WRITE(6,'(3A)') ' %XRWRIT-ERR: reciprocal space object ',
     &  WD(1:WDLEN),' undefined.'
      CALL WRNDIE(-5,'XRWRIT',' object undefined.')
      ELSEIF (NOBJ.GE.XSFNUM) THEN
      CALL WRNDIE(-5,'XRWRIT',
     &  'Exceeded number of objects to be written.')
      ELSE
      NOBJ=NOBJ+1
      IOBJ(NOBJ)=II
      END IF
      END IF
C
      IF (.NOT.OK) THEN
      CALL CHKEND('WRITE>',DONE)
      END IF
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.QQSELE) THEN
C make default selection ( all )
      DO REFLCT=1,XRNREF
      QSELE(REFLCT)=.TRUE.
      END DO
      END IF
C
      IF (NSTAR .GT. 0) THEN
        IF (NOBJ .GT. 0) THEN
          CALL WRNDIE(-5, 'XRWRIT',
     &                'STAR and non-STAR output cannot be mixed.')
        ELSE
          CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
          IF (.NOT.ERROR) THEN
            STARBUF = CALLHP(NSTAR, 80)
            LSTARBUF = ALLHP(INTEG4(NSTAR))
            CALL XWRSTAR(OUNIT, C80HEAP(STARBUF), HEAP(LSTARBUF),
     &                   CWSHEAP(STARLBL), HEAP(STAROBJ),
     &                   C4HEAP(STARVAL),
     &                   HEAP(STARTRN), HEAP(STARTRPI), HEAP(STARTRPS),
     &                   NSTAR,
     &                   XSFTYPE, HPSF,
     &                   XRNREF, XRH, XRK, XRL, QSELE)
            CALL CFREHP(STARBUF, NSTAR, 80)
            CALL FREHP(LSTARBUF, INTEG4(NSTAR))
            CALL VCLOSE(OUNIT,'KEEP',ERROR)
          END IF
        END IF
C
      ELSE
C if no objects are explicitly specified then write
C all objects
      IF (NOBJ.EQ.0) THEN
C
C write all declared and allocated user objects
      DO I=1,XSFNUM
      IF (HPSF(I).NE.0) THEN
CCC      QTOUCH=.FALSE.
C
C check if reciprocal space object is not equal to zero
CCC      CALL XQTOUCH(QTOUCH,XRNREF,XSFNUM,HPSF,
CCC     &                   XSFTYPE,I)
CCC      IF (QTOUCH) THEN
      NOBJ=NOBJ+1
      IOBJ(NOBJ)=I
CCC      END IF
      END IF
      END DO
C
      END IF
C
      IF (NOBJ.GT.0) THEN
C
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
C
C write the number of reflection written to file
      II=0
      DO REFLCT=1,XRNREF
      IF (QSELE(REFLCT)) II=II+1
      END DO
      WRITE(OUNIT,'(A,I10)') ' NREFlection=',II
C
      IF (QHERM) THEN
      WRITE(OUNIT,'(A)') ' ANOMalous=FALSe { equiv. to HERMitian=TRUE}'
      ELSE
      WRITE(OUNIT,'(A)') ' ANOMalous=TRUE { equiv. to HERMitian=FALSe}'
      END IF
C
      DO I=1,NOBJ
      WRITE(OUNIT,'(5A)') ' DECLare NAME=',XSFNAM(IOBJ(I)),
     &   '   DOMAin=RECIprocal   TYPE=',XSFTYPE(IOBJ(I)),' END'
      END DO
C
C get object group list
      CALL XGLIST(XSFNUM,XSFNAM,XSFTYPE,XSFGNAM,
     &                  XSFGTYP,XSFGORD,NGROUP,GRPOBJ,
     &                  GRPPTR)
C
      DO IGROUP=1,NGROUP
C
C check if this group is included in the output list
      FOUND=.FALSE.
      DO COUNTER=GRPPTR(IGROUP)+1,GRPPTR(IGROUP+1)
      DO I=1,NOBJ
      IF (IOBJ(I).EQ.GRPOBJ(COUNTER)) THEN
      FOUND=.TRUE.
      END IF
      END DO
      END DO
C
      IF (FOUND) THEN
      WRITE(OUNIT,'(3A)') ' GROUp ',
     &      ' TYPE=',XSFGTYP(GRPOBJ(GRPPTR(IGROUP)+1))
      DO COUNTER=GRPPTR(IGROUP)+1,GRPPTR(IGROUP+1)
C
C check if this group is included in the output list
      FOUND=.FALSE.
      DO I=1,NOBJ
      IF (IOBJ(I).EQ.GRPOBJ(COUNTER)) THEN
      FOUND=.TRUE.
      END IF
      END DO
      IF (.NOT.FOUND) THEN
      WRITE(6,'(2A)')
     &  ' %XWRITE-ERR: must write whole group for object ',
     &    XSFNAM(GRPOBJ(COUNTER))
      CALL WRNDIE(-5,'XWRITE','specify all objects of group.')
C
      ELSE
      WRITE(OUNIT,'(2A)') '     OBJEct=',XSFNAM(GRPOBJ(COUNTER))
      END IF
C
      END DO
      WRITE(OUNIT,'(A)') ' END'
      END IF
C
      END DO
C
C
      DO REFLCT=1,XRNREF
C
C
      IF (QSELE(REFLCT)) THEN
      LINE(1:6)=' INDE'
      WRITE(LINE(7:11),'(I5)') XRH(REFLCT)
      WRITE(LINE(12:16),'(I5)') XRK(REFLCT)
      WRITE(LINE(17:21),'(I5)') XRL(REFLCT)
      IS=22
C
      DO I=1,NOBJ
C
      IF (XSFTYPE(IOBJ(I)).EQ.'COMP') THEN
      CALL XCOPY(CTEMP,1,HEAP(HPSF(IOBJ(I))),REFLCT)
      CALL XPHASE(CTEMP,AFCALC,PFCALC)
      LENG=LEN(XSFNAM(IOBJ(I)))
      CALL TRIMM(XSFNAM(IOBJ(I)),LENG)
      PFCALC=PFCALC/RAD
      IF (ABS(AFCALC).GT.1.0D+5.OR.ABS(PFCALC).GT.1.0D+5) THEN
      CALL XWRDMP(IS,33+LENG,IL,LINE,MAXLEN,OUNIT)
      WRITE(LINE(IS:IL),'(3A,E16.6,E16.6)')
     &           ' ',XSFNAM(IOBJ(I))(1:LENG),'=',AFCALC,PFCALC
      ELSE
      CALL XWRDMP(IS,21+LENG,IL,LINE,MAXLEN,OUNIT)
      WRITE(LINE(IS:IL),'(3A,F10.3,F10.3)')
     &           ' ',XSFNAM(IOBJ(I))(1:LENG),'=',AFCALC,PFCALC
      END IF
      IS=IL+1
C
      ELSEIF (XSFTYPE(IOBJ(I)).EQ.'REAL') THEN
      CALL XCOPYR(TEMP,1,HEAP(HPSF(IOBJ(I))),REFLCT)
      LENG=LEN(XSFNAM(IOBJ(I)))
      CALL TRIMM(XSFNAM(IOBJ(I)),LENG)
      IF (ABS(TEMP).GT.1.0D+5) THEN
      CALL XWRDMP(IS,17+LENG,IL,LINE,MAXLEN,OUNIT)
      WRITE(LINE(IS:IL),'(3A,E16.6)')
     &           ' ',XSFNAM(IOBJ(I))(1:LENG),'=',TEMP
      ELSE
      CALL XWRDMP(IS,11+LENG,IL,LINE,MAXLEN,OUNIT)
      WRITE(LINE(IS:IL),'(3A,F10.3)')
     &           ' ',XSFNAM(IOBJ(I))(1:LENG),'=',TEMP
      END IF
C
      IS=IL+1
C
      ELSEIF (XSFTYPE(IOBJ(I)).EQ.'INTE') THEN
      CALL XCOPYI(ITEMP,1,HEAP(HPSF(IOBJ(I))),REFLCT)
      LENG=LEN(XSFNAM(IOBJ(I)))
      CALL TRIMM(XSFNAM(IOBJ(I)),LENG)
      CALL XWRDMP(IS,11+LENG,IL,LINE,MAXLEN,OUNIT)
      WRITE(LINE(IS:IL),'(3A,I10)')
     &           ' ',XSFNAM(IOBJ(I))(1:LENG),'=',ITEMP
      IS=IL+1
      END IF
C
      END DO
C
      WRITE(OUNIT,'(A)') LINE(1:IS-1)
C
      END IF
      END DO
      CALL VCLOSE(OUNIT,'KEEP',ERROR)
      END IF
C
      END IF
C
      END IF
C
      CALL FREESTARTR(HEAP(STARTRN),
     &                HEAP(STARTRPI), HEAP(STARTRPS), NSTAR)
C
      IF (MSTAR .GT. 0) THEN
        CALL CFREHP(STARLBL, MSTAR, WORD_SIZE)
        CALL  FREHP(STAROBJ,  INTEG4(MSTAR))
        CALL CFREHP(STARVAL, MSTAR, 4)
        CALL  FREHP(STARTRN,  INTEG4(MSTAR))
        CALL  FREHP(STARTRPI, INTEG4(MSTAR))
        CALL  FREHP(STARTRPS, INTEG4(MSTAR))
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XWRDMP(IS,LENGHT,IL,LINE,MAXLEN,OUNIT)
C
C Routine dumps line if it exceeds the maximum length
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER IS,LENGHT,IL
      CHARACTER*(*) LINE
      INTEGER MAXLEN, OUNIT
C begin
      IL=IS+LENGHT
      IF (IL.GT.MAXLEN) THEN
      WRITE(OUNIT,'(A)') LINE(1:IS-1)
      IL=IL-IS+19
      IS=19
      LINE(1:18)=' '
      END IF
      RETURN
      END
C=======================================================================
      SUBROUTINE XREFDEL(XRMREF,XRNREF,HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,
     &           HPSF,HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C
C Routine deletes selected reflections
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XRMREF, XRNREF
      INTEGER HPTSEL, HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      INTEGER NQSELE
C pointer
      INTEGER QSELE, MAP2
C begin
C
      NQSELE=XRNREF
      QSELE=ALLHP(ILOGIC(NQSELE))
      MAP2=ALLHP(INTEG4(NQSELE))
      CALL XREFDE2(XRMREF,XRNREF,
     &           HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,
     &           HPSF,HPMULT,HPTYPE,
     &           XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,HEAP(QSELE),HEAP(MAP2),
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
      CALL FREHP(MAP2,INTEG4(NQSELE))
      CALL FREHP(QSELE,ILOGIC(NQSELE))
      RETURN
      END
C======================================================================
      SUBROUTINE XREFDE2(XRMREF,XRNREF,
     &           HPTSEL,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,
     &           HPSF,HPMULT,HPTYPE,
     &           XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QSELE,MAP2,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C
C
C See routine XREFDEF above
C =========================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XRMREF, XRNREF
      INTEGER HPTSEL, HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT, HPTYPE
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QSELE(*)
      INTEGER MAP2(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      LOGICAL QQSELE
      INTEGER N, I, REFLCT
C pointer
      INTEGER IPTR
C begin
C
C defaults
      QQSELE=.FALSE.
C
C
C parsing
      CALL PUSEND('DELEte>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DELEte>')
      CALL MISCOM('DELEte>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-delete')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,QSELE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
      QQSELE=.TRUE.
C
C
C============================================================
      ELSE
      CALL CHKEND('DELEte>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (QQSELE.AND.XRNREF.GT.0) THEN
C
      N=0
      DO REFLCT=1,XRNREF
      IF (.NOT.QSELE(REFLCT)) THEN
      N=N+1
      MAP2(N)=REFLCT
      END IF
      END DO
C
C now we throw out multiple entries by mapping all data
      IPTR=ALLHP(INTEG4(N))
      CALL AINDX4(MAP2,HEAP(HPH),N,HEAP(IPTR))
      CALL AINDX4(MAP2,HEAP(HPK),N,HEAP(IPTR))
      CALL AINDX4(MAP2,HEAP(HPL),N,HEAP(IPTR))
      CALL AINDX4(MAP2,HEAP(HPMULT),N,HEAP(IPTR))
      CALL AINDX4(MAP2,HEAP(HPTYPE),N,HEAP(IPTR))
      CALL AINDX4(MAP2,HEAP(HPTSEL),N,HEAP(IPTR))
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
      WRITE(6,'(A,I8,A)') ' XREFDE2: ',XRNREF-N,
     &  ' reflections have been deleted.'
      XRNREF=N
C
      END IF
C
      RETURN
      END
C
