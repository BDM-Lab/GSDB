C======================================================================
      SUBROUTINE XSCALE(XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,
     &           HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRCELL,XRVOL)
C
C Scaling routine for structure factors.
C
C Authors: Axel T. Brunger and J.-S. Jiang
C ========================================
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
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*), XRTR(3,3)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
C pointer
      INTEGER QSELE, INDX
C begin
C allocate space for the selection
C
      QSELE=ALLHP(ILOGIC(XRNREF))
      INDX=ALLHP(INTEG4(XRNREF))
      CALL XSCALE2(XRMREF,XRNREF,HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,HEAP(QSELE),HEAP(INDX),
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
      CALL FREHP(INDX,INTEG4(XRNREF))
      CALL FREHP(QSELE,ILOGIC(XRNREF))
      RETURN
      END
C======================================================================
      SUBROUTINE XSCALE2(XRMREF,XRNREF,
     &           HPH,HPK,HPL,
     &           XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           HPMULT,HPTYPE,
     &           XRTR,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &           XRSYMM,XRITSY,QSELE,INDX,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,XRCELL,XRVOL)
C
C See routine XSCALE above
C
C Authors: Axel T. Brunger and J.-S. Jiang
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      DOUBLE PRECISION XRTR(3,3)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QSELE(*)
      INTEGER INDX(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION XRCELL(6), XRVOL
C local
      DOUBLE PRECISION FFK
      LOGICAL XRFFKQ
      INTEGER REFLCT, NSET, NNSELE, J
      DOUBLE PRECISION KSCAL(20), BSCAL(20), TEMP
      CHARACTER*6 SET(20), CTEMP, PROMPT
      DOUBLE COMPLEX DUCOMP
      LOGICAL ERR, OK, COND
      INTEGER ISET, MAXSET, TGSET
      INTEGER HPSFSET(20), NNSET(20), TYPESF(20), ISFN(20)
      CHARACTER*4 MODE, UNIF
      LOGICAL QUPDA
      INTEGER NCYC, DIAG, FIXMOD
      DOUBLE PRECISION EPS, KSMIN, BFMIN, BFMAX, KRES
      DOUBLE PRECISION KINI, BINI
      DOUBLE PRECISION FFK1, FFK2
      LOGICAL QFFK1, QFFK2, QTLOW
      DOUBLE PRECISION RVAL, TARG
      INTEGER SILENT
      LOGICAL QANISO, QISO
      DOUBLE PRECISION BTENS(120)
      CHARACTER*4 RESTRC
C incorporate with VSTACK
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N, VSTACK, XEDNI
C parameters
      DOUBLE PRECISION MARK, ONE, ZERO, BLARGE
      PARAMETER (MARK=-9999.D0, ONE=1.0D0, ZERO=0.0D0, BLARGE=500.0D0)
C begin
      FFK=ZERO
      XRFFKQ=.TRUE.
C
C defaults for VSTACK - not use always
      QVSTACK=.FALSE.
      VLEVEL=1
      VMAX=2
      N=2
      XEDNI=ALLHP(INTEG4(N))
      VSTACK=ALLHP(ICPLX8(VMAX*N))
C
C make all defaults
      DO REFLCT=1,XRNREF
      QSELE(REFLCT)=.TRUE.
      INDX(REFLCT)=REFLCT
      END DO
      FFK1=ZERO
      FFK2=ZERO
      QFFK1=.TRUE.
      QFFK2=.TRUE.
      QTLOW=.FALSE.
C
C modified, ATB, 12/02/08
      SILENT=5
      QANISO=.FALSE.
      QISO=.TRUE.
      RESTRC='ALL '
C
      MAXSET=0
      DO J=1,6
      BTENS(J)=ZERO
      END DO
      DO J=7,120
      BTENS(J)=MARK
      END DO
      DO NSET=1,20
      KSCAL(NSET)=MARK
      BSCAL(NSET)=MARK
      SET(NSET)=' '
      HPSFSET(NSET)=0
      NNSET(NSET)=0
      TYPESF(NSET)=0
      ISFN(NSET)=0
      END DO
      KSCAL(1)=-ONE
      BSCAL(1)=ZERO
      MODE='TARG'
      QUPDA=.TRUE.
      NCYC=30
      DIAG=0
      EPS=0.001D0
      KSMIN=ZERO
      BFMIN=-BLARGE
      BFMAX=BLARGE
      FIXMOD=0
      UNIF=' '
      KINI=ONE
      BINI=ZERO
      KRES=ZERO
C
      ERR=.FALSE.
      ISET=0
C
C parsing
      CALL PUSEND('MULTiscale>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MULTiscale>')
      CALL MISCOM('MULTiscale>',USED)
      IF (.NOT.USED) THEN
C======================================================================
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-multiscale')
C
C======================================================================
      ELSE IF (WD(1:1).EQ.'?') THEN
      WRITE(6,'(A,A)') ' |---------------------------',
     & ' MULTiscale parameters ---------------------------'
      WRITE(6,'(A,I4,A,I3,A,F9.4)')
     & ' | least-squares: NCYC=',NCYC,'   DIAG=',DIAG,'   EPS=',EPS
      WRITE(6,'(A,F8.3,A,F8.3)')
     & ' | initials: KINI=',KINI,'   BINI=',BINI
      WRITE(6,'(A,F8.3,A,F8.3,A,F8.3)')
     & ' | restraints: KSMIn=',KSMIN,
     & '   BMIN=',BFMIN,'   BMAX=',BFMAX
      IF (UNIF.EQ.'K') WRITE(6,'(A)')
     & ' | uniform scale factor applied'
      IF (UNIF.EQ.'B') WRITE(6,'(A)')
     & ' | uniform B factor applied'
      IF (.NOT.XRFFKQ) WRITE(6,'(A,F8.3)')
     & ' | the overall scale factor is fixed FFK=',FFK
      IF (KRES.GT.RSMALL) THEN
      WRITE(6,'(A,F8.3)')
     & ' | two resolution-dependent scales are applied.  RESK=',
     & SQRT(ONE/KRES)
      IF (.NOT.QFFK1) WRITE(6,'(A,F8.3)')
     & ' | overall scale for low resolution is fixed FFK1=',FFK1
      IF (.NOT.QFFK2) WRITE(6,'(A,F8.3)')
     & ' | overall scale for high resolution is fixed FFK2=',FFK2
      END IF
      IF (QANISO) THEN
      WRITE(6,'(A,A)') ' | RESTriction= ',RESTRC
      WRITE(6,'(A,A)') ' | MULTiscale:  <i>    SET<i>',
     & '      K<i>          B<i>_tensor'
      ELSE
      WRITE(6,'(A,A)') ' | MULTiscale:  <i>    SET<i>',
     & '      K<i>          B<i>'
      END IF
      DO J=1,ISET
      IF (HPSFSET(J).NE.0) THEN
      IF (QANISO) THEN
      WRITE(6,'(A,I2,3A,F10.4/6(A,F8.2))')
     &  ' |            ',NNSET(J),
     &  '      ',SET(NNSET(J)),'   ', KSCAL(NNSET(J)),
     &  ' | B11=',BTENS(NNSET(J)*6-5),
     &  ' B22=',BTENS(NNSET(J)*6-4),
     &  ' B33=',BTENS(NNSET(J)*6-3),
     &  ' B12=',BTENS(NNSET(J)*6-2),
     &  ' B13=',BTENS(NNSET(J)*6-1),
     &  ' B23=',BTENS(NNSET(J)*6)
      ELSE
      WRITE(6,'(A,I2,3A,F10.4,A,F10.4)') ' |            ',NNSET(J),
     &   '      ',SET(NNSET(J)),'   ',
     &   KSCAL(NNSET(J)),'   ',BSCAL(NNSET(J))
      END IF
      END IF
      END DO
      WRITE(6,'(A,A)') ' |-------------------------------',
     & '-----------------------------------------------'
C======================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL XFSELE(XRTR,XRMREF,XRNREF,HPH,HPK,HPL,
     &    XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &    HPMULT,HPTYPE,
     &    QHERM,XRNSYM,XRMSYM,XRSYTH,
     &    XRSYMM,XRITSY,QSELE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &    XRCELL,XRVOL)
C======================================================================
      ELSE IF (WD(1:3).EQ.'SET') THEN
      OK=WDLEN.GE.4
      IF (OK) NSET=DECODI(WD(4:WDLEN),WDLEN-3,OK)
      IF (.NOT.OK.OR.NSET.LE.0.OR.NSET.GE.20) THEN
      ERR=.TRUE.
      WRITE(6,'(A)') ' %SCALE-ERR: incorrect definition of SET<i>'
      END IF
      PROMPT=WD(1:WDLEN)//'='
      CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD(PROMPT)
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') PROMPT,SET(NSET)
      ELSE
      CALL XSCALSET(NSET,SET,WD,WDLEN,ISET,NNSET,HPSFSET,TYPESF,
     &              ISFN,COND,XSFNUM,XSFNAM,XSFTYPE,HPSF)
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XSCALE-ERR: reciprocal space object ',
     &                 WD(1:WDLEN),' is undeclared.'
      CALL WRNDIE(-5,'XSCALE','undeclared object.')
      END IF
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'ANIS') THEN
      CALL NEXTLO('ANISotropic=',QANISO)
      ELSE IF (WD(1:WDLEN).EQ.'FFK') THEN
      CALL NEXTF('FFK=',FFK)
      IF (FFK.GT.ZERO) THEN
      XRFFKQ=.FALSE.
      ELSE
      XRFFKQ=.TRUE.
      END IF
      ELSE IF (WD(1:4).EQ.'ISOT') THEN
      CALL NEXTLO('ISOTropic=',QISO)
      ELSE IF (WD(1:4).EQ.'REST') THEN
      CALL NEXTA4('RESTriction=',RESTRC)
C======================================================================
      ELSE IF (WD(1:1).EQ.'K'.AND.(WDLEN.EQ.2.OR.WDLEN.EQ.3)) THEN
      OK=WDLEN.EQ.2.OR.WDLEN.EQ.3
      IF (OK) NSET=DECODI(WD(2:WDLEN),WDLEN-1,OK)
      IF (OK.AND.(NSET.GT.0.OR.NSET.LE.20)) THEN
      CALL NEXTF(WD(1:WDLEN)//'=',KSCAL(NSET))
      END IF
C======================================================================
      ELSE IF (WD(1:1).EQ.'B'.AND.(WDLEN.EQ.2.OR.WDLEN.EQ.3)) THEN
      IF (OK) NSET=DECODI(WD(2:WDLEN),WDLEN-1,OK)
      IF (OK.AND.(NSET.GT.0.OR.NSET.LE.20)) THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BSCAL(NSET))
      END IF
C======================================================================
      ELSE IF (WD(1:1).EQ.'B'.AND.(INDEX(WD(1:WDLEN),'_').NE.0)) THEN
      J=INDEX(WD(1:WDLEN),'_')
      OK=J.GT.2.AND.QANISO
      IF (OK) NSET=DECODI(WD(2:J-1),J-2,OK)
      IF (OK.AND.(NSET.GT.0.OR.NSET.LE.20)) THEN
      IF (WD(J+1:WDLEN).EQ.'11') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6-5))
      ELSE IF (WD(J+1:WDLEN).EQ.'22') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6-4))
      ELSE IF (WD(J+1:WDLEN).EQ.'33') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6-3))
      ELSE IF (WD(J+1:WDLEN).EQ.'12') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6-2))
      ELSE IF (WD(J+1:WDLEN).EQ.'13') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6-1))
      ELSE IF (WD(J+1:WDLEN).EQ.'23') THEN
      CALL NEXTF(WD(1:WDLEN)//'=',BTENS(NSET*6))
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('scale-mode=',MODE)
      IF (MODE.EQ.'TLOW') QTLOW=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UPDA') THEN
      CALL NEXTLO('update-parameters=',QUPDA)
C======================================================================
      ELSE IF (WD(1:4).EQ.'NCYC') THEN
      CALL NEXTI('cycle limits=',NCYC)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DIAG') THEN
      CALL NEXTI('diagonal-approximation=',DIAG)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'EPS ') THEN
      CALL NEXTF('convergence-terminitor=',EPS)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'KSMI') THEN
      CALL NEXTF('minimum-scale-factor=',KSMIN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BFMI'.OR.WD(1:4).EQ.'BMIN') THEN
      CALL NEXTF('minimum-B-factor=',BFMIN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BFMA'.OR.WD(1:4).EQ.'BMAX') THEN
      CALL NEXTF('maximum-B-factor=',BFMAX)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FIXM') THEN
      CALL NEXTI('one-uniform-parameter=',FIXMOD)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'UNIF') THEN
      CALL NEXTA4('one-uniform-parameter=',UNIF)
      IF (UNIF(1:1).EQ.'K') THEN
      FIXMOD=3
      ELSE IF (UNIF(1:1).EQ.'B') THEN
      FIXMOD=4
      ELSE
      FIXMOD=0
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'KINI') THEN
      CALL NEXTF('initial-scale-factor=',KINI)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'BINI') THEN
      CALL NEXTF('initial-B-factor=',BINI)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RESK') THEN
      CALL NEXTF('scaling-resolution-boundary=',KRES)
      IF (KRES.GT.RSMALL) THEN
      KRES=ONE/(KRES*KRES)
      ELSE
      KRES=ZERO
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FFK1') THEN
      CALL NEXTF('FFK1=',FFK1)
      IF (ABS(FFK1).LT.RSMALL) THEN
      QFFK1=.TRUE.
      ELSE
      QFFK1=.FALSE.
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FFK2') THEN
      CALL NEXTF('FFK2=',FFK2)
      IF (ABS(FFK2).LT.RSMALL) THEN
CCC modification ATB 4/27/08
      QFFK2=.TRUE.
      ELSE
      QFFK2=.FALSE.
      END IF
C=====================================================================
      ELSE
      CALL CHKEND('MULTiscale>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C finish parsing
C
      IF (.NOT.QISO.AND..NOT.QANISO) THEN
      WRITE(6,'(A)')
     & ' %XSCALE-ERR: Inconsistent options: QISO=FALSE and ',
     & '              QANISO=FALSE.        QISO set to TRUE.'
      QISO=.TRUE.
      END IF
C
C
C {{{{{{{{{{{{{{{{{{{{{{{{{{{{{
C total number of data sets were selected
      MAXSET=ISET
      IF (MAXSET.GT.1) THEN
C
C modification, 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,A,I4)') ' XSCALE:',
     & ' total number of data sets were selected ',MAXSET
      END IF
      ELSE
      WRITE(6,'(A,A,I4)') ' %XSCALE-err:',
     & ' at least two data sets are required ',MAXSET
      ERR=.TRUE.
      END IF
C which one is the target data set ('observed')? the first K<i>=-1
      TGSET=0
      DO ISET=MAXSET,1,-1
      IF (KSCAL(NNSET(ISET)).EQ.-ONE) TGSET=ISET
      END DO
      IF (TGSET.EQ.0) THEN
      WRITE(6,'(A)') ' %XSCALE-ERR: no target data set is specified'
      ERR=.TRUE.
      ELSE
C
C modification, 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN      
      WRITE(6,'(A,I2,A,A,A)') ' XSCALE: set',NNSET(TGSET),' ',
     & SET(NNSET(TGSET)),' is specified as the target data set'
      END IF
      END IF
C }}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
C
C check data selection
      NNSELE=0
      DO REFLCT=1,XRNREF
      IF (QSELE(REFLCT)) THEN
      NNSELE=NNSELE+1
      INDX(NNSELE)=REFLCT
      END IF
      END DO
      DO REFLCT=NNSELE+1,XRNREF
      INDX(REFLCT)=0
      END DO
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A,I9,A)')
     & ' Total of ',NNSELE,' structure factor elements were selected.'
      END IF
      TEMP=NNSELE
      CALL DECLAR('SELECT','DP',' ',DUCOMP,TEMP)
C
C
      IF (.NOT.ERR) THEN
C
      IF (QANISO) THEN
      CALL XSCALIT(KSCAL,BTENS,NNSELE,INDX,MARK,
     &             HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &             MAXSET,TGSET,NNSET,HPSFSET,TYPESF,SET,
     &             MODE,NCYC,DIAG,EPS,KSMIN,BFMIN,BFMAX,
     &             FIXMOD,KRES,KINI,BINI,QUPDA,FFK,RVAL,TARG,
     &             XRFFKQ,FFK1,FFK2,QFFK1,QFFK2,QTLOW,
     &             QVSTACK,VLEVEL,VMAX,HEAP(VSTACK),N,HEAP(XEDNI),
     &             SILENT,QANISO,QISO,RESTRC,XRTR,XRCELL,XRVOL)
      ELSE
      CALL XSCALIT(KSCAL,BSCAL,NNSELE,INDX,MARK,
     &             HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &             MAXSET,TGSET,NNSET,HPSFSET,TYPESF,SET,
     &             MODE,NCYC,DIAG,EPS,KSMIN,BFMIN,BFMAX,
     &             FIXMOD,KRES,KINI,BINI,QUPDA,FFK,RVAL,TARG,
     &             XRFFKQ,FFK1,FFK2,QFFK1,QFFK2,QTLOW,
     &             QVSTACK,VLEVEL,VMAX,HEAP(VSTACK),N,HEAP(XEDNI),
     &             SILENT,QANISO,QISO,RESTRC,XRTR,XRCELL,XRVOL)
      END IF
C
C print information and declare symbols.
C
C modification 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      IF (QANISO) THEN
      WRITE(6,'(A,A)') ' MULTiscale:  <i>    SET<i>',
     & '      K<i>          B<i>_tensor'
      ELSE
      WRITE(6,'(A,A)') ' MULTiscale:  <i>    SET<i>',
     & '      K<i>          B<i>'
      END IF
      END IF
      DO ISET=1,MAXSET
      NSET=NNSET(ISET)
      CALL ENCODI(NSET,WD,WDMAX,WDLEN)
C
C modification, 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      IF (QANISO) THEN
      WRITE(6,'(A,I2,3A,F10.4/6(A,F8.3))')
     &  '              ',NSET,
     &  '      ',SET(NSET),'   ', KSCAL(NSET),
     &  ' B11=',BTENS(NSET*6-5),
     &  ' B22=',BTENS(NSET*6-4),
     &  ' B33=',BTENS(NSET*6-3),
     &  ' B12=',BTENS(NSET*6-2),
     &  ' B13=',BTENS(NSET*6-1),
     &  ' B23=',BTENS(NSET*6)
      ELSE
      WRITE(6,'(A,I2,3A,F10.4,A,F10.4)') '             ',NSET,
     &   '      ',SET(NSET),'   ',
     &   KSCAL(NSET),'   ',BSCAL(NSET)
      END IF
      END IF
C
      TEMP=KSCAL(NSET)
      CTEMP='K'//WD(1:WDLEN)
      CALL DECLAR(CTEMP(1:WDLEN+1),'DP',' ',DUCOMP,TEMP)
      IF (QANISO) THEN
      TEMP=BTENS(NSET*6-5)
      CTEMP='B'//WD(1:WDLEN)//'_11'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      TEMP=BTENS(NSET*6-4)
      CTEMP='B'//WD(1:WDLEN)//'_22'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      TEMP=BTENS(NSET*6-3)
      CTEMP='B'//WD(1:WDLEN)//'_33'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      TEMP=BTENS(NSET*6-2)
      CTEMP='B'//WD(1:WDLEN)//'_12'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      TEMP=BTENS(NSET*6-1)
      CTEMP='B'//WD(1:WDLEN)//'_13'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      TEMP=BTENS(NSET*6)
      CTEMP='B'//WD(1:WDLEN)//'_23'
      CALL DECLAR(CTEMP(1:WDLEN+4),'DP',' ',DUCOMP,TEMP)
      ELSE
      TEMP=BSCAL(NSET)
      CTEMP='B'//WD(1:WDLEN)
      CALL DECLAR(CTEMP(1:WDLEN+1),'DP',' ',DUCOMP,TEMP)
      END IF
      END DO
C
C modification 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      IF (QANISO) THEN
      WRITE(6,'(A/A)')
     & ' MULTiscale: K<i> is stored in symbols $K<i>, and B<i>_tensor',
     & ' MULTiscale: are stored in $B<i>_11, $B<i>_22, $B<i>_33, etc.'
      ELSE
      WRITE(6,'(A)')
     & ' MULTiscale: K<i>, B<i> are stored in symbols $K<i>, $B<i>'
      END IF
      END IF
C
      CALL DECLAR( 'XSCFFK', 'DP', ' ', DUCOMP, FFK )
C
C modification 12/02/08, ATB
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     &  ' MULTiscale: the overall scale is stored in symbol $XSCFFK'
      END IF
      IF (KRES.GT.RSMALL) THEN
      CALL DECLAR( 'XSCFFK1', 'DP', ' ', DUCOMP, FFK1 )
      IF (WRNLEV.GT.5) WRITE(6,'(A,F8.2,A)')
     &  ' MULTiscale: $XSCFFK1 is for low resolution d > ',
     &  SQRT(ONE/KRES),'(A)'
      CALL DECLAR( 'XSCFFK2', 'DP', ' ', DUCOMP, FFK2 )
      IF (WRNLEV.GT.5) WRITE(6,'(A,F8.2,A)')
     &  ' MULTiscale: $XSCFFK2 is for high resolution d <= ',
     &  SQRT(ONE/KRES),'(A)'
      END IF
C
      ELSE
      CALL WRNDIE(-5,'MULTiscale','There were some errors')
      END IF
C
C free space for VSTACK
      CALL FREHP(VSTACK,ICPLX8(VMAX*N))
      CALL FREHP(XEDNI,INTEG4(N))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCALSET(NSET,SET,WD,WDLEN,ISET,NNSET,HPSFSET,
     &                TYPESF,ISFN,COND,XSFNUM,XSFNAM,XSFTYPE,HPSF)
C
C Get the operand name, heap and type.
C Type=1 complex, type=2 real and type=3 integer.
C
C Authors: Axel T. Brunger and J.-S. Jiang
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER NSET
      CHARACTER*(*) SET(*)
      CHARACTER*(*) WD
      INTEGER WDLEN, ISET
      INTEGER NNSET(*), HPSFSET(*), TYPESF(*), ISFN(*)
      LOGICAL COND
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
C local
      INTEGER I
C begin
C
      SET(NSET)=WD(1:WDLEN)
C
C check reciprocal space objects
      COND=.FALSE.
      DO I=1,XSFNUM
      IF (WD(1:WDLEN).EQ.XSFNAM(I)) THEN
      COND=.TRUE.
      IF (HPSF(I).EQ.0) THEN
      WRITE(6,'(3A)') ' %XSCALE-ERR: reciprocal space object ',
     & WD(1:WDLEN),' undefined.'
      CALL WRNDIE(-5,'XSCALE','object undefined.')
      ELSE
C
      ISET=ISET+1
      HPSFSET(ISET)=HPSF(I)
      ISFN(ISET)=I
      NNSET(ISET)=NSET
C
      IF (XSFTYPE(I).EQ.'COMP') THEN
      TYPESF(ISET)=1
      ELSEIF (XSFTYPE(I).EQ.'REAL'.OR.XSFTYPE(I).EQ.'PHAS') THEN
      TYPESF(ISET)=2
      ELSEIF (XSFTYPE(I).EQ.'INTE') THEN
      TYPESF(ISET)=3
      END IF
C
      END IF
      END IF
      END DO
C
C
      RETURN
      END
C======================================================================
