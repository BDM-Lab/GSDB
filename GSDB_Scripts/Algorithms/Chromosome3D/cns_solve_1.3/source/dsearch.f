C============================================================================
      SUBROUTINE DRSEAR
C
C Direct rotation search
C
C Authors: Axel T. Brunger and Warren DeLano
C ==========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'expression.inc'
C local
      DOUBLE PRECISION T1MIN, T1MAX, T2MIN, T2MAX, T3MIN, T3MAX
      DOUBLE PRECISION DELTA
      INTEGER DUMMY
      DOUBLE PRECISION  EPS
      INTEGER OUNIT, LUNIT, NLIST
C
      INTEGER MAXNR
      DOUBLE PRECISION RAVE, RSIGMA, RMIN, RMAX
C
      CHARACTER*(WORD_SIZE) LFILE
C
C
      CHARACTER*4 ROTMOD
C
      DOUBLE PRECISION SECT, SECS
C
      LOGICAL QDEBUG, QSELF, QFORM
C
      INTEGER N2DIM
C
      CHARACTER*(WDMAX) RFSYMB
      INTEGER RFSYML
C pointer
      INTEGER N1DIM, N3DIM
      INTEGER NR, T1, T2, T3, RFCORR
      INTEGER XTEMP,YTEMP,ZTEMP
C symmetry handling
      INTEGER TSNSYM,XROLD
C parameters
      DOUBLE PRECISION ZERO, TWO, THREE, T1P4
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, THREE=3.0D0, T1P4=1.4D0)
C begin
C
C default values
      TSNSYM=1
      OFILE=' '
      LFILE=' '
      OUNIT=0
      LUNIT=0
      T1MIN=ZERO
      T1MAX=ZERO
      T2MIN=ZERO
      T2MAX=ZERO
      T3MIN=ZERO
      T3MAX=ZERO
      T3MIN=ZERO
      T3MAX=ZERO
      T1MIN=ZERO
      T1MAX=ZERO
      DELTA=ZERO
      ROTMOD='EULE'
      NLIST=100
      EPS=0.2D0
      QDEBUG=.FALSE.
      QSELF=.FALSE.
      QFORM=.TRUE.
      RFSYMB=' '
      RFSYML=1
C
C parsing
      CALL PUSEND('DRSEARCH>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('DRSEARCH>')
      CALL MISCOM('DRSEARCH>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-search-direct')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMB') THEN
      CALL NEXTST('SYMBol=',WD)
      CALL COPYST(RFSYMB,WDMAX,RFSYML,WD,WDLEN)
      IF ( RFSYML.GT.0 ) THEN
         CALL TOUPPER(RFSYMB)
      ELSE
         RFSYMB=' '
         RFSYML=1
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMM') THEN
      XROLD=XRNSYM
      XRNSYM=TSNSYM
      CALL NEXTEX('SYMMetry=',WDD,WDMAX,WDDLEN)
      CALL XRSYPA(WDD,WDDLEN,XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
      TSNSYM=XRNSYM
      XRNSYM=XROLD
C=====================================================================
C dummy parse asymmetry command (is not stored in spacegroup asymmetric unit)
C
      ELSE IF (WD(1:4).EQ.'ASYM') THEN
      CALL XASYMM(RPNMX,RPNN3,RPNX,RPN3,RPNL3,RPNDB3,RPNMLT3,
     &            RPNTYP3,RPNDOM3,RPNLEV3,DUMMY)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'LIST') THEN
      CALL NEXTFI('LIST=',LFILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'EPSI') THEN
      CALL NEXTF('EPSIlon=',EPS)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NLIS') THEN
      CALL NEXTI('NLISt=',NLIST)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DEBU') THEN
      CALL NEXTLO('DEBUg=',QDEBUG)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTLO('FORMatted=',QFORM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T1MI') THEN
      CALL NEXTF('T1MIn= ',T1MIN)
      ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T1MA') THEN
      CALL NEXTF('T1MAx= ',T1MAX)
      ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T2MI') THEN
      CALL NEXTF('T2MIn= ',T2MIN)
      IF (ROTMOD.NE.'LATT') ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T2MA') THEN
      CALL NEXTF('T2MAx= ',T2MAX)
      IF (ROTMOD.NE.'LATT') ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T3MI') THEN
      CALL NEXTF('T3MIn= ',T3MIN)
      ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'T3MA') THEN
      CALL NEXTF('T3MAx= ',T3MAX)
      ROTMOD='EULE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TMMI') THEN
      CALL NEXTF('TMMIn= ',T3MIN)
      ROTMOD='LATT'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TMMA') THEN
      CALL NEXTF('TMMAx= ',T3MAX)
      ROTMOD='LATT'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TPMI') THEN
      CALL NEXTF('TPMIn= ',T1MIN)
      ROTMOD='LATT'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TPMA') THEN
      CALL NEXTF('TPMAx= ',T1MAX)
      ROTMOD='LATT'
C=====================================================================
      ELSE IF (WD(1:5).EQ.'PSIMI') THEN
      CALL NEXTF('PSIMIn= ',T1MIN)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:5).EQ.'PSIMA') THEN
      CALL NEXTF('PSIMAx= ',T1MAX)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:5).EQ.'PHIMI') THEN
      CALL NEXTF('PHIMIn= ',T2MIN)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:5).EQ.'PHIMA') THEN
      CALL NEXTF('PHIMAx= ',T2MAX)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:7).EQ.'KAPPAMI') THEN
      CALL NEXTF('KAPPAMIn= ',T3MIN)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:7).EQ.'KAPPAMA') THEN
      CALL NEXTF('KAPPAMAx= ',T3MAX)
      ROTMOD='SPHE'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DELT') THEN
      CALL NEXTF('DELTa= ',DELTA)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' -----------rotation sear',
     & 'ch-parameters---------------------------------------'
      IF (ROTMOD.EQ.'LATT') THEN
      WRITE(6,'(A)')
     & ' | The quasi-orthogonal Eulerian angle sampling method is used'
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | TPMIn=',T1MIN,' TPMAx=',T1MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | TMMIn=',T3MIN,' TMMAx=',T3MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | T2MIn=',T2MIN,' T2MAx=',T2MAX
      ELSE IF (ROTMOD.EQ.'EULE') THEN
      WRITE(6,'(A)')
     & ' | The equidistant Eulerian angle sampling method is used'
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | T1MIn=',T1MIN,' T1MAx=',T1MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | T2MIn=',T2MIN,' T2MAx=',T2MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | T3MIn=',T3MIN,' T3MAx=',T3MAX
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(6,'(A)')
     & ' | The spherical polar angle sampling method is used'
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | PSIMIn=',T1MIN,' PSIMAx=',T1MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | PHIMIn=',T2MIN,' PHIMAx=',T2MAX
      WRITE(6,'(A,F12.5,A,F12.5)')
     & ' | KAPPAMIn=',T3MIN,' KAPPAMAx=',T3MAX
      END IF
      WRITE(6,'(A,F12.5)')
     & ' | DELTa=',DELTA
      WRITE(6,'(A,I8)')
     & ' | NLISt=',NLIST
      WRITE(6,'(A,F5.2)')
     & ' | EPSIlon=',EPS
      WRITE(6,'(A)')
     & ' | Symmetry operators for cluster analysis'
      XROLD=XRNSYM
      XRNSYM=TSNSYM
      CALL XSYPRI(XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
      XRNSYM=XROLD
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE
      CALL CHKEND('DRSEARCH>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (XRNSYM.NE.1) THEN
      CALL WRNDIE(-5,'DSEARCH','P1 spacegroup must be selected.')
      ELSE
C
C assign output file if required
      IF (OFILE.NE.' ') THEN
      CALL ASSFIL(OFILE,OUNIT,'WRITE','FORMATTED',ERROR)
      IF (ERROR) GOTO 9999
      END IF
C
C assign list output file if required
      IF (LFILE.NE.' ') THEN
      CALL ASSFIL(LFILE,LUNIT,'WRITE','FORMATTED',ERROR)
      IF (ERROR) GOTO 9999
      END IF
C
C
      IF(NATOM.GT.0) THEN
C determine the number of angles to be searched
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL MAKTHE('SCAN',NR,MAXNR,ROTMOD,QSELF,
     & T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     & DELTA,0,0,0,0,N2DIM,0)
      MAXNR=NR
C
C allocate space for the theta's to be searched and
C the rotation function values
      RFCORR=ALLHP(IREAL8(MAXNR))
      T1=ALLHP(IREAL8(MAXNR))
      T2=ALLHP(IREAL8(MAXNR))
      T3=ALLHP(IREAL8(MAXNR))
      N1DIM=ALLHP(INTEG4(N2DIM))
      N3DIM=ALLHP(INTEG4(N2DIM))
      XTEMP=ALLHP(IREAL8(NATOM))
      YTEMP=ALLHP(IREAL8(NATOM))
      ZTEMP=ALLHP(IREAL8(NATOM))
C
C make the list of the angles to be searched
      CALL MAKTHE('STOR',NR,MAXNR,ROTMOD,QSELF,
     & T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     & DELTA,HEAP(T1),HEAP(T2),HEAP(T3),
     & HEAP(N1DIM),N2DIM,HEAP(N3DIM))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' RSEARC: CPU-Time:  generation of angle grid=',SECT
      END IF
      WRITE(6,'(A,I8)')
     & ' DRSEAR: Number of rotations to be searched: ',NR
C carry out rotation search -- store values of RF search in RFCORR
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL RFDIRE(NR,HEAP(T1),HEAP(T2),
     & HEAP(T3),HEAP(RFCORR),
     & HEAP(XTEMP),HEAP(YTEMP),HEAP(ZTEMP),
     & RAVE,RSIGMA,RMIN,RMAX,ROTMOD)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' DRSEAR: CPU-Time:  rotation search=',SECT
      END IF
C
C switch on symmetry operators for cluster analysis
      XROLD=XRNSYM
      XRNSYM=TSNSYM
C
C write rotation function / print summary and statistics
      CALL RFWRIT(OUNIT,NR,HEAP(T1),HEAP(T2),HEAP(T3),HEAP(RFCORR),
     & RAVE,RSIGMA,RMIN,RMAX,
     & ROTMOD,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     & HEAP(N1DIM),N2DIM,HEAP(N3DIM))
C
C number of highest peaks to be analysed can't be larger than NR
      NLIST=MIN(NLIST,NR)
C
C make a cluster analysis of the NLIST highest peaks and
C mark redundant peaks
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL RFCLUS(LUNIT,EPS,
     & NR,HEAP(T1),HEAP(T2),HEAP(T3),HEAP(RFCORR),
     & NLIST,ROTMOD,QDEBUG,QSELF,
     & T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &             RFSYMB,RFSYML)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' DRSEAR: CPU-Time:  cluster analysis=',SECT
      END IF
C
C restore original symmetry (should be identity operator)
      XRNSYM=XROLD
C
C free up space for rotation function values
      CALL FREHP(N3DIM,INTEG4(N2DIM))
      CALL FREHP(N1DIM,INTEG4(N2DIM))
      CALL FREHP(T3,IREAL8(MAXNR))
      CALL FREHP(T2,IREAL8(MAXNR))
      CALL FREHP(T1,IREAL8(MAXNR))
      CALL FREHP(RFCORR,IREAL8(MAXNR))
      CALL FREHP(XTEMP,IREAL8(NATOM))
      CALL FREHP(YTEMP,IREAL8(NATOM))
      CALL FREHP(ZTEMP,IREAL8(NATOM))
C
      END IF
      END IF
C
9999  CONTINUE
      RETURN
      END
C
      SUBROUTINE RFDIRE(NR,T1,T2,T3,RFCORR,
     & XTEMP,YTEMP,ZTEMP,
     & RAVE,RSIGMA,RMIN,RMAX,ROTMOD)
C
C Perform direct rotation search
C
C Authors: Warren DeLano and Axel T. Brunger
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'timer.inc'
      INTEGER NR
      DOUBLE PRECISION T1(*),T2(*),T3(*)
      DOUBLE PRECISION RFCORR(*),XTEMP(*),YTEMP(*),ZTEMP(*)
      DOUBLE PRECISION RAVE,RSIGMA,RMIN,RMAX
      CHARACTER*4 ROTMOD
C local
      INTEGER I, J
      DOUBLE PRECISION ROT(3,3), RAVE2
      DOUBLE PRECISION AXIS(3), Q(4)
      DOUBLE PRECISION CORR
      DOUBLE COMPLEX DBCOMP
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0, ONE=1.0)
C begin
C
C initialize statistics for rotation search
      RAVE=ZERO
      RAVE2=ZERO
      RMAX=-ONE
      RMIN=+ONE
C
      DO I=1,NATOM
      XTEMP(I)=X(I)
      YTEMP(I)=Y(I)
      ZTEMP(I)=Z(I)
      END DO
C
C loop over all grid points of rotation search
      DO J=1,NR
C
C compute unitary rotation matrix at this grid point
      CALL ROTMAT(ROT,T1(J),T2(J),T3(J),Q,AXIS,ROTMOD)
C
      DO I=1,NATOM
      X(I)=ROT(1,1)*XTEMP(I)+ROT(1,2)*YTEMP(I)+ROT(1,3)*ZTEMP(I)
      Y(I)=ROT(2,1)*XTEMP(I)+ROT(2,2)*YTEMP(I)+ROT(2,3)*ZTEMP(I)
      Z(I)=ROT(3,1)*XTEMP(I)+ROT(3,2)*YTEMP(I)+ROT(3,3)*ZTEMP(I)
      END DO
C
C compute ASSOciated atom selections and target function,
C but no derivatives and no monitor function.
      CALL XCALCS(.FALSE.,.FALSE.,XRE,
     &           MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &           XRTR,XRINTR,XRSCAL,XCVTEST,
     &           TRPNMX,TRPNX,TRPN,TRPNL,TRPNN,TRPNDB,
     &           TRPNMLT,TRPNLEV,TRPNTYP,TRPNDOM,TDEPTH,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           1,1,' ',1,0,DBCOMP,1,1,' ',' ',1,
     &           SRPNMX,SRPNX,SRPN,SRPNL,SRPNN,SRPNDB,
     &           SRPNMLT,SRPNLEV,SRPNTYP,SRPNDOM,SDEPTH,
     &           CRPNMX,CRPNX,CRPN,CRPNL,CRPNN,CRPNDB,
     &           CRPNMLT,CRPNLEV,CRPNTYP,CRPNDOM,CDEPTH,
     &           XRNREF,XRMREF,HPH,HPK,HPL,HPMULT,HPTYPE,
     &           HPTSEL,QHERM,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,
     &           XRSYGP, XRSYIV,
     &           XRCELL,XRVOL,
     &           HPANOMFLAG,HPATOF,HPINDF,XNAMEAS,XASSOC,
     &           QASELE,XRNATF,QFFT,
     &           HPDX,HPDY,HPDZ,HPDT,HPDQ,
     &           XSFMX,XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &           XSFGNAM,XSFGTYP,XSFGORD,
     &           XRSN,XRSM,XRSA,XRSB,XRSC,XRFP,XRFDP,QLOOK,
     &           MAPR,IHPFLAG,HPFLAG,MAXHPFLAG,
     &           XRRED,XRREUP)
      XRQCHK=.TRUE.
C
C compute correlation (negative for sorting)
      CORR=ONE-(XRE/XRSCAL)
C
      IF (WRNLEV.GE.10) THEN
      IF (ROTMOD.EQ.'EULE') THEN
      WRITE(6,'(3(A,F9.3),A,F9.5)') ' DRSEAR: T1=',T1(J),
     & ' T2=',T2(J),' T3=',T3(J),' CORR=', CORR
      ELSE IF(ROTMOD.EQ.'SPHE') THEN
      WRITE(6,'(3(A,F9.3),A,F9.5)') ' DRSEAR: PSI=',T1(J),
     & ' PHI=',T2(J),' KAPPA=',T3(J),' CORR=', CORR
      ELSE IF(ROTMOD.EQ.'LATT') THEN
      WRITE(6,'(3(A,F9.3),A,F9.5)') ' DRSEAR: TP=',T1(J),
     & ' T2=',T2(J),' TM=',T3(J),' CORR=', CORR
      END IF
      END IF
C
C update statistics of rotation search
      RAVE=RAVE+CORR
      RAVE2=RAVE2+CORR**2
      RMAX=MAX(RMAX,CORR)
      RMIN=MIN(RMIN,CORR)
C
C store negative of value for sorting
      RFCORR(J)=-CORR
C
      END DO
C
      DO I=1,NATOM
      X(I)=XTEMP(I)
      Y(I)=YTEMP(I)
      Z(I)=ZTEMP(I)
      END DO
C
C compute average and sigma of rotation function
      RAVE=RAVE/NR
      RSIGMA=SQRT(MAX(ZERO,RAVE2/NR-RAVE**2))
C
      RETURN
      END
C============================================================================
