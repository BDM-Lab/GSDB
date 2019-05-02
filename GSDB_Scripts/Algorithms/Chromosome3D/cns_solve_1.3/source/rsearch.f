      SUBROUTINE RSEARC
C
C Rotation search of two Patterson maps
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'expression.inc'
C local
      DOUBLE PRECISION T1MIN, T1MAX, T2MIN, T2MAX, T3MIN, T3MAX
      DOUBLE PRECISION DELTA
      LOGICAL QCORR
      INTEGER DUMMY
      DOUBLE PRECISION VLOW, VHIGH, THRES, TEMP, EPS
      INTEGER OUNIT, LUNIT, NLIST, NPEAKS
C
      INTEGER MPEAKS
C
      INTEGER MAXNR
      DOUBLE PRECISION RAVE, RSIGMA, RMIN, RMAX
C
      CHARACTER*(WORD_SIZE) LFILE
C
      CHARACTER*(WORD_SIZE) P1FILE
      INTEGER P1NA,P1AMIN,P1AMAX,P1NB,P1BMIN,P1BMAX,P1NC,P1CMIN,P1CMAX
      INTEGER P1DIM, P1UNIT
      DOUBLE PRECISION P1CELL(9), P1TR(3,3), P1INTR(3,3), P1VOL
C
      CHARACTER*(WORD_SIZE) P2FILE
      INTEGER P2NA,P2AMIN,P2AMAX,P2NB,P2BMIN,P2BMAX,P2NC,P2CMIN,P2CMAX
      INTEGER P2DIM, P2UNIT
      DOUBLE PRECISION P2CELL(9), P2TR(3,3), P2INTR(3,3), P2VOL
C
      CHARACTER*4 STARG, ROTMOD
C
      INTEGER P2TDIM
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
      INTEGER PEAKS, PEAKX, PEAKY, PEAKZ, N1DIM, N3DIM
      INTEGER NR, T1, T2, T3, RFCORR, P1MAP, P2MAP
      INTEGER P2TMAP, WORKA, WORKB
C parameters
      DOUBLE PRECISION ZERO, TWO, THREE, T1P4
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, THREE=3.0D0, T1P4=1.4D0)
C begin
C
C default values
      OFILE=' '
      LFILE=' '
      P1FILE=' '
      P2FILE=' '
      OUNIT=0
      LUNIT=0
      P1UNIT=0
      P2UNIT=0
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
      QCORR=.FALSE.
      VLOW=3.0D0
      VHIGH=20.0D0
      NPEAKS=5000
      NLIST=100
      THRES=500.0D0
      EPS=0.2D0
      QDEBUG=.FALSE.
      QSELF=.FALSE.
      QFORM=.TRUE.
      RFSYMB=' '
      RFSYML=1
C
C parsing
      CALL PUSEND('RSEARCH>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('RSEARCH>')
      CALL MISCOM('RSEARCH>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-search-rotation')
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
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'LIST') THEN
      CALL NEXTFI('LIST=',LFILE)
      ELSE IF (WD(1:4).EQ.'P1IN') THEN
      CALL NEXTFI('P1INput=',P1FILE)
      ELSE IF (WD(1:4).EQ.'P2IN') THEN
      CALL NEXTFI('P2INput=',P2FILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'RANG') THEN
      CALL NEXTF('lower-cutoff=',VLOW)
      CALL NEXTF('upper-cutoff=',VHIGH)
      IF (VLOW.GT.VHIGH) THEN
      TEMP=VLOW
      VLOW=VHIGH
      VHIGH=TEMP
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'THRE') THEN
      CALL NEXTF('THREshold=',THRES)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NPEA') THEN
      CALL NEXTI('NPEAks=',NPEAKS)
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
      ELSE IF (WD(1:4).EQ.'SELF') THEN
      CALL NEXTLO('SELF_symmetry=',QSELF)
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
      ELSE IF (WD(1:4).EQ.'TARG') THEN
      IF (QCORR) THEN
      STARG='CORR'
      ELSE
      STARG='PROD'
      END IF
      CALL NEXTA4('TARGet=',STARG)
      IF (STARG.EQ.'CORR') THEN
      QCORR=.TRUE.
      ELSE IF (STARG.EQ.'PROD') THEN
      QCORR=.FALSE.
      ELSE
      CALL DSPERR('RSEARCH>','only CORRelation or PRODuct allowed')
      END IF
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
      WRITE(6,'(A,F12.5,A,L1,A,L1)')
     & ' | DELTa=',DELTA,' SELF_symmetry=',QSELF,
     & ' FORMatted=',QFORM
      WRITE(6,'(A,F10.5,A,F10.5)')
     & ' | vector-RANGe-of-P1= ',VLOW,' to ',VHIGH
      WRITE(6,'(A,F12.6)')
     & ' | THREshold-of-peaks-of-P1=',THRES
      IF (QCORR) THEN
      WRITE(6,'(A)')
     & ' | TARGet=CORRelation'
      ELSE
      WRITE(6,'(A)')
     & ' | TARGet=PRODuct'
      END IF
      WRITE(6,'(A,I8)')
     & ' | NPEAks-of-P1=',NPEAKS
      WRITE(6,'(A,I8)')
     & ' | NLISt=',NLIST
      WRITE(6,'(A,F5.2)')
     & ' | EPSIlon=',EPS
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SYMM') THEN
      CALL NEXTEX('SYMMetry=',WDD,WDMAX,WDDLEN)
      CALL XRSYPA(WDD,WDDLEN,XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
C=====================================================================
C dummy parse asymmetry command (is not stored in spacegroup asymmetric unit)
C
      ELSE IF (WD(1:4).EQ.'ASYM') THEN
      CALL XASYMM(RPNMX,RPNN3,RPNX,RPN3,RPNL3,RPNDB3,RPNMLT3,
     &            RPNTYP3,RPNDOM3,RPNLEV3,DUMMY)
C=====================================================================
      ELSE
      CALL CHKEND('RSEARCH>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
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
      IF (P1FILE.NE.' '.AND.P2FILE.NE.' ') THEN
C
C allocate space for "permanent" peak list of map P1
      MPEAKS=NPEAKS
      PEAKS=ALLHP(IREAL8(MPEAKS))
      PEAKX=ALLHP(IREAL8(MPEAKS))
      PEAKY=ALLHP(IREAL8(MPEAKS))
      PEAKZ=ALLHP(IREAL8(MPEAKS))
C
C read header of map P1
      CALL RMAP1(P1FILE,P1UNIT,ERROR,EOF,
     &           P1NA,P1AMIN,P1AMAX,P1NB,P1BMIN,P1BMAX,
     &           P1NC,P1CMIN,P1CMAX,P1CELL,QDEBUG,QFORM)
      IF (ERROR.OR.EOF) GOTO 9999
C
C allocate space for map P1
      P1DIM=(P1AMAX-P1AMIN+1)*(P1BMAX-P1BMIN+1)*(P1CMAX-P1CMIN+1)
      P1MAP=ALLHP(IREAL8(P1DIM))
C
C read map P1
      CALL RMAP2(P1UNIT,ERROR,EOF,
     &           P1AMIN,P1AMAX,P1BMIN,P1BMAX,P1CMIN,P1CMAX,
     &           HEAP(P1MAP),QFORM)
      IF (ERROR.OR.EOF) GOTO 9999
C
C get fractional<->orthogonal transformation matrices
C (P1INTR transforms fractional to orthogonal A coordinates)
      CALL XRFRAC(P1CELL,P1TR,P1INTR,P1VOL)
C
C convert transformation matrices into the index space of the map
C that is, P1INTR will transform the map index into orthogonal
C A coordinates.
      CALL CONIND(P1TR,P1INTR,P1NA,P1NB,P1NC)
C
C make a list of NPEAKS highest peaks of P1 with distances
C from the origin that fall within the range VLOW and VHIGH
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL SELPEK(VLOW,VHIGH,THRES,NPEAKS,P1CELL,P1NA,P1NB,P1NC,
     &            P1AMAX,P1AMIN,P1BMAX,P1BMIN,P1CMAX,P1CMIN,
     &            HEAP(P1MAP),P1TR,P1INTR,
     &            HEAP(PEAKS),HEAP(PEAKX),HEAP(PEAKY),HEAP(PEAKZ),
     &            QDEBUG)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' RSEARC: CPU-Time:  P1 peak selection=',SECT
      END IF
C
C free up space for map P1
      CALL FREHP(P1MAP,IREAL8(P1DIM))
C
C determine the number of angles to be searched
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL MAKTHE('SCAN',NR,MAXNR,ROTMOD,QSELF,
     &     T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &     DELTA,0,0,0,0,N2DIM,0)
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
C
C make the list of the angles to be searched
      CALL MAKTHE('STOR',NR,MAXNR,ROTMOD,QSELF,
     &     T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &     DELTA,HEAP(T1),HEAP(T2),HEAP(T3),
     &     HEAP(N1DIM),N2DIM,HEAP(N3DIM))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' RSEARC: CPU-Time:  generation of angle grid=',SECT
      END IF
C
C read header of map P2
      CALL RMAP1(P2FILE,P2UNIT,ERROR,EOF,
     &           P2NA,P2AMIN,P2AMAX,P2NB,P2BMIN,P2BMAX,
     &           P2NC,P2CMIN,P2CMAX,P2CELL,QDEBUG,QFORM)
      IF (ERROR.OR.EOF) GOTO 9999
C
C allocate permanent space for map P2
      P2DIM=(P2NA+1)*(P2NB+1)*(P2NC+1)
      P2MAP=ALLHP(IREAL8(P2DIM))
C
C check whether map has to be "shuffled"
      IF (P2AMAX.NE.P2NA.OR.P2AMIN.NE.0.OR.
     &    P2BMAX.NE.P2NB.OR.P2BMIN.NE.0.OR.
     &    P2CMAX.NE.P2NC.OR.P2CMIN.NE.0 ) THEN
C
      WRITE(6,'(A)') ' RSEARC: shuffling map 2'
C
C allocate temporary (work) space for map P2
      P2TDIM=(P2AMAX-P2AMIN+1)*(P2BMAX-P2BMIN+1)*(P2CMAX-P2CMIN+1)
      P2TMAP=ALLHP(IREAL8(P2TDIM))
      WORKA=ALLHP(INTEG4(P2AMAX-P2AMIN+1))
      WORKB=ALLHP(INTEG4(P2BMAX-P2BMIN+1))
C
C read map P2
      CALL RMAP2(P2UNIT,ERROR,EOF,
     &           P2AMIN,P2AMAX,P2BMIN,P2BMAX,P2CMIN,P2CMAX,
     &           HEAP(P2TMAP),QFORM)
      IF (ERROR.OR.EOF) GOTO 9999
C
C shuffle the map P2 -- we want a complete primary unit cell
C with indices running from 1 to NA, 1 to NB, 1 to NC.
      CALL MAPSHF(P2NA,P2NB,P2NC,
     &            P2AMIN,P2AMAX,P2BMIN,P2BMAX,P2CMIN,P2CMAX,
     &            HEAP(P2TMAP),HEAP(P2MAP),HEAP(WORKA),HEAP(WORKB),
     &            QDEBUG)
C
C free up temporary space for map P2
      CALL FREHP(WORKB,INTEG4(P2BMAX-P2BMIN+1))
      CALL FREHP(WORKA,INTEG4(P2AMAX-P2AMIN+1))
      CALL FREHP(P2TMAP,IREAL8(P2TDIM))
      ELSE
C
C read map P2 directly
      CALL RMAP2(P2UNIT,ERROR,EOF,
     &           P2AMIN,P2AMAX,P2BMIN,P2BMAX,P2CMIN,P2CMAX,
     &           HEAP(P2MAP),QFORM)
      END IF
C
C get fractional<->orthogonal transformation matrices
C (P2INTR transforms fractional to orthogonal coordinates)
      CALL XRFRAC(P2CELL,P2TR,P2INTR,P2VOL)
C
C convert transformation matrices into the index space of the map
C that is, P1INTR will transform the map index into orthogonal
C A coordinates.
      CALL CONIND(P2TR,P2INTR,P2NA,P2NB,P2NC)
C
C carry out rotation search -- store values of RF search in RFCORR
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL XPATR2(NR,HEAP(T1),HEAP(T2),HEAP(T3),HEAP(RFCORR),
     &     NPEAKS,HEAP(PEAKS),HEAP(PEAKX),HEAP(PEAKY),HEAP(PEAKZ),
     &     P2NA,P2NB,P2NC,HEAP(P2MAP),P2TR,
     &     RAVE,RSIGMA,RMIN,RMAX,QCORR,ROTMOD)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' RSEARC: CPU-Time:  rotation search=',SECT
      END IF
C
C free up space for permanent MAP P2
      CALL FREHP(P2MAP,IREAL8(P2DIM))
C
C free up space for "permanent" peak list
      CALL FREHP(PEAKZ,IREAL8(MPEAKS))
      CALL FREHP(PEAKY,IREAL8(MPEAKS))
      CALL FREHP(PEAKX,IREAL8(MPEAKS))
      CALL FREHP(PEAKS,IREAL8(MPEAKS))
C
C write rotation function / print summary and statistics
      CALL RFWRIT(OUNIT,NR,HEAP(T1),HEAP(T2),HEAP(T3),HEAP(RFCORR),
     &     RAVE,RSIGMA,RMIN,RMAX,
     &     ROTMOD,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &     HEAP(N1DIM),N2DIM,HEAP(N3DIM))
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
     &            NR,HEAP(T1),HEAP(T2),HEAP(T3),HEAP(RFCORR),
     &            NLIST,ROTMOD,QDEBUG,QSELF,
     &            T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &            RFSYMB,RFSYML)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' RSEARC: CPU-Time:  cluster analysis=',SECT
      END IF
C
C free up space for rotation function values
      CALL FREHP(N3DIM,INTEG4(N2DIM))
      CALL FREHP(N1DIM,INTEG4(N2DIM))
      CALL FREHP(T3,IREAL8(MAXNR))
      CALL FREHP(T2,IREAL8(MAXNR))
      CALL FREHP(T1,IREAL8(MAXNR))
      CALL FREHP(RFCORR,IREAL8(MAXNR))
C
C error label
      END IF
9999  CONTINUE
      RETURN
      END
C=====================================================================
      SUBROUTINE SELPEK(VLOW,VHIGH,THRES,NPEAKS,CELL,NA,NB,NC,
     &                  AMAX,AMIN,BMAX,BMIN,CMAX,CMIN,
     &                  MAP,TR,INTR,
     &                  PEAKS,PEAKX,PEAKY,PEAKZ,QDEBUG)
C
C Front-end routine for peak selection routine SELPE2
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      DOUBLE PRECISION VLOW, VHIGH, THRES
      INTEGER NPEAKS
      DOUBLE PRECISION CELL(9)
      INTEGER NA, NB, NC
      INTEGER AMAX, AMIN, BMAX, BMIN, CMAX, CMIN
      DOUBLE PRECISION MAP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      DOUBLE PRECISION TR(3,3), INTR(3,3)
      DOUBLE PRECISION PEAKS(*), PEAKX(*), PEAKY(*), PEAKZ(*)
      LOGICAL QDEBUG
C local
      INTEGER TPDIM, TPEAKS, TPEAKX, TPEAKY, TPEAKZ, PERM
      INTEGER A, B, C
      INTEGER AMAXP, AMINP, BMAXP, BMINP, CMAXP, CMINP
      DOUBLE PRECISION GRIDX, GRIDY, GRIDZ, DIST2, VLOW2, VHIGH2
C begin
C
C Determine AMAXP, AMAXN, BMAXP, AMINP, BMINP, CMINP of a
C box isomorphous to the unitcell geometry which surrounds
C the sphere with radius VHIGH.  Then the box is reduced to
C fit into the extent of the map AMAX, BMAX, CMAX,
C AMIN, BMIN, CMIN.   Note: CELL(7), CELL(8), CELL(9) contain
C the norms of a*, b*, c* .
      AMAXP=MIN(AMAX,INT(VHIGH*CELL(7)*NA+RSMALL)+1)
      BMAXP=MIN(BMAX,INT(VHIGH*CELL(8)*NB+RSMALL)+1)
      CMAXP=MIN(CMAX,INT(VHIGH*CELL(9)*NC+RSMALL)+1)
      AMINP=MAX(AMIN,-INT(VHIGH*CELL(7)*NA+RSMALL)-1)
      BMINP=MAX(BMIN,-INT(VHIGH*CELL(8)*NB+RSMALL)-1)
      CMINP=MAX(CMIN,-INT(VHIGH*CELL(9)*NC+RSMALL)-1)
C
      TPDIM=0
C
C determine the number of selected peaks.  We need this number in order
C to allocate space from the heap.  Therefore, we have to search through
C all peaks again after the heap space has been allocated.  This sounds
C pretty inefficient, but is the only memory-efficient method since
C FORTRAN doesn't support dynamic expansion of lists or arrays.
C
C compute squared distances for selection criterium
      VLOW2=VLOW**2
      VHIGH2=VHIGH**2
C
C loop through all peaks of the reduced box
      DO C=CMINP,CMAXP
      DO B=BMINP,BMAXP
      DO A=AMINP,AMAXP
C
C get orthogonal A coordinate for the grid point
      GRIDX=INTR(1,1)*A+INTR(1,2)*B+INTR(1,3)*C
      GRIDY=INTR(2,1)*A+INTR(2,2)*B+INTR(2,3)*C
      GRIDZ=INTR(3,1)*A+INTR(3,2)*B+INTR(3,3)*C
C
C check whether the peak value of the grid point is within the
C specified range and check whether the peak value of the
C grid point is greater than threshold
      DIST2=GRIDX**2+GRIDY**2+GRIDZ**2
      IF (VLOW2.LE.DIST2.AND.DIST2.LE.VHIGH2.AND.MAP(A,B,C).GT.THRES)
     &   THEN
      TPDIM=TPDIM+1
      END IF
      END DO
      END DO
      END DO
C
C allocate space for temporary peak list
      TPEAKS=ALLHP(IREAL8(TPDIM))
      TPEAKX=ALLHP(IREAL8(TPDIM))
      TPEAKY=ALLHP(IREAL8(TPDIM))
      TPEAKZ=ALLHP(IREAL8(TPDIM))
      PERM=ALLHP(INTEG4(TPDIM))
C
C call routine that actually carries out the peak selection
      CALL SELPE2(VLOW,VHIGH,THRES,NPEAKS,
     &            AMAX,AMIN,BMAX,BMIN,CMAX,CMIN,
     &            MAP,TR,INTR,
     &            AMAXP,AMINP,BMAXP,BMINP,CMAXP,CMINP,
     &            TPDIM,HEAP(TPEAKS),HEAP(TPEAKX),HEAP(TPEAKY),
     &            HEAP(TPEAKZ),PEAKS,PEAKX,PEAKY,PEAKZ,HEAP(PERM),
     &            QDEBUG)
C
C free up space of for the temporary peak list
      CALL FREHP(PERM,INTEG4(TPDIM))
      CALL FREHP(TPEAKZ,IREAL8(TPDIM))
      CALL FREHP(TPEAKY,IREAL8(TPDIM))
      CALL FREHP(TPEAKX,IREAL8(TPDIM))
      CALL FREHP(TPEAKS,IREAL8(TPDIM))
      RETURN
      END
C============================================================================
      SUBROUTINE SELPE2(VLOW,VHIGH,THRES,NPEAKS,
     &                  AMAX,AMIN,BMAX,BMIN,CMAX,CMIN,
     &                  MAP,TR,INTR,
     &                  AMAXP,AMINP,BMAXP,BMINP,CMAXP,CMINP,
     &                  TPDIM,TPEAKS,TPEAKX,TPEAKY,TPEAKZ,
     &                  PEAKS,PEAKX,PEAKY,PEAKZ,PERM,QDEBUG)
C
C Makes a list of NPEAKS highest peaks of map larger than
C THRES with distances from the origin that fall within
C the range VLOW and VHIGH.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION VLOW, VHIGH, THRES
      INTEGER NPEAKS
      INTEGER AMAX, AMIN, BMAX, BMIN, CMAX, CMIN
      DOUBLE PRECISION MAP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      DOUBLE PRECISION TR(3,3), INTR(3,3)
      INTEGER AMAXP, AMINP, BMAXP, BMINP, CMAXP, CMINP
      INTEGER TPDIM
      DOUBLE PRECISION TPEAKS(*), TPEAKX(*), TPEAKY(*), TPEAKZ(*)
      DOUBLE PRECISION PEAKS(*), PEAKX(*), PEAKY(*), PEAKZ(*)
      INTEGER PERM(*)
      LOGICAL QDEBUG
C external
      EXTERNAL ORDRRR
C local
      INTEGER A, B, C, IPEAK, I, IRANGE
      DOUBLE PRECISION VLOW2, VHIGH2, GRIDX, GRIDY, GRIDZ, DIST2
C begin
C
C IPEAK is the number of selected peaks
      IPEAK=0
C
C IRANGE is the number of peaks within the vector length range
      IRANGE=0
C
C compute squared distances for selection criterium
      VLOW2=VLOW**2
      VHIGH2=VHIGH**2
C
C loop through all peaks of the reduced box
      DO C=CMINP,CMAXP
      DO B=BMINP,BMAXP
      DO A=AMINP,AMAXP
C
C get orthogonal A coordinate for the grid point
      GRIDX=INTR(1,1)*A+INTR(1,2)*B+INTR(1,3)*C
      GRIDY=INTR(2,1)*A+INTR(2,2)*B+INTR(2,3)*C
      GRIDZ=INTR(3,1)*A+INTR(3,2)*B+INTR(3,3)*C
C
C check whether the peak value of the grid point is within the
C specified range and check whether the peak value of the
C grid point is greater than threshold
      DIST2=GRIDX**2+GRIDY**2+GRIDZ**2
      IF (VLOW2.LE.DIST2.AND.DIST2.LE.VHIGH2) IRANGE=IRANGE+1
      IF (VLOW2.LE.DIST2.AND.DIST2.LE.VHIGH2.AND.MAP(A,B,C).GT.THRES)
     &   THEN
      IPEAK=IPEAK+1
C
C store negative grid point (for the sake of the sorting routine)
      TPEAKS(IPEAK)=-MAP(A,B,C)
      TPEAKX(IPEAK)=GRIDX
      TPEAKY(IPEAK)=GRIDY
      TPEAKZ(IPEAK)=GRIDZ
      END IF
      END DO
      END DO
      END DO
C
      IF (IPEAK.GT.TPDIM) THEN
      CALL WRNDIE(-5,'SELPE2','Fatal coding error: IPEAK > TPDIM')
      END IF
C
      WRITE(6,'(A,I8,A,F10.5,/,A,F10.5,A)')
     &' XPATRT: ',IRANGE,
     & ' peaks of P1 have coordinates between',
     &  VLOW,' and ',VHIGH,' A relative to the origin.'
      WRITE(6,'(A,I8,A,G14.6)')
     &' XPATRT: of those ',IPEAK,
     & ' peaks have values that are above',THRES
C
C determine permutation PERM of the temporary peak list
C (the first element of the permutation will be the smallest)
      CALL SORTP(IPEAK,PERM,ORDRRR,TPEAKS,0,0,0,0,0,0,0)
C
C copy the NPEAKS highest peaks into the permanent peak list
      NPEAKS=MIN(IPEAK,NPEAKS)
      DO I=1,NPEAKS
C
C need the minus sign (see above)
      PEAKS(I)=-TPEAKS(PERM(I))
      PEAKX(I)=TPEAKX(PERM(I))
      PEAKY(I)=TPEAKY(PERM(I))
      PEAKZ(I)=TPEAKZ(PERM(I))
      END DO
C
      WRITE(6,'(A,I8,A)')
     &' XPATRT: of those the ',NPEAKS,' highest peaks will be used'
C
      IF (QDEBUG) THEN
      WRITE(6,'(A)') ' P1-TR'
      WRITE(6,'(3F10.5)') TR
      WRITE(6,'(A)') ' P1-INTR'
      WRITE(6,'(3F10.5)') INTR
      WRITE(6,'(A)') ' P1 map in reduced box'
      WRITE(6,'(6I6)') AMAXP,AMINP,BMAXP,BMINP,CMAXP,CMINP
      WRITE(6,'(A)') ' Writing z=0 section of map P1'
      DO B=BMINP,BMAXP
      WRITE(6,'(6E12.5)') (MAP(A,B,MAX(CMINP,0)),A=AMINP,AMAXP)
      END DO
C
      WRITE(6,'(A)') ' List of peaks: x,y,z, value'
      DO I=1,NPEAKS
C
C get index the grid point
      A=INT(TR(1,1)*PEAKX(I)+TR(1,2)*PEAKY(I)+TR(1,3)*PEAKZ(I)
     &     +RSMALL+10000)-10000
      B=INT(TR(2,1)*PEAKX(I)+TR(2,2)*PEAKY(I)+TR(2,3)*PEAKZ(I)
     &     +RSMALL+10000)-10000
      C=INT(TR(3,1)*PEAKX(I)+TR(3,2)*PEAKY(I)+TR(3,3)*PEAKZ(I)
     &     +RSMALL+10000)-10000
      WRITE(6,'(3I10,3F12.5,G14.6)')
     &     A,B,C,PEAKX(I),PEAKY(I),PEAKZ(I),PEAKS(I)
      END DO
      END IF
C
      RETURN
      END
C============================================================================
      SUBROUTINE XPATR2(NR,T1,T2,T3,RFCORR,
     &           NPEAKS,PEAKS,PEAKX,PEAKY,PEAKZ,
     &           P2NA,P2NB,P2NC,P2MAP,P2TR,
     &           RAVE,RSIGMA,RMIN,RMAX,QCORR,ROTMOD)
C
C Routine rotates set of selected peaks P1 with respect to
C map P2 and computes correlation function between P1 and P2.
C A eight-point interpolation is carried out for the non-integral
C grid points.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER NR
      DOUBLE PRECISION T1(*), T2(*), T3(*), RFCORR(*)
      INTEGER NPEAKS
      DOUBLE PRECISION PEAKS(*), PEAKX(*), PEAKY(*), PEAKZ(*)
      INTEGER P2NA, P2NB, P2NC
      DOUBLE PRECISION P2MAP(0:P2NA,0:P2NB,0:P2NC), P2TR(3,3)
      DOUBLE PRECISION RAVE, RSIGMA, RMIN, RMAX
      LOGICAL QCORR
      CHARACTER*4 ROTMOD
C local
      INTEGER IT, I, J, K
      DOUBLE PRECISION ROT(3,3), RND(3,3), RAVE2
      DOUBLE PRECISION DX, DY, DZ, AXIS(3), Q(4)
      DOUBLE PRECISION CORR, INTERP
      DOUBLE PRECISION CI, CJ, CIJ, CII, CJJ
      INTEGER NX, NY, NZ
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
      IF (NPEAKS.GT.0.AND.NR.GT.0) THEN
C
C precompute parameters needed for the rotation correlation function
      CJ=ZERO
      CJJ=ZERO
      DO I=1,NPEAKS
      CJ=CJ+PEAKS(I)
      CJJ=CJJ+PEAKS(I)**2
      END DO
C
C loop over all grid points of rotation search
      DO IT=1,NR
C
C compute unitary rotation matrix at this grid point
      CALL ROTMAT(ROT,T1(IT),T2(IT),T3(IT),Q,AXIS,ROTMOD)
C
C compute the matrix product of the rotation matrix ROT
C and the matrix that converts orthogonal coordinates into
C index space of map P2
      DO I=1,3
      DO K=1,3
      RND(I,K)=ZERO
      DO J=1,3
      RND(I,K)=RND(I,K)+P2TR(I,J)*ROT(J,K)
      END DO
      END DO
      END DO
C
C initialize correlation coefficients
      CI=ZERO
      CII=ZERO
      CIJ=ZERO
C
C loop through all selected peaks of map P1
      DO I=1,NPEAKS
C
C apply combined rotation / transformation matrix to
C the position of the grid point of P1.
      DX=RND(1,1)*PEAKX(I)+RND(1,2)*PEAKY(I)+RND(1,3)*PEAKZ(I)
      DY=RND(2,1)*PEAKX(I)+RND(2,2)*PEAKY(I)+RND(2,3)*PEAKZ(I)
      DZ=RND(3,1)*PEAKX(I)+RND(3,2)*PEAKY(I)+RND(3,3)*PEAKZ(I)
C Get the index position of the lowest corner of the grid box of
C map P2 surrounding the peak of P1
      NX=-10000+INT(10000+DX+RSMALL)
      NY=-10000+INT(10000+DY+RSMALL)
      NZ=-10000+INT(10000+DZ+RSMALL)
C
C compute difference of the P2 grid point to the position of the peak in
C P1
      DX=DX-NX
      DY=DY-NY
      DZ=DZ-NZ
C
C map index back into box if necessary (NOTE: NX, NY, NZ
C are always well defined since the extent of the map P2
C goes from 0 to NX, 0 to NY, 0 to NZ. (see routine MAPSHF))
      NX=MOD(NX +10000*P2NA,P2NA)
      NY=MOD(NY +10000*P2NB,P2NB)
      NZ=MOD(NZ +10000*P2NC,P2NC)
C
C perform linear interpolation using an eight point interpolation
C formula
      INTERP= P2MAP(NX,NY,NZ)*(1-DX)*(1-DY)*(1-DZ)
     &       +P2MAP(NX+1,NY,NZ)*DX*(1-DY)*(1-DZ)
     &       +P2MAP(NX,NY+1,NZ)*(1-DX)*DY*(1-DZ)
     &       +P2MAP(NX+1,NY+1,NZ)*DX*DY*(1-DZ)
     &       +P2MAP(NX,NY,NZ+1)*(1-DX)*(1-DY)*DZ
     &       +P2MAP(NX+1,NY,NZ+1)*DX*(1-DY)*DZ
     &       +P2MAP(NX,NY+1,NZ+1)*(1-DX)*DY*DZ
     &       +P2MAP(NX+1,NY+1,NZ+1)*DX*DY*DZ
C
C accumulate correlation function
      CI=CI+INTERP
      CII=CII+INTERP**2
      CIJ=CIJ+PEAKS(I)*INTERP
C
C end loop over all selected peaks
      END DO
C
      IF (QCORR) THEN
C
C compute rotation correlation function
      CII=(CII-CI**2/NPEAKS)*(CJJ-CJ**2/NPEAKS)
      IF (ABS(CII).GT.RSMALL) THEN
      CORR=(CIJ - CI*CJ/NPEAKS ) / SQRT(CII)
      ELSE
      CORR=ZERO
      END IF
C
      ELSE
C
C compute rotation product function
      CORR=CIJ/NPEAKS
      END IF
C
C update statistics of rotation search
      RAVE=RAVE+CORR
      RAVE2=RAVE2+CORR**2
      RMAX=MAX(RMAX,CORR)
      RMIN=MIN(RMIN,CORR)
C
C store negative correlation for this set of rotation angles
C (the negative value is used for the sake of the sorting routine)
      RFCORR(IT)=-CORR
      END DO
C
C compute average and sigma of rotation function
      RAVE=RAVE/NR
      RSIGMA=SQRT(MAX(ZERO,RAVE2/NR-RAVE**2))
      END IF
C
      RETURN
      END
C
      SUBROUTINE CONIND(TR,INTR,NA,NB,NC)
C
C convert transformation matrices into the index space of the map
C that is, INTR will transform the map index into orthogonal
C A coordinates.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION TR(3,3), INTR(3,3)
      INTEGER NA, NB, NC
C local
      INTEGER J
C begin
      DO J=1,3
      TR(1,J)=NA*TR(1,J)
      END DO
      DO J=1,3
      TR(2,J)=NB*TR(2,J)
      END DO
      DO J=1,3
      TR(3,J)=NC*TR(3,J)
      END DO
C
      DO J=1,3
      INTR(J,1)=INTR(J,1)/NA
      END DO
      DO J=1,3
      INTR(J,2)=INTR(J,2)/NB
      END DO
      DO J=1,3
      INTR(J,3)=INTR(J,3)/NC
      END DO
C
      RETURN
      END
C============================================================================
      SUBROUTINE MAPSHF(NA,NB,NC,
     &            AMIN,AMAX,BMIN,BMAX,CMIN,CMAX,
     &            TMAP,MAP,WORKA,WORKB,QDEBUG)
C
C MAP "shuffle" routine.
C
C Routine copies map TMAP with extent (AMIN--AMAX;BMIN--BMAX;
C CMIN--CMAX) into map MAP with extent (0--NA;0--NB;0--NC).
C MAP will then represent the primary unit cell INCLUDING all
C planes and edges at the boundary of the unit cell.  The operation is
C carried out by projecting the indices of TMAP into the primary
C unit cell.  The inclusion of all planes and edges makes
C interpolations at the boundary of the unit cell more efficient
C (see routine XPATR2).
C
C WORKA is a work array of dimension AMAX-AMIN+1
C WORKB is a work array of dimension BMAX-BMIN+1
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER NA, NB, NC, AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
      DOUBLE PRECISION TMAP(AMIN:AMAX,BMIN:BMAX,CMIN:CMAX)
      DOUBLE PRECISION MAP(0:NA,0:NB,0:NC)
      INTEGER WORKA(AMIN:AMAX), WORKB(BMIN:BMAX)
      LOGICAL QDEBUG
C local
      INTEGER A, B, C, AT, BT, CT
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0)
C begin
C
C initialize map MAP
      DO AT=0,NA
      DO BT=0,NB
      DO CT=0,NC
      MAP(AT,BT,CT)=ZERO
      END DO
      END DO
      END DO
C
C pre-compute mod's for inner-most loops
      DO A=AMIN,AMAX
      WORKA(A)=MOD(A +10000*NA,NA)
      END DO
      DO B=BMIN,BMAX
      WORKB(B)=MOD(B +10000*NB,NB)
      END DO
C
      DO C=CMIN,CMAX
      CT=MOD(C +10000*NC,NC)
      DO B=BMIN,BMAX
      BT=WORKB(B)
      DO A=AMIN,AMAX
      AT=WORKA(A)
      MAP(AT,BT,CT)=TMAP(A,B,C)
      END DO
      END DO
      END DO
C
C now we have to generate the planes and edges that are at the
C boundary of the unit cell
      DO AT=0,NA-1
      DO BT=0,NB-1
      MAP(AT,BT,NC)=MAP(AT,BT,0)
      END DO
      END DO
      DO AT=0,NA-1
      DO CT=0,NC-1
      MAP(AT,NB,CT)=MAP(AT,0,CT)
      END DO
      END DO
      DO BT=0,NB-1
      DO CT=0,NC-1
      MAP(NA,BT,CT)=MAP(0,BT,CT)
      END DO
      END DO
C
      DO AT=0,NA-1
      MAP(AT,NB,NC)=MAP(AT,0,0)
      END DO
      DO BT=0,NB-1
      MAP(NA,BT,NC)=MAP(0,BT,0)
      END DO
      DO CT=0,NC-1
      MAP(NA,NB,CT)=MAP(0,0,CT)
      END DO
C
      MAP(NA,NB,NC)=MAP(0,0,0)
C
C
      IF (QDEBUG) THEN
      WRITE(6,'(A)') ' P2 map after shuffling'
      WRITE(6,'(6I6)') NA,NB,NC
      WRITE(6,'(A)') ' Writing z=0 section of map P2'
      DO BT=0,NB
      WRITE(6,'(6E12.5)') (MAP(AT,BT,0),AT=0,NA)
      END DO
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE RFWRIT(OUNIT,NR,T1,T2,T3,RFCORR,
     &           RAVE,RSIGMA,RMIN,RMAX,
     &           ROTMOD,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &           N1DIM,N2DIM,N3DIM)
C
C Writes/prints information about rotation search
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER OUNIT, NR
      DOUBLE PRECISION T1(*), T2(*), T3(*), RFCORR(*)
      DOUBLE PRECISION RAVE, RSIGMA, RMIN, RMAX
      CHARACTER*4 ROTMOD
      DOUBLE PRECISION T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX
      INTEGER N1DIM(*), N2DIM, N3DIM(*)
C local
      INTEGER I, IT1, IT2, IT3
C begin
C write complete rotation function as a Mathematica data structure
C to the specified file
      IF (OUNIT.GT.0) THEN
C
      WRITE(OUNIT,'(A)') ' (* Rotation Function *)'
C
C write statistics about R-function
      WRITE(OUNIT,'(A,I7,A,4(/,A,F8.4,A))')
     & '   rnumber=',NR,';',
     & '   rave=',RAVE,';',
     & '   rsigma=',RSIGMA,';',
     & '   rmax=',RMAX,';',
     & '   rmin=',RMIN,';'
C
C write variable for each level
      IF (ROTMOD.EQ.'LATT') THEN
      WRITE(OUNIT,'(A)') ' var={ "theta2", "theta-", "theta+" };'
      ELSE IF (ROTMOD.EQ.'EULE') THEN
      WRITE(OUNIT,'(A)') ' var={ "theta2", "theta3", "theta1" };'
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(OUNIT,'(A)') ' var={ "phi", "kappa", "psi" };'
      END IF
C
C write lower bounds for each level
      WRITE(OUNIT,'(3(A,F10.5),A)')
     &        ' min={ ',T2MIN,' , ',T3MIN,' , ',T1MIN,' };'
C
C write upper bounds for each level
      WRITE(OUNIT,'(3(A,F10.5),A)')
     &        ' max={ ',T2MAX,' , ',T3MAX,' , ',T1MAX,' };'
C
C write RF 3-d-matrix
      WRITE(OUNIT,'(A)') ' rf={'
      I=0
      DO IT2=1,N2DIM
      WRITE(OUNIT,'(A)') ' {'
      DO IT3=1,N3DIM(IT2)
      WRITE(OUNIT,'(A)') ' {'
      WRITE(OUNIT,'(6(F11.4,A))')
     &  (-RFCORR(IT1),',',IT1=I+1,I+N1DIM(IT2)-1),
     &   -RFCORR(I+N1DIM(IT2))
      I=I+N1DIM(IT2)
      IF (IT3.LT.N3DIM(IT2)) THEN
      WRITE(OUNIT,'(A)') ' },'
      ELSE
      WRITE(OUNIT,'(A)') ' }'
      END IF
      END DO
      IF (IT2.LT.N2DIM) THEN
      WRITE(OUNIT,'(A)') ' },'
      ELSE
      WRITE(OUNIT,'(A)') ' }'
      END IF
      END DO
      WRITE(OUNIT,'(A)') ' };'
      IF (I.NE.NR) CALL WRNDIE(-5,'RFWRIT','Fatal Coding Error')
      END IF
C
C print summary of rotation search to standard output
      WRITE(6,'(A,/,A,I6,4(/,A,F7.3))')
     & ' Summary of rotation search:',
     & '   number of grid points=',NR,
     & '   mean of R-function=',RAVE,
     & '   sigma of R-function=',RSIGMA,
     & '   maximum of R-function=',RMAX,
     & '   minimum of R-function=',RMIN
      RETURN
      END
C
      SUBROUTINE RFCLUS(LUNIT,EPS,NR,T1,T2,T3,RFCORR,NLIST,ROTMOD,
     &           QDEBUG,QSELF,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &           RFSYMB,RFSYML)
C
C Make a cluster analysis of the NLIST highest peaks.  Patterson
C symmetry is used.  Redundant peaks are marked.
C
C Front-end routine to RFCLU2
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER LUNIT
      DOUBLE PRECISION EPS
      INTEGER NR
      DOUBLE PRECISION T1(*), T2(*), T3(*), RFCORR(*)
      INTEGER NLIST
      CHARACTER*4 ROTMOD
      LOGICAL QDEBUG, QSELF
      DOUBLE PRECISION T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX
      CHARACTER*(*) RFSYMB
      INTEGER RFSYML
C local
      INTEGER ROT, SYMM, PRIO, PERM2, NEAR, REDUND, PERM, SYOP1, SYOP2
C begin
      PERM=ALLHP(INTEG4(NR))
      ROT=ALLHP(IREAL8(NLIST*9))
      SYMM=ALLHP(IREAL8(NLIST*9))
      PRIO=ALLHP(INTEG4(NLIST))
      PERM2=ALLHP(INTEG4(NLIST))
      NEAR=ALLHP(ILOGIC(NLIST))
      REDUND=ALLHP(ILOGIC(NLIST))
      SYOP1=ALLHP(INTEG4(NLIST))
      SYOP2=ALLHP(INTEG4(NLIST))
      CALL RFCLU2(LUNIT,EPS,
     &            NR,T1,T2,T3,RFCORR,NLIST,ROTMOD,HEAP(PERM),
     &            HEAP(ROT),HEAP(SYMM),HEAP(PRIO),HEAP(PERM2),
     &            HEAP(NEAR),HEAP(REDUND),HEAP(SYOP1),HEAP(SYOP2),
     &            QDEBUG,QSELF,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &            RFSYMB,RFSYML)
      CALL FREHP(SYOP2,INTEG4(NLIST))
      CALL FREHP(SYOP1,INTEG4(NLIST))
      CALL FREHP(REDUND,ILOGIC(NLIST))
      CALL FREHP(NEAR,ILOGIC(NLIST))
      CALL FREHP(PERM2,INTEG4(NLIST))
      CALL FREHP(PRIO,INTEG4(NLIST))
      CALL FREHP(SYMM,IREAL8(NLIST*9))
      CALL FREHP(ROT,IREAL8(NLIST*9))
      CALL FREHP(PERM,INTEG4(NR))
      RETURN
      END
C============================================================================
      SUBROUTINE RFCLU2(LUNIT,EPS,NR,T1,T2,T3,RFCORR,NLIST,ROTMOD,
     &             PERM,ROT,SYMM,PRIO,PERM2,NEAR,REDUND,SYOP1,SYOP2,
     &             QDEBUG,QSELF,T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &             RFSYMB,RFSYML)
C
C Make a cluster analysis of the NLIST highest peaks.  Patterson
C symmetry is used.  Redundant peaks are marked.  Peaks are
C printed to standard output as well as printed to LUNIT if
C LUNIT is assigned.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'symbol.inc'
      INCLUDE 'comand.inc'
      INTEGER LUNIT
      DOUBLE PRECISION EPS
      INTEGER NR
      DOUBLE PRECISION T1(*), T2(*), T3(*), RFCORR(*)
      INTEGER NLIST
      CHARACTER*4 ROTMOD
      INTEGER PERM(*)
      DOUBLE PRECISION ROT(3,3,NLIST), SYMM(NLIST,3,3)
      INTEGER PRIO(NLIST), PERM2(NLIST)
      LOGICAL NEAR(NLIST), REDUND(NLIST)
      INTEGER SYOP1(NLIST), SYOP2(NLIST)
      LOGICAL QDEBUG, QSELF
      DOUBLE PRECISION T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX
      CHARACTER*(*) RFSYMB
      INTEGER RFSYML
C local
      DOUBLE PRECISION DIFF, EPS2, AXIS(3), TT1, TT2, TT3, TTP, TTM
      INTEGER ISYM1, ISYM2, NSYM1, DET
      INTEGER U, V, X, Y, I, J, II, III
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC, AOP(3,3), BOP(3,3), Q(4)
      CHARACTER BUF*256
C external
      EXTERNAL ORDRRR, ORDER
C parameter
      DOUBLE PRECISION ZERO, ONE, SMALL, T360, TEN, TWO, SMALLR
      PARAMETER (ZERO=0.0, ONE=1.0, SMALL=0.0001, T360=360.D0)
      PARAMETER (TEN=10.D0, TWO=2.0D0, SMALLR=1.0D-6)
C begin
CSGI-specific optimization
C*$* OPTIMIZE(4)
C
C sort R-correlation function
      CALL SORTP(NR,PERM,ORDRRR,RFCORR,0,0,0,0,0,0,0)
C
C store highest value in symbol $RESULT
      DBPREC=-RFCORR(PERM(1))
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, DBPREC )
C
C precompute the square of epsilon for the matrix difference test
      EPS2=EPS**2
C
C loop through NLIST highest peaks (see routine RFWRIT)
      DO I=1,NLIST
C
C initialize priority list and redundancy indicator
      PRIO(I)=I
      REDUND(I)=.FALSE.
      SYOP1(I)=1
      SYOP2(I)=1
C
C compute unitary rotation matrix at this set of angles
      CALL ROTMAT(ROT(1,1,I),T1(PERM(I)),T2(PERM(I)),T3(PERM(I)),
     &     Q,AXIS,ROTMOD)
      END DO
C
C loop through all pairs of the NLIST peaks (triangular search)
      DO I=1,NLIST
C
C only look at peaks that are isolated so far.
      IF (PRIO(I).EQ.I) THEN
C
C loop over all crystallographic symmetry operators for map P1
      IF (QSELF) THEN
      NSYM1=XRNSYM
      ELSE
      NSYM1=1
      END IF
      DO ISYM1=1,NSYM1
C
C only use rotations (ATB 2/8/94)
      DET=(XRSYMM(ISYM1,1,1)*XRSYMM(ISYM1,2,2)
     &     -XRSYMM(ISYM1,1,2)*XRSYMM(ISYM1,2,1))*XRSYMM(ISYM1,3,3)
     &   +(XRSYMM(ISYM1,2,1)*XRSYMM(ISYM1,3,2)
     &     -XRSYMM(ISYM1,2,2)*XRSYMM(ISYM1,3,1))*XRSYMM(ISYM1,1,3)
     &   +(XRSYMM(ISYM1,3,1)*XRSYMM(ISYM1,1,2)
     &     -XRSYMM(ISYM1,3,2)*XRSYMM(ISYM1,1,1))*XRSYMM(ISYM1,2,3)
      IF (DET.EQ.1) THEN
C
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      BOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      BOP(U,V)=BOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
C
C loop over all crystallographic symmetry operators for map P2
      DO ISYM2=1,XRNSYM
C
C only use rotations (ATB 2/8/94)
      DET=(XRSYMM(ISYM2,1,1)*XRSYMM(ISYM2,2,2)
     &     -XRSYMM(ISYM2,1,2)*XRSYMM(ISYM2,2,1))*XRSYMM(ISYM2,3,3)
     &   +(XRSYMM(ISYM2,2,1)*XRSYMM(ISYM2,3,2)
     &     -XRSYMM(ISYM2,2,2)*XRSYMM(ISYM2,3,1))*XRSYMM(ISYM2,1,3)
     &   +(XRSYMM(ISYM2,3,1)*XRSYMM(ISYM2,1,2)
     &     -XRSYMM(ISYM2,3,2)*XRSYMM(ISYM2,1,1))*XRSYMM(ISYM2,2,3)
      IF (DET.EQ.1) THEN
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
C
C multiply P2 and P1 symmetry operators ( orthogonal A frame ) and
C rotation operator for each peak, i.e.,  SYMM=AOP*ROT*BOP
      DO U=1,3
      DO V=1,3
      DO J=I+1,NLIST
      SYMM(J,U,V)=ZERO
      END DO
      DO X=1,3
      DO Y=1,3
      DO J=I+1,NLIST
      SYMM(J,U,V)=SYMM(J,U,V)+AOP(U,X)*ROT(X,Y,J)*BOP(Y,V)
      END DO
      END DO
      END DO
      END DO
      END DO
C
      DO J=I+1,NLIST
C
C check whether rotation matrices are close
      DIFF= (SYMM(J,1,1)-ROT(1,1,I))**2
     &     +(SYMM(J,1,2)-ROT(1,2,I))**2
     &     +(SYMM(J,1,3)-ROT(1,3,I))**2
     &     +(SYMM(J,2,1)-ROT(2,1,I))**2
     &     +(SYMM(J,2,2)-ROT(2,2,I))**2
     &     +(SYMM(J,2,3)-ROT(2,3,I))**2
     &     +(SYMM(J,3,1)-ROT(3,1,I))**2
     &     +(SYMM(J,3,2)-ROT(3,2,I))**2
     &     +(SYMM(J,3,3)-ROT(3,3,I))**2
      REDUND(J)=REDUND(J).OR.DIFF.LT.SMALL
      NEAR(J)=DIFF.LT.EPS2
      END DO
C
C adjust priorities of peaks that are close to peak i
      DO J=I+1,NLIST
C if j is isolated so far --> set the priority of j to the priority of i
      IF (NEAR(J).AND.PRIO(J).EQ.J) THEN
      PRIO(J)=PRIO(I)
      SYOP1(J)=ISYM1
      SYOP2(J)=ISYM2
      END IF
      END DO
C
C end do symmetry operators
      END IF
      END DO
      END IF
      END DO
C
      END IF
C
C end do i loop
      END DO
C
C sort the priorities
      CALL SORTP(NLIST,PERM2,ORDER,PRIO,1,0,0,0,0,0,0)
C
C print the clustered and sorted list
      WRITE(6,'(A,F5.2,A,I6,A)')
     & ' Clustered list (EPSIlon=',EPS,') of first ',NLIST,
     & ' highest peaks'
      WRITE(6,'(A)') ' (the [] indicate a symmetry operation)'
      IF (ROTMOD.EQ.'EULE'.OR.ROTMOD.EQ.'LATT') THEN
      WRITE(6,'(A)')
     & '     no.  priority     theta1, theta2, theta3          RF'
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(6,'(A)')
     & '     no.  priority     psi,    phi,    kappa           RF'
      END IF
      WRITE(6,'(A)')
     & ' ==========================================================='
      II=PRIO(PERM2(1))
      III=1
      DO I=1,NLIST
      TT1=T1(PERM(PERM2(I)))
      TT2=T2(PERM(PERM2(I)))
      TT3=T3(PERM(PERM2(I)))
C
C convert Lattman angles into regular t1, t2, t3 angles
      IF (ROTMOD.EQ.'LATT') THEN
      TTP=TT1
      TTM=TT3
      TT1=(TTP+TTM)/TWO
      TT1=MOD(TT1+TEN*T360+SMALLR,T360)
      TT3=(TTP-TTM)/TWO
      TT3=MOD(TT3+TEN*T360+SMALLR,T360)
      END IF
C
      IF (PRIO(PERM2(I)).NE.II) THEN
      III=III+1
      II=PRIO(PERM2(I))
      WRITE(6,'(A)') ' '
      END IF
      IF (SYOP1(PERM2(I)).NE.1.OR.SYOP2(PERM2(I)).NE.1) THEN
      WRITE(6,'(A,I6,A,I6,A,I3,A,I3,A,3F8.3,A,F8.4)')
     & ' ',PERM2(I),
     & ' ',III,' [',SYOP1(PERM2(I)),',',SYOP2(PERM2(I)),'] (',
     & TT1,TT2,TT3,')    ',-RFCORR(PERM(PERM2(I)))
      ELSE IF (.NOT.REDUND(PERM2(I))) THEN
      WRITE(6,'(A,I6,A,I6,A,3F8.3,A,F8.4)')
     & ' ',PERM2(I),' ',III,'           (',
     & TT1,TT2,TT3,')    ',-RFCORR(PERM(PERM2(I)))
      END IF
      END DO
C
C write a summary to standard output
      WRITE(6,'(A,F5.2,A)')
     & ' List of the highest peaks of the clusters (EPSIlon=',
     &  EPS,')'
      IF (ROTMOD.EQ.'EULE'.OR.ROTMOD.EQ.'LATT') THEN
      WRITE(6,'(A)')
     & '     no.  priority     theta1, theta2, theta3          RF'
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(6,'(A)')
     & '     no.  priority     psi,    phi,    kappa           RF'
      END IF
      WRITE(6,'(A)')
     & ' ==========================================================='
C
C write a summary to the list file if requested
      IF (LUNIT.GT.0) THEN
      IF (ROTMOD.EQ.'EULE'.OR.ROTMOD.EQ.'LATT') THEN
      WRITE(LUNIT,'(A,F5.2,A)')
     & ' ! index, theta1, theta2, theta3, RF-function (EPSIlon=',
     &  EPS,')'
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(LUNIT,'(A,F5.2,A)')
     & ' ! index, psi, phi, kappa, RF-function (',
     &  EPS,')'
      END IF
      END IF
C
      II=0
      III=0
      DO I=1,NLIST
      IF (PRIO(PERM2(I)).NE.II) THEN
      III=III+1
      II=PRIO(PERM2(I))
      TT1=T1(PERM(PERM2(I)))
      TT2=T2(PERM(PERM2(I)))
      TT3=T3(PERM(PERM2(I)))
C
C map angles into asymmetric unit
      CALL ROTASY(QSELF,ROTMOD,TT1,TT2,TT3,
     &            T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX)
C
C convert Lattman angles into regular t1, t2, t3 angles
      IF (ROTMOD.EQ.'LATT') THEN
      TTP=TT1
      TTM=TT3
      TT1=(TTP+TTM)/TWO
      TT1=MOD(TT1+TEN*T360+SMALLR,T360)
      TT3=(TTP-TTM)/TWO
      TT3=MOD(TT3+TEN*T360+SMALLR,T360)
      END IF
C
      WRITE(6,'(A,I6,A,I6,A,3F8.3,A,F8.4)')
     & ' ',PERM2(I),' ',III,'       (',
     & TT1,TT2,TT3,')    ',-RFCORR(PERM(PERM2(I)))
C
      IF (LUNIT.GT.0) THEN
      WRITE(LUNIT,'(1X,I6,2X,3F8.3,2X,F8.4)')
     & PERM2(I),TT1,TT2,TT3,-RFCORR(PERM(PERM2(I)))
      END IF
C
      IF (RFSYMB(1:RFSYML).NE.' ') THEN
      IF (ROTMOD.EQ.'EULE'.OR.ROTMOD.EQ.'LATT') THEN
      CALL ENCODI(III,WD,WDMAX,WDLEN)
      DBPREC=-RFCORR(PERM(PERM2(I)))
      BUF=RFSYMB(1:RFSYML)//'_RFVALUE_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT1
      BUF=RFSYMB(1:RFSYML)//'_THETA1_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT2
      BUF=RFSYMB(1:RFSYML)//'_THETA2_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT3
      BUF=RFSYMB(1:RFSYML)//'_THETA3_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      CALL ENCODI(III,WD,WDMAX,WDLEN)
      DBPREC=-RFCORR(PERM(PERM2(I)))
      BUF=RFSYMB(1:RFSYML)//'_RFVALUE_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT1
      BUF=RFSYMB(1:RFSYML)//'_PSI_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT2
      BUF=RFSYMB(1:RFSYML)//'_PHI_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      DBPREC=TT3
      BUF=RFSYMB(1:RFSYML)//'_KAPPA_'//WD(1:WDLEN)
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      END IF
      END IF
C
      END IF
      END DO
C
      IF (RFSYMB(1:RFSYML).NE.' ') THEN
      DBPREC=III
      BUF=RFSYMB(1:RFSYML)//'_NLIST'
      CALL DECLAR( BUF, 'DP', ' ', DBCOMP, DBPREC )
      END IF
C
CSGI-specific optimization
C*$* OPTIMIZE(5)
      RETURN
      END
C============================================================================
      SUBROUTINE MAKTHE(CMODE,NR,MAXNR,ROTMOD,QSELF,
     &                  T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX,
     &                  DELTA,T1,T2,T3,N1DIM,N2DIM,N3DIM)
C
C Makes the list of the rotation angles T1, T2, T3
C to be searched.  In mode CMODE='SCAN' the number of angles
C to be searched is determined *without* storing them.  This
C mode is used for purposes of efficient dynamic storage allocation.
C
C Organization of the list
C    for (t2=t2min,t2max,t2delta)
C       for (t3=t3min,t3max,t3delta(t2))
C          for (t1=t1min,t1max,t1delta(t2))
C
C that is, t2 is the first level, t3 the second, and t1 the third.
C In the case of the Lattman angles, the stepsize for the second and
C third level is a function of t2, whereas for the equidistant modes
C the stepsize for the second and third level is a constant.
C
C The number of elements in the first, second, and third levels
C are returned in N2DIM, N3DIM(i2), N1DIM(i2).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      CHARACTER*4 CMODE
      INTEGER NR, MAXNR
      CHARACTER*4 ROTMOD
      LOGICAL QSELF
      DOUBLE PRECISION T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX
      DOUBLE PRECISION DELTA, T1(*), T2(*), T3(*)
      INTEGER N1DIM(*), N2DIM, N3DIM(*)
C local
      INTEGER N1, N3
      INTEGER IT1, IT2, IT3
      DOUBLE PRECISION DELTA1, DELTA2, DELTA3
      DOUBLE PRECISION THETA1, THETA2, THETA3, TEMP
C parameter
      DOUBLE PRECISION ZERO, TWO, RAD, TEN, T360
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, RAD=PI/180.0D0)
      PARAMETER (TEN=10.0D0, T360=360.0D0)
C begin
C
C check consistency of sampling intervals
      IF (T1MIN.GT.T1MAX) THEN
      TEMP=T1MIN
      T1MIN=T1MAX
      T1MAX=TEMP
      END IF
      IF (T2MIN.GT.T2MAX) THEN
      TEMP=T2MIN
      T2MIN=T2MAX
      T2MAX=TEMP
      END IF
      IF (T3MIN.GT.T3MAX) THEN
      TEMP=T3MIN
      T3MIN=T3MAX
      T3MAX=TEMP
      END IF
C
C
      IF (ROTMOD.EQ.'LATT') THEN
C
C use optimal sampling method by Lattman
C
C determine number of intervals for theta2
      IF (DELTA.GT.RSMALL) THEN
      N2DIM=INT((T2MAX-T2MIN)/DELTA+RSMALL) +1
      IF (N2DIM.GT.1) THEN
      DELTA2=(T2MAX-T2MIN)/(N2DIM-1)
      ELSE
      DELTA2=ZERO
      END IF
      ELSE
      N2DIM=1
      END IF
C
C initialize number of grid points
      NR=0
C
C loop over all theta2 grid points of rotation search
      DO IT2=0,N2DIM-1
      THETA2=T2MIN+IT2*DELTA2
C
C determine theta+ grid
      DELTA1=ABS(COS(RAD*THETA2/TWO))
      IF (DELTA1.GT.RSMALL) THEN
      DELTA1=DELTA/DELTA1
      N1=INT((T1MAX-T1MIN)/DELTA1+RSMALL) +1
      IF (N1.GT.1) THEN
      DELTA1=(T1MAX-T1MIN)/(N1-1)
      ELSE
      DELTA1=ZERO
      END IF
      ELSE
      N1=1
      DELTA1=ZERO
      END IF
C
C determine theta- grid
      DELTA3=ABS(SIN(RAD*THETA2/TWO))
      IF (DELTA3.GT.RSMALL) THEN
      DELTA3=DELTA/DELTA3
      N3=INT((T3MAX-T3MIN)/DELTA3+RSMALL) +1
      IF (N3.GT.1) THEN
      DELTA3=(T3MAX-T3MIN)/(N3-1)
      ELSE
      DELTA3=ZERO
      END IF
      ELSE
      N3=1
      DELTA3=ZERO
      END IF
C
      IF (CMODE.EQ.'STOR') THEN
      N1DIM(IT2+1)=N1
      N3DIM(IT2+1)=N3
      END IF
C
C loop through all theta- values
      DO IT3=0,N3-1
      THETA3=T3MIN+IT3*DELTA3
C
C loop through all theta+ values
      DO IT1=0,N1-1
      THETA1=T1MIN+IT1*DELTA1
C
C increment number of grid points
      IF (CMODE.EQ.'STOR'.AND.NR.GE.MAXNR) THEN
      CALL WRNDIE(-5,'MAKTHE','NR > MAXNR fatal coding error')
      ELSE
      NR=NR+1
      END IF
C
C store angles of this grid point
      IF (CMODE.EQ.'STOR') THEN
      T1(NR)=THETA1
      T2(NR)=THETA2
      T3(NR)=THETA3
      END IF
      END DO
      END DO
      END DO
C
      ELSE
C
C equidistant sampling method
C
C determine number of intervals for theta2
      IF (DELTA.GT.RSMALL) THEN
      N2DIM=INT((T2MAX-T2MIN)/DELTA+RSMALL) +1
      IF (N2DIM.GT.1) THEN
      DELTA2=(T2MAX-T2MIN)/(N2DIM-1)
      ELSE
      DELTA2=ZERO
      END IF
      ELSE
      N2DIM=1
      END IF
C
C initialize number of grid points
      NR=0
C
C loop over all grid points
      DO IT2=0,N2DIM-1
C
C determine number of intervals for theta1
      IF (DELTA.GT.RSMALL) THEN
      N1=INT((T1MAX-T1MIN)/DELTA+RSMALL) +1
      IF (N1.GT.1) THEN
      DELTA1=(T1MAX-T1MIN)/(N1-1)
      ELSE
      DELTA1=ZERO
      END IF
      ELSE
      N1=1
      END IF
C
C determine number of intervals for theta3
      IF (DELTA.GT.RSMALL) THEN
      N3=INT((T3MAX-T3MIN)/DELTA+RSMALL) +1
      IF (N3.GT.1) THEN
      DELTA3=(T3MAX-T3MIN)/(N3-1)
      ELSE
      DELTA3=ZERO
      END IF
      ELSE
      N3=1
      END IF
C
      IF (CMODE.EQ.'STOR') THEN
      N1DIM(IT2+1)=N1
      N3DIM(IT2+1)=N3
      END IF
C
      DO IT3=0,N3-1
      DO IT1=0,N1-1
C
C compute rotation angles at this grid point
      THETA1=T1MIN+IT1*DELTA1
      THETA2=T2MIN+IT2*DELTA2
      THETA3=T3MIN+IT3*DELTA3
C
C increment number of grid points
      IF (CMODE.EQ.'STOR'.AND.NR.GE.MAXNR) THEN
      CALL WRNDIE(-5,'MAKTHE','NR > MAXNR fatal coding error')
      ELSE
      NR=NR+1
      END IF
C
C compute rotation angles at this grid point
      IF (CMODE.EQ.'STOR') THEN
      T1(NR)=THETA1
      T2(NR)=THETA2
      T3(NR)=THETA3
      END IF
      END DO
      END DO
      END DO
C
      END IF
C
      RETURN
      END
C============================================================================
      SUBROUTINE ROTASY(QSELF,ROTMOD,THETA1,THETA2,THETA3,
     &            T1MIN,T1MAX,T2MIN,T2MAX,T3MIN,T3MAX)
C
C This routine returns the asymmetric unit equilvalent of the
C rotation THETA1, THETA2, THETA3.  P1
C is the Patterson function that is rotated, P2
C is the second Patterson function.  If QSELF is false we assume that P1 has
C P(1) symmetry, whereas the symmetry of P2 is specified by the
C symmetry operators XRSYMM.  If QSELF is true we assume the same symmetry
C for both P1 and P2 as specified by the symmetry operators XRSYMM.
C
C The routine should be utilized as follows:
C Input is theta1, theta2, theta3, QSELF.  Output is theta1, theta2,
C theta2.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'consta.inc'
      LOGICAL QSELF
      CHARACTER*4 ROTMOD
      DOUBLE PRECISION THETA1, THETA2, THETA3
      DOUBLE PRECISION T1MIN, T1MAX, T2MIN, T2MAX, T3MIN, T3MAX
C local
      LOGICAL COND
      INTEGER ISYM1, ISYM2, NSYM1, I, U, V, X, Y, DET
      DOUBLE PRECISION AXIS(3), ROT(3,3), SYMM(3,3), AOP(3,3), BOP(3,3)
      DOUBLE PRECISION TS1, TS2, TS3, TTS1, TTS2, TTS3, Q(4)
C parameter
      DOUBLE PRECISION ZERO, R180, R360, R720, TEN, SMALL, ONE
      PARAMETER (ZERO=0.0D0, R360=360.D0, R720=720.D0, TEN=10.D0)
      PARAMETER (R180=180.D0, SMALL=1.0D-3, ONE=1.0D0)
C begin
C asymmetric unit check
C
C compute rotation matrix
      CALL ROTMAT(ROT,THETA1,THETA2,THETA3,Q,AXIS,ROTMOD)
C
C loop over all symmetry operators of map P1
      IF (QSELF) THEN
      NSYM1=XRNSYM
      ELSE
      NSYM1=1
      END IF
      DO ISYM1=1,NSYM1
C
C only use rotations (ATB 2/8/94)
      DET=(XRSYMM(ISYM1,1,1)*XRSYMM(ISYM1,2,2)
     &     -XRSYMM(ISYM1,1,2)*XRSYMM(ISYM1,2,1))*XRSYMM(ISYM1,3,3)
     &   +(XRSYMM(ISYM1,2,1)*XRSYMM(ISYM1,3,2)
     &     -XRSYMM(ISYM1,2,2)*XRSYMM(ISYM1,3,1))*XRSYMM(ISYM1,1,3)
     &   +(XRSYMM(ISYM1,3,1)*XRSYMM(ISYM1,1,2)
     &     -XRSYMM(ISYM1,3,2)*XRSYMM(ISYM1,1,1))*XRSYMM(ISYM1,2,3)
      IF (DET.EQ.1) THEN
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      BOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      BOP(U,V)=BOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM1,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
C
C loop over all symmetry operators of map P2
      DO ISYM2=1,XRNSYM
C
C only use rotations (ATB 2/8/94)
      DET=(XRSYMM(ISYM2,1,1)*XRSYMM(ISYM2,2,2)
     &     -XRSYMM(ISYM2,1,2)*XRSYMM(ISYM2,2,1))*XRSYMM(ISYM2,3,3)
     &   +(XRSYMM(ISYM2,2,1)*XRSYMM(ISYM2,3,2)
     &     -XRSYMM(ISYM2,2,2)*XRSYMM(ISYM2,3,1))*XRSYMM(ISYM2,1,3)
     &   +(XRSYMM(ISYM2,3,1)*XRSYMM(ISYM2,1,2)
     &     -XRSYMM(ISYM2,3,2)*XRSYMM(ISYM2,1,1))*XRSYMM(ISYM2,2,3)
      IF (DET.EQ.1) THEN
C
C compute symmetry operator in real space
      DO U=1,3
      DO V=1,3
      AOP(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      AOP(U,V)=AOP(U,V)+XRINTR(U,X)*XRSYMM(ISYM2,X,Y)*XRTR(Y,V)
      END DO
      END DO
      END DO
      END DO
C
C compute symmetry related rotation matrix
C i.e.,  SYMM=AOP*ROT*BOP
      DO U=1,3
      DO V=1,3
      SYMM(U,V)=ZERO
      DO X=1,3
      DO Y=1,3
      SYMM(U,V)=SYMM(U,V)+AOP(U,X)*ROT(X,Y)*BOP(Y,V)
      END DO
      END DO
      END DO
      END DO
C
C compute angles corresponding to symmetry related rotation matrix
      CALL MATROT(SYMM,TTS1,TTS2,TTS3,Q,AXIS,ROTMOD)
C
C
      IF (ROTMOD.EQ.'LATT') THEN
C
C Lattman angles asymmetric unit check
C
C loop through redundant operations
      DO I=1,4
      TS1=TTS1
      TS2=TTS2
      TS3=TTS3
      IF (I.EQ.2.OR.I.EQ.4) THEN
C redundant operation (t+,t2,t- -> t+ + 2pi, -t2,t-)
      TS1=TS1+R360
      TS2=-TS2
      TS3=TS3
      END IF
C
      IF (I.EQ.3.OR.I.EQ.4) THEN
C redundant operation (t+,t2,t- -> t+ + 2pi, t2, t- + 2pi)
      TS1=TS1+R360
      TS2=TS2
      TS3=TS3+R360
      END IF
C
C check whether the symmetry related rotation ts1,ts2,ts3 is *less*
C than the current rotation t1,t2,t3.  The order is defined as follows
C (ts1,ts2,ts3) < (t1,t2,t3) if and only if ts2<t2 or if they are equal
C ts1<t1 or if they are equal ts3<t3.
C
C Since the t2 angle is periodic by 360 intervals and the t-,t+ angles
C are periodic by 720 intervals we want the smallest possible values
C greater than the minimum value of the specified ranges.
      TS2=TS2-INT((TS2-T2MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
      TS1=TS1-INT((TS1-T1MIN+TEN*R720)/R720+SMALL)*R720+TEN*R720
      TS3=TS3-INT((TS3-T3MIN+TEN*R720)/R720+SMALL)*R720+TEN*R720
C
C we've found a symmetry-related peak if it is inside the specified
C intervals and if it satisfies the order condition
      COND= TS1.LE.T1MAX+SMALL.AND.
     &      TS2.LE.T2MAX+SMALL.AND.
     &      TS3.LE.T3MAX+SMALL
      IF (COND) THEN
      COND=TS2.LT.THETA2-SMALL
      IF (.NOT.COND.AND.TS2.LT.THETA2+SMALL) THEN
      COND=TS3.LT.THETA3-SMALL
      IF (.NOT.COND.AND.TS3.LT.THETA3+SMALL) COND=TS1.LT.THETA1-SMALL
      END IF
      END IF
      IF (COND) THEN
C we want to replace the old theta1, theta2, theta3
C with the new values
      THETA1=TS1
      THETA2=TS2
      THETA3=TS3
      END IF
      END DO
C
      ELSE IF (ROTMOD.EQ.'EULE') THEN
C
C Euler angles asymmetric unit check
C
C loop through redundant operations
      DO I=1,2
      TS1=TTS1
      TS2=TTS2
      TS3=TTS3
      IF (I.EQ.2) THEN
C redundant operation (t1,t2,t3 -> t1 + pi, -t2,t3+pi)
      TS1=TS1+R180
      TS2=-TS2
      TS3=TS3+R180
      END IF
C
C check whether the symmetry related rotation ts1,ts2,ts3 is *less*
C than the current rotation t1,t2,t3.  The order is defined as follows
C (ts1,ts2,ts3) < (t1,t2,t3) if and only if ts2<t2 or if they are equal
C ts1<t1 or if they are equal ts3<t3.
C
C Since the t1, t2, t3  angles are periodic by 360 intervals
C we want the smallest possible values greater than the minimum
C value of the specified ranges.
      TS1=TS1-INT((TS1-T1MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
      TS2=TS2-INT((TS2-T2MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
      TS3=TS3-INT((TS3-T3MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
C
C we've found a symmetry-related peak if it is inside the specified
C intervals and if it satisfies the order condition
      COND= TS1.LE.T1MAX+SMALL.AND.
     &      TS2.LE.T2MAX+SMALL.AND.
     &      TS3.LE.T3MAX+SMALL
      IF (COND) THEN
      COND=TS2.LT.THETA2-SMALL
      IF (.NOT.COND.AND.TS2.LT.THETA2+SMALL) THEN
      COND=TS3.LT.THETA3-SMALL
      IF (.NOT.COND.AND.TS3.LT.THETA3+SMALL) COND=TS1.LT.THETA1-SMALL
      END IF
      END IF
      IF (COND) THEN
C we want to replace the old theta1, theta2, theta3
C with the new values
      THETA1=TS1
      THETA2=TS2
      THETA3=TS3
      END IF
      END DO
C
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
C
C Spherical polar angles asymmetric unit check
C
C loop through redundant operations
      DO I=1,4
      TS1=TTS1
      TS2=TTS2
      TS3=TTS3
      IF (I.EQ.2.OR.I.EQ.4) THEN
C redundant operation (psi,phi,kappa -> pi-psi, phi+pi, -kappa)
      TS1=R180-TS1
      TS2=TS2+R180
      TS3=-TS3
      END IF
C
      IF (I.EQ.3.OR.I.EQ.4) THEN
C redundant operation (psi,phi,kappa -> psi+pi, phi, -kappa)
      TS1=TS1+R180
      TS2=TS2
      TS3=-TS3
      END IF
C
C check whether the symmetry related rotation ts1,ts2,ts3 is *less*
C than the current rotation t1,t2,t3.  The order is defined as follows
C (ts1,ts2,ts3) < (t1,t2,t3) if and only if ts2<t2 or if they are equal
C ts1<t1 or if they are equal ts3<t3.
C
C Since the psi, phi, kappa angles are periodic by 360 intervals
C we want the smallest possible values greater than the minimum
C value of the specified ranges.
      TS1=TS1-INT((TS1-T1MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
      TS2=TS2-INT((TS2-T2MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
      TS3=TS3-INT((TS3-T3MIN+TEN*R360)/R360+SMALL)*R360+TEN*R360
C
C we've found a symmetry-related peak if it is inside the specified
C intervals and if it satisfies the order condition
      COND= TS1.LE.T1MAX+SMALL.AND.
     &      TS2.LE.T2MAX+SMALL.AND.
     &      TS3.LE.T3MAX+SMALL
      IF (COND) THEN
      COND=TS2.LT.THETA2-SMALL
      IF (.NOT.COND.AND.TS2.LT.THETA2+SMALL) THEN
      COND=TS3.LT.THETA3-SMALL
      IF (.NOT.COND.AND.TS3.LT.THETA3+SMALL) COND=TS1.LT.THETA1-SMALL
      END IF
      END IF
      IF (COND) THEN
C we want to replace the old theta1, theta2, theta3
C with the new values
      THETA1=TS1
      THETA2=TS2
      THETA3=TS3
      END IF
      END DO
      END IF
      END IF
      END DO
      END IF
      END DO
C
      RETURN
      END
C============================================================================
