C
C============================================================================
C
      SUBROUTINE SCLUST
C
C Perform cluster analysis and sorting of peak list
C
C Authors: Axel T. Brunger and Paul D. Adams
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
      INCLUDE 'cluster.inc'
      INCLUDE 'expression.inc'
C local
      DOUBLE PRECISION T1TEMP, T2TEMP, T3TEMP, RFTEMP
      INTEGER DUMMY
      DOUBLE PRECISION  EPS
      INTEGER LUNIT, NLIST
C
      CHARACTER*(WORD_SIZE) LFILE
      CHARACTER*4 ROTMOD
C
      DOUBLE PRECISION SECT, SECS
      LOGICAL QDEBUG, QSELF, QFORM, QCLUST
C
      CHARACTER*(WDMAX) RFSYMB
      INTEGER RFSYML
C symmetry handling
      INTEGER TSNSYM,XROLD
C parameters
      DOUBLE PRECISION ZERO, TWO, THREE, T1P4
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, THREE=3.0D0, T1P4=1.4D0)
C begin
C
C default values
      TSNSYM=1
      LFILE=' '
      LUNIT=0
      ROTMOD='EULE'
      NLIST=100
      EPS=0.2D0
      QDEBUG=.FALSE.
      QSELF=.FALSE.
      QFORM=.TRUE.
      RFSYMB=' '
      RFSYML=1
      QCLUST=.FALSE.
C
C allocate space for the peak list
      IF (PEAKALLOC.EQ.0) THEN
         PEAKALLOC=PEAKALLOC+100
         HPRF=ALLHP(IREAL8(PEAKALLOC))
         HPT1=ALLHP(IREAL8(PEAKALLOC))
         HPT2=ALLHP(IREAL8(PEAKALLOC))
         HPT3=ALLHP(IREAL8(PEAKALLOC))
      END IF
C
C parsing
      CALL PUSEND('CLUSter>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('CLUSter>')
      CALL MISCOM('CLUSter>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-search-cluster')
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
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODE= ',ROTMOD)
      IF (ROTMOD.NE.'LATT'.AND.
     &    ROTMOD.NE.'EULE'.AND.
     &    ROTMOD.NE.'SPHE') THEN
      WRITE(6,'(A)')'SCLUST: Unknown rotation mode, using EULEr'
      ROTMOD='EULE'
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' -----------cluster analys',
     & 'is-parameters---------------------------------------'
      IF (ROTMOD.EQ.'LATT') THEN
      WRITE(6,'(A)')
     & ' | The quasi-orthogonal Eulerian angle method is used'
      ELSE IF (ROTMOD.EQ.'EULE') THEN
      WRITE(6,'(A)')
     & ' | The equidistant Eulerian angle method is used'
      ELSE IF (ROTMOD.EQ.'SPHE') THEN
      WRITE(6,'(A)')
     & ' | The spherical polar angle method is used'
      END IF
      WRITE(6,'(A,I8)')
     & ' | NLISt=',NLIST
      WRITE(6,'(A,F5.2)')
     & ' | EPSIlon=',EPS
      WRITE(6,'(A,L1)')
     & ' | CLUSter=',QCLUST
      WRITE(6,'(A)')
     & ' | Symmetry operators for cluster analysis'
      XROLD=XRNSYM
      XRNSYM=TSNSYM
      CALL XSYPRI(XRMSYM,XRNSYM,XRSYMM,XRITSY,XRSYTH)
      XRNSYM=XROLD
      WRITE(6,'(2A)') ' ---------------------------',
     & '----------------------------------------------------'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PEAK') THEN
      CALL NEXTF('T1=',T1TEMP)
      CALL NEXTF('T2=',T2TEMP)
      CALL NEXTF('T3=',T3TEMP)
      CALL NEXTF('RF=',RFTEMP)
      NPEAK=NPEAK+1
      IF (NPEAK.GT.PEAKALLOC) THEN
      HPT1=REAHP(HPT1,IREAL8(PEAKALLOC),IREAL8(PEAKALLOC+100))
      HPT2=REAHP(HPT2,IREAL8(PEAKALLOC),IREAL8(PEAKALLOC+100))
      HPT3=REAHP(HPT3,IREAL8(PEAKALLOC),IREAL8(PEAKALLOC+100))
      HPRF=REAHP(HPRF,IREAL8(PEAKALLOC),IREAL8(PEAKALLOC+100))
      PEAKALLOC=PEAKALLOC+100
      END IF
      CALL ADDPEAK(T1TEMP,T2TEMP,T3TEMP,RFTEMP,NPEAK,
     &             HEAP(HPT1),HEAP(HPT2),HEAP(HPT3),HEAP(HPRF))
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CLUS') THEN
      CALL NEXTLO('CLUSter=',QCLUST)
C=====================================================================
      ELSE
      CALL CHKEND('CLUSter>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (XRNSYM.NE.1) THEN
      CALL WRNDIE(-5,'CLUSter','P1 spacegroup must be selected.')
      ELSE
C
C assign list output file if required
      IF (LFILE.NE.' ') THEN
      CALL ASSFIL(LFILE,LUNIT,'WRITE','FORMATTED',ERROR)
      IF (ERROR) GOTO 9999
      END IF
C
C switch on symmetry operators for cluster analysis
      XROLD=XRNSYM
      XRNSYM=TSNSYM
C
C number of highest peaks to be analysed can't be larger than NPEAK
      NLIST=MIN(NLIST,NPEAK)
C
C make a cluster analysis of the NLIST highest peaks and
C mark redundant peaks
      IF (QCLUST) THEN
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECT)
      END IF
      CALL RFCLUS(LUNIT,EPS,
     & NPEAK,HEAP(HPT1),HEAP(HPT2),HEAP(HPT3),HEAP(HPRF),
     & NLIST,ROTMOD,QDEBUG,QSELF,
     & ZERO,ZERO,ZERO,ZERO,ZERO,ZERO,
     &             RFSYMB,RFSYML)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      SECT=SECS-SECT
      WRITE(6,'(A,F10.4)')
     & ' DRSEAR: CPU-Time:  cluster analysis=',SECT
      END IF
      END IF
C
C restore original symmetry (should be identity operator)
      XRNSYM=XROLD
C
      END IF
C
9999  CONTINUE
      RETURN
      END
C
C ===========================================================================
C
      SUBROUTINE ADDPEAK(PEAKT1,PEAKT2,PEAKT3,PEAKRF,NPEAK,T1,T2,T3,RF)
C
      IMPLICIT NONE
C
      INTEGER NPEAK
      DOUBLE PRECISION PEAKT1, PEAKT2, PEAKT3, PEAKRF
      DOUBLE PRECISION T1(*), T2(*), T3(*), RF(*)
C
      T1(NPEAK)=PEAKT1
      T2(NPEAK)=PEAKT2
      T3(NPEAK)=PEAKT3
C
C make peak height negative for later sorting
      RF(NPEAK)=-PEAKRF
C
      RETURN
      END
C
C ===========================================================================
C
      SUBROUTINE CLUSINI
C
      IMPLICIT NONE
C
      INCLUDE 'cluster.inc'
C
C initial space for peak list
      PEAKALLOC=0
C
C initial number of peaks
      NPEAK=0
C
      RETURN
      END
C
C ===========================================================================
C
      SUBROUTINE CLUSFRE
C
      IMPLICIT NONE
C
      INCLUDE 'cluster.inc'
      INCLUDE 'funct.inc'
C
C free up space for rotation function values
      IF (PEAKALLOC.GT.0) THEN
         CALL FREHP(HPT3,IREAL8(PEAKALLOC))
         CALL FREHP(HPT2,IREAL8(PEAKALLOC))
         CALL FREHP(HPT1,IREAL8(PEAKALLOC))
         CALL FREHP(HPRF,IREAL8(PEAKALLOC))
      END IF
C
      RETURN
      END
C
C ===========================================================================
C
