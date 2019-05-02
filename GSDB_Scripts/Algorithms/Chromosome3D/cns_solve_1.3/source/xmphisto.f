      SUBROUTINE XMHISTO(NA,NB,NC,NAP,NBP,NCPL,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &  HPRHOMA,NRHO,IRHO,NMASK,XRMAP)
C
C Computes electron density map histogram and histogram matching.
C
C References:
C
C Zhang, K.Y.J., and Main, P. (1990), Acta Cryst. A46, 41-46.
C "Histogram Matching as a New Density Modification Technique
C for Phase Refinement and Extension of Protein Molecules".
C
C Main, P. (1990), Acta Cryst. A46, 507-509.
C "A Formula for Electron Density Histograms for Equal-Atom Structures".
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'ncs.inc'
      INTEGER NA, NB, NC, NAP, NBP, NCPL, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      INTEGER ADEPTH, ARPNMX, ARPNN, ARPNX
      CHARACTER*(*) ARPN(4,*)
      INTEGER ARPNL(4,*)
      DOUBLE COMPLEX ARPNDB(4,*)
      INTEGER ARPNMLT(*)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
      LOGICAL XRMAP
C local
C pointer
      INTEGER FLAGS, XRVDW, XL, YL, ZL, WDATA
C begin
      FLAGS=ALLHP(INTEG4(NATOM))
      XRVDW=ALLHP(IREAL8(NATOM*XNNSYM))
      XL=ALLHP(IREAL8(NATOM*XNNSYM))
      YL=ALLHP(IREAL8(NATOM*XNNSYM))
      ZL=ALLHP(IREAL8(NATOM*XNNSYM))
      WDATA=ALLHP(IREAL8(10000))
C
      CALL XMHIST2(NATOM,HEAP(FLAGS),HEAP(XRVDW),HEAP(XL),
     &   HEAP(YL),HEAP(ZL),MAXCN,CNBVR,LOOKUP,NA,NB,NC,NAP,NBP,
     &   NCPL,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,HEAP(WDATA))
C
      CALL FREHP(WDATA,IREAL8(10000))
      CALL FREHP(ZL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(YL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(XL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(XRVDW,IREAL8(NATOM*XNNSYM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHIST2(NATOM,FLAGS,XRVDW,XL,
     &   YL,ZL,MAXCN,CNBVR,LOOKUP,NA,NB,NC,NAP,NBP,
     &   NCP,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,WDATA)
C
C See the above
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'timer.inc'
      INTEGER FLAGS(*), NATOM
      DOUBLE PRECISION XRVDW(*), XL(*), YL(*), ZL(*)
      INTEGER MAXCN
      DOUBLE PRECISION CNBVR(MAXCN,MAXCN)
      INTEGER LOOKUP(*)
      INTEGER NA, NB, NC, NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRCELL(9)
      INTEGER ADEPTH, ARPNMX, ARPNN, ARPNX
      CHARACTER*(*) ARPN(4,*)
      INTEGER ARPNL(4,*)
      DOUBLE COMPLEX ARPNDB(4,*)
      INTEGER ARPNMLT(*)
      DOUBLE PRECISION MAPR
      INTEGER XRHONUM
      CHARACTER*(*) XRHONAM(*)
      INTEGER HPRRHO(*), HPIRHO(*)
      INTEGER HPRHOMA, IRHO, NRHO, NMASK
      LOGICAL XRMAP
      DOUBLE PRECISION WDATA(*)
C
C local
      DOUBLE PRECISION DSUM, DAVE, DMIN, DMAX, DRMS, DNORM
      DOUBLE PRECISION DSLOT, DMINC, DMAXC
      DOUBLE PRECISION DSLOTE, DAVEE, DRMSE, HSCALE
      DOUBLE PRECISION LOWFRAC, CUTOFF
      INTEGER MBINS, MBINSE, I
      CHARACTER*4 MODEL
      INTEGER IFROM
      DOUBLE PRECISION SECS, SECS2
      LOGICAL COND, ERR
      LOGICAL QMBINS, QDSLOT, QDMIN, QDMAX, QDATA, QMATCH, QHRESO
      INTEGER EXUNIT, EXLEN
      CHARACTER*(WORD_SIZE) EXFILE
      DOUBLE PRECISION HRESO
C pointers
      INTEGER HPRMAP, HPIMAP
      INTEGER MPACK
      INTEGER NHIST, HIST, HISTE, HSLOT, HSCAL, HSHIF
C number of grid points
      INTEGER NNSELE, NNHIST, NPROT
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, SIX, R100, HALF, DSMALL
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, SIX=6.0D0)
      PARAMETER (R100=100.D0, HALF=0.5D0, DSMALL=1.0D-4)
C
C begin
C default values
      ERR=.FALSE.
      IFROM=1
      MODEL=' '
      LOWFRAC=ZERO
      MBINS=100
      QMBINS=.FALSE.
      QDSLOT=.FALSE.
      QDMIN=.FALSE.
      QDMAX=.FALSE.
      QDATA=.FALSE.
      QMATCH=.FALSE.
      QHRESO=.FALSE.
      EXUNIT=0
      EXFILE=' '
      EXLEN=0
      HRESO=ZERO
      DSLOT=ZERO
      DSLOTE=ZERO
C
C make/heap default map selection (MAP1)
      IF (HPRHOMA.EQ.0) THEN
      CALL WRNDIE(-5,'XMHIST2','Map is undefined.')
      ELSE
      MPACK=ALLHP(INTEG4(NRHO))
      CALL XMPPCK0(HEAP(MPACK),NNSELE,
     & MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,HEAP(HPRHOMA))
C
C parsing
      CALL PUSEND('HISTO>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('HISTO>')
      CALL MISCOM('HISTO>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-histogram')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODEl=',MODEL)
      IF (MODEL.NE.'PRED'.AND.MODEL.NE.'EXPE'
     &   .AND.MODEL.NE.'pred'.AND.MODEL.NE.'expe') THEN
      CALL DSPERR('HISTO','unknown model for MODEL')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'FROM=',XRHONAM(IFROM)
      ELSE
C
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      IFROM=I
      END IF
      END DO
C
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XMHISTO-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMHISTO','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE'.OR.WD(1:1).EQ.'(') THEN
      CALL XMPSELE(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &      QHERM,XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &      HPRHOMA,NRHO,NMASK,
     &      XRNSYM,HEAP(MPACK),NNSELE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SOLC'.OR.WD(1:4).EQ.'LOWF') THEN
      CALL NEXTF('LOWFraction=',LOWFRAC)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MBIN') THEN
      IF (QMATCH) THEN
      CALL NEXTI('MBINS=',MBINSE)
      IF (MBINS.NE.MBINSE) MBINS=MBINSE
      ELSE
      CALL NEXTI('MBINS=',MBINS)
      END IF
      IF (MBINS.GT.0.AND.MBINS.LE.10000) QMBINS=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SLOT') THEN
      IF (QMATCH) THEN
      CALL NEXTF('SLOT=',DSLOTE)
      IF (DSLOTE.GT.ZERO) QDSLOT=.TRUE.
      ELSE
      CALL NEXTF('SLOT=',DSLOT)
      IF (DSLOT.GT.ZERO) QDSLOT=.TRUE.
      END IF
C=====================================================================
      ELSE IF (WD(1:6).EQ.'RHOMIN') THEN
      CALL NEXTF('RHOMIN cutoff=',DMINC)
      QDMIN=.TRUE.
C=====================================================================
      ELSE IF (WD(1:6).EQ.'RHOMAX') THEN
      CALL NEXTF('RHOMAX cutoff=',DMAXC)
      QDMAX=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'HRES') THEN
      CALL NEXTF('HRESolution=',HRESO)
      IF (HRESO.GT.DSMALL) QHRESO=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'EXPO'.OR.WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('EXPOrt filename=',EXFILE)
      IF (EXFILE.NE.' ') THEN
      EXLEN=LEN(EXFILE)
      CALL TRIMM(EXFILE,EXLEN)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DATA') THEN
      IF (QMBINS.AND.MBINS.LE.10000) THEN
      I=1
      DO WHILE (I.LE.MBINS)
      CALL NEXTF('DATA=',WDATA(I))
      I=I+1
      END DO
      QDATA=.TRUE.
      ELSE
      WRITE(6,'(A,A)') ' %XMHISTO-ERR:',
     & ' MBINS is not specified for reading data!'
      QDATA=.FALSE.
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MATC') THEN
      QMATCH=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     & ' histo ---------------------------------------'
      WRITE(6,'(A,4A,I4)') ' | FROM = ',
     & XRHONAM(IFROM),'  MODEl = ',MODEL, '   MBINS = ',MBINS
      IF (QDSLOT.AND.QMATCH) THEN
      WRITE(6,'(A,F10.6)')
     & ' | SLOT width = ', DSLOTE
      ELSE IF (QDSLOT.AND..NOT.QMATCH) THEN
      WRITE(6,'(A,F10.6)')
     & ' | SLOT width = ', DSLOT
      END IF
      IF (QDMIN) WRITE(6,'(A,F8.4)')
     & ' | RHOMIN low density bound = ', DMINC
      IF (QDMAX) WRITE(6,'(A,F8.4)')
     & ' | RHOMAX high density bound = ', DMAXC
      IF (LOWFRAC.GT.ZERO) WRITE(6,'(A,F7.3)')
     & ' | LOWFraction = ', LOWFRAC
      IF (HRESO.GT.ZERO) WRITE(6,'(A,F7.3)')
     & ' | HRESolution = ', HRESO
      IF (EXFILE.NE.' ') WRITE(6,'(A,A)')
     & ' | EXPort file=',EXFILE(1:EXLEN)
      WRITE(6,'(2A)') ' |--------------------------',
     & '--------------------------------------------'
      ELSE
      CALL CHKEND('HISTOgram>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (QMATCH) WRITE(6,'(A,A,A)') ' XMHISTO: CAUTION: overwrite to ',
     & XRHONAM(IFROM),' after histogram matching.'
C
C assign export file if required
      IF (EXFILE.NE.' ') THEN
      CALL ASSFIL(EXFILE,EXUNIT,'WRITE','FORMATTED',ERROR)
      END IF
C
C initial the timer
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C===== SECTION 1 =========================================
C===== Get the Heap Allocation for the Map.
C
      IF (.NOT.ERR.AND.HPRRHO(IFROM).EQ.0) THEN
      WRITE(6,'(3A)') ' %XMHISTO-ERR: real space object ',
     & XRHONAM(IFROM),' undefined.'
      CALL WRNDIE(-5,'XMHISTO','object undefined.')
      ELSEIF (.NOT.ERR) THEN
      HPRMAP=HPRRHO(IFROM)
      HPIMAP=HPIRHO(IFROM)
C
      IF (HPRMAP.NE.0.AND.HPRHOMA.NE.0) THEN
C
C
C===== SECTION 2 =========================================
C===== Statistical Information on the Map
C
C initialize statistical properties
      CALL XMPSTIN(DSUM,DAVE,DMIN,DMAX,DNORM,DRMS)
C
C evaluate statistical properties
      CALL XMPSTAT(HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   DSUM,DMIN,DMAX,DNORM)
C
C write out and declare variables
      CALL XMPSTOU(NNSELE,DSUM,DAVE,DMIN,DMAX,DNORM,DRMS)
C
C
C===== SECTION 3 =========================================
C===== Histogram Settings
C
C defaults from the map (trim the accuracy within DSMALL)
      DMAX=FLOAT(NINT(DMAX/DSMALL))*DSMALL
      DMIN=FLOAT(NINT(DMIN/DSMALL))*DSMALL
      WRITE(6,'(A,A,2F9.4)') ' XMHISTO: (default from map)',
     & ' RHOMIN and RHOMAX =',DMIN,DMAX
      IF ((DMAX-DSMALL).LT.DMIN) THEN
      ERR=.TRUE.
      WRITE(6,'(A)') ' %XMHISTO-ERR:  a complete flat map.'
C      CALL WRNDIE(-5,'XMHISTO',' a complete flat map.')
      END IF
C histogram slot width
      IF (QMATCH.OR..NOT.QDSLOT) THEN
      DSLOT=(DMAX-DMIN)/MBINS
      IF (.NOT.(LOWFRAC.GT.ZERO))
     & DSLOT=FLOAT(NINT(DSLOT/DSMALL))*DSMALL
      WRITE(6,'(A,A,F10.6)') ' XMHISTO: (default from map)',
     & ' SLOT width =',DSLOT
      END IF
C minimum density cutoff
      IF (QDMIN) THEN
      IF (QMATCH) THEN
      WRITE(6,'(A,A,A,F9.4)') ' XMHISTO:',
     & ' histogram matching starts at minimum density',
     & ' RHOMIN=',DMINC
      DMIN=DMINC
      ELSE IF (ABS(DMINC-DMIN).GT.DSMALL) THEN
      WRITE(6,'(A,A,F9.4)') ' XMHISTO: minimum density bound',
     & ' set for histogram RHOMIN=',DMINC
      DMIN=DMINC
      END IF
      END IF
C maximum density cutoff
      IF (QDMAX) THEN
      IF (QMATCH) THEN
      WRITE(6,'(A,A,A,F9.4)') ' XMHISTO:',
     & ' histogram matching ends   at maximum density',
     & ' RHOMAX=',DMAXC
      DMAX=DMAXC
      ELSE IF (ABS(DMAXC-DMAX).GT.DSMALL) THEN
      WRITE(6,'(A,A,F9.4)') ' XMHISTO: maximum density bound',
     & ' set for histogram RHOMAX=',DMAXC
      DMAX=DMAXC
      END IF
      END IF
C modify histogram SLOT and MBINS
      IF (QMATCH.AND.DSLOTE.NE.DSLOT) THEN
      WRITE(6,'(A,A,F10.6)') ' XMHISTO: histogram matching',
     & ' must use the same width SLOT=',DSLOTE
      DSLOT=DSLOTE
      END IF
      IF (DSLOT.LE.DSMALL*DSMALL) THEN
      ERR=.TRUE.
      WRITE(6,'(A,A,F9.5)') ' %XMHISTO-ERR: the SLOT width',
     & ' is too small ',DSLOT
C      CALL WRNDIE(-5,'XMHISTO',' SLOT is too small.')
      END IF
      IF (.NOT.ERR.AND..NOT.QMATCH) THEN
      MBINSE=INT((DMAX-DMIN)/DSLOT+HALF)
      IF (MBINSE.NE.MBINS) THEN
      WRITE(6,'(A,A)') ' XMHISTO: MBINS is modified',
     & ' according to actual histogram setting'
      MBINS=MBINSE
      END IF
      END IF
      WRITE(6,'(A,I5,A,F9.4)')
     & ' XMHISTO: the number of slots MBINS=',MBINS,
     & ' and width SLOT=',DSLOT
      IF (MBINS.GT.10000) THEN
      ERR=.TRUE.
      WRITE(6,'(A)') '%XMHISTO-ERR: MBINS exceeds 10000.'
C      CALL WRNDIE(-5,'XMHISTO',' MBINS exceeds 10000.')
      END IF
C check resolution specification
      IF (QMATCH.AND..NOT.QHRESO)
     & WRITE(6,'(A)') ' %XMHISTO-WRN: resolution is not specified!'
C
C
      IF (.NOT.ERR) THEN
C
C allocate heap for histogram slots
      NHIST=ALLHP(INTEG4(MBINS+1))
      HIST=ALLHP(IREAL8(MBINS+1))
      HISTE=ALLHP(IREAL8(MBINS+1))
      HSLOT=ALLHP(IREAL8(MBINS+1))
      HSCAL=ALLHP(IREAL8(MBINS+1))
      HSHIF=ALLHP(IREAL8(MBINS+1))
C
C===== SECTION 4 =========================================
C===== Evaluate Histogram from Map
C
C initialize histogram slots
      CALL XMHISIN(MBINS,HEAP(NHIST),HEAP(HIST))
C
C accumulate histogram
      CALL XMHISTA(MBINS,HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   DMIN,DSLOT,HEAP(NHIST))
C
C evaluate histogram in percentage
      CALL XMHISTC(1,MBINS,HEAP(NHIST),HEAP(HIST),NNSELE)
C
C===== SECTION 5 =========================================
C===== Histogram Matching
C
C copy experimental histogram from the buffer data array
      IF (QMATCH.AND.QDATA) THEN
      CALL COPYR8(WDATA,HEAP(HISTE),MBINS)
      END IF
C
C compute the expected histogram according to Main's formula
      IF (MODEL.EQ.'PRED') THEN
      WRITE(6,'(A)') ' XMHISTO: computing expected histogram....'
      CALL XMHTHEO(MBINS,HEAP(HISTE),HRESO,DSLOT,DMIN)
      IF (.NOT.QMATCH) THEN
      CALL COPYR8(HEAP(HISTE),HEAP(HIST),MBINS)
      CALL XMHISTC(-1,MBINS,HEAP(NHIST),HEAP(HISTE),NNSELE)
      END IF
      END IF
C
C estimate the mean and r.m.s. of the expected map
      IF (QMATCH.OR.MODEL.EQ.'PRED') THEN
      CALL XMHISTE(MBINS,HEAP(HISTE),DSLOT,DMIN,DAVEE,DRMSE)
      WRITE(6,'(A,A,F8.4)') ' XMHISTO:',
     & ' the estimated mean of the expected map =',DAVEE
      WRITE(6,'(A,A,F8.4)') ' XMHISTO:',
     & ' the estimated r.m.s. of the expected map =',DRMSE
C the scale factor (on r.m.s.) between expected and observed maps
      HSCALE=DRMSE/DRMS
      WRITE(6,'(A,A,F8.5)') ' XMHISTO:',
     & ' the estimated scale (on r.m.s.) between two maps =',HSCALE
      ELSE
      HSCALE=ZERO
      END IF
C
C apply the scale factor to the observed map
      IF (QMATCH.AND.ABS(HSCALE-ONE).GT.ONE/R100) THEN
      WRITE(6,'(A,A)') ' XMHISTO: apply the overall scale factor',
     & ' to the observed map'
      CALL XMHSCAL(MBINS,HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   DAVE,HSCALE)
C and update histogram
      CALL XMHISIN(MBINS,HEAP(NHIST),HEAP(HIST))
      CALL XMHISTA(MBINS,HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   DMIN,DSLOT,HEAP(NHIST))
      CALL XMHISTC(1,MBINS,HEAP(NHIST),HEAP(HIST),NNSELE)
      END IF
C
C evaluate the transformed solts, individual scale factors and shifts
      IF (QMATCH) THEN
      CALL XMHMATC(MBINS,MBINS,HEAP(HIST),HEAP(HISTE),HEAP(HSLOT),
     &             HEAP(HSCAL),HEAP(HSHIF),DSLOT,DMIN)
      END IF
C
C apply the individual scale factors and shifts to the map
      IF (QMATCH) THEN
      CALL XMHTRAN(MBINS,HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   HEAP(HSCAL),HEAP(HSHIF),DSLOT,DMIN)
C and update histogram
      CALL XMHISIN(MBINS,HEAP(NHIST),HEAP(HIST))
      CALL XMHISTA(MBINS,HEAP(MPACK),
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   DMIN,DSLOT,HEAP(NHIST))
      CALL XMHISTC(1,MBINS,HEAP(NHIST),HEAP(HIST),NNSELE)
      END IF
C
C===== SECTION 6 =========================================
C===== Print Histogram
      CALL XMHISTP(MBINS,HEAP(NHIST),HEAP(HIST),
     &             DMIN,DSLOT,NNHIST,NNSELE)
C
C===== SECTION 7 =========================================
C===== Determine the Cutoff for the Boundary of Protein/Solvent
C
      IF (LOWFRAC.GT.ZERO) THEN
      CALL XMHICUT(MBINS,HEAP(NHIST),NNSELE,DMIN,DSLOT,
     &                   LOWFRAC,CUTOFF,NPROT)
      END IF
C
C===== SECTION 8 =========================================
C===== Export Histogram
C
      IF (EXFILE.NE.' ') THEN
      IF (.NOT.QHRESO)
     & WRITE(6,'(A)') ' %XMHIST-WRN: resolution is not specified!'
      CALL XMHEXPO(EXUNIT,EXFILE,EXLEN,HRESO,
     &             MBINS,DSLOT,DMIN,DMAX,HEAP(HIST))
      IF (EXUNIT.NE.0) CALL VCLOSE(EXUNIT,'KEEP',ERROR)
      END IF
C
C=========================
C deallocate space for heaps
      CALL FREHP(HSHIF,IREAL8(MBINS+1))
      CALL FREHP(HSCAL,IREAL8(MBINS+1))
      CALL FREHP(HSLOT,IREAL8(MBINS+1))
      CALL FREHP(HISTE,IREAL8(MBINS+1))
      CALL FREHP(HIST,IREAL8(MBINS+1))
      CALL FREHP(NHIST,INTEG4(MBINS+1))
C
      END IF
C=========================
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' HISTO2: CPU-time:  histogram = ',SECS2-SECS
      END IF
C
C
C
      END IF
      END IF
C
C free heap
      IF (HPRHOMA.NE.0) CALL FREHP(MPACK,INTEG4(NRHO))
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMPSTIN(DSUM,DAVE,DMIN,DMAX,DNORM,DRMS)
C
C Initializes statistical properties.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION DSUM,DAVE,DMIN,DMAX,DNORM,DRMS
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      DSUM=ZERO
      DAVE=ZERO
      DMIN=R4BIG
      DMAX=-R4BIG
      DNORM=ZERO
      DRMS=ZERO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMPSTAT(MPACK,
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   RRHO,IRHO,
     &   DSUM,DMIN,DMAX,DNORM)
C
C Computes statistical properties on the selected set
C (N map elements from START to STOP)
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION DSUM, DMIN, DMAX, DNORM
C local
      INTEGER A, B, C
      DOUBLE PRECISION DENS
C begin
C
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (MPACK(A,B,C).EQ.1) THEN
      DENS=RRHO(A,B,C)
      DSUM=DSUM+DENS
      DMIN=MIN(DMIN,DENS)
      DMAX=MAX(DMAX,DENS)
      DNORM=DNORM+DENS*DENS
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMPSTOU(N,DSUM,DAVE,DMIN,DMAX,DNORM,DRMS)
C
C Writes out statistical properties.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER N
      DOUBLE PRECISION DSUM, DAVE, DMIN, DMAX, DNORM, DRMS
C local
      DOUBLE PRECISION DBPREC
C begin
C
      IF (N.NE.0) THEN
C
C number of elements
C average minimum maximum
      DAVE=DSUM/N
      WRITE(6,'(A,F8.4,A,F8.4,A,F8.4)')
     & ' XMPST: average = ',DAVE,
     &       '  minimum = ',DMIN,
     &       '  maximum = ',DMAX
C norm and standard deviation (r.m.s.) or sigma
      DRMS=SQRT(DNORM/N-DAVE*DAVE)
      DBPREC=SQRT(DNORM/N)
      WRITE(6,'(A,F8.4,A,F8.4)')
     & ' XMPST: r.m.s.  = ',DRMS,'     norm = ',DBPREC
C
      ELSE
C
      CALL WRNDIE(-5,'XMPST','Number of selected elements is zero.')
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHISIN(MBINS,NHIST,HIST)
C
C Initialize histogram slots.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS, NHIST(*)
      DOUBLE PRECISION HIST(*)
C local
      INTEGER I
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      DO I=1,MBINS+1
      NHIST(I)=0
      HIST(I)=ZERO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHISTA(MBINS,MPACK,
     &   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &   RRHO,IRHO,DMIN,DSLOT,NHIST)
C
C Projects density into histogram slots.
C (operates on N map elements from START to STOP)
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION DMIN, DSLOT
      INTEGER NHIST(*)
C local
      INTEGER A, B, C, J
      DOUBLE PRECISION DENS
C parameter
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
C begin
C
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (MPACK(A,B,C).EQ.1) THEN
      DENS=RRHO(A,B,C)
      IF (DENS.GE.DMIN) THEN
      J=MAX(1,NINT((DENS-DMIN)/DSLOT))
      IF (J.LE.MBINS) NHIST(J)=NHIST(J)+1
      END IF
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHISTC(MODE,MBINS,NHIST,HIST,NPTS)
C
C Converts histogram (always in percentage).
C MODE=1, from points to percentages
C MODE=-1, from percentages to points
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MODE, MBINS, NHIST(*)
      DOUBLE PRECISION HIST(*)
      INTEGER NPTS
C local
      INTEGER I
C parameters
      DOUBLE PRECISION HUNDRED
      PARAMETER (HUNDRED=100.0D0)
C begin
C
      IF (MODE.EQ.1) THEN
      DO I=1,MBINS
      HIST(I)=NHIST(I)*HUNDRED/NPTS
      END DO
      ELSE IF (MODE.EQ.-1) THEN
      DO I=1,MBINS
      NHIST(I)=HIST(I)/HUNDRED*NPTS
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHISTP(MBINS,NHIST,HIST,DMIN,DSLOT,NNHIST,NPTS)
C
C Prints histogram.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      INTEGER MBINS, NHIST(*)
      DOUBLE PRECISION HIST(*)
      DOUBLE PRECISION DMIN, DSLOT
      INTEGER NNHIST, NPTS
C local
      INTEGER I, ANHIST, GISNA
      DOUBLE PRECISION DENS, DENS1, GISPA
      CHARACTER*30 LINE
      INTEGER LLINE
C parameters
      DOUBLE PRECISION ZERO, ONE, HUNDRED
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HUNDRED=100.0D0)
C begin
C
C total elements
      NNHIST=0
      DO I=1,MBINS
      NNHIST=NNHIST+NHIST(I)
      END DO
      WRITE(6,'(A,I9)')
     & ' XMHISTP: total grid-points within histogram slots',NNHIST
      ANHIST=NPTS-NNHIST
      IF (ANHIST.NE.0) WRITE(6,'(A,I6)')
     & ' XMHISTP: grid-points outside the histogram',ANHIST
C
      IF (.NOT.(PUNIT.EQ.6.AND.WRNLEV.EQ.0)) THEN
C
C print headings
      LINE=' Density Histogram '
      LLINE=19
      WRITE(PUNIT,'(3A)')
     &  ' {* ================== ',LINE(1:LLINE),
     &  ' ================== *}'
      WRITE(PUNIT,'(A,A)')
     &  ' {* density-range     points   percent',
     &  '    acc. pts  acc. pcnt *}'
C
C print histogram
      DO I=1,MBINS
      DENS=DMIN+(I-1)*DSLOT
      DENS1=DENS+DSLOT
      ANHIST=ANHIST+NHIST(I)
      GISNA=NPTS-ANHIST+NHIST(I)
      GISPA=GISNA*HUNDRED/NPTS
      WRITE(PUNIT,'(A,2F8.3,I10,F10.3,I12,F10.3)')
     & '  ',DENS,DENS1,NHIST(I),HIST(I),
     & GISNA,GISPA
      END DO
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHICUT(MBINS,NHIST,NALL,DMIN,DSLOT,
     &                   LOWFRAC,CUTOFF,NPROT)
C
C Determines the cutoff such that LOWFRAC % of the
C selected density points is below the cutoff.
C The cutoff is stored in the symbol $CUTOFF
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS, NALL
      INTEGER NHIST(*)
      DOUBLE PRECISION DMIN, DSLOT
      DOUBLE PRECISION LOWFRAC, CUTOFF
      INTEGER NPROT
C local
      INTEGER I, N
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBPREC
C parameter
      DOUBLE PRECISION HALF, ONE, HUNDRED
      PARAMETER (HALF=0.5D0,ONE=1.0D0,HUNDRED=100.0D0)
C begin
C
      NPROT=INT((ONE-LOWFRAC)*NALL+HALF)
C
      N=0
      I=MBINS+1
      DO WHILE (I.GT.1.AND.N.LT.NPROT)
      I=I-1
      N=N+NHIST(I)
      END DO
      CUTOFF=DMIN+I*DSLOT
      NPROT=N
      WRITE(6,'(A,F8.3,A,F10.5)')
     & ' XMHICUT: ',LOWFRAC*HUNDRED,
     & '% of the selected density points are below cutoff=',
     & CUTOFF
      CALL DECLAR( 'RESULT','DP',' ',DBCOMP,CUTOFF)
      CALL DECLAR( 'CUTOFF','DP',' ',DBCOMP,CUTOFF)
      DBPREC=NPROT
      CALL DECLAR( 'NPROT','DP',' ',DBCOMP,DBPREC)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHEXPO(EXUNIT,EXFILE,EXLEN,HRESO,
     &                   MBINS,DSLOT,DMIN,DMAX,HIST)
C
C Exports histogram (always in percentage).
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INTEGER EXUNIT, EXLEN
      CHARACTER*(*) EXFILE
      INTEGER MBINS
      DOUBLE PRECISION HIST(*)
      DOUBLE PRECISION HRESO, DSLOT, DMIN, DMAX
C local
      INTEGER I
C parameters
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C begin
C
C headings
      WRITE(EXUNIT,'(A,A)') ' REMARK  histogram data file: ',
     & EXFILE(1:EXLEN)
      WRITE(EXUNIT,'(A,F9.4,A)') '  hreso=',HRESO,
     & '          {* histogram at resolution *}'
      WRITE(EXUNIT,'(A,I5,A)') '  mbins=',MBINS,
     & '              {* the number of slots *}'
      WRITE(EXUNIT,'(A,F9.4,A)') '  slot=',DSLOT,
     & '           {* the width of the slot *}'
      WRITE(EXUNIT,'(A,F9.4,A)') '  rhomin=',DMIN,
     & '         {* the minimum density bound *}'
      WRITE(EXUNIT,'(A,F9.4,A)') '  rhomax=',DMAX,
     & '         {* the maximum density bound *}'
C
C data block
      WRITE(EXUNIT,'(A,A,I4,A)') '  data=    {%}    ',
     & '         {* the total number of data =',MBINS,' (MBINS) *}'
      DO I=1,MBINS
      WRITE(EXUNIT,'(7X,F10.6)') HIST(I)
      END DO
C
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHMATC(MBINS,MBINSE,HIST,HISTE,HSLOT,
     &                   HSCAL,HSHIF,DSLOT,DMIN)
C
C Matches the observed histogram to the expected histogram.
C The transformed slots should satisfy the following condition
C
C   Integral [P'(rho') d(rho')] = Integral [P(rho) d(rho)]
C
C And computes the individual transform scale factors and shifts
C from the observed histogram to the expected histogram.
C
C     Rho' = scale * rho + shift
C          Scale = d(rho') / d(rho)
C          Shift = rho' - scale * rho
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS, MBINSE
      DOUBLE PRECISION HIST(*), HISTE(*)
      DOUBLE PRECISION HSLOT(*), HSCAL(*), HSHIF(*)
      DOUBLE PRECISION DSLOT, DMIN
C local
      INTEGER I, J, DHS
      DOUBLE PRECISION ACC, ACC2, HSUM
C parameter
      DOUBLE PRECISION HALF, ZERO, EPS, RSMALL
      PARAMETER (HALF=0.5D0, ZERO=0.0D0, EPS=0.1D-4, RSMALL=0.1D-9)
C begin
C
C define the accuracy of a histogram (in percentage)
      ACC=EPS
      ACC2=EPS*HALF
C assume the increment of subscripts for histogram
      DHS=1
C
C initialization
      HSLOT(1)=ZERO
      J=0
      HSUM=ACC
C
C loop over all histogram slots
      DO I=1,MBINS
      HSUM=HSUM+HIST(I)
C
      DO WHILE (HSUM.GE.ACC2.AND.J.LE.MBINSE)
      J=J+DHS
      IF (J.LE.MBINSE) HSUM=HSUM-HISTE(J)
      END DO
C
      IF (J.LE.MBINSE) THEN
      IF (HISTE(J).LT.RSMALL) THEN
      HSLOT(I+1)=FLOAT(J)+HSUM
      ELSE
      HSLOT(I+1)=FLOAT(J)+HSUM/HISTE(J)
      END IF
      ELSE IF (HSLOT(I).LT.MBINSE) THEN
      HSLOT(I+1)=MBINSE
      ELSE
      HSLOT(I+1)=HSLOT(I)+DHS
      END IF
C
      END DO
C
C compute individual scale factors and shifts
      DO I=1,MBINS
      HSCAL(I)=HSLOT(I+1)-HSLOT(I)
      HSHIF(I)=HSLOT(I+1)*DSLOT+DMIN-HSCAL(I)*(FLOAT(I)*DSLOT+DMIN)
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHTRAN(MBINS,MPACK,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &           RRHO,IRHO,HSCAL,HSHIF,DSLOT,DMIN)
C
C Apply the individual scale factors and shifts to map
C NOTE: modified electron densities overwrite to the original map.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION HSCAL(*), HSHIF(*)
      DOUBLE PRECISION DSLOT, DMIN
C local
      INTEGER A, B, C, J
      DOUBLE PRECISION DENS
C parameter
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
C begin
C
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (MPACK(A,B,C).EQ.1) THEN
      DENS=RRHO(A,B,C)
      J=MIN(MBINS,MAX(1,NINT((DENS-DMIN)/DSLOT)))
      RRHO(A,B,C)=HSCAL(J)*DENS+HSHIF(J)
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHSCAL(MBINS,MPACK,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,QHERM,
     &           RRHO,IRHO,DAVE,HSCALE)
C
C Apply the overall scale factor so that the r.m.s. of the map
C agrees with that of the desired histogram
C
C NOTE: keep the same mean of the map
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      LOGICAL QHERM
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MPACK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION DAVE, HSCALE
C local
      INTEGER A, B, C
      DOUBLE PRECISION DENS
C begin
C
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (MPACK(A,B,C).EQ.1) THEN
      DENS=RRHO(A,B,C)
      RRHO(A,B,C)=HSCALE*(DENS-DAVE)+DAVE
      END IF
      END DO
      END DO
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHISTE(MBINS,HISTE,DSLOT,DMIN,DAVEE,DRMSE)
C
C Estimate the mean and r.m.s. of the map that corresponds
C to the expected histogram
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS
      DOUBLE PRECISION HISTE(*)
      DOUBLE PRECISION DSLOT, DMIN, DAVEE, DRMSE
C local
      INTEGER I
      DOUBLE PRECISION HSUME, TEMP
C parameter
      DOUBLE PRECISION ZERO, HALF, RSMALL
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, RSMALL=0.1D-6)
C begin
      HSUME=ZERO
      DAVEE=ZERO
      DRMSE=ZERO
      DO I=1,MBINS
      TEMP=FLOAT(I)*DSLOT+DMIN
      HSUME=HSUME+HISTE(I)
      DAVEE=DAVEE+HISTE(I)*TEMP
      DRMSE=DRMSE+HISTE(I)*TEMP*TEMP
      END DO
      IF (HSUME.GT.RSMALL) THEN
      DAVEE=DAVEE/HSUME
      DRMSE=SQRT(DRMSE/HSUME-DAVEE*DAVEE)
      ELSE
      WRITE(6,'(A)') ' %XMHISTE-ERR: the expected hsitogram is zero'
      CALL WRNDIE(-5,'XMHISTE',' histogram is zero.')
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE INTRPOLAT(N,X,Y,U,V)
C
C Interpolation at u for v=f(u)
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER N
      DOUBLE PRECISION X(N), Y(N)
      DOUBLE PRECISION U, V
C local
      INTEGER I, J
      DOUBLE PRECISION TEMP
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      V=ZERO
      DO I=1,N
      TEMP=ONE
      DO J=1,N
      IF (J.NE.I) THEN
      TEMP=TEMP*(U-X(J))/(X(I)-X(J))
      END IF
      END DO
      V=V+TEMP*Y(I)
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMHTHEO(MBINS,HISTE,HRESO,DSLOT,DMIN)
C
C Computes predict (expected) histogram according to Main's formula
C Resolution limits: 4.5 to 0.9 angstrom.
C Modification:  the curve is smooth by averaging over Eq.(7a).
C
C Reference:
C
C Main, P. (1990), Acta Cryst. A46, 507-509.
C "A Formula for Electron Density Histograms for Equal-Atom Structures".
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MBINS
      DOUBLE PRECISION HISTE(*)
      DOUBLE PRECISION HRESO, DSLOT, DMIN
C local
      INTEGER I
      INTEGER NTAB1, NTAB2
      DOUBLE PRECISION T1RE(4), T1RM(4), T1S2(4), T1RA(4)
      DOUBLE PRECISION T1R1(4), T1R2(4), T1R0(4)
      DOUBLE PRECISION T2RE(5), T2RM(5), T2S2(5), T2RA(5)
      DOUBLE PRECISION T2R1(5), T2R2(5), T2R0(5)
      DOUBLE PRECISION RHOM, S22, AA, RHO1, RHO2, RHO0
      DOUBLE PRECISION Y1, Y2, SLG, R12, S1, S2, S0
      DOUBLE PRECISION SA, SB, SC, SD
      DOUBLE PRECISION DENS, P1, P2, P3
C parameter
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, THREE, FOUR, R100
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0, R100=100.0D0)
C data block
      DATA NTAB1/4/, NTAB2/5/
C Table 1
      DATA T1RE/4.5D0, 3.5D0, 2.8D0, 2.2D0/
      DATA T1RM/0.119D0, 0.036D0, -0.028D0, -0.042D0/
      DATA T1S2/0.065D0, 0.071D0, 0.084D0, 0.114D0/
      DATA T1RA/0.396D0, 0.378D0, 0.322D0, 0.225D0/
      DATA T1R1/0.167D0, 0.098D0, 0.059D0, 0.122D0/
      DATA T1R2/0.845D0, 1.200D0, 1.243D0, 1.283D0/
      DATA T1R0/1.25D0, 1.52D0, 1.93D0, 2.73D0/
C Table 2
      DATA T2RE/2.2D0, 1.8D0, 1.4D0, 1.1D0, 0.9D0/
      DATA T2RM/-0.035D0, -0.022D0, 0.010D0, 0.011D0, 0.029D0/
      DATA T2S2/0.121D0, 0.134D0, 0.142D0, 0.211D0, 0.176D0/
      DATA T2RA/0.230D0, 0.139D0, 0.074D0, 0.039D0, 0.020D0/
      DATA T2R1/0.275D0, 0.395D0, 0.535D0, 0.375D0, 0.356D0/
      DATA T2R2/1.382D0, 1.540D0, 1.843D0, 1.382D0, 1.463D0/
      DATA T2R0/2.73D0, 3.64D0, 5.50D0, 8.09D0, 11.45D0/
C begin
C
      IF (HRESO.LE.T1RE(1).AND.HRESO.GE.T2RE(5)) THEN
C
C interpolations
      IF (HRESO.GE.T1RE(4)) THEN
C at Table 1
      CALL INTRPOLAT(NTAB1,T1RE,T1RM,HRESO,RHOM)
      CALL INTRPOLAT(NTAB1,T1RE,T1S2,HRESO,S22)
      CALL INTRPOLAT(NTAB1,T1RE,T1RA,HRESO,AA)
      CALL INTRPOLAT(NTAB1,T1RE,T1R1,HRESO,RHO1)
      CALL INTRPOLAT(NTAB1,T1RE,T1R2,HRESO,RHO2)
      CALL INTRPOLAT(NTAB1,T1RE,T1R0,HRESO,RHO0)
C at Table 2
      ELSE
      CALL INTRPOLAT(NTAB2,T2RE,T2RM,HRESO,RHOM)
      CALL INTRPOLAT(NTAB2,T2RE,T2S2,HRESO,S22)
      CALL INTRPOLAT(NTAB2,T2RE,T2RA,HRESO,AA)
      CALL INTRPOLAT(NTAB2,T2RE,T2R1,HRESO,RHO1)
      CALL INTRPOLAT(NTAB2,T2RE,T2R2,HRESO,RHO2)
      CALL INTRPOLAT(NTAB2,T2RE,T2R0,HRESO,RHO0)
      END IF
C
C pre-compute coefficients
      Y1=EXP(-(RHO1-RHOM)*(RHO1-RHOM)/S22)
      SLG=SQRT(LOG(RHO0/RHO2))
      Y2=AA/RHO2*SLG
      S1=-(RHO1-RHOM)/S22/TWO*Y1
      S2=-AA/(RHO2*RHO2)*(SLG+HALF/SLG)
      R12=RHO1-RHO2
      S0=(Y1-Y2)/R12
      IF (HRESO.GE.T1RE(4)) THEN
      SA=(S1+S2-TWO*S0)/(R12*R12)
      ELSE
      SA=(TWO*S0+S2-SQRT(THREE*S2*(FOUR*S0-S2)))/(TWO*R12*R12)
      END IF
      SB=(S0-S2)/R12-(RHO1+TWO*RHO2)*SA
      SC=S2-THREE*SA*RHO2*RHO2-TWO*SB*RHO2
      SD=Y2-((SA*RHO2+SB)*RHO2+SC)*RHO2
C
C cumpute expected histogram
C average with P1 to smooth the curve
      DO I=1,MBINS
      DENS=(FLOAT(I)-HALF)*DSLOT+DMIN
      P1=EXP(-(DENS-RHOM)*(DENS-RHOM)/S22)
      IF (DENS.LT.RHO1) THEN
      HISTE(I)=P1
      ELSE IF (DENS.LT.RHO2) THEN
      P2=((SA*DENS+SB)*DENS+SC)*DENS+SD
      HISTE(I)=(P1+P2)/TWO
      ELSE IF (DENS.LT.RHO0) THEN
      P3=AA/DENS*SQRT(LOG(RHO0/DENS))
      HISTE(I)=(P1+P3)/TWO
      ELSE
      HISTE(I)=ZERO
      END IF
      END DO
C
C normalize to percentage
      AA=ZERO
      DO I=1,MBINS
      AA=AA+HISTE(I)
      END DO
      S1=R100/AA
      DO I=1,MBINS
      HISTE(I)=S1*HISTE(I)
      END DO
C
      ELSE
      WRITE(6,'(A,A,F8.4,A,F8.4)') ' %XMHTHEO-ERR:',
     & ' expected histogram avaliable only from resolution',
     & T1RE(1),' to',T2RE(5)
      WRITE(6,'(A,A,F8.4,A,F8.4)') ' %XMHTHEO-ERR:',
     & ' check the resolution (HRESO) ',HRESO
      CALL WRNDIE(-5,'XMHTHEO',
     & ' no expected histogram at such resolution')
      END IF
C
      RETURN
      END
C
