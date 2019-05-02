      SUBROUTINE XMASK(NA,NB,NC,NAP,NBP,NCPL,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,XRHONUM,XRHONAM,HPRRHO,HPIRHO,HPRHOMA,NRHO,
     &  IRHO,NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &  HPH,HPK,HPL,
     &  XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &  HPMULT,HPTYPE,XRSYGP,XRSYIV)
C
C Creates density masks based on atomic positions.
C
C Front-end routine for XMASK2
C
C Author: Axel T. Brunger
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
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C local
C pointer
      INTEGER FLAGS, XRVDW, XL, YL, ZL
C begin
      FLAGS=ALLHP(INTEG4(NATOM))
      XRVDW=ALLHP(IREAL8(NATOM*XNNSYM))
      XL=ALLHP(IREAL8(NATOM*XNNSYM))
      YL=ALLHP(IREAL8(NATOM*XNNSYM))
      ZL=ALLHP(IREAL8(NATOM*XNNSYM))
C
      CALL XMASK2(NATOM,HEAP(FLAGS),HEAP(XRVDW),HEAP(XL),
     &   HEAP(YL),HEAP(ZL),MAXCN,CNBVR,LOOKUP,NA,NB,NC,NAP,NBP,
     &   NCPL,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &   HPH,HPK,HPL,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &   HPMULT,HPTYPE,XRSYGP,XRSYIV)
C
      CALL FREHP(ZL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(YL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(XL,IREAL8(NATOM*XNNSYM))
      CALL FREHP(XRVDW,IREAL8(NATOM*XNNSYM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMASK2(NATOM,FLAGS,XRVDW,XL,
     &   YL,ZL,MAXCN,CNBVR,LOOKUP,NA,NB,NC,NAP,NBP,
     &   NCP,NAPP,NBPP,NCPP,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRCELL,ADEPTH,
     &   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,MAPR,
     &   XRHONUM,XRHONAM,HPRRHO,HPIRHO,
     &   HPRHOMA,NRHO,IRHO,
     &   NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &   HPH,HPK,HPL,
     &   XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &   HPMULT,HPTYPE,XRSYGP,XRSYIV)
C
C Parsing routine for computing a mask.
C
C Author: Axel T. Brunger
C =======================
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
      INCLUDE 'consta.inc'
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
      DOUBLE PRECISION XRTR(3,3), XRINTR(3,3), XRVOL
      INTEGER XRMREF, XRNREF
      INTEGER HPH, HPK, HPL
      INTEGER XSFNUM
      CHARACTER*(*) XSFNAM(*), XSFTYPE(*)
      INTEGER HPSF(*)
      INTEGER HPMULT
      INTEGER HPTYPE
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
C local
      DOUBLE PRECISION SOLRAD, SOLDEN, SRAD, PROVAL, VINSIDE, VOUTSIDE
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      INTEGER IAT, LNATOM, LLNATOM, INSYM, I
      CHARACTER*4 SMODE
      INTEGER ISTO
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION SECS, SECS2
      INTEGER NSHELL, LENG
      DOUBLE PRECISION THICK, AVERVDW
      LOGICAL COND, ERR, QAVER
C pointer
      INTEGER HPRMAP, HPIMAP, VOLUME, XSYM, XX, YY, ZZ
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, SIX, R100, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, SIX=6.0D0)
      PARAMETER (R100=100.D0, HALF=0.5D0)
C begin
C
C default values
      ERR=.FALSE.
      PROVAL=ZERO
      SOLDEN=ONE
      SOLRAD=1.4D0
      SRAD=ZERO
      ISTO=1
      LENG=LEN(XRHONAM(ISTO))
      CALL TRIMM(XRHONAM(ISTO),LENG)
      SMODE='VDW'
      NSHELL=1
      THICK=ONE
      QAVER=.FALSE.
      XMIN=ZERO
      XMAX=ONE
      YMIN=ZERO
      YMAX=ONE
      ZMIN=ZERO
      ZMAX=ONE
C
C
C default selection: none
      LNATOM=0
C
C
      CALL PUSEND('MASK>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MASK>')
      CALL MISCOM('MASK>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-mask')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'AVER') THEN
      CALL NEXTLO('AVERaging-mode=',QAVER)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NSHE'.OR.WD(1:4).EQ.'PART') THEN
      CALL NEXTI('PARTitions=',NSHELL)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'XMIN') THEN
      CALL NEXTF('XMIN=',XMIN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'YMIN') THEN
      CALL NEXTF('XMIN=',YMIN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ZMIN') THEN
      CALL NEXTF('XMIN=',ZMIN)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'XMAX') THEN
      CALL NEXTF('XMIN=',XMAX)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'YMAX') THEN
      CALL NEXTF('XMIN=',YMAX)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ZMAX') THEN
      CALL NEXTF('XMIN=',ZMAX)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'THIC') THEN
      CALL NEXTF('THICkness=',THICK)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(FLAGS,LNATOM,X,Y,Z,.TRUE.)
      CALL MAKIND(FLAGS,NATOM,LNATOM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PROB'.OR.WD(1:4).EQ.'SOLR') THEN
      CALL NEXTF('solvent-radius=',SOLRAD)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('MODE=',SMODE)
      IF (SMODE.NE.'VDW'.AND.SMODE.NE.'SIGM'.AND.SMODE.NE.'PROB'
     &   .AND.SMODE.NE.'FBOX') THEN
      CALL DSPERR('XMASKX','unknown value for MODE.')
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SHRI') THEN
      CALL NEXTF('shrink-radius=',SRAD)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TO') THEN
      CALL NEXTWD('TO=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('TO=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'TO=',XRHONAM(ISTO)
      ELSE
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      ISTO=I
      LENG=LEN(XRHONAM(ISTO))
      CALL TRIMM(XRHONAM(ISTO),LENG)
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XMASK-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XMASK','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     & 'mask---------------------------------------'
      WRITE(6,'(5A,L1)') ' |  TO = ',XRHONAM(ISTO)(1:LENG),
     &   '  MODE=',SMODE,' AVERaging_mode=',QAVER
      IF (SMODE.EQ.'VDW'.OR.SMODE.EQ.'SIGM'.OR.SMODE.EQ.'PROB') THEN
      WRITE(6,'(A,F7.4)') ' | SHRINk = ',SRAD
      WRITE(6,'(A,F7.4,A,A,I6,A,F7.4,A)') ' |  PROBe=',SOLRAD,' A',
     & ' PARTition=',NSHELL,' THICkness=',THICK,' A'
      ELSEIF (SMODE.EQ.'FBOX') THEN
      WRITE(6,'(A,F7.4)')
     &  ' |  XMIN=',XMIN,' XMAX=',XMAX,
     &  ' |  YMIN=',YMIN,' YMAX=',YMAX,
     &  ' |  ZMIN=',ZMIN,' ZMAX=',ZMAX
      END IF
      WRITE(6,'(2A)') ' |--------------------------',
     & '---------------------------------------------------'
      ELSE
      CALL CHKEND('MASK>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (LNATOM.EQ.0) THEN
      CALL WRNDIE(-5,'XMASK','zero atoms selected.')
      ELSEIF (.NOT.ERR) THEN
C
      IF (SMODE.EQ.'VDW'.OR.SMODE.EQ.'SIGM'.OR.SMODE.EQ.'PROB') THEN
C
C gather coordinates and radii of atoms for all selected
C atoms
      DO IAT=1,LNATOM
      XL(IAT)=X(FLAGS(IAT))
      YL(IAT)=Y(FLAGS(IAT))
      ZL(IAT)=Z(FLAGS(IAT))
      END DO
C
C Make sure that lookup tables for van der Waals radii are set
      IF (UPNBLK) THEN
      UPNBLK=.FALSE.
      CALL NBUPDA
      END IF
C
      IF (SMODE.EQ.'VDW') THEN
C
C use half the vdw radius
      DO IAT=1,LNATOM
      XRVDW(IAT)=HALF*CNBVR(LOOKUP(FLAGS(IAT)),LOOKUP(FLAGS(IAT)))
     &  +SOLRAD
      END DO
      ELSEIF (SMODE.EQ.'SIGM') THEN
C
C use half the sigma value
      DO IAT=1,LNATOM
      XRVDW(IAT)=HALF*CNBVR(LOOKUP(FLAGS(IAT)),LOOKUP(FLAGS(IAT)))
     &  *TWO**(-ONE/SIX)+SOLRAD
      END DO
C
      ELSEIF (SMODE.EQ.'PROB') THEN
C
C use probe radius
      DO IAT=1,LNATOM
      XRVDW(IAT)=SOLRAD
      END DO
      END IF
C
      IF (LNATOM.LE.0) THEN
      CALL WRNDIE(-5,'XMASK','Zero atoms selected.')
      ELSE
      WRITE(6,'(A,I8,A)') ' XMASK: ',LNATOM,
     &  ' atoms have been selected for mask calculation.'
C
C apply non-crystallographic symmetry operations
      IF (XNNSYM.GT.1) THEN
      LLNATOM=0
C loop over non-crystallographic symmetry operators
      DO INSYM=1,XNNSYM
C apply the NCS-operator to the atomic coordinates.
      DO IAT=1,LNATOM
      LLNATOM=LLNATOM+1
      XL(LLNATOM)=NCSOP(INSYM,1,1)*XL(IAT) + NCSOP(INSYM,1,2)*YL(IAT)
     &          + NCSOP(INSYM,1,3)*ZL(IAT) + NCSOP(INSYM,1,4)
      YL(LLNATOM)=NCSOP(INSYM,2,1)*XL(IAT) + NCSOP(INSYM,2,2)*YL(IAT)
     &          + NCSOP(INSYM,2,3)*ZL(IAT) + NCSOP(INSYM,2,4)
      ZL(LLNATOM)=NCSOP(INSYM,3,1)*XL(IAT) + NCSOP(INSYM,3,2)*YL(IAT)
     &          + NCSOP(INSYM,3,3)*ZL(IAT) + NCSOP(INSYM,3,4)
      XRVDW(LLNATOM)=XRVDW(IAT)
      END DO
      END DO
      LNATOM=LLNATOM
      WRITE(6,'(A)') ' XMASK: applying non-crystallographic symmetry.'
      END IF
C
      IF (QAVER) THEN
      XSYM=ALLHP(INTEG4(LNATOM))
      XX=ALLHP(INTEG4(LNATOM))
      YY=ALLHP(INTEG4(LNATOM))
      ZZ=ALLHP(INTEG4(LNATOM))
      ELSE
      XSYM=1
C
C set to dummy pointers since arrays are not used  (ATB 01/07/07)
      XX=1
      YY=1
      ZZ=1
      END IF
C
      END IF
      END IF
C
C define the asymmetric unit if required
      IF (XRMAP) CALL XASUDEF(MAPR,NA,NB,NC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     &                   MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                   XRCELL,ADEPTH,
     &                   ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,ARPNDB,ARPNMLT,
     &                   HPRHOMA,NRHO,NMASK,XRSYGP,XRSYIV)
      XRMAP=.FALSE.
C
C
C get the pointer for the map.  Allocate space if required.
      IF (HPRRHO(ISTO).EQ.0) THEN
      CALL XMAPAL(HPRRHO(ISTO),HPIRHO(ISTO),QHERM,NRHO,IRHO)
      END IF
      HPRMAP=HPRRHO(ISTO)
      HPIMAP=HPIRHO(ISTO)
C
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      IF (SMODE.EQ.'VDW'.OR.SMODE.EQ.'SIGM'.OR.SMODE.EQ.'PROB') THEN
      IF (LNATOM.GT.0) THEN
      CALL XMASKX(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   HEAP(HPRMAP),HEAP(HPIMAP),HEAP(HPRHOMA),
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,SOLRAD,SOLDEN,PROVAL,SRAD,LNATOM,XL,YL,ZL,XRVDW,
     &   VINSIDE,VOUTSIDE,NSHELL,THICK,XRSYGP,XRSYIV,QAVER,
     &   HEAP(XSYM),HEAP(XX),HEAP(YY),HEAP(ZZ))
      END IF
C
      ELSEIF (SMODE.EQ.'FBOX') THEN
      CALL XFBOX(XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &           NA,NB,NC,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           HEAP(HPRMAP),HEAP(HPIMAP),HEAP(HPRHOMA),
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,
     &           XRCELL,XRTR,XRINTR,
     &           VINSIDE,VOUTSIDE,XRSYGP,XRSYIV,
     &           QAVER)
      END IF
C
      IF (SMODE.EQ.'VDW'.OR.SMODE.EQ.'SIGM'.OR.SMODE.EQ.'PROB') THEN
      IF (QAVER) THEN
      CALL FREHP(ZZ,INTEG4(LNATOM))
      CALL FREHP(YY,INTEG4(LNATOM))
      CALL FREHP(XX,INTEG4(LNATOM))
      CALL FREHP(XSYM,INTEG4(LNATOM))
      END IF
C
C compute average exclusion radius
      AVERVDW=ZERO
      DO IAT=1,LNATOM
      AVERVDW=AVERVDW+XRVDW(IAT)
      END DO
      AVERVDW=AVERVDW/LNATOM-SOLRAD
      WRITE(6,'(A,F8.4,A)')
     & ' XMASK: average mask radius around selected atoms',
     & AVERVDW,' A'
      WRITE(6,'(A,F8.4,A,F8.4)')
     & ' XMASK: probe radius=',SOLRAD,' shrink radius=',SRAD
      CALL DECLAR( 'AVERADIUS','DP',' ',DBCOMP,AVERVDW)
      AVERVDW=AVERVDW+SOLRAD
C
      IF (NSHELL.GT.1) THEN
      VOLUME=ALLHP(IREAL8(NSHELL+1))
      CALL XMASKVOL(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &              HEAP(HPRHOMA),HEAP(HPRMAP),XRNSYM,AVERVDW,SRAD,
     &              THICK,NSHELL,VINSIDE,VOUTSIDE,HEAP(VOLUME),
     &              XRHONAM(ISTO)(1:LENG))
      CALL FREHP(VOLUME,IREAL8(NSHELL+1))
      ELSE
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume inside mask=',
     &  VINSIDE*R100,'% (',XRHONAM(ISTO)(1:LENG),'<=0)'
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume outside mask=',
     &  VOUTSIDE*R100,'% (',XRHONAM(ISTO)(1:LENG),'=1)'
      END IF
C
      ELSEIF (SMODE.EQ.'FBOX') THEN
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume inside mask=',
     &  VINSIDE*R100,'% (',XRHONAM(ISTO)(1:LENG),'<=0)'
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume outside mask=',
     &  VOUTSIDE*R100,'% (',XRHONAM(ISTO)(1:LENG),'=1)'
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MASKX2: CPU-time:  mask-calc=',SECS2-SECS
      END IF
C
C
C declare symbols (unitcell volume inside and outside of mask )
      CALL DECLAR( 'INSIDE','DP',' ',DBCOMP,VINSIDE)
      CALL DECLAR( 'OUTSIDE','DP',' ',DBCOMP,VOUTSIDE)
      END IF
C
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XMASKX(NA,NB,NC,MAASY,
     &   MBASY,MCASY,NAASY,NBASY,NCASY,
     &   RRHO,IRHO,RHOMASK,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,SOLRAD,SOLDEN,PROVAL,SRAD,NATOM,X,Y,Z,XRVDW,
     &   VINSIDE,VOUTSIDE,NSHELL,THICK,XRSYGP,XRSYIV,QAVER,XSYM,
     &   XX,YY,ZZ)
C
C Front-end routine for XMASKX2.  Determines the size
C of the temporary mask that surrounds the selected
C atoms plus their radii.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*), XRTR(3,3), XRINTR(3,3)
      LOGICAL QHERM
      DOUBLE PRECISION SOLRAD, SOLDEN, PROVAL, SRAD
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XRVDW(*)
      DOUBLE PRECISION VINSIDE, VOUTSIDE
      INTEGER NSHELL
      DOUBLE PRECISION THICK
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      LOGICAL QAVER
      INTEGER XSYM(*), XX(*), YY(*), ZZ(*)
C local
      INTEGER CUSHNA, CUSHNB, CUSHNC, CUSHMA, CUSHMB, CUSHMC
      INTEGER IAT, NN, I, J, SYM, AS, BS, CS, AA, BB, CC
      LOGICAL COND
      DOUBLE PRECISION DX, DY, DZ
      INTEGER AX, BX, CX, NABC(3)
      DOUBLE PRECISION ABC(3,3), ABCINV(3,3), FTAS, FTBS, FTCS
      DOUBLE PRECISION ADD
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C pointer
      INTEGER MASK, XSYMASK
C begin
C add to cushion for the radial shell model
      ADD=(NSHELL-1)*THICK
C
C compute sublattice transformation matrices
      NABC(1)=NA
      NABC(2)=NB
      NABC(3)=NC
      DO I=1,3
      DO J=1,3
      ABC(I,J)=NABC(I)*XRTR(I,J)
      END DO
      END DO
      DO I=1,3
      DO J=1,3
      ABCINV(I,J)=XRINTR(I,J)/NABC(J)
      END DO
      END DO
C
C compute norm of A*, B*, C* in FFT grid box
      FTAS=XRCELL(7)*NA
      FTBS=XRCELL(8)*NB
      FTCS=XRCELL(9)*NC
C
C initialize extent of cushion mask
      CUSHNA=NAASY
      CUSHNB=NBASY
      CUSHNC=NCASY
      CUSHMA=MAASY
      CUSHMB=MBASY
      CUSHMC=MCASY
C
      DO IAT=1,NATOM
C
C initialize symmetry and lattice translation info
      IF (QAVER) THEN
      XSYM(IAT)=1
      XX(IAT)=0
      YY(IAT)=0
      ZZ(IAT)=0
      END IF
C
C determine indices of lattice point closest to atomic position x,y,z
      DX= ABC(1,1)*(X(IAT))
     &   +ABC(1,2)*(Y(IAT))
     &   +ABC(1,3)*(Z(IAT))
      AX=NINT(DX)
      DY= ABC(2,1)*(X(IAT))
     &   +ABC(2,2)*(Y(IAT))
     &   +ABC(2,3)*(Z(IAT))
      BX=NINT(DY)
      DZ= ABC(3,1)*(X(IAT))
     &   +ABC(3,2)*(Y(IAT))
     &   +ABC(3,3)*(Z(IAT))
      CX=NINT(DZ)
C
C map AX, BX, CX into asymmetric unit
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*AX+XRSYMM(SYM,1,2)*BX+XRSYMM(SYM,1,3)*CX
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*AX+XRSYMM(SYM,2,2)*BX+XRSYMM(SYM,2,3)*CX
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*AX+XRSYMM(SYM,3,2)*BX+XRSYMM(SYM,3,3)*CX
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
C
C map coordinates into asymmetric unit.  The result
C will overwrite the X(IAT), Y(IAT), Z(IAT) array!
      X(IAT)=XRSYMM(SYM,1,1)*DX+XRSYMM(SYM,1,2)*DY+
     &         XRSYMM(SYM,1,3)*DZ
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH + AS-AA
      Y(IAT)=XRSYMM(SYM,2,1)*DX+XRSYMM(SYM,2,2)*DY+
     &         XRSYMM(SYM,2,3)*DZ
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH + BS-BB
      Z(IAT)=XRSYMM(SYM,3,1)*DX+XRSYMM(SYM,3,2)*DY+
     &         XRSYMM(SYM,3,3)*DZ
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH + CS-CC
C
C update symmetry and lattice translation info
      IF (QAVER) THEN
      XSYM(IAT)=SYM
      XX(IAT)=(AS-AA)/NA
      YY(IAT)=(BS-BB)/NB
      ZZ(IAT)=(CS-CC)/NC
      END IF
C
      DX=X(IAT)
      DY=Y(IAT)
      DZ=Z(IAT)
C
C copy index
      AX=NINT(DX)
      BX=NINT(DY)
      CX=NINT(DZ)
      ELSE
      CALL WRNDIE(-5,'MASKX','Internal Error No. 3')
      END IF
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
      CUSHNA=MAX(CUSHNA,INT((XRVDW(IAT)+ADD)*FTAS-DX+RSMALL)+AX)
      CUSHNB=MAX(CUSHNB,INT((XRVDW(IAT)+ADD)*FTBS-DY+RSMALL)+BX)
      CUSHNC=MAX(CUSHNC,INT((XRVDW(IAT)+ADD)*FTCS-DZ+RSMALL)+CX)
      CUSHMA=MIN(CUSHMA,-INT((XRVDW(IAT)+ADD)*FTAS+DX+RSMALL)+AX)
      CUSHMB=MIN(CUSHMB,-INT((XRVDW(IAT)+ADD)*FTBS+DY+RSMALL)+BX)
      CUSHMC=MIN(CUSHMC,-INT((XRVDW(IAT)+ADD)*FTCS+DZ+RSMALL)+CX)
      END DO
C
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(6(A,I6))')
     &  '  CUSHA=',CUSHMA,',...,',CUSHNA,
     &  '  CUSHB=',CUSHMB,',...,',CUSHNB,
     &  '  CUSHC=',CUSHMC,',...,',CUSHNC
      END IF
C
C allocate space for temporary mask
      NN=(-CUSHMA+CUSHNA+1)*
     &   (-CUSHMB+CUSHNB+1)*
     &   (-CUSHMC+CUSHNC+1)
      MASK=ALLHP(INTEG4(NN))
      IF (QAVER) THEN
      XSYMASK=ALLHP(INTEG4(NN))
      ELSE
      XSYMASK=1
      END IF
C
      CALL XMASKX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RHOMASK,HEAP(MASK),
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,SOLDEN,PROVAL,SRAD,ABC,
     &           ABCINV,NATOM,X,Y,Z,XRVDW,FTAS,FTBS,FTCS,
     &           VINSIDE,VOUTSIDE,NSHELL,THICK,ADD,XRSYGP,XRSYIV,
     &           QAVER,XSYM,XX,YY,ZZ,HEAP(XSYMASK))
C
      IF (QAVER) THEN
      CALL FREHP(XSYMASK,INTEG4(NN))
      END IF
C
      CALL FREHP(MASK,INTEG4(NN))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XMASKX2(NA,NB,NC,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RHOMASK,MASK,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,SOLDEN,PROVAL,SRAD,ABC,
     &           ABCINV,NATOM,X,Y,Z,XRVDW,FTAS,FTBS,FTCS,
     &           VINSIDE,VOUTSIDE,NSHELL,THICK,ADD,XRSYGP,XRSYIV,
     &           QAVER,XSYM,XX,YY,ZZ,XSYMASK)
C
C Computes a mask of the atomic model specified by coordinates
C (X,Y,Z, NATOM) on a lattice given by i*A + j*B + k*C
C where i,j,k are integers and (A,B,C)=ABCINV.  Grid points
C within XRVDW plus SOLRAD around any atom are set to PROVAL
C Grid points outside atoms are set to SOLDEN
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER CUSHMC, CUSHNC, CUSHMB, CUSHNB, CUSHMA, CUSHNA
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MASK(CUSHMA:CUSHNA,CUSHMB:CUSHNB,CUSHMC:CUSHNC)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION SOLDEN, PROVAL, SRAD, ABC(3,3)
      DOUBLE PRECISION ABCINV(3,3)
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XRVDW(*)
      DOUBLE PRECISION FTAS, FTBS, FTCS
      DOUBLE PRECISION VINSIDE, VOUTSIDE
      INTEGER NSHELL
      DOUBLE PRECISION THICK, ADD
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      LOGICAL QAVER
      INTEGER XSYM(*), XX(*), YY(*), ZZ(*)
      INTEGER XSYMASK(CUSHMA:CUSHNA,CUSHMB:CUSHNB,CUSHMC:CUSHNC)
C local
      INTEGER A, B, C, AX, BX, CX, IAT, AA, BB, CC
      INTEGER AS, BS, CS, ASS, BSS, CSS
      DOUBLE PRECISION RCUT2, DX, DY, DZ, DELX, DELY, DELZ
      DOUBLE PRECISION DELXB, DELYB, DELZB, DELXC, DELYC, DELZC
      DOUBLE PRECISION INSIDE, OUTSIDE
      INTEGER AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN, SYM
      LOGICAL COND, ERR
      DOUBLE PRECISION TEMP, THICKINV
      INTEGER ISHELL
      INTEGER SYMMOP, LTX, LTY, LTZ
      INTEGER LSYMMOP, LLTX, LLTY, LLTZ
      REAL KEY
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      ERR=.FALSE.
C
      IF (THICK.GT.RSMALL) THEN
      THICKINV=ONE/THICK
      ELSE
      THICKINV=ZERO
      END IF
C
C initialize mask
      DO C=CUSHMC,CUSHNC
      DO B=CUSHMB,CUSHNB
      DO A=CUSHMA,CUSHNA
      MASK(A,B,C)=MAX(1,NSHELL)
C
C initialize symmetry and lattice translation info matrix
      IF (QAVER) THEN
      XSYMASK(A,B,C)=-9999
      END IF
C
      END DO
      END DO
      END DO
C
C for all atoms set all grid points within XRVDW to 1.
      DO IAT=1,NATOM
C
      RCUT2=XRVDW(IAT)**2
C
C determine indices of lattice point closest to atomic position x,y,z
      DX=X(IAT)
      AX=NINT(DX)
      DY=Y(IAT)
      BX=NINT(DY)
      DZ=Z(IAT)
      CX=NINT(DZ)
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
C convert this difference into real space
      DELX=ABCINV(1,3)*DZ+ABCINV(1,2)*DY+ABCINV(1,1)*DX
      DELY=ABCINV(2,3)*DZ+ABCINV(2,2)*DY+ABCINV(2,1)*DX
      DELZ=ABCINV(3,3)*DZ+ABCINV(3,2)*DY+ABCINV(3,1)*DX
C
C determine AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN of a
C box isomorphous to the unitcell geometry which surrounds
C the sphere with the atomic radius.  Compute maximum box.
      AMAXP=INT((XRVDW(IAT)+ADD)*FTAS-DX+RSMALL)
      BMAXP=INT((XRVDW(IAT)+ADD)*FTBS-DY+RSMALL)
      CMAXP=INT((XRVDW(IAT)+ADD)*FTCS-DZ+RSMALL)
      AMAXN=INT((XRVDW(IAT)+ADD)*FTAS+DX+RSMALL)
      BMAXN=INT((XRVDW(IAT)+ADD)*FTBS+DY+RSMALL)
      CMAXN=INT((XRVDW(IAT)+ADD)*FTCS+DZ+RSMALL)
C
      DO C=-CMAXN,CMAXP
      CC=C+CX
      DELXC=DELX+ABCINV(1,3)*C
      DELYC=DELY+ABCINV(2,3)*C
      DELZC=DELZ+ABCINV(3,3)*C
C
      DO B=-BMAXN,BMAXP
      BB=B+BX
      DELXB=DELXC+ABCINV(1,2)*B
      DELYB=DELYC+ABCINV(2,2)*B
      DELZB=DELZC+ABCINV(3,2)*B
C
      DO A=-AMAXN,AMAXP
      AA=A+AX
C
      IF (CC.LT.CUSHMC.OR.CC.GT.CUSHNC.OR.
     &    BB.LT.CUSHMB.OR.BB.GT.CUSHNB.OR.
     &    AA.LT.CUSHMA.OR.AA.GT.CUSHNA) THEN
      WRITE(6,'(A,I3)') ' aa,bb,cc=',AA,BB,CC
      WRITE(6,'(A,I3)') ' a,b,c=',A,B,C
      CALL WRNDIE(-5,'MASKX','Internal Error No. 1')
      END IF
C
      TEMP=(ABCINV(1,1)*A+DELXB)**2
     &    +(ABCINV(2,1)*A+DELYB)**2
     &    +(ABCINV(3,1)*A+DELZB)**2
      IF (RCUT2.GT.TEMP) THEN
      ISHELL=0
      ELSE
      ISHELL=MAX(1,MIN(NSHELL,INT((SQRT(TEMP)-XRVDW(IAT))*THICKINV)+1))
      END IF
C
      MASK(AA,BB,CC)=MIN(MASK(AA,BB,CC),ISHELL)
C
C define symmetry and lattice translation info
      IF (QAVER) THEN
C
C pack the symmetry operator and lattice translation info
      CALL MSKM2A(XSYM(IAT),XX(IAT),YY(IAT),ZZ(IAT),KEY,
     &            XRMSYM)
      XSYMASK(AA,BB,CC)=KEY
C
      END IF
C
      END DO
C
      END DO
      END DO
C
      END DO
C
C
C
C map all grid points outside asymmetric unit into
C asymmetric unit
C =================================================
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
      IF (MASK(AA,BB,CC).GE.0) THEN
C
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (.NOT.COND) THEN
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      ASS=XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BSS=XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CSS=XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
      AS=MOD(ASS+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BSS+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CSS+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      IF (COND) THEN
C
      IF (MASK(AA,BB,CC).LT.MASK(AS,BS,CS)) THEN
      MASK(AS,BS,CS)=MASK(AA,BB,CC)
C
C update symmetry and lattice translation info
C (product of operator corresponding to AA,BB,CC and SYM)
CCC modification ATB 4/27/08
      IF (QAVER) THEN
      IF (XSYMASK(AA,BB,CC).NE.-9999) THEN
C
C unpack matrix to get operator and lattice translation
C corresponding to AA, BB, CC
      CALL MSKA2M(REAL(XSYMASK(AA,BB,CC)),LSYMMOP,LLTX,LLTY,LLTZ,
     &            XRMSYM)
C
C id of product of symmetry operators
      SYMMOP=XRSYGP(SYM,LSYMMOP)
C
      LTX=XRSYMM(SYM,1,1)*XRSYMM(LSYMMOP,1,4)
     &   +XRSYMM(SYM,1,2)*XRSYMM(LSYMMOP,2,4)
     &   +XRSYMM(SYM,1,3)*XRSYMM(LSYMMOP,3,4)
     &   +XRSYMM(SYM,1,4)
     &   -XRSYMM(SYMMOP,1,4)
      LTY=XRSYMM(SYM,2,1)*XRSYMM(LSYMMOP,1,4)
     &   +XRSYMM(SYM,2,2)*XRSYMM(LSYMMOP,2,4)
     &   +XRSYMM(SYM,2,3)*XRSYMM(LSYMMOP,3,4)
     &   +XRSYMM(SYM,2,4)
     &   -XRSYMM(SYMMOP,2,4)
      LTZ=XRSYMM(SYM,3,1)*XRSYMM(LSYMMOP,1,4)
     &   +XRSYMM(SYM,3,2)*XRSYMM(LSYMMOP,2,4)
     &   +XRSYMM(SYM,3,3)*XRSYMM(LSYMMOP,3,4)
     &   +XRSYMM(SYM,3,4)
     &   -XRSYMM(SYMMOP,3,4)
C
C concatenated lattice translation
      LTX=LTX/XRSYTH+(AS-ASS)/NA+XRSYMM(SYM,1,1)*LLTX
     &    +XRSYMM(SYM,1,2)*LLTY+XRSYMM(SYM,1,3)*LLTZ
      LTY=LTY/XRSYTH+(BS-BSS)/NB+XRSYMM(SYM,2,1)*LLTX
     &    +XRSYMM(SYM,2,2)*LLTY+XRSYMM(SYM,2,3)*LLTZ
      LTZ=LTZ/XRSYTH+(CS-CSS)/NC+XRSYMM(SYM,3,1)*LLTX
     &    +XRSYMM(SYM,3,2)*LLTY+XRSYMM(SYM,3,3)*LLTZ
C
C repack the symmetry operator and lattice translation info
      CALL MSKM2A(SYMMOP,LTX,LTY,LTZ,KEY,XRMSYM)
      XSYMASK(AS,BS,CS)=KEY
C
      END IF
      END IF
C
      END IF
      END IF
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'MASKX','Internal error no. 2')
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
C now we have to shrink the mask if this is required
C shrink the first shell only if the radial shell model
C ==================================================
      IF (SRAD.GT.RSMALL) THEN
      RCUT2=SRAD**2
C
      AMAXP=INT(SRAD*FTAS+RSMALL)
      BMAXP=INT(SRAD*FTBS+RSMALL)
      CMAXP=INT(SRAD*FTCS+RSMALL)
      AMAXN=INT(SRAD*FTAS+RSMALL)
      BMAXN=INT(SRAD*FTBS+RSMALL)
      CMAXN=INT(SRAD*FTCS+RSMALL)
C
C set all grid points outside asymmetric unit to 0
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (.NOT.COND) MASK(AA,BB,CC)=0
      END DO
      END DO
      END DO
C
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
C
C look at all points inside asymmetric unit and outside mask
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (COND.AND.MASK(AA,BB,CC).EQ.1) THEN
C
C loop over box
      DO C=-CMAXN+CC,CMAXP+CC
      IF (C.GE.CUSHMC.AND.C.LE.CUSHNC) THEN
      CX=C-CC
      DELXC=ABCINV(1,3)*CX
      DELYC=ABCINV(2,3)*CX
      DELZC=ABCINV(3,3)*CX
C
      DO B=-BMAXN+BB,BMAXP+BB
      IF (B.GE.CUSHMB.AND.B.LE.CUSHNB) THEN
      BX=B-BB
      DELXB=DELXC+ABCINV(1,2)*BX
      DELYB=DELYC+ABCINV(2,2)*BX
      DELZB=DELZC+ABCINV(3,2)*BX
C
      DO A=-AMAXN+AA,AMAXP+AA
      IF (A.GE.CUSHMA.AND.A.LE.CUSHNA) THEN
      AX=A-AA
      IF (RCUT2.GT.
     &  (ABCINV(1,1)*AX+DELXB)**2
     & +(ABCINV(2,1)*AX+DELYB)**2
     & +(ABCINV(3,1)*AX+DELZB)**2.AND.MASK(A,B,C).EQ.0) THEN
      MASK(A,B,C)=-1
      END IF
      END IF
      END DO
      END IF
      END DO
      END IF
      END DO
C
      END IF
C
      END DO
      END DO
      END DO
C
C second pass: map all grid points marked (-1) into asymmetric
C unit
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
      IF (MASK(AA,BB,CC).EQ.-1) THEN
C
      COND=(AA.GE.MAASY.AND.AA.LE.NAASY.AND.
     &      BB.GE.MBASY.AND.BB.LE.NBASY.AND.
     &      CC.GE.MCASY.AND.CC.LE.NCASY)
      IF (COND) COND=RHOMASK(AA,BB,CC).GT.0
      IF (COND) THEN
      MASK(AA,BB,CC)=1
      ELSE
C
C apply all symmetry operators to grid point until
C the symmetry-related grid point falls into the
C asymmetric unit
C
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
      AS=MOD(XRSYMM(SYM,1,1)*AA+XRSYMM(SYM,1,2)*BB+XRSYMM(SYM,1,3)*CC
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH+10000*NA-MAASY,NA) +MAASY
      BS=MOD(XRSYMM(SYM,2,1)*AA+XRSYMM(SYM,2,2)*BB+XRSYMM(SYM,2,3)*CC
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH+10000*NB-MBASY,NB) +MBASY
      CS=MOD(XRSYMM(SYM,3,1)*AA+XRSYMM(SYM,3,2)*BB+XRSYMM(SYM,3,3)*CC
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=RHOMASK(AS,BS,CS).GT.0
      IF (COND) MASK(AS,BS,CS)=1
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'MASKX','Internal error no. 3')
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
      END IF
C
C now we have to map the mask into RHO
C ====================================
C
C initialize RRHO, IRHO
      IF (QHERM) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (MASK(A,B,C).EQ.0) THEN
      IF (QAVER) THEN
      RRHO(A,B,C)=XSYMASK(A,B,C)
      IF (XSYMASK(A,B,C).EQ.-9999) THEN
      CALL WRNDIE(-5,'MASKX2','Internal error no. 6')
      END IF
      ELSE
      RRHO(A,B,C)=PROVAL
      END IF
      ELSE
      RRHO(A,B,C)=SOLDEN*MASK(A,B,C)
      END IF
      END IF
      END DO
      END DO
      END DO
      ELSE
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (MASK(A,B,C).EQ.0) THEN
      IF (QAVER) THEN
      RRHO(A,B,C)=XSYMASK(A,B,C)
      IF (XSYMASK(A,B,C).EQ.-9999) THEN
      CALL WRNDIE(-5,'MASKX2','Internal error no. 7')
      END IF
      ELSE
      RRHO(A,B,C)=PROVAL
      END IF
      ELSE
      RRHO(A,B,C)=SOLDEN*MASK(A,B,C)
      END IF
      IRHO(A,B,C)=ZERO
      END IF
      END DO
      END DO
      END DO
      END IF
C
C
C compute percentage of volume inside and outside of mask
      INSIDE=ZERO
      OUTSIDE=ZERO
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (ABS(RRHO(A,B,C)-SOLDEN).GE.R4SMAL) THEN
      INSIDE=INSIDE+ONE/RHOMASK(A,B,C)
      ELSEIF (ABS(RRHO(A,B,C)-SOLDEN*MASK(A,B,C)).LT.R4SMAL) THEN
      OUTSIDE=OUTSIDE+ONE/RHOMASK(A,B,C)
      END IF
      END IF
      END DO
      END DO
      END DO
C
      VINSIDE=XRNSYM*INSIDE/(NC*NB*NA)
      VOUTSIDE=XRNSYM*OUTSIDE/(NC*NB*NA)
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XMASKVOL(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           RHOMASK,RHO,XRNSYM,SOLRAD,SRAD,THICK,NSHELL,
     &           VINSIDE,VOUTSIDE,VOLUME,STO)
C
C Computes percentage volume for each mask shell
C
C Authors: J.-S. Jiang and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM
      DOUBLE PRECISION SOLRAD, SRAD, THICK
      INTEGER NSHELL
      DOUBLE PRECISION VINSIDE, VOUTSIDE
      DOUBLE PRECISION VOLUME(*)
      CHARACTER*(*) STO
C local
      INTEGER A, B, C, I, I1, NS
      DOUBLE PRECISION SOLRD1, SOLRD2, OUTSIDE
      DOUBLE PRECISION ZERO, ONE, R100
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, R100=100.0D0)
C begin
      NS=NSHELL+1
C
C initialization
      DO I=1,NS
      VOLUME(I)=ZERO
      END DO
C
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      I=NINT(RHO(A,B,C))+1
      IF (I.LT.1) I=1
      IF (I.GT.NS) I=NS
      VOLUME(I)=VOLUME(I)+ONE/RHOMASK(A,B,C)
      END IF
      END DO
      END DO
      END DO
C
      VINSIDE=XRNSYM*VOLUME(1)/(NC*NB*NA)
      VOLUME(1)=VINSIDE
      OUTSIDE=ZERO
      DO I=2,NS
      OUTSIDE=OUTSIDE+VOLUME(I)
      VOLUME(I)=XRNSYM*VOLUME(I)/(NC*NB*NA)
      END DO
      VOUTSIDE=XRNSYM*OUTSIDE/(NC*NB*NA)
C
      DO I=1,NS
      VOLUME(I)=VOLUME(I)*R100
      END DO
C
C some outputs
      WRITE(6,'(A,I3,A,F7.2,A,F8.3,3A)') ' XMASK:  shell=',0,
     & '   from=   0     to',SOLRAD-SRAD,
     & ' A,  volume=',VOLUME(1),'% (',STO,'=0)'
      DO I=2,NS
      I1=I-1
      SOLRD1=SOLRAD+(I1-1)*THICK
      SOLRD2=SOLRD1+THICK
      IF (SRAD.GT.RSMALL.AND.I1.EQ.1) SOLRD1=SOLRAD-SRAD
      IF (I.EQ.NS) THEN
      WRITE(6,'(A,I3,A,F7.2,3A,F8.3,3A,I3,A)') ' XMASK:  shell=',I1,
     & '   from=',SOLRD1, '  to',' infin. A, ',
     & ' volume=',VOLUME(I),'% (',STO,'=',I1,')'
      ELSE
      WRITE(6,'(A,I3,2(A,F7.2),A,F8.3,3A,I3,A)') ' XMASK:  shell=',
     & I1,
     & '   from=',SOLRD1, '  to',SOLRD2,
     & ' A,  volume=',VOLUME(I),'% (',STO,'=',I1,')'
      END IF
      END DO
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume inside mask=',
     &  VINSIDE*R100,'% (',STO,'<=0)'
      WRITE(6,'(A,F10.4,3A)') ' XMASK: volume outside mask=',
     &  VOUTSIDE*R100,'% (',STO,'>0)'
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XFBOX(XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &           NA,NB,NC,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           RRHO,IRHO,RHOMASK,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,
     &           XRCELL,XRTR,XRINTR,
     &           VINSIDE,VOUTSIDE,XRSYGP,XRSYIV,
     &           QAVER)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*), XRTR(3,3), XRINTR(3,3)
      LOGICAL QHERM
      DOUBLE PRECISION VINSIDE, VOUTSIDE
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      LOGICAL QAVER
C local
      INTEGER NN
C pointer
      INTEGER MASK, XSYMASK
C
C allocate space for temporary mask
      NN=(-MAASY+NAASY+1)*
     &   (-MBASY+NBASY+1)*
     &   (-MCASY+NCASY+1)
      MASK=ALLHP(INTEG4(NN))
      IF (QAVER) THEN
      XSYMASK=ALLHP(INTEG4(NN))
      ELSE
      XSYMASK=1
      END IF
C
      CALL XFBOX2(XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &           NA,NB,NC,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           RRHO,IRHO,RHOMASK,HEAP(MASK),
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,
     &           XRCELL,XRTR,XRINTR,
     &           VINSIDE,VOUTSIDE,XRSYGP,XRSYIV,
     &           QAVER,HEAP(XSYMASK))
C
      IF (QAVER) THEN
      CALL FREHP(XSYMASK,INTEG4(NN))
      END IF
C
      CALL FREHP(MASK,INTEG4(NN))
      RETURN
      END
C=====================================================================
      SUBROUTINE XFBOX2(XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX,
     &           NA,NB,NC,
     &           MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           RRHO,IRHO,RHOMASK,MASK,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,
     &           XRCELL,XRTR,XRINTR,
     &           VINSIDE,VOUTSIDE,XRSYGP,XRSYIV,
     &           QAVER,XSYMASK)
C
C Computes a mask for specified box in fractional coordinates
C The mask will be set to values <=0 inside the box and to 1
C outside.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*), XRTR(3,3), XRINTR(3,3)
      LOGICAL QHERM
      DOUBLE PRECISION VINSIDE, VOUTSIDE
      INTEGER XRSYGP(XRMSYM,XRMSYM), XRSYIV(XRMSYM)
      LOGICAL QAVER
      INTEGER XSYMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
C local
      DOUBLE PRECISION FTAS, FTBS, FTCS
      DOUBLE PRECISION ABC(3,3)
      DOUBLE PRECISION ABCINV(3,3)
      DOUBLE PRECISION INSIDE, OUTSIDE
      INTEGER CUSHMC, CUSHNC, CUSHMB, CUSHNB, CUSHMA, CUSHNA
      INTEGER I, J
      INTEGER A, B, C, AX, BX, CX, SYM, AA, BB, CC
      INTEGER AS, BS, CS, XX, YY, ZZ, NABC(3)
      REAL KEY
      LOGICAL COND
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C compute sublattice transformation matrices
      NABC(1)=NA
      NABC(2)=NB
      NABC(3)=NC
      DO I=1,3
      DO J=1,3
      ABC(I,J)=NABC(I)*XRTR(I,J)
      END DO
      END DO
      DO I=1,3
      DO J=1,3
      ABCINV(I,J)=XRINTR(I,J)/NABC(J)
      END DO
      END DO
C
C compute norm of A*, B*, C* in FFT grid box
      FTAS=XRCELL(7)*NA
      FTBS=XRCELL(8)*NB
      FTCS=XRCELL(9)*NC
C initialize mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      MASK(A,B,C)=1
C initialize symmetry and lattice translation info matrix
      IF (QAVER) THEN
      XSYMASK(A,B,C)=-9999
      END IF
      END DO
      END DO
      END DO
C
      CUSHNA=NINT(XMAX*NA)
      CUSHNB=NINT(YMAX*NB)
      CUSHNC=NINT(ZMAX*NC)
      CUSHMA=NINT(XMIN*NA)
      CUSHMB=NINT(YMIN*NB)
      CUSHMC=NINT(ZMIN*NC)
C
      DO CX=CUSHMC,CUSHNC
      DO BX=CUSHMB,CUSHNB
      DO AX=CUSHMA,CUSHNA
C
C map AX, BX, CX into asymmetric unit
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*AX+XRSYMM(SYM,1,2)*BX+XRSYMM(SYM,1,3)*CX
     &       +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*AX+XRSYMM(SYM,2,2)*BX+XRSYMM(SYM,2,3)*CX
     &       +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*AX+XRSYMM(SYM,3,2)*BX+XRSYMM(SYM,3,3)*CX
     &       +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA) +MAASY
      BS=MOD(BB+10000*NB-MBASY,NB) +MBASY
      CS=MOD(CC+10000*NC-MCASY,NC) +MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND)  COND=RHOMASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      SYM=SYM-1
      XX=(AS-AA)/NA
      YY=(BS-BB)/NB
      ZZ=(CS-CC)/NC
      MASK(AS,BS,CS)=0
C define symmetry and lattice translation info
      IF (QAVER) THEN
      CALL MSKM2A(SYM,XX,YY,ZZ,KEY,XRMSYM)
      XSYMASK(AS,BS,CS)=KEY
      END IF
      END IF
      END DO
      END DO
      END DO
C
C initialize RRHO, IRHO
      IF (QHERM) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (MASK(A,B,C).EQ.0) THEN
      IF (QAVER) THEN
      RRHO(A,B,C)=XSYMASK(A,B,C)
      IF (XSYMASK(A,B,C).EQ.-9999) THEN
      CALL WRNDIE(-5,'MASKX2','Internal error no. 6')
      END IF
      ELSE
      RRHO(A,B,C)=ZERO
      END IF
      ELSE
      RRHO(A,B,C)=ONE
      END IF
      END IF
      END DO
      END DO
      END DO
      ELSE
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (MASK(A,B,C).EQ.0) THEN
      IF (QAVER) THEN
      RRHO(A,B,C)=XSYMASK(A,B,C)
      IF (XSYMASK(A,B,C).EQ.-9999) THEN
      CALL WRNDIE(-5,'MASKX2','Internal error no. 7')
      END IF
      ELSE
      RRHO(A,B,C)=ZERO
      END IF
      ELSE
      RRHO(A,B,C)=ONE
      END IF
      IRHO(A,B,C)=ZERO
      END IF
      END DO
      END DO
      END DO
      END IF
C
C
C compute percentage of volume inside and outside of mask
      INSIDE=ZERO
      OUTSIDE=ZERO
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).GT.0) THEN
      IF (ABS(RRHO(A,B,C)-ONE).GE.R4SMAL) THEN
      INSIDE=INSIDE+ONE/RHOMASK(A,B,C)
      ELSEIF (ABS(RRHO(A,B,C)-MASK(A,B,C)).LT.R4SMAL) THEN
      OUTSIDE=OUTSIDE+ONE/RHOMASK(A,B,C)
      END IF
      END IF
      END DO
      END DO
      END DO
C
      VINSIDE=XRNSYM*INSIDE/(NC*NB*NA)
      VOUTSIDE=XRNSYM*OUTSIDE/(NC*NB*NA)
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE MSKM2A(SYM,XTRANS,YTRANS,ZTRANS,KEY,XRMSYM)
C
C given the symmetry operator, lattice translations in x, y and z
C return the remapping key from ASU lattice point to original molecule
C
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C
      INCLUDE 'xmask.inc'
C
      INTEGER SYM,XTRANS,YTRANS,ZTRANS
      REAL KEY
      INTEGER XRMSYM
C local
      LOGICAL ERROR
C
      ERROR=.FALSE.
C
      IF (SYM.LE.0) THEN
         CALL WRNDIE(-5,'MSKM2A',
     &     'Internal error: symmetry operator <= zero')
      ELSE IF (SYM.GT.XRMSYM) THEN
      CALL WRNDIE(-5,'MSKM2A',
     &     'Internal error: symmetry operator > MXSYM')
      END IF
C
      ERROR=(ABS(XTRANS).GT.MXTRAN.OR.
     &       ABS(YTRANS).GT.MXTRAN.OR.
     &       ABS(ZTRANS).GT.MXTRAN)
C
      IF (ERROR) THEN
      WRITE(6,'(2A,/,2A)')
     & ' %MSKM2A-ERR: Lattice point translated by',
     & ' more than four lattice units.',
     & '               Translate molecule closer to the',
     & ' origin (0,0,0) and rerun job.'
      CALL WRNDIE(-5,'MSKM2A','Lattice point too far from origin.')
      END IF
C
C the symmetry operator is incremented by one so that the reserved
C flag -1 can be used for reading/writing molecular averaging masks
      KEY=-(((XTRANS+MXTRAN)*PMX)+((YTRANS+MXTRAN)*PMY)
     &     +((ZTRANS+MXTRAN)*PMZ)+SYM+1)
C
      RETURN
      END
C======================================================================
      SUBROUTINE MSKA2M(KEY,SYM,XTRANS,YTRANS,ZTRANS,XRMSYM)
C
C given the remapping key from ASU lattice point to original molecule
C return the symmetry operator, lattice translations in x, y and z
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C
      INCLUDE 'xmask.inc'
C
      INTEGER SYM,XTRANS,YTRANS,ZTRANS
      REAL KEY
      INTEGER XRMSYM
C local
      INTEGER PACK
      LOGICAL ERROR
C
      ERROR=.FALSE.
C
C the key -1 is reserved, and signifies the identity operator and
C zero lattice translations - needed for reading and writing
C molecular averaging masks
      IF (KEY.EQ.-1.OR.KEY.EQ.0) THEN
         SYM=1
         XTRANS=0
         YTRANS=0
         ZTRANS=0
      ELSE IF (KEY.LT.-1) THEN
         PACK = NINT(-KEY)
         XTRANS = INT((PACK)/PMX)-MXTRAN
         YTRANS = INT((MOD(PACK,PMX))/PMY)-MXTRAN
         ZTRANS = INT((MOD(PACK,PMY))/PMZ)-MXTRAN
         SYM = INT(MOD(PACK,PMZ))-1
      ELSE
         CALL WRNDIE(-5,'MSKA2M',
     &     'Internal error: remapping key >=0')
      END IF
C
      IF (SYM.LE.0) THEN
         CALL WRNDIE(-5,'MSKA2M',
     &     'Internal error: mask returned symmetry operator <= zero')
      ELSE IF (SYM.GT.XRMSYM) THEN
      CALL WRNDIE(-5,'MSKA2M',
     &     'Internal error: mask returned symmetry operator > MXSYM')
      END IF
C
      ERROR=(ABS(XTRANS).GT.MXTRAN.OR.
     &       ABS(YTRANS).GT.MXTRAN.OR.
     &       ABS(ZTRANS).GT.MXTRAN)
C
      IF (ERROR) THEN
      WRITE(6,'(2A,/,2A)')
     & ' %MSKA2M-ERR: Lattice point translated by',
     & ' more than four lattice units.',
     & '               Translate molecule closer to the',
     & ' origin (0,0,0) and rerun job.'
      CALL WRNDIE(-5,'MSKA2M','Lattice point too far from origin.')
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE MSKDIM(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     XRNSYM,XRSYTH,XRSYMM,XRITSY,XRMSYM,XRTR,XRINTR,RHOMASK,MASK,
     &     RTMAT,MATRIX,CENTER,CUSHION,AMIN,AMAX,BMIN,BMAX,CMIN,CMAX)
C
C given a mask and rotation matrix R and translation T
C return the limits of the grid space that will cover the transformed
C mask
C
C Author: Paul Adams
C
      IMPLICIT NONE
      INCLUDE 'xmask.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
C I/O
      INTEGER NA, NB, NC
      INTEGER MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER XRNSYM, XRSYTH, XRMSYM
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      DOUBLE PRECISION RTMAT(XRNSYM,3,4),XRTR(3,3),XRINTR(3,3)
      DOUBLE PRECISION MATRIX(MXMAVG,3,4),CENTER(3)
      DOUBLE PRECISION CUSHION
      INTEGER AMIN,AMAX,BMIN,BMAX,CMIN,CMAX
C local
      INTEGER A,B,C,AN1,AN2,BN1,BN2,CN1,CN2,POINTS
      DOUBLE PRECISION AS,BS,CS,AP,BP,CP
      INTEGER LAMIN,LAMAX,LBMIN,LBMAX,LCMIN,LCMAX
      INTEGER MSO,LTX,LTY,LTZ
      INTEGER ACUSH,BCUSH,CCUSH
      DOUBLE PRECISION CENA,CENB,CENC
      DOUBLE PRECISION ATRANS,BTRANS,CTRANS
C parameter
      INTEGER IZERO,ONE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0,IZERO=0,ONE=1)
C
C initialise
      ACUSH=NINT((CUSHION*XRTR(1,1))*DBLE(NA))
      BCUSH=NINT((CUSHION*XRTR(2,2))*DBLE(NB))
      CCUSH=NINT((CUSHION*XRTR(3,3))*DBLE(NC))
      CENA=CENTER(1)*XRTR(1,1)*NA
      CENB=CENTER(2)*XRTR(2,2)*NB
      CENC=CENTER(3)*XRTR(3,3)*NC
      AMIN=9999
      AMAX=-9999
      BMIN=9999
      BMAX=-9999
      CMIN=9999
      CMAX=-9999
C
C set up the realspace transformation matrices
      CALL PREMAT(NA,NB,NC,XRNSYM,XRINTR,MATRIX,ONE,
     &            XRTR,RTMAT,'NONI',XRNSYM,XRMSYM,
     &            XRSYTH,XRSYMM,XRITSY)
C
      POINTS=0
C
C loop over points in primary mask
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
C
      IF (RHOMASK(A,B,C).GT.IZERO.AND.MASK(A,B,C).LE.ZERO) THEN
C
      POINTS=POINTS+1
C
C get the symmetry operation which was applied to map this grid point
C into the ASU
      CALL MSKA2M(MASK(A,B,C),MSO,LTX,LTY,LTZ,XRMSYM)
C
C map point from ASU to original lattice point (ie around molecule)
      ATRANS=A-(LTX*NA)-(XRSYMM(MSO,1,4)*NA)/XRSYTH
      BTRANS=B-(LTY*NB)-(XRSYMM(MSO,2,4)*NB)/XRSYTH
      CTRANS=C-(LTZ*NC)-(XRSYMM(MSO,3,4)*NC)/XRSYTH
      AS=XRITSY(MSO,1,1)*ATRANS+XRITSY(MSO,2,1)*BTRANS
     &     +XRITSY(MSO,3,1)*CTRANS
      BS=XRITSY(MSO,1,2)*ATRANS+XRITSY(MSO,2,2)*BTRANS
     &     +XRITSY(MSO,3,2)*CTRANS
      CS=XRITSY(MSO,1,3)*ATRANS+XRITSY(MSO,2,3)*BTRANS
     &     +XRITSY(MSO,3,3)*CTRANS
C
C apply the centre of rotation correction
      ATRANS=AS-CENA
      BTRANS=BS-CENB
      CTRANS=CS-CENC
C
C apply the transformation to off-grid space
      AP=RTMAT(1,1,1)*ATRANS+RTMAT(1,1,2)*BTRANS
     &  +RTMAT(1,1,3)*CTRANS+RTMAT(1,1,4)
      BP=RTMAT(1,2,1)*ATRANS+RTMAT(1,2,2)*BTRANS
     &  +RTMAT(1,2,3)*CTRANS+RTMAT(1,2,4)
      CP=RTMAT(1,3,1)*ATRANS+RTMAT(1,3,2)*BTRANS
     &  +RTMAT(1,3,3)*CTRANS+RTMAT(1,3,4)
C
C apply center of rotation correction
      AP=AP+CENA
      BP=BP+CENB
      CP=CP+CENC
C
C the limits of the nearest grid points
      AN1=INT(AP)
      AN2=AN1+INT(SIGN(1.0D0,AP))
      BN1=INT(BP)
      BN2=BN1+INT(SIGN(1.0D0,BP))
      CN1=INT(CP)
      CN2=CN1+INT(SIGN(1.0D0,CP))
C
C find the min and max of these limits
      LAMIN=MIN(AN1,AN2)
      LAMAX=MAX(AN1,AN2)
      LBMIN=MIN(BN1,BN2)
      LBMAX=MAX(BN1,BN2)
      LCMIN=MIN(CN1,CN2)
      LCMAX=MAX(CN1,CN2)
C
C if these are the best so far store them
      IF (LAMIN.LT.AMIN) AMIN=LAMIN
      IF (LAMAX.GT.AMAX) AMAX=LAMAX
      IF (LBMIN.LT.BMIN) BMIN=LBMIN
      IF (LBMAX.GT.BMAX) BMAX=LBMAX
      IF (LCMIN.LT.CMIN) CMIN=LCMIN
      IF (LCMAX.GT.CMAX) CMAX=LCMAX
C
      END IF
C
      END DO
      END DO
      END DO
C
C check to see if we actually visited any points
      IF (POINTS.LE.0) THEN
        CALL WRNDIE(-5,'MSKDIM',
     &     'no points in mask - check mask')
      END IF
C
C add the cushion around the molecule
      AMIN=AMIN-ACUSH
      AMAX=AMAX+ACUSH
      BMIN=BMIN-BCUSH
      BMAX=BMAX+BCUSH
      CMIN=CMIN-CCUSH
      CMAX=CMAX+CCUSH
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE PREMAT(NA,NB,NC,NSYM,A,R,MOL,B,MATRIX,INVFLG,
     &                  XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C
C Author: Paul Adams
C
C precompute the rotation matrix and translation vector to be applied
C to a set of grid points
C
C must account for the symmetry operator which was originally used to
C map the grid points into the asymmetric unit (this is stored in the
C mask). This requires that we calculate matrices and translations for
C each symmetry operator. We use the inverse of the symmetry operator
C and its translational vector.
C
C for both the non-inverse and inverse case will compute:
C
C     Rotation = G*B*Rr*A*Rs*F
C
C for the non-inverse case will compute:
C
C     Translation = G*B*Rr*A*Rs*Ts + G*B*Tr
C
C for the inverse case will compute:
C
C     Translation = G*B*Rr*A*Rs*Ts + G*B*Rr*Tr
C
C where F is the conversion from grid space to fractional (diagonal)
C       A is the orthogonalisation matrix
C       Rr is the real space rotation matrix
C       Tr is the real space translation vector
C       Rs is the inverse of the symmetry operator rotation matrix
C       Ts is the inverse of the symmetry operator translation vector
C       B is the fractionalisation matrix
C       G is the conversion to grid space from fractional (diagonal)
C
C in reality G and F are NA,NB,NC and 1/NA,1/NB,1/NC respectively
C Rr and Tr are both stored in R(MOL,..,..)
C Rs is taken from the transpose of XRITSY
C Ts is taken from XRSYMM
C the resulting rotation matrix and translation are both stored in MATRIX
C
C N.B. if the inverse flag is set then the inverse of Rr and Tr
C      are used
C
C I/O
      IMPLICIT NONE
C
      INCLUDE 'xmask.inc'
      INCLUDE 'timer.inc'
C
      CHARACTER*4 INVFLG
      INTEGER NA,NB,NC,NSYM,MOL
      DOUBLE PRECISION A(3,3),R(MXMAVG,3,4),B(3,3)
      DOUBLE PRECISION MATRIX(NSYM,3,4)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER I,J,K,S
      DOUBLE PRECISION F(3),G(3),DET
      DOUBLE PRECISION TEMP1(3,3),TEMP2(3,3),TRANS1(3),TRANS2(3)
      DOUBLE PRECISION RSYTH
C parameter
      DOUBLE PRECISION ZERO,ONE,SMALL
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-3)
C
      RSYTH=XRSYTH
C
      F(1)=ONE/NA
      F(2)=ONE/NB
      F(3)=ONE/NC
      G(1)=NA
      G(2)=NB
      G(3)=NC
C
C calculate the matrix product
C if inverse then use the transpose of R
C
      DO S=1,NSYM
C
C G*B
      DO I=1,3
         DO J=1,3
            TEMP1(I,J)=G(I)*B(I,J)
         END DO
      END DO
C
C if not the inverse real space operator then:
C G*B*Tr
      IF (INVFLG.NE.'INVE') THEN
         DO I=1,3
            TRANS1(I)=ZERO
            DO K=1,3
               TRANS1(I)=TRANS1(I)+(TEMP1(I,K)*R(MOL,K,4))
            END DO
         END DO
      END IF
C
C G*B*Rr
C if inverse use the transpose of real space rotation matrix
      DO I=1,3
         DO J=1,3
            TEMP2(I,J)=ZERO
            DO K=1,3
               IF (INVFLG.EQ.'INVE') THEN
                  TEMP2(I,J)=TEMP2(I,J)+(TEMP1(I,K)*R(MOL,J,K))
               ELSE
                  TEMP2(I,J)=TEMP2(I,J)+(TEMP1(I,K)*R(MOL,K,J))
               END IF
            END DO
         END DO
      END DO
C
C if inverse real space operator then:
C G*B*Rr*Tr
C use the inverse (negative) of real space translation
      IF (INVFLG.EQ.'INVE') THEN
         DO I=1,3
            TRANS1(I)=ZERO
            DO K=1,3
               TRANS1(I)=TRANS1(I)+(TEMP2(I,K)*(-R(MOL,K,4)))
            END DO
         END DO
      END IF
C
C G*B*Rr*A
      DO I=1,3
         DO J=1,3
            TEMP1(I,J)=ZERO
            DO K=1,3
               TEMP1(I,J)=TEMP1(I,J)+(TEMP2(I,K)*A(K,J))
            END DO
         END DO
      END DO
C
C G*B*Rr*A*Rs
C note transpose of inverse of symmetry rotation matrix (is stored as
C transpose)
      DO I=1,3
         DO J=1,3
            TEMP2(I,J)=ZERO
            DO K=1,3
               TEMP2(I,J)=TEMP2(I,J)+(TEMP1(I,K)*XRITSY(S,J,K))
            END DO
         END DO
      END DO
C
C G*B*Rr*A*Rs*Ts
C note inverse of symmetry translation component
      DO I=1,3
         TRANS2(I)=ZERO
         DO K=1,3
            TRANS2(I)=TRANS2(I)+(TEMP2(I,K)*(-XRSYMM(S,K,4)/RSYTH))
         END DO
      END DO
C
C G*B*Rr*A*Rs*F
      DO I=1,3
         DO J=1,3
            MATRIX(S,I,J)=TEMP2(I,J)*F(J)
         END DO
      END DO
C
C sum translational components and insert into MATRIX array
      DO I=1,3
         MATRIX(S,I,4)=TRANS1(I)+TRANS2(I)
      END DO
C
C check determinant of matrix
      DET=(MATRIX(S,1,1)*MATRIX(S,2,2)-MATRIX(S,1,2)
     &    *MATRIX(S,2,1))*MATRIX(S,3,3)
     &   +(MATRIX(S,2,1)*MATRIX(S,3,2)-MATRIX(S,2,2)
     &    *MATRIX(S,3,1))*MATRIX(S,1,3)
     &   +(MATRIX(S,3,1)*MATRIX(S,1,2)-MATRIX(S,3,2)
     &    *MATRIX(S,1,1))*MATRIX(S,2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      WRITE(6,'(A,F8.5)') ' %XMAVRT-ERR: invalid determinant =',DET
      CALL WRNDIE(-5,'MAPAVER',
     &            'cannot precompute transformation matrix.')
      END IF
C
C if debugging is set high enough, print out the precomputed matrices
      IF (WRNLEV.GT.10) THEN
      WRITE(6,'(2A)') ' |-----------------------',
     &     '--------------------------------|'
      WRITE(6,'(A,3F10.4,A)')
     & ' | Matrix = ( ',(MATRIX(S,1,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' |            ',(MATRIX(S,2,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' |            ',(MATRIX(S,3,J),J=1,3),' )           |'
      WRITE(6,'(A,3F10.4,A)')
     & ' | Vector = ( ',(MATRIX(S,J,4),J=1,3),' )           |'
      WRITE(6,'(2A)') ' |-----------------------',
     &     '--------------------------------|'
      ENDIF
C
      END DO
C
      RETURN
      END
C
C======================================================================
C
