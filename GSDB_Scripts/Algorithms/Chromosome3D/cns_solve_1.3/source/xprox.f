      SUBROUTINE XPROX(NA,NB,NC,NAP,NBP,NCPL,NAPP,NBPP,NCPP,MAASY,
     &  MBASY,MCASY,NAASY,NBASY,NCASY,XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &  XRITSY,QHERM,XRCELL,ADEPTH,ARPNMX,ARPNN,ARPNX,ARPN,ARPNL,
     &  ARPNDB,ARPNMLT,MAPR,XRHONUM,XRHONAM,HPRRHO,HPIRHO,HPRHOMA,NRHO,
     &  IRHO,NMASK,XRMAP,XRTR,XRINTR,XRVOL,XRMREF,XRNREF,
     &  HPH,HPK,HPL,
     &  XSFNUM,XSFNAM,XSFTYPE,HPSF,
     &  HPMULT,HPTYPE,XRSYGP,XRSYIV)
C
C Creates distance and property maps using selected atomic positions.
C
C Front-end routine for XPROX2
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
      CALL XPROX2(HEAP(FLAGS),HEAP(XRVDW),HEAP(XL),
     &   HEAP(YL),HEAP(ZL),NA,NB,NC,NAP,NBP,
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
      SUBROUTINE XPROX2(FLAGS,XRVDW,XL,
     &   YL,ZL,NA,NB,NC,NAP,NBP,
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
C Parsing routine for distance and property maps.
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
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INTEGER FLAGS(*)
      DOUBLE PRECISION XRVDW(*), XL(*), YL(*), ZL(*)
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
      INTEGER XRSYGP(XRMSYM,XRMSYM),XRSYIV(XRMSYM)
C local
      INTEGER IAT, LNATOM, LLNATOM, INSYM, I
      CHARACTER*8 SFROM
      INTEGER ISTO, ISTO2, LENG, LENG2
      DOUBLE PRECISION SECS, SECS2, CUTOFF
      LOGICAL COND, ERR
C pointer
      INTEGER HPRMAP, HPIMAP, HPRMAP2, HPIMAP2
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, SIX, R100, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, SIX=6.0D0)
      PARAMETER (R100=100.D0, FOUR=4.0D0)
C begin
C
C default values
      ERR=.FALSE.
      ISTO=1
      LENG=LEN(XRHONAM(ISTO))
      CALL TRIMM(XRHONAM(ISTO),LENG)
      ISTO2=2
      LENG2=LEN(XRHONAM(ISTO2))
      CALL TRIMM(XRHONAM(ISTO2),LENG2)
C
      CUTOFF=FOUR
C
      SFROM='WMAIN'
C
C default selection: none
      LNATOM=0
C
      CALL PUSEND('PROXimity>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PROXimity>')
      CALL MISCOM('PROXimity>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-xray-proximity')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(FLAGS,LNATOM,X,Y,Z,.TRUE.)
      CALL MAKIND(FLAGS,NATOM,LNATOM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CUTO') THEN
      CALL NEXTF('CUTOff=',CUTOFF)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FROM') THEN
      CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('FROM=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'FROM=',SFROM
      ELSE
      SFROM=WD(1:WDLEN)
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DIST') THEN
      CALL NEXTWD('DISTance_map=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('DISTance_map=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'DISTance_map=',XRHONAM(ISTO)
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
      WRITE(6,'(3A)') ' %XPROX-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XPROX','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'PROP') THEN
      CALL NEXTWD('PROPerty_map=')
      IF (WD(1:4).EQ.'=   ') CALL NEXTWD('PROPerty_map=')
      IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(1X,2A)') 'PROPerty_map=',XRHONAM(ISTO2)
      ELSE
      COND=.FALSE.
      DO I=1,XRHONUM
      IF (WD(1:WDLEN).EQ.XRHONAM(I)) THEN
      COND=.TRUE.
      ISTO2=I
      LENG2=LEN(XRHONAM(ISTO2))
      CALL TRIMM(XRHONAM(ISTO2),LENG2)
      END IF
      END DO
      IF (.NOT.COND) THEN
      WRITE(6,'(3A)') ' %XPROX-ERR: real space object ',
     & WD(1:WDLEN),' not declared.'
      CALL WRNDIE(-5,'XPROX','object undeclared.')
      ERR=.TRUE.
      END IF
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)') ' |-----------------------',
     & 'proximity---------------------------------------------'
      WRITE(6,'(2A)')
     & ' |  DISTance_map = ',XRHONAM(ISTO)(1:LENG)
      WRITE(6,'(2A)')
     & ' |  PROPerty_map = ',XRHONAM(ISTO2)(1:LENG2)
      WRITE(6,'(2A)')
     & ' |  FROM = ',SFROM
      WRITE(6,'(A,F10.4)') ' |  CUTOff=',CUTOFF
      WRITE(6,'(2A)') ' |--------------------------',
     & '---------------------------------------------------'
      ELSE
      CALL CHKEND('PROXimity>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (LNATOM.EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','zero atoms selected.')
      ELSEIF (.NOT.ERR) THEN
C
C gather coordinates and properties for all selected atoms
      DO IAT=1,LNATOM
      XL(IAT)=X(FLAGS(IAT))
      YL(IAT)=Y(FLAGS(IAT))
      ZL(IAT)=Z(FLAGS(IAT))
      END DO
C
      IF (SFROM(1:4).EQ.'X   ') THEN
      CALL XPROMAP(FLAGS,X,LNATOM,XRVDW)
      ELSEIF (SFROM(1:4).EQ.'Y   ') THEN
      CALL XPROMAP(FLAGS,Y,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'Z   ') THEN
      CALL XPROMAP(FLAGS,Z,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'WMAI') THEN
      CALL XPROMAP(FLAGS,WMAIN,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'QMAI') THEN
      CALL XPROMAP(FLAGS,QMAIN,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'B   ') THEN
      CALL XPROMAP(FLAGS,WMAIN,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'Q   ') THEN
      CALL XPROMAP(FLAGS,QMAIN,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'XCOM') THEN
      CALL XPROMAP(FLAGS,XCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'YCOM') THEN
      CALL XPROMAP(FLAGS,YCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'ZCOM') THEN
      CALL XPROMAP(FLAGS,ZCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'WCOM') THEN
      CALL XPROMAP(FLAGS,WCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'BCOM') THEN
      CALL XPROMAP(FLAGS,WCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'QCOM') THEN
      CALL XPROMAP(FLAGS,QCOMP,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'REFX') THEN
      CALL XPROMAP(FLAGS,REFX,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'REFY') THEN
      CALL XPROMAP(FLAGS,REFY,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'REFZ') THEN
      CALL XPROMAP(FLAGS,REFZ,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'RMSD') THEN
      CALL XPROMAP(FLAGS,RMSD,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'HARM') THEN
      CALL XPROMAP(FLAGS,KCNSTR,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'MASS') THEN
      CALL XPROMAP(FLAGS,AMASS,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'CHAR') THEN
      CALL XPROMAP(FLAGS,CG,LNATOM,XRVDW)
      ELSE IF (SFROM(1:6).EQ.'STORE1') THEN
      IF (PTRSTO(1).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE1 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(1)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE2') THEN
      IF (PTRSTO(2).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE2 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(2)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE3') THEN
      IF (PTRSTO(3).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE3 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(3)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE4') THEN
      IF (PTRSTO(4).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE4 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(4)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE5') THEN
      IF (PTRSTO(5).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE5 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(5)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE6') THEN
      IF (PTRSTO(6).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE6 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(6)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE7') THEN
      IF (PTRSTO(7).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE7 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(7)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE8') THEN
      IF (PTRSTO(8).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE8 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(8)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:6).EQ.'STORE9') THEN
      IF (PTRSTO(9).EQ.0) THEN
      CALL WRNDIE(-5,'XPROX','STORE9 undefined')
      ELSE
      CALL XPROMAP(FLAGS,HEAP(PTRSTO(9)),LNATOM,XRVDW)
      END IF
      ELSE IF (SFROM(1:4).EQ.'DX  ') THEN
      CALL XPROMAP(FLAGS,DX,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'DY  ') THEN
      CALL XPROMAP(FLAGS,DY,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'DZ  ') THEN
      CALL XPROMAP(FLAGS,DZ,LNATOM,XRVDW)
      ELSE IF (SFROM(1:4).EQ.'VX  ') THEN
      DO I=1,NATOM
      XV(I)=XV(I)/TIMFAC
      END DO
      CALL XPROMAP(FLAGS,XV,LNATOM,XRVDW)
      DO I=1,NATOM
      XV(I)=XV(I)*TIMFAC
      END DO
      ELSE IF (SFROM(1:4).EQ.'VY  ') THEN
      DO I=1,NATOM
      YV(I)=YV(I)/TIMFAC
      END DO
      CALL XPROMAP(FLAGS,YV,LNATOM,XRVDW)
      DO I=1,NATOM
      YV(I)=YV(I)*TIMFAC
      END DO
      ELSE IF (SFROM(1:4).EQ.'VZ  ') THEN
      DO I=1,NATOM
      ZV(I)=ZV(I)/TIMFAC
      END DO
      CALL XPROMAP(FLAGS,ZV,LNATOM,XRVDW)
      DO I=1,NATOM
      ZV(I)=ZV(I)*TIMFAC
      END DO
      ELSE IF (SFROM(1:4).EQ.'FBET') THEN
      CALL XPROMAP(FLAGS,FBETA,LNATOM,XRVDW)
      ELSE
      CALL WRNDIE(-5,'XPROX',' unknown atom property.')
      END IF
C
      IF (LNATOM.LE.0) THEN
      CALL WRNDIE(-5,'XPROX','Zero atoms selected.')
      ELSE
      WRITE(6,'(A,I8,A)') ' XPROX: ',LNATOM,
     &  ' atoms have been selected for proximity calculation.'
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
      WRITE(6,'(A)') ' XPROX: applying non-crystallographic symmetry.'
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
      IF (HPRRHO(ISTO).EQ.0) THEN
      CALL XMAPAL(HPRRHO(ISTO),HPIRHO(ISTO),QHERM,NRHO,IRHO)
      END IF
      HPRMAP=HPRRHO(ISTO)
      HPIMAP=HPIRHO(ISTO)
C
      IF (HPRRHO(ISTO2).EQ.0) THEN
      CALL XMAPAL(HPRRHO(ISTO2),HPIRHO(ISTO2),QHERM,NRHO,IRHO)
      END IF
      HPRMAP2=HPRRHO(ISTO2)
      HPIMAP2=HPIRHO(ISTO2)
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      CALL XPROXX(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   HEAP(HPRMAP),HEAP(HPIMAP),
     &   HEAP(HPRMAP2),HEAP(HPIMAP2),HEAP(HPRHOMA),
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,LNATOM,XL,YL,ZL,XRVDW,CUTOFF)
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' XPROX2: CPU-time:  proximity-calc=',SECS2-SECS
      END IF
C
      END IF
      END IF
C
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XPROXX(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &   RRHO,IRHO,RRHO2,IRHO2,RHOMASK,
     &   XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,XRCELL,XRTR,XRINTR,
     &   QHERM,NATOM,X,Y,Z,XRVDW,CUTOFF)
C
C Front-end routine for XPROXX2.  Determines the size
C of the cushion to surround the asymmetric unit.
C Atoms are mapped into the asymmetric unit.
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
      REAL RRHO2(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO2(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION XRCELL(*), XRTR(3,3), XRINTR(3,3)
      LOGICAL QHERM
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XRVDW(*), CUTOFF
C local
      INTEGER CUSHNA, CUSHNB, CUSHNC, CUSHMA, CUSHMB, CUSHMC
      INTEGER IAT, NN, I, J, SYM, AS, BS, CS, AA, BB, CC
      LOGICAL COND
      DOUBLE PRECISION DX, DY, DZ
      INTEGER AX, BX, CX, NABC(3)
      DOUBLE PRECISION ABC(3,3), ABCINV(3,3), FTAS, FTBS, FTCS
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C pointer
      INTEGER DIST, PROP
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
C
C initialize extent of cushion
      CUSHNA=NAASY
      CUSHNB=NBASY
      CUSHNC=NCASY
      CUSHMA=MAASY
      CUSHMB=MBASY
      CUSHMC=MCASY
C
      DO IAT=1,NATOM
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
      DX=X(IAT)
      DY=Y(IAT)
      DZ=Z(IAT)
C
C copy index
      AX=NINT(DX)
      BX=NINT(DY)
      CX=NINT(DZ)
      ELSE
      CALL WRNDIE(-5,'XPROXX','Internal Error No. 3')
      END IF
C
C determine the difference between lattice point AX, AY, AZ and
C actual position
      DX=AX-DX
      DY=BX-DY
      DZ=CX-DZ
C
      CUSHNA=MAX(CUSHNA,INT(CUTOFF*FTAS-DX+RSMALL)+AX)
      CUSHNB=MAX(CUSHNB,INT(CUTOFF*FTBS-DY+RSMALL)+BX)
      CUSHNC=MAX(CUSHNC,INT(CUTOFF*FTCS-DZ+RSMALL)+CX)
      CUSHMA=MIN(CUSHMA,-INT(CUTOFF*FTAS+DX+RSMALL)+AX)
      CUSHMB=MIN(CUSHMB,-INT(CUTOFF*FTBS+DY+RSMALL)+BX)
      CUSHMC=MIN(CUSHMC,-INT(CUTOFF*FTCS+DZ+RSMALL)+CX)
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
C allocate space for temporary distance and property matrices
      NN=(-CUSHMA+CUSHNA+1)*
     &   (-CUSHMB+CUSHNB+1)*
     &   (-CUSHMC+CUSHNC+1)
      DIST=ALLHP(IREAL4(NN))
      PROP=ALLHP(IREAL4(NN))
C
      CALL XPROXX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RRHO2,IRHO2,RHOMASK,HEAP(DIST),HEAP(PROP),
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,ABC,
     &           ABCINV,NATOM,X,Y,Z,XRVDW,FTAS,FTBS,FTCS,CUTOFF)
C
      CALL FREHP(PROP,IREAL4(NN))
      CALL FREHP(DIST,IREAL4(NN))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XPROXX2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &           CUSHMC,CUSHNC,CUSHMB,CUSHNB,CUSHMA,CUSHNA,
     &           RRHO,IRHO,RRHO2,IRHO2,RHOMASK,DIST,PROP,
     &           XRNSYM,XRMSYM,XRSYTH,XRSYMM,
     &           XRITSY,QHERM,ABC,
     &           ABCINV,NATOM,X,Y,Z,XRVDW,FTAS,FTBS,FTCS,CUTOFF)
C
C Computes distance and property maps of selected atoms
C specified by coordinates
C (X,Y,Z, NATOM) on a lattice given by i*A + j*B + k*C
C where i,j,k are integers and (A,B,C)=ABCINV.
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
      REAL RRHO2(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO2(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL DIST(CUSHMA:CUSHNA,CUSHMB:CUSHNB,CUSHMC:CUSHNC)
      REAL PROP(CUSHMA:CUSHNA,CUSHMB:CUSHNB,CUSHMC:CUSHNC)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION ABC(3,3)
      DOUBLE PRECISION ABCINV(3,3)
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XRVDW(*)
      DOUBLE PRECISION FTAS, FTBS, FTCS, CUTOFF
C local
      INTEGER A, B, C, AX, BX, CX, IAT, AA, BB, CC
      INTEGER AS, BS, CS
      DOUBLE PRECISION RCUT2, DX, DY, DZ, DELX, DELY, DELZ
      DOUBLE PRECISION DELXB, DELYB, DELZB, DELXC, DELYC, DELZC
      INTEGER AMAXP, AMAXN, BMAXP, BMAXN, CMAXP, CMAXN, SYM
      LOGICAL COND
      DOUBLE PRECISION TEMP
C parameters
      DOUBLE PRECISION ZERO, ONE
      REAL SZERO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, SZERO=0.0)
C begin
C
      RCUT2=CUTOFF**2
C
C initialize mask
      DO C=CUSHMC,CUSHNC
      DO B=CUSHMB,CUSHNB
      DO A=CUSHMA,CUSHNA
      DIST(A,B,C)=R4BIG
      PROP(A,B,C)=ZERO
      END DO
      END DO
      END DO
C
C loop over all atoms
      DO IAT=1,NATOM
C
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
      AMAXP=INT(CUTOFF*FTAS-DX+RSMALL)
      BMAXP=INT(CUTOFF*FTBS-DY+RSMALL)
      CMAXP=INT(CUTOFF*FTCS-DZ+RSMALL)
      AMAXN=INT(CUTOFF*FTAS+DX+RSMALL)
      BMAXN=INT(CUTOFF*FTBS+DY+RSMALL)
      CMAXN=INT(CUTOFF*FTCS+DZ+RSMALL)
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
      CALL WRNDIE(-5,'XPROXX2','Internal Error No. 1')
      END IF
C
      TEMP=(ABCINV(1,1)*A+DELXB)**2
     &    +(ABCINV(2,1)*A+DELYB)**2
     &    +(ABCINV(3,1)*A+DELZB)**2
      IF (RCUT2.GT.TEMP) THEN
C
C distance is within cutoff
      IF (DIST(AA,BB,CC).GT.TEMP) THEN
C
C we've found an atom that is closer to grid point
C AA, BB, CC
      DIST(AA,BB,CC)=TEMP
      PROP(AA,BB,CC)=XRVDW(IAT)
      END IF
      END IF
C
      END DO
C
      END DO
      END DO
C
      END DO
C
C map all grid points outside asymmetric unit into
C asymmetric unit
C =================================================
      DO CC=CUSHMC,CUSHNC
      DO BB=CUSHMB,CUSHNB
      DO AA=CUSHMA,CUSHNA
      IF (DIST(AA,BB,CC).LT.R4BIG) THEN
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
      IF (COND) THEN
      IF (DIST(AS,BS,CS).GT.DIST(AA,BB,CC)) THEN
C we've found an atom that is closer to grid point
C AS, BS, CS
      DIST(AS,BS,CS)=DIST(AA,BB,CC)
      PROP(AS,BS,CS)=PROP(AA,BB,CC)
      END IF
      END IF
      SYM=SYM+1
      END DO
C
      IF (.NOT.COND) CALL WRNDIE(-5,'XPROXX2','Internal error no. 2')
C
      END IF
      END IF
      END DO
      END DO
      END DO
C
C
C now we have to map the distance and property maps into RHO and RHO2
C ===================================================================
C
C compute the square root to get the distance.
C
C initialize RRHO, IRHO
      IF (QHERM) THEN
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (DIST(A,B,C).LT.R4BIG) THEN
      RRHO(A,B,C)=SQRT(MAX(SZERO,DIST(A,B,C)))
      RRHO2(A,B,C)=PROP(A,B,C)
      ELSE
      RRHO(A,B,C)=ZERO
      RRHO2(A,B,C)=ZERO
      END IF
      END DO
      END DO
      END DO
      ELSE
      DO C=MCASY,NCASY
      DO B=MBASY,NBASY
      DO A=MAASY,NAASY
      IF (DIST(A,B,C).LT.R4BIG) THEN
      RRHO(A,B,C)=SQRT(MAX(SZERO,DIST(A,B,C)))
      RRHO2(A,B,C)=PROP(A,B,C)
      ELSE
      RRHO(A,B,C)=ZERO
      RRHO2(A,B,C)=ZERO
      END IF
      IRHO(A,B,C)=ZERO
      IRHO2(A,B,C)=ZERO
      END DO
      END DO
      END DO
      END IF
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XPROMAP(INVMAP,ARRAY,DIM,WORK)
C
C maps ARRAY into WORK using INVMAP
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INTEGER DIM, INVMAP(*)
      DOUBLE PRECISION ARRAY(*), WORK(*)
C local
      INTEGER I
C begin
      DO I=1,DIM
      WORK(I)=ARRAY(INVMAP(I))
      END DO
      RETURN
      END
