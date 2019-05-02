C
C=====================================================================
C
      SUBROUTINE CALIBR(ISLCT,CALPTR,DIAGON)
C
C Set up distance calibration.
C
C Author: Michael Nilges, EMBL, 1994.
C
      IMPLICIT NONE
C
C global
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'noe.inc'
C
C input/ output
C   ISLCT,           flags for selected atoms
      INTEGER ISLCT(*),CALPTR(*)
      DOUBLE PRECISION DIAGON(*)
C
C local
      DOUBLE PRECISION RTEMP
      INTEGER NISLCT, ICAL, JCAL, I, J
      LOGICAL QMATCH,QCALIB
      DOUBLE PRECISION CALSHU,CALSHL
      CHARACTER*4 CALMOD,ACTION
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      QCALIB=.TRUE.
      CALMOD='ANIS'
      CALSHL=ZERO
      CALSHU=9999.0D0
      ACTION='DEFA'
      CALL PUSEND('CALIBR>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('CALIBR>')
C
C[
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-aria-calibrate')
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALCNT=0
      CALL FILL4(CALPTR,NATOM,0)
      CALL FILLR8(DIAGON, NATOM, ONE)
      CALEXP=SIX
      CALDIS=SIX
C
      ELSE IF (WD(1:4).EQ.'DIAGON') THEN
      CALL NEXTF('DIAGONal=',RTEMP)
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      WRITE (6,'(A)') ' %CALIBR-ERR: method not implemented. no effect.'
C[
      IF ((NISLCT.GT.0).AND.(RTEMP.GT.RSMALL)) THEN
      CALL MAKIND(ISLCT,NATOM,NISLCT)
      DO I=1, NISLCT
      DIAGON(ISLCT(I))=RTEMP
      END DO
      END IF
C]

C
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      CALL NEXTA4('Calibration-mode=',CALMOD)
C
      ELSE IF (WD(1:4).EQ.'EXPO') THEN
      CALL NEXTF('Calibration-exponent=',CALEXP)
C
      ELSE IF (WD(1:4).EQ.'CUTO') THEN
      CALL NEXTF('Calibration-cutoff=',CALDIS)
C
      ELSE IF (WD(1:4).EQ.'SHEL') THEN
      CALL NEXTF('Lower_distance=',CALSHL)
      CALL NEXTF('Upper_distance=',CALSHU)
C
      ELSE IF (WD(1:4).EQ.'VECT') THEN
C select two groups of atoms
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL CHKSEL
     &  (ISLCT,NISLCT,NATOM,CALPTR,CALCNT,ICAL,MAXCAL,QMATCH)
C if a new group has been created initialize the arrays
C[
      IF (ICAL.GT.0) THEN
C[
      IF (.NOT.QMATCH) THEN
      DO I=1,CALCNT
      DREFER(ICAL,I)=ZERO
      DREFER(I,ICAL)=ZERO
      VREFER(ICAL,I)=ZERO
      VREFER(I,ICAL)=ZERO
      END DO
      END IF
C]
      ELSE
      WRITE (6,'(A)') ' CALIBR: ill-formed group pointer'
      END IF
C]
      CALL SELCTA(ISLCT,NISLCT,X,Y,Z,.TRUE.)
      CALL CHKSEL
     &  (ISLCT,NISLCT,NATOM,CALPTR,CALCNT,JCAL,MAXCAL,QMATCH)
C if a new group has been created initialize the arrays
C[
      IF (JCAL.GT.0) THEN
C[
      IF (.NOT.QMATCH) THEN
      DO J=1,CALCNT
      DREFER(JCAL,J)=ZERO
      DREFER(J,JCAL)=ZERO
      VREFER(JCAL,J)=ZERO
      VREFER(J,JCAL)=ZERO
      END DO
      END IF
C]
      ELSE
      WRITE (6,'(A)') ' CALIBR: ill-formed group pointer'
      END IF
C]
C
      ELSE IF (WD(1:4).EQ.'DREF') THEN
      CALL NEXTF('Reference-distance=',RTEMP)
C[
      IF ((ICAL.GT.0).AND.(JCAL.GT.0)) THEN
C[
      IF (RTEMP.GT.RSMALL) THEN
      RTEMP=RTEMP**(-CALEXP)
      END IF
C]
      DREFER(ICAL,JCAL)=RTEMP
      DREFER(JCAL,ICAL)=RTEMP
      ELSE
      WRITE (6,'(A)') ' CALIBR: ill-formed group pointer'
      END IF
C]
      ELSE IF (WD(1:4).EQ.'VREF') THEN
      CALL NEXTF('Reference-intensity=',RTEMP)
C[
      IF ((ICAL.GT.0).AND.(JCAL.GT.0)) THEN
      VREFER(ICAL,JCAL)=RTEMP
      VREFER(JCAL,ICAL)=RTEMP
      ELSE
      WRITE (6,'(A)') ' CALIBR: ill-formed group pointer'
      END IF
C]
      ELSE IF (WD(1:4).EQ.'AUTO') THEN
C[
      IF ((ICAL.GT.0).AND.(JCAL.GT.0)) THEN
      DREFER(ICAL,JCAL)=ZERO
      DREFER(JCAL,ICAL)=ZERO
      VREFER(ICAL,JCAL)=ZERO
      VREFER(JCAL,ICAL)=ZERO
      ELSE
      WRITE (6,'(A)') ' CALIBR: ill-formed group pointer'
      END IF
C]
      ELSE IF (WD(1:4).EQ.'ERRS') THEN
      ACTION='ERRS'
      ELSE
      CALL CHKEND('CALIBR>',DONE)
      END IF
C]
      END DO
      DONE= .FALSE.
C[
      IF (CALCNT.EQ.0) THEN
      WRITE(6,'(A)')
     &  ' CALIBR: no calibration groups defined. uniform calibration'
      CALL FILL4(CALPTR,NATOM,1)
      END IF
C]
C check first if any coordinates are in the database
      IF (HPNMAT.EQ.0) THEN
C get distances for current coordinates
      CALL ALLDSP (HEAP(HPNHGL),
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      HPNMAT=ALLHP(IREAL8(MATDIM))
      HPNMA2=ALLHP(INTEG4(MATDIM))
      CALL FILLR8(HEAP(HPNMAT),MATDIM,ZERO)
      CALL FILL4(HEAP(HPNMA2),MATDIM,1)
      CALL FILLR8(HEAP(HPNRAV),NOEMAX,ZERO)
      CALL AMBACC(HEAP(HPNHGL),HEAP(HPNMAT),'ACCU',
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
      AMBSTR=1
      END IF
C
C now set errors or calibrate
C[
      IF (ACTION.EQ.'ERRS') THEN
      CALL NOERRS(HEAP(HPNHGL),HEAP(HPNMAT),CALSHL,CALSHU,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
C
      ELSE IF (CALMOD.EQ.'VOLU') THEN
C
C do volume dependent uniform calibration
      CALL CALLIN(CALMOD,HEAP(HPNHGL),HEAP(HPNCAL),
     &            HEAP(HPNMAT),DIAGON,CALSHL,CALSHU,
     &            HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &            HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &            HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &            HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &            HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))

C
      ELSE
C
C set up calibraion matrix
C now do following: go through DREFER and VREFER table.
C if DREFER=0 get DREFER from structures (hpnmat)
C if VREFER=0 get VREFER from input average (hpnvol,hpnmat)
C if auto VREFER and DREFER are 0
C
C[
      CALL CALFIL(CALMOD,HEAP(HPNHGL),HEAP(HPNCAL),
     &           HEAP(HPNMAT),DIAGON,CALSHL,CALSHU,
     &           HEAP(HPNORR),HEAP(HPNIPR),HEAP(HPNILS),HEAP(HPNJPR),
     &           HEAP(HPNJLS),HEAP(HPNCND),HEAP(HPNRAV),HEAP(HPNRRV),
     &           HEAP(HPNDIS),HEAP(HPNLOW),HEAP(HPNHIG),HEAP(HPNVOL),
     &           HEAP(HPNPP1),HEAP(HPNPP2),HEAP(HPNWGH),HEAP(HPNVIO),
     &           HEAP(HPNNSP),HEAP(HPNCV), HEAP(HPNPID))
C
      WRITE (6,'(A,I5)') ' CAILBR: Number of calibration groups: ',
     & CALCNT
      WRITE (6,'(A)') ' CAILBR: reference distances '
      DO I=1,CALCNT
      WRITE (6,'(10E12.3)')
     &   (MAX(RSMALL,DREFER(I,J))**(-ONE/CALEXP), J=1,CALCNT)
      END DO
      WRITE (6,'(A)') ' CAILBR: reference volumes '
      DO I=1,CALCNT
      WRITE (6,'(10E12.3)') (VREFER(I,J), J=1,CALCNT)
      END DO
      WRITE (6,'(A)') ' CAILBR: calibration factors c=d^-6/v '
      DO I=1,CALCNT
      WRITE (6,'(10E12.3)')
     &   (DREFER(I,J)/MAX(RSMALL,VREFER(I,J)), J=1,CALCNT)
      END DO
      WRITE (6,'(A)') ' CAILBR: isotropic calibration factors c=d^-6/v '
      DO I=1,CALCNT
      WRITE (6,'(10E12.3)') (CALFAC(I,J), J=1,CALCNT)
      END DO
C
      END IF
C]
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE CALFIL(MODE,NOEHGL,CALPTR,NOEMAT,DIAGON,
     &           CALSHL,CALSHU,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
C
      CHARACTER*4 MODE
      INTEGER NOEHGL(*), CALPTR(*)
      DOUBLE PRECISION NOEMAT(*),DIAGON(*),CALSHL,CALSHU
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SQERR, SUMDIS
      INTEGER I, J, N, L, NSTART, NSTOP, II, JJ
      INTEGER CALNUM
      INTEGER IAT, JAT, NMAT, IGRP, JGRP, NORR
      DOUBLE PRECISION RN, RNT, NISO, WEIGHT, SCALE,RN6
      DOUBLE PRECISION CALISO(MAXCAL)
      LOGICAL QFILL, QDIST
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
C
C check which elements of drefer and vrefer are not defined.
C initialize dcount and vcount to -1 if defined, 0 otherwise.
      QFILL=.FALSE.
      DO I=1,CALCNT
      DO J=1,CALCNT
C[
      IF (DREFER(I,J).LT.RSMALL) THEN
      DCOUNT(I,J)=0
      WRITE(6,'(A,I3,I3,A)')
     &  ' CALIBR: reference distance for ', I,J, ' set to average '
      QFILL=.TRUE.
      ELSE
      DCOUNT(I,J)=-1
      WRITE(6,'(A,I3,I3,A,E12.3)')
     &  ' CALIBR: reference distance for ', I,J, ' set to ', DREFER(I,J)
      END IF
C]
C[
      IF (VREFER(I,J).LT.RSMALL) THEN
      VCOUNT(I,J)=0
      WRITE(6,'(A,I3,I3,A)')
     &  ' CALIBR: reference volume for ', I,J, ' set to average '
      QFILL=.TRUE.
      ELSE
      VCOUNT(I,J)=-1
      WRITE(6,'(A,I3,I3,A,E12.3)')
     &  ' CALIBR: reference volume for ', I,J, ' set to ', VREFER(I,J)
      END IF
C]
      END DO
      END DO
C
      NMAT=0
      CALNUM=0
      DO N=NSTART,NSTOP
C loop over all pairs of atoms belonging to restraint N
C get r-6 summed distance
      L=0
      RN=ZERO
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
C noemat contains r-6 distances
      RN=RN+NOEMAT(NMAT+L)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
C
      RNT=RN
      RN6=RNT**(-SIXTH)
      NOERAV(N)=RN6
      QDIST=((RN6.LT.CALSHU).AND.(RN6.GE.CALSHL))
C[
      IF (QDIST) THEN
C[
      IF ((NOECV(N).NE.0)) THEN
      CALNUM=CALNUM+1
C
C
C calculate a scale factor for the volumes from the size of the diagonal peaks.
C we REQUIRE that the second selection contains the heteronuclear bound
C hydrogens.
      SCALE=ZERO
      L=0
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
      SCALE=SCALE+NOEMAT(NMAT+L)/RNT*DIAGON(JAT)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
C[]
      SCALE=ONE
C
C
C now calculate average distances and volumes for each
C element in vrefer and drefer that is not defined.
C the volumes of an ambiguous peak are scaled by the
C distances. r-6 sum distance is temporarily stored in
C noerrv.
      L=0
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
      RN=(NOEMAT(NMAT+L))
      IGRP=CALPTR(IAT)
      JGRP=CALPTR(JAT)
C
C test if the distance is smaller than the overall cutoff
C[
      IF ((RN.GT.CALDIS**(-SIX))
     &     .AND.((NOEVOL(N).GT.RSMALL)
     &        .OR.(VCOUNT(IGRP,JGRP).LT.0))) THEN
C[
      IF (DCOUNT(IGRP,JGRP).GE.0) THEN
      DREFER(IGRP,JGRP)=DREFER(IGRP,JGRP)+RN
      DCOUNT(IGRP,JGRP)=DCOUNT(IGRP,JGRP)+RN/RNT
      DREFER(JGRP,IGRP)=DREFER(IGRP,JGRP)
      DCOUNT(JGRP,IGRP)=DCOUNT(IGRP,JGRP)
      END IF
C]
      END IF
C]
C
C volumes are accumulated if either the distance criterion is
C satisfied, or if the reference distance is set explicitely.
C The reference volume is scaled by the ratio of the -sixth
C powers of the r6summed distance for an ambiguous peak, and
C the individual distance.
C[
      IF (VCOUNT(IGRP,JGRP).GE.0) THEN
C[
      IF ((NOEVOL(N).GT.RSMALL)
     &   .AND.((RN.GT.CALDIS**(-SIX))
     &        .OR.(DCOUNT(IGRP,JGRP).LT.0))) THEN
      VREFER(IGRP,JGRP)=VREFER(IGRP,JGRP)
     &  +RN/RNT*NOEVOL(N)/SCALE
      VCOUNT(IGRP,JGRP)=VCOUNT(IGRP,JGRP)+RN/RNT
      VREFER(JGRP,IGRP)=VREFER(IGRP,JGRP)
      VCOUNT(JGRP,IGRP)=VCOUNT(IGRP,JGRP)
      END IF
C]
      END IF
C]
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
C
      END IF
C]
      END IF
C]
      NMAT=NMAT+L
      END DO
      WRITE(6,'(A,I6,A)')
     &  'NOECAL: ',CALNUM, ' distances in calibration shell'
C}
      DO I=1,CALCNT
      DO J=1,CALCNT
      DREFER(I,J)=(DREFER(I,J)/(MAX(ONE,DCOUNT(I,J))))
      VREFER(I,J)=(VREFER(I,J)/(MAX(ONE,VCOUNT(I,J))))
      WRITE (6,'(2F15.5)') DCOUNT(I,J), VCOUNT(I,J)
C[
      IF (VREFER(I,J).GT.RSMALL) THEN
      CALFAC(I,J)=DREFER(I,J)/VREFER(I,J)
      ELSE
      CALFAC(I,J)=ZERO
      END IF
C]
      END DO
      END DO
C
C now make the calibration factors "isotropic", i.e.
C C_ab^2 = C_aa.C_bb
C for that, we search for fixed values in each row of the matrix
C i.e. vcount and dcount < 0.
C if we find some, we average over the fixed ones only.
C if we find none, we average over the whole row.
C[
      IF (MODE.EQ.'ISOT') THEN
      DO I=1,CALCNT
      CALISO(I)=ZERO
      NISO=ZERO
      DO J=1,CALCNT
C[
      IF ((VCOUNT(I,J).LT.0).AND.(DCOUNT(I,J).LT.0)) THEN
C[
      IF ((VCOUNT(J,J).LT.0).AND.(DCOUNT(J,J).LT.0)) THEN
      CALISO(I)=CALISO(I)+(CALFAC(I,J)**2/CALFAC(J,J))
      NISO=NISO+ONE
      END IF
C]
      END IF
C]
      END DO
C done the fixed calibration factors. if there are none...
C TODO: there is a bug somewhere here!!
C[
      IF (NISO.LE.ZERO) THEN
      DO J=1,CALCNT
C error propagation is treated very casually here...
      WEIGHT
     &  =MIN(MAX(DCOUNT(I,J),VCOUNT(I,J)),MAX(DCOUNT(J,J),VCOUNT(J,J)))
C[
      IF (CALFAC(J,J).GT.RSMALL) THEN
      CALISO(I)=
     &   CALISO(I)+WEIGHT*(CALFAC(I,J)**2/CALFAC(J,J))
      ELSE
      WEIGHT=ZERO
      END IF
C]
      NISO=NISO+WEIGHT
      END DO
      END IF
C]
C[
      IF (NISO.GT.RSMALL) THEN
      CALISO(I)=CALISO(I)/NISO
      ELSE
      CALISO(I)=ZERO
      END IF
C]
      END DO
C
C now we have to re-define the calibration factors.
      DO I=1,CALCNT
      DO J=1,CALCNT
      CALFAC(I,J)=SQRT(CALISO(I)*CALISO(J))
      END DO
      END DO
C
      END IF
C]
C start over again and apply calibration factors and calculate mean error
      CALNUM=0
      SQERR=ZERO
      SUMDIS=ZERO
      NMAT=0
      DO N=NSTART,NSTOP
C
C calculate a scale factor for the volumes from the size of the diagonal peaks.
C we REQUIRE that the second selection contains the heteronuclear bound
C hydrogens.
      SCALE=ZERO
      L=0
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
      SCALE=SCALE+NOEMAT(NMAT+L)/RNT*DIAGON(JAT)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
C[]
      SCALE=ONE
C
C apply calibration factor to distances in distance window
      RN6=NOERAV(N)
      RNT=RN6**(-SIX)
      QDIST=((RN6.LT.CALSHU).AND.(RN6.GE.CALSHL))
C[
      IF (QDIST) THEN
      NOERRV(N)=ZERO
      L=0
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      IGRP=CALPTR(IAT)
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      JGRP=CALPTR(JAT)
      L=L+1
C noemat contains r-6 distances
      RN=NOEMAT(NMAT+L)
      IF (NOEVOL(N).GT.RSMALL) THEN
      NOERRV(N)=NOERRV(N)+(RN/RNT)
     &   *CALFAC(IGRP,JGRP)*NOEVOL(N)/SCALE
      END IF
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
      IF (NOEVOL(N).GT.RSMALL) THEN
      NOEDIS(N)=NOERRV(N)**(-ONE/CALEXP)
      END IF
      END IF
C]
      NMAT=NMAT+L
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE CHKSEL
     &  (ISLCT,NISLCT,NATOM,CALPTR,CALCNT,ICAL,MAXCAL,QMATCH)
C
C
C Author: Michael Nilges, EMBL, 1994.
C
      IMPLICIT NONE
C
C input/ output
C   ISLCT:           flags for selected atoms
C   NISLCT:          number of selected atoms
C   CALPTR:          pointers to correlation time group
C   CALCNT:          lookup table for correlation times
      INTEGER ISLCT(*),NISLCT,NATOM,CALPTR(*),CALCNT,ICAL,MAXCAL
      LOGICAL QMATCH
C
C local
C   ICNT,IATOM:      loop variables (group number,atom number)
      INTEGER ICNT,IATOM,I
      LOGICAL QOVERL,QIDENT,QDISJU
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
C check if the groups just defined are disjunct with or identical to
C any existing groups
      QOVERL=.FALSE.
      QMATCH=.FALSE.
      DO ICNT=1,CALCNT
C
      QIDENT=.FALSE.
      QDISJU=.FALSE.
      DO IATOM=1,NATOM
      QIDENT=(QIDENT.OR.
     &       ((ISLCT(IATOM).EQ.1).AND.(CALPTR(IATOM).EQ.ICNT)))
      QDISJU=(QDISJU.OR.
     &       ((ISLCT(IATOM).EQ.1).AND.(CALPTR(IATOM).NE.ICNT)))
      QMATCH=(QMATCH.OR.QIDENT)
      QOVERL=(QOVERL.OR.(QIDENT.AND.QDISJU))
C[
      IF (QIDENT) THEN
      ICAL=ICNT
      END IF
C]
      END DO
C
      END DO
C[
      IF (QOVERL) THEN
      CALL WRNDIE(-1,'CHKSEL',
     &  'groups have to be disjunct or identical')
      ELSE IF (.NOT.QMATCH) THEN
C[
      IF (CALCNT+1.GT.MAXCAL) THEN
      CALL WRNDIE(-1,'CHKSEL',
     &  'maximum numer of groups exceeded ==> recompile program.')
      ELSE
      CALCNT=CALCNT+1
      ICAL=CALCNT
C now fill the CALPTR array with pointers
      DO I=1,NATOM
C[
      IF (ISLCT(I).GT.0) THEN
      CALPTR(I)=ICAL
      END IF
C]
      END DO
C
      END IF
C]
      END IF
C]
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NOERRS(NOEHGL,NOEMAT,CALSHL,CALSHU,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
C
      INTEGER NOEHGL(*)
      DOUBLE PRECISION NOEMAT(*), CALSHL,CALSHU
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SQERR, SUMDIS, SQERRM
      INTEGER N, L, NSTART, NSTOP, II, JJ
      INTEGER CALNUM
      INTEGER IAT, JAT, NMAT, NORR
      DOUBLE PRECISION RN, RNT
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
      CALNUM=0
      SQERR=ZERO
      SUMDIS=ZERO
      DO N=NSTART,NSTOP
C loop over all pairs of atoms belonging to restraint N
C get r-6 summed distance
      L=0
      RN=ZERO
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
C noemat contains r-6 distances
      RN=RN+NOEMAT(NMAT+L)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
      RNT=RN
      NOERAV(N)=RNT**(-SIXTH)
C[
      IF((NOERAV(N).LT.CALSHU).AND.(NOERAV(N).GE.CALSHL)) THEN
C[
      IF (NOECV(N).NE.0) THEN
      SQERR=SQERR+((NOEDIS(N)-NOERAV(N))/NOEDIS(N))**2
      SUMDIS=SUMDIS+NOEDIS(N)**2
      CALNUM=CALNUM+1
      END IF
C]
      END IF
C]
      NMAT=NMAT+L
      END DO
C
      SQERR=SQRT(MAX(RSMALL,SQERR/SUMDIS))
      WRITE(6,'(A,e15.5)')
     &  'NOECAL: rms error: ', SQERR
      SQERRM=SQERR
C
C set error estimates to distance times mean error
      DO N=NSTART,NSTOP
C[
      IF((NOERAV(N).LT.CALSHU).AND.(NOERAV(N).GE.CALSHL)) THEN
      NOELOW(N)=NOEDIS(N)*NOEDIS(N)*SQERRM
      NOEHIG(N)=NOEDIS(N)*NOEDIS(N)*SQERR
      END IF
C]
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE CALLIN(MODE,
     &           NOEHGL,CALPTR,NOEMAT,DIAGON,CALSHL,CALSHU,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
C
      CHARACTER*4 MODE
      INTEGER NOEHGL(*), CALPTR(*)
      DOUBLE PRECISION DIAGON(*), NOEMAT(*), CALSHL,CALSHU
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SUMVOL, SUMDIS
      INTEGER N, L, NSTART, NSTOP, II, JJ
      INTEGER CALNUM
      INTEGER IAT, JAT, NMAT, NORR
      DOUBLE PRECISION RN, RNT, RN6
      DOUBLE PRECISION ACOEFF, BCOEFF, AERR, BERR, CHI2
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
C
      NMAT=0
      CALNUM=0
      SUMDIS=ZERO
      SUMVOL=ZERO
      DO N=NSTART,NSTOP
C loop over all pairs of atoms belonging to restraint N
C get r-6 summed distance
      L=0
      RN=ZERO
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
C noemat contains r-6 distances
      RN=RN+NOEMAT(NMAT+L)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
C
      RNT=RN
      RN6=RNT**(-SIXTH)
      NOERAV(N)=RN6
C[
      IF ((NOECV(N).NE.0)) THEN
      CALNUM=CALNUM+1
      NOERAV(N)=RNT/NOEVOL(N)
      SUMDIS=SUMDIS+NOERAV(N)
      SUMVOL=SUMVOL+NOEVOL(N)
      END IF
C]
      NMAT=NMAT+L
      END DO
C
      WRITE(6,'(A,I6,A)')
     &  ' CALIBR: ',CALNUM, ' distances in calibration shell'
      WRITE(6,'(A,E15.5)')
     &  ' CALIBR: average calibration factor: ', SUMDIS/CALNUM
      WRITE(6,'(A,E15.5)')
     &  ' CALIBR: average volume: ', SUMVOL/CALNUM
C
      CALL LINFIT(0,NOENUM,NOEVOL,NOERAV,NOECV,
     &     ACOEFF,BCOEFF,AERR,BERR,CHI2)
      WRITE (6,'(A,E15.5,A,E15.5,A)')
     & ' CALIBR: calibration factors: ', ACOEFF,
     & ' + ', BCOEFF, ' * Volume.'
      WRITE (6,'(A,E15.5,E15.5,A,E15.5)')
     & ' CALIBR: errors: ', AERR, BERR, ' Chi2 ', CHI2
C
C start over again to apply calibration factors
      DO N=NSTART,NSTOP
      NOEDIS(N)=(ACOEFF+BCOEFF*NOEVOL(N))*NOEVOL(N)
      NOEDIS(N)=NOEDIS(N)**(-SIXTH)
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NOERRL(NOEHGL,NOEMAT,CALSHL,CALSHU,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
C
      INTEGER NOEHGL(*)
      DOUBLE PRECISION NOEMAT(*), CALSHL,CALSHU
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C local
      DOUBLE PRECISION SQERR, SUMDIS
      INTEGER N, L, NSTART, NSTOP, II, JJ
      INTEGER CALNUM
      INTEGER IAT, JAT, NMAT, NORR
      DOUBLE PRECISION RN, RN6
      DOUBLE PRECISION ACOEFF, BCOEFF, AERR, BERR, CHI2
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
C
      NMAT=0
      CALNUM=0
      SQERR=ZERO
      SUMDIS=ZERO
      DO N=NSTART,NSTOP
C loop over all pairs of atoms belonging to restraint N
C get r-6 summed distance
      L=0
      RN=ZERO
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
C noemat contains r-6 distances
      RN=RN+NOEMAT(NMAT+L)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
      RN6=RN**(-SIXTH)
      NOERAV(N)=-(NOEDIS(N)-RN6)/NOEDIS(N)
      WRITE (6,'(3E15.5)') NOERAV(N),NOEDIS(N),RN6
      NMAT=NMAT+L
      END DO
C
      CALL LINFIT(-1,NOENUM,NOEDIS,NOERAV,NOECV,
     &     ACOEFF,BCOEFF,AERR,BERR,CHI2)
      WRITE (6,'(A,E15.5,A,E15.5,A)')
     & ' CALIBR: relative error: ', ACOEFF,
     & ' + ', BCOEFF, ' * Distance.'
      WRITE (6,'(A,E15.5,E15.5,A,E15.5)')
     & ' CALIBR: errors: ', AERR, BERR, ' Chi2: ', CHI2
      NOELOW(N)=NOEDIS(N)*(ACOEFF+BCOEFF*NOEDIS(N))
      CALL LINFIT(+1,NOENUM,NOEDIS,NOERAV,NOECV,
     &     ACOEFF,BCOEFF,AERR,BERR,CHI2)
      WRITE (6,'(A,E15.5,A,E15.5,A)')
     & ' CALIBR: relative error: ', ACOEFF,
     & ' + ', BCOEFF, ' * Distance.'
      WRITE (6,'(A,E15.5,E15.5,A,E15.5)')
     & ' CALIBR: errors: ', AERR, BERR, ' Chi2: ', CHI2
      NOEHIG(N)=NOEDIS(N)*(ACOEFF+BCOEFF*NOEDIS(N))
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE NOERRLold(NOEHGL,NOEMAT,CALSHL,CALSHU,
     &           NOEORR,NOEIPR,NOEILS,NOEJPR,NOEJLS,NOECND,
     &           NOERAV,NOERRV,NOEDIS,NOELOW,NOEHIG,NOEVOL,
     &           NOEPP1,NOEPP2,NOEWGH,NOEVIO,NOESTP,NOECV,NOEPID)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'numbers.inc'
C
C global NOE arrays on HEAP
C restraint and atom pointers
      INTEGER NOEORR(*),NOEIPR(*),NOEILS(*),NOEJPR(*),NOEJLS(*)
C classes
      INTEGER NOECND(*)
C averages
      DOUBLE PRECISION NOERAV(*),NOERRV(*)
C target distance and errors
      DOUBLE PRECISION NOEDIS(*),NOELOW(*),NOEHIG(*)
C volume, chemical shifts
      DOUBLE PRECISION NOEVOL(*),NOEPP1(*),NOEPP2(*)
C weights
      DOUBLE PRECISION NOEWGH(*)
C number of violations
      INTEGER NOEVIO(*)
C time average pointer, test set, restraint number
      INTEGER NOESTP(*), NOECV(*), NOEPID(*)
C
      INTEGER NOEHGL(*)
      DOUBLE PRECISION NOEMAT(*), CALSHL,CALSHU
C local
      DOUBLE PRECISION SQERR, SUMDIS
      INTEGER N, L, NSTART, NSTOP, II, JJ
      INTEGER CALNUM
      INTEGER IAT, JAT, NMAT, NORR
      DOUBLE PRECISION RN, RN6
      DOUBLE PRECISION ACOEFF, BCOEFF, AERR, BERR, CHI2
C functions
      LOGICAL FRSTEL
      EXTERNAL FRSTEL
C
C begin
C
      NSTART=1
      NSTOP=NOENUM
      NMAT=0
      CALNUM=0
      SQERR=ZERO
      SUMDIS=ZERO
      DO N=NSTART,NSTOP
C loop over all pairs of atoms belonging to restraint N
C get r-6 summed distance
      L=0
      RN=ZERO
      DO NORR=NOEORR(N)+1,NOEORR(N+1)
      DO II=NOEIPR(NORR)+1,NOEIPR(NORR+1)
      IAT=NOEILS(II)
C[
      IF (FRSTEL(IAT,NOEHGL)) THEN
      DO JJ=NOEJPR(NORR)+1,NOEJPR(NORR+1)
      JAT=NOEJLS(JJ)
C[
      IF (FRSTEL(JAT,NOEHGL)) THEN
      L=L+1
C noemat contains r-6 distances
      RN=RN+NOEMAT(NMAT+L)
      END IF
C]
      END DO
      END IF
C]
      END DO
      END DO
      RN6=RN**(-SIXTH)
      NOERAV(N)=ABS(NOEDIS(N)-RN6)/NOEDIS(N)
      NMAT=NMAT+L
      END DO
C
      CALL LINFIT(0,NOENUM,NOEDIS,NOERAV,NOECV,
     &     ACOEFF,BCOEFF,AERR,BERR,CHI2)
      WRITE (6,'(A,E15.5,A,E15.5,A)')
     & ' CALIBR: relative error: ', ACOEFF,
     & ' + ', BCOEFF, ' * Distance.'
      WRITE (6,'(A,E15.5,E15.5,A,E15.5)')
     & ' CALIBR: errors: ', AERR, BERR, ' Chi2: ', CHI2
C
C set error estimates to distance times mean error
      DO N=NSTART,NSTOP
      NOELOW(N)=NOEDIS(N)*(ACOEFF+BCOEFF*NOEDIS(N))
      NOEHIG(N)=NOELOW(N)
      END DO
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE LINFIT(S,N,ABSZ,ORDI,FLAGS,
     &                  ACOEFF,BCOEFF,AERR,BERR,CHI2)
C
      IMPLICIT NONE
C
      INCLUDE 'numbers.inc'
C
      INTEGER S,N,FLAGS(N)
      DOUBLE PRECISION ABSZ(N),ORDI(N),ACOEFF,BCOEFF,AERR,BERR,CHI2
      DOUBLE PRECISION SABSZ,SABSZW,SORDI,ST2,T,SWEIGHT
      INTEGER I
      SABSZ=ZERO
      SORDI=ZERO
      ST2=ZERO
      BCOEFF=ZERO
      SWEIGHT=ZERO
      DO I=1,N
      IF ((FLAGS(I).NE.0).AND.(S*ORDI(I).GE.ZERO)) THEN
      SWEIGHT=SWEIGHT+ONE
      SABSZ=SABSZ+ABSZ(I)
      SORDI=SORDI+ORDI(I)
      END IF
      END DO
      SABSZW=SABSZ/SWEIGHT
      DO I=1,N
      IF ((FLAGS(I).NE.0).AND.(S*ORDI(I).GE.ZERO)) THEN
      T=ABSZ(I)-SABSZW
      ST2=ST2+T*T
      BCOEFF=BCOEFF+T*ORDI(I)
      END IF
      END DO
      BCOEFF=BCOEFF/ST2
      ACOEFF=(SORDI-BCOEFF*SABSZ)/SWEIGHT
      AERR=SQRT(((ONE+SABSZ*SABSZ)/(SWEIGHT*ST2))/SWEIGHT)
      BERR=SQRT(ONE/ST2)
      CHI2=ZERO
      DO I=1,N
      IF ((FLAGS(I).NE.0).AND.(S*ORDI(I).GE.ZERO)) THEN
      CHI2=CHI2+(ORDI(I)-ACOEFF-BCOEFF*ABSZ(I))**2
      END IF
      END DO
      AERR=AERR*SQRT(CHI2/(SWEIGHT-TWO))
      BERR=BERR*SQRT(CHI2/(SWEIGHT-TWO))
      RETURN
      END
C
C=====================================================================
C
