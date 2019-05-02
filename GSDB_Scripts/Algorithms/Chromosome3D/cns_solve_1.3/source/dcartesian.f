      SUBROUTINE DCART
C
C Cartesian molecular dynamics - using the Verlet Leapfrog algorithm
C
C Paul Adams and Axel T. Brunger - June 1997
C
      IMPLICIT NONE
C
C IO
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
C
C local
      CHARACTER*10 FORM
      CHARACTER*(COMMAX) CFILE
      INTEGER I, NSTEP, NPRINT, CMPERIOD, NTRFRQ
      INTEGER CRDIND(MAXA), CRDNUM, NSAVC, IUNCRD, STLEN
      DOUBLE PRECISION TSTEP, TTEMP, OFFSET, SCALE
      LOGICAL TCOUPLING, VSCALING, CMREMOVE, QFORM, QFIRSC
C
C pointer
      INTEGER APOINT
C
C defaults
      NSTEP=100
      NPRINT=1
      CMREMOVE=.TRUE.
      CMPERIOD=0
      TSTEP=0.001D0/TIMFAC
      TTEMP=300.0D0
      TCOUPLING=.FALSE.
      VSCALING=.FALSE.
      QFORM=.FALSE.
      OFFSET=800.0D0
      SCALE=10000.0D0
      FORM='12Z6'
      CFILE=' '
      NSAVC=0
      QFIRSC=.TRUE.
      IUNCRD=0
C
C set default selection for trajectory to be all atoms
      CRDNUM=0
      DO I=1,NATOM
      CRDNUM=CRDNUM+1
      CRDIND(I)=1
      END DO
C
C parsing
      CALL PUSEND('Cartesian Dynamics>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('Cartesian Dynamics>')
      CALL MISCOM('Cartesian Dynamics>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-dynamics-cartesian')
C
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('NSTEp=',NSTEP)
      ELSE IF (WD(1:4).EQ.'TIME') THEN
      CALL NEXTF('TIMEstep=',TSTEP)
      TSTEP=TSTEP/TIMFAC
      ELSE IF (WD(1:4).EQ.'TEMP') THEN
      CALL NEXTF('TEMPerature=',TTEMP)
      ELSE IF (WD(1:4).EQ.'TCOU') THEN
      CALL NEXTLO('TCOUpling=',TCOUPLING)
      IF (TCOUPLING) THEN
      WRITE(6,'(A)')
     &   ' DCART: temperature coupling (TCOUpling) enabled'
      IF (VSCALING) THEN
      VSCALING=.FALSE.
      WRITE(6,'(A)')
     &   ' DCART: velocity rescaling (VSCAling) disabled'
      END IF
      END IF
      ELSE IF (WD(1:4).EQ.'VSCA') THEN
      CALL NEXTLO('VSCAling=',VSCALING)
      IF (VSCALING) THEN
      WRITE(6,'(A)')
     &   ' DCART: velocity rescaling (VSCAling) enabled'
      IF (TCOUPLING) THEN
      TCOUPLING=.FALSE.
      WRITE(6,'(A)')
     &   ' DCART: temperature coupling (TCOUpling) disabled'
      END IF
      END IF
      ELSE IF (WD(1:4).EQ.'NPRI') THEN
      CALL NEXTI('NPRInt=',NPRINT)
      ELSE IF (WD(1:4).EQ.'CMRE') THEN
      CALL NEXTLO('CMREmove=',CMREMOVE)
      ELSE IF (WD(1:4).EQ.'CMPE') THEN
      CALL NEXTI('CMPEriodic=',CMPERIOD)
      ELSE IF (WD(1:4).EQ.'NTRF') THEN
      CALL NEXTI('NTRFrq=',NTRFRQ)
      IF (NTRFRQ.LE.0) CMREMOVE=.FALSE.
      IF (NTRFRQ.GE.1) CMREMOVE=.TRUE.
      ELSE IF (WD(1:5).EQ.'NSAVC') THEN
      CALL NEXTI('NSAVC=',NSAVC)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTST('FORMat=',FORM)
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
      CALL NEXTF('SCALe=',SCALE)
      ELSE IF (WD(1:4).EQ.'OFFS') THEN
      CALL NEXTF('OFFSet=',OFFSET)
      ELSE IF (WD(1:4).EQ.'TRAJ') THEN
      CALL NEXTFI('TRAJectory-file=',CFILE)
      ELSE IF (WD(1:4).EQ.'CSEL') THEN
      CALL SELCTA(CRDIND,CRDNUM,X,Y,Z,.TRUE.)
      ELSE IF (WD(1:4).EQ.'ASCI') THEN
      CALL NEXTLO('ASCIi=',QFORM)
      ELSE IF (WD(1:1).EQ.'?') THEN
      WRITE(6,'(A)')
     &' -----------------------------------------------------------'
      WRITE(6,'(A,I10,A,F8.5,A,I10,A)')
     & ' | NSTEp= ',NSTEP,' TIMEstep= ',TSTEP*TIMFAC,
     & ' NPRInt= ',NPRINT,' |'
      IF (CMREMOVE) THEN
      WRITE(6,'(A,I10,A)')
     & ' | CMREmove=true  CMPEriodic= ',CMPERIOD,
     & '                   |'
      ELSE
      WRITE(6,'(A,I10,A)')
     & ' | CMREmove=false  CMPEriodic= ',CMPERIOD,
     & '                  |'
      END IF
      IF (TCOUPLING) THEN
      WRITE(6,'(A,F10.4,A,I10,A)')
     & ' | TCOUpling=true  TEMPerature= ',TTEMP,'                 |'
      END IF
      IF (VSCALING) THEN
      WRITE(6,'(A,F10.4,A,I10,A)')
     & ' | VSCAling=true  TEMPerature= ',TTEMP,'                  |'
      END IF
      IF (CFILE.NE.' ') THEN
      IF (QFORM) THEN
      WRITE(6,'(A,I10,A,A)')
     & ' | NSAVC=',NSAVC,' ASCIi=true',
     & '                             |'
      ELSE
      WRITE(6,'(A,I10,A,A)')
     & ' | NSAVC=',NSAVC,' ASCIi=false',
     & '                            |'
      END IF
      WRITE(6,'(3A,F10.4,A,F10.4,A)')
     & ' | FORMat= ',FORM,' OFFSet= ',OFFSET,' SCALE= ',SCALE,' |'
      END IF
      WRITE(6,'(A)')
     &' -----------------------------------------------------------'
      IF (CFILE.NE.' ') THEN
      WRITE(6,'(A)')
     & ' Trajectory will be written to file:'
      STLEN=COMMAX
      CALL TRIMM(CFILE,STLEN)
      WRITE(6,'(A,A)')'   ',CFILE(1:STLEN)
      END IF
C
      ELSE
      CALL CHKEND('Cartesian Dynamics>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C open trajectory files
      IF (CFILE.NE.' ') THEN
      IF (QFORM) THEN
      CALL ASSFIL(CFILE,IUNCRD,'WRITE','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(CFILE,IUNCRD,'WRITE','UNFORMATTED',ERROR)
      END IF
      END IF
C
C set up the trajectory index lists
      CALL MAKIND(CRDIND,NATOM,CRDNUM)
C
      APOINT=ALLHP(INTEG4(NATOM))
C
      CALL DCART2(NSTEP,TSTEP,TTEMP,TCOUPLING,VSCALING,NPRINT,CMREMOVE,
     &            CMPERIOD,CRDNUM,CRDIND,QFIRSC,NSAVC,IUNCRD,QFORM,
     &            FORM,SCALE,OFFSET,HEAP(APOINT))
C
      CALL FREHP(APOINT,INTEG4(NATOM))
C
C close trajectory file
      IF (IUNCRD.GT.0) CALL VCLOSE(IUNCRD,'KEEP',ERROR)
C
      RETURN
      END
C
C ========================================================================
C
      SUBROUTINE DCART2(NSTEP,TSTEP,TTEMP,TCOUPLING,VSCALING,NPRINT,
     &                  CMREMOVE,CMPERIOD,CRDNUM,CRDIND,QFIRSC,
     &                  NSAVC,IUNCRD,QFORM,FORM,SCALE,OFFSET,APOINT)
C
C this routine implements the Verlet Leapfrog algorithm for dynamics
C
C given:
C
C    x (coordinates at time t)
C     t
C
C    v       (velocities at time t-dt/2)
C     t-dt/2
C
C    F (forces at time t)
C     t
C
C then new velocities at t+dt/2 are:
C
C    v       =  v       + (F dt)/m
C     t+dt/2     t-dt/2     t
C
C if velocity rescaling is requested:
C
C    v'      = [ v       ] sqrt[T      /T      ]
C     t+dt/2      t+dt/2         target  t+dt/2
C
C if temperature coupling (method of Berendsen) is requested:
C
C    v'      = [ v       ] sqrt( 1 + dt f [ ( T      /T       ) - 1 ] )
C     t+dt/2      t+dt/2                 b     target  t+dt/2
C
C if no temperature control is requested:
C
C    v'      = v
C     t+dt/2    t+dt/2
C
C then new coordinates at t+dt are:
C
C    x     =  x  + [ v'      ] dt
C     t+dt     t      t+dt/2
C
C
C Paul Adams - June 1997
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
C
C IO
      CHARACTER*10 FORM
      INTEGER NSTEP, NPRINT, CMPERIOD, CRDNUM, NSAVC, IUNCRD
      DOUBLE PRECISION TSTEP, TTEMP, SCALE, OFFSET
      LOGICAL TCOUPLING, VSCALING, CMREMOVE, QFIRSC, QFORM
      INTEGER APOINT(*), CRDIND(*)
C
C local
      INTEGER I, J, FREE, CYCLE
      DOUBLE PRECISION FACTOR, CTEMP
      DOUBLE PRECISION XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
C
C parameter
      DOUBLE PRECISION ZERO, HALF, ONE, TWO
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO=2.0D0)
C
C begin
C
      J=1
      DO I=1,NATOM
         IF (IMOVE(I).EQ.0) THEN
            APOINT(J)=I
            J=J+1
         END IF
      END DO
      FREE=J-1
C
C check masses before anything else is done
      DO I=1,FREE
      IF (AMASS(APOINT(I)).LE.ZERO) THEN
         CALL WRNDIE(-5,'DCART2',
     &   'there are atomic masses <= 0')
      END IF
      IF (AMASS(APOINT(I)).GE.9999) THEN
         CALL WRNDIE(-5,'DCART2',
     &   'there are very large (>= 9999) or undefined atomic masses')
      END IF
      END DO
C
      IF (FREE.GT.0) THEN
C
C do centre of mass motion removal if CMREMOVE is true
      IF (CMREMOVE) THEN
      CALL STOPRT(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.)
      END IF
C
C get current energies
      CALL ENERGY
C
C calculate the temperature and kinetic energy
      CALL GETKIN(XV,YV,ZV)
C
C store current temperature in temporary variable
      CTEMP=RENR(SSTEMP)
C
C do velocity rescaling if required
      IF (VSCALING) THEN
      IF (CTEMP.LE.RSMALL) THEN
      FACTOR=ONE
      ELSE
      FACTOR=SQRT(TTEMP/CTEMP)
      END IF
      DO I=1,FREE
      XV(APOINT(I))=XV(APOINT(I))*FACTOR
      YV(APOINT(I))=YV(APOINT(I))*FACTOR
      ZV(APOINT(I))=ZV(APOINT(I))*FACTOR
      END DO
      END IF
C
C do temperature coupling if required
      IF (TCOUPLING) THEN
      DO I=1,FREE
      IF (CTEMP.LE.RSMALL) THEN
      FACTOR=ONE
      ELSE
      FACTOR=SQRT( ONE + ( TIMFAC*TSTEP*FBETA(APOINT(I)) *
     &           ( (TTEMP/CTEMP ) - ONE ) ) )
      END IF
      XV(APOINT(I))=XV(APOINT(I))*FACTOR
      YV(APOINT(I))=YV(APOINT(I))*FACTOR
      ZV(APOINT(I))=ZV(APOINT(I))*FACTOR
      END DO
      END IF
C
C calculate the new temperature and kinetic energy
      CALL GETKIN(XV,YV,ZV)
C
      WRITE(6,'(2A)')
     & ' -------------------------- Cartesian dynamics start ',
     & '---------------------------'
      CALL PRINTD(RENR)
C
C start integration loop
      DO CYCLE=1,NSTEP
C
C do centre of mass motion removal if CMPERIOD is greater than 0
      IF (CMPERIOD.GT.0) THEN
      IF (MOD(CYCLE,CMPERIOD).EQ.0) THEN
      CALL STOPRT(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.)
      END IF
      END IF
C
C calculate velocities at t+dt/2
      DO I=1,FREE
      FACTOR=TSTEP/AMASS(APOINT(I))
      XV(APOINT(I))=XV(APOINT(I)) -
     &              DX(APOINT(I)) * FACTOR
      YV(APOINT(I))=YV(APOINT(I)) -
     &              DY(APOINT(I)) * FACTOR
      ZV(APOINT(I))=ZV(APOINT(I)) -
     &              DZ(APOINT(I)) * FACTOR
      END DO
C
C calculate the temperature and kinetic energy from new velocities
      CALL GETKIN(XV,YV,ZV)
C
C store current temperature in temporary variable
      CTEMP=RENR(SSTEMP)
C
C do velocity rescaling if required
      IF (VSCALING) THEN
      IF (CTEMP.LE.RSMALL) THEN
      FACTOR=ONE
      ELSE
      FACTOR=SQRT(TTEMP/CTEMP)
      END IF
      DO I=1,FREE
      XV(APOINT(I))=XV(APOINT(I))*FACTOR
      YV(APOINT(I))=YV(APOINT(I))*FACTOR
      ZV(APOINT(I))=ZV(APOINT(I))*FACTOR
      END DO
      END IF
C
C do temperature coupling if required
      IF (TCOUPLING) THEN
      DO I=1,FREE
      IF (CTEMP.LE.RSMALL) THEN
      FACTOR=ONE
      ELSE
      FACTOR=SQRT( ONE + ( TIMFAC*TSTEP*FBETA(APOINT(I)) *
     &           ( (TTEMP/CTEMP ) - ONE ) ) )
      END IF
      XV(APOINT(I))=XV(APOINT(I))*FACTOR
      YV(APOINT(I))=YV(APOINT(I))*FACTOR
      ZV(APOINT(I))=ZV(APOINT(I))*FACTOR
      END DO
      END IF
C
C do the integration to get coordinates at t+dt
      DO I=1,FREE
C
      X(APOINT(I))=(  X(APOINT(I)) +
     &               XV(APOINT(I)) * TSTEP )
      Y(APOINT(I))=(  Y(APOINT(I)) +
     &               YV(APOINT(I)) * TSTEP )
      Z(APOINT(I))=(  Z(APOINT(I)) +
     &               ZV(APOINT(I)) * TSTEP )
C
      END DO
C
C get new energies and forces on atoms DX,DY,DZ
      CALL ENERGY
C
C calculate the temperature and kinetic energy - from current velocities
      CALL GETKIN(XV,YV,ZV)
C
C Coordinate write options
      IF (NSAVC.GT.0.AND.IUNCRD.GT.0) THEN
      IF (MOD(CYCLE,NSAVC).EQ.0) THEN
      CALL WRITTC('CORD',NATOM,X,Y,Z,CRDNUM,CRDIND,CYCLE,QFIRSC,
     &            TSTEP,NSAVC,IUNCRD,QFORM,FORM,SCALE,OFFSET)
      END IF
      END IF
C
      IF (NPRINT.GT.0) THEN
C
C always print information for last step
      IF (CYCLE.EQ.NSTEP) THEN
      WRITE(6,'(A,I9,A,F12.5,A)')
     & ' ----------------- final step=',NSTEP,' at ',
     & TIMFAC*TSTEP*NSTEP,' ps ---------------------'
      CALL PRINTD(RENR)
C
C print information if at correct number of cycles
      ELSEIF (MOD(CYCLE,NPRINT).EQ.0) THEN
      WRITE(6,'(A,I9,A,F12.5,A)')
     & ' -------------------- step=',CYCLE,' at ',
     & TIMFAC*TSTEP*CYCLE,' ps ------------------------'
      CALL PRINTD(RENR)
      END IF
      END IF
C
C end integration loop
      END DO
C
C print current centre of mass motion information
      CALL  CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.,
     &             XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
C
      END IF
C
      RETURN
      END
C ========================================================================
      SUBROUTINE GETKIN(XV,YV,ZV)
C
C calculate the kinetic energy and temperature given velocities
C masses and information about moving atoms are taken from common block data
C
C Paul Adams - June 1997
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
C
C IO
      DOUBLE PRECISION XV(*), YV(*), ZV(*)
C
C local
      INTEGER I, NDEGF
      DOUBLE COMPLEX DCVAL
C
C begin
      RENR(SSTOTK)=0.0D0
      NDEGF=0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      NDEGF=NDEGF+3
      RENR(SSTOTK)=RENR(SSTOTK)+AMASS(I)*(XV(I)**2+YV(I)**2+ZV(I)**2)
      END IF
      END DO
C
      RENR(SSTOTK)=RENR(SSTOTK)*0.5D0
      RENR(SSTOTE)=RENR(SSENER)+RENR(SSTOTK)
      RENR(SSTEMP)=2.0D0*RENR(SSTOTK)/(NDEGF*KBOLTZ)
C
      CALL DECLAR( ANER(SSTOTK), 'DP', ' ', DCVAL, RENR(SSTOTK) )
      CALL DECLAR( ANER(SSTOTE), 'DP', ' ', DCVAL, RENR(SSTOTE) )
      CALL DECLAR( ANER(SSTEMP), 'DP', ' ', DCVAL, RENR(SSTEMP) )
C
      RETURN
      END
C=================================================================
      SUBROUTINE STOPRT(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,QPR)
C
C removes CM translational and rotational velocity from the
C system by inverting the inertia tensor
C
      IMPLICIT NONE
      INTEGER NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XV(*), YV(*), ZV(*)
      DOUBLE PRECISION AMASS(*)
      INTEGER IMOVE(*)
      LOGICAL QPR
C local
      INTEGER I
      DOUBLE PRECISION TCM(3,3),LH(3),MH(3),D
      DOUBLE PRECISION XX, XY, XZ, YY, YZ, ZZ, XI, YI, ZI
      DOUBLE PRECISION XCM, YCM, ZCM, VXCM, VYCM, VZCM, AXCM, AYCM, AZCM
      DOUBLE PRECISION OXCM, OYCM, OZCM
C begin
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,
     1      AMASS,IMOVE,QPR,XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
      XX = 0.0D0
      XY = 0.0D0
      XZ = 0.0D0
      YY = 0.0D0
      YZ = 0.0D0
      ZZ = 0.0D0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      XI = X(I) - XCM
      YI = Y(I) - YCM
      ZI = Z(I) - ZCM
      XX = XX + XI*XI*AMASS(I)
      XY = XY + XI*YI*AMASS(I)
      XZ = XZ + XI*ZI*AMASS(I)
      YY = YY + YI*YI*AMASS(I)
      YZ = YZ + YI*ZI*AMASS(I)
      ZZ = ZZ + ZI*ZI*AMASS(I)
      END IF
      END DO
      TCM(1,1) = YY + ZZ
      TCM(2,1) = -XY
      TCM(3,1) = -XZ
      TCM(1,2) = -XY
      TCM(2,2) = XX + ZZ
      TCM(3,2) = -YZ
      TCM(1,3) = -XZ
      TCM(2,3) = -YZ
      TCM(3,3) = XX + YY
C
C invert the inertia tensor TCM
      CALL MINVV(TCM,3,D,LH,MH)
C
C get angular velocity OXCM, OYCM, OZCM
      OXCM = AXCM*TCM(1,1) + AYCM*TCM(1,2) +  AZCM*TCM(1,3)
      OYCM = AXCM*TCM(2,1) + AYCM*TCM(2,2) +  AZCM*TCM(2,3)
      OZCM = AXCM*TCM(3,1) + AYCM*TCM(3,2) +  AZCM*TCM(3,3)
C
C remove CM translational and rotational motion from velocities
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      XI = X(I) - XCM
      YI = Y(I) - YCM
      ZI = Z(I) - ZCM
      XV(I) = XV(I) - VXCM - OYCM*ZI + OZCM*YI
      YV(I) = YV(I) - VYCM - OZCM*XI + OXCM*ZI
      ZV(I) = ZV(I) - VZCM - OXCM*YI + OYCM*XI
      END IF
      END DO
C
      WRITE(6,1000)
1000  FORMAT(' ROTSTP: CM translation and rotation removed.')
      IF (QPR) THEN
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,
     1      AMASS,IMOVE,QPR,XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
      END IF
C
      RETURN
      END
C=================================================================
      SUBROUTINE MINVV (A,N,D,L,M)
C
CCCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, APR. 1979 CCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE MINVV (A,N,D,L,M)                                      C
C                                                                      C
COMMENT   MINVV INVERTS THE N*N MATRIX A, USING THE STANDARD GAUSS-     C
C     JORDAN METHOD. THE DETERMINANT IS ALSO CALCULATED. IT HAS BEEN   C
C     COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.               C
C                                                                      C
C     A(1..N,1..N) = MATRIX TO BE INVERTED; DELIVERED WITH THE INVERSE C
C     N = ORDER OF MATRIX A                                            C
C     D = DELIVERED WITH THE DETERMINANT OF A                          C
C     L,M(1..N) = DUMMY ARRAYS                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT NONE
C input/ouput
      DOUBLE PRECISION A(*)
      INTEGER N
      DOUBLE PRECISION D
      INTEGER L(*), M(*)
C local
      DOUBLE PRECISION BIGA, HOLD
      INTEGER NK, K, I, J, IJ, KJ, JQ, JR, KI, JI
      INTEGER KK, IZ, JP, JK, IK
C begin
C
C*****SEARCH FOR LARGEST ELEMENT
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
      IF ( ABS(BIGA)- ABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C*****INTERCHANGE ROWS
      J=L(K)
      IF (J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C
C*****INTERCHANGE COLUMNS
   35 I=M(K)
      IF (I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
C
C*****DIVIDE COLUMN BY MINUS PIVOT
   45 IF (BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF (I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C*****REDUCE MATRIX
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF (I-K) 60,65,60
   60 IF (J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C*****DIVIDE ROW BY PIVOT
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF (J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C*****PRODUCT OF PIVOTS
      D=D*BIGA
C
C*****REPLACE PIVOT BY RECIPROCAL
      A(KK)=1.0D0/BIGA
   80 CONTINUE
C
C*****FINAL ROW AND COLUMN INTERCHANGE
      K=N
  100 K=(K-1)
      IF (K) 150,150,105
  105 I=L(K)
      IF (I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF (J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GOTO 100
  150 RETURN
      END
C=================================================================
      SUBROUTINE GAUSS (AM,SD,V,IGAUSS)
C
CCCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, MAR. 1979 CCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
COMMENT   GAUSS WILL SUPPLY A NORMALLY DISTRIBUTED RANDOM NUMBER WITH  C
C     A GIVEN MEAN AND STANDARD DEVIATION. IT USES 12 UNIFORM RANDOM   C
C     NUMBERS TO COMPUTE NORMAL RANDOM NUMBERS BY THE CENTRAL LIMIT    C
C     THEOREM (IBM 360 SCI.SUBR.PACKAGE).                              C
C                                                                      C
C     AM = THE DESIRED MEAN OF THE NORMAL DISTRIBUTION                 C
C     SD = THE DESIRED STANDARD DEVIATION OF THE NORMAL DISTRIBUTION   C
C     V = VALUE OF THE COMPUTED NORMAL RANDOM VARIABLE                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT NONE
C input/output
      INCLUDE 'funct.inc'
      DOUBLE PRECISION AM, SD, V
      INTEGER IGAUSS
C local
      DOUBLE PRECISION Y, A
      INTEGER I
C begin
      IF (IGAUSS.LE.0) THEN
C
C random distribution of negative & positive values of SD in VEL
      V = SD
      CALL GGUBFS(Y)
      IF(Y.LE.0.5D0) V = -SD
      ELSE
      A=0.0D0
      DO I=1,12
      CALL GGUBFS(Y)
      A=A+Y
      END DO
      V=(A-6.0D0)*SD+AM
      END IF
      RETURN
      END
C
C=================================================================
      SUBROUTINE PRINTD(RENRL)
C
C prints kinetic energy, temperature and partial energies
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION RENRL(*)
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
C local
      CHARACTER*20 ST1, ST2, ST3
      INTEGER ST1LEN, ST2LEN, ST3LEN
C parameters
      DOUBLE PRECISION MAXNUM, MINNUM
      PARAMETER (MAXNUM=999999.D0)
      PARAMETER (MINNUM=-99999.D0)
C begin
      IF (RENRL(SSTOTE).GT.MAXNUM.OR.RENRL(SSTOTE).LT.MINNUM) THEN
      WRITE(ST1,'(E10.2)') RENRL(SSTOTE)
      ELSE
      WRITE(ST1,'(F10.3)') RENRL(SSTOTE)
      END IF
      ST1LEN=11
      CALL TRIML(ST1,ST1LEN)
      IF (RENRL(SSTOTK).GT.MAXNUM.OR.RENRL(SSTOTK).LT.MINNUM) THEN
      WRITE(ST2,'(E10.2)') RENRL(SSTOTK)
      ELSE
      WRITE(ST2,'(F10.3)') RENRL(SSTOTK)
      END IF
      ST2LEN=11
      CALL TRIML(ST2,ST2LEN)
      IF (RENRL(SSTEMP).GT.MAXNUM.OR.RENRL(SSTEMP).LT.MINNUM) THEN
      WRITE(ST3,'(E10.2)') RENRL(SSTEMP)
      ELSE
      WRITE(ST3,'(F10.3)') RENRL(SSTEMP)
      END IF
      ST3LEN=11
      CALL TRIML(ST3,ST3LEN)
      WRITE(6,'(7A)')
     & ' | E(kin)+E(total)=',ST1(1:11),'     E(kin)=',
     & ST2(1:11),'   temperature=',ST3(1:11),'|'
      CALL PRINTE(RENRL)
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      RETURN
      END
C=================================================================
      SUBROUTINE CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,QPR,
     &                  XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
C
C Determines the center of mass, center of mass motion, and angular
C momentum of all non-fixed atoms.
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      INTEGER   NATOM
      DOUBLE PRECISION X(*), Y(*), Z(*), XV(*), YV(*), ZV(*)
      DOUBLE PRECISION AMASS(*)
      INTEGER IMOVE(*)
      LOGICAL   QPR
      DOUBLE PRECISION XCM,YCM,ZCM,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
C local
      DOUBLE PRECISION EKCM, TMASS
      INTEGER I
C begin
      VXCM = 0.0D0
      VYCM = 0.0D0
      VZCM = 0.0D0
      AXCM = 0.0D0
      AYCM = 0.0D0
      AZCM = 0.0D0
      XCM = 0.0D0
      YCM = 0.0D0
      ZCM = 0.0D0
      TMASS = 0.0D0
      IF (NATOM.GT.0) THEN
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      TMASS = TMASS + AMASS(I)
      VXCM = VXCM + XV(I)*AMASS(I)
      VYCM = VYCM + YV(I)*AMASS(I)
      VZCM = VZCM + ZV(I)*AMASS(I)
      XCM = XCM + X(I)*AMASS(I)
      YCM = YCM + Y(I)*AMASS(I)
      ZCM = ZCM + Z(I)*AMASS(I)
      AXCM=AXCM+(Y(I)*ZV(I)-Z(I)*YV(I))*AMASS(I)
      AYCM=AYCM+(Z(I)*XV(I)-X(I)*ZV(I))*AMASS(I)
      AZCM=AZCM+(X(I)*YV(I)-Y(I)*XV(I))*AMASS(I)
      END IF
      END DO
      AXCM = AXCM - (YCM*VZCM - ZCM*VYCM)/TMASS
      AYCM = AYCM - (ZCM*VXCM - XCM*VZCM)/TMASS
      AZCM = AZCM - (XCM*VYCM - YCM*VXCM)/TMASS
      XCM = XCM/TMASS
      YCM = YCM/TMASS
      ZCM = ZCM/TMASS
      VXCM = VXCM/TMASS
      VYCM = VYCM/TMASS
      VZCM = VZCM/TMASS
      END IF
      EKCM = (VXCM**2+VYCM**2+VZCM**2)*TMASS*0.5D0
C
      IF (QPR) THEN
C
C note the time conversion for velocities and momenta !
      WRITE(6,1000) XCM,YCM,ZCM,VXCM/TIMFAC,VYCM/TIMFAC,VZCM/TIMFAC,
     2             AXCM/TIMFAC,AYCM/TIMFAC,AZCM/TIMFAC,EKCM
1000  FORMAT(' CENMAS: Information about center of free masses',/,
     2       '         position [A]          :',3F13.5,/,
     3       '         velocity [A/ps]       :',3F13.5,/,
     4       '         ang. mom. [amu A/ps]  :',3F13.5,/,
     5       '         kin. ener. [Kcal/mol] :',F13.5)
      END IF
      RETURN
      END
