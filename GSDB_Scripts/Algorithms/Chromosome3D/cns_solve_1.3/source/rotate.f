      SUBROUTINE ROTMAT(ROT,T1,T2,T3,Q,AXIS,MODE)
C
C Routine computes unitary rotation matrix ROT using Eulerian angles
C (MODE="EULE"), Lattman angles (MODE="LATT"), spherical polar angles
C (MODE="SPHE"), a rotation about the specified axis (MODE="AXIS"),
C or quaternions (MODE='QUAT').
C
C Input:
C    MODE specifies angle mode
C    T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
C    T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
C    T1,T2,T3 are psi (incl. vs. y), phi (azimuthal angle, that is,
C             the angle between the x-axis and the projection of the
C             axis into the x,z plane ), kappa for MODE="SPHE"
C    T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
C  Note: all rotations are counter-clockwise
C
C Output:
C    ROT(3,3) contains the rotation matrix.  Should be applied as
C    r'(i)=sum_j ROT(i,j)*r(j)
C
C    AXIS contains the normalized input AXIS vector and T1,T2,T3 will
C    contain the spherical polar angles for MODE="SPHE"
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ROT(3,3), T1, T2, T3, Q(0:3), AXIS(3)
      CHARACTER*4 MODE
C local
      DOUBLE PRECISION S1, S2, S3, C1, C2, C3, C3C, S1SQ, NN, TP, TM
      DOUBLE PRECISION TT1, TT3
C parameter
      DOUBLE PRECISION RAD, ZERO, ONE, TWO, SMALL
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (SMALL=1.0D-7)
C begin
C=================================================================
      IF (MODE.EQ.'EULE'.OR.MODE.EQ.'LATT') THEN
C
C Computes rotation matrix corresponding to the three
C Eulerian angles theta1=T1, theta2=T2, theta3=T3.  This uses
C the Rossmann and Blow convention, i.e., a theta1 rotation
C around the Z axis, followed by
C a theta2 rotation around the moved X axis, and followed by a
C theta3 rotation around the moved Z axis.  These angles are positive
C if they are anti-clockwise when looking along the relevant axis.
C
      IF (MODE.EQ.'LATT') THEN
C In "Lattman" mode, theta+=theta1+theta3=T1, theta2=T2,
C theta-=theta1-theta3=T3
      TP=T1
      TM=T3
      TT1=(TP+TM)/TWO
      TT3=(TP-TM)/TWO
      ELSE
      TT1=T1
      TT3=T3
      END IF
C
      S1=SIN(RAD*TT1)
      S2=SIN(RAD*T2)
      S3=SIN(RAD*TT3)
      C1=COS(RAD*TT1)
      C2=COS(RAD*T2)
      C3=COS(RAD*TT3)
      ROT(1,1)=-S1*C2*S3+C1*C3
      ROT(1,2)=C1*C2*S3+S1*C3
      ROT(1,3)=S2*S3
      ROT(2,1)=-S1*C2*C3-C1*S3
      ROT(2,2)=C1*C2*C3-S1*S3
      ROT(2,3)=S2*C3
      ROT(3,1)=S1*S2
      ROT(3,2)=-C1*S2
      ROT(3,3)=C2
C=================================================================
      ELSE IF (MODE.EQ.'SPHE'.OR.MODE.EQ.'AXIS') THEN
C
      IF (MODE.EQ.'AXIS') THEN
C
C In axis mode we obtain the psi, phi spherical polar angles from
C the AXIS vector.
C The rotation kappa corresonds to an anticlockwise rotation (t3) about
C the specified axis. The axis is desribed as a vector (AXIS).
      NN=SQRT(AXIS(1)**2+AXIS(2)**2+AXIS(3)**2)
      IF (NN.LT.RSMALL) THEN
      CALL WRNDIE(-5,'ROTMAT','rotation vector has zero length')
      T1=ZERO
      T2=ZERO
      T3=ZERO
      ELSE
      AXIS(1)=AXIS(1)/NN
      AXIS(2)=AXIS(2)/NN
      AXIS(3)=AXIS(3)/NN
      T1=ACOS(AXIS(2))/RAD
      IF (AXIS(1)**2+AXIS(3)**2.LT.SMALL) THEN
      T2=ZERO
      ELSE
      T2=ACOS(MAX(-ONE,MIN(ONE,AXIS(1)/SQRT(AXIS(1)**2+AXIS(3)**2))))
     &    /RAD
      IF (AXIS(3).GT.ZERO) T2=-T2
      END IF
      END IF
      END IF
C
C compute rotation matrix corresponding to the three
C spherical polar angles psi=t1, phi=t2 and kappa=t3.  The rotation is
C described by specification of the direction of an axis through
C phi (azimutal angle between the x axis and the projection of the
C rotation axis on the x-z plane) and psi (inclination versus y axis).
C The angle kappa specifies the rotation around the specified axis.
C The kappa angle is anti-clockwise when looking along the rotation axis.
C The phi angle is anti-clockwise when looking along y.
      S1=SIN(RAD*T1)
      S2=SIN(RAD*T2)
      S3=SIN(RAD*T3)
      C1=COS(RAD*T1)
      C2=COS(RAD*T2)
      C3=COS(RAD*T3)
      S1SQ=S1*S1
      C3C=ONE-C3
      ROT(1,1) = C3           + S1SQ*C2*C2*C3C
      ROT(1,2) = S1*C1*C2*C3C   - S1*S2*S3
      ROT(1,3) =-S1SQ*C2*S2*C3C - C1*S3
      ROT(2,1) = S1*C1*C2*C3C   + S1*S2*S3
      ROT(2,2) = C3           + C1*C1*C3C
      ROT(2,3) =-S1*C1*S2*C3C   + S1*C2*S3
      ROT(3,1) =-S1SQ*S2*C2*C3C + C1*S3
      ROT(3,2) =-S1*C1*S2*C3C   - S1*C2*S3
      ROT(3,3) = C3           + S1SQ*S2*S2*C3C
C================================================================
      ELSE IF (MODE.EQ.'QUAT') THEN
      ROT(1,1)=Q(0)*Q(0)+Q(1)*Q(1)-Q(2)*Q(2)-Q(3)*Q(3)
      ROT(2,1)=TWO*(Q(1)*Q(2)-Q(0)*Q(3))
      ROT(3,1)=TWO*(Q(1)*Q(3)+Q(0)*Q(2))
      ROT(1,2)=TWO*(Q(1)*Q(2)+Q(0)*Q(3))
      ROT(2,2)=Q(0)*Q(0)-Q(1)*Q(1)+Q(2)*Q(2)-Q(3)*Q(3)
      ROT(3,2)=TWO*(Q(2)*Q(3)-Q(0)*Q(1))
      ROT(1,3)=TWO*(Q(1)*Q(3)-Q(0)*Q(2))
      ROT(2,3)=TWO*(Q(2)*Q(3)+Q(0)*Q(1))
      ROT(3,3)=Q(0)*Q(0)-Q(1)*Q(1)-Q(2)*Q(2)+Q(3)*Q(3)
C=================================================================
      END IF
C
      RETURN
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MATROT(ROT,T1,T2,T3,Q,AXIS,MODE)
C
C Routine computes Eulerian angles (MODE="EULE"), Lattman angles
C (MODE="LATT"), spherical polar angles (MODE="SPHE"),  rotation
C axis and angle (MODE="AXIS"), or quaternions
C from unitary matrix ROT.  The angular
C definitions are described in routine ROTMAT.
C
C Restrictions:
C    Matrix ROT has to be unitary, i.e., det(ROT)=+1 otherwise a fatal
C    warning will be issued.
C
C Conventions:
C   In Eulerian angle mode T2 can be forced without restriction of generality
C   to be located in the interval 0<= T2 <= pi (this is a consequence of
C   the identity operation t1,t2,t3 -> pi+t1,-t2,pi+t3 in Eulerian angle space).
C   T1 is forced into 0 <= t1 < 2*pi and T3 is forced into
C   0<=T3< 2*pi.  If T2 is 0 or PI we set T3 to zero without restriction
C   of generality.
C
C   Lattman angles are forced into the intervals (0 <= theta- <= 2*pi,
C   0<= theta2 <= pi, 0 <= theta+ < 4*pi).
C
C   For spherical polar angles we force psi (=t1) into 0<= t1 < pi and
C   phi (=t2) into 0<= t2 < pi, kappa (=t3) into 0 <= t3 < 2*pi
C   without restriction of generality.  If kappa is equal to zero,
C   all angles will be zero.  If psi is equal to zero, phi will be set to zero.
C
C   Spherical polar angles are used to compute the rotation axis in AXIS mode.
C   Without restriction of generality, kappa is forced into 0 <= kappa <= pi.
C
C Input:
C    MODE specifies angle mode
C    ROT(3,3) contains the rotation matrix.  The matrix is defines
C    the rotation according to r'(i)=sum_j ROT(i,j)*r(j)
C
C Ouput:
C    T1,T2,T3 are theta1 (z), theta2 (x'), theta3 (z') for MODE="EULE"
C    T1,T2,T3 are theta+, theta2, theta- for MODE="LATT"
C    T1,T2,T3 are psi (incl. vs. y), phi (azimuthal), kappa for MODE="SPHE"
C             and MODE="AXIS"
C    T3, AXIS(3) are kappa and a 3-D vector specifying the axis for MODE="AXIS"
C Note: all rotations are counter-clockwise
C
C Author: Axel T. Brunger
C =======================
C Fixed by J.-S. Jiang, 4/27/95
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ROT(3,3), T1, T2, T3, Q(0:3), AXIS(3)
      CHARACTER*4 MODE
C local
      DOUBLE PRECISION S1, S2, S3, C1, C2, C3, DET, C12, TP,TM
      DOUBLE PRECISION ROT2(3,3), NORM
      LOGICAL COND
      INTEGER I, J
C parameter
      DOUBLE PRECISION RAD, ZERO, ONE, TWO, SMALL, R180, R360, SMALLR
      DOUBLE PRECISION QUART
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (SMALL=1.0D-2, R180=180.0D0, R360=360.0D0)
      PARAMETER (SMALLR=RSMALL*100.D0, QUART=0.25D0)
C
C test to make sure that the matrix is unitary (allow
C mirror inversions as well)
      DET=(ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3)
     &   +(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3)
     &   +(ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)
C
      IF (ABS(DET-ONE).GT.SMALL) THEN
      CALL WRNDIE(-5,'MATROT','matrix ROT is not unitary')
      ELSE
C==================================================================
      IF (MODE.EQ.'EULE'.OR.MODE.EQ.'LATT') THEN
C
C first, let's get t2:
      C2=MAX(-ONE,MIN(ONE,ROT(3,3)))
      T2=ACOS(C2)/RAD
C
C second, let's compute SIN(T2) (OE positive)
C      S2=SQRT(MAX(ZERO,ONE-ROT(3,3)**2))
      S2=SQRT((ROT(3,1)**2+ROT(3,2)**2+ROT(1,3)**2+ROT(2,3)**2)/TWO)
C
C if SIN(T2) not equal 0 then T1 and T3 are uniquely determined
      IF (ABS(ABS(C2)-ONE).GT.SMALLR) THEN
      S3=ROT(1,3)/S2
      C3=ROT(2,3)/S2
      C3=MAX(-ONE,MIN(ONE,C3))
      S3=MAX(-ONE,MIN(ONE,S3))
      T3=ATAN2(S3,C3)+TWO*TWO*PI
      T3=MOD(T3,TWO*PI)/RAD
C
      S1=ROT(3,1)/S2
      C1=-ROT(3,2)/S2
      C1=MAX(-ONE,MIN(ONE,C1))
      S1=MAX(-ONE,MIN(ONE,S1))
      T1=ATAN2(S1,C1)+TWO*TWO*PI
      T1=MOD(T1,TWO*PI)/RAD
C
      ELSE
C
C without restriction T3 can be set to zero
      T3=ZERO
      C1=ROT(1,1)
      S1=ROT(1,2)
      C1=MAX(-ONE,MIN(ONE,C1))
      S1=MAX(-ONE,MIN(ONE,S1))
      T1=ATAN2(S1,C1)+TWO*TWO*PI
      T1=MOD(T1,TWO*PI)/RAD
C
      END IF
C
      IF (MODE.EQ.'LATT') THEN
C
C in "Lattman" mode compute theta+, theta- from theta1, theta3
      TP=T1+T3
      TM=T1-T3
      T1=TP
      T3=TM
C
C we know that 0 <= theta1 < 2*pi and 0 <= theta3 < 2*pi.  We have to
C project this into the Lattman space asymmetric unit (0<= theta+ < 4*pi,
C 0 <= theta- < 2*pi).
      IF (T3.LT.ZERO) THEN
      IF (T1.LT.R360) THEN
      T1=T1+R360
      T3=T3+R360
      ELSE
      T1=T1-R360
      T3=T3+R360
      END IF
      END IF
C
      END IF
C==========================================================
      ELSE IF (MODE.EQ.'SPHE'.OR.MODE.EQ.'AXIS') THEN
C
C first lets get COS (kappa)
      C3=(ROT(1,1)+ROT(2,2)+ROT(3,3)-ONE)/TWO
      C3=MAX(-ONE,MIN(ONE,C3))
      S3=SQRT( (ROT(2,1)-ROT(1,2))**2
     &        +(ROT(3,1)-ROT(1,3))**2
     &        +(ROT(2,3)-ROT(3,2))**2 )/TWO
      S3=MAX(-ONE,MIN(ONE,S3))
      T3=ATAN2(S3,C3)+TWO*TWO*PI
      T3=MOD(T3,TWO*PI)
C
C special case COS(kappa)=+1
      IF (ABS(C3-ONE).LT.SMALLR) THEN
      T1=ZERO
      T2=ZERO
      T3=ZERO
      ELSE
C
C determine COS^2(psi)
      C12=MAX(ZERO,(ROT(2,2)-C3)/(ONE-C3))
      C1=SQRT(C12)
      S1=SQRT(MAX(ZERO,((ROT(1,1)-C3)+(ROT(3,3)-C3))/(ONE-C3)))
C
C special case COS(psi)=+-1
      IF (ABS(C12-ONE).LT.SMALLR) THEN
      T1=ZERO
      T2=ZERO
      IF (ROT(3,1).LT.ZERO) T3=TWO*PI-T3
      ELSE
C
C determine phi:
C special case COS(kappa)=-1
      IF (ABS(C3+ONE).LT.SMALLR) THEN
      T3=PI
      C2=SQRT(MAX(ZERO,(ROT(1,1)+ONE)/(TWO*(ONE-C12))))
      S2=SQRT(MAX(ZERO,(ROT(3,3)+ONE)/(TWO*(ONE-C12))))
      IF (ROT(1,3).GT.ZERO) C2=-C2
      T2=ATAN2(S2,C2)+TWO*TWO*PI
      T2=MOD(T2,TWO*PI)
      ELSE
C
C otherwise determine phi by computing TAN(phi)
      C2=SQRT(MAX(ZERO,(ROT(1,1)-C3)/(ONE-C3)))
      S2=SQRT(MAX(ZERO,(ROT(3,3)-C3)/(ONE-C3)))
      IF ((-ROT(1,3)-ROT(3,1))/(ONE-C3).LT.ZERO) C2=-C2
      T2=ATAN2(S2,C2)+TWO*TWO*PI
      T2=MOD(T2,TWO*PI)
      END IF
C
C now determine signs of C1 and kappa.  We have a special case
C for COS(phi)=+-1
      IF (ABS(ABS(COS(T2))-ONE).LT.SMALL) THEN
      IF ((ROT(2,3)-ROT(3,2))/COS(T2).LT.ZERO) T3=TWO*PI-T3
      IF ((ROT(1,2)+ROT(2,1))/COS(T2).LT.ZERO) C1=-C1
      T1=ATAN2(S1,C1)+TWO*TWO*PI
      T1=MOD(T1,TWO*PI)
      ELSE
C determine sign of COS(psi) and then determine psi
      C1=SIGN(ONE,-ROT(2,3)-ROT(3,2))*C1
      T1=ATAN2(S1,C1)+TWO*TWO*PI
      T1=MOD(T1,TWO*PI)
      IF (ROT(1,2)-ROT(2,1).GT.ZERO) T3=TWO*PI-T3
      END IF
C
      END IF
      END IF
C
C convert all angles into degrees.
      T1=T1/RAD
      T2=T2/RAD
      T3=T3/RAD
C
      IF (MODE.EQ.'AXIS') THEN
C
      AXIS(2)=COS(T1*RAD)
      AXIS(1)=COS(T2*RAD)*SQRT(MAX(ZERO,ONE-AXIS(2)**2))
      AXIS(3)=-SQRT(MAX(ZERO,ONE-AXIS(1)**2-AXIS(2)**2))
C
C without restriction of generality we force T3 in 0 <= T3 <= pi
      IF (T3.GT.R180) THEN
      T3=TWO*R180-T3
      AXIS(1)=-AXIS(1)
      AXIS(2)=-AXIS(2)
      AXIS(3)=-AXIS(3)
      END IF
      END IF
C
      ELSEIF (MODE.EQ.'QUAT') THEN
      Q(0)=QUART*(ONE+ROT(1,1)+ROT(2,2)+ROT(3,3))
      Q(1)=QUART*(ONE+ROT(1,1)-ROT(2,2)-ROT(3,3))
      Q(2)=QUART*(ONE-ROT(1,1)+ROT(2,2)-ROT(3,3))
      Q(3)=QUART*(ONE-ROT(1,1)-ROT(2,2)+ROT(3,3))
      IF (ABS(Q(0)).GT.R4SMAL) THEN
      Q(0)=SQRT(Q(0))
      Q(1)=QUART*(ROT(2,3)-ROT(3,2))/Q(0)
      Q(2)=QUART*(ROT(3,1)-ROT(1,3))/Q(0)
      Q(3)=QUART*(ROT(1,2)-ROT(2,1))/Q(0)
      ELSE IF (ABS(Q(1)).GT.R4SMAL) THEN
      Q(1)=SQRT(Q(1))
      Q(0)=QUART*(ROT(2,3)-ROT(3,2))/Q(1)
      Q(2)=QUART*(ROT(1,2)+ROT(2,1))/Q(1)
      Q(3)=QUART*(ROT(1,3)+ROT(3,1))/Q(1)
      ELSE IF (ABS(Q(2)).GT.R4SMAL) THEN
      Q(2)=SQRT(Q(2))
      Q(0)=QUART*(ROT(3,1)-ROT(1,3))/Q(2)
      Q(1)=QUART*(ROT(1,2)+ROT(2,1))/Q(2)
      Q(3)=QUART*(ROT(2,3)+ROT(3,2))/Q(2)
      ELSE IF (ABS(Q(3)).GT.R4SMAL) THEN
      Q(3)=SQRT(Q(3))
      Q(0)=QUART*(ROT(1,2)-ROT(2,1))/Q(3)
      Q(1)=QUART*(ROT(1,3)+ROT(3,1))/Q(3)
      Q(2)=QUART*(ROT(2,3)+ROT(3,2))/Q(3)
      ELSE
      CALL WRNDIE(-5,'MATROT','quaternions cannot be defined.')
      END IF
      NORM=SQRT(Q(0)**2+Q(1)**2+Q(2)**2+Q(3)**2)
      Q(0)=Q(0)/NORM
      Q(1)=Q(1)/NORM
      Q(2)=Q(2)/NORM
      Q(3)=Q(3)/NORM
      END IF
C
C now back-compute the matrix as an internal consistency check
      CALL ROTMAT(ROT2,T1,T2,T3,Q,AXIS,MODE)
C
      COND=
     &     ABS(ROT(1,1)-ROT2(1,1)).GT.SMALL
     & .OR.ABS(ROT(1,2)-ROT2(1,2)).GT.SMALL
     & .OR.ABS(ROT(1,3)-ROT2(1,3)).GT.SMALL
     & .OR.ABS(ROT(2,1)-ROT2(2,1)).GT.SMALL
     & .OR.ABS(ROT(2,2)-ROT2(2,2)).GT.SMALL
     & .OR.ABS(ROT(2,3)-ROT2(2,3)).GT.SMALL
     & .OR.ABS(ROT(3,1)-ROT2(3,1)).GT.SMALL
     & .OR.ABS(ROT(3,2)-ROT2(3,2)).GT.SMALL
     & .OR.ABS(ROT(3,3)-ROT2(3,3)).GT.SMALL
      IF (COND) THEN
      WRITE(6,'(A,3G12.4)')
     &   ' %ROTMAT-ERR: inconsistent T1,T2,T3=',T1, T2, T3
      WRITE(6,'(/A,3(/3F12.6))')
     &  ' ROT =',((ROT(I,J),J=1,3),I=1,3)
      WRITE(6,'(/A,3(/3F12.6))')
     &  ' ROT2 =',((ROT2(I,J),J=1,3),I=1,3)
      CALL WRNDIE(-5,'MATROT','Error in internal consistency check')
      END IF
C
      END IF
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE XRFRAC(XRCELL,XRTR,XRINTR,XRVOL)
C
C Computes orthogonal to fractional transformation matrix XRTR and
C its inverse.  Similar to routine ORTHO.FOR in PROTEIN program
C package (W. Steigemann, MPI Biochemie, FRG).
C
C New coordinates are obtained by r'(i)=sum_j matrix(i,j)*r(j)
C where matrix is XRTR and XRINTR for orthogonal to fractional and
C fractional to orthogonal, respectively.
C
C The convention to setup the matrices is as follows: x same direction
C as a, y is in (a,b) plane.
C
C Author: Axel T. Brunger
C =======================
C Modification: Jian-Sheng Jiang
C (1) declare the unit cell dimensions in $XRCELL
C (2) declare orthogonal <--> fractional matrixes $XRTR and $XRINTR
C (3) declare reciprocal unit cell dimensions a*, b*, c* in $ABCS,
C
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
      DOUBLE PRECISION XRCELL(9), XRTR(3,3), XRINTR(3,3), XRVOL
C local
      INTEGER I, J
      DOUBLE PRECISION CABG(3), SABG(3), CABGS(3), ABCS(3), SABGS1
      DOUBLE COMPLEX DBCOMP
C parameters
      DOUBLE PRECISION RAD, ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, RAD=PI/180.0D0)
C begin
C
      DO I=1,3
      DO J=1,3
      XRTR(I,J)=ZERO
      XRINTR(I,J)=ZERO
      END DO
      END DO
C
      DO I=1,3
      CABG(I)=COS(XRCELL(I+3)*RAD)
      SABG(I)=SIN(XRCELL(I+3)*RAD)
      END DO
      CABGS(1)=(CABG(2)*CABG(3)-CABG(1))/(SABG(2)*SABG(3))
      CABGS(2)=(CABG(3)*CABG(1)-CABG(2))/(SABG(3)*SABG(1))
      CABGS(3)=(CABG(1)*CABG(2)-CABG(3))/(SABG(1)*SABG(2))
      XRVOL=XRCELL(1)*XRCELL(2)*XRCELL(3)*
     &                   SQRT(ONE+TWO*CABG(1)*CABG(2)*CABG(3)
     &                  -CABG(1)**2-CABG(2)**2-CABG(3)**2)
      ABCS(1)=XRCELL(2)*XRCELL(3)*SABG(1)/XRVOL
      ABCS(2)=XRCELL(1)*XRCELL(3)*SABG(2)/XRVOL
      ABCS(3)=XRCELL(1)*XRCELL(2)*SABG(3)/XRVOL
      SABGS1=SQRT(ONE-CABGS(1)**2)
C
C cartesian to fractional conversion
      XRTR(1,1)=ONE/XRCELL(1)
      XRTR(1,2)=-CABG(3)/(SABG(3)*XRCELL(1))
      XRTR(1,3)=-(CABG(3)*SABG(2)*CABGS(1)+CABG(2)*SABG(3))/
     &           (SABG(2)*SABGS1*SABG(3)*XRCELL(1))
      XRTR(2,2)=ONE/(SABG(3)*XRCELL(2))
      XRTR(2,3)=CABGS(1)/(SABGS1*SABG(3)*XRCELL(2))
      XRTR(3,3)=ONE/(SABG(2)*SABGS1*XRCELL(3))
C
C fractional to cartesian
      XRINTR(1,1)= XRCELL(1)
      XRINTR(1,2)= CABG(3)*XRCELL(2)
      XRINTR(1,3)= CABG(2)*XRCELL(3)
      XRINTR(2,2)= SABG(3)*XRCELL(2)
      XRINTR(2,3)=-SABG(2)*CABGS(1)*XRCELL(3)
      XRINTR(3,3)=SABG(2)*SABGS1*XRCELL(3)
C
C compute norm of A*, B*, C*
      XRCELL(7)=SQRT(XRTR(1,1)**2 + XRTR(1,2)**2 +XRTR(1,3)**2)
      XRCELL(8)=SQRT(XRTR(2,1)**2 + XRTR(2,2)**2 +XRTR(2,3)**2)
      XRCELL(9)=SQRT(XRTR(3,1)**2 + XRTR(3,2)**2 +XRTR(3,3)**2)
C
      CALL DECLAR('VOLUME','DP',' ',DBCOMP,XRVOL)
C
C declare the unit cell dimensions and matrixes - JSJ
      CALL DECLAR( 'XRCELL_1', 'DP', ' ', DBCOMP, XRCELL(1) )
      CALL DECLAR( 'XRCELL_2', 'DP', ' ', DBCOMP, XRCELL(2) )
      CALL DECLAR( 'XRCELL_3', 'DP', ' ', DBCOMP, XRCELL(3) )
      CALL DECLAR( 'XRCELL_4', 'DP', ' ', DBCOMP, XRCELL(4) )
      CALL DECLAR( 'XRCELL_5', 'DP', ' ', DBCOMP, XRCELL(5) )
      CALL DECLAR( 'XRCELL_6', 'DP', ' ', DBCOMP, XRCELL(6) )
      CALL DECLAR( 'XRTR_1_1', 'DP', ' ', DBCOMP, XRTR(1,1) )
      CALL DECLAR( 'XRTR_1_2', 'DP', ' ', DBCOMP, XRTR(1,2) )
      CALL DECLAR( 'XRTR_1_3', 'DP', ' ', DBCOMP, XRTR(1,3) )
      CALL DECLAR( 'XRTR_2_1', 'DP', ' ', DBCOMP, XRTR(2,1) )
      CALL DECLAR( 'XRTR_2_2', 'DP', ' ', DBCOMP, XRTR(2,2) )
      CALL DECLAR( 'XRTR_2_3', 'DP', ' ', DBCOMP, XRTR(2,3) )
      CALL DECLAR( 'XRTR_3_1', 'DP', ' ', DBCOMP, XRTR(3,1) )
      CALL DECLAR( 'XRTR_3_2', 'DP', ' ', DBCOMP, XRTR(3,2) )
      CALL DECLAR( 'XRTR_3_3', 'DP', ' ', DBCOMP, XRTR(3,3) )
      CALL DECLAR( 'XRINTR_1_1', 'DP', ' ', DBCOMP, XRINTR(1,1) )
      CALL DECLAR( 'XRINTR_1_2', 'DP', ' ', DBCOMP, XRINTR(1,2) )
      CALL DECLAR( 'XRINTR_1_3', 'DP', ' ', DBCOMP, XRINTR(1,3) )
      CALL DECLAR( 'XRINTR_2_1', 'DP', ' ', DBCOMP, XRINTR(2,1) )
      CALL DECLAR( 'XRINTR_2_2', 'DP', ' ', DBCOMP, XRINTR(2,2) )
      CALL DECLAR( 'XRINTR_2_3', 'DP', ' ', DBCOMP, XRINTR(2,3) )
      CALL DECLAR( 'XRINTR_3_1', 'DP', ' ', DBCOMP, XRINTR(3,1) )
      CALL DECLAR( 'XRINTR_3_2', 'DP', ' ', DBCOMP, XRINTR(3,2) )
      CALL DECLAR( 'XRINTR_3_3', 'DP', ' ', DBCOMP, XRINTR(3,3) )
C declare reciprocal unit cell dimensions
      CALL DECLAR( 'ASTAR', 'DP', ' ', DBCOMP, ABCS(1) )
      CALL DECLAR( 'BSTAR', 'DP', ' ', DBCOMP, ABCS(2) )
      CALL DECLAR( 'CSTAR', 'DP', ' ', DBCOMP, ABCS(3) )
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MATPAR(MODE,ROT)
C
C Routine parses matrix inputs:
C
C    MATRix <vector> <vector> <vector>  ! direct input of the matrix
C    EULEr <vector>  ! Eulerian angles T1, T2, T3
C    LATTman <vector> ! Lattman's angles T+, T2, T-
C    SPHErical < vector>  ! spherical polar angles psi, phi, kappa
C    AXIS <vector> <real>  !  axis vector and angle kappa
C    QUATernion <real> <real> <real> <real>! q0, q1, q2, q3 quaternions
C
C Remarks:
C    Eulerian angle convention: theta1 (z), theta2 (x'), theta3 (z')
C    Lattman's angles: T+=T1+t3, T-=T1-T3
C    Spherical polars:  psi (incl. vs. y), phi (azimuthal angle,
C             that is, the angle between the x-axis and the
C             projection of the axis into the x,z plane ), kappa
C             (rotation around the axis)
C Note: all rotations are counter-clockwise
C
C Output:
C    ROT(3,3) contains the matrix.  Should be applied as
C    r'(i)=sum_j ROT(i,j)*r(j)
C
C Author: Axel T. Brunger
C =======================
C
C     IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      CHARACTER*4 MODE
      DOUBLE PRECISION ROT(3,3)
C local
      DOUBLE PRECISION EU(3), R(3), Q(0:3), NORM
      INTEGER I
C begin
C
C====================================================================
      IF (MODE(1:4).EQ.'MATR') THEN
C
      CALL NEXTWD('MATRix-first-row=')
      CALL SAVEWD
      CALL NEXTVF('MATRix-first-row=',EU)
      DO I=1,3
      ROT(1,I)=EU(I)
      END DO
      CALL NEXTVF('MATRix-second-row=',EU)
      DO I=1,3
      ROT(2,I)=EU(I)
      END DO
      CALL NEXTVF('MATRix-third-row=',EU)
      DO I=1,3
      ROT(3,I)=EU(I)
      END DO
C====================================================================
      ELSE IF (MODE(1:4).EQ.'EULE') THEN
      CALL NEXTVF('EULErian-angles=',EU)
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'EULE')
C====================================================================
      ELSE IF (WD(1:4).EQ.'LATT') THEN
      CALL NEXTVF('Lattman-angles=',EU)
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'LATT')
C====================================================================
      ELSE IF (WD(1:4).EQ.'SPHE') THEN
      CALL NEXTVF('SPHErical-polar-angles=',EU)
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'SPHE')
C====================================================================
      ELSE IF (WD(1:4).EQ.'AXIS') THEN
      CALL NEXTVF('AXIS-vector=',R)
      CALL NEXTF('AXIS-angle=',EU(3))
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'AXIS')
C====================================================================
      ELSE IF (WD(1:4).EQ.'QUAT') THEN
      CALL NEXTF('Quaternion q0=',Q(0))
      CALL NEXTF('Quaternion q1=',Q(1))
      CALL NEXTF('Quaternion q2=',Q(2))
      CALL NEXTF('Quaternion q3=',Q(3))
      NORM=SQRT(Q(0)**2+Q(1)**2+Q(2)**2+Q(3)**2)
      Q(0)=Q(0)/NORM
      Q(1)=Q(1)/NORM
      Q(2)=Q(2)/NORM
      Q(3)=Q(3)/NORM
      WRITE(6,'(A,4F10.4/)')
     & ' Normalized quaternions: ',Q(0), Q(1), Q(2), Q(3)
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'QUAT')
C====================================================================
      END IF
C
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MATPRI(ROT)
C
C Prints matrix stored in ROT(3,3).
C If it is a rotation matrix then
C Eulerian, and spherical polar angles are printed.
C
      IMPLICIT NONE
C Author: Axel T. Brunger
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ROT(3,3)
C local
      DOUBLE PRECISION EU(3), R(3), TEMP, Q(0:3)
      INTEGER I, J, K
      LOGICAL QROT
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C check that it is a rotation matrix:
      QROT=.TRUE.
      DO I=1,3
      DO J=1,3
      TEMP=ZERO
      DO K=1,3
      TEMP=TEMP+ROT(I,K)*ROT(J,K)
      END DO
      IF (I.NE.J.AND.ABS(TEMP).GT.R4SMAL) THEN
      QROT=.FALSE.
      END IF
      IF (I.EQ.J.AND.ABS(TEMP-ONE).GT.R4SMAL) THEN
      QROT=.FALSE.
      END IF
      END DO
      END DO
      IF (QROT) THEN
      WRITE(6,'(/A,3(/3F12.6))')
     &  ' Rotation matrix =',((ROT(I,J),J=1,3),I=1,3)
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,R,'EULE')
      WRITE(6,'(A,3F10.4)')
     & ' Corresp. Eulerian angles (theta1,theta2,theta3) ',
     &  EU(1),EU(2),EU(3)
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,R,'SPHE')
      WRITE(6,'(A,3F10.4)')
     & ' Corresp. spherical polar angles (psi,phi,kappa) ',
     &  EU(1),EU(2),EU(3)
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,R,'AXIS')
      WRITE(6,'(A,F10.4,A,3F10.4)')
     & ' Corresp. rotation angle ',EU(3),' about axis   ',
     &  R(1),R(2),R(3)
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,R,'QUAT')
      WRITE(6,'(A,4F10.4/)')
     & ' Corresp. quaternions ',Q(0), Q(1), Q(2), Q(3)
      ELSE
      WRITE(6,'(/A,3(/3F12.6))')
     &  ' Matrix =',((ROT(I,J),J=1,3),I=1,3)
      END IF
      RETURN
      END
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MATDCL(ROT)
C
C Declares symbols that describe the rotation matrix
C stored in ROT(3,3).
C If it is a rotation matrix then
C Eulerian, and spherical polar angles are printed.
C
C Author: Axel T. Brunger
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ROT(3,3)
C local
      DOUBLE PRECISION EU(3), TEMP, Q(0:3), AXIS(3)
      INTEGER I, J, K
      LOGICAL QROT
      DOUBLE COMPLEX DBCOMP
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C check that it is a rotation matrix:
      QROT=.TRUE.
      DO I=1,3
      DO J=1,3
      TEMP=ZERO
      DO K=1,3
      TEMP=TEMP+ROT(I,K)*ROT(J,K)
      END DO
      IF (I.NE.J.AND.ABS(TEMP).GT.R4SMAL) THEN
      QROT=.FALSE.
      END IF
      IF (I.EQ.J.AND.ABS(TEMP-ONE).GT.R4SMAL) THEN
      QROT=.FALSE.
      END IF
      END DO
      END DO
      IF (QROT) THEN
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,AXIS,'EULE')
      CALL DECLAR( 'THETA1', 'DP', ' ', DBCOMP, EU(1) )
      CALL DECLAR( 'THETA2', 'DP', ' ', DBCOMP, EU(2) )
      CALL DECLAR( 'THETA3', 'DP', ' ', DBCOMP, EU(3) )
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,AXIS,'AXIS')
      CALL DECLAR( 'AXIS_X', 'DP', ' ', DBCOMP, AXIS(1) )
      CALL DECLAR( 'AXIS_Y', 'DP', ' ', DBCOMP, AXIS(2) )
      CALL DECLAR( 'AXIS_Z', 'DP', ' ', DBCOMP, AXIS(3) )
      CALL DECLAR( 'AXIS_KAPPA', 'DP', ' ', DBCOMP, EU(3) )
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,AXIS,'SPHE')
      CALL DECLAR( 'PSI', 'DP', ' ', DBCOMP, EU(1) )
      CALL DECLAR( 'PHI', 'DP', ' ', DBCOMP, EU(2) )
      CALL DECLAR( 'KAPPA', 'DP', ' ', DBCOMP, EU(3) )
      CALL MATROT(ROT,EU(1),EU(2),EU(3),Q,AXIS,'QUAT')
      CALL DECLAR( 'Q0', 'DP', ' ', DBCOMP, Q(0) )
      CALL DECLAR( 'Q1', 'DP', ' ', DBCOMP, Q(1) )
      CALL DECLAR( 'Q2', 'DP', ' ', DBCOMP, Q(2) )
      CALL DECLAR( 'Q3', 'DP', ' ', DBCOMP, Q(3) )
      END IF
      CALL DECLAR( 'ROT_1_1', 'DP', ' ', DBCOMP, ROT(1,1) )
      CALL DECLAR( 'ROT_1_2', 'DP', ' ', DBCOMP, ROT(1,2) )
      CALL DECLAR( 'ROT_1_3', 'DP', ' ', DBCOMP, ROT(1,3) )
      CALL DECLAR( 'ROT_2_1', 'DP', ' ', DBCOMP, ROT(2,1) )
      CALL DECLAR( 'ROT_2_2', 'DP', ' ', DBCOMP, ROT(2,2) )
      CALL DECLAR( 'ROT_2_3', 'DP', ' ', DBCOMP, ROT(2,3) )
      CALL DECLAR( 'ROT_3_1', 'DP', ' ', DBCOMP, ROT(3,1) )
      CALL DECLAR( 'ROT_3_2', 'DP', ' ', DBCOMP, ROT(3,2) )
      CALL DECLAR( 'ROT_3_3', 'DP', ' ', DBCOMP, ROT(3,3) )
      RETURN
      END
