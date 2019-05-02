      SUBROUTINE TORMD3(NSTEP,LIST,POINTR,QFIRSC,IUNCRD,DT,
     @     CRDNUM,CRDIND,NSAVC,NPRINT,QFORM,FORM,SCALE,OFFSET,
     @     CMREMOVE,CMPERIOD,GRPSEQ,NIN,NOUT,SIZE,XCM,YCM,ZCM,
     @     COEFFI,COEFFO,MM,YY,Q0,Q1,Q2,Q3,XB,YB,ZB,FRICT,QQ,B1,B2,
     @     DD,KK,LL,YYDOT,Q0O,Q1O,Q2O,Q3O,YYO,RKQ,RKQDOT,CHNLEN,
     @     QDOT,QDDOT,QDOTO,QO,Q,HHX,HHY,HHZ,SSX,SSY,SSZ,
     @     NCHAINS,TCOUPL,VSCALE,XCMO,
     @     YCMO,ZCMO,RKCRD,RKVEL)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'dtorsion.inc'
C
C I/O
      INTEGER NSTEP,LIST(*),POINTR(*)
      INTEGER NCHAINS(MAXTREE),GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER NOUT(MAXTREE,MAXCHN,*),SIZE(*),ISTEP,NPRINT,IGP,III
      INTEGER POINTER,CRDNUM,CRDIND(*),J
      INTEGER NIN(MAXTREE,MAXCHN,*),CHNLEN(MAXTREE,*)
      INTEGER IUNCRD,NSAVC,IIIM,CMPERIOD
      INTEGER ICHAIN,NNOUT,I,ITREE
      LOGICAL CMREMOVE
C
      LOGICAL TCOUPL
      LOGICAL VSCALE
C
      LOGICAL QFORM,QFIRSC,ERROR
C
      DOUBLE PRECISION XCMO(*),YCMO(*),ZCMO(*)
      DOUBLE PRECISION XCM(*),YCM(*),ZCM(*)
      DOUBLE PRECISION COEFFI(MAXTREE,MAXCHN,*),COEFFO(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION YY(6,NGP),Q0(*),Q1(*),Q2(*),Q3(*),XB(*)
      DOUBLE PRECISION MM(6,6,NGP),DTT,DT,FRICT(*),QQ(6,NGP)
      DOUBLE PRECISION B1(6,6,NGP),B2(6,NGP),KK(6,6,NGP)
      DOUBLE PRECISION DD(6,NGP),LL(6,NGP),YYDOT(6,NGP)
      DOUBLE PRECISION Q0O(*),Q1O(*),Q2O(*),SNORM
      DOUBLE PRECISION YYO(6,NGP),QVEL(0:3),QACC(0:3),K(7)
      DOUBLE PRECISION QSQR,NORM
      DOUBLE PRECISION XREL,YREL,ZREL
      DOUBLE PRECISION SCALE,OFFSET,QDOT(*)
      DOUBLE PRECISION Q3O(*),YB(*),ZB(*),HHX(*),HHY(*),HHZ(*)
      DOUBLE PRECISION SSX(*),SSY(*),SSZ(*),HNORM,CMX,CMY,CMZ,AXCM
      DOUBLE PRECISION AYCM,AZCM,VXCM,VYCM,VZCM
      DOUBLE PRECISION QDDOT(*),QDOTO(*)
      DOUBLE PRECISION RKCRD(7,MAXTREE),RKVEL(7,MAXTREE)
      DOUBLE PRECISION RKQ(*),RKQDOT(*),QO(*),Q(*)
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3,E0,E1,E2,E3
      DOUBLE PRECISION TEMP4,TEMP5,TEMP6,Q02,Q12,Q22,Q32
      DOUBLE PRECISION A(3,3), AO(3,3), DA(3,3)
C
      DOUBLE PRECISION HALF, EIGHTH, TWO, QUART, ONE, SIXTH, ZERO
      PARAMETER (HALF=0.5D0,EIGHTH=0.125D0,TWO=2.0D0,QUART=0.25D0)
      PARAMETER (ONE=1.0D0,SIXTH=1.0D0/6.0D0,ZERO=0.0D0)
C
      CHARACTER*10 FORM
C
      DO I=1,NATOM
      FRICT(I) = ZERO
      END DO
C
      DO I=1,NGP
      Q(I) = ZERO
      END DO
C
      DTT=DT/TIMFAC
C
C stop COM motion if CMREMOVE is true
      IF (CMREMOVE) THEN
      CALL RTSTOP(LIST,POINTR,NCHAINS,CHNLEN,YY,
     @     GRPSEQ,XCM,YCM,ZCM,COEFFO)
      END IF
C
      DO ISTEP=1,NSTEP
C
      IF (WRNLEV.GE.10) THEN
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.,
     &     CMX,CMY,CMZ,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
      END IF
C
C stop COM motion periodically if CMPERIOD is greater than 0
      IF (CMPERIOD.GT.0) THEN
      IF (MOD(ISTEP,CMPERIOD).EQ.0) THEN
      CALL RTSTOP(LIST,POINTR,NCHAINS,CHNLEN,YY,
     @     GRPSEQ,XCM,YCM,ZCM,COEFFO)
      END IF
      END IF
C
      CALL ENERGY
C
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
      IF (ISTEP.EQ.1) THEN
        IF (VSCALE) THEN
          CALL TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
          CALL TDYNKIN(NDEGF,XV,YV,ZV)
         ENDIF
      WRITE(6,'(A,I8)')
     &  ' TORMD3: number of degrees of freedom=',NDEGF
      WRITE(6,'(A,A)')
     @        ' -------------------------- Initial Conditions ',
     @        '---------------------------------'
      CALL PRINTD(RENR)
      END IF
C
C Now do T-coupling if required
      IF (TCOUPL) THEN
      IF (RENR(SSTEMP).GT.RSMALL) THEN
      DO I=1,NATOM
      FRICT(I) = FBETA(I)*TIMFAC*AMASS(I)*(ONE-TBATH/RENR(SSTEMP))
      END DO
      END IF
      END IF
C
C
C Compute velocity transform matrices, forces, etc ...
      CALL SWEEPOUT1(XCM,YCM,ZCM,FRICT,YY,QQ,MM,
     @     B1,B2,DD,GRPSEQ,SIZE,LIST,POINTR,NIN,
     @     NOUT,COEFFI,COEFFO,NCHAINS,CHNLEN,QDOT)
C
C Do inward recursion to get the KK and LL matrices
C  - these will allow calculation of the base body acceleration
      CALL SWEEPIN(MM,KK,B1,B2,DD,QQ,LL,GRPSEQ,
     @     NCHAINS,CHNLEN)
C
C Get base acceleration, and map it out the tree, including
C  relative accelerations along the way
      CALL SWEEPOUT2(GRPSEQ,KK,MM,B1,B2,DD,
     @     QQ,LL,YYDOT,NCHAINS,CHNLEN,QDDOT)
C
C Now take trial steps and update coordinates, velocities.
C Also compute the intermediate Runge Kutta predictions
      DO ITREE=1,NTREE
C
C First do the base
      III=GRPSEQ(ITREE,1,1)
C
C Save the initial coordinates and velocities
      XCMO(ITREE)=XCM(III)
      YCMO(ITREE)=YCM(III)
      ZCMO(ITREE)=ZCM(III)
      Q0O(III)=Q0(III)
      Q1O(III)=Q1(III)
      Q2O(III)=Q2(III)
      Q3O(III)=Q3(III)
C
      DO J=1,6
      YYO(J,III)=YY(J,III)
      END DO
C
C Set the Runge Kutta step accumulator to zero
      DO J=1,7
      RKCRD(J,ITREE)=ZERO
      RKVEL(J,ITREE)=ZERO
      END DO
C
C Transform angular terms to quaternion variables
      CALL AV2Q(Q0(III),Q1(III),Q2(III),Q3(III),YY(1,III),QVEL)
C
      CALL AA2Q(Q0(III),Q1(III),Q2(III),Q3(III),
     &          YY(1,III),YYDOT(1,III),QVEL,QACC)
C
C COMPUTE R-K TRIAL 1
      K(1)=DTT*YYDOT(1,III)
      K(2)=DTT*YYDOT(2,III)
      K(3)=DTT*YYDOT(3,III)
      K(4)=DTT*QACC(0)
      K(5)=DTT*QACC(1)
      K(6)=DTT*QACC(2)
      K(7)=DTT*QACC(3)
C
C Update the Runge Kutta final predictor
      DO J=1,7
         RKCRD(J,ITREE)=RKCRD(J,ITREE)+K(J)
         RKVEL(J,ITREE)=RKVEL(J,ITREE)+K(J)
      END DO
C
C Step to a midpoint
      XCM(III)=XCMO(ITREE)+HALF*DTT*YYO(1,III)+
     @     EIGHTH*DTT*K(1)
      YCM(III)=YCMO(ITREE)+HALF*DTT*YYO(2,III)+
     @     EIGHTH*DTT*K(2)
      ZCM(III)=ZCMO(ITREE)+HALF*DTT*YYO(3,III)+
     @     EIGHTH*DTT*K(3)
      Q0(III)=Q0O(III)+HALF*DTT*QVEL(0)+EIGHTH*DTT*K(4)
      Q1(III)=Q1O(III)+HALF*DTT*QVEL(1)+EIGHTH*DTT*K(5)
      Q2(III)=Q2O(III)+HALF*DTT*QVEL(2)+EIGHTH*DTT*K(6)
      Q3(III)=Q3O(III)+HALF*DTT*QVEL(3)+EIGHTH*DTT*K(7)
C
C Normalize
      QSQR=Q0(III)*Q0(III)+Q1(III)*Q1(III)+Q2(III)*
     @     Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
      YY(1,III)=YYO(1,III)+HALF*K(1)
      YY(2,III)=YYO(2,III)+HALF*K(2)
      YY(3,III)=YYO(3,III)+HALF*K(3)
C
      QVEL(0)=QVEL(0)+HALF*K(4)
      QVEL(1)=QVEL(1)+HALF*K(5)
      QVEL(2)=QVEL(2)+HALF*K(6)
      QVEL(3)=QVEL(3)+HALF*K(7)
C
      YY(4,III)=TWO*(-QVEL(0)*Q1(III)+QVEL(1)*Q0(III)-QVEL(2)*
     @     Q3(III)+QVEL(3)*Q2(III))
      YY(5,III)=TWO*(-QVEL(0)*Q2(III)+QVEL(2)*Q0(III)+QVEL(1)*
     @     Q3(III)-QVEL(3)*Q1(III))
      YY(6,III)=TWO*(-QVEL(1)*Q2(III)+QVEL(2)*Q1(III)-QVEL(0)*
     @     Q3(III)+QVEL(3)*Q0(III))
C
C Update atom coordinates, velocities
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
C Get rotation matrix which describes the step
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL1(A,AO,DA)
C
C Go on to all other bodies
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
C
C Save the initial coordinates and velocities
      Q0O(III)=Q0(III)
      Q1O(III)=Q1(III)
      Q2O(III)=Q2(III)
      Q3O(III)=Q3(III)
C
      QO(III)=Q(III)
      QDOTO(III)=QDOT(III)
C
      DO J=1,6
      YYO(J,III)=YY(J,III)
      END DO
C
C Set the Runge Kutta step accumulator to zero
      RKQ(III)=ZERO
      RKQDOT(III)=ZERO
C
C Compute R-K trial 1
      K(1)=DTT*QDDOT(III)
C
C Update the Runge Kutta final predictor
      RKQ(III)=RKQ(III)+K(1)
      RKQDOT(III)=RKQDOT(III)+K(1)
C
C Step to a midpoint
      Q(III)=QO(III)+HALF*DTT*QDOTO(III)+EIGHTH*DTT*K(1)
      QDOT(III)=QDOTO(III)+HALF*K(1)
C
C Rotate body about the old bond direction
      HNORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      HNORM=SQRT(HNORM)
      E0=COS(HALF*(Q(III)-QO(III)))
      E1=(HHX(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E2=(HHY(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E3=(HHZ(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      QSQR=E0**2+E1**2+E2**2+E3**2
      NORM=SQRT(QSQR)
      E0=E0/NORM
      E1=E1/NORM
      E2=E2/NORM
      E3=E3/NORM
C
      CALL TQ2ROT(E0,E1,E2,E3,AO)
C
      TEMP1=SSX(III)
      TEMP2=SSY(III)
      TEMP3=SSZ(III)
      SNORM=TEMP1**2+TEMP2**2+TEMP3**2
      SNORM=SQRT(SNORM)
      TEMP4=AO(1,1)*TEMP1+AO(1,2)*TEMP2+AO(1,3)*TEMP3
      TEMP5=AO(2,1)*TEMP1+AO(2,2)*TEMP2+AO(2,3)*TEMP3
      TEMP6=AO(3,1)*TEMP1+AO(3,2)*TEMP2+AO(3,3)*TEMP3
C
C Now rotate HH, SS by DA, the inboard body rotation
      TEMP1=DA(1,1)*TEMP4+DA(1,2)*TEMP5+DA(1,3)*TEMP6
      TEMP2=DA(2,1)*TEMP4+DA(2,2)*TEMP5+DA(2,3)*TEMP6
      TEMP3=DA(3,1)*TEMP4+DA(3,2)*TEMP5+DA(3,3)*TEMP6
      NORM=TEMP1**2+TEMP2**2+TEMP3**2
      NORM=SQRT(NORM)
      TEMP1=TEMP1*(SNORM/NORM)
      TEMP2=TEMP2*(SNORM/NORM)
      TEMP3=TEMP3*(SNORM/NORM)
      TEMP4=DA(1,1)*HHX(III)+DA(1,2)*HHY(III)+DA(1,3)*HHZ(III)
      TEMP5=DA(2,1)*HHX(III)+DA(2,2)*HHY(III)+DA(2,3)*HHZ(III)
      TEMP6=DA(3,1)*HHX(III)+DA(3,2)*HHY(III)+DA(3,3)*HHZ(III)
      NORM=TEMP4**2+TEMP5**2+TEMP6**2
      NORM=SQRT(NORM)
      TEMP4=TEMP4*(HNORM/NORM)
      TEMP5=TEMP5*(HNORM/NORM)
      TEMP6=TEMP6*(HNORM/NORM)
C
C Get center of mass position
      NNOUT=NOUT(ITREE,ICHAIN,IGP-1)
      XCM(III)=X(NNOUT)+TEMP4-TEMP1
      YCM(III)=Y(NNOUT)+TEMP5-TEMP2
      ZCM(III)=Z(NNOUT)+TEMP6-TEMP3
C
C Compute new DA, where DA represents the overall rotation
C   applied to this body, ie DA = DA_{inboard} * AO_{relative}
      CALL ROTMUL2(DA,AO,A)
C
      DO J=1,3
         DO I=1,3
            DA(I,J)=A(I,J)
         END DO
      END DO
C
C Get new quaternions to describe the body frame orientation
C   The new matrix from body frame to lab frame is A = DA * AO^T_{lab->body}
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),AO)
C
      CALL ROTMUL3(DA,AO,A)
C
      Q02=QUART*(ONE+A(1,1)+A(2,2)+A(3,3))
      Q12=QUART*(ONE+A(1,1)-A(2,2)-A(3,3))
      Q22=QUART*(ONE-A(1,1)+A(2,2)-A(3,3))
      Q32=QUART*(ONE-A(1,1)-A(2,2)+A(3,3))
      IF (ABS(Q02).GT.R4SMAL) THEN
      Q0(III)=SQRT(Q02)
      Q1(III)=QUART*(A(2,3)-A(3,2))/Q0(III)
      Q2(III)=QUART*(A(3,1)-A(1,3))/Q0(III)
      Q3(III)=QUART*(A(1,2)-A(2,1))/Q0(III)
      ELSE IF (ABS(Q12).GT.R4SMAL) THEN
      Q1(III)=SQRT(Q12)
      Q0(III)=QUART*(A(2,3)-A(3,2))/Q1(III)
      Q2(III)=QUART*(A(1,2)+A(2,1))/Q1(III)
      Q3(III)=QUART*(A(1,3)+A(3,1))/Q1(III)
      ELSE IF (ABS(Q22).GT.R4SMAL) THEN
      Q2(III)=SQRT(Q22)
      Q0(III)=QUART*(A(3,1)-A(1,3))/Q2(III)
      Q1(III)=QUART*(A(1,2)+A(2,1))/Q2(III)
      Q3(III)=QUART*(A(2,3)+A(3,2))/Q2(III)
      ELSE IF (ABS(Q32).GT.R4SMAL) THEN
      Q3(III)=SQRT(Q32)
      Q0(III)=QUART*(A(1,2)-A(2,1))/Q3(III)
      Q1(III)=QUART*(A(1,3)+A(3,1))/Q3(III)
      Q2(III)=QUART*(A(2,3)+A(3,2))/Q3(III)
      ELSE
      WRITE(6,'(A)') 'ERROR: matrix is not a rotation matrix'
      END IF
      QSQR= Q0(III)*Q0(III)+Q1(III)*Q1(III)
     &     +Q2(III)*Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
C Update coordinates
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      END DO
C
C Get body velocities
      YY(1,III)=YY(1,IIIM)+(-ZCM(IIIM)+ZCM(III))*
     @     YY(5,IIIM)+(YCM(IIIM)-YCM(III))*YY(6,IIIM)+
     @     ((TEMP6*TEMP2-TEMP5*TEMP3)/HNORM)*QDOT(III)
      YY(2,III)=YY(2,IIIM)+(ZCM(IIIM)-ZCM(III))*
     @     YY(4,IIIM)+(-XCM(IIIM)+XCM(III))*YY(6,IIIM)+
     @     ((TEMP4*TEMP3-TEMP6*TEMP1)/HNORM)*QDOT(III)
      YY(3,III)=YY(3,IIIM)+(-YCM(IIIM)+YCM(III))*
     @     YY(4,IIIM)+(XCM(IIIM)-XCM(III))*YY(5,IIIM)+
     @     ((TEMP5*TEMP1-TEMP4*TEMP2)/HNORM)*QDOT(III)
      YY(4,III)=YY(4,IIIM)+QDOT(III)*TEMP4/HNORM
      YY(5,III)=YY(5,IIIM)+QDOT(III)*TEMP5/HNORM
      YY(6,III)=YY(6,IIIM)+QDOT(III)*TEMP6/HNORM
C
C Update atom velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      END DO
C
C Put in a special test in case of a single atom chain tip
      IF ((SIZE(III).EQ.0).AND.
     @     COEFFO(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)).EQ.ONE) THEN
      IIIM=III
      III=GRPSEQ(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)+1)
      POINTER=POINTR(III)+1
      X(LIST(POINTER))=XCM(IIIM)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(IIIM)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(IIIM)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(IIIM)
      YREL=Y(LIST(POINTER))-YCM(IIIM)
      ZREL=Z(LIST(POINTER))-ZCM(IIIM)
      XV(LIST(POINTER))=YY(1,IIIM)+YY(5,IIIM)*ZREL
     @     -YY(6,IIIM)*YREL
      YV(LIST(POINTER))=YY(2,IIIM)+YY(6,IIIM)*XREL
     @     -YY(4,IIIM)*ZREL
      ZV(LIST(POINTER))=YY(3,IIIM)+YY(4,IIIM)*YREL
     @     -YY(5,IIIM)*XREL
      END IF
C
C Recover the DA rotation matrix corresponding
C   to the base of the next chain
      IF (ICHAIN.NE.NCHAINS(ITREE)) THEN
      III=GRPSEQ(ITREE,ICHAIN+1,1)
C
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      CALL ROTMUL1(A,AO,DA)
C
      END IF
C
      END DO
C
      END DO
C
      CALL ENERGY
C
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
C
C Now do T-coupling if required
      IF (TCOUPL) THEN
      IF (RENR(SSTEMP).GT.RSMALL) THEN
      DO I=1,NATOM
      FRICT(I) = FBETA(I)*TIMFAC*AMASS(I)*(ONE-TBATH/RENR(SSTEMP))
      END DO
      END IF
      END IF
C
C      IF (VSCALE) THEN
C        CALL TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
C      END IF
C
C Compute velocity transform matrices, forces, etc ...
      CALL SWEEPOUT1(XCM,YCM,ZCM,FRICT,YY,QQ,MM,
     @     B1,B2,DD,GRPSEQ,SIZE,LIST,POINTR,NIN,
     @     NOUT,COEFFI,COEFFO,NCHAINS,CHNLEN,QDOT)
C
C Do inward recursion to get the KK and LL matrices
C  - these will allow calculation of the base body acceleration
      CALL SWEEPIN(MM,KK,B1,B2,DD,QQ,LL,GRPSEQ,
     @     NCHAINS,CHNLEN)
C
C Get base acceleration, and map it out the tree, including
C  relative accelerations along the way
      CALL SWEEPOUT2(GRPSEQ,KK,MM,B1,B2,DD,
     @     QQ,LL,YYDOT,NCHAINS,CHNLEN,QDDOT)
C
C Now take trial steps and update coordinates, velocities.
C Also compute the intermediate Runge Kutta predictions
      DO ITREE=1,NTREE
C
C First do the base
      III=GRPSEQ(ITREE,1,1)
C
C Transform angular terms to quaternion variables
      CALL AV2Q(Q0(III),Q1(III),Q2(III),Q3(III),YY(1,III),QVEL)
C
      CALL AA2Q(Q0(III),Q1(III),Q2(III),Q3(III),
     &          YY(1,III),YYDOT(1,III),QVEL,QACC)
C
C Compute R-K trial 2
      K(1)=DTT*YYDOT(1,III)
      K(2)=DTT*YYDOT(2,III)
      K(3)=DTT*YYDOT(3,III)
      K(4)=DTT*QACC(0)
      K(5)=DTT*QACC(1)
      K(6)=DTT*QACC(2)
      K(7)=DTT*QACC(3)
C
C Update the Runge Kutta final predictor
      DO J=1,7
         RKCRD(J,ITREE)=RKCRD(J,ITREE)+K(J)
         RKVEL(J,ITREE)=RKVEL(J,ITREE)+TWO*K(J)
      END DO
C
C Recompute the initial angular velocities -
C   we need them for the runge kutta integration
      CALL AV2Q(Q0O(III),Q1O(III),Q2O(III),Q3O(III),YYO(1,III),QVEL)
C
C Now we make the second step - only the velocities change this time
      YY(1,III)=YYO(1,III)+HALF*K(1)
      YY(2,III)=YYO(2,III)+HALF*K(2)
      YY(3,III)=YYO(3,III)+HALF*K(3)
C
      QVEL(0)=QVEL(0)+HALF*K(4)
      QVEL(1)=QVEL(1)+HALF*K(5)
      QVEL(2)=QVEL(2)+HALF*K(6)
      QVEL(3)=QVEL(3)+HALF*K(7)
C
      YY(4,III)=TWO*(-QVEL(0)*Q1(III)+QVEL(1)*Q0(III)-QVEL(2)*
     @     Q3(III)+QVEL(3)*Q2(III))
      YY(5,III)=TWO*(-QVEL(0)*Q2(III)+QVEL(2)*Q0(III)+QVEL(1)*
     @     Q3(III)-QVEL(3)*Q1(III))
      YY(6,III)=TWO*(-QVEL(1)*Q2(III)+QVEL(2)*Q1(III)-QVEL(0)*
     @     Q3(III)+QVEL(3)*Q0(III))
C
C Update atom velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
C Get rotation matrix which describes the step
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL1(A,AO,DA)
C
C Now go on to other bodies
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
C
C Compute R-K trial 2
      K(1)=DTT*QDDOT(III)
C
C Update the Runge Kutta final predictor
      RKQ(III)=RKQ(III)+K(1)
      RKQDOT(III)=RKQDOT(III)+TWO*K(1)
C
C Now we make the second step - only the velocities change this time
      QDOT(III)=QDOTO(III)+HALF*K(1)
C
C Recompute the new HH, SS vectors
      HNORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      HNORM=SQRT(HNORM)
      E0=COS(HALF*(Q(III)-QO(III)))
      E1=(HHX(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E2=(HHY(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E3=(HHZ(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      QSQR=E0**2+E1**2+E2**2+E3**2
      NORM=SQRT(QSQR)
      E0=E0/NORM
      E1=E1/NORM
      E2=E2/NORM
      E3=E3/NORM
C
      CALL TQ2ROT(E0,E1,E2,E3,AO)
C
      TEMP1=SSX(III)
      TEMP2=SSY(III)
      TEMP3=SSZ(III)
      SNORM=TEMP1**2+TEMP2**2+TEMP3**2
      SNORM=SQRT(SNORM)
      TEMP4=AO(1,1)*TEMP1+AO(1,2)*TEMP2+AO(1,3)*TEMP3
      TEMP5=AO(2,1)*TEMP1+AO(2,2)*TEMP2+AO(2,3)*TEMP3
      TEMP6=AO(3,1)*TEMP1+AO(3,2)*TEMP2+AO(3,3)*TEMP3
C
C Now rotate by DA, the inboard body rotation
      TEMP1=DA(1,1)*TEMP4+DA(1,2)*TEMP5+DA(1,3)*TEMP6
      TEMP2=DA(2,1)*TEMP4+DA(2,2)*TEMP5+DA(2,3)*TEMP6
      TEMP3=DA(3,1)*TEMP4+DA(3,2)*TEMP5+DA(3,3)*TEMP6
      NORM=TEMP1**2+TEMP2**2+TEMP3**2
      NORM=SQRT(NORM)
      TEMP1=TEMP1*(SNORM/NORM)
      TEMP2=TEMP2*(SNORM/NORM)
      TEMP3=TEMP3*(SNORM/NORM)
      TEMP4=DA(1,1)*HHX(III)+DA(1,2)*HHY(III)+DA(1,3)*HHZ(III)
      TEMP5=DA(2,1)*HHX(III)+DA(2,2)*HHY(III)+DA(2,3)*HHZ(III)
      TEMP6=DA(3,1)*HHX(III)+DA(3,2)*HHY(III)+DA(3,3)*HHZ(III)
      NORM=TEMP4**2+TEMP5**2+TEMP6**2
      NORM=SQRT(NORM)
      TEMP4=TEMP4*(HNORM/NORM)
      TEMP5=TEMP5*(HNORM/NORM)
      TEMP6=TEMP6*(HNORM/NORM)
C
C Compute new DA, where DA represents the overall rotation
C   applied to this body, ie DA = DA_{inboard} * AO_{relative}
      CALL ROTMUL2(DA,AO,A)
C
      DO J=1,3
         DO I=1,3
            DA(I,J)=A(I,J)
         END DO
      END DO
C
C Get body velocities
      YY(1,III)=YY(1,IIIM)+(-ZCM(IIIM)+ZCM(III))*
     @     YY(5,IIIM)+(YCM(IIIM)-YCM(III))*YY(6,IIIM)+
     @     ((TEMP6*TEMP2-TEMP5*TEMP3)/HNORM)*QDOT(III)
      YY(2,III)=YY(2,IIIM)+(ZCM(IIIM)-ZCM(III))*
     @     YY(4,IIIM)+(-XCM(IIIM)+XCM(III))*YY(6,IIIM)+
     @     ((TEMP4*TEMP3-TEMP6*TEMP1)/HNORM)*QDOT(III)
      YY(3,III)=YY(3,IIIM)+(-YCM(IIIM)+YCM(III))*
     @     YY(4,IIIM)+(XCM(IIIM)-XCM(III))*YY(5,IIIM)+
     @     ((TEMP5*TEMP1-TEMP4*TEMP2)/HNORM)*QDOT(III)
      YY(4,III)=YY(4,IIIM)+QDOT(III)*TEMP4/HNORM
      YY(5,III)=YY(5,IIIM)+QDOT(III)*TEMP5/HNORM
      YY(6,III)=YY(6,IIIM)+QDOT(III)*TEMP6/HNORM
C
C Update atom velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      END DO
C
C Put in a special test in case of a single atom chain tip
      IF ((SIZE(III).EQ.0).AND.
     @     COEFFO(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)).EQ.ONE) THEN
      IIIM=III
      III=GRPSEQ(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)+1)
      POINTER=POINTR(III)+1
      XREL=X(LIST(POINTER))-XCM(IIIM)
      YREL=Y(LIST(POINTER))-YCM(IIIM)
      ZREL=Z(LIST(POINTER))-ZCM(IIIM)
      XV(LIST(POINTER))=YY(1,IIIM)+YY(5,IIIM)*ZREL
     @     -YY(6,IIIM)*YREL
      YV(LIST(POINTER))=YY(2,IIIM)+YY(6,IIIM)*XREL
     @     -YY(4,IIIM)*ZREL
      ZV(LIST(POINTER))=YY(3,IIIM)+YY(4,IIIM)*YREL
     @     -YY(5,IIIM)*XREL
      END IF
C
C
C Recover the DA rotation matrix corresponding
C   to the base of the next chain
      IF (ICHAIN.NE.NCHAINS(ITREE)) THEN
      III=GRPSEQ(ITREE,ICHAIN+1,1)
C
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      CALL ROTMUL1(A,AO,DA)
C
      END IF
C
      END DO
C
      END DO
C
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
C
C Now do T-coupling if required
      IF (TCOUPL) THEN
      IF (RENR(SSTEMP).GT.RSMALL) THEN
      DO I=1,NATOM
      FRICT(I) = FBETA(I)*TIMFAC*AMASS(I)*(ONE-TBATH/RENR(SSTEMP))
      END DO
      END IF
      END IF
C
C      IF (VSCALE) THEN
C        CALL TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
C      END IF
C
C Compute velocity transform matrices, forces, etc ...
      CALL SWEEPOUT1(XCM,YCM,ZCM,FRICT,YY,QQ,MM,
     @     B1,B2,DD,GRPSEQ,SIZE,LIST,POINTR,NIN,
     @     NOUT,COEFFI,COEFFO,NCHAINS,CHNLEN,QDOT)
C
C Do inward recursion to get the KK and LL matrices
C  - these will allow calculation of the base body acceleration
      CALL SWEEPIN(MM,KK,B1,B2,DD,QQ,LL,GRPSEQ,
     @     NCHAINS,CHNLEN)
C
C Get base acceleration, and map it out the tree, including
C  relative accelerations along the way
      CALL SWEEPOUT2(GRPSEQ,KK,MM,B1,B2,DD,
     @     QQ,LL,YYDOT,NCHAINS,CHNLEN,QDDOT)
C
C Now take trial steps and update coordinates, velocities.
C Also compute the intermediate Runge Kutta predictions
      DO ITREE=1,NTREE
C
C First do base
      III=GRPSEQ(ITREE,1,1)
C
      CALL AV2Q(Q0(III),Q1(III),Q2(III),Q3(III),YY(1,III),QVEL)
C
      CALL AA2Q(Q0(III),Q1(III),Q2(III),Q3(III),
     &          YY(1,III),YYDOT(1,III),QVEL,QACC)
C
C Compute R-K trial 3
      K(1)=DTT*YYDOT(1,III)
      K(2)=DTT*YYDOT(2,III)
      K(3)=DTT*YYDOT(3,III)
      K(4)=DTT*QACC(0)
      K(5)=DTT*QACC(1)
      K(6)=DTT*QACC(2)
      K(7)=DTT*QACC(3)
C
C Update the Runge Kutta final predictor
      DO J=1,7
         RKCRD(J,ITREE)=RKCRD(J,ITREE)+K(J)
         RKVEL(J,ITREE)=RKVEL(J,ITREE)+TWO*K(J)
      END DO
C
C Step to an endpoint
C
C Recompute the initial angular velocities -
C   we need them for the runge kutta integration
      CALL AV2Q(Q0O(III),Q1O(III),Q2O(III),Q3O(III),YYO(1,III),QVEL)
C
      XCM(III)=XCMO(ITREE)+DTT*YYO(1,III)+HALF*DTT*K(1)
      YCM(III)=YCMO(ITREE)+DTT*YYO(2,III)+HALF*DTT*K(2)
      ZCM(III)=ZCMO(ITREE)+DTT*YYO(3,III)+HALF*DTT*K(3)
      Q0(III)=Q0O(III)+DTT*QVEL(0)+HALF*DTT*K(4)
      Q1(III)=Q1O(III)+DTT*QVEL(1)+HALF*DTT*K(5)
      Q2(III)=Q2O(III)+DTT*QVEL(2)+HALF*DTT*K(6)
      Q3(III)=Q3O(III)+DTT*QVEL(3)+HALF*DTT*K(7)
C
C Normalize
      QSQR=Q0(III)*Q0(III)+Q1(III)*Q1(III)+Q2(III)*
     @     Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
      YY(1,III)=YYO(1,III)+K(1)
      YY(2,III)=YYO(2,III)+K(2)
      YY(3,III)=YYO(3,III)+K(3)
C
      QVEL(0)=QVEL(0)+K(4)
      QVEL(1)=QVEL(1)+K(5)
      QVEL(2)=QVEL(2)+K(6)
      QVEL(3)=QVEL(3)+K(7)
C
      YY(4,III)=TWO*(-QVEL(0)*Q1(III)+QVEL(1)*Q0(III)-QVEL(2)*
     @     Q3(III)+QVEL(3)*Q2(III))
      YY(5,III)=TWO*(-QVEL(0)*Q2(III)+QVEL(2)*Q0(III)+QVEL(1)*
     @     Q3(III)-QVEL(3)*Q1(III))
      YY(6,III)=TWO*(-QVEL(1)*Q2(III)+QVEL(2)*Q1(III)-QVEL(0)*
     @     Q3(III)+QVEL(3)*Q0(III))
C
C Update atom coordinates, velocities
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
C Get rotation matrix which describes the step
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL1(A,AO,DA)
C
C Now go on to all other bodies
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
C
C Compute R-K trial 3
      K(1)=DTT*QDDOT(III)
C
C Update the Runge Kutta final predictor
      RKQ(III)=RKQ(III)+K(1)
      RKQDOT(III)=RKQDOT(III)+TWO*K(1)
C
C Take a step to an endpoint
      Q(III)=QO(III)+DTT*QDOTO(III)+HALF*DTT*K(1)
      QDOT(III)=QDOTO(III)+K(1)
C
C Rotate body about the old bond direction
      HNORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      HNORM=SQRT(HNORM)
      E0=COS(HALF*(Q(III)-QO(III)))
      E1=(HHX(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E2=(HHY(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E3=(HHZ(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      QSQR=E0**2+E1**2+E2**2+E3**2
      NORM=SQRT(QSQR)
      E0=E0/NORM
      E1=E1/NORM
      E2=E2/NORM
      E3=E3/NORM
C
      CALL TQ2ROT(E0,E1,E2,E3,AO)
C
      TEMP1=SSX(III)
      TEMP2=SSY(III)
      TEMP3=SSZ(III)
      SNORM=TEMP1**2+TEMP2**2+TEMP3**2
      SNORM=SQRT(SNORM)
      TEMP4=AO(1,1)*TEMP1+AO(1,2)*TEMP2+AO(1,3)*TEMP3
      TEMP5=AO(2,1)*TEMP1+AO(2,2)*TEMP2+AO(2,3)*TEMP3
      TEMP6=AO(3,1)*TEMP1+AO(3,2)*TEMP2+AO(3,3)*TEMP3
C
C Now rotate HH, SS by DA, the inboard body rotation
      TEMP1=DA(1,1)*TEMP4+DA(1,2)*TEMP5+DA(1,3)*TEMP6
      TEMP2=DA(2,1)*TEMP4+DA(2,2)*TEMP5+DA(2,3)*TEMP6
      TEMP3=DA(3,1)*TEMP4+DA(3,2)*TEMP5+DA(3,3)*TEMP6
      NORM=TEMP1**2+TEMP2**2+TEMP3**2
      NORM=SQRT(NORM)
      TEMP1=TEMP1*(SNORM/NORM)
      TEMP2=TEMP2*(SNORM/NORM)
      TEMP3=TEMP3*(SNORM/NORM)
      TEMP4=DA(1,1)*HHX(III)+DA(1,2)*HHY(III)+DA(1,3)*HHZ(III)
      TEMP5=DA(2,1)*HHX(III)+DA(2,2)*HHY(III)+DA(2,3)*HHZ(III)
      TEMP6=DA(3,1)*HHX(III)+DA(3,2)*HHY(III)+DA(3,3)*HHZ(III)
      NORM=TEMP4**2+TEMP5**2+TEMP6**2
      NORM=SQRT(NORM)
      TEMP4=TEMP4*(HNORM/NORM)
      TEMP5=TEMP5*(HNORM/NORM)
      TEMP6=TEMP6*(HNORM/NORM)
C
C Get center of mass position
      NNOUT=NOUT(ITREE,ICHAIN,IGP-1)
      XCM(III)=X(NNOUT)+TEMP4-TEMP1
      YCM(III)=Y(NNOUT)+TEMP5-TEMP2
      ZCM(III)=Z(NNOUT)+TEMP6-TEMP3
C
C Compute new DA, where DA represents the overall rotation
C   applied to this body, ie DA = DA_{inboard} * AO_{relative}
      CALL ROTMUL2(DA,AO,A)
C
      DO J=1,3
         DO I=1,3
            DA(I,J)=A(I,J)
         END DO
      END DO
C
C Get new quaternions to describe the body frame orientation
C   The new matrix from body frame to lab frame is A = DA * AO^T_{lab->body}
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL3(DA,AO,A)
C
      Q02=QUART*(ONE+A(1,1)+A(2,2)+A(3,3))
      Q12=QUART*(ONE+A(1,1)-A(2,2)-A(3,3))
      Q22=QUART*(ONE-A(1,1)+A(2,2)-A(3,3))
      Q32=QUART*(ONE-A(1,1)-A(2,2)+A(3,3))
      IF (ABS(Q02).GT.R4SMAL) THEN
      Q0(III)=SQRT(Q02)
      Q1(III)=QUART*(A(2,3)-A(3,2))/Q0(III)
      Q2(III)=QUART*(A(3,1)-A(1,3))/Q0(III)
      Q3(III)=QUART*(A(1,2)-A(2,1))/Q0(III)
      ELSE IF (ABS(Q12).GT.R4SMAL) THEN
      Q1(III)=SQRT(Q12)
      Q0(III)=QUART*(A(2,3)-A(3,2))/Q1(III)
      Q2(III)=QUART*(A(1,2)+A(2,1))/Q1(III)
      Q3(III)=QUART*(A(1,3)+A(3,1))/Q1(III)
      ELSE IF (ABS(Q22).GT.R4SMAL) THEN
      Q2(III)=SQRT(Q22)
      Q0(III)=QUART*(A(3,1)-A(1,3))/Q2(III)
      Q1(III)=QUART*(A(1,2)+A(2,1))/Q2(III)
      Q3(III)=QUART*(A(2,3)+A(3,2))/Q2(III)
      ELSE IF (ABS(Q32).GT.R4SMAL) THEN
      Q3(III)=SQRT(Q32)
      Q0(III)=QUART*(A(1,2)-A(2,1))/Q3(III)
      Q1(III)=QUART*(A(1,3)+A(3,1))/Q3(III)
      Q2(III)=QUART*(A(2,3)+A(3,2))/Q3(III)
      ELSE
      WRITE(6,'(A)') 'ERROR: matrix is not a rotation matrix'
      END IF
      QSQR= Q0(III)*Q0(III)+Q1(III)*Q1(III)
     &     +Q2(III)*Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
C Update coordinates
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      END DO
C
C Get body velocities
      YY(1,III)=YY(1,IIIM)+(-ZCM(IIIM)+ZCM(III))*
     @     YY(5,IIIM)+(YCM(IIIM)-YCM(III))*YY(6,IIIM)+
     @     ((TEMP6*TEMP2-TEMP5*TEMP3)/HNORM)*QDOT(III)
      YY(2,III)=YY(2,IIIM)+(ZCM(IIIM)-ZCM(III))*
     @     YY(4,IIIM)+(-XCM(IIIM)+XCM(III))*YY(6,IIIM)+
     @     ((TEMP4*TEMP3-TEMP6*TEMP1)/HNORM)*QDOT(III)
      YY(3,III)=YY(3,IIIM)+(-YCM(IIIM)+YCM(III))*
     @     YY(4,IIIM)+(XCM(IIIM)-XCM(III))*YY(5,IIIM)+
     @     ((TEMP5*TEMP1-TEMP4*TEMP2)/HNORM)*QDOT(III)
      YY(4,III)=YY(4,IIIM)+QDOT(III)*TEMP4/HNORM
      YY(5,III)=YY(5,IIIM)+QDOT(III)*TEMP5/HNORM
      YY(6,III)=YY(6,IIIM)+QDOT(III)*TEMP6/HNORM
C
C Update atom velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      END DO
C
C Put in a special test in case of a single atom chain tip
      IF ((SIZE(III).EQ.0).AND.
     @     COEFFO(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)).EQ.ONE) THEN
      IIIM=III
      III=GRPSEQ(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)+1)
      POINTER=POINTR(III)+1
      X(LIST(POINTER))=XCM(IIIM)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(IIIM)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(IIIM)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(IIIM)
      YREL=Y(LIST(POINTER))-YCM(IIIM)
      ZREL=Z(LIST(POINTER))-ZCM(IIIM)
      XV(LIST(POINTER))=YY(1,IIIM)+YY(5,IIIM)*ZREL
     @     -YY(6,IIIM)*YREL
      YV(LIST(POINTER))=YY(2,IIIM)+YY(6,IIIM)*XREL
     @     -YY(4,IIIM)*ZREL
      ZV(LIST(POINTER))=YY(3,IIIM)+YY(4,IIIM)*YREL
     @     -YY(5,IIIM)*XREL
      END IF
C
C Recover the DA rotation matrix corresponding
C   to the base of the next chain
      IF (ICHAIN.NE.NCHAINS(ITREE)) THEN
      III=GRPSEQ(ITREE,ICHAIN+1,1)
C
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      CALL ROTMUL1(A,AO,DA)
C
      END IF
C
      END DO
C
      END DO
C
      CALL ENERGY
C
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
C
C Now do T-coupling if required
      IF (TCOUPL) THEN
      IF (RENR(SSTEMP).GT.RSMALL) THEN
      DO I=1,NATOM
      FRICT(I) = FBETA(I)*TIMFAC*AMASS(I)*(ONE-TBATH/RENR(SSTEMP))
      END DO
      END IF
      END IF
C
C      IF (VSCALE) THEN
C        CALL TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
C      END IF
C
C Compute velocity transform matrices, forces, etc ...
      CALL SWEEPOUT1(XCM,YCM,ZCM,FRICT,YY,QQ,MM,
     @     B1,B2,DD,GRPSEQ,SIZE,LIST,POINTR,NIN,
     @     NOUT,COEFFI,COEFFO,NCHAINS,CHNLEN,QDOT)
C
C Do inward recursion to get the KK and LL matrices
C  - these will allow calculation of the base body acceleration
      CALL SWEEPIN(MM,KK,B1,B2,DD,QQ,LL,GRPSEQ,
     @     NCHAINS,CHNLEN)
C
C Get base acceleration, and map it out the tree, including
C  relative accelerations along the way
      CALL SWEEPOUT2(GRPSEQ,KK,MM,B1,B2,DD,
     @     QQ,LL,YYDOT,NCHAINS,CHNLEN,QDDOT)
C
C Now take trial steps and update coordinates, velocities.
C Also compute the intermediate Runge Kutta predictions
      DO ITREE=1,NTREE
C
C First do base
      III=GRPSEQ(ITREE,1,1)
C
      CALL AV2Q(Q0(III),Q1(III),Q2(III),Q3(III),YY(1,III),QVEL)
C
      CALL AA2Q(Q0(III),Q1(III),Q2(III),Q3(III),
     &          YY(1,III),YYDOT(1,III),QVEL,QACC)
C
C Compute R-K trial 4
      K(1)=DTT*YYDOT(1,III)
      K(2)=DTT*YYDOT(2,III)
      K(3)=DTT*YYDOT(3,III)
      K(4)=DTT*QACC(0)
      K(5)=DTT*QACC(1)
      K(6)=DTT*QACC(2)
      K(7)=DTT*QACC(3)
C
C Update the Runge Kutta final predictor
      DO J=1,7
         RKVEL(J,ITREE)=RKVEL(J,ITREE)+K(J)
      END DO
C
C Step to the correct endpoint
C
C Recompute the initial angular velocities -
C   we need them for the runge kutta integration
      CALL AV2Q(Q0O(III),Q1O(III),Q2O(III),Q3O(III),YYO(1,III),QVEL)
C
      XCM(III)=XCMO(ITREE)+DTT*(YYO(1,III)+SIXTH*RKCRD(1,ITREE))
      YCM(III)=YCMO(ITREE)+DTT*(YYO(2,III)+SIXTH*RKCRD(2,ITREE))
      ZCM(III)=ZCMO(ITREE)+DTT*(YYO(3,III)+SIXTH*RKCRD(3,ITREE))
C
      Q0(III)=Q0O(III)+DTT*(QVEL(0)+SIXTH*RKCRD(4,ITREE))
      Q1(III)=Q1O(III)+DTT*(QVEL(1)+SIXTH*RKCRD(5,ITREE))
      Q2(III)=Q2O(III)+DTT*(QVEL(2)+SIXTH*RKCRD(6,ITREE))
      Q3(III)=Q3O(III)+DTT*(QVEL(3)+SIXTH*RKCRD(7,ITREE))
C
C Normalize
      QSQR=Q0(III)*Q0(III)+Q1(III)*Q1(III)+Q2(III)*
     @     Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
      YY(1,III)=YYO(1,III)+SIXTH*RKVEL(1,ITREE)
      YY(2,III)=YYO(2,III)+SIXTH*RKVEL(2,ITREE)
      YY(3,III)=YYO(3,III)+SIXTH*RKVEL(3,ITREE)
C
      QVEL(0)=QVEL(0)+SIXTH*RKVEL(4,ITREE)
      QVEL(1)=QVEL(1)+SIXTH*RKVEL(5,ITREE)
      QVEL(2)=QVEL(2)+SIXTH*RKVEL(6,ITREE)
      QVEL(3)=QVEL(3)+SIXTH*RKVEL(7,ITREE)
C
      YY(4,III)=TWO*(-QVEL(0)*Q1(III)+QVEL(1)*Q0(III)-QVEL(2)*
     @     Q3(III)+QVEL(3)*Q2(III))
      YY(5,III)=TWO*(-QVEL(0)*Q2(III)+QVEL(2)*Q0(III)+QVEL(1)*
     @     Q3(III)-QVEL(3)*Q1(III))
      YY(6,III)=TWO*(-QVEL(1)*Q2(III)+QVEL(2)*Q1(III)-QVEL(0)*
     @     Q3(III)+QVEL(3)*Q0(III))
C
C Update atom coordinates, velocities
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
C Get rotation matrix which describes the step
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL1(A,AO,DA)
C
C Now go on to other bodies
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
C
C Compute R-K trial 4
      K(1)=DTT*QDDOT(III)
C
C Update the Runge Kutta final predictor
      RKQDOT(III)=RKQDOT(III)+K(1)
C
C Step to the correct endpoint
      Q(III)=QO(III)+DTT*(QDOTO(III)+SIXTH*RKQ(III))
      QDOT(III)=QDOTO(III)+SIXTH*RKQDOT(III)
C
C Rotate body about the old bond direction
      HNORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      HNORM=SQRT(HNORM)
      E0=COS(HALF*(Q(III)-QO(III)))
      E1=(HHX(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E2=(HHY(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      E3=(HHZ(III)/HNORM)*SIN(HALF*(Q(III)-QO(III)))
      QSQR=E0**2+E1**2+E2**2+E3**2
      NORM=SQRT(QSQR)
      E0=E0/NORM
      E1=E1/NORM
      E2=E2/NORM
      E3=E3/NORM
C
      CALL TQ2ROT(E0,E1,E2,E3,AO)
C
      TEMP1=SSX(III)
      TEMP2=SSY(III)
      TEMP3=SSZ(III)
      SNORM=TEMP1**2+TEMP2**2+TEMP3**2
      SNORM=SQRT(SNORM)
      SSX(III)=AO(1,1)*TEMP1+AO(1,2)*TEMP2+AO(1,3)*TEMP3
      SSY(III)=AO(2,1)*TEMP1+AO(2,2)*TEMP2+AO(2,3)*TEMP3
      SSZ(III)=AO(3,1)*TEMP1+AO(3,2)*TEMP2+AO(3,3)*TEMP3
C
C Now rotate HH, SS by DA, the inboard body rotation
      TEMP1=SSX(III)
      TEMP2=SSY(III)
      TEMP3=SSZ(III)
      SSX(III)=DA(1,1)*TEMP1+DA(1,2)*TEMP2+DA(1,3)*TEMP3
      SSY(III)=DA(2,1)*TEMP1+DA(2,2)*TEMP2+DA(2,3)*TEMP3
      SSZ(III)=DA(3,1)*TEMP1+DA(3,2)*TEMP2+DA(3,3)*TEMP3
      NORM=SSX(III)**2+SSY(III)**2+SSZ(III)**2
      NORM=SQRT(NORM)
      SSX(III)=SSX(III)*(SNORM/NORM)
      SSY(III)=SSY(III)*(SNORM/NORM)
      SSZ(III)=SSZ(III)*(SNORM/NORM)
      TEMP1=HHX(III)
      TEMP2=HHY(III)
      TEMP3=HHZ(III)
      HHX(III)=DA(1,1)*TEMP1+DA(1,2)*TEMP2+DA(1,3)*TEMP3
      HHY(III)=DA(2,1)*TEMP1+DA(2,2)*TEMP2+DA(2,3)*TEMP3
      HHZ(III)=DA(3,1)*TEMP1+DA(3,2)*TEMP2+DA(3,3)*TEMP3
      NORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      NORM=SQRT(NORM)
      HHX(III)=HHX(III)*(HNORM/NORM)
      HHY(III)=HHY(III)*(HNORM/NORM)
      HHZ(III)=HHZ(III)*(HNORM/NORM)
C
C Get center of mass position
      NNOUT=NOUT(ITREE,ICHAIN,IGP-1)
      XCM(III)=X(NNOUT)+HHX(III)-SSX(III)
      YCM(III)=Y(NNOUT)+HHY(III)-SSY(III)
      ZCM(III)=Z(NNOUT)+HHZ(III)-SSZ(III)
C
C Compute new DA, where DA represents the overall rotation
C   applied to this body, ie DA = DA_{inboard} * AO_{relative}
      CALL ROTMUL2(DA,AO,A)
C
      DO J=1,3
         DO I=1,3
            DA(I,J)=A(I,J)
         END DO
      END DO
C
C Get new quaternions to describe the body frame orientation
C   The new matrix from body frame to lab frame is A = DA * AO^T_{lab->body}
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL ROTMUL3(DA,AO,A)
C
      Q02=QUART*(ONE+A(1,1)+A(2,2)+A(3,3))
      Q12=QUART*(ONE+A(1,1)-A(2,2)-A(3,3))
      Q22=QUART*(ONE-A(1,1)+A(2,2)-A(3,3))
      Q32=QUART*(ONE-A(1,1)-A(2,2)+A(3,3))
      IF (ABS(Q02).GT.R4SMAL) THEN
      Q0(III)=SQRT(Q02)
      Q1(III)=QUART*(A(2,3)-A(3,2))/Q0(III)
      Q2(III)=QUART*(A(3,1)-A(1,3))/Q0(III)
      Q3(III)=QUART*(A(1,2)-A(2,1))/Q0(III)
      ELSE IF (ABS(Q12).GT.R4SMAL) THEN
      Q1(III)=SQRT(Q12)
      Q0(III)=QUART*(A(2,3)-A(3,2))/Q1(III)
      Q2(III)=QUART*(A(1,2)+A(2,1))/Q1(III)
      Q3(III)=QUART*(A(1,3)+A(3,1))/Q1(III)
      ELSE IF (ABS(Q22).GT.R4SMAL) THEN
      Q2(III)=SQRT(Q22)
      Q0(III)=QUART*(A(3,1)-A(1,3))/Q2(III)
      Q1(III)=QUART*(A(1,2)+A(2,1))/Q2(III)
      Q3(III)=QUART*(A(2,3)+A(3,2))/Q2(III)
      ELSE IF (ABS(Q32).GT.R4SMAL) THEN
      Q3(III)=SQRT(Q32)
      Q0(III)=QUART*(A(1,2)-A(2,1))/Q3(III)
      Q1(III)=QUART*(A(1,3)+A(3,1))/Q3(III)
      Q2(III)=QUART*(A(2,3)+A(3,2))/Q3(III)
      ELSE
      WRITE(6,'(A)') 'ERROR: matrix is not a rotation matrix'
      END IF
      QSQR= Q0(III)*Q0(III)+Q1(III)*Q1(III)
     &     +Q2(III)*Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
C
C Update coordinates
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      END DO
C
C Get body velocities
      YY(1,III)=YY(1,IIIM)+(-ZCM(IIIM)+ZCM(III))*
     @     YY(5,IIIM)+(YCM(IIIM)-YCM(III))*YY(6,IIIM)+
     @     ((HHZ(III)*SSY(III)-HHY(III)*SSZ(III))/HNORM)*QDOT(III)
      YY(2,III)=YY(2,IIIM)+(ZCM(IIIM)-ZCM(III))*
     @     YY(4,IIIM)+(-XCM(IIIM)+XCM(III))*YY(6,IIIM)+
     @     ((HHX(III)*SSZ(III)-HHZ(III)*SSX(III))/HNORM)*QDOT(III)
      YY(3,III)=YY(3,IIIM)+(-YCM(IIIM)+YCM(III))*
     @     YY(4,IIIM)+(XCM(IIIM)-XCM(III))*YY(5,IIIM)+
     @     ((HHY(III)*SSX(III)-HHX(III)*SSY(III))/HNORM)*QDOT(III)
      YY(4,III)=YY(4,IIIM)+QDOT(III)*HHX(III)/HNORM
      YY(5,III)=YY(5,IIIM)+QDOT(III)*HHY(III)/HNORM
      YY(6,III)=YY(6,IIIM)+QDOT(III)*HHZ(III)/HNORM
C
C Update atom velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      END DO
C
C Put in a special test in case of a single atom chain tip
      IF ((SIZE(III).EQ.0).AND.
     @     COEFFO(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)).EQ.ONE) THEN
      IIIM=III
      III=GRPSEQ(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)+1)
      POINTER=POINTR(III)+1
      X(LIST(POINTER))=XCM(IIIM)+XB(LIST(POINTER))*A(1,1)
     @     +YB(LIST(POINTER))*A(2,1)+ZB(LIST(POINTER))*A(3,1)
      Y(LIST(POINTER))=YCM(IIIM)+XB(LIST(POINTER))*A(1,2)
     @     +YB(LIST(POINTER))*A(2,2)+ZB(LIST(POINTER))*A(3,2)
      Z(LIST(POINTER))=ZCM(IIIM)+XB(LIST(POINTER))*A(1,3)
     @     +YB(LIST(POINTER))*A(2,3)+ZB(LIST(POINTER))*A(3,3)
      XREL=X(LIST(POINTER))-XCM(IIIM)
      YREL=Y(LIST(POINTER))-YCM(IIIM)
      ZREL=Z(LIST(POINTER))-ZCM(IIIM)
      XV(LIST(POINTER))=YY(1,IIIM)+YY(5,IIIM)*ZREL
     @     -YY(6,IIIM)*YREL
      YV(LIST(POINTER))=YY(2,IIIM)+YY(6,IIIM)*XREL
     @     -YY(4,IIIM)*ZREL
      ZV(LIST(POINTER))=YY(3,IIIM)+YY(4,IIIM)*YREL
     @     -YY(5,IIIM)*XREL
      END IF
C
C Recover the DA rotation matrix corresponding
C   to the base of the next chain
      IF (ICHAIN.NE.NCHAINS(ITREE)) THEN
      III=GRPSEQ(ITREE,ICHAIN+1,1)
C
      CALL Q2ROT(Q0O(III),Q1O(III),Q2O(III),Q3O(III),AO)
C
      CALL Q2ROT(Q0(III),Q1(III),Q2(III),Q3(III),A)
C
      CALL ROTMUL1(A,AO,DA)
C
      END IF
C
      END DO
C
      END DO
C
C Coordinate write options
      IF (NSAVC.GT.0.AND.IUNCRD.GT.0) THEN
      IF (MOD(ISTEP,NSAVC).EQ.0) THEN
      CALL WRITTC('CORD',NATOM,X,Y,Z,CRDNUM,CRDIND,ISTEP,QFIRSC,
     @     DTT,NSAVC,IUNCRD,QFORM,FORM,SCALE,OFFSET)
      END IF
      END IF
C
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
      IF (VSCALE) THEN
        CALL TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
      END IF
      CALL TDYNKIN(NDEGF,XV,YV,ZV)
C
      IF (NPRINT.GT.0) THEN
      IF (MOD(ISTEP,NPRINT).EQ.0) THEN
      WRITE(6,'(A,I9,A,F12.5,A)')
     @     ' --------------- step=',ISTEP,' at ',
     @     TIMFAC*DTT*ISTEP,' ps -----------------------------'
      CALL PRINTD(RENR)
      IF (WRNLEV.GE.10) THEN
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.,
     @     CMX,CMY,CMZ,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
      END IF
      END IF
      END IF
C
      END DO
C
      IF (IUNCRD.GT.0) CALL VCLOSE(IUNCRD,'KEEP',ERROR)
C
      RETURN
      END
C===============================================================
      SUBROUTINE RTSTOP(LIST,POINTR,NCHAINS,CHNLEN,
     @     YY,GRPSEQ,XCM,YCM,ZCM,COEFFO)
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'timer.inc'
C
      INTEGER NCHAINS(*)
      INTEGER CHNLEN(MAXTREE,*),GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER LIST(*),POINTR(*)
      DOUBLE PRECISION YY(6,NGP)
      DOUBLE PRECISION XCM(*),YCM(*),ZCM(*),COEFFO(MAXTREE,MAXCHN,*)
C
C Local
      INTEGER I,ICHAIN,IGP,III,IIIM,POINTER,ITREE
      DOUBLE PRECISION TCM(3,3),LH(3),MH(3),D,ONE
      DOUBLE PRECISION IXX,IXY,IXZ,IYY,IYZ,IZZ,XI,YI,ZI
      DOUBLE PRECISION CMX,CMY,CMZ,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM
      DOUBLE PRECISION OXCM,OYCM,OZCM,XREL,YREL,ZREL
      PARAMETER (ONE=1.0D0)
C
C
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.,
     @     CMX,CMY,CMZ,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
C
C
      IXX = 0.0D0
      IXY = 0.0D0
      IXZ = 0.0D0
      IYY = 0.0D0
      IYZ = 0.0D0
      IZZ = 0.0D0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) THEN
      XI = X(I) - CMX
      YI = Y(I) - CMY
      ZI = Z(I) - CMZ
      IXX = IXX + XI*XI*AMASS(I)
      IXY = IXY + XI*YI*AMASS(I)
      IXZ = IXZ + XI*ZI*AMASS(I)
      IYY = IYY + YI*YI*AMASS(I)
      IYZ = IYZ + YI*ZI*AMASS(I)
      IZZ = IZZ + ZI*ZI*AMASS(I)
      END IF
      END DO
      TCM(1,1) = IYY + IZZ
      TCM(2,1) = -IXY
      TCM(3,1) = -IXZ
      TCM(1,2) = -IXY
      TCM(2,2) = IXX + IZZ
      TCM(3,2) = -IYZ
      TCM(1,3) = -IXZ
      TCM(2,3) = -IYZ
      TCM(3,3) = IXX + IYY
C
C invert the inertia tensor TCM
      CALL MINVV(TCM,3,D,LH,MH)
C
C get angular velocity OXCM, OYCM, OZCM
      OXCM = AXCM*TCM(1,1) + AYCM*TCM(1,2) +  AZCM*TCM(1,3)
      OYCM = AXCM*TCM(2,1) + AYCM*TCM(2,2) +  AZCM*TCM(2,3)
      OZCM = AXCM*TCM(3,1) + AYCM*TCM(3,2) +  AZCM*TCM(3,3)
C
C Loop over trees
      DO ITREE=1,NTREE
C
C remove CM translational and rotational motion from velocities
      III=GRPSEQ(ITREE,1,1)
      YY(1,III)=YY(1,III)-VXCM-OYCM*(ZCM(III)-CMZ)
     @     +OZCM*(YCM(III)-CMY)
      YY(2,III)=YY(2,III)-VYCM-OZCM*(XCM(III)-CMX)
     @     +OXCM*(ZCM(III)-CMZ)
      YY(3,III)=YY(3,III)-VZCM-OXCM*(YCM(III)-CMY)
     @     +OYCM*(XCM(III)-CMX)
      YY(4,III)=YY(4,III)-OXCM
      YY(5,III)=YY(5,III)-OYCM
      YY(6,III)=YY(6,III)-OZCM
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      YY(1,III)=YY(1,III)-VXCM-OYCM*(ZCM(III)-CMZ)
     @     +OZCM*(YCM(III)-CMY)
      YY(2,III)=YY(2,III)-VYCM-OZCM*(XCM(III)-CMX)
     @     +OXCM*(ZCM(III)-CMZ)
      YY(3,III)=YY(3,III)-VZCM-OXCM*(YCM(III)-CMY)
     @     +OYCM*(XCM(III)-CMX)
      YY(4,III)=YY(4,III)-OXCM
      YY(5,III)=YY(5,III)-OYCM
      YY(6,III)=YY(6,III)-OZCM
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      END DO
C
      IF (((POINTR(III+1)-POINTR(III)-1).EQ.0).AND.
     @     COEFFO(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)).EQ.ONE) THEN
      IIIM=III
      III=GRPSEQ(ITREE,ICHAIN,CHNLEN(ITREE,ICHAIN)+1)
      POINTER=POINTR(III)+1
      XREL=X(LIST(POINTER))-XCM(IIIM)
      YREL=Y(LIST(POINTER))-YCM(IIIM)
      ZREL=Z(LIST(POINTER))-ZCM(IIIM)
      XV(LIST(POINTER))=YY(1,IIIM)+YY(5,IIIM)*ZREL-YY(6,IIIM)*YREL
      YV(LIST(POINTER))=YY(2,IIIM)+YY(6,IIIM)*XREL-YY(4,IIIM)*ZREL
      ZV(LIST(POINTER))=YY(3,IIIM)+YY(4,IIIM)*YREL-YY(5,IIIM)*XREL
      END IF
C
      END DO
      END DO
C
      WRITE(6,1000)
 1000 FORMAT(' ROTSTP: CM translation and rotation removed.')
      CALL CENMAS(NATOM,X,Y,Z,XV,YV,ZV,AMASS,IMOVE,.TRUE.,
     1     CMX,CMY,CMZ,VXCM,VYCM,VZCM,AXCM,AYCM,AZCM)
C
      RETURN
      END
C===============================================================
      SUBROUTINE SWEEPOUT1(XCM,YCM,ZCM,FRICT,YY,QQ,MM,B1,
     @     B2,DD,GRPSEQ,SIZE,LIST,POINTR,
     @     NIN,NOUT,COEFFI,COEFFO,NCHAINS,CHNLEN,QDOT)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
C I/O
      INTEGER SIZE(*),LIST(*),POINTR(*),NIN(MAXTREE,MAXCHN,*)
      INTEGER NOUT(MAXTREE,MAXCHN,*),GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER NCHAINS(*),CHNLEN(MAXTREE,*)
      DOUBLE PRECISION XCM(*),YCM(*),ZCM(*),FRICT(*)
      DOUBLE PRECISION YY(6,NGP)
      DOUBLE PRECISION MM(6,6,NGP),B1(6,6,NGP)
      DOUBLE PRECISION COEFFI(MAXTREE,MAXCHN,*),QQ(6,NGP)
      DOUBLE PRECISION COEFFO(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION B2(6,NGP),DD(6,NGP),QDOT(*)
C
C Local
      INTEGER III,IIIM,IGP,POINTER,J,K,ICHAIN
      INTEGER NNIN,NNINP,NNOUT,ITREE
      DOUBLE PRECISION A(6),ZERO,XREL,YREL,ZREL,XRELP,YRELP,ZRELP
      DOUBLE PRECISION HHX,HHY,HHZ,SSX,SSY,SSZ,HNORM,ONE
      DOUBLE PRECISION B2DOT(6),B1DOT(6,6),HALF,COEFIN,COEFOU
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
C
C Loop over all trees
      DO ITREE=1,NTREE
C Do Base first
      III=GRPSEQ(ITREE,1,1)
      COEFIN=COEFFI(ITREE,1,1)
      COEFOU=COEFFO(ITREE,1,1)
C
C A is the inertia tensor. It's a 3X3 tensor.  The
C diagonal elements are numbered 1, 3, and 6
      A(1)=ZERO
      A(2)=ZERO
      A(3)=ZERO
      A(4)=ZERO
      A(5)=ZERO
      A(6)=ZERO
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      A(1)=A(1)+AMASS(LIST(POINTER))*(YREL**2+ZREL**2)
      A(3)=A(3)+AMASS(LIST(POINTER))*(XREL**2+ZREL**2)
      A(6)=A(6)+AMASS(LIST(POINTER))*(XREL**2+YREL**2)
      A(2)=A(2)-AMASS(LIST(POINTER))*XREL*YREL
      A(4)=A(4)-AMASS(LIST(POINTER))*ZREL*XREL
      A(5)=A(5)-AMASS(LIST(POINTER))*YREL*ZREL
      END DO
C
      MM(4,4,III)=A(1)
      MM(4,5,III)=A(2)
      MM(4,6,III)=A(4)
      MM(5,4,III)=A(2)
      MM(5,5,III)=A(3)
      MM(5,6,III)=A(5)
      MM(6,4,III)=A(4)
      MM(6,5,III)=A(5)
      MM(6,6,III)=A(6)
C
C Get generalized force
      QQ(4,III)=-A(4)*YY(4,III)*YY(5,III)-A(5)*
     @     (YY(5,III)**2-YY(6,III)**2)+A(2)*YY(4,III)
     @     *YY(6,III)-(A(6)-A(3))*YY(5,III)*YY(6,III)
      QQ(5,III)=+A(4)*(YY(4,III)**2-YY(6,III)**2)+A(5)*
     @     YY(4,III)*YY(5,III)-(A(1)-A(6))*YY(4,III)*YY(6,III)
     @     -A(2)*YY(5,III)*YY(6,III)
      QQ(6,III)=-A(2)*(YY(4,III)**2-YY(5,III)**2)-A(5)*
     @     YY(4,III)*YY(6,III)+(A(1)-A(3))*YY(4,III)*YY(5,III)
     @     +A(4)*YY(5,III)*YY(6,III)
C
      QQ(1,III)=ZERO
      QQ(2,III)=ZERO
      QQ(3,III)=ZERO
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      QQ(1,III)=QQ(1,III)-DX(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     XV(LIST(POINTER))
      QQ(2,III)=QQ(2,III)-DY(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     YV(LIST(POINTER))
      QQ(3,III)=QQ(3,III)-DZ(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     ZV(LIST(POINTER))
      QQ(4,III)=QQ(4,III)-YREL*(DZ(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     ZV(LIST(POINTER)))+ZREL*(DY(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*YV(LIST(POINTER)))
      QQ(5,III)=QQ(5,III)-ZREL*(DX(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     XV(LIST(POINTER)))+XREL*(DZ(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*ZV(LIST(POINTER)))
      QQ(6,III)=QQ(6,III)-XREL*(DY(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     YV(LIST(POINTER)))+YREL*(DX(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*XV(LIST(POINTER)))
      END DO
C
C Go on to other bodies
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
C
      NNIN=NIN(ITREE,ICHAIN,IGP)
      NNINP=NIN(ITREE,ICHAIN,IGP+1)
      NNOUT=NOUT(ITREE,ICHAIN,IGP-1)
      COEFIN=COEFFI(ITREE,ICHAIN,IGP)
      COEFOU=COEFFO(ITREE,ICHAIN,IGP)
C
      III=GRPSEQ(ITREE,ICHAIN,IGP)
C
C Get B1, B2, and their time derivatives; also get DD and QDOT
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
      DO J=1,6
         DO K=1,6
            B1(K,J,III)=ZERO
         END DO
         B1(J,J,III)=ONE
      END DO
      B1(1,5,III)=-ZCM(IIIM)+ZCM(III)
      B1(1,6,III)=YCM(IIIM)-YCM(III)
      B1(2,4,III)=ZCM(IIIM)-ZCM(III)
      B1(2,6,III)=-XCM(IIIM)+XCM(III)
      B1(3,4,III)=-YCM(IIIM)+YCM(III)
      B1(3,5,III)=XCM(IIIM)-XCM(III)
C
C Update the HH and SS vectors
      HHX=X(NNIN)-X(NNOUT)
      HHY=Y(NNIN)-Y(NNOUT)
      HHZ=Z(NNIN)-Z(NNOUT)
      HNORM=HHX**2+HHY**2+HHZ**2
      HNORM=SQRT(HNORM)
      HHX=HHX/HNORM
      HHY=HHY/HNORM
      HHZ=HHZ/HNORM
      SSX=X(NNIN)-XCM(III)
      SSY=Y(NNIN)-YCM(III)
      SSZ=Z(NNIN)-ZCM(III)
C
      B2(1,III)=HHZ*SSY-HHY*SSZ
      B2(2,III)=-HHZ*SSX+HHX*SSZ
      B2(3,III)=HHY*SSX-HHX*SSY
      B2(4,III)=HHX
      B2(5,III)=HHY
      B2(6,III)=HHZ
C
C then the B time derivatives
      DO J=1,6
         DO K=1,6
            B1DOT(K,J)=ZERO
         END DO
      END DO
      B1DOT(1,5)=-YY(3,IIIM)+YY(3,III)
      B1DOT(1,6)=YY(2,IIIM)-YY(2,III)
      B1DOT(2,4)=YY(3,IIIM)-YY(3,III)
      B1DOT(2,6)=-YY(1,IIIM)+YY(1,III)
      B1DOT(3,4)=-YY(2,IIIM)+YY(2,III)
      B1DOT(3,5)=YY(1,IIIM)-YY(1,III)
C
      B2DOT(1)=(HHY*SSY+HHZ*SSZ)*
     @     (YY(4,IIIM)-YY(4,III))-HHX*(SSY*YY(5,IIIM)
     @     +SSZ*YY(6,IIIM))+SSX*(HHY*YY(5,III)
     @     +HHZ*YY(6,III))
      B2DOT(2)=(HHX*SSX+HHZ*SSZ)*
     @     (YY(5,IIIM)-YY(5,III))-HHY*(SSX*YY(4,IIIM)
     @     +SSZ*YY(6,IIIM))+SSY*(HHX*YY(4,III)
     @     +HHZ*YY(6,III))
      B2DOT(3)=(HHX*SSX+HHY*SSY)*
     @     (YY(6,IIIM)-YY(6,III))-HHZ*(SSX*YY(4,IIIM)
     @     +SSY*YY(5,IIIM))+SSZ*(HHX*YY(4,III)
     @     +HHY*YY(5,III))
      B2DOT(4)=HHZ*YY(5,IIIM)-HHY*YY(6,IIIM)
      B2DOT(5)=-HHZ*YY(4,IIIM)+HHX*YY(6,IIIM)
      B2DOT(6)=HHY*YY(4,IIIM)-HHX*YY(5,IIIM)
C
C get DD vector
      DO J=1,6
      DD(J,III)=ZERO
      DO K=1,6
      DD(J,III)=DD(J,III)+B1DOT(J,K)*YY(K,IIIM)
      END DO
      DD(J,III)=DD(J,III)+B2DOT(J)*QDOT(III)
      END DO
C
      IF (SIZE(III).GT.0) THEN
C
      A(1)=ZERO
      A(2)=ZERO
      A(3)=ZERO
      A(4)=ZERO
      A(5)=ZERO
      A(6)=ZERO
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      IF (LIST(POINTER).NE.NNIN) THEN
      A(1)=A(1)+AMASS(LIST(POINTER))*(YREL**2+ZREL**2)
      A(3)=A(3)+AMASS(LIST(POINTER))*(XREL**2+ZREL**2)
      A(6)=A(6)+AMASS(LIST(POINTER))*(XREL**2+YREL**2)
      A(2)=A(2)-AMASS(LIST(POINTER))*XREL*YREL
      A(4)=A(4)-AMASS(LIST(POINTER))*ZREL*XREL
      A(5)=A(5)-AMASS(LIST(POINTER))*YREL*ZREL
      ELSE
      A(1)=A(1)+AMASS(LIST(POINTER))*(YREL**2+ZREL**2)*COEFIN
      A(3)=A(3)+AMASS(LIST(POINTER))*(XREL**2+ZREL**2)*COEFIN
      A(6)=A(6)+AMASS(LIST(POINTER))*(XREL**2+YREL**2)*COEFIN
      A(2)=A(2)-AMASS(LIST(POINTER))*XREL*YREL*COEFIN
      A(4)=A(4)-AMASS(LIST(POINTER))*ZREL*XREL*COEFIN
      A(5)=A(5)-AMASS(LIST(POINTER))*YREL*ZREL*COEFIN
      END IF
      END DO
C
      ELSE
C
      XREL=X(NNIN)-XCM(III)
      YREL=Y(NNIN)-YCM(III)
      ZREL=Z(NNIN)-ZCM(III)
      XRELP=X(NNINP)-XCM(III)
      YRELP=Y(NNINP)-YCM(III)
      ZRELP=Z(NNINP)-ZCM(III)
      A(1)=AMASS(NNIN)*(YREL**2+ZREL**2)*COEFIN+
     @     AMASS(NNINP)*(YRELP**2+ZRELP**2)*COEFOU
      A(3)=AMASS(NNIN)*(XREL**2+ZREL**2)*COEFIN+
     @     AMASS(NNINP)*(XRELP**2+ZRELP**2)*COEFOU
      A(6)=AMASS(NNIN)*(XREL**2+YREL**2)*COEFIN+
     @     AMASS(NNINP)*(XRELP**2+YRELP**2)*COEFOU
      A(2)=-AMASS(NNIN)*XREL*YREL*COEFIN-
     @     AMASS(NNINP)*XRELP*YRELP*COEFOU
      A(4)=-AMASS(NNIN)*ZREL*XREL*COEFIN-
     @     AMASS(NNINP)*ZRELP*XRELP*COEFOU
      A(5)=-AMASS(NNIN)*YREL*ZREL*COEFIN-
     @     AMASS(NNINP)*YRELP*ZRELP*COEFOU
C
      END IF
C
      MM(4,4,III)=A(1)
      MM(4,5,III)=A(2)
      MM(4,6,III)=A(4)
      MM(5,4,III)=A(2)
      MM(5,5,III)=A(3)
      MM(5,6,III)=A(5)
      MM(6,4,III)=A(4)
      MM(6,5,III)=A(5)
      MM(6,6,III)=A(6)
C
C Get generalized force
      QQ(4,III)=-A(4)*YY(4,III)*YY(5,III)-A(5)*
     @     (YY(5,III)**2-YY(6,III)**2)+A(2)*YY(4,III)
     @     *YY(6,III)-(A(6)-A(3))*YY(5,III)*YY(6,III)
      QQ(5,III)=+A(4)*(YY(4,III)**2-YY(6,III)**2)+A(5)*
     @     YY(4,III)*YY(5,III)-(A(1)-A(6))*YY(4,III)*YY(6,III)
     @     -A(2)*YY(5,III)*YY(6,III)
      QQ(6,III)=-A(2)*(YY(4,III)**2-YY(5,III)**2)-A(5)*
     @     YY(4,III)*YY(6,III)+(A(1)-A(3))*YY(4,III)*YY(5,III)
     @     +A(4)*YY(5,III)*YY(6,III)
C
      QQ(1,III)=ZERO
      QQ(2,III)=ZERO
      QQ(3,III)=ZERO
      IF (SIZE(III).GT.0) THEN
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      IF (LIST(POINTER).NE.NNIN) THEN
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      QQ(1,III)=QQ(1,III)-DX(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     XV(LIST(POINTER))
      QQ(2,III)=QQ(2,III)-DY(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     YV(LIST(POINTER))
      QQ(3,III)=QQ(3,III)-DZ(LIST(POINTER))-FRICT(LIST(POINTER))*
     @     ZV(LIST(POINTER))
      QQ(4,III)=QQ(4,III)-YREL*(DZ(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     ZV(LIST(POINTER)))+ZREL*(DY(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*YV(LIST(POINTER)))
      QQ(5,III)=QQ(5,III)-ZREL*(DX(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     XV(LIST(POINTER)))+XREL*(DZ(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*ZV(LIST(POINTER)))
      QQ(6,III)=QQ(6,III)-XREL*(DY(LIST(POINTER))+FRICT(LIST(POINTER))*
     @     YV(LIST(POINTER)))+YREL*(DX(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*XV(LIST(POINTER)))
      ELSE
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      QQ(1,III)=QQ(1,III)-COEFIN*(DX(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*XV(LIST(POINTER)))
      QQ(2,III)=QQ(2,III)-COEFIN*(DY(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*YV(LIST(POINTER)))
      QQ(3,III)=QQ(3,III)-COEFIN*(DZ(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*ZV(LIST(POINTER)))
      QQ(4,III)=QQ(4,III)-YREL*COEFIN*(DZ(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*ZV(LIST(POINTER)))+
     @     ZREL*COEFIN*(DY(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*YV(LIST(POINTER)))
      QQ(5,III)=QQ(5,III)-ZREL*COEFIN*(DX(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*XV(LIST(POINTER)))+
     @     XREL*COEFIN*(DZ(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*ZV(LIST(POINTER)))
      QQ(6,III)=QQ(6,III)-XREL*COEFIN*(DY(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*YV(LIST(POINTER)))+
     @     YREL*COEFIN*(DX(LIST(POINTER))+
     @     FRICT(LIST(POINTER))*XV(LIST(POINTER)))
      END IF
      END DO
C
C small group
      ELSE
      XREL=X(NNIN)-XCM(III)
      YREL=Y(NNIN)-YCM(III)
      ZREL=Z(NNIN)-ZCM(III)
      XRELP=X(NNINP)-XCM(III)
      YRELP=Y(NNINP)-YCM(III)
      ZRELP=Z(NNINP)-ZCM(III)
C
      QQ(1,III)=-DX(NNIN)*COEFIN-DX(NNINP)*COEFOU-FRICT(NNIN)*
     @     COEFIN*XV(NNIN)-FRICT(NNINP)*COEFOU*XV(NNINP)
      QQ(2,III)=-DY(NNIN)*COEFIN-DY(NNINP)*COEFOU-FRICT(NNIN)*
     @     COEFIN*YV(NNIN)-FRICT(NNINP)*COEFOU*YV(NNINP)
      QQ(3,III)=-DZ(NNIN)*COEFIN-DZ(NNINP)*COEFOU-FRICT(NNIN)*
     @     COEFIN*ZV(NNIN)-FRICT(NNINP)*COEFOU*ZV(NNINP)
C
      QQ(4,III)=QQ(4,III)-YREL*COEFIN*(DZ(NNIN)+ZV(NNIN)*
     @     FRICT(NNIN))+ZREL*COEFIN*(DY(NNIN)+
     @     FRICT(NNIN)*YV(NNIN))-YRELP*COEFOU*(DZ(NNINP)+
     @     ZV(NNINP)*FRICT(NNINP))+ZRELP*COEFOU*
     @     (DY(NNINP)+FRICT(NNINP)*YV(NNINP))
      QQ(5,III)=QQ(5,III)-ZREL*COEFIN*(DX(NNIN)+XV(NNIN)*
     @     FRICT(NNIN))+XREL*COEFIN*(DZ(NNIN)+
     @     FRICT(NNIN)*ZV(NNIN))-ZRELP*COEFOU*(DX(NNINP)+
     @     XV(NNINP)*FRICT(NNINP))+XRELP*COEFOU*
     @     (DZ(NNINP)+FRICT(NNINP)*ZV(NNINP))
      QQ(6,III)=QQ(6,III)-XREL*COEFIN*(DY(NNIN)+YV(NNIN)*
     @     FRICT(NNIN))+YREL*COEFIN*(DX(NNIN)+
     @     FRICT(NNIN)*XV(NNIN))-XRELP*COEFOU*(DY(NNINP)+
     @     YV(NNINP)*FRICT(NNINP))+YRELP*COEFOU*
     @     (DX(NNINP)+FRICT(NNINP)*XV(NNINP))
      END IF
C
      END DO
      END DO
C
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE SWEEPIN(MM,KK,B1,B2,DD,QQ,LL,
     @     GRPSEQ,NCHAINS,CHNLEN)
C
C modified by Paul Adams August-98 to optimize operations
C
      IMPLICIT NONE
C
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
C I/O
      INTEGER NCHAINS(*),GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER CHNLEN(MAXTREE,*)
      DOUBLE PRECISION MM(6,6,NGP),KK(6,6,NGP)
      DOUBLE PRECISION QQ(6,NGP),LL(6,NGP),B1(6,6,NGP)
      DOUBLE PRECISION B2(6,NGP),DD(6,NGP)
C
C Local
      INTEGER J,K,IGP,IDUM,IDUMP,III,IIIP,L,ICHAIN,JCHAIN
      INTEGER ITREE
      DOUBLE PRECISION TEMP1,SCRATCH(6),TEMP2(6),TEMP5(6,6)
      DOUBLE PRECISION TEMP3(6,6),TEMP4(6),TEMP,SCRATCH1(6,6)
C parameter
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
C
C
C First set everything to zero
      DO J=1,NGP
         DO K=1,6
            LL(K,J)=ZERO
            DO L=1,6
               KK(L,K,J)=ZERO
            END DO
         END DO
      END DO
C
C Do inward recursion to get base acceleration
      DO ITREE=1,NTREE
      DO ICHAIN=1,NCHAINS(ITREE)
      JCHAIN=NCHAINS(ITREE)-ICHAIN+1
C
      DO IGP=1,CHNLEN(ITREE,JCHAIN)-1
C
      IDUMP=CHNLEN(ITREE,JCHAIN)-IGP+1
      IDUM =CHNLEN(ITREE,JCHAIN)-IGP
      III=GRPSEQ(ITREE,JCHAIN,IDUM)
      IIIP=GRPSEQ(ITREE,JCHAIN,IDUMP)
C
C Assemble useful matrix products
C
C TEMP5 = (M + K)
      DO J=1,6
         TEMP5(J,1)=MM(J,1,IIIP)+KK(J,1,IIIP)
         TEMP5(J,2)=MM(J,2,IIIP)+KK(J,2,IIIP)
         TEMP5(J,3)=MM(J,3,IIIP)+KK(J,3,IIIP)
         TEMP5(J,4)=MM(J,4,IIIP)+KK(J,4,IIIP)
         TEMP5(J,5)=MM(J,5,IIIP)+KK(J,5,IIIP)
         TEMP5(J,6)=MM(J,6,IIIP)+KK(J,6,IIIP)
      END DO
C
C TEMP1 = [B2^T * (M + K) * B2 ]^(-1)
      TEMP1=ZERO
      DO J=1,6
      SCRATCH(J)=TEMP5(J,1)*B2(1,IIIP) +
     &           TEMP5(J,2)*B2(2,IIIP) +
     &           TEMP5(J,3)*B2(3,IIIP) +
     &           TEMP5(J,4)*B2(4,IIIP) +
     &           TEMP5(J,5)*B2(5,IIIP) +
     &           TEMP5(J,6)*B2(6,IIIP)
      TEMP1=TEMP1+B2(J,IIIP)*SCRATCH(J)
      END DO
      TEMP1=ONE/TEMP1
C
C TEMP2 = (M + K) * D - (Q + L)
      DO J=1,6
      TEMP2(J)=( TEMP5(J,1)*DD(1,IIIP) +
     &           TEMP5(J,2)*DD(2,IIIP) +
     &           TEMP5(J,3)*DD(3,IIIP) +
     &           TEMP5(J,4)*DD(4,IIIP) +
     &           TEMP5(J,5)*DD(5,IIIP) +
     &           TEMP5(J,6)*DD(6,IIIP) ) -
     &            QQ(J,IIIP) - LL(J,IIIP)
      END DO
C
C TEMP3 = (M + K) * B1
      DO K=1,6
      DO J=1,6
      TEMP3(J,K)=TEMP5(J,1)*B1(1,K,IIIP) +
     &           TEMP5(J,2)*B1(2,K,IIIP) +
     &           TEMP5(J,3)*B1(3,K,IIIP) +
     &           TEMP5(J,4)*B1(4,K,IIIP) +
     &           TEMP5(J,5)*B1(5,K,IIIP) +
     &           TEMP5(J,6)*B1(6,K,IIIP)
      END DO
      END DO
C
C TEMP4 = B1^T * (M + K) * B2
      DO J=1,6
      TEMP4(J)=(B1(1,J,IIIP)*SCRATCH(1)) +
     &         (B1(2,J,IIIP)*SCRATCH(2)) +
     &         (B1(3,J,IIIP)*SCRATCH(3)) +
     &         (B1(4,J,IIIP)*SCRATCH(4)) +
     &         (B1(5,J,IIIP)*SCRATCH(5)) +
     &         (B1(6,J,IIIP)*SCRATCH(6))
      END DO
C
C Get the KK's and LL's
      DO J=1,6
      SCRATCH(J)=( (B2(1,IIIP)*TEMP3(1,J)) +
     &             (B2(2,IIIP)*TEMP3(2,J)) +
     &             (B2(3,IIIP)*TEMP3(3,J)) +
     &             (B2(4,IIIP)*TEMP3(4,J)) +
     &             (B2(5,IIIP)*TEMP3(5,J)) +
     &             (B2(6,IIIP)*TEMP3(6,J)) ) * TEMP1
      END DO
C
      DO J=1,6
      DO K=1,6
      KK(J,K,III)=KK(J,K,III) + B1(1,J,IIIP)*TEMP3(1,K)
     &                        + B1(2,J,IIIP)*TEMP3(2,K)
     &                        + B1(3,J,IIIP)*TEMP3(3,K)
     &                        + B1(4,J,IIIP)*TEMP3(4,K)
     &                        + B1(5,J,IIIP)*TEMP3(5,K)
     &                        + B1(6,J,IIIP)*TEMP3(6,K)
      END DO
      END DO
C
      DO J=1,6
      SCRATCH1(J,1)=TEMP4(J)*SCRATCH(1)
      SCRATCH1(J,2)=TEMP4(J)*SCRATCH(2)
      SCRATCH1(J,3)=TEMP4(J)*SCRATCH(3)
      SCRATCH1(J,4)=TEMP4(J)*SCRATCH(4)
      SCRATCH1(J,5)=TEMP4(J)*SCRATCH(5)
      SCRATCH1(J,6)=TEMP4(J)*SCRATCH(6)
      END DO
C
      DO J=1,6
      KK(J,1,III)=KK(J,1,III)-SCRATCH1(J,1)
      KK(J,2,III)=KK(J,2,III)-SCRATCH1(J,2)
      KK(J,3,III)=KK(J,3,III)-SCRATCH1(J,3)
      KK(J,4,III)=KK(J,4,III)-SCRATCH1(J,4)
      KK(J,5,III)=KK(J,5,III)-SCRATCH1(J,5)
      KK(J,6,III)=KK(J,6,III)-SCRATCH1(J,6)
      END DO
C
      TEMP=( B2(1,IIIP)*TEMP2(1) +
     &       B2(2,IIIP)*TEMP2(2) +
     &       B2(3,IIIP)*TEMP2(3) +
     &       B2(4,IIIP)*TEMP2(4) +
     &       B2(5,IIIP)*TEMP2(5) +
     &       B2(6,IIIP)*TEMP2(6) ) * TEMP1
C
      DO J=1,6
         LL(J,III)=LL(J,III) - B1(1,J,IIIP)*TEMP2(1)
     &                       - B1(2,J,IIIP)*TEMP2(2)
     &                       - B1(3,J,IIIP)*TEMP2(3)
     &                       - B1(4,J,IIIP)*TEMP2(4)
     &                       - B1(5,J,IIIP)*TEMP2(5)
     &                       - B1(6,J,IIIP)*TEMP2(6)
     &                       + TEMP4(J)*TEMP
      END DO
C
C End inward recursion
      END DO
      END DO
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE SWEEPOUT2(GRPSEQ,KK,MM,B1,
     @     B2,DD,QQ,LL,YYDOT,NCHAINS,CHNLEN,QDDOT)
C
C modified by Paul Adams August-98 to optimize operations
C
      IMPLICIT NONE
C
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
C I/O
      INTEGER NCHAINS(*),GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER CHNLEN(MAXTREE,*)
      DOUBLE PRECISION KK(6,6,NGP),MM(6,6,NGP)
      DOUBLE PRECISION DD(6,NGP),QQ(6,NGP),QDDOT(*)
      DOUBLE PRECISION B1(6,6,NGP),B2(6,NGP)
      DOUBLE PRECISION LL(6,NGP),YYDOT(6,NGP)
C
C Local
      INTEGER IGP,III,J,K,IIIM,ICHAIN,ITREE
      DOUBLE PRECISION TEMP3(6,6),LDUM(6),MDUM(6),DET,TEMP1
      DOUBLE PRECISION SCRATCH(6),TEMP,TEMP5(6,6)
C
      DOUBLE PRECISION ONE, ZERO, SMALL
      PARAMETER (ONE=1.0D0,ZERO=0.0D0,SMALL=1.0D-10)
C
C Loop over all trees
      DO ITREE=1,NTREE
C
C Get base acceleration, and map outwards
      III=GRPSEQ(ITREE,1,1)
      DO K=1,6
      DO J=1,6
      TEMP3(J,K)=KK(J,K,III)+MM(J,K,III)
      END DO
      END DO
C
      CALL MINVV(TEMP3,6,DET,LDUM,MDUM)
C
      IF (DET.LT.SMALL) THEN
      CALL WRNDIE(-5,'SWEEPOUT2','THERE IS A ZERO DETERMINANT')
      END IF
C
      DO J=1,6
      YYDOT(J,III)=TEMP3(J,1)*(LL(1,III)+QQ(1,III)) +
     &             TEMP3(J,2)*(LL(2,III)+QQ(2,III)) +
     &             TEMP3(J,3)*(LL(3,III)+QQ(3,III)) +
     &             TEMP3(J,4)*(LL(4,III)+QQ(4,III)) +
     &             TEMP3(J,5)*(LL(5,III)+QQ(5,III)) +
     &             TEMP3(J,6)*(LL(6,III)+QQ(6,III))
      END DO
C
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM=GRPSEQ(ITREE,ICHAIN,IGP-1)
C
C COMPUTE QDDOT
C
C TEMP5 = (M + K)
      DO J=1,6
         TEMP5(J,1)=MM(J,1,III)+KK(J,1,III)
         TEMP5(J,2)=MM(J,2,III)+KK(J,2,III)
         TEMP5(J,3)=MM(J,3,III)+KK(J,3,III)
         TEMP5(J,4)=MM(J,4,III)+KK(J,4,III)
         TEMP5(J,5)=MM(J,5,III)+KK(J,5,III)
         TEMP5(J,6)=MM(J,6,III)+KK(J,6,III)
      END DO
C
C TEMP1 = [B2^T * (M + K) * B2 ]^(-1)
      TEMP1=ZERO
      DO J=1,6
      TEMP1=TEMP1+( TEMP5(J,1)*B2(1,III) +
     &              TEMP5(J,2)*B2(2,III) +
     &              TEMP5(J,3)*B2(3,III) +
     &              TEMP5(J,4)*B2(4,III) +
     &              TEMP5(J,5)*B2(5,III) +
     &              TEMP5(J,6)*B2(6,III) ) *
     &               B2(J,III)
      END DO
      TEMP1=ONE/TEMP1
C
C TEMP = B_2^T * [ (M + K) * D - (Q + L) ]
      TEMP=ZERO
      DO J=1,6
      TEMP=TEMP+( ( TEMP5(J,1)*DD(1,III) +
     &              TEMP5(J,2)*DD(2,III) +
     &              TEMP5(J,3)*DD(3,III) +
     &              TEMP5(J,4)*DD(4,III) +
     &              TEMP5(J,5)*DD(5,III) +
     &              TEMP5(J,6)*DD(6,III) ) -
     &               QQ(J,III)-LL(J,III) ) * B2(J,III)
      END DO
C
C SCRATCH = B_2^T * (M + K) * B1
      DO J=1,6
      DO K=1,6
      TEMP3(J,K)=TEMP5(J,1)*B1(1,K,III) +
     &           TEMP5(J,2)*B1(2,K,III) +
     &           TEMP5(J,3)*B1(3,K,III) +
     &           TEMP5(J,4)*B1(4,K,III) +
     &           TEMP5(J,5)*B1(5,K,III) +
     &           TEMP5(J,6)*B1(6,K,III)
      END DO
      END DO
C
      DO J=1,6
      SCRATCH(J)=B2(1,III)*TEMP3(1,J) +
     &           B2(2,III)*TEMP3(2,J) +
     &           B2(3,III)*TEMP3(3,J) +
     &           B2(4,III)*TEMP3(4,J) +
     &           B2(5,III)*TEMP3(5,J) +
     &           B2(6,III)*TEMP3(6,J)
      END DO
C
      QDDOT(III)=-TEMP1*( SCRATCH(1)*YYDOT(1,IIIM) +
     &                    SCRATCH(2)*YYDOT(2,IIIM) +
     &                    SCRATCH(3)*YYDOT(3,IIIM) +
     &                    SCRATCH(4)*YYDOT(4,IIIM) +
     &                    SCRATCH(5)*YYDOT(5,IIIM) +
     &                    SCRATCH(6)*YYDOT(6,IIIM) + TEMP )
C
      DO J=1,6
      YYDOT(J,III)=B1(J,1,III)*YYDOT(1,IIIM) +
     &             B1(J,2,III)*YYDOT(2,IIIM) +
     &             B1(J,3,III)*YYDOT(3,IIIM) +
     &             B1(J,4,III)*YYDOT(4,IIIM) +
     &             B1(J,5,III)*YYDOT(5,IIIM) +
     &             B1(J,6,III)*YYDOT(6,IIIM) +
     &             B2(J,III)*QDDOT(III)+DD(J,III)
      END DO
C
      END DO
      END DO
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE TDYNKIN(NDEGF,XV,YV,ZV)
C
C routine computes kinetic energy and temperature of the system
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INTEGER NDEGF
      DOUBLE PRECISION XV(*), YV(*), ZV(*)
C local
      INTEGER I
      DOUBLE COMPLEX DCVAL
C begin
      RENR(SSTOTK)=0.0D0
      DO I=1,NATOM
      IF ( IMOVE(I).EQ.0 ) THEN
      RENR(SSTOTK)=RENR(SSTOTK)+AMASS(I)*(XV(I)**2+YV(I)**2+ZV(I)**2)
      END IF
      END DO
      RENR(SSTOTK)=RENR(SSTOTK)*0.5D0
      RENR(SSTOTE)=RENR(SSENER)+RENR(SSTOTK)
      RENR(SSTEMP)=2.0D0*RENR(SSTOTK)/(NDEGF*KBOLTZ)
      CALL DECLAR( ANER(SSTOTK), 'DP', ' ', DCVAL, RENR(SSTOTK) )
      CALL DECLAR( ANER(SSTOTE), 'DP', ' ', DCVAL, RENR(SSTOTE) )
      CALL DECLAR( ANER(SSTEMP), 'DP', ' ', DCVAL, RENR(SSTEMP) )
      RETURN
      END
C===============================================================
      SUBROUTINE TSCVEL(NCHAINS,CHNLEN,GRPSEQ,YY,QDOT)
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'deriv.inc'
C
C i/o
      INTEGER NCHAINS(MAXTREE),CHNLEN(MAXTREE,*)
      INTEGER GRPSEQ(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION YY(6,NGP),QDOT(*)
C
C
C local
C
      INTEGER IGP,III,J
      INTEGER ICHAIN,ITREE
      DOUBLE PRECISION SCVEL
C
C
C calculate the scaling term (scvel)
C so that velocity(t) * scvel = velocity(tbath)
C
C t (= SSTEMP in RENR array) is calculated at
C each step in the subroutine TDYNKIN
C
C
C
      SCVEL=SQRT(TBATH/RENR(SSTEMP))
C
C
C
C apply scale to body velocities
C
C base first
C
      DO ITREE=1,NTREE
        III=GRPSEQ(ITREE,1,1)
        DO J=1,6
          YY(J,III)=SCVEL*YY(J,III)
        END DO
C
C now all groups
C
        DO ICHAIN=1,NCHAINS(ITREE)
         DO IGP=2,CHNLEN(ITREE,ICHAIN)
           III=GRPSEQ(ITREE,ICHAIN,IGP)
           QDOT(III)=SCVEL*QDOT(III)
           DO J=1,6
             YY(J,III)=SCVEL*YY(J,III)
           END DO
         END DO
        END DO
      END DO
C
C also scale the atom velocites
C
      DO J=1,NATOM
        XV(J)=XV(J)*SCVEL
        YV(J)=YV(J)*SCVEL
        ZV(J)=ZV(J)*SCVEL
      END DO
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE Q2ROT(Q0,Q1,Q2,Q3,MAT)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION Q0, Q1, Q2, Q3, MAT(3,3)
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C begin
      MAT(1,1)=Q0*Q0+Q1*Q1-Q2*Q2-Q3*Q3
      MAT(2,1)=TWO*(Q1*Q2-Q0*Q3)
      MAT(3,1)=TWO*(Q1*Q3+Q0*Q2)
      MAT(1,2)=TWO*(Q1*Q2+Q0*Q3)
      MAT(2,2)=Q0*Q0-Q1*Q1+Q2*Q2-Q3*Q3
      MAT(3,2)=TWO*(Q2*Q3-Q0*Q1)
      MAT(1,3)=TWO*(Q1*Q3-Q0*Q2)
      MAT(2,3)=TWO*(Q2*Q3+Q0*Q1)
      MAT(3,3)=Q0*Q0-Q1*Q1-Q2*Q2+Q3*Q3
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE TQ2ROT(Q0,Q1,Q2,Q3,MAT)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION Q0, Q1, Q2, Q3, MAT(3,3)
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C begin
      MAT(1,1)=Q0*Q0+Q1*Q1-Q2*Q2-Q3*Q3
      MAT(1,2)=TWO*(Q1*Q2-Q0*Q3)
      MAT(1,3)=TWO*(Q1*Q3+Q0*Q2)
      MAT(2,1)=TWO*(Q1*Q2+Q0*Q3)
      MAT(2,2)=Q0*Q0-Q1*Q1+Q2*Q2-Q3*Q3
      MAT(2,3)=TWO*(Q2*Q3-Q0*Q1)
      MAT(3,1)=TWO*(Q1*Q3-Q0*Q2)
      MAT(3,2)=TWO*(Q2*Q3+Q0*Q1)
      MAT(3,3)=Q0*Q0-Q1*Q1-Q2*Q2+Q3*Q3
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE ROTMUL1(A,B,C)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
C begin
      C(1,1)=A(1,1)*B(1,1)+A(2,1)*B(2,1)+A(3,1)*B(3,1)
      C(2,1)=A(1,2)*B(1,1)+A(2,2)*B(2,1)+A(3,2)*B(3,1)
      C(3,1)=A(1,3)*B(1,1)+A(2,3)*B(2,1)+A(3,3)*B(3,1)
      C(1,2)=A(1,1)*B(1,2)+A(2,1)*B(2,2)+A(3,1)*B(3,2)
      C(2,2)=A(1,2)*B(1,2)+A(2,2)*B(2,2)+A(3,2)*B(3,2)
      C(3,2)=A(1,3)*B(1,2)+A(2,3)*B(2,2)+A(3,3)*B(3,2)
      C(1,3)=A(1,1)*B(1,3)+A(2,1)*B(2,3)+A(3,1)*B(3,3)
      C(2,3)=A(1,2)*B(1,3)+A(2,2)*B(2,3)+A(3,2)*B(3,3)
      C(3,3)=A(1,3)*B(1,3)+A(2,3)*B(2,3)+A(3,3)*B(3,3)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE ROTMUL2(A,B,C)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
C begin
      C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
      C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
      C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
      C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
      C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
      C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
      C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)
      C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
      C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE ROTMUL3(A,B,C)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION A(3,3), B(3,3), C(3,3)
C begin
      C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(1,2)+A(1,3)*B(1,3)
      C(1,2)=A(2,1)*B(1,1)+A(2,2)*B(1,2)+A(2,3)*B(1,3)
      C(1,3)=A(3,1)*B(1,1)+A(3,2)*B(1,2)+A(3,3)*B(1,3)
      C(2,1)=A(1,1)*B(2,1)+A(1,2)*B(2,2)+A(1,3)*B(2,3)
      C(2,2)=A(2,1)*B(2,1)+A(2,2)*B(2,2)+A(2,3)*B(2,3)
      C(2,3)=A(3,1)*B(2,1)+A(3,2)*B(2,2)+A(3,3)*B(2,3)
      C(3,1)=A(1,1)*B(3,1)+A(1,2)*B(3,2)+A(1,3)*B(3,3)
      C(3,2)=A(2,1)*B(3,1)+A(2,2)*B(3,2)+A(2,3)*B(3,3)
      C(3,3)=A(3,1)*B(3,1)+A(3,2)*B(3,2)+A(3,3)*B(3,3)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE AV2Q(Q0,Q1,Q2,Q3,YY,QVEL)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION Q0, Q1, Q2, Q3, YY(6)
      DOUBLE PRECISION QVEL(0:3)
C parameter
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
C begin
      QVEL(0)=HALF*(-Q1*YY(4)-Q2*YY(5)-Q3*YY(6))
      QVEL(1)=HALF*( Q0*YY(4)+Q3*YY(5)-Q2*YY(6))
      QVEL(2)=HALF*(-Q3*YY(4)+Q0*YY(5)+Q1*YY(6))
      QVEL(3)=HALF*( Q2*YY(4)-Q1*YY(5)+Q0*YY(6))
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE AA2Q(Q0,Q1,Q2,Q3,YY,YYDOT,QVEL,QACC)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION Q0, Q1, Q2, Q3, YY(6), YYDOT(6)
      DOUBLE PRECISION QVEL(0:3), QACC(0:3)
C parameter
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
C begin
      QACC(0)=HALF*(-QVEL(1)*YY(4)-Q1*YYDOT(4) -
     @               QVEL(2)*YY(5)-Q2*YYDOT(5) -
     @               QVEL(3)*YY(6)-Q3*YYDOT(6))
      QACC(1)=HALF*( QVEL(0)*YY(4)+Q0*YYDOT(4) +
     @               QVEL(3)*YY(5)+Q3*YYDOT(5) -
     @               QVEL(2)*YY(6)-Q2*YYDOT(6))
      QACC(2)=HALF*(-QVEL(3)*YY(4)-Q3*YYDOT(4) +
     @               QVEL(0)*YY(5)+Q0*YYDOT(5) +
     @               QVEL(1)*YY(6)+Q1*YYDOT(6))
      QACC(3)=HALF*( QVEL(2)*YY(4)+Q2*YYDOT(4) -
     @               QVEL(1)*YY(5)-Q1*YYDOT(5) +
     @               QVEL(0)*YY(6)+Q0*YYDOT(6))
C
      RETURN
      END
C
C =====================================================================
C
