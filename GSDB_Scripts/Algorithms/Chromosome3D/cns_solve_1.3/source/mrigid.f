C
C======================================================================
C
      SUBROUTINE MINIMI
C
C parses minimizer commands
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      CALL NEXTWD('MINImize-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-minimize')
C
      ELSE IF (WD(1:4).EQ.'POWE') THEN
      CALL POWELL
      ELSE IF (WD(1:4).EQ.'RIGI') THEN
      CALL MRIGID
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      ELSE IF (WD(1:5).EQ.'LBFGS') THEN
      CALL LBFGS
C=====================================================================
C #endif
C=====================================================================
      ELSE
      CALL DSPERR('MINIMIZE','unknown qualifier')
      END IF
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE MRIGID
C
C Parser routine for rigid body minimization.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'mrigid.inc'
C local
      INTEGER FLAGS
C begin
      MAXNGR=MAX(NATOM,1)
      IPOINT=ALLHP(INTEG4(MAXNGR+1))
      HPLIST=ALLHP(INTEG4(NATOM))
      FLAGS=ALLHP(INTEG4(NATOM))
      CALL MRIGI2(HEAP(FLAGS),HEAP(HPLIST),HEAP(IPOINT))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      CALL FREHP(HPLIST,INTEG4(NATOM))
      CALL FREHP(IPOINT,INTEG4(MAXNGR+1))
      RETURN
      END
C
      SUBROUTINE MRIGI2(FLAGS,LIST,POINTR)
C
C Parser routine for rigid body minimization.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mrigid.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER FLAGS(*), LIST(*), POINTR(*)
C local
      INTEGER NSTEP, I, II, IER, NFLAGS
      INTEGER ROTTRA, GRADNT, WORK, FINITE
      DOUBLE PRECISION STEP, TOLER, TARGET
      EXTERNAL RTFUNC
      LOGICAL DEFAUL, QDEBUG
C begin
C defaults
      QTRANS=.TRUE.
      QDEBUG=.FALSE.
      DEFAUL=.TRUE.
C make a default selection consisting of all free known atoms, the first
C GROUp statement will overwrite this default selection.
      NGROUP=1
      POINTR(1)=0
      POINTR(2)=0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0.AND.INITIA(I,X,Y,Z)) THEN
      POINTR(2)=POINTR(2)+1
      LIST(POINTR(2))=I
      END IF
      END DO
C
      NSTEP=50
      STEP=1.0D0
      TOLER=0.01D0
      NPRINT=1
      NCALLS=0
C
      CALL PUSEND('RIGID>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('RIGID>')
      CALL MISCOM('RIGID>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-minimize-rigid')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(A,I6,A,G12.4,A,G12.4)')
     &' NSTEp=',NSTEP,' DROP=',STEP,' TOLErance=',TOLER
      WRITE(6,'(A,I5)')
     & ' Number of rigid groups=',NGROUP
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('NSTEp=',NSTEP)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NPRI') THEN
      CALL NEXTI('NPRInt=',NPRINT)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'STEP'.OR.WD(1:4).EQ.'DROP') THEN
      CALL NEXTF('DROP=',STEP)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TOLE') THEN
      CALL NEXTF('TOLErance=',TOLER)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAN') THEN
      CALL NEXTLO('TRANslate=',QTRANS)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
C
      IF (DEFAUL) THEN
C
C erase default selection
      NGROUP=0
      DEFAUL=.FALSE.
      END IF
C
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
C
C make sure that the selections are disjoint.
      DO I=1,POINTR(NGROUP+1)
      IF (FLAGS(LIST(I)).EQ.1) THEN
      WRITE(6,'(2A,1X,A,1X,A,1X,A)')
     & ' MRIGID-ERR: selection overlap for atom ',
     & SEGID(LIST(I)),RESID(LIST(I)),RES(LIST(I)),TYPE(LIST(I))
      WRITE(6,'(A)')
     & ' The GROUp statements select at least one atom to be ',
     & ' in two or more groups.  Please check your GROUp ',
     & ' selections with respect to overlapping definitions.'
      CALL WRNDIE(-1,'MRIGID',
     & 'GROUp definition overlap --> change GROUp selections')
      FLAGS(LIST(I))=0
      END IF
      END DO
C
C exclude unknown atoms
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      FLAGS(I)=0
      END IF
      END DO
C
      CALL MAKIND(FLAGS,NATOM,NFLAGS)
C
      IF (NFLAGS.GT.0) THEN
C
C make new group
      IF (NGROUP.GE.MAXNGR) THEN
      CALL WRNDIE(-5,'MRIGID',
     & 'MAXNGR (max. no. of groups) exceeded')
      ELSE
      NGROUP=NGROUP+1
      END IF
C
C ok, now that we know that everything is ok we can copy the
C list of atom pointers
      II=POINTR(NGROUP)
      DO I=1,NFLAGS
      II=II+1
      LIST(II)=FLAGS(I)
      END DO
      POINTR(NGROUP+1)=II
      END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'DEBU') THEN
      QDEBUG=.TRUE.
      RPRINT=.FALSE.
      ROTTRA=ALLHP(IREAL8(NGROUP*6))
      GRADNT=ALLHP(IREAL8(NGROUP*6))
      FINITE=ALLHP(IREAL8(NGROUP*6))
      CALL RGDDBG(HEAP(ROTTRA),HEAP(GRADNT),HEAP(FINITE),NGROUP,
     &            POINTR,LIST,QTRANS)
      CALL FREHP(FINITE,IREAL8(NGROUP*6))
      CALL FREHP(GRADNT,IREAL8(NGROUP*6))
      CALL FREHP(ROTTRA,IREAL8(NGROUP*6))
C=====================================================================
      ELSE
      CALL CHKEND('OPTIMIZE>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (POINTR(NGROUP+1).GT.0) THEN
C
      IF (.NOT.QDEBUG) THEN
      TOLER=(TOLER**2)*6
      RPRINT=.TRUE.
C
C allocate space for the parameters, the gradient and the work array
      ROTTRA=ALLHP(IREAL8(NGROUP*6))
      GRADNT=ALLHP(IREAL8(NGROUP*6))
      WORK=ALLHP(IREAL8(NGROUP*36))
      CALL RGDINI(HEAP(ROTTRA),NGROUP)
      IF (QTRANS) THEN
      CALL ZXCGR(RTFUNC,NGROUP*6,TOLER,NSTEP,STEP,HEAP(ROTTRA),
     &           HEAP(GRADNT),TARGET,HEAP(WORK),IER)
      ELSE
      CALL ZXCGR(RTFUNC,NGROUP*6-3,TOLER,NSTEP,STEP,HEAP(ROTTRA),
     &           HEAP(GRADNT),TARGET,HEAP(WORK),IER)
      END IF
      IF (IER.EQ.0) THEN
      WRITE(6,'(A)') ' ZXCGR: gradient converged'
      ELSE IF (IER.EQ.129) THEN
      WRITE(6,'(A)') ' ZXCGR: Line search terminated'
      ELSE IF (IER.EQ.130) THEN
      WRITE(6,'(A)') ' %ZXCGR-ERR: Search direction uphill'
      ELSE IF (IER.EQ.131) THEN
      WRITE(6,'(A)') ' ZXCGR: NSTEP limit reached'
      ELSE IF (IER.EQ.132) THEN
      WRITE(6,'(A)') ' %ZXCGR-ERR: failure to reduce E'
      END IF
C
      IF (IER.EQ.0.OR.IER.EQ.131.OR.IER.EQ.129.OR.IER.EQ.131) THEN
C
C apply the rotation/translation to the main coordinates
      CALL RTFUNC(0,HEAP(ROTTRA),TARGET,HEAP(GRADNT))
      WRITE(6,'(A)') ' RIGID: main coordinates set to best minimum'
      END IF
      CALL FREHP(WORK,IREAL8(NGROUP*36))
      CALL FREHP(GRADNT,IREAL8(NGROUP*6))
      CALL FREHP(ROTTRA,IREAL8(NGROUP*6))
      END IF
C
      END IF
C
      RETURN
      END
C
      SUBROUTINE RGDINI(ROTTRA,NGROUP)
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION ROTTRA(*)
      INTEGER NGROUP
C local
      INTEGER I
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      DO I=1,NGROUP*6
      ROTTRA(I)=ZERO
      END DO
      RETURN
      END
C
      SUBROUTINE RGDDBG(ROTTRA,GRADNT,FINITE,NGROUP,POINTR,LIST,QTRANS)
C
C Tests first derivatives of grouped rigid body minimizer
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C I/O
      DOUBLE PRECISION ROTTRA(*), GRADNT(*), FINITE(*)
      INTEGER NGROUP, POINTR(*), LIST(*)
      LOGICAL QTRANS
C local
      DOUBLE PRECISION DELTA, TARGET
      INTEGER IG, OFFSET, I, II
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0)
C begin
      WRITE(6,'(A)') ' List of group selections'
      DO I=1,NGROUP
      WRITE(6,'(A,I5,A)') ' ------ GROUP ',I,' -----------'
      DO II=POINTR(I)+1,POINTR(I+1)
      WRITE(6,'(1X,A,1X,A,1X,A,1X,A)')
     &    SEGID(LIST(II)),RESID(LIST(II)),RES(LIST(II)),TYPE(LIST(II))
      END DO
      END DO
      WRITE(6,'(A)') ' Performing finite difference test'
C
      DO IG=1,NGROUP
      OFFSET=(IG-1)*6
      CALL NEXTF('thetaz=',ROTTRA(1+OFFSET))
      CALL NEXTF('thetay=',ROTTRA(2+OFFSET))
      CALL NEXTF('thetax=',ROTTRA(3+OFFSET))
      IF (QTRANS.OR.IG.LT.NGROUP) THEN
      CALL NEXTF('transx=',ROTTRA(4+OFFSET))
      CALL NEXTF('transy=',ROTTRA(5+OFFSET))
      CALL NEXTF('transz=',ROTTRA(6+OFFSET))
      ELSE
      ROTTRA(4+OFFSET)=0
      ROTTRA(5+OFFSET)=0
      ROTTRA(6+OFFSET)=0
      END IF
      END DO
      CALL NEXTF('delta=',DELTA)
C
      DO I=1,6*NGROUP
      ROTTRA(I)=ROTTRA(I)+DELTA
      CALL RTFUNC(NGROUP*6,ROTTRA,TARGET,GRADNT)
      FINITE(I)=TARGET
      ROTTRA(I)=ROTTRA(I)-TWO*DELTA
      CALL RTFUNC(NGROUP*6,ROTTRA,TARGET,GRADNT)
      ROTTRA(I)=ROTTRA(I)+DELTA
      FINITE(I)=FINITE(I)-TARGET
      FINITE(I)=FINITE(I)/(TWO*DELTA)
      END DO
      DO IG=1,NGROUP
      OFFSET=(IG-1)*6
      WRITE(6,'(A,I6,/,A,6F11.4)')
     &' Group ',IG,
     & ' fin.dif.=',(FINITE(I+OFFSET),I=1,6)
      CALL RTFUNC(NGROUP*6,ROTTRA,TARGET,GRADNT)
      WRITE(6,'(A,6F11.4)')
     & ' analyt. =',(GRADNT(I+OFFSET),I=1,6)
      END DO
      RETURN
      END
C
      SUBROUTINE RTFUNC(DIM,ROTTRA,TARGET,GRADNT)
C
C Target function for conjugent gradient minimizer.
C
C Subroutine applies rotation, translation according to specified
C ROTTRA (thetaz, thetay, thetax, trans-x, trans-y, trans-z) to
C current coordinates, computes energy, and derivatives with regard
C to (thetaz, thetay, thetaz, trans-x, trans-y, trans-z).  Energy
C is returned in TARGET, gradient is returned in GRADNT.
C
C NOTE: DIM is used as a flag to indicate whether to store the
C rotated/translated coordinates into main coordinate set.  Normally,
C we don't want to store them because they are recomputed each time.
C The minimizer is calling RTFUNC with DIM > 0.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mrigid.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INTEGER DIM
      DOUBLE PRECISION ROTTRA(*), TARGET, GRADNT(*)
C local
      INTEGER XCENT, YCENT, ZCENT
C begin
      XCENT=ALLHP(IREAL8(NGROUP))
      YCENT=ALLHP(IREAL8(NGROUP))
      ZCENT=ALLHP(IREAL8(NGROUP))
      CALL RTFUN2(DIM,ROTTRA,TARGET,GRADNT,HEAP(HPLIST),
     &           HEAP(IPOINT),NGROUP,HEAP(XCENT),HEAP(YCENT),
     &           HEAP(ZCENT),RPRINT,NCALLS,NPRINT,QTRANS)
      CALL FREHP(ZCENT,IREAL8(NGROUP))
      CALL FREHP(YCENT,IREAL8(NGROUP))
      CALL FREHP(XCENT,IREAL8(NGROUP))
      RETURN
      END
      SUBROUTINE RTFUN2(DIM,ROTTRA,TARGET,GRADNT,LIST,POINTR,
     &           NGROUP,XCENT,YCENT,ZCENT,RPRINT,NCALLS,NPRINT,QTRANS)
C
C Target function for conjugent gradient minimizer.  See above
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER DIM
      DOUBLE PRECISION ROTTRA(*), TARGET, GRADNT(*)
      INTEGER LIST(*), POINTR(*), NGROUP
      DOUBLE PRECISION XCENT(*), YCENT(*), ZCENT(*)
      LOGICAL RPRINT
      INTEGER NCALLS, NPRINT
      LOGICAL QTRANS
C local
      DOUBLE PRECISION XNEW, YNEW, ZNEW
      DOUBLE PRECISION S1, S2, S3, C1, C2, C3, R(3,3)
      INTEGER I, NCENT, OFFSET
      INTEGER IX, IY, IZ, IG
C parameters
      DOUBLE PRECISION RAD, ZERO, TWO
      PARAMETER (RAD=PI/180.0D0, ZERO=0.0D0, TWO=2.0D0)
C begin
C
C copy current cartesian coordinates
      IF (DIM.GT.0) THEN
      IX=ALLHP(IREAL8(NATOM))
      IY=ALLHP(IREAL8(NATOM))
      IZ=ALLHP(IREAL8(NATOM))
      CALL COPYCO(NATOM,X,Y,Z,HEAP(IX),HEAP(IY),HEAP(IZ))
      END IF
C
C loop over all groups
      DO IG=1,NGROUP
      OFFSET=(IG-1)*6
C
C compute rotation matrix
      S1=SIN(RAD*ROTTRA(1+OFFSET))
      S2=SIN(RAD*ROTTRA(2+OFFSET))
      S3=SIN(RAD*ROTTRA(3+OFFSET))
      C1=COS(RAD*ROTTRA(1+OFFSET))
      C2=COS(RAD*ROTTRA(2+OFFSET))
      C3=COS(RAD*ROTTRA(3+OFFSET))
C
      R(1,1)=C1*C2
      R(1,2)=S1*C3-C1*S2*S3
      R(1,3)=S1*S3+C1*S2*C3
      R(2,1)=-S1*C2
      R(2,2)=C1*C3+S1*S2*S3
      R(2,3)=C1*S3-S1*S2*C3
      R(3,1)=-S2
      R(3,2)=-C2*S3
      R(3,3)=C2*C3
C
C determine geometric center of atoms of group
      XCENT(IG)=ZERO
      YCENT(IG)=ZERO
      ZCENT(IG)=ZERO
      NCENT=0
      DO I=POINTR(IG)+1,POINTR(IG+1)
      XCENT(IG)=XCENT(IG)+X(LIST(I))
      YCENT(IG)=YCENT(IG)+Y(LIST(I))
      ZCENT(IG)=ZCENT(IG)+Z(LIST(I))
      NCENT=NCENT+1
      END DO
      XCENT(IG)=XCENT(IG)/NCENT
      YCENT(IG)=YCENT(IG)/NCENT
      ZCENT(IG)=ZCENT(IG)/NCENT
C
C transform coordinates of group
      DO I=POINTR(IG)+1,POINTR(IG+1)
      XNEW=X(LIST(I))-XCENT(IG)
      YNEW=Y(LIST(I))-YCENT(IG)
      ZNEW=Z(LIST(I))-ZCENT(IG)
      X(LIST(I))=R(1,1)*XNEW+R(1,2)*YNEW+R(1,3)*ZNEW +XCENT(IG)
      Y(LIST(I))=R(2,1)*XNEW+R(2,2)*YNEW+R(2,3)*ZNEW +YCENT(IG)
      Z(LIST(I))=R(3,1)*XNEW+R(3,2)*YNEW+R(3,3)*ZNEW +ZCENT(IG)
C
C add the translation for all groups except the last when QTRANS is
C FALSE
      IF (QTRANS.OR.IG.LT.NGROUP) THEN
      X(LIST(I))=X(LIST(I))+ROTTRA(4+OFFSET)
      Y(LIST(I))=Y(LIST(I))+ROTTRA(5+OFFSET)
      Z(LIST(I))=Z(LIST(I))+ROTTRA(6+OFFSET)
      END IF
      END DO
C
      END DO
C
C compute energy
      CALL ENERGY
C
C target is simply the total energy
      TARGET=RENR(SSENER)
C
C copy original coordinates back
      IF (DIM.GT.0) THEN
      CALL COPYCO(NATOM,HEAP(IX),HEAP(IY),HEAP(IZ),X,Y,Z)
      CALL FREHP(IX,IREAL8(NATOM))
      CALL FREHP(IY,IREAL8(NATOM))
      CALL FREHP(IZ,IREAL8(NATOM))
      END IF
C
C loop over all groups
      DO IG=1,NGROUP
      OFFSET=(IG-1)*6
C
C initialize gradient
      DO I=1,6
      GRADNT(I+OFFSET)=ZERO
      END DO
C
C gradient for the X, Y, Z translations
C compute the gradient for all groups except the last when QTRANS is
C FALSE
      IF (QTRANS.OR.IG.LT.NGROUP) THEN
      DO I=POINTR(IG)+1,POINTR(IG+1)
      GRADNT(4+OFFSET)=GRADNT(4+OFFSET)+DX(LIST(I))
      GRADNT(5+OFFSET)=GRADNT(5+OFFSET)+DY(LIST(I))
      GRADNT(6+OFFSET)=GRADNT(6+OFFSET)+DZ(LIST(I))
      END DO
      END IF
C
      S1=SIN(RAD*ROTTRA(1+OFFSET))
      S2=SIN(RAD*ROTTRA(2+OFFSET))
      S3=SIN(RAD*ROTTRA(3+OFFSET))
      C1=COS(RAD*ROTTRA(1+OFFSET))
      C2=COS(RAD*ROTTRA(2+OFFSET))
      C3=COS(RAD*ROTTRA(3+OFFSET))
C
C gradient for thetaz
      R(1,1)=-S1*C2
      R(1,2)=C1*C3+S1*S2*S3
      R(1,3)=C1*S3-S1*S2*C3
      R(2,1)=-C1*C2
      R(2,2)=-S1*C3+C1*S2*S3
      R(2,3)=-S1*S3-C1*S2*C3
      R(3,1)=ZERO
      R(3,2)=ZERO
      R(3,3)=ZERO
      DO I=POINTR(IG)+1,POINTR(IG+1)
      GRADNT(1+OFFSET)=GRADNT(1+OFFSET)
     & +DX(LIST(I))*
     &   (R(1,1)*(X(LIST(I))-XCENT(IG))+R(1,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(1,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DY(LIST(I))*
     &   (R(2,1)*(X(LIST(I))-XCENT(IG))+R(2,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(2,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DZ(LIST(I))*
     &   (R(3,1)*(X(LIST(I))-XCENT(IG))+R(3,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(3,3)*(Z(LIST(I))-ZCENT(IG)))
      END DO
      GRADNT(1+OFFSET)=GRADNT(1+OFFSET)*RAD
C
C gradient for thetay
      R(1,1)=-C1*S2
      R(1,2)=-C1*C2*S3
      R(1,3)=+C1*C2*C3
      R(2,1)=+S1*S2
      R(2,2)=S1*C2*S3
      R(2,3)=-S1*C2*C3
      R(3,1)=-C2
      R(3,2)=+S2*S3
      R(3,3)=-S2*C3
      DO I=POINTR(IG)+1,POINTR(IG+1)
      GRADNT(2+OFFSET)=GRADNT(2+OFFSET)
     & +DX(LIST(I))*
     &   (R(1,1)*(X(LIST(I))-XCENT(IG))+R(1,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(1,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DY(LIST(I))*
     &   (R(2,1)*(X(LIST(I))-XCENT(IG))+R(2,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(2,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DZ(LIST(I))*
     &   (R(3,1)*(X(LIST(I))-XCENT(IG))+R(3,2)
     &  *(Y(LIST(I))-YCENT(IG))+R(3,3)*(Z(LIST(I))-ZCENT(IG)))
      END DO
      GRADNT(2+OFFSET)=GRADNT(2+OFFSET)*RAD
C
C gradient for thetax
      R(1,1)=ZERO
      R(1,2)=-S1*S3-C1*S2*C3
      R(1,3)=S1*C3-C1*S2*S3
      R(2,1)=ZERO
      R(2,2)=-C1*S3+S1*S2*C3
      R(2,3)=C1*C3+S1*S2*S3
      R(3,1)=ZERO
      R(3,2)=-C2*C3
      R(3,3)=-C2*S3
      DO I=POINTR(IG)+1,POINTR(IG+1)
      GRADNT(3+OFFSET)=GRADNT(3+OFFSET)
     & +DX(LIST(I))*
     &   (R(1,1)*(X(LIST(I))-XCENT(IG))+R(1,2)
     &   *(Y(LIST(I))-YCENT(IG))+R(1,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DY(LIST(I))*
     &   (R(2,1)*(X(LIST(I))-XCENT(IG))+R(2,2)
     &   *(Y(LIST(I))-YCENT(IG))+R(2,3)*(Z(LIST(I))-ZCENT(IG)))
     & +DZ(LIST(I))*
     &   (R(3,1)*(X(LIST(I))-XCENT(IG))+R(3,2)
     &   *(Y(LIST(I))-YCENT(IG))+R(3,3)*(Z(LIST(I))-ZCENT(IG)))
      END DO
      GRADNT(3+OFFSET)=GRADNT(3+OFFSET)*RAD
      END DO
C
      IF (RPRINT) THEN
      NCALLS=NCALLS+1
      IF (MOD(NCALLS,NPRINT).EQ.0) THEN
      WRITE(6,'(A,A)')
          WRITE(6,'(A,I6,A)')
     & ' --------------- cycle=',NCALLS,
     & ' --------------------------------------------------'
      IF (WRNLEV.GE.5) THEN
      DO IG=1,NGROUP
      OFFSET=(IG-1)*6
      WRITE(6,'(A,I5,A,6F8.2,A)')
     & ' | group=',IG,
     &' rot/tran=(',ROTTRA(1+OFFSET),ROTTRA(2+OFFSET),
     &   ROTTRA(3+OFFSET),ROTTRA(4+OFFSET),ROTTRA(5+OFFSET),
     &   ROTTRA(6+OFFSET),')     |'
      END DO
      END IF
      CALL PRINTE(RENR)
      WRITE(6,'(2A)') ' --------------------------------------------',
     & '-----------------------------------'
      END IF
      END IF
C
      RETURN
      END
C
C   IMSL ROUTINE NAME   - ZXCGR
C
C-----------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   LATEST REVISION     - NOVEMBER 1, 1979
C
C   PURPOSE             - A CONJUGATE GRADIENT ALGORITHM FOR FINDING
C                           THE MINIMUM OF A FUNCTION OF N VARIABLES
C
C   USAGE               - CALL ZXCGR (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,
C                           IER)
C
C   ARGUMENTS    FUNCT  - A USER SUPPLIED SUBROUTINE WHICH CALCULATES
C                           THE OBJECTIVE FUNCTION AND ITS GRADIENT
C                           FOR GIVEN PARAMETER VALUES
C                           X(1),X(2),...,X(N).
C                           THE CALLING SEQUENCE HAS THE FOLLOWING FORM
C                           CALL FUNCT (N,X,F,G)
C                           WHERE X AND G ARE VECTORS OF LENGTH N.
C                           THE SCALAR F IS FOR THE OBJECTIVE FUNCTION.
C                           G(1), G(2), ..., G(N) ARE FOR THE COMPONENTS
C                           OF THE GRADIENT OF F.
C                           FUNCT MUST APPEAR IN AN EXTERNAL STATEMENT
C                           IN THE CALLING PROGRAM. FUNCT MUST NOT
C                           ALTER THE VALUES OF X(I),I=1,...,N OR N.
C                N      - THE NUMBER OF PARAMETERS OF THE OBJECTIVE
C                           FUNCTION. (INPUT) (I.E.,THE LENGTH OF X)
C                ACC    - CONVERGENCE CRITERION. (INPUT)
C                           THE CALCULATION ENDS WHEN THE SUM OF SQUARES
C                           OF THE COMPONENTS OF G IS LESS THAN ACC.
C                MAXFN  - MAXIMUM NUMBER OF FUNCTION EVALUATIONS (I.E.,
C                           CALLS TO SUBROUTINE FUNCT) ALLOWED. (INPUT)
C                DFPRED - A ROUGH ESTIMATE OF THE EXPECTED REDUCTION
C                           IN F, WHICH IS USED TO DETERMINE THE SIZE
C                           OF THE INITIAL CHANGE TO X. (INPUT)
C                           NOTE THAT DFPRED IS THE EXPECTED REDUCTION
C                           ITSELF, AND DOES NOT DEPEND ON ANY RATIOS.
C                           A BAD VALUE OF DFPRED CAUSES AN ERROR
C                           MESSAGE, WITH IER=129, AND A RETURN ON THE
C                           FIRST ITERATION. (SEE THE DESCRIPTION OF
C                           IER BELOW)
C                X      - VECTOR OF LENGTH N CONTAINING PARAMETER
C                           VALUES.
C                         ON INPUT, X MUST CONTAIN THE INITIAL
C                           PARAMETER ESTIMATES.
C                         ON OUTPUT, X CONTAINS THE FINAL PARAMETER
C                           ESTIMATES AS DETERMINED BY ZXCGR.
C                G      - A VECTOR OF LENGTH N CONTAINING THE
C                           COMPONENTS OF THE GRADIENT OF F AT THE
C                           FINAL PARAMETER ESTIMATES. (OUTPUT)
C                F      - A SCALAR CONTAINING THE VALUE OF THE FUNCTION
C                           AT THE FINAL PARAMETER ESTIMATES. (OUTPUT)
C                W      - WORK VECTOR OF LENGTH 6*N.
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         IER = 0 IMPLIES THAT CONVERGENCE WAS
C                           ACHIEVED AND NO ERRORS OCCURRED.
C                         TERMINAL ERROR
C                           IER = 129 IMPLIES THAT THE LINE SEARCH OF
C                             AN INTEGRATION WAS ABANDONED. THIS
C                             ERROR MAY BE CAUSED BY AN ERROR IN THE
C                             GRADIENT.
C                           IER = 130 IMPLIES THAT THE CALCULATION
C                             CANNOT CONTINUE BECAUSE THE SEARCH
C                             DIRECTION IS UPHILL.
C                           IER = 131 IMPLIES THAT THE ITERATION WAS
C                             TERMINATED BECAUSE MAXFN WAS EXCEEDED.
C                           IER = 132 IMPLIES THAT THE CALCULATION
C                             WAS TERMINATED BECAUSE TWO CONSECUTIVE
C                             ITERATIONS FAILED TO REDUCE F.
C
C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
C                       - SINGLE/H36,H48,H60
C
C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
C
C   REMARKS  1.  THE ROUTINE INCLUDES NO THOROUGH CHECKS ON THE PART
C                OF THE USER PROGRAM THAT CALCULATES THE DERIVATIVES
C                OF THE OBJECTIVE FUNCTION. THEREFORE, BECAUSE
C                DERIVATIVE CALCULATION IS A FREQUENT SOURCE OF
C                ERROR, THE USER SHOULD VERIFY INDEPENDENTLY THE
C                CORRECTNESS OF THE DERIVATIVES THAT ARE GIVEN TO
C                THE ROUTINE.
C            2.  BECAUSE OF THE CLOSE RELATION BETWEEN THE CONJUGATE
C                GRADIENT METHOD AND THE METHOD OF STEEPEST DESCENTS,
C                IT IS VERY HELPFUL TO CHOOSE THE SCALE OF THE
C                VARIABLES IN A WAY THAT BALANCES THE MAGNITUDES OF
C                THE COMPONENTS OF A TYPICAL DERIVATE VECTOR. IT
C                CAN BE PARTICULARLY INEFFICIENT IF A FEW COMPONENTS
C                OF THE GRADIENT ARE MUCH LARGER THAN THE REST.
C            3.  IF THE VALUE OF THE PARAMETER ACC IN THE ARGUMENT
C                LIST OF THE ROUTINE IS SET TO ZERO, THEN THE
C                SUBROUTINE WILL CONTINUE ITS CALCULATION UNTIL IT
C                STOPS REDUCING THE OBJECTIVE FUNCTION. IN THIS CASE
C                THE USUAL BEHAVIOUR IS THAT CHANGES IN THE
C                OBJECTIVE FUNCTION BECOME DOMINATED BY COMPUTER
C                ROUNDING ERRORS BEFORE PRECISION IS LOST IN THE
C                GRADIENT VECTOR. THEREFORE, BECAUSE THE POINT OF
C                VIEW HAS BEEN TAKEN THAT THE USER REQUIRES THE
C                LEAST POSSIBLE VALUE OF THE FUNCTION, A VALUE OF
C                THE OBJECTIVE FUNCTION THAT IS SMALL DUE TO
C                COMPUTER ROUNDING ERRORS CAN PREVENT FURTHER
C                PROGRESS. HENCE THE PRECISION IN THE FINAL VALUES
C                OF THE VARIABLES MAY BE ONLY ABOUT HALF THE
C                NUMBER OF SIGNIFICANT DIGITS IN THE COMPUTER
C                ARITHMETIC, BUT THE LEAST VALUE OF F IS USUALLY
C                FOUND TO QUITE HIGH ACCURACY.
C
C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.
C
C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
C                           EXPRESSED OR IMPLIED, IS APPLICABLE.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE ZXCGR  (FUNCT,N,ACC,MAXFN,DFPRED,X,G,F,W,IER)
C                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N,MAXFN,IER
      DOUBLE PRECISION   ACC,DFPRED,X(N),G(N),F,W(*)
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      EXTERNAL FUNCT
      INTEGER            MAXLIN,MXFCON,I,IGINIT,IGOPT,IRETRY,IRSDG
      INTEGER            IRSDX,ITERC,ITERFM,ITERRS,IXOPT,NCALLS,NFBEG
      INTEGER            NFOPT
      DOUBLE PRECISION   BETA,DDSPLN,DFPR,FCH,FINIT,FMIN,GAMDEN,GAMA
      DOUBLE PRECISION   GINIT,GMIN,GNEW,GSPLN,GSQRD,SBOUND,STEP,STEPCH
      DOUBLE PRECISION   STMIN,SUM,WORK
      PARAMETER (MAXLIN=5,MXFCON=2)
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
C                                  THE WORKING SPACE ARRAY IS SPLIT
C                                    INTO SIX VECTORS OF LENGTH N. THE
C                                    FIRST PART IS USED FOR THE SEARCH
C                                    DIRECTION OF AN ITERATION. THE
C                                    SECOND AND THIRD PARTS CONTAIN THE
C                                    INFORMATION THAT IS REQUIRED BY
C                                    THE CONJUGACY CONDITIONS OF THE
C                                    RESTART PROCEDURE. THE FOURTH PART
C                                    CONTAINS THE GRADIENT AT THE START
C                                    OF AN ITERATION. THE FIFTH PART
C                                    CONTAINS THE PARAMETERS THAT GIVE
C                                    THE LEAST CALCULATED VALUE OF F.
C                                    THE SIXTH PART CONTAINS THE
C                                    GRADIENT VECTOR WHERE F IS LEAST.
      IRSDX = N
      IRSDG = IRSDX+N
      IGINIT = IRSDG+N
      IXOPT = IGINIT+N
      IGOPT = IXOPT+N
C                                  SET SOME PARAMETERS TO BEGIN THE
C                                    CALCULATION. ITERC AND
C                                    NCALLS COUNT THE NUMBER OF
C                                    ITERATIONS AND CALLS OF FUNCT.
C                                    ITERFM IS THE NUMBER OF THE MOST
C                                    RECENT ITERATION THAT DECREASES F.
      ITERC = 0
      NCALLS = 0
      ITERFM = ITERC
C                                  CALL SUBROUTINE FUNCT. LET THE
C                                    INITIAL SEARCH DIRECTION BE MINUS
C                                    THE GRADIENT VECTOR. USUALLY THE
C                                    PARAMETER ITERRS GIVES THE
C                                    ITERATION NUMBER OF THE MOST
C                                    RECENT RESTART, BUT IT IS SET TO
C                                    ZERO WHEN THE STEEPEST DESCENT
C                                    DIRECTION IS USED.
    5 NCALLS = NCALLS+1
      CALL FUNCT (N,X,F,G)
      IF (NCALLS.GE.2) GO TO 20
   10 DO 15 I=1,N
   15 W(I) = -G(I)
      ITERRS = 0
      IF (ITERC.GT.0) GO TO 80
C                                  SET SUM TO G SQUARED. GMIN AND GNEW
C                                    ARE THE OLD AND THE NEW
C                                    DIRECTIONAL DERIVATIVES ALONG THE
C                                    CURRENT SEARCH DIRECTION. LET FCH
C                                    BE THE DIFFERENCE BETWEEN F AND
C                                    THE PREVIOUS BEST VALUE OF THE
C                                    OBJECTIVE FUNCTION.
   20 GNEW = 0.0D0
      SUM = 0.0D0
      DO 25 I=1,N
         GNEW = GNEW+W(I)*G(I)
   25 SUM = SUM+G(I)**2
      IF (NCALLS.EQ.1) GO TO 35
      FCH = F-FMIN
C                                  STORE THE VALUES OF X, F AND G, IF
C                                    THEY ARE THE BEST THAT HAVE BEEN
C                                    CALCULATED SO FAR, AND NOTE G
C                                    SQUARED AND THE VALUE OF NCALLS.
C                                    TEST FOR CONVERGENCE.
      IF (FCH) 35,30,50
   30 IF (GNEW/GMIN.LT.-1.0D0) GO TO 45
   35 FMIN = F
      GSQRD = SUM
      NFOPT = NCALLS
      DO 40 I=1,N
         W(IXOPT+I) = X(I)
   40 W(IGOPT+I) = G(I)
   45 IF (SUM.LE.ACC) GO TO 9005
C                                  TEST IF THE VALUE OF MAXFN ALLOWS
C                                    ANOTHER CALL OF FUNCT.
CCC   50 IF (NCALLS.NE.MAXFN) GO TO 55
C modification ATB 11/14/96
   50 IF (NCALLS.LT.MAXFN) GO TO 55
      IER = 131
      GO TO 9000
   55 IF (NCALLS.GT.1) GO TO 100
C                                  SET DFPR TO THE ESTIMATE OF THE
C                                    REDUCTION IN F GIVEN IN THE
C                                    ARGUMENT LIST, IN ORDER THAT THE
C                                    INITIAL CHANGE TO THE PARAMETERS
C                                    IS OF A SUITABLE SIZE. THE VALUE
C                                    OF STMIN IS USUALLY THE
C                                    STEP-LENGTH OF THE MOST RECENT
C                                    LINE SEARCH THAT GIVES THE LEAST
C                                    CALCULATED VALUE OF F.
      DFPR = DFPRED
      STMIN = DFPRED/GSQRD
C                                  BEGIN THE ITERATION
   80 ITERC = ITERC+1
C                                  STORE THE INITIAL FUNCTION VALUE AND
C                                    GRADIENT, CALCULATE THE INITIAL
C                                    DIRECTIONAL DERIVATIVE, AND BRANCH
C                                    IF ITS VALUE IS NOT NEGATIVE. SET
C                                    SBOUND TO MINUS ONE TO INDICATE
C                                    THAT A BOUND ON THE STEP IS NOT
C                                    KNOWN YET, AND SET NFBEG TO THE
C                                    CURRENT VALUE OF NCALLS. THE
C                                    PARAMETER IRETRY SHOWS THE NUMBER
C                                    OF ATTEMPTS AT SATISFYING THE BETA
C                                    CONDITION.
      FINIT = F
      GINIT = 0.0D0
      DO 85 I=1,N
         W(IGINIT+I) = G(I)
   85 GINIT = GINIT+W(I)*G(I)
      IF (GINIT.GE.0.0D0) GO TO 165
      GMIN = GINIT
      SBOUND = -1.0D0
      NFBEG = NCALLS
      IRETRY = -1
C                                  SET STEPCH SO THAT THE INITIAL
C                                    STEP-LENGTH IS CONSISTENT WITH THE
C                                    PREDICTED REDUCTION IN F, SUBJECT
C                                    TO THE CONDITION THAT IT DOES NOT
C                                    EXCEED THE STEP-LENGTH OF THE
C                                    PREVIOUS ITERATION. LET STMIN BE
C                                    THE STEP TO THE LEAST CALCULATED
C                                    VALUE OF F.
      STEPCH = MIN(STMIN,ABS(DFPR/GINIT))
      STMIN = 0.0D0
C                                  CALL SUBROUTINE FUNCT AT THE VALUE
C                                    OF X THAT IS DEFINED BY THE NEW
C                                    CHANGE TO THE STEP-LENGTH, AND LET
C                                    THE NEW STEP-LENGTH BE STEP. THE
C                                    VARIABLE WORK IS USED AS WORK
C                                    SPACE.
   90 STEP = STMIN+STEPCH
      WORK = 0.0D0
      DO 95 I=1,N
         X(I) = W(IXOPT+I)+STEPCH*W(I)
   95 WORK = MAX(WORK,ABS(X(I)-W(IXOPT+I)))
      IF (WORK.GT.0.0D0) GO TO 5
C                                  TERMINATE THE LINE SEARCH IF STEPCH
C                                    IS EFFECTIVELY ZERO.
      IF (NCALLS.GT.NFBEG+1) GO TO 115
      IF (ABS(GMIN/GINIT)-0.2D0) 170,170,115
C                                  LET SPLN BE THE QUADRATIC SPLINE
C                                    THAT INTERPOLATES THE CALCULATED
C                                    FUNCTION VALUES AND DIRECTIONAL
C                                    DERIVATIVES AT THE POINTS STMIN
C                                    AND STEP OF THE LINE SEARCH, WHERE
C                                    THE KNOT OF THE SPLINE IS AT
C                                    0.5*(STMIN+STEP). REVISE STMIN,
C                                    GMIN AND SBOUND, AND SET DDSPLN TO
C                                    THE SECOND DERIVATIVE OF SPLN AT
C                                    THE NEW STMIN. HOWEVER, IF FCH IS
C                                    ZERO, IT IS ASSUMED THAT THE
C                                    MAXIMUM ACCURACY IS ALMOST
C                                    ACHIEVED, SO DDSPLN IS CALCULATED
C                                    USING ONLY THE CHANGE IN THE
C                                    GRADIENT.
  100 WORK = (FCH+FCH)/STEPCH-GNEW-GMIN
      DDSPLN = (GNEW-GMIN)/STEPCH
      IF (NCALLS.GT.NFOPT) SBOUND = STEP
      IF (NCALLS.GT.NFOPT) GO TO 105
      IF (GMIN*GNEW.LE.0.0D0) SBOUND = STMIN
      STMIN = STEP
      GMIN = GNEW
      STEPCH = -STEPCH
  105 IF (FCH.NE.0.0D0) DDSPLN = DDSPLN+(WORK+WORK)/STEPCH
C
C                                  TEST FOR CONVERGENCE OF THE LINE
C                                    SEARCH, BUT FORCE AT LEAST TWO
C                                    STEPS TO BE TAKEN IN ORDER NOT TO
C                                    LOSE QUADRATIC TERMINATION.
      IF (GMIN.EQ.0.0D0) GO TO 170
      IF (NCALLS.LE.NFBEG+1) GO TO 120
      IF (ABS(GMIN/GINIT).LE.0.2D0) GO TO 170
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MAXLIN.
  110 IF (NCALLS.LT.NFOPT+MAXLIN) GO TO 120
  115 IER = 129
      GO TO 170
C                                  SET STEPCH TO THE GREATEST CHANGE TO
C                                    THE CURRENT VALUE OF STMIN THAT IS
C                                    ALLOWED BY THE BOUND ON THE LINE
C                                    SEARCH. SET GSPLN TO THE GRADIENT
C                                    OF THE QUADRATIC SPLINE AT
C                                    (STMIN+STEPCH). HENCE CALCULATE
C                                    THE VALUE OF STEPCH THAT MINIMIZES
C                                    THE SPLINE FUNCTION, AND THEN
C                                    OBTAIN THE NEW FUNCTION AND
C                                    GRADIENT VECTOR, FOR THIS VALUE OF
C                                    THE CHANGE TO THE STEP-LENGTH.
  120 STEPCH = 0.5D0*(SBOUND-STMIN)
      IF (SBOUND.LT.-0.5D0) STEPCH = 9.0D0*STMIN
      GSPLN = GMIN+STEPCH*DDSPLN
      IF (GMIN*GSPLN.LT.0.0D0) STEPCH = STEPCH*GMIN/(GMIN-GSPLN)
      GO TO 90
C                                  CALCULATE THE VALUE OF BETA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION.
  125 SUM = 0.0D0
      DO 130 I=1,N
  130 SUM = SUM+G(I)*W(IGINIT+I)
      BETA = (GSQRD-SUM)/(GMIN-GINIT)
C                                  TEST THAT THE NEW SEARCH DIRECTION
C                                    CAN BE MADE DOWNHILL. IF IT
C                                    CANNOT, THEN MAKE ONE ATTEMPT TO
C                                    IMPROVE THE ACCURACY OF THE LINE
C                                    SEARCH.
      IF (ABS(BETA*GMIN).LE.0.2D0*GSQRD) GO TO 135
      IRETRY = IRETRY+1
      IF (IRETRY.LE.0) GO TO 110
C                                  APPLY THE TEST THAT DEPENDS ON THE
C                                    PARAMETER MXFCON.
C                                    SET DFPR TO THE PREDICTED
C                                    REDUCTION IN F ON THE NEXT
C                                    ITERATION.
  135 IF (F.LT.FINIT) ITERFM = ITERC
      IF (ITERC.LT.ITERFM+MXFCON) GO TO 140
      IER = 132
      GO TO 9000
  140 DFPR = STMIN*GINIT
C                                  BRANCH IF A RESTART PROCEDURE IS
C                                    REQUIRED DUE TO THE ITERATION
C                                    NUMBER OR DUE TO THE SCALAR
C                                    PRODUCT OF CONSECUTIVE GRADIENTS.
      IF (IRETRY.GT.0) GO TO 10
      IF (ITERRS.EQ.0) GO TO 155
      IF (ITERC-ITERRS.GE.N) GO TO 155
      IF (ABS(SUM).GE.0.2D0*GSQRD) GO TO 155
C                                  CALCULATE THE VALUE OF GAMA THAT
C                                    OCCURS IN THE NEW SEARCH
C                                    DIRECTION, AND SET SUM TO A SCALAR
C                                    PRODUCT FOR THE TEST BELOW. THE
C                                    VALUE OF GAMDEN IS SET BY THE
C                                    RESTART PROCEDURE.
      GAMA = 0.0D0
      SUM = 0.0D0
      DO 145 I=1,N
         GAMA = GAMA+G(I)*W(IRSDG+I)
  145 SUM = SUM+G(I)*W(IRSDX+I)
      GAMA = GAMA/GAMDEN
C                                  RESTART IF THE NEW SEARCH DIRECTION
C                                    IS NOT SUFFICIENTLY DOWNHILL.
C
      IF (ABS(BETA*GMIN+GAMA*SUM).GE.0.2D0*GSQRD) GO TO 155
C
C                                  CALCULATE THE NEW SEARCH DIRECTION.
      DO 150 I=1,N
  150 W(I) = -G(I)+BETA*W(I)+GAMA*W(IRSDX+I)
      GO TO 80
C                                  APPLY THE RESTART PROCEDURE.
  155 GAMDEN = GMIN-GINIT
      DO 160 I=1,N
         W(IRSDX+I) = W(I)
         W(IRSDG+I) = G(I)-W(IGINIT+I)
  160 W(I) = -G(I)+BETA*W(I)
      ITERRS = ITERC
      GO TO 80
C                                  SET IER TO INDICATE THAT THE SEARCH
C                                    DIRECTION IS UPHILL.
  165 IER = 130
C                                  ENSURE THAT F, X AND G ARE OPTIMAL.
  170 IF (NCALLS.EQ.NFOPT) GO TO 180
      F = FMIN
      DO 175 I=1,N
         X(I) = W(IXOPT+I)
  175 G(I) = W(IGOPT+I)
  180 IF (IER.EQ.0) GO TO 125
 9000 CONTINUE
C
 9005 RETURN
      END
