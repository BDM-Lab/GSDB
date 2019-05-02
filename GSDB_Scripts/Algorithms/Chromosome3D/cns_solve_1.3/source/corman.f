      SUBROUTINE CORMAN
C
C This routine handles all of the coordinate read and
C manipulation options
C
C For syntax see main parsing loop "HELP"
C
C Remarks: all operations are performed on the main coordinate set
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER ISLCT
C begin
      ISLCT=ALLHP(INTEG4(NATOM))
      CALL CORMA2(HEAP(ISLCT))
      CALL FREHP(ISLCT,INTEG4(NATOM))
      RETURN
      END
C
      SUBROUTINE CORMA2(ISLCT)
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'xcrystal.inc'
      INTEGER ISLCT(*)
C local
      INTEGER NMISS, NSELCT, I, J, K, ISCR
      DOUBLE PRECISION XCM, YCM, ZCM, WW, TMASS
      DOUBLE PRECISION XN, YN, ZN, FACT, RMST, RMSTT, AMASSV, AMASST
      DOUBLE PRECISION RMSV, RG, AR, TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION R(3), DIST, CR(3), XNN, YNN, ZNN, EU(3)
      DOUBLE PRECISION ROT(3,3), PHI, A, Q(4)
      DOUBLE PRECISION SYFR(3,4), SY(3,4), XX, YY, ZZ, RTH
      LOGICAL LRMS, LNORO, LMASS, QDIST, QUL, QCHAIN
      CHARACTER*4 SDISP, ACTION, SREFE, MODE
      DOUBLE COMPLEX DBCOMP
      INTEGER CRNSYM, CRSYMM(1,3,4), CRITSY(1,3,3)
C parameters
      DOUBLE PRECISION ANUM, RAD, ZERO, NINETY, ONE
      PARAMETER (ANUM=9999.0D0, RAD=PI/180.0D0, ZERO=0.0D0,NINETY=90.D0)
      PARAMETER (ONE=1.0D0)
C
C begin
C
C defaults
      NMISS=0
      ACTION='READ'
      SDISP='MAIN'
      SREFE='MAIN'
      CALL FILL4(ISLCT,NATOM,1)
      NSELCT=NATOM
      LNORO=.FALSE.
      LMASS=.FALSE.
      QDIST=.FALSE.
      QCHAIN=.TRUE.
      PHI=ZERO
      DIST=ZERO
      FACT=ZERO
      QUL=.FALSE.
      DO I=1,3
      R(I)=ZERO
      END DO
      DO I=1,3
      CR(I)=ZERO
      END DO
C
      CALL PUSEND('COOR>')
      DO WHILE (.NOT.DONE)
      IF (ACTION.NE.'READ') THEN
      CALL NEXTWD('COOR-'//ACTION//'>')
      CALL MISCOM('COOR-'//ACTION//'>',USED)
      ELSE
      CALL NEXTWD('COOR>')
      CALL MISCOM('COOR>',USED)
      END IF
      IF (.NOT.USED) THEN
      IF (ACTION.NE.'READ') THEN
      ELSE IF (SDISP.EQ.'MAIN') THEN
      CALL CREAD(X,Y,Z,WMAIN,QMAIN,ISLCT,QCHAIN)
      ELSE IF (SDISP.EQ.'COMP') THEN
      CALL CREAD(XCOMP,YCOMP,ZCOMP,WCOMP,QCOMP,ISLCT,QCHAIN)
      ELSE IF (SDISP.EQ.'REFE') THEN
      CALL CREAD(REFX,REFY,REFZ,KCNSTR,KCNSTR,ISLCT,QCHAIN)
      ELSE IF (SDISP.EQ.'DERI') THEN
      CALL CREAD(DX,DY,DZ,FBETA,FBETA,ISLCT,QCHAIN)
      END IF
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-coordinate')
C
      ELSE IF (WD(1:4).EQ.'DISP'.AND.ACTION.EQ.'READ') THEN
      CALL NEXTA4('DISPosition=',SDISP)
      IF (.NOT.(SDISP.EQ.'MAIN'.OR.SDISP.EQ.'COMP'
     &   .OR.SDISP.EQ.'REFE'.OR.SDISP.EQ.'DERI')) THEN
      CALL DSPERR('COOR','unknown value for disposition')
      END IF
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(3A)') ' COOR: coordinates read into ',SDISP,' set'
      END IF
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(ISLCT,NSELCT,X,Y,Z,.TRUE.)
      IF (NSELCT.LT.NATOM.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' COOR: using atom subset.'
      END IF
C
      ELSE IF (WD(1:4).EQ.'CONV') THEN
      CALL NEXTLO('CONVert-chainid-to-segid=',QCHAIN)
C
      ELSE IF (WD(1:4).EQ.'READ') THEN
      ACTION='READ'
      ELSE IF (WD(1:4).EQ.'REFE') THEN
      CALL NEXTA4('REFErence=',SREFE)
      ELSE IF (WD(1:4).EQ.'VECT') THEN
      CALL NEXTVF('VECTor=',R)
      ELSE IF (WD(1:4).EQ.'CENT') THEN
      CALL NEXTVF('CENTer-of-rotation=',CR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'MATR'.OR.WD(1:4).EQ.'EULE'.OR.
     &         WD(1:4).EQ.'LATT'.OR.WD(1:4).EQ.'SPHE'.OR.
     &         WD(1:4).EQ.'AXIS'.OR.WD(1:4).EQ.'QUAT') THEN
      MODE=WD(1:4)
      CALL MATPAR(MODE,ROT)
C
C set the flag that disables the old VECTOR ... ANGLE statement
      QUL=.TRUE.
C
C write some info
      IF (WRNLEV.GE.5) CALL MATPRI(ROT)
      CALL MATDCL(ROT)
C==================================================================
C
      ELSE IF (WD(1:4).EQ.'ANGL') THEN
      CALL NEXTF('ANGLe=',PHI)
      ELSE IF (WD(1:4).EQ.'DIST') THEN
      CALL NEXTF('DISTance=',DIST)
      QDIST=.TRUE.
      ELSE IF (WD(1:4).EQ.'FACT') THEN
      CALL NEXTF('FACTor=',FACT)
      ELSE IF (WD(1:4).EQ.'MASS') THEN
      CALL NEXTLO('MASS=',LMASS)
      ELSE IF (WD(1:3).EQ.'LSQ') THEN
      LNORO=.NOT.LNORO
      CALL NEXTLO('LSQ=',LNORO)
      LNORO=.NOT.LNORO
C
      ELSE IF (WD(1:4).EQ.'SWAP'.AND.ACTION.EQ.'READ') THEN
      ACTION='SWAP'
      ELSE IF (WD(1:4).EQ.'INIT'.AND.ACTION.EQ.'READ') THEN
      ACTION='INIT'
      ELSE IF (WD(1:4).EQ.'COPY'.AND.ACTION.EQ.'READ') THEN
      ACTION='COPY'
      ELSE IF (WD(1:4).EQ.'ROTA'.AND.ACTION.EQ.'READ') THEN
      ACTION='ROTA'
      ELSE IF (WD(1:4).EQ.'SYMM'.AND.ACTION.EQ.'READ') THEN
      ACTION='SYMM'
      CALL NEXTEX('SYMMetry=',WDD,WDMAX,WDDLEN)
      CRNSYM=0
      CALL XRSYPA(WDD,WDDLEN,1,CRNSYM,CRSYMM,CRITSY,XRSYTH)
C
      IF (CRNSYM.GT.0) THEN
C
C concatenate symmetry operator and orthogonal-to-fractional operator
      DO I=1,3
      DO J=1,3
      SYFR(I,J)=ZERO
      DO K=1,3
      SYFR(I,J)=SYFR(I,J)+CRSYMM(CRNSYM,I,K)*XRTR(K,J)
      END DO
      END DO
      END DO
      RTH=ONE/XRSYTH
      DO I=1,3
      SYFR(I,4)=CRSYMM(CRNSYM,I,4)*RTH
      END DO
C
C apply fractional-to-orthogonal operator
      DO I=1,3
      DO J=1,4
      SY(I,J)=ZERO
      DO K=1,3
      SY(I,J)=SY(I,J)+XRINTR(I,K)*SYFR(K,J)
      END DO
      END DO
      END DO
C
      ELSE
      DO I=1,3
      DO J=1,4
      SY(I,J)=ZERO
      END DO
      END DO
      DO J=1,3
      SY(J,J)=ONE
      END DO
      END IF
C
      ELSE IF (WD(1:4).EQ.'TRAN'.AND.ACTION.EQ.'READ') THEN
      ACTION='TRAN'
      ELSE IF (WD(1:4).EQ.'ORIE'.AND.ACTION.EQ.'READ') THEN
      ACTION='ORIE'
      ELSE IF (WD(1:3).EQ.'FIT'.AND.ACTION.EQ.'READ') THEN
      ACTION='FIT'
      ELSE IF (WD(1:3).EQ.'RMS'.AND.ACTION.EQ.'READ') THEN
      ACTION='RMS'
      ELSE IF (WD(1:4).EQ.'ORTH'.AND.ACTION.EQ.'READ') THEN
      ACTION='ORTH'
      ELSE IF (WD(1:4).EQ.'FRAC'.AND.ACTION.EQ.'READ') THEN
      ACTION='FRAC'
      ELSE IF (WD(1:4).EQ.'RGYR'.AND.ACTION.EQ.'READ') THEN
      ACTION='RGYR'
      ELSE
      IF (ACTION.NE.'READ') THEN
      CALL CHKEND('COOR-'//ACTION//'>',DONE)
      ELSE
      CALL CHKEND('COOR>',DONE)
      END IF
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (ACTION.EQ.'SWAP') THEN
C
C swaps the comparison coordinates with the main coordinates
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z).OR..NOT.INITIA(I,XCOMP,YCOMP,ZCOMP)) THEN
      NMISS=NMISS+1
      END IF
      A=Y(I)
      Y(I)=YCOMP(I)
      YCOMP(I)=A
      A=Z(I)
      Z(I)=ZCOMP(I)
      ZCOMP(I)=A
      A=X(I)
      X(I)=XCOMP(I)
      XCOMP(I)=A
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' COOR: selected coordinates swaped (main<->comp)'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'INIT') THEN
C
C initializes the selected main coordinates
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      X(I)=ANUM
      Y(I)=ANUM
      Z(I)=ANUM
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' COOR: selected main coordinates initialized'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'COPY') THEN
C
C copies main coordinate set into comparsion set
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF(.NOT.INITIA(I,X,Y,Z)) NMISS=NMISS+1
      XCOMP(I)=X(I)
      YCOMP(I)=Y(I)
      ZCOMP(I)=Z(I)
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)') ' COOR: selected main coordinates copied to comp'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'ORIE') THEN
C
C orient selected atoms in space
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1.AND..NOT.INITIA(I,X,Y,Z)) THEN
      NSELCT=NSELCT-1
      ISLCT(I)=0
      NMISS=NMISS+1
      END IF
      END DO
C
      IF (NSELCT.GT.0) THEN
      ISCR=ALLHP(INTEG4(NATOM+NATOM))
      LRMS=.FALSE.
      CALL ORINTC(NATOM,X,Y,Z,XCOMP,YCOMP,ZCOMP,AMASS,LMASS,.FALSE.,
     2     HEAP(ISCR),ISLCT,LNORO)
      CALL FREHP(ISCR,INTEG4(NATOM+NATOM))
      END IF
C
      END IF
C
      ELSE IF (ACTION.EQ.'FIT') THEN
C
C make lsq fit of main and comp coordinate set
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z).OR..NOT.INITIA(I,XCOMP,YCOMP,ZCOMP)) THEN
      NSELCT=NSELCT-1
      ISLCT(I)=0
      NMISS=NMISS+1
      END IF
      END IF
      END DO
C
      IF (NSELCT.GT.0) THEN
      ISCR=ALLHP(INTEG4(NATOM+NATOM))
      LRMS=.FALSE.
      CALL ORINTC(NATOM,X,Y,Z,XCOMP,YCOMP,ZCOMP,AMASS,LMASS,.TRUE.,
     &     HEAP(ISCR),ISLCT,LNORO)
      CALL FREHP(ISCR,INTEG4(NATOM+NATOM))
      END IF
C
      END IF
C
      ELSE IF (ACTION.EQ.'ROTA') THEN
C
C find rotation matrix, if not explicitly specified by vector/angle mode
      IF (.NOT.QUL) THEN
      EU(3)=PHI
      CALL ROTMAT(ROT,EU(1),EU(2),EU(3),Q,R,'AXIS')
C
C print out some information
      IF (WRNLEV.GE.5) CALL MATPRI(ROT)
      CALL MATDCL(ROT)
      END IF
C
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE
      XNN=X(I)-CR(1)
      YNN=Y(I)-CR(2)
      ZNN=Z(I)-CR(3)
      XN=ROT(1,1)*XNN+ROT(1,2)*YNN+ROT(1,3)*ZNN
      YN=ROT(2,1)*XNN+ROT(2,2)*YNN+ROT(2,3)*ZNN
      ZN=ROT(3,1)*XNN+ROT(3,2)*YNN+ROT(3,3)*ZNN
      X(I)=XN+CR(1)
      Y(I)=YN+CR(2)
      Z(I)=ZN+CR(3)
      END IF
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,3F12.6)')
     &  ' COOR: rotation center =',CR(1),CR(2),CR(3)
      WRITE(6,'(A)') ' COOR: selected coordinates rotated'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'SYMM') THEN
C
C apply symmetry operator to selected coordinates
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE
      XX=SY(1,1)*X(I) +SY(1,2)*Y(I)
     &         +SY(1,3)*Z(I) +SY(1,4)
      YY=SY(2,1)*X(I) +SY(2,2)*Y(I)
     &         +SY(2,3)*Z(I) +SY(2,4)
      ZZ=SY(3,1)*X(I) +SY(3,2)*Y(I)
     &         +SY(3,3)*Z(I) +SY(3,4)
      X(I)=XX
      Y(I)=YY
      Z(I)=ZZ
      END IF
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' COOR: applied symmetry operator to selected coordinates'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'TRAN') THEN
C
C translate selected coordinates
      IF (QDIST) THEN
      AR=(R(1)**2+R(2)**2+R(3)**2)
      IF (AR.GT.RSMALL) THEN
      DIST=DIST/SQRT(AR)
      ELSE
      DIST=ZERO
      END IF
      DO I=1,3
      R(I)=R(I)*DIST
      END DO
      END IF
C
      IF (NSELCT.GT.0) THEN
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE
      X(I)=X(I)+R(1)
      Y(I)=Y(I)+R(2)
      Z(I)=Z(I)+R(3)
      END IF
      END IF
      END DO
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,3F12.6,A)')
     &    ' COOR: translation vector =( ',(R(I),I=1,3),' )'
      WRITE(6,'(A)') ' COOR: selected coordinates translated'
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'RMS') THEN
C
C compute RMS differences for individual atoms
      IF (NSELCT.GT.0) THEN
      AMASST=ZERO
      RMST=ZERO
      DO I=1,NATOM
      IF (ISLCT(I).EQ.1) THEN
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE
      IF (.NOT.INITIA(I,XCOMP,YCOMP,ZCOMP)) THEN
      NMISS=NMISS+1
      ELSE
      IF (LMASS) THEN
      AMASSV=AMASS(I)
      ELSE
      AMASSV=ONE
      END IF
      AMASST=AMASST+AMASSV
      RMSTT=(X(I)-XCOMP(I))**2
      RMSTT=RMSTT+(Y(I)-YCOMP(I))**2
      RMSTT=RMSTT+(Z(I)-ZCOMP(I))**2
      RMST=RMSTT*AMASSV+RMST
      RMSD(I)=SQRT(RMSTT)
      END IF
      END IF
      END IF
      END DO
      IF (AMASST.GT.ZERO) THEN
      RMSV=SQRT(RMST/AMASST)
      ELSE
      RMSV=ZERO
      END IF
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' COOR: ATTENTION: atom rms differences copied to RMSD array. '
      WRITE(6,'(A,F12.6,A,F12.4,A,F12.4)')
     & ' COOR: rms= ',RMSV,', square sum= ',RMST,
     &              ', denominator= ',AMASST
      END IF
      CALL DECLAR( 'RESULT', 'DP', ' ', DBCOMP, RMSV )
      CALL DECLAR( 'RMS', 'DP', ' ', DBCOMP, RMSV )
      END IF
C
      ELSE IF (ACTION.EQ.'ORTH') THEN
C
      IF (NSELCT.GT.0) THEN
      NMISS=0
C
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE IF (ISLCT(I).EQ.1) THEN
      TEMPX=X(I)
      TEMPY=Y(I)
      TEMPZ=Z(I)
      X(I)=XRINTR(1,1)*TEMPX +XRINTR(1,2)*TEMPY +XRINTR(1,3)*TEMPZ
      Y(I)=XRINTR(2,1)*TEMPX +XRINTR(2,2)*TEMPY +XRINTR(2,3)*TEMPZ
      Z(I)=XRINTR(3,1)*TEMPX +XRINTR(3,2)*TEMPY +XRINTR(3,3)*TEMPZ
      END IF
      END DO
C
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,/,3(1X,F12.6),/,3(1X,F12.6),/,3(1X,F12.6))')
     & ' ORTHO: transformation matrix =',((XRINTR(I,J),J=1,3),I=1,3)
      END IF
C
      END IF
C
      ELSE IF (ACTION.EQ.'FRAC') THEN
C
      IF (NSELCT.GT.0) THEN
      NMISS=0
C
      DO I=1,NATOM
      IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE IF (ISLCT(I).EQ.1) THEN
      TEMPX=X(I)
      TEMPY=Y(I)
      TEMPZ=Z(I)
      X(I)=XRTR(1,1)*TEMPX +XRTR(1,2)*TEMPY +XRTR(1,3)*TEMPZ
      Y(I)=XRTR(2,1)*TEMPX +XRTR(2,2)*TEMPY +XRTR(2,3)*TEMPZ
      Z(I)=XRTR(3,1)*TEMPX +XRTR(3,2)*TEMPY +XRTR(3,3)*TEMPZ
      END IF
      END DO
C
      IF (WRNLEV.GE.5) THEN
CCC modification ATB 4/27/08
      WRITE(6,'(A,/,3(1X,F12.6),/,3(1X,F12.6),/,3(1X,F12.6))')
     & ' FRACT: transformation matrix =',((XRTR(I,J),J=1,3),I=1,3)
      END IF
      END IF
C
      ELSE IF (ACTION.EQ.'RGYR') THEN
C
C compute radius of gyration
      IF (NSELCT.GT.0) THEN
C
C Center-of-mass:
      XCM=ZERO
      YCM=ZERO
      ZCM=ZERO
      TMASS=ZERO
      NMISS=0
C
      DO I=1,NATOM
      IF (ISLCT(I).NE.1) THEN
      ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
      NMISS=NMISS+1
      ELSE
      IF (LMASS) THEN
      WW=AMASS(I)
      ELSE
      WW=ONE
      END IF
      WW=WW-FACT
      XCM=WW*X(I)+XCM
      YCM=WW*Y(I)+YCM
      ZCM=WW*Z(I)+ZCM
      TMASS=WW+TMASS
      END IF
      END DO
C
      IF (TMASS .LE. RSMALL) THEN
      WRITE(6,'(A,F12.5)') ' %COOR-RGYR-ERR: Net "mass"=',TMASS
      ELSE
      XCM=XCM/TMASS
      YCM=YCM/TMASS
      ZCM=ZCM/TMASS
C
C Radius of gyration:
      RG=ZERO
      DO I=1,NATOM
      IF (ISLCT(I).NE.1) THEN
      ELSE IF (.NOT.INITIA(I,X,Y,Z)) THEN
      ELSE
      IF (LMASS) THEN
      WW=AMASS(I)
      ELSE
      WW=ONE
      END IF
      WW=WW-FACT
      RG=RG+WW*((X(I)-XCM)**2+(Y(I)-YCM)**2+(Z(I)-ZCM)**2)
      END IF
      END DO
      AR=ABS(RG/TMASS)
C
C Compute an RG with the same sign as RG/TMASS:
      RG=SQRT(AR)*RG/(TMASS*AR)
C
      WRITE(6,'(A,F12.5,A,F12.5,/,A,F12.5,A,F12.5,A,F12.5,A)')
     & ' COOR: Radius of gyration=',RG,' Net "mass"=',TMASS,
     & ' COOR: center-of-"mass" = (' ,XCM,' ',YCM,' ',ZCM,' )'
C
      CALL DECLAR('RG', 'DP', ' ', (1.0D0, 0.0D0), RG)
      CALL DECLAR('XCM', 'DP', ' ', (1.0D0, 0.0D0), XCM)
      CALL DECLAR('YCM', 'DP', ' ', (1.0D0, 0.0D0), YCM)
      CALL DECLAR('ZCM', 'DP', ' ', (1.0D0, 0.0D0), ZCM)
C
      END IF
      END IF
      END IF
C
C
      IF (NMISS.GT.0) WRITE(6,22) NMISS
22    FORMAT(' %COOR-ERR: ',I5,' missing coordinates encountered')
C
      RETURN
      END
