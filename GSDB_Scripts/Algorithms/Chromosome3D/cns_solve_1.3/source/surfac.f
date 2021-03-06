      SUBROUTINE SURFAC
C
C For syntax see main parsing loop "HELP"
C
C Keyword ACCEss means accessible surface, whereas CONTact means
C contact area.
C
C Input are vdw radii of atoms and atomic positions.
C Output is in RMSD array.
C
C Authors: B.K. Lee, F.M. Richards, Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'coord.inc'
C local
      INTEGER FLAGS, SRAD, SX, SY, SZ
C begin
      FLAGS=ALLHP(INTEG4(NATOM))
      SRAD=ALLHP(IREAL8(NATOM))
      SX=ALLHP(IREAL8(NATOM))
      SY=ALLHP(IREAL8(NATOM))
      SZ=ALLHP(IREAL8(NATOM))
      CALL SURFA1(HEAP(FLAGS),HEAP(SRAD),
     &            HEAP(SX),HEAP(SY),HEAP(SZ))
      CALL FREHP(SZ,IREAL8(NATOM))
      CALL FREHP(SY,IREAL8(NATOM))
      CALL FREHP(SX,IREAL8(NATOM))
      CALL FREHP(SRAD,IREAL8(NATOM))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      RETURN
      END
C
      SUBROUTINE SURFA1(FLAGS,SRAD,SX,SY,SZ)
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'cnst.inc'
      INTEGER FLAGS(*)
      DOUBLE PRECISION SX(*),SY(*),SZ(*),SRAD(*)
C local
      DOUBLE PRECISION RH2O, NP
      INTEGER I, J, SNATOM, NSELCT, MODE, DIM
      INTEGER TAG, TAG1, INZ, KN, ZLB, ZUB, ZR1, YR1, XR1, RAD1
      INTEGER RSEC2, AREA, RSEC, ARCI, ARCF
      CHARACTER*4 SMODE, SRADI
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, SIX, HALF
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,SIX=6.0D0,HALF=0.5D0)
C   IICT gives the maximum number of contacts for each sphere.
      INTEGER MICT
      PARAMETER (MICT=15)
C
C begin
C
C defaults
      CALL FILL4(FLAGS,NATOM,1)
      NSELCT=NATOM
      RH2O=1.6D0
      NP=0.05D0
      MODE=1
      SRADI='SIGM'
C
C parsing
      CALL PUSEND('SURFAC>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('SURFAC>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-surface')
C
      ELSE IF (WD(1:4).EQ.'SELE') THEN
      CALL SELCTA(FLAGS,NSELCT,X,Y,Z,.TRUE.)
C======================================================================
      ELSE IF (WD(1:4).EQ.'MODE') THEN
      IF (MODE.EQ.1) THEN
      SMODE='ACCE'
      ELSE
      SMODE='CONT'
      END IF
      CALL NEXTA4('MODE=',SMODE)
      IF (SMODE.EQ.'CONT') THEN
      MODE=2
      ELSE IF (SMODE.EQ.'ACCE') THEN
      MODE=1
      ELSE
      CALL DSPERR('SURFAC','unknown value for MODE.')
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'RADI') THEN
      CALL NEXTA4('RADIus',SRADI)
      IF (SRADI.NE.'VDW'.AND.SRADI.NE.'SIGM') THEN
      CALL DSPERR('SURFAC','unknown value for RADIus.')
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'RH2O') THEN
      CALL NEXTF('RH2O=',RH2O)
C======================================================================
      ELSE IF (WD(1:4).EQ.'ACCU') THEN
      CALL NEXTF('ACCUracy=',NP)
      ELSE
      CALL CHKEND('SURFAC>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (NSELCT.GT.0) THEN
C
C
C Make sure that lookup tables for van der Waals radii are set
      IF (UPNBLK) THEN
      UPNBLK=.FALSE.
      CALL NBUPDA
      END IF
C
C map all coordinates according to selected subset
      SNATOM=0
      DO I=1,NATOM
      IF (FLAGS(I).EQ.1) THEN
      SNATOM=SNATOM+1
      SX(SNATOM)=X(I)
      SY(SNATOM)=Y(I)
      SZ(SNATOM)=Z(I)
C
      IF (SRADI.EQ.'SIGM') THEN
C RAD is set to half the sigma value
      SRAD(SNATOM)=HALF*CNBVR(LOOKUP(I),LOOKUP(I))*TWO**(-ONE/SIX)
      ELSE
C RAD is set to half the vdw value
      SRAD(SNATOM)=HALF*CNBVR(LOOKUP(I),LOOKUP(I))
      END IF
      END IF
      END DO
      IF (SRADI.EQ.'SIGM') THEN
      WRITE(6,'(A)')
     &    ' SURFAC: half the sigma value being used for atomic radii.'
      ELSE
      WRITE(6,'(A)')
     &    ' SURFAC: half the vdw value being used for atomic radii.'
      END IF
C
      DIM=MAX(SNATOM,MICT)
      TAG=ALLHP(INTEG4(DIM))
      TAG1=ALLHP(INTEG4(DIM))
      INZ=ALLHP(INTEG4(DIM))
      KN=ALLHP(INTEG4(DIM))
      ZLB=ALLHP(IREAL8(DIM))
      ZUB=ALLHP(IREAL8(DIM))
      ZR1=ALLHP(IREAL8(DIM))
      YR1=ALLHP(IREAL8(DIM))
      XR1=ALLHP(IREAL8(DIM))
      RAD1=ALLHP(IREAL8(DIM))
      RSEC2=ALLHP(IREAL8(DIM))
      AREA=ALLHP(IREAL8(DIM))
      RSEC=ALLHP(IREAL8(DIM))
      ARCI=ALLHP(IREAL8(MICT*DIM))
      ARCF=ALLHP(IREAL8(MICT*DIM))
C
      CALL SURFA2(SNATOM,MICT,HEAP(TAG),HEAP(TAG1),HEAP(INZ),
     2            HEAP(KN),HEAP(ZLB),HEAP(ZUB),
     3            HEAP(ZR1),HEAP(YR1),HEAP(XR1),HEAP(RAD1),
     4            HEAP(RSEC2),HEAP(AREA),HEAP(RSEC),
     5            HEAP(ARCI),HEAP(ARCF),SX,SY,SZ,SRAD,
     6            NP,RH2O,RMSD,MODE)
C
      CALL FREHP(ARCF,IREAL8(MICT*DIM))
      CALL FREHP(ARCI,IREAL8(MICT*DIM))
      CALL FREHP(RSEC,IREAL8(DIM))
      CALL FREHP(AREA,IREAL8(DIM))
      CALL FREHP(RSEC2,IREAL8(DIM))
      CALL FREHP(RAD1,IREAL8(DIM))
      CALL FREHP(XR1,IREAL8(DIM))
      CALL FREHP(YR1,IREAL8(DIM))
      CALL FREHP(ZR1,IREAL8(DIM))
      CALL FREHP(ZUB,IREAL8(DIM))
      CALL FREHP(ZLB,IREAL8(DIM))
      CALL FREHP(KN,INTEG4(DIM))
      CALL FREHP(INZ,INTEG4(DIM))
      CALL FREHP(TAG1,INTEG4(DIM))
      CALL FREHP(TAG,INTEG4(DIM))
C
C finally map the RMSD (=RMSD) array back.
      J=SNATOM
      DO I=NATOM,1,-1
      IF (FLAGS(I).EQ.1) THEN
      RMSD(I)=RMSD(J)
      J=J-1
      ELSE
      RMSD(I)=ZERO
      END IF
      END DO
C
      END IF
      RETURN
      END
C
      SUBROUTINE SURFA2(NATOM,ICT,TAG,TAG1,INZ,KN,ZLB,
     &                  ZUB,ZR1,YR1,XR1,RAD1,RSEC2,AREA,RSEC,
     &                  ARCI,ARCF,XR,YR,ZR,RAD,P,RH2O,ACCESS,MODE)
C
C     ACCESS: CALCULATE ACCESSIBLE CONTACT SURFACE AREA FOR A GROUP OF
C     ATOMS. THE ACCESSIBLE AREA FOR A GIVEN ATOM IS CALCULATED BY THE
C     FORMULA,
C
C         (ARCSUM) X (ATOM RADIUS+PROBE RADIUS) X (DELTAZ)
C
C     NUMERICAL INTEGRATION IS CARRIED OUT OVER Z. IN EACH Z-SECTION,
C     THE ARCSUM FOR A GIVEN ATOM IS THE ARCLENGTH OF THE CIRCLE
C     (INTERSECTION OF THE ATOM SPHERE WITH THE Z-SECTION) THAT IS NOT
C     INTERIOR TO ANY OTHER ATOM CIRCLES IN THE SAME Z-SECTION.
C
      IMPLICIT NONE
C input/ouput
      INTEGER NATOM, ICT, TAG(*), TAG1(*), INZ(*), KN(*)
      DOUBLE PRECISION ZLB(*), ZUB(*), XR1(*), YR1(*), ZR1(*), RAD1(*)
      DOUBLE PRECISION RSEC2(*), AREA(*), RSEC(*), ARCI(ICT,*)
      DOUBLE PRECISION ARCF(ICT,*)
      DOUBLE PRECISION XR(*), YR(*), ZR(*), RAD(*), P, RH2O, ACCESS(*)
      INTEGER MODE
C local
      INTEGER ANO, I, IDUM, N, K, J, J1, J2, J3, M, IEND, ITAB, ICT1
      INTEGER IANO, ICNT1, ICNT2, IFLAG
      DOUBLE PRECISION PIE, RMIN, PIEX2, ZUBMAX, RINC, ZOR, HZRES
      DOUBLE PRECISION ZRES, CUTOFF, ZNEXT, ZGRID, A, DX, DY
      DOUBLE PRECISION D, B, Q, D2
C begin
      PIE=ACOS(-1.0D0)
      PIEX2=2.0D0*PIE
      ICT1=ICT-1
      ANO=NATOM
      RMIN=10000.D0
C
      DO 7 I=1,ANO
C
C     CALCULATE LOWEST ZBOUND FOR EACH ATOM AND SORT ATOMS FROM LOW TO
C     HIGH ON ZLB
C
      ZLB(I)=ZR(I)-RAD(I)
    7 CONTINUE
      CALL SORTAG(ZLB,ANO,TAG)
      DO 6 I=1,ANO
      J=TAG(I)
      TAG1(J)=I
      XR1(I)=XR(J)
      YR1(I)=YR(J)
      ZR1(I)=ZR(J)
      RAD1(I)=RAD(J)
      IF(RAD1(I).LT.RMIN)RMIN=RAD1(I)
    6 CONTINUE
      ZUBMAX=0.0D0
C
C     THE RADIUS OF AN ATOM SPHERE = ATOM RADIUS + PROBE RADIUS
C
      RINC=RH2O
      ZOR=ZLB(1)-RINC
      DO 13 I=1,ANO
      RAD1(I)=RAD1(I)+RINC
      ZR1(I)=ZR1(I)-ZOR
      ZUB(I)=ZR1(I)+RAD1(I)
      ZLB(I)=ZR1(I)-RAD1(I)
      IF(ZUB(I).GT.ZUBMAX)ZUBMAX=ZUB(I)
      AREA(I)=0.0D0
   13 CONTINUE
C
C     Z RESOLUTION DETERMINED
C
      HZRES=(RMIN+RH2O)*P
      ZRES=2.0D0*HZRES
      CUTOFF=ZRES/100.D0
      IANO=1
      J1=1
      IFLAG=0
      ICNT1=0
      ICNT2=0
      ZNEXT=ZRES
      IEND=ZUBMAX/ZRES
C
C     SECTION ATOM SPHERES PERPENDICULAR TO THE Z AXIS
C
      DO 9 I=1,IEND
      ITAB=0
      ZGRID=ZNEXT
      ZNEXT=ZGRID+ZRES
      DO 10 N=IANO,ANO
C
C     THE UPPER AND LOWER Z BOUNDS OF AN ATOM ARE USED TO DETERMINE
C     WHETHER THIS Z-SECTION CUTS THE ATOM SPHERE
C
      IF(ZUB(N).LE.ZGRID) GO TO 21
      IF(ZLB(N).GE.ZGRID) GO TO 30
      KN(N)=0
      DO 34 K=1,ICT
      ARCI(K,N)=0.0D0
   34 CONTINUE
C
C     COUNT AND RECORD SPHERES CUT BY SECTION
C
      ITAB=ITAB+1
      INZ(ITAB)=N
C
C     FIND RADIUS OF CIRCLE LOCUS
C
      A=RAD1(N)**2-(ZGRID-ZR1(N))**2
      RSEC2(N)=A
      RSEC(N)=SQRT(A)
   21 IF(IFLAG.EQ.1)GO TO 10
C
C     FIND 1ST SPHERE CUT BY SECTION I,SKIP ATOMS PREDEEDING IT IN LIST
C     FOR NEXT SECTION
C
      IF(ZUB(N).GT.ZNEXT) GO TO 32
      J1=J1+1
      GO TO 10
   32 IFLAG=1
   10 CONTINUE
   30 IANO=J1
      IFLAG=0
C
C     ZERO, ONE, OR MORE CIRCLES ON SECTION REQUIRE DIFFERENT PROCESSING
C
      IF(ITAB.LT.1) GO TO 9
      IF(ITAB.EQ.1) GO TO 28
      J3=ITAB-1
C
C     FIND INTERSECTIONS OF CIRCLES IN SECTION CALARC CALLED TO FIND
C     INITIAL AND FINAL ANGLES OF INTERSECTION OF CIRCLES. IF ARCI AND
C     ARCF ARRAYS ARE FILLED, REDUCE IS CALLED. THE CURRENT 1ST INDEX OF
C     THESE ARRAYS FOR THE ATOM INDICATED BY 2ND INDEX IS STORED IN THE
C     ARRAY KN, IF KN FOR ANY ATOM IS 10000, THEN NO AREA REMAINS
C     ACCESSIBLE FOR THAT ATOM ON THIS SECTION OR THE ATOM IS NOT OF
C     INTEREST(KNSET=1).
C
      DO 11 K=1,J3
      N=INZ(K)
      J2=K+1
      DO 22 J=J2,ITAB
      M=INZ(J)
      A=RSEC(M)+RSEC(N)
      DX=XR1(M)-XR1(N)
      IF(ABS(DX).GE.A) GO TO 22
      DY=YR1(M)-YR1(N)
      IF(ABS(DY).GE.A) GO TO 22
      D2=DY**2+DX**2
      D=SQRT(D2)
      IF(D.GE.A)GO TO 22
      IF(KN(N).LT.ICT1)GO TO 4
      IF (KN(N).LE.ICT) THEN
      CALL REDUCE(N,0,ICT,ARCI,ARCF,TAG,TAG1,ZR1,RAD1,AREA,
     2            INZ,KN,ICNT1,ICNT2,
     3            PIE,PIEX2,ICT1,ZRES,RH2O,ZGRID,HZRES,CUTOFF)
      END IF
    4 KN(N)=KN(N)+1
      IF(KN(M).LT.ICT1)GO TO 5
      IF (KN(M).LE.ICT) THEN
      CALL REDUCE(M,0,ICT,ARCI,ARCF,TAG,TAG1,ZR1,RAD1,AREA,
     2            INZ,KN,ICNT1,ICNT2,
     3            PIE,PIEX2,ICT1,ZRES,RH2O,ZGRID,HZRES,CUTOFF)
      END IF
    5 KN(M)=KN(M)+1
C
C     DO THE CIRCLES INTERSECT, OR IS ONE COMPLETELY INSIDE THE OTHER?
C
      B=RSEC(M)-RSEC(N)
      IF(D.GT.ABS(B))GO TO 20
      IF(B.GT.0.0D0)GO TO 12
      KN(M)=10000
      KN(N)=KN(N)-1
      GO TO 22
   12 KN(N)=10000
      KN(M)=KN(M)-1
      GO TO 22
C
C     IF THE CIRCLES INTERSECT, FIND THE POINTS OF INTERSECTION
C
   20 Q=RSEC2(M)-RSEC2(N)
      D=2.0D0*D
      IF(KN(M).GT.ICT)GO TO 45
      CALL CALARC(M,1,ICT,D2,Q,D,RSEC,ARCI,ARCF,KN,DY,DX,PIE,PIEX2)
   45 IF(KN(N).GT.ICT)GO TO 22
      CALL CALARC(N,-1,ICT,D2,Q,D,RSEC,ARCI,ARCF,KN,DY,DX,PIE,PIEX2)
   22 CONTINUE
   11 CONTINUE
C
C     FIND THE ACCESSIBLE CONTACT SURFACE AREA FOR ALL THE SPHERES
C     INTERSECTING THIS SECTION
C
   28 CALL REDUCE(IDUM,ITAB,ICT,ARCI,ARCF,TAG,TAG1,ZR1,RAD1,
     2            AREA,INZ,KN,ICNT1,ICNT2,
     3            PIE,PIEX2,ICT1,ZRES,RH2O,ZGRID,HZRES,CUTOFF)
    9 CONTINUE
C
C     OUTPUT OPERATION PARAMETERS
C
      IF (MODE.EQ.1) THEN
      WRITE(6,9000) ' SURFAC: ACCEssible surface area'
      ELSE
      WRITE(6,9000) ' SURFAC: CONTact area'
      END IF
      WRITE(6,9000)
     1       ' SURFAC: ACCUracy=',P,' RH2O=',RH2O,
     2       '         Z-grid=',ZRES,' number-of-Z-sections=',IEND,
     3       '         measures-of-arc=',ICNT1,' and',ICNT2
9000  FORMAT(A,F5.2,A,F6.2,/,A,F6.2,A,I5,/,A,I6,A,I6)
C
      DO 31 J=1,ANO
      I=TAG1(J)
C
C     SCALE AREA TO VDW SHELL IF NECESSARY
C
      IF (MODE.EQ.2) THEN
      ACCESS(J)=AREA(I)*((RAD1(I)-RH2O)/RAD1(I))**2
      ELSE
      ACCESS(J)=AREA(I)
      END IF
   31 CONTINUE
      RETURN
      END
C
      SUBROUTINE CALARC(M,ISIGN,ICT,D2,Q,D,RSEC,ARCI,ARCF,KN,DY,DX,
     2                  PIE,PIEX2)
C
C     INITIAL AND FINAL ARC ENDPOINTS ARE FOUND FOR A REFERENCE CIRCLE
C     INTERSECTED BY ANOTHER CIRCLE CONTAINED IN THE SAME PLANE. THE
C     INITIAL ENDPOINT OF THE ENCLOSED ARC IS STORED IN ARCI, AND THE
C     FINAL ARC IN ARCF
C
      IMPLICIT NONE
C input/ouput
      INTEGER M
      INTEGER ISIGN
      INTEGER ICT
      DOUBLE PRECISION D2, Q, D, RSEC(*), ARCI(ICT,*), ARCF(ICT,*)
      INTEGER KN(*)
      DOUBLE PRECISION DY, DX, PIE, PIEX2
C local
      INTEGER K1
      DOUBLE PRECISION ARG, ALPHA1, BETA1, TI, TF
C begin
      K1=KN(M)
C
C     LAW OF COSINES
C
      ARG=(D2+Q*ISIGN)/(D*RSEC(M))
      ALPHA1=ACOS(ARG)
C
C     ALPHA1 IS THE ANGLE BETWEEN A LINE CONTAINING A POINT OF
C     INTERSECTION AND THE REFERENCE CIRCLE CENTER AND THE LINE
C     CONTAINING BOTH CIRCLE CENTERS
C
      BETA1=ATAN2(DY,DX)
C
C     BETA1 IS THE ANGLE BETWEEN THE LINE CONTAINING BOTH CIRCLE CENTERS
C     AND THE X-AXIS
C
      IF(ISIGN.EQ.1)BETA1=BETA1+PIE
      TI=BETA1-ALPHA1
      TF=BETA1+ALPHA1
      IF(TI.LT.0.0D0)TI=TI+PIEX2
      IF(TF.GT.PIEX2)TF=TF-PIEX2
      IF(TF.LT.0.0D0)TF=TF+PIEX2
      IF(TF.GE.TI)GO TO 3
C
C     IF THE ARC CROSSES ZERO, THEN BREAK IT INTO TWO SEGMENTS. THE
C     FIRST ENDS AT 2XPI AND THE SECOND BEGINS AT ZERO
C
      ARCF(K1+1,M)=TF
      ARCF(K1,M)=PIEX2
      KN(M)=KN(M)+1
      GO TO 2
    3 ARCF(K1,M)=TF
    2 ARCI(K1,M)=TI
      RETURN
      END
C
      SUBROUTINE REDUCE(N,ITAB,ICT,ARCI,ARCF,TAG,TAG1,ZR1,RAD1,
     2                  AREA,INZ,KN,ICNT1,ICNT2,
     3                  PIE,PIEX2,ICT1,ZRES,RH2O,ZGRID,HZRES,CUTOFF)
C
C     1) ITAB=0, REMOVE DEGENERACIES IN ARCI AND ARCF FOR A SINGLE ATOM
C     SINCE THESE ARRAYS ARE FILLED, OR 2) ITAB>0, CALCULATE THE
C     ACCESSIBLE SURFACE AREA FOR THE N ATOMS IN THIS SECTION
C
      IMPLICIT NONE
C input/output
      INTEGER N, ITAB, ICT
      DOUBLE PRECISION    ARCI(ICT,*),ARCF(ICT,*)
      INTEGER TAG(*),TAG1(*)
      DOUBLE PRECISION    ZR1(*),RAD1(*),AREA(*)
      INTEGER INZ(*), KN(*), ICNT1, ICNT2
      DOUBLE PRECISION    PIE, PIEX2
      INTEGER ICT1
      DOUBLE PRECISION    ZRES, RH2O, ZGRID, HZRES, CUTOFF
C local
      INTEGER II, K1, JJ, K, M, I, I1, J
      DOUBLE PRECISION    ARCSUM, T, A, PAREA
C begin
      DO 23 II=1,MAX(1,ITAB)
C
C     N IS PASSED IN THE SUBROUTINE ARGUMENT LIST FOR 1), OR IN INZ FOR
C     2)
C
      IF(ITAB.NE.0)N=INZ(II)
      K1=KN(N)
      IF(K1.GT.ICT)GO TO 23
C
C     IF K1>ICT, THIS ATOM IS NOT OF INTEREST OR IT HAS NO AREA
C     REMAINING IS THIS CIRCLE INTERSECTED BY OTHERS?
C
      IF(K1.NE.0) GO TO 19
C
C     THERE IS ONLY ONE CIRCLE IN THIS SECTION
C
      ARCSUM=PIEX2
      GO TO 25
C
C     THE ARC ENDPOINTS ARE SORTED ON THE VALUE OF THE INITIAL ARC
C     ENDPOINT
C
   19 CALL SORTAG(ARCI(1,N),K1,TAG)
C
C     CALCULATE THE ACCESSIBLE AREA
C
      ARCSUM=ARCI(1,N)
      JJ=TAG(1)
      T=ARCF(JJ,N)
      IF(K1.EQ.1) GO TO 8
      DO 27 K=2,MAX(2,K1)
      IF(T.GE.ARCI(K,N))GO TO 39
      ARCSUM=ARCSUM+ARCI(K,N)-T
      JJ=TAG(K)
      T=ARCF(JJ,N)
      GO TO 27
   39 M=TAG(K)
      IF(ARCF(M,N).GT.T)T=ARCF(M,N)
   27 CONTINUE
    8 ARCSUM=ARCSUM+PIEX2-T
   25 A=ZR1(N)-ZGRID
      A=RAD1(N)-ABS(A)
C
C     THE AREA IS EQUAL TO THE ACCESSIBLE ARC LENGTH X THE SECTION
C     THICKNESS, CORRECTED IF IT IS THE FIRST OR LAST SECTION, X THE
C     RADIUS OF THE SPHERE
C
      PAREA=ARCSUM*(HZRES+MIN(A,HZRES))*RAD1(N)
C
C     ADD THE ACCESSIBLE AREA FOR THIS ATOM IN THIS SECTION TO THE AREA
C     FOR THIS ATOM FOR ALL THE SECTION ENCOUNTERED THUS FAR
C
      AREA(N)=AREA(N)+PAREA
   23 CONTINUE
C
C     IF THIS WAS THE FINAL ATOM FOR THIS Z SECTION RETURN TO ACCESS
C
      IF(ITAB.GT.0)GO TO 1
C
C     THE ARCI AND ARCF ARRAYS WERE FILLED, DOES SIGNIFICANT AREA
C     REMAIN?
C
      IF(PAREA.GT.CUTOFF)GO TO 47
C
C     NO, THIS ATOM IS USED ONLY IN CALCS FOR OTHER ATOMS, UNTIL NEXT Z
C     SECTION
C
      KN(N)=10000
      ICNT1=ICNT1+1
      GO TO 1
C
C     YES, REMOVE DEGENERACIES FROM ARCI AND ARCF (COMBINE OVERLAPPING
C     ARCS)
C
   47 JJ=TAG(1)
      T=ARCF(JJ,N)
      I=1
      DO 2 K=2,MAX(2,K1)
      IF(T.GE.ARCI(K,N))GO TO 3
      I=I+1
      IF(I.EQ.ICT1) GO TO 4
      I1=I-1
      DO 7 J=K,MAX(K,K1)
      IF(TAG(J).EQ.I1)GO TO 6
    7 CONTINUE
      GO TO 9
    6 ARCF(JJ,N)=ARCF(I1,N)
      TAG(J)=JJ
    9 ARCF(I1,N)=T
      ARCI(I,N)=ARCI(K,N)
      JJ=TAG(K)
      T=ARCF(JJ,N)
      GO TO 2
    3 M=TAG(K)
      IF(ARCF(M,N).GT.T)T=ARCF(M,N)
    2 CONTINUE
      ARCF(I,N)=T
      KN(N)=I
      I=I+1
      DO 5 K=I,MAX(K1,I)
      ARCI(K,N)=0.0D0
    5 CONTINUE
C
C     REMOVE THE PARTIAL AREA ADDED, THIS CIRCLE MAY HAVE MORE
C     INTERSECTIONS
C
      AREA(N)=AREA(N)-PAREA
      ICNT2=ICNT2+1
      GO TO 1
C
C     NO DEGENERACIES WERE FOUND IN THE OVERLAPS, START OVER WITH LARGER
C     ICT
C
    4 WRITE(6,15) TAG1(N),ICT
   15 FORMAT('0ERROR IN SURFA1 - ATOM ',I5,
     2       ' HAS MORE THAN ',I3,' CONTACTS')
      CALL DIE
    1 RETURN
      END
C
      SUBROUTINE SORTAG(A,N,TAG)
C
      IMPLICIT NONE
      INTEGER N,TAG(*)
      DOUBLE PRECISION A(N)
C
      DOUBLE PRECISION T, TT
      INTEGER I, L, M, J, K, IJ
      INTEGER TG, IU(16), IL(16)
C
      DO 1  I=1,N
      TAG(I)=I
    1 CONTINUE
      M=1
      I=1
      J=N
  5   IF(I.GE.J) GO TO 70
 10   K=I
      IJ=(J+I)/2
      T=A(IJ)
      IF(A(I).LE.T) GO TO 20
      A(IJ)= A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
 20   L=J
      IF(A(J).GE.T) GO TO 40
      A(IJ)=A(J)
      A(J)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(J)
      TAG(J)=TG
      IF(A(I).LE.T) GO TO 40
      A(IJ)=A(I)
      A(I)=T
      T=A(IJ)
      TG=TAG(IJ)
      TAG(IJ)=TAG(I)
      TAG(I)=TG
      GO TO 40
 30   A(L)=A(K)
      A(K)=TT
      TG=TAG(L)
      TAG(L)=TAG(K)
      TAG(K)=TG
 40   L=L-1
      IF(A(L).GT.T) GO TO 40
      TT=A(L)
 50   K=K+1
      IF(A(K).LT.T) GO TO 50
      IF(K.LE.L) GO TO 30
      IF(L-I.LE.J-K) GO TO 60
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 80
 60   IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 80
 70   M=M-1
      IF(M.EQ.0) RETURN
      I=IL(M)
      J=IU(M)
 80   IF(J-I.GE.1) GO TO 10
      IF(I.EQ.1) GO TO 5
      I=I-1
 90   I=I+1
      IF(I.EQ.J) GO TO 70
      T=A(I+1)
      IF(A(I).LE.T) GO TO 90
      TG=TAG(I+1)
      K=I
 100  A(K+1)=A(K)
      TAG(K+1)=TAG(K)
      K=K-1
      IF(T.LT.A(K)) GO TO 100
      A(K+1)=T
      TAG(K+1)=TG
      GO TO 90
      END
