      PROGRAM RHOMAP
C     --------------
C     
C---- Modified to read "ZYX" formatted CNS maps
C---- 4/15/91 SRH

      DIMENSION CL(50),IORGX(3),INCX(3),NG(3),
     $     NP(3),ZAPP(12),IALPHA(3)
      real*8 rho(400,400), cellDP(6)
C     COMMON /A/ RHO(400,400)
      common /A/ rho
      COMMON /LIMIT/ PLTMAX
      COMMON /SKEW / SB,CB,BETA
      common / lines / contourLineWidth, borderLineWidth
      CHARACTER*80 LABELS,LABELX,LABELY
      CHARACTER*80 MAPFILE, PLOTFILE
      CHARACTER*82 NAME
      character*80 title(400)
      character*3 zyx
      CHARACTER*2 NUMBER
      character*1 ialpha
      INTEGER STLEN

      data contourLineWidth, borderLineWidth / 0.50, 0.50 /


      idisk = 26

      LABELX = ' '
      LABELY = ' '
      LABELS = ' '
      READ (5,'(A)') MAPFILE
      READ (5,'(A)') PLOTFILE
      READ (5,*) SCALE,PLTMAX
      READ (5,*) CTRMIN,CTRINT,ISP,AXMOD
      READ (5,*) KK,KMIN,KINC,KSKIP
      READ (5,*) LINWID,ITICK
C     IF (ITICK.LE.0) ITICK=0
      LAND=0
      IF (LINWID.LT.0) THEN
         LAND=1
         LINWID=ABS(LINWID)
      ENDIF
      IF (KMIN.LE.0) KMIN = 1
      IF (KINC.LE.0) KINC = 1
      IF (KSKIP.LE.0) KSKIP = 1
      NLB=80
      READ (5,'(A)') LABELS
      WRITE (6,80000) SCALE,PLTMAX,CTRMIN,CTRINT,ISP,AXMOD
80000 FORMAT (1H ,'Input scale ..................... ',F10.4/
     $        1H ,'Maximum plot size ............... ',F10.2/
     $        1H ,'Minimum contur value ............ ',F10.2/
     $        1H ,'Contur increment ................ ',F10.2/
     $        1H ,'Spline constant ................. ',I10/
     $        1H ,'Axmod constant .................. ',F10.1/)
      WRITE (6,80100) KK,KMIN,KINC,KSKIP
80100 FORMAT (1H ,'Number of sections .............. ',I10/
     $        1H ,'First section to plot ........... ',I10/
     $        1H ,'Section increment ............... ',I10/
     $        1H ,'Sections per plot page .......... ',I10/)
      WRITE (6,80200) LABELS(1:NLB)
80200 FORMAT (1H ,'Plot title input ................ ',A//)
      IF (ITICK.EQ.1) THEN
         IF (NLB.GT.68) NLB = 68
         LABELS(NLB+1:NLB+13) = ' SECTION     '
         NLB = NLB+13
      END IF
CCC     READ(10) IDCOMP,IORGX,INCX,NG,NP,ISUM1,ISUM2,ISUM3,
CCC     $         IALPHA,IORD,            JSUM1,JSUM2,JSUM3,ZAPP

      mxtitl = 400
      do j=1,mxtitl
         title(j) = ' '
      end do

      STLEN=80
      CALL TRIML(MAPFILE,STLEN)
      CALL TRIMM(MAPFILE,STLEN)
      OPEN (UNIT=10,FILE=MAPFILE(1:STLEN),
     &     FORM='FORMATTED',STATUS='OLD')

C
C MODIFICATION 7/12/2006 ATB
      OPEN (UNIT=10,FILE=MAPFILE(1:STLEN),
     &     FORM='FORMATTED',STATUS='OLD')
      read(10,'(/I8)') ntitle
      read(10,'(A)')(title(j)(1:80), j=1,ntitle)
      read(10,'(9I8)') nxPts, nxStart, nxEnd, nyPts, nyStart,
     $     nyEnd, nzPts, nzStart, nzEnd
      read(10,'(6E12.5)') (cellDP(i), i=1,6)
      read(10,'(A)') zyx
C
C END MODIFICATION ATB

      zapp(3) = cellDP(1)
      zapp(4) = cellDP(2)
      zapp(5) = cellDP(3)
      zapp(6) = cellDP(4)
      zapp(7) = cellDP(5)
      zapp(8) = cellDP(6)

      iorgx(1) = nzStart
      iorgx(2) = nxStart
      iorgx(3) = nyStart
      incx(1) = 1
      incx(2) = 1
      incx(3) = 1
      ng(1) = nzPts
      ng(2) = nyPts
      ng(3) = nxPts
      np(1) = nzEnd - nzStart + 1
      np(2) = nyEnd - nyStart + 1
      np(3) = nxEnd - nxStart + 1
      isum1 = 3
      isum2 = 2
      isum3 = 1
      ialpha(1) = 'Z'
      ialpha(2) = 'Y'
      ialpha(3) = 'X'

      BETA=ZAPP(11-(ISUM2+ISUM3))
      BE90=(BETA-90)*ATAN(1.0)/45.0
      CB=COS(BE90)
      SB=SIN(BE90)
C     
C---- ASSUME UNIT STEPS AND A RECTANGULAR NET
      LBX = 16
      WRITE (LABELX(1:LBX),80300) IALPHA(3),NG(3)
      WRITE (LABELY(1:LBX),80300) IALPHA(2),NG(2)
80300 FORMAT (A4,' IN ',I4,' THS')
      IF(KK.EQ.0) KK=NP(1)
      M = NP(2)
      N = NP(3)
      lw = n
      XA = IORGX(3)
      YA = IORGX(2)
      XB = IORGX(3) + NP(3) - 1
      YB = IORGX(2) + NP(2) - 1
      SX = ZAPP(ISUM3+2) / FLOAT(NG(3)) * SCALE
      SY = ZAPP(ISUM2+2) / FLOAT(NG(2)) * SCALE
      XG = (XB-XA) * SX
      YG = (YB-YA) * SY
      XO = (YB-YA) * SB * SX
      XL = XG + ABS(XO)
      YL = YG * CB
      IF (XL.GT.PLTMAX .OR. YL.GT.PLTMAX) THEN
         PLOTSCALE = PLTMAX / MAX (XL,YL)
         XO = XO * PLOTSCALE
         XG = XG * PLOTSCALE
         YG = YG * PLOTSCALE
         XL = XL * PLOTSCALE
         YL = YL * PLOTSCALE
         WRITE (6,999) PLOTSCALE
 999     FORMAT (' Input scale reduced by a factor of',F8.5)
      END IF
C     CALL PLOTS (DUMMY,DUMMY,1)

C---- Set up output file name
      iFileIndx = 1
      STLEN=80
      CALL TRIML(PLOTFILE,STLEN)
      CALL TRIMM(PLOTFILE,STLEN)
      OPEN (UNIT=25,FILE=PLOTFILE(1:STLEN),
     &     FORM='FORMATTED',STATUS='UNKNOWN')
      write(25, '(a)') '%! This is a PostScript file!'
      write(25, *) '/in {72 mul} def'
      write(25, *) '/m {moveto} def'
      write(25, *) '/l {lineto} def'
      write(25, '(a,/,a,/a)') 
     @     '/Helvetica findfont', '15 scalefont', 'setfont'
C     
C---- SKIP SECTIONS UNTIL START SECTION REACHED
      DO 10 K = 1,KMIN-1
C     DO 10 I = 1,M
C     10	READ (10)
         read(10,*)               ! map section number
 10      read(10,*)               ! 2-D map section
C     
C---- PLOT SECTIONS
         XOFF=XO
         IF (XO.LT.0.0) XOFF=0.0
         KSKIP = KINC*KSKIP
              DO 50 K = 0,KK-1
C     DO 10000 I=1,M
C     READ(10) IXX1,IXX2,IXX3,LW,(RHO(I,J),J=1,LW)
C
C MODIFICATION 7/12/1006 ATB
                 read(10,'(I8)') ixx1
                 read(10,'(6E12.5)') ((rho(i,j), j=1,lw), i=1,m)
C END MODIFICATION
                 ixx2 = 0
                 ixx3 = 0
C     10000	CONTINUE
                 IF (MOD(K,KINC).NE.0) GOTO 50
                 KKPLT = MOD(K,KSKIP)
                 IF (KKPLT.EQ.0) THEN
C     IF (K.NE.0) CALL PLOT (0.,0.,-999)
                 if (k .ne. 0) then
                    write(25, '(a,/,a)') 'stroke', 'showpage'
                    close(25)
                    iFileIndx = iFileIndx + 1
                    if (iFileIndx .lt. 10) then
                       write(number(1:1),'(i1)') iFileIndx
 1010                  format(i1)
                       name=NUMBER(1:1)//PLOTFILE(1:STLEN)
                       open(unit=25, file=name,
     &                      form='formatted', status='unknown')
                    else
                       write(number(1:2),'(i2)') iFileIndx
 1011                  format(i2)
                       name=NUMBER(1:2)//PLOTFILE(1:STLEN)
                       open(unit=25, file=name,
     &                      form='formatted', status='unknown')
                    end if
                    write(25, '(a)') '%! This is a PostScript file!'
                    write(25, *) '/in {72 mul} def'
                    write(25, *) '/m {moveto} def'
                    write(25, *) '/l {lineto} def'
                    write(25, '(a,/,a,/,a)') 
     @            '/Helvetica findfont', '15 scalefont', 'setfont'
                 end if
C     CALL PLOT (1.0+XOFF,1.0,-3)
C     IF (LAND.EQ.0) CALL PORTR
C     CALL NEWPEN (LINWID)
                    IF (ITICK.EQ.1) WRITE (LABELS(NLB-3:NLB),'(I4)') IXX1
                    if (land .ne. 0) then
                       write(25, *) '8.5 in  0.0 in translate'
                       write(25, *) '90 rotate'
                    end if
                    write(25, '(f6.3,a,f6.3,a)') 
     @                   1.0+xoff, ' in ', 1.0, ' in translate'
                    write(25, '(f4.2,a)') contourLineWidth, 
     $                          ' setlinewidth'
                    XZR = -XO
                    YZR =  YL+0.5
                    HGH = 0.15
                    IF (ITICK.NE.0)
     $                 CALL PSTRNG (LABELS,NLB,XZR,YZR,HGH)
                    XZR =  0.00
                    YZR = -0.40
                    IF (ITICK.EQ.1) 
     $                 CALL PSTRNG (LABELX, LBX,XZR,YZR,HGH)
                    XZR = -0.25
                    YZR =  0.00
                    IF (ITICK.EQ.1)
     $                  CALL PSTRNG (LABELY,-LBX,XZR,YZR,HGH)
                    write(25, '(a,/,a,/a)') 
     @                   '/Courier findfont', '9 scalefont', 'setfont'
                    YZR =  YL+3.6
                    DO J=1,NTITLE
                       XZR = -XO
                       YZR =  YZR-0.12
                       HGH = 0.15
                       CALL PSTRNG (TITLE(J)(9:80),NLB,XZR,YZR,HGH)
                    END DO
                 END IF
                 RHOMAX = 0.      
                 DO 20 I=1,M
                    DO 20 J=1,LW
                       IF(RHO(I,J).GT.RHOMAX) RHOMAX = RHO(I,J)
 20                 CONTINUE
                    CL(1) = CTRMIN
                    DO 30 I=2,50
                       CL(I) = CL(I-1) + CTRINT
                       IF(CL(I).GT.RHOMAX) GO TO 40
 30                 CONTINUE
 40                 NCL = I
                    IF (NCL.GT.50) NCL=50
                    CALL CONTUR
     $           (M,N,ISP,XA,XB,YA,YB,XG,YG,NCL,CL,KKPLT,AXMOD,ITICK)
                    WRITE (6,81000) 
     $                K,KKPLT,M,N,IXX1,IXX2,IXX3,LW,XA,XB,YA,YB,XG,YG
81000               FORMAT (1H ,4I10/1H ,4I10/1H ,6F10.2)
 50              CONTINUE
C     CALL PLOT (0.,0.,999)
              write(25, *) 'showpage'
              close(25)
              close(10)	
              END
      SUBROUTINE CONTUR
     $     (M,N,ISP,XA,XB,YA,YB,XG,YG,NCL,CL,KKPLT,AXMOD,ITICK)
      real*8 Z(400,400)
C     COMMON /A/ Z(400,400)
      common /A/ z
      COMMON/INDEX/MROW,NCOL,MMROW,NNCOL
      COMMON/XYBNDS/XMIN,XMAX,YMIN,YMAX,XSIZE,YSIZE,
     1     HX,HY,XS,XSS,YS,YSS,FXA,FYA
      COMMON/CTRLVL/NLVLS,NLV,CLEVEL(50)
      COMMON/CAVIN/IDIM,DUM(10035)
      common / lines / contourLineWidth, borderLineWidth
      DIMENSION CL(50)
C     Z(I,J) IS THE ORDINATE AT POINT X(J), Y(I)
C     MXN IS THE SIZE OF THE CALCULATED X-Y GRID
C     MMXNN IS THE SIZE OF THE EXPANDED(BY INTERPOLATION) X-Y GRID
C     XA,XB,YA,YB ARE THE MINIMUM AND MAXIMUM VALUES
C     OF X AND Y,
C     XG IS THE WIDTH OF THE GRAPH IN inES.
C     YG IS THE HEIGHT OF THE GRAPH IN inES.
C     NCL IS THE NO. OF CONTOUR LEVELS
C     CL(I) ARE THE CONTOUR LEVELS
C     THE X(I) ARE ASSUMED TO BE EQUALLY SPACED, AND
C     LIKEWISE, THE Y(I).
C     FX IS THE FUNCTION TO BE PLOTTED ALONG THE X-AXIS.
C     FY IS THE FUNCTION TO BE PLOTTED ALONG THE Y-AXIS.
      MROW=M
      NCOL=N
      XMIN=XA
      XMAX=XB
      YMIN=YA
      YMAX=YB
      XSIZE=XG
      YSIZE=YG
      NLVLS=NCL
      DO 1 I=1,NLVLS
 1       CLEVEL (I)=CL(I)
 2       CALL SPLINE(ISP)
         HX=(XMAX-XMIN)/FLOAT(NCOL-1)
         HY=(YMAX-YMIN)/FLOAT(MROW-1)
         XS=(XMAX-XMIN)/FLOAT(NNCOL-1)
         YS=(YMAX-YMIN)/FLOAT(MMROW-1)
         FXA=XMIN
         FYA=YMIN
         XSS=XG/(XMAX-XMIN)
         YSS=YG/(YMAX-YMIN)
         DO3NLV=1,NLVLS
 3       CALL SCAN
         IF (KKPLT.EQ.0) then
	    write(25, '(f4.2,a)') borderLineWidth, ' setlinewidth'
	    CALL LABEL (AXMOD,ITICK)
	    write(25, '(f4.2,a)') contourLineWidth, ' setlinewidth'
         end if
         END
C     
      SUBROUTINE SCAN
C     AM IS THE MATRIX TO BE CONTOURED. MT AND NT ARE ITS X AND Y DIMENSIONS
C     CL(NLV) IS THE CONTOUR LEVEL.
C     THE N  (X,Y)  VALUES OF ONE CONTOUR LINE ARE PLOTTED WHEN
C     THEY ARE AVAILABLE.
C     COMMON /A/ AM(160000)
      real*8 am(160000)
      common /A/ am
      COMMON /CTRLVL/ NCL,NLV,CL(50)
      COMMON/INDEX/DUM(2),MT,NT
      COMMON/CAVIN/DIM, IX,IY,IDX,IDY,ISS,
     1     NP,N,CV,IS,IS0,IX0,IY0,DCP,
     2     INX(8),INY(8),REC(2000),X(4003),Y(4003)
      INTEGER REC,DIM
      DATA INX/-1,-1,0,1,1,1,0,-1/
      DATA INY/0,1,1,1,0,-1,-1,-1/
      DATA DIM /400/
      NP=0
      ISS=0
      CV=CL(NLV)
      MT1=MT-1
      NT1=NT-1
      DO 110 I=1,MT1
         IF(AM(I)  -CV)55,110,110
 55      IF(AM(I+1)  -CV)110,57,57
 57      IX0=I+1
         IX=I+1
         IY0=1
         IY=1
         IS0=1
         IS=1
         IDX=-1
         IDY=0
         CALL TRACE
 110  CONTINUE
      J=MT-DIM
      DO20I=1,NT1
      J=J+DIM
      IF(AM(J)-CV)15,20,20
 15   IF(AM(J+DIM)-CV)20,17,17
 17   IX0=MT
      IX=MT
      IY0=I+1
      IY=I+1
      IDX=0
      IDY=-1
      IS0=7
      IS=7
      CALL TRACE
 20   CONTINUE
      J=MT+NT1*DIM+1
      DO30I=1,MT1
      J=J-1
      IF(AM(J)-CV)25,30,30
 25   IF(AM(J-1)-CV)30,27,27
 27   IX0=MT-I
      IX=MT-I
      IY0=NT
      IY=NT
      IDX=1
      IDY=0
      IS0=5
      IS=5
      CALL TRACE
 30   CONTINUE
      J=NT*DIM+1
      DO40I=1,NT1
      J=J-DIM
      IF(AM(J)-CV)35,40,40
 35   IF(AM(J-DIM)-CV)40,37,37
 37   IX0=1
      IX=1
      IY0=NT-I
      IY=NT-I
      IDX=0
      IDY=1
      IS0=3
      IS=3
      CALL TRACE
 40   CONTINUE
      ISS=1
      L=0
      DO13J=2,NT1
      L=L+DIM
      DO10I=1,MT1
      L=L+1
      IF(AM(L)-CV)5,10,10
    5 IF(AM(L+1)-CV)10,7,7
    7 K=L+1
      DO 9 ID = 1,NP
         IF(REC(ID)-K)9,10,9
    9 CONTINUE
      IX0=I+1
      IX=I+1
      IY0=J
      IY=J
      IDX=-1
      IDY=0
      IS0=1
      IS=1
      CALL TRACE
 10   CONTINUE
 13   L=L-MT1
      END

      SUBROUTINE TRACE
      real*8 am(160000)
C     COMMON /A/ AM(160000)
      common /A/ am
      COMMON /INDEX/ DUM(2),MT,NT
      COMMON/XYBNDS/XA,XB,YA,YB,XSIZE,YSIZE,HX,HY,
     1     XS,XSS,YS,YSS,FXA,FYA
      COMMON/CAVIN/DIM, IX,IY,IDX,IDY,ISS,
     1     NP,N,CV,IS,IS0,IX0,IY0,DCP,
     2     INX(8),INY(8),REC(2000),X(4003),Y(4003)
      COMMON/CTRLVL/NCL,NLV,CL(50)
      COMMON /SKEW / SB,CB,BETA
      INTEGER REC,DIM
      N=0
      JY=DIM*(IY-1)+IX
      MY=DIM*IDY+IDX+JY
    2 N=N+1
      IF(N-4000) 3,3,32
    3 IF(IDX)5,4,6
    4 X(N)=FLOAT(IY-1)+FLOAT(IDY)*(AM(JY)-CV)/(AM(JY)-AM(DIM*IDY+JY))
      Y(N)=FLOAT(IX-1)
      GOTO7
    5 NP=NP+1
      REC(NP)=JY
    6 Y(N)=FLOAT(IX-1)+FLOAT(IDX)*(AM(JY)-CV)/(AM(JY)-AM(JY+IDX))
      X(N)=FLOAT(IY-1)
    7 IS=IS+1
    8 IF(IS-8)10,10,9
    9 IS=IS-8
 10   IDX=INX(IS)
      IDY=INY(IS)
      IX2=IX+IDX
      IY2=IY+IDY
      IR=IDX*IDY
 11   IF(ISS.EQ.0) GO TO 15
 13   IF(IS.NE.IS0.OR.IY.NE.IY0.OR.IX.NE.IX0) GO TO 16
 14   N=N+1
      X(N)=X(1)
      Y(N)=Y(1)
      GOTO73
 15   IF(IX2.EQ.0.OR.IX2.GT.MT.OR.IY2.EQ.0.OR.IY2.GT.NT) GO TO 73
 16   MY=DIM*IDY+IDX+JY
      IF(IR)19,17,20
 17   IF(CV-AM(MY))18,18,2
 18   IX=IX2
      IY=IY2
 81   IS=IS+5
      JY=MY
      GOTO8
 19   KY=JY+IDX
      LY=MY-IDX
      GOTO21
 20   KY=MY-IDY
      LY=JY+IDY
 21   DCP=(AM(JY)+AM(KY)+AM(LY)+AM(MY))*.25
      IF(CV-DCP)23,23,22
 22   CALL GETPT(JY)
      GOTO7
 23   IF(IR)24,25,25
 24   IX=IX2
      IDX=-IDX
      CALL GETPT(KY)
      IY=IY2
      IDY=-IDY
      GOTO26
 25   IY=IY2
      IDY=-IDY
      CALL GETPT(KY)
      IX=IX2
      IDX=-IDX
 26   IF(CV-AM(MY))81,81,28
 28   CALLGETPT(MY)
      IF(IR)29,30,30
 29   IX=IX+IDX
      IDX=-IDX
      GOTO31
 30   IY=IY+IDY
      IDY=-IDY
 31   IF(CV-AM(LY))33,33,34
 33   IS=IS-1
      JY=LY
      GOTO10
 34   CALLGETPT(LY)
      IF(IR)35,36,36
 35   IY=IY+IDY
      GOTO7
 36   IX=IX+IDX
      GOTO7
 32   PRINT103,CV
 73   DO74I=1,N
      X(I)=XSS*((X(I)*XS+XA)-FXA)
 74   Y(I)=YSS*((Y(I)*YS+YA)-FYA)
      DO I = 1,N
         X(I) = X(I) - Y(I) * SB
         Y(I) = Y(I) * CB
      END DO
C     CALL PLOT(X(1), Y(1), 3)
      write(25, *) 'newpath'
      write(25, '(f6.3,a,f6.3,a)') x(1), ' in ', y(1), ' in m'

      DO75I=2,N
C     75 CALL PLOT(X(I), Y(I), 2)
 75   write(25, '(f6.3,a,f6.3,a)') x(i), ' in ', y(i), ' in l'
      write(25, *) 'stroke'
      RETURN
 103  FORMAT(1H0,23HA CONTOUR LINE AT LEVEL,F10.5,
     1     41H WAS TERMINATED BECAUSE IT CONTAINED MORE,
     2     23H THAN 1600 PLOT POINTS.)
      END

      SUBROUTINE GET PT(J)
      real*8 am(160000)
C     COMMON /A/ AM(160000)
      common /A/ am
      COMMON/CAVIN/DIM, IX,IY,IDX,IDY,ISS,
     1     NP,N,CV,IS,IS0,IX0,IY0,DCP,
     2     INX(8),INY(8),REC(2000),X(4003),Y(4003)
      N=N+1
      B=AM(J)-DCP
      IF(B.NE.0) GO TO 2
    1 V=.5
      GOTO3
    2 V=.5*(AM(J)-CV)/B
    3 Y(N)=FLOAT(IX-1)+FLOAT(IDX)*V
      X(N)=FLOAT(IY-1)+FLOAT(IDY)*V
      END

      SUBROUTINE LABEL(AXMOD,ITICK)
      COMMON/XYBNDS/XA,XB,YA,YB,XSIZE,YSIZE,HX,HY,
     1     XS,XSS,YS,YSS,FXA,FYA
      COMMON/INDEX/M,N,MM,NN
      COMMON /SKEW / SB,CB,BETA
      common / lines / contourLineWidth, borderLineWidth
      DIMENSION XX(1001),YY(1001),NX(1001),NY(1001)
      DATA INIT,JX,JY /0,0,0/
      IF(INIT.EQ.1) GO TO 21
      INIT = 1
      XG = XSIZE
      YG = YSIZE
C---- SET UP X-AXIS LABELLING
      X=XA
      DO 10 I=1,N
         IF(AMOD(X,AXMOD).NE.0) GO TO 10
         JX=JX+1
         XX(JX) = XSS*(X-FXA)
         NX(JX)=X
 10      X=X+HX
C---- SET UP Y-AXIS LABELLING
         Y=YA
         DO 20 I=1,M
            IF(AMOD(Y,AXMOD).NE.0) GO TO 20
            JY=JY+1
            YY(JY) = YSS*(Y-FYA)
            NY(JY)=Y
 20         Y=Y+HY
C---- MAKE A BOUNDARY
 21         XO = -YG*SB
            XP=XG+XO
            YP=YG*CB
C     30 CALL PLOT(0.,0.,3)
C     CALL PLOT(XG,0.,2)
C     CALL PLOT(XP,YP,2)
C     CALL PLOT(XO,YP,2)
C     CALL PLOT(0.,0.,2)
C BEGIN MODIFICATION ( 0 -> 0.0) ATB 7/12/2006
            write(25, *) 'newpath'
            write(25, '(f6.3,a,f6.3,a)') 0.0, ' in ', 0.0, ' in m'
            write(25, '(f6.3,a,f6.3,a)') xg, ' in ', 0.0, ' in l'
            write(25, '(f6.3,a,f6.3,a)') xp, ' in ', yp, ' in l'
            write(25, '(f6.3,a,f6.3,a)') xo, ' in ', yp, ' in l'
            write(25, '(f6.3,a,f6.3,a)') 0.0, ' in ', 0.0, ' in l'
            write(25, *) 'stroke'
            IF (ITICK.EQ.0) RETURN
            write(25, *) 'newpath'
C---- TICK-MARK LOWER X-AXIS
            DO 35 I=1,JX
C     CALL PLOT(XX(I),0.,3)
C     CALL PLOT(XX(I),-0.1,2)
               write(25, '(f6.3,a,f6.3,a)') xx(i), ' in ', 0.0, ' in m'
               write(25, '(f6.3,a,f6.3,a)') xx(i), ' in ', -0.1, ' in l'
C END MODIFICATION
 35         CONTINUE
C---- TICK-MARK UPPER X-AXIS
            DO 45 I=1,JX
C     CALL PLOT(XX(I)+XO,YP,3)
C     CALL PLOT(XX(I)+XO,YP+0.1,2)
          write(25, '(f6.3,a,f6.3,a)') xx(i)+xo, ' in ', yp, ' in m'
          write(25, '(f6.3,a,f6.3,a)') xx(i)+xo, ' in ', yp+0.1, ' in l'
 45         CONTINUE
C---- TICK-MARK RIGHT Y-AXIS
            DO 40 I=1,JY
               XP=XG-YY(I)*SB
               YP=YY(I)*CB
C     CALL PLOT(XP,YP,3)
C     CALL PLOT(XP+0.1,YP,2)
               write(25, '(f6.3,a,f6.3,a)') xp, ' in ', yp, ' in m'
               write(25, '(f6.3,a,f6.3,a)') xp+0.1, ' in ', yp, ' in l'
 40         CONTINUE
C---- TICK-MARK LEFT Y-AXIS
            DO 50 I=1,JY
               XP=-YY(I)*SB
               YP= YY(I)*CB
C     CALL PLOT(XP,YP,3)
C     CALL PLOT(XP-0.1,YP,2)
               write(25, '(f6.3,a,f6.3,a)') xp, ' in ', yp, ' in m'
               write(25, '(f6.3,a,f6.3,a)') xp-0.1, ' in ', yp, ' in l'
               write(25, *) 'stroke'
 50         CONTINUE

            RETURN
            END
C     
      SUBROUTINE SPLINE(ISP)
      real*8 am(400,400)
C     COMMON /A/ AM(400,400)
      common /A/ am
      COMMON/INDEX/M,N,MM,NN
      COMMON /SPLIN/ Q(1001),Z(1001),R(1001)
      DATA INIT /0/
C---- INITIALIZE
      IF(INIT.EQ.1) GO TO 8
      INIT = 1
      GO TO 4
    2 ISP = ISP - 1
    4 MM = M + ISP*(M-1)
      NN = N + ISP*(N-1)
      IF(MM.GT.400.OR.NN.GT.400) GO TO 2
      ISPP = ISP + 1
      SIP = ISPP
      WRITE( 6,6) M,N,MM,NN,ISPP
 6    FORMAT(21H1THE ORIGINAL GRID OF,I3,3H BY,I3,21H HAS BEEN EXPANDED
     $TO,I3,3H BY,I3,44H USING SPLINE INTERPOLATION INCREMENTS OF 1/ ,I1
     $,18H THE ORIGINAL GRID)
C---- INTERPOLATE FOR ROWS
    8 IF(ISP.EQ.0) RETURN
      LN=N-1
      IM=M+1
      DO 30 IT=1,M
         I=IM-IT
         DO 10 J=1,N
            Q(J) = AM(I,J)
 10      CONTINUE
         CALL SPLICN(N)
         II = I+ISP*(I-1)
         DO 20 J=1,LN
            JJ = J+ISP*(J-1)
            AM(II,JJ) = Q(J)
            DO 20 K=1,ISP
               RAT = FLOAT(K)/SIP
               QAT = 1. - RAT
               JK=JJ+K
               AM(II,JK) = QAT*(Z(J)*QAT**2+R(J)) + 
     $                     RAT*(Z(J+1)*RAT**2+R(J+1))
 20         CONTINUE
            AM(II,NN) = Q(N)
 30      CONTINUE
C---- INTERPOLATE FOR COLUMNS
         LM=M-1
         DO 60 JJ=1,NN
            DO 40 I=1,M
               II = I+ISP*(I-1)
               Q(I) = AM(II,JJ)
 40         CONTINUE
            CALL SPLICN(M)
            DO 50 I=1,LM
               II = I+ISP*(I-1)
               DO 50 K=1,ISP
                  RAT = FLOAT(K)/SIP
                  QAT = 1. - RAT
                  IK=II+K
                  AM(IK,JJ) = QAT*(Z(I)*QAT**2+R(I)) + 
     $                        RAT*(Z(I+1)*RAT**2+R(I+1))
 50            CONTINUE
 60         CONTINUE
            RETURN
            END
C     
      SUBROUTINE SPLICN(M)
      COMMON /SPLIN/ Y(1001),Z(1001),R(1001)
      DIMENSION A(1001,3),B(1001),E(1001)
      EQUIVALENCE (E,Z)
      MM=M-1
      DO 2 K=1,MM
 2       E(K) = Y(K+1)-Y(K)
         DO 3 K=2,MM
 3          B(K)=E(K)-E(K-1)
            A(1,2)=-2.
            A(1,3)=1.
            A(2,2)=1.
            A(2,3)=0.
            DO 4 K=3,MM
               A(K,2)=2./3. - A(K-1,3)/6.
               B(K)=B(K)-B(K-1)/6.
               A(K,3)=1./(6.*A(K,2))
 4             B(K)=B(K)/A(K,2)
               A(M,1)=2. + A(M-2,3)
               A(M,2) = -A(M,1)*A(M-1,3) -1.
               B(M)=B(M-2)-A(M,1)*B(M-1)
               Z(M)=B(M)/A(M,2)
               MN=M-2
               DO 6 I=1,MN
                  K=M-I
 6                Z(K)=B(K)-A(K,3)*Z(K+1)
                  Z(1)=-A(1,2)*Z(2)-A(1,3)*Z(3)
                  DO 7 K=1,M
                     Z(K) = Z(K)/6.
                     R(K) = Y(K) - Z(K)
 7                CONTINUE
                  RETURN
                  END
      SUBROUTINE PSTRNG (LABEL,NL,X0,Y0,H)
C     ----------------------------------------
C     
C---- PLOT A CHARACTER STRING OF NL CHARCATERS
C     
      CHARACTER LABEL*(*)
      COMMON /LIMIT/ XXMAX
      COMMON /SKEW / SB,CB,BETA
      NLL = ABS(NL)
      ANGLE=0.0
      IF (NL.LT.0) ANGLE=BETA
      IF (H*NLL.GT.XXMAX) H = XXMAX/NLL
C     CALL FONT(2)
C     CALL SYMBOL
C     1 (X0,Y0,H,LABEL,ANGLE,NLL)
      write(25, '(f6.3,a,f6.3,a)') x0, ' in ', y0, ' in m'
      write(25, '(f7.1,a)') angle, ' rotate'
      write(25, '(a,a80,a)') '(', label,  ') show'
      write(25, '(f7.1,a)') -angle, ' rotate'
      RETURN
      END
C     
      SUBROUTINE TRIMM(ST,STLEN)
C     routine removes trailing blanks in string ST
C     
C     Axel T. Brunger
C     ===============
C     
      IMPLICIT NONE
C     input/output
      CHARACTER*(*) ST
      INTEGER STLEN
C     begin
      IF (STLEN.GT.0) THEN
         DO WHILE (STLEN.GE.1.AND.ST(STLEN:STLEN).EQ.' ')
            STLEN=STLEN-1
         END DO
      END IF
      RETURN
      END
C     
      SUBROUTINE TRIML(ST,STLEN)
C     routine removes leading blanks in string ST and fills the
C     remainder on the right site with blanks
C     
C     Axel T. Brunger
C     ===============
C     
      IMPLICIT NONE
C     input/output
      CHARACTER*(*) ST
      INTEGER STLEN
C     local
      INTEGER IS, I, STOLD
      CHARACTER*1 TEMP
C     begin
      IF (STLEN.GT.0) THEN
         STOLD=STLEN
         IS=1
         DO WHILE (IS.LE.STLEN.AND.ST(IS:IS).EQ.' ')
            IS=IS+1
         END DO
         STLEN=STLEN-IS+1
         IF (STLEN.LT.STOLD) THEN
            IS=IS-1
C     
C     The following code avoids a self-assignment ( ST(..) = ST(..) )
C     and therefore makes the CFT 1.13 compiler on the CRAY happy.
C     As soon as this is fixed (in the compiler) this inefficient loop
C     should be replaced !
            DO I=1,STLEN
               TEMP=ST(I+IS:I+IS)
               ST(I:I)=TEMP
            END DO
            ST(STLEN+1:STOLD)=' '
         END IF
      END IF
      RETURN
      END
C
