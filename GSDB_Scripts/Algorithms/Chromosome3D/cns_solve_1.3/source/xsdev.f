! OpenMP version Kay Diederichs 2006
C======================================================================
      SUBROUTINE XMDOSDEV(VLEVEL,VMAX,VSTACK,N,INDEXA,INDEXB,INDEXC,
     &                    MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &                    NA,NB,NC,MASK,ARGS)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER NAASY, NBASY, NCASY, MAASY, MBASY, MCASY
      INTEGER NA, NB, NC
      INTEGER MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER INDEXA(*), INDEXB(*), INDEXC(*)
C parameters
      DOUBLE PRECISION RDEF, MDEF
      PARAMETER (RDEF=3.0D0, MDEF=0.0D0)
C local
      DOUBLE PRECISION RADIUS, MEAN
      INTEGER NAMIN, NAMAX, NBMIN, NBMAX, NCMIN, NCMAX
      INTEGER XAMIN, XAMAX, XBMIN, XBMAX, XCMIN, XCMAX
C pointers
      INTEGER SMASK, EXMAP
C begin
C
C get radius if an argument
      IF (ABS(DBLE(ARGS(1))).LT.RSMALL) THEN
         RADIUS=RDEF
      ELSE
         RADIUS=ABS(DBLE(ARGS(1)))
      END IF
C
C get mean if an argument
      IF (ABS(DBLE(ARGS(2))).LT.RSMALL) THEN
         MEAN=MDEF
      ELSE
         MEAN=ABS(DBLE(ARGS(2)))
      END IF
C
C calculate size of sphere mask and allocate space
      CALL SPHFRAC(RADIUS,NAMIN,NAMAX,NBMIN,NBMAX,NCMIN,NCMAX)
      SMASK=ALLHP(INTEG4(2*(NBMAX-NBMIN+1)*
     &                     (NCMAX-NCMIN+1)))
C
C calculate space needed for expanded map and allocate space
      XAMAX=NAASY+NAMAX
      XAMIN=MAASY-ABS(NAMIN)
      XBMAX=NBASY+NBMAX
      XBMIN=MBASY-ABS(NBMIN)
      XCMAX=NCASY+NCMAX
      XCMIN=MCASY-ABS(NCMIN)
      EXMAP=ALLHP(IREAL4((XAMAX-XAMIN+1)*
     &                   (XBMAX-XBMIN+1)*
     &                   (XCMAX-XCMIN+1)))
C
C fill expanded map
      CALL XMSDEVEX(VLEVEL,VMAX,VSTACK,N,INDEXA,INDEXB,INDEXC,MASK,
     &              HEAP(EXMAP),XAMIN,XAMAX,XBMIN,XBMAX,XCMIN,XCMAX)
C
C calculate standard deviation of points within radius at each point in map
      CALL XMDOSDEV2(VLEVEL,VMAX,VSTACK,N,INDEXA,INDEXB,INDEXC,MASK,
     &               HEAP(EXMAP),XAMIN,XAMAX,XBMIN,XBMAX,XCMIN,XCMAX,
     &               HEAP(SMASK),NAMIN,NAMAX,NBMIN,NBMAX,NCMIN,NCMAX,
     &               RADIUS,MEAN)
C
C free up space for expanded map and sphere mask
      CALL FREHP(EXMAP,IREAL4((XAMAX-XAMIN+1)*
     &                        (XBMAX-XBMIN+1)*
     &                        (XCMAX-XCMIN+1)))
      CALL FREHP(SMASK,INTEG4(2*(NBMAX-NBMIN+1)*
     &                          (NCMAX-NCMIN+1)))
C
      RETURN
      END
C======================================================================
      SUBROUTINE SPHFRAC(RADIUS,NAMIN,NAMAX,NBMIN,NBMAX,NCMIN,NCMAX)
C
C calculate grid point limits that cover the sphere of given radius
C return NAMIN, NAMAX, NBMIN, NBMAX, NCMIN, NCMAX
C
C Author: Axel T. Brunger and Paul Adams
C ======================================
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
C I/O
      DOUBLE PRECISION RADIUS
      INTEGER NAMIN, NAMAX, NBMIN, NBMAX, NCMIN, NCMAX
C local
      INTEGER IAT
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ
      DOUBLE PRECISION XL(8), YL(8), ZL(8)
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C set up fake coordinate list of cube edges
      XL(1)=-RADIUS
      YL(1)=-RADIUS
      ZL(1)=-RADIUS
      XL(2)= RADIUS
      YL(2)=-RADIUS
      ZL(2)=-RADIUS
      XL(3)=-RADIUS
      YL(3)= RADIUS
      ZL(3)=-RADIUS
      XL(4)= RADIUS
      YL(4)= RADIUS
      ZL(4)=-RADIUS
      XL(5)=-RADIUS
      YL(5)=-RADIUS
      ZL(5)= RADIUS
      XL(6)= RADIUS
      YL(6)=-RADIUS
      ZL(6)= RADIUS
      XL(7)=-RADIUS
      YL(7)= RADIUS
      ZL(7)= RADIUS
      XL(8)= RADIUS
      YL(8)= RADIUS
      ZL(8)= RADIUS
C
C fractionalize cube edge coordinates
      DO IAT=1,8
         TEMPX=XL(IAT)
         TEMPY=YL(IAT)
         TEMPZ=ZL(IAT)
         XL(IAT)=XRTR(1,1)*TEMPX +XRTR(1,2)*TEMPY +XRTR(1,3)*TEMPZ
         YL(IAT)=XRTR(2,1)*TEMPX +XRTR(2,2)*TEMPY +XRTR(2,3)*TEMPZ
         ZL(IAT)=XRTR(3,1)*TEMPX +XRTR(3,2)*TEMPY +XRTR(3,3)*TEMPZ
      END DO
C
      DO IAT=1,8
         XL(IAT) = INT(XL(IAT)*NA) + INT(SIGN(ONE,XL(IAT)))
         YL(IAT) = INT(YL(IAT)*NB) + INT(SIGN(ONE,YL(IAT)))
         ZL(IAT) = INT(ZL(IAT)*NC) + INT(SIGN(ONE,ZL(IAT)))
      END DO
C
C determine min and max
      NAMIN=XL(1)
      NAMAX=XL(1)
      NBMIN=YL(1)
      NBMAX=YL(1)
      NCMIN=ZL(1)
      NCMAX=ZL(1)
      DO IAT=2,8
         NAMIN=MIN(NAMIN,NINT(XL(IAT)))
         NAMAX=MAX(NAMAX,NINT(XL(IAT)))
         NBMIN=MIN(NBMIN,NINT(YL(IAT)))
         NBMAX=MAX(NBMAX,NINT(YL(IAT)))
         NCMIN=MIN(NCMIN,NINT(ZL(IAT)))
         NCMAX=MAX(NCMAX,NINT(ZL(IAT)))
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMSDEVEX(VLEVEL,VMAX,VSTACK,N,
     &                    INDEXA,INDEXB,INDEXC,MASK,
     &                    EXMAP,XAMIN,XAMAX,XBMIN,XBMAX,XCMIN,XCMAX)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XAMIN, XAMAX, XBMIN, XBMAX, XCMIN, XCMAX
      REAL EXMAP(XAMIN:XAMAX,XBMIN:XBMAX,XCMIN:XCMAX)
      INTEGER INDEXA(*), INDEXB(*), INDEXC(*)
C local
      INTEGER A, B, C, AA, BB, CC, AS, BS, CS, SYM, I
      LOGICAL COND
C begin
C
C fill expanded map with blank tags
C
      DO C=XCMIN,XCMAX
      DO B=XBMIN,XBMAX
      DO A=XAMIN,XAMAX
         EXMAP(A,B,C)=-9999
      END DO
      END DO
      END DO
C
C fill expanded map all map points passed from original map
C
      DO I=1,N
         EXMAP(INDEXA(I),INDEXB(I),INDEXC(I))=DBLE(VSTACK(I,VLEVEL))
      END DO
C
C go through all points in expanded map.
C if outside the ASU defined by mask then try to fill by symmetry
C
!$omp parallel do default(none)
!$omp& private(c,b,a,cond,sym,aa,bb,cc,as,bs,cs)
!$omp& shared(XCMIN,XCMAX,XBMIN,XBMAX,XAMIN,XAMAX,MAASY,NAASY,MBASY)
!$omp& shared(NBASY,MCASY,NCASY,mask,XRNSYM,XRSYMM,na,nb,nc)
!$omp& shared(exmap)
      DO C=XCMIN,XCMAX
      DO B=XBMIN,XBMAX
      DO A=XAMIN,XAMAX
      COND=.FALSE.
      IF ((A.GE.MAASY.AND.A.LE.NAASY).AND.
     &    (B.GE.MBASY.AND.B.LE.NBASY).AND.
     &    (C.GE.MCASY.AND.C.LE.NCASY)) THEN
      IF (MASK(A,B,C).LE.0) THEN
      COND=.TRUE.
      END IF
      ELSE
      COND=.TRUE.
      END IF
C
      IF (COND) THEN
      SYM=1
      COND=.FALSE.
      DO WHILE (SYM.LE.XRNSYM.AND..NOT.COND)
C
      AA=XRSYMM(SYM,1,1)*A+XRSYMM(SYM,1,2)*B+XRSYMM(SYM,1,3)*C
     &        +(XRSYMM(SYM,1,4)*NA)/XRSYTH
      BB=XRSYMM(SYM,2,1)*A+XRSYMM(SYM,2,2)*B+XRSYMM(SYM,2,3)*C
     &        +(XRSYMM(SYM,2,4)*NB)/XRSYTH
      CC=XRSYMM(SYM,3,1)*A+XRSYMM(SYM,3,2)*B+XRSYMM(SYM,3,3)*C
     &        +(XRSYMM(SYM,3,4)*NC)/XRSYTH
C
      AS=MOD(AA+10000*NA-MAASY,NA)+MAASY
      BS=MOD(BB+10000*NB-MBASY,NB)+MBASY
      CS=MOD(CC+10000*NC-MCASY,NC)+MCASY
      COND=(AS.GE.MAASY.AND.AS.LE.NAASY.AND.
     &      BS.GE.MBASY.AND.BS.LE.NBASY.AND.
     &      CS.GE.MCASY.AND.CS.LE.NCASY)
      IF (COND) COND=MASK(AS,BS,CS).GT.0
      SYM=SYM+1
      END DO
C
      IF (COND) THEN
      EXMAP(A,B,C)=EXMAP(AS,BS,CS)
      ELSE
      CALL WRNDIE(-5,'XMDOSDEV2','point could not be remapped')
      END IF
C
      END IF
C
      END DO
      END DO
      END DO
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XMDOSDEV2(VLEVEL,VMAX,VSTACK,N,
     &                     INDEXA,INDEXB,INDEXC,MASK,
     &                     EXMAP,XAMIN,XAMAX,XBMIN,XBMAX,XCMIN,XCMAX,
     &                     SMASK,NAMIN,NAMAX,NBMIN,NBMAX,NCMIN,NCMAX,
     &                     RADIUS,MEAN)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C I/O
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER MASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XAMIN, XAMAX, XBMIN, XBMAX, XCMIN, XCMAX
      REAL EXMAP(XAMIN:XAMAX,XBMIN:XBMAX,XCMIN:XCMAX)
      INTEGER NAMIN,NAMAX,NBMIN,NBMAX,NCMIN,NCMAX
      INTEGER SMASK(2,NBMIN:NBMAX,NCMIN:NCMAX)
      INTEGER INDEXA(*), INDEXB(*), INDEXC(*)
      DOUBLE PRECISION RADIUS, MEAN
C local
      INTEGER A, B, C, I, NUM_POINT, END, BEGIN
      INTEGER INDA, INDB, INDC
      DOUBLE PRECISION ATEMP, BTEMP, CTEMP, XDIST, YDIST, ZDIST
      DOUBLE PRECISION RADIUS2, LENGTH, SDEV
      LOGICAL START
C begin
      RADIUS2=RADIUS*RADIUS
C
C first fill sphere mask
      DO C=NCMIN,NCMAX
      CTEMP=C/DBLE(NC)
      DO B=NBMIN,NBMAX
      BTEMP=B/DBLE(NB)
      START=.FALSE.
      SMASK(1,B,C)=-9999
      SMASK(2,B,C)=-9999
      DO A=NAMIN,NAMAX
      ATEMP=A/DBLE(NA)
      XDIST=XRINTR(1,1)*ATEMP +XRINTR(1,2)*BTEMP +XRINTR(1,3)*CTEMP
      YDIST=XRINTR(2,1)*ATEMP +XRINTR(2,2)*BTEMP +XRINTR(2,3)*CTEMP
      ZDIST=XRINTR(3,1)*ATEMP +XRINTR(3,2)*BTEMP +XRINTR(3,3)*CTEMP
      LENGTH=(XDIST*XDIST)+(YDIST*YDIST)+(ZDIST*ZDIST)
      IF (LENGTH.LE.RADIUS2) THEN
         IF (.NOT.START) THEN
            SMASK(1,B,C)=A
            START=.TRUE.
            SMASK(2,B,C)=A
         ELSE
            SMASK(2,B,C)=A
         ENDIF
      ENDIF
      END DO
      END DO
      END DO
C
C loop over all selected points in map
!$omp parallel do default(none)
!$omp& private(i,sdev,num_point,c,indc,b,indb,begin,end,a,inda)
!$omp& shared(n,ncmin,ncmax,indexc,nbmin,nbmax,indexb,smask,indexa)
!$omp& shared(exmap,mean,vstack,vlevel)
      DO I=1,N
C
      SDEV=ZERO
      NUM_POINT=0
C
C accumulate standard deviation from mean (passed from outside)
      DO C=NCMIN,NCMAX
      INDC=INDEXC(I)+C
      DO B=NBMIN,NBMAX
      INDB=INDEXB(I)+B
      BEGIN=SMASK(1,B,C)
      END=SMASK(2,B,C)
      IF (BEGIN.NE.-9999) THEN
         DO A=BEGIN,END
            INDA=INDEXA(I)+A
            IF (EXMAP(INDA,INDB,INDC).NE.-9999) THEN
               SDEV=SDEV+(EXMAP(INDA,INDB,INDC)-MEAN)**2
               NUM_POINT=NUM_POINT+1
            END IF
         END DO
      END IF
      END DO
      END DO
C
C calculate standard deviation (RMS)
      SDEV=SQRT(SDEV/(NUM_POINT-1))
C
C put standard deviation in original map point
      VSTACK(I,VLEVEL)=DCMPLX(SDEV,ZERO)
C
      END DO
C
      RETURN
      END
C======================================================================
