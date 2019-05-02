      SUBROUTINE MODES
C
C Parser routine for MODES
C
C Author: Luke Rice and Axel T. Brunger
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
C
      INTEGER FLAGS,POINTR,HPLIST
C
C begin
      HPLIST=ALLHP(INTEG4(NATOM))
      FLAGS=ALLHP(INTEG4(NATOM))
      POINTR=ALLHP(INTEG4(MAX(NATOM,1)+1))
      CALL MODES2(HEAP(FLAGS),HEAP(HPLIST),HEAP(POINTR))
      CALL FREHP(POINTR,INTEG4(MAX(NATOM,1)+1))
      CALL FREHP(FLAGS,INTEG4(NATOM))
      CALL FREHP(HPLIST,INTEG4(NATOM))
      RETURN
      END
C
      SUBROUTINE MODES2(FLAGS,LIST,POINTR)
C
C For syntax see main parsing loop "HELP"
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'coordc.inc'
      INTEGER FLAGS(*),LIST(*),POINTR(*)
C local
      INTEGER MXLIST, FLIST, XL, YL, ZL
      INTEGER XSUM, YSUM, ZSUM, i, nflags
      INTEGER NGROUP, EUNIT, FUNIT, cunit, sigma, nmod(25)
      INTEGER BEGIN, SKIP, STOP,NUNIT,ii,munit
      double precision t
      LOGICAL QFORM, qprint,defaul
      CHARACTER*4 frunit
      CHARACTER*6 output
      PARAMETER (MXLIST=20)
C dynamic
      FLIST=ALLHP(INTEG4(MXLIST))
C
C defaults
      OFILE='OUTPUT'
C
C Make a default selection consisting of all atoms
      defaul=.true.
      NGROUP=1
      POINTR(1)=0
      POINTR(2)=0
      DO I=1,NATOM
      POINTR(2)=POINTR(2)+1
      LIST(POINTR(2))=I
      END DO
C
      qprint=.false.
      EUNIT = 0
      FUNIT = 6
      MUNIT = 0
      CUNIT = 0
      nmod(1)=10
      t=298.
      frunit='AKMA'
      CALL DYNSET('INIT',USED,BEGIN,SKIP,STOP,NUNIT,
     1                  MXLIST,HEAP(FLIST),QFORM)
C
C parsing
      CALL PUSEND('MODES>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('MODES>')
      CALL DYNSET('PARSE',USED,BEGIN,SKIP,STOP,NUNIT,
     1                  MXLIST,HEAP(FLIST),QFORM)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-modes')
C
C==================================================================
      ELSE IF (WD(1:4).EQ.'EOUT') THEN
      CALL NEXTFI('EOUTPut=',OFILE)
      CALL ASSFIL(OFILE,EUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'FOUT') THEN
      CALL NEXTFI('FOUTput=',OFILE)
      CALL ASSFIL(OFILE,FUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'MOUT') THEN
      CALL NEXTFI('MOUTput=',OFILE)
      CALL ASSFIL(OFILE,MUNIT,'WRITE','FORMATTED',ERROR)
C==================================================================
      else if (wd(1:4).eq.'COUT') then
      call nextfi('COUTput=',OFILE)
      call assfil(ofile,cunit,'WRITE','FORMATTED',ERROR)
C==================================================================
      ELSE IF (WD(1:4).EQ.'GROU') THEN
C
      if (defaul) then
      NGROUP=0
      defaul=.false.
      end if
C
      CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
C
C Make sure the selections are disjoint
      DO I=1,POINTR(NGROUP+1)
      IF (FLAGS(LIST(I)).EQ.1) THEN
      WRITE(6,'(2A,1X,A,1X,A,1X,A)')
     &' MODES-ERR: selection overlap for atom ',
     &    SEGID(LIST(I)),RESID(LIST(I)),RES(LIST(I)),TYPE(LIST(I))
      WRITE(6,'(A)')
     &' The GROUp statements select at least one atom to be ',
     &' in two or more groups.  Please check your GROUp ',
     &' selections with respect to overlapping definitions.'
      CALL WRNDIE(-1,'MODES',
     &'GROUp definition overlap --> change GROUp selections')
      FLAGS(LIST(I))=0
      END IF
      END DO
C
      CALL MAKIND(FLAGS,NATOM,NFLAGS)
C
      IF (NFLAGS.GT.0) THEN
         IF (NGROUP.GE.NATOM) THEN
         CALL WRNDIE(-5,'MODES','NGROUP exceeds NATOM')
         ELSE
         NGROUP=NGROUP+1
         END IF
C
C Now get pointers
         II=POINTR(NGROUP)
         DO I=1,NFLAGS
         II=II+1
         LIST(II)=FLAGS(I)
         END DO
         POINTR(NGROUP+1)=II
      END IF
C==================================================================
      else if (wd(1:4).eq.'QPRI') then
      call nextlo('QPRInt=',qprint)
C==================================================================
      else if (wd(1:4).eq.'NMOD') then
      call nexti('NMODes=',nmod(ngroup))
C==================================================================
      ELSE if (wd(1:4).eq.'TEMP') then
      call nextf('TEMPerature=',T)
C==================================================================
      else if (wd(1:4).eq.'FQUN') then
      call nextst('FQUNit=',frunit)
C==================================================================
      else if (wd(1:4).eq.'?   ') then
         if (qprint) then
            output='TRUE'
         else
            output='FALSE'
         end if
      WRITE(6,9000)
     @ 'MODES: Reference Temperature=', T,
     @ 'MODES: Frequency Output Units=', FRUNIT
9000  FORMAT(A,F6.2,/,
     @       A,A)
C==================================================================
      else
      CALL CHKEND('MODES>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
C
C
      IF (.NOT.ERROR) THEN
C
C dynamic
      XSUM=ALLHP(IREAL8(NATOM))
      YSUM=ALLHP(IREAL8(NATOM))
      ZSUM=ALLHP(IREAL8(NATOM))
      XL=ALLHP(IREAL8(NATOM))
      YL=ALLHP(IREAL8(NATOM))
      ZL=ALLHP(IREAL8(NATOM))
      SIGMA=ALLHP(IREAL8(9*Natom**2))
C
      CALL MODES3(EUNIT,FUNIT,cunit,munit,NATOM,HEAP(XSUM),
     1            HEAP(YSUM),HEAP(ZSUM),MXLIST,
     2            HEAP(FLIST),BEGIN,SKIP,STOP,NUNIT,
     3            HEAP(XL),HEAP(YL),HEAP(ZL),HEAP(SIGMA),QFORM,T,
     4            nmod,
     5            LIST,POINTR,NGROUP,frunit,amass)
C
      CALL FREHP(SIGMA,IREAL8(9*Natom**2))
      CALL FREHP(ZL,IREAL8(NATOM))
      CALL FREHP(YL,IREAL8(NATOM))
      CALL FREHP(XL,IREAL8(NATOM))
      CALL FREHP(ZSUM,IREAL8(NATOM))
      CALL FREHP(YSUM,IREAL8(NATOM))
      CALL FREHP(XSUM,IREAL8(NATOM))
      END IF
C
      CALL FREHP(FLIST, INTEG4(MXLIST))
C
      RETURN
      END
C
      SUBROUTINE MODES3(EUNIT,FUNIT,cunit,munit,NATOM,XSUM,
     1                  YSUM,ZSUM,MXLIST,FLIST,BEGIN,
     2                  SKIP,STOP,NUNIT,XL,YL,ZL,SIGMA,QFORM,
     3                  T,nmod,LIST,
     4                  POINTR,NGROUP,frunit,amass)
C
C computes covariances of the spatial atom displacements of
C a dynamics trajectory for selected pairs of atoms:
C
C      ab         a      a    b     b
C sigma   =  E[ (R  - E[R ])(R - E[R ]) ]
C      jk         j      j    k     k
C
C and computes quasi-normal modes.
C
C Author: Luke Rice
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER EUNIT,FUNIT,cunit,munit
      INTEGER NATOM
      DOUBLE PRECISION XSUM(*),YSUM(*),ZSUM(*),T
      INTEGER MXLIST, FLIST(*), BEGIN, SKIP, STOP, NUNIT
      INTEGER LIST(*),POINTR(*),NGROUP, nmod(*)
      DOUBLE PRECISION XL(*), YL(*), ZL(*),amass(*)
      DOUBLE PRECISION SIGMA(3*Natom,3*Natom)
      LOGICAL QFORM
      character*4 frunit
C local
      LOGICAL START
      INTEGER NCOORD, I, I1, I2,pointer,ii, NATMS
      INTEGER  L, M, EVECT, EVAL, GPCOV
      INTEGER ISTEP, IUNIT,ii1,ii2
      DOUBLE PRECISION DELTA,two
      double precision c,zero,one
      CHARACTER*4 HDR
      LOGICAL TDONE, TEOF, ERROR
      parameter (two=2.0d0,c=2.9979D10,zero=0.0d0,one=1.0d0)
C
C initialize
      NCOORD=0
      START=.TRUE.
      DO I=1,NGROUP
         DO POINTER=POINTR(I)+1,POINTR(I+1)
            II=LIST(POINTER)
            XSUM(II)=0.0D0
            YSUM(II)=0.0D0
            ZSUM(II)=0.0D0
         END DO
C
C initialize covariance matrix
         DO I1=POINTR(I)+1,POINTR(I+1)
            II1=LIST(I1)
            DO I2=POINTR(I)+1,POINTR(I+1)
               II2=LIST(I2)
               SIGMA(3*II1-2,3*II2-2)=ZERO
               SIGMA(3*II1-2,3*II2-1)=ZERO
               SIGMA(3*II1-2,3*II2)=ZERO
               SIGMA(3*II1-1,3*II2-2)=ZERO
               SIGMA(3*II1-1,3*II2-1)=ZERO
               SIGMA(3*II1-1,3*II2)=ZERO
               SIGMA(3*II1,3*II2-2)=ZERO
               SIGMA(3*II1,3*II2-1)=ZERO
               SIGMA(3*II1,3*II2)=ZERO
            END DO
         END DO
      END DO
C
C begin analysis - read trajectory frames
      IF (NUNIT.EQ.0) THEN
         CALL WRNDIE(-5,'MODES','PLEASE SPECIFY A TRAJECTORY FILE')
      END IF
C
      TDONE = .FALSE.
      TEOF = .FALSE.
      ERROR = .FALSE.
      DO WHILE (.NOT. (TDONE.OR.TEOF.OR.ERROR))
         CALL READDC(NATOM,XL,YL,ZL,
     1           START,TDONE,TEOF,ERROR,
     2           NUNIT,FLIST,IUNIT,BEGIN,STOP,SKIP,
     3           ISTEP,DELTA,HDR,QFORM)
         NCOORD=NCOORD+1
C
C accumulate average
         DO I=1,NGROUP
            DO POINTER=POINTR(I)+1,POINTR(I+1)
               II=LIST(POINTER)
               XSUM(II)=XSUM(II)+XL(II)
               YSUM(II)=YSUM(II)+YL(II)
               ZSUM(II)=ZSUM(II)+ZL(II)
            END DO
C
C now add in the contribution for this coordinate set
            DO I1=POINTR(I)+1,POINTR(I+1)
               II1=LIST(I1)
               DO I2=POINTR(I)+1,POINTR(I+1)
                  II2=LIST(I2)
C                                        sigmaXiXj
                  sigma(3*II1-2,3*II2-2)=sigma(3*II1-2,3*II2-2)+
     @                     xl(ii1)*xl(ii2)
C                                        sigmaXiYj
                  sigma(3*II1-2,3*II2-1)=sigma(3*II1-2,3*II2-1)+
     @                     xl(ii1)*yl(ii2)
C                                        sigmaXiZj
                  sigma(3*II1-2,3*II2)=sigma(3*II1-2,3*II2)+
     @                     xl(ii1)*zl(ii2)
C                                        sigmaYiXj
                  sigma(3*II1-1,3*II2-2)=sigma(3*II1-1,3*II2-2)+
     @                     yl(ii1)*xl(ii2)
C                                        sigmaYiYj
                  sigma(3*II1-1,3*II2-1)=sigma(3*II1-1,3*II2-1)+
     @                     yl(ii1)*yl(ii2)
C                                        sigmaYiZj
                  sigma(3*II1-1,3*II2)=sigma(3*II1-1,3*II2)+
     @                     yl(ii1)*zl(ii2)
C                                        sigmaZiXj
                  sigma(3*II1,3*II2-2)=sigma(3*II1,3*II2-2)+
     @                     zl(ii1)*xl(ii2)
C                                        sigmaZiYj
                  sigma(3*II1,3*II2-1)=sigma(3*II1,3*II2-1)+
     @                     zl(ii1)*yl(ii2)
C                                        sigmaZiZj
                  sigma(3*II1,3*II2)=sigma(3*II1,3*II2)+
     @                     zl(ii1)*zl(ii2)
               END DO
            END DO
         END DO
      END DO
C collect averages
      IF ((TDONE).AND.(NCOORD.GT.0)) THEN
         DO I=1,NGROUP
            DO POINTER=POINTR(I)+1,POINTR(I+1)
               II=LIST(POINTER)
               XSUM(II)=XSUM(II)/NCOORD
               YSUM(II)=YSUM(II)/NCOORD
               ZSUM(II)=ZSUM(II)/NCOORD
            END DO
C
C complete the covariance matrix
            DO I1=POINTR(I)+1,POINTR(I+1)
               II1=LIST(I1)
               DO I2=POINTR(I)+1,POINTR(I+1)
                  II2=LIST(I2)
                  sigma(3*II1-2,3*II2-2)=sigma(3*II1-2,3*II2-2)/NCOORD-
     @                 xsum(ii1)*xsum(ii2)
                  sigma(3*II1-2,3*II2-1)=sigma(3*II1-2,3*II2-1)/NCOORD-
     @                 xsum(ii1)*ysum(ii2)
                  sigma(3*II1-2,3*II2)=sigma(3*II1-2,3*II2)/NCOORD-
     @                 xsum(ii1)*zsum(ii2)
                  sigma(3*II1-1,3*II2-2)=sigma(3*II1-1,3*II2-2)/NCOORD-
     @                 ysum(ii1)*xsum(ii2)
                  sigma(3*II1-1,3*II2-1)=sigma(3*II1-1,3*II2-1)/NCOORD-
     @                 ysum(ii1)*ysum(ii2)
                  sigma(3*II1-1,3*II2)=sigma(3*II1-1,3*II2)/NCOORD-
     @                 ysum(ii1)*zsum(ii2)
                  sigma(3*II1,3*II2-2)=sigma(3*II1,3*II2-2)/NCOORD-
     @                 zsum(ii1)*xsum(ii2)
                  sigma(3*II1,3*II2-1)=sigma(3*II1,3*II2-1)/NCOORD-
     @                 zsum(ii1)*ysum(ii2)
                  sigma(3*II1,3*II2)=sigma(3*II1,3*II2)/NCOORD-
     @                 zsum(ii1)*zsum(ii2)
               END DO
            END DO
C
C now get covariances for each group
C
            NATMS = POINTR(I+1) - POINTR(I)
            GPCOV = ALLHP(IREAL8(9*NATMS**2))
            EVECT = ALLHP(IREAL8(9*NATMS**2))
            EVAL  = ALLHP(IREAL8(3*NATMS))
            L     = ALLHP(IREAL8(3*NATMS))
            M     = ALLHP(IREAL8(3*NATMS))
            CALL GPMODE(LIST,POINTR,I,NATMS,NATOM,AMASS,HEAP(GPCOV),
     @                  SIGMA,NMOD,HEAP(EVECT),HEAP(EVAL),HEAP(L),
     @                  HEAP(M),T,EUNIT,FUNIT,FRUNIT)
            CALL FREHP(M, IREAL8(3*NATMS))
            CALL FREHP(L, IREAL8(3*NATMS))
            CALL FREHP(EVAL, IREAL8(3*NATMS))
            CALL FREHP(EVECT, IREAL8(9*NATMS**2))
            CALL FREHP(GPCOV, IREAL8(9*NATMS**2))
C
C write out nominal coordinates
            IF (CUNIT.NE.0) THEN
               do i1=pointr(I)+1,pointr(I+1)
                  ii1=list(i1)
                  write(cunit,'(1x,e14.8)') xsum(ii1)
                  write(cunit,'(1x,e14.8)') ysum(ii1)
                  write(cunit,'(1x,e14.8)') zsum(ii1)
               end do
            END IF
C
C write out group mode information
            IF (MUNIT.NE.0) THEN
               write(munit,'(1x,i5)') nmod(i)
            END IF
         END DO
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE GPMODE(LIST,POINTR,I,NATMS,NATOM,AMASS,GPCOV,SIGMA,
     @                  NMOD,EVECT,EVAL,L,M,T,EUNIT,FUNIT,FRUNIT)
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'consta.inc'
      INTEGER NATMS, LIST(*), POINTR(*), I, NMOD(*)
      INTEGER FUNIT, EUNIT, NATOM
      DOUBLE PRECISION GPCOV(3*NATMS,3*NATMS), SIGMA(3*NATOM,3*NATOM)
      DOUBLE PRECISION EVECT(3*NATMS,3*NATMS), EVAL(*), AMASS(*)
      DOUBLE PRECISION L(*), M(*), T
      CHARACTER*4 FRUNIT
C
C Local
      INTEGER N1, N2, I1, I2, II1, II2
      INTEGER ISTART, ICOUNT
      DOUBLE PRECISION M1, M2, ONE, DET, C, OUT
      LOGICAL QFIRST
      PARAMETER (ONE=1.0D0, C=2.9979D10)
C
C get covariance for group
      N1=0
      DO I1=POINTR(I)+1,POINTR(I+1)
         II1=LIST(I1)
         N1=N1+1
         N2=0
         DO I2=POINTR(I)+1,POINTR(I+1)
            II2=LIST(I2)
            N2=N2+1
            GPCOV(3*N1-2,3*N2-2)=SIGMA(3*II1-2,3*II2-2)
            GPCOV(3*N1-2,3*N2-1)=SIGMA(3*II1-2,3*II2-1)
            GPCOV(3*N1-2,3*N2)=SIGMA(3*II1-2,3*II2)
            GPCOV(3*N1-1,3*N2-2)=SIGMA(3*II1-1,3*II2-2)
            GPCOV(3*N1-1,3*N2-1)=SIGMA(3*II1-1,3*II2-1)
            GPCOV(3*N1-1,3*N2)=SIGMA(3*II1-1,3*II2)
            GPCOV(3*N1,3*N2-2)=SIGMA(3*II1,3*II2-2)
            GPCOV(3*N1,3*N2-1)=SIGMA(3*II1,3*II2-1)
            GPCOV(3*N1,3*N2)=SIGMA(3*II1,3*II2)
         END DO
      END DO
C
C invert the covariance matrix
      call minvv(GPCOV,3*NATMS,det,l,m)
C
C now multiply by kT to get the effective force matrix
C also do the symmetric mass-weighting
      N1=0
      do i1=POINTR(I)+1,POINTR(I+1)
         II1=LIST(I1)
         N1=N1+1
         M1=ONE/SQRT(AMASS(II1))
         N2=0
         do i2=POINTR(I)+1,POINTR(I+1)
            II2=LIST(I2)
            N2=N2+1
            M2=ONE/SQRT(AMASS(II2))
            GPCOV(3*N1-2,3*N2-2)=KBOLTZ*T*GPCOV(3*N1-2,3*N2-2)*M1*M2
            GPCOV(3*N1-2,3*N2-1)=KBOLTZ*T*GPCOV(3*N1-2,3*N2-1)*M1*M2
            GPCOV(3*N1-2,3*N2)=KBOLTZ*T*GPCOV(3*N1-2,3*N2)*M1*M2
            GPCOV(3*N1-1,3*N2-2)=KBOLTZ*T*GPCOV(3*N1-1,3*N2-2)*M1*M2
            GPCOV(3*N1-1,3*N2-1)=KBOLTZ*T*GPCOV(3*N1-1,3*N2-1)*M1*M2
            GPCOV(3*N1-1,3*N2)=KBOLTZ*T*GPCOV(3*N1-1,3*N2)*M1*M2
            GPCOV(3*N1,3*N2-2)=KBOLTZ*T*GPCOV(3*N1,3*N2-2)*M1*M2
            GPCOV(3*N1,3*N2-1)=KBOLTZ*T*GPCOV(3*N1,3*N2-1)*M1*M2
            GPCOV(3*N1,3*N2)=KBOLTZ*T*GPCOV(3*N1,3*N2)*M1*M2
         end do
      end do
C
C diagonalize to get eigenvectors and eiganvalues
      call diagsq(3*N1,3*N1,GPCOV,evect,eval)
C
C write out eigenvectors and eigenvalues after parsing out
C all modes corresponding to non-oscillatory motions
C write the frequencies to FUNIT
C
      qfirst=.false.
      do i1=1,3*N1
         if (.not.qfirst) then
            if (eval(i1).gt.0) then
               istart=i1
               icount=1
               qfirst=.true.
               WRITE(6,'(3A)')
     @              '---------------------------------',
     @              ' MODE SUMMARY ',
     @              '---------------------------------'
               WRITE(6,'(a,i5,a,i5)')
     @              ' MODES: There are ',istart-1,
     @              ' negative eigenvalues for group ',I
               WRITE(6,'(A,i5,A)')
     @              ' MODES: ',3*N1+1-ISTART,
     @              ' positive eigenvalues remain.'
               WRITE(6,'(3A)')
     @              '---------------------------------',
     @              '--------------',
     @              '---------------------------------'
            end if
         end if
         if ((qfirst).and.(icount.le.nmod(i))) then
            icount=icount+1
            if (frunit.eq.'AKMA') then
               out=sqrt(eval(i1))
            else if (frunit.eq.'PS  ') then
               out=sqrt(eval(i1))/timfac
            else if (frunit.eq.'CM  ') then
               out=sqrt(eval(i1))*(1/(2*pi*c*timfac*1.0D-12))
            end if
               write(funit,'(1x,e14.8)') out
         end if
      end do
C
C write out modeshapes to EUNIT (if necessary)
C take care of the mass normalization first
      IF (EUNIT.NE.0) THEN
         DO I2=istart,istart+nmod(i)-1
            N1=0
            DO I1=POINTR(I)+1,POINTR(I+1)
               II1=LIST(I1)
               N1=N1+1
               M1=ONE/SQRT(AMASS(II1))
               WRITE(EUNIT,'(1X,e14.8)') M1*EVECT(3*N1-2,I2)
               WRITE(EUNIT,'(1X,e14.8)') M1*EVECT(3*N1-1,I2)
               WRITE(EUNIT,'(1X,e14.8)') M1*EVECT(3*N1,I2)
            END DO
         END DO
      END IF
C
      RETURN
      END
C
