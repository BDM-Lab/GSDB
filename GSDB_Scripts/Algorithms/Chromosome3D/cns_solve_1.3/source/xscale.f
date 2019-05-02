! OpenMP version KD 1/2007
C======================================================================
      SUBROUTINE XSCALIT(KSCAL,BSCAL,NNSELE,INDX,MARK,
     &                   XRH,XRK,XRL,
     &                   MAXSET,TGSET,NNSET,HPSFSET,TYPESF,SET,
     &                   MODE,NCYC,DIAG,EPS,KSMIN,BFMIN,BFMAX,
     &                   FIXMOD,KRES,KINI,BINI,QUPDA,FFK,RVAL,TARG,
     &                   XRFFKQ,FFK1,FFK2,QFFK1,QFFK2,QTLOW,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   SILENT,QANISO,QISO,RESTRC,XRTR,XRCELL,XRVOL)
C
C Performes least-squares minimization.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'timer.inc'
      DOUBLE PRECISION KSCAL(*), BSCAL(*)
      INTEGER NNSELE
      INTEGER INDX(*)
      DOUBLE PRECISION MARK
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MAXSET, TGSET
      INTEGER NNSET(*), HPSFSET(*), TYPESF(*)
      CHARACTER*(*) SET(*)
      CHARACTER*4 MODE
      INTEGER NCYC, DIAG, FIXMOD
      DOUBLE PRECISION EPS, KSMIN, BFMIN, BFMAX, KRES
      DOUBLE PRECISION KINI, BINI
      LOGICAL QUPDA
      DOUBLE PRECISION FFK
      LOGICAL XRFFKQ
      DOUBLE PRECISION FFK1, FFK2
      LOGICAL QFFK1, QFFK2, QTLOW
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      INTEGER SILENT
      LOGICAL QANISO, QISO
      CHARACTER*4 RESTRC
      DOUBLE PRECISION XRTR(3,3), XRCELL(6), XRVOL
C pointers
      INTEGER ARG
      INTEGER PARAM, GDERI, PARAM0, FIX, PWORK
      INTEGER AA, A0, H
C local
      INTEGER M, M1, MPERSET
      DOUBLE PRECISION FOSUM, FO2SUM, FOSUM1, FOSUM2
      DOUBLE PRECISION RVAL, TARG, RLOW, TLOW, RHIG, THIG
      LOGICAL ERROR
      INTEGER IERR
      DOUBLE PRECISION ABCS(3), CABGS(3)
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C begin
      ERROR=.FALSE.
      IERR=0
C
      IF (QANISO) THEN
C number of parameters per set for anisotropic is 7
      MPERSET=7
C compute reciprocal unit cell dimensions
      CALL XRABCS(ABCS,CABGS,XRCELL,XRVOL)
      ELSE
      MPERSET=2
      END IF
C
C number of parameters
      M=MPERSET*MAXSET
      M1=M+1
C
C allocate HEAP space for parameter arraies
      ARG=ALLHP(ICPLX8(MAXSET))
      PARAM=ALLHP(IREAL8(M))
      GDERI=ALLHP(IREAL8(M1))
      PARAM0=ALLHP(IREAL8(M))
      FIX=ALLHP(INTEG4(M))
      PWORK=ALLHP(INTEG4(M1))
      AA=ALLHP(IREAL8(M*M1))
      A0=ALLHP(IREAL8(M*M1))
      H=ALLHP(IREAL8(M1))
C
C
C if KSCAL, BSCAL .eq. MARK means that they need to be refined!!!
C install least-squares parameters
      CALL XSCPARM(MAXSET,NNSET,KSCAL,BSCAL,MARK,HEAP(PARAM),HEAP(FIX),
     &             KINI,BINI,MPERSET,QANISO,QISO,RESTRC,CABGS)
C
C optimal all parameters for least-squares targets
      IF (MODE.EQ.'TARG'.OR.MODE.EQ.'TLOW') THEN
      CALL XSCLSQ(NCYC,DIAG,M,M1,EPS,
     &            NNSELE,INDX,XRH,XRK,XRL,HEAP(ARG),
     &            HEAP(PARAM),HEAP(GDERI),HEAP(PARAM0),
     &            HEAP(FIX),HEAP(PWORK),HEAP(AA),HEAP(A0),HEAP(H),
     &            FIXMOD,FFK,XRFFKQ,KRES,FFK1,FFK2,QFFK1,QFFK2,
     &            MAXSET,TGSET,NNSET,HPSFSET,TYPESF,SET,
     &            KSMIN,BFMIN,BFMAX,ERROR,IERR,QTLOW,
     &            QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &            SILENT,QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,XRTR)
      ELSE
      ERROR=.TRUE.
      WRITE(6,'(A,A)') ' %XSCALE-err: no such target: ',MODE
      END IF
C
C final overall scale factor, target, R value and paramters
      IF (.NOT.ERROR) THEN
C
C update the overall scale factor
      CALL XSCFFKO(FFK,XRFFKQ,FOSUM,FO2SUM,ERROR,
     &             ZERO,FFK1,FFK2,QFFK1,QFFK2,FOSUM1,FOSUM2,
     &             HEAP(ARG),HEAP(PARAM),
     &             XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
C
C calculate the target and R-factor
      CALL XSCTARG(FFK,FOSUM,FO2SUM,ZERO,FFK1,FFK2,
     &             FOSUM1,FOSUM2,
     &             HEAP(ARG),HEAP(PARAM),
     &             XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             RVAL,TARG,RLOW,TLOW,RHIG,THIG,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
C
C output
      IF (WRNLEV.GE.SILENT) THEN
      WRITE(6,'(A,A,A,A)') ' XSCALE:',
     & ' ---------------------- Final ','-------------',
     & '-----------------------'
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCALE: ',
     &   'Target=',TARG,'  Scale= ',FFK,'  R-factor= ',RVAL
      IF (IERR.EQ.1) THEN
      WRITE(6,'(A)') ' XSCALE: stop by the limited number of cycles.'
      ELSE
      WRITE(6,'(A)') ' XSCALE: normal termination - convergence.'
      END IF
      END IF
C
C transfer parameters from PARAM to KSCAL and BSCAL
      CALL XSCPARU(MAXSET,NNSET,KSCAL,BSCAL,HEAP(PARAM),
     &             HEAP(FIX),TGSET,FFK,ERROR,QUPDA,SILENT,
     &             MPERSET)
C
C report error status
      ELSE
      WRITE(6,'(A)') ' %XSCALE: aborting termination - errors.'
      IF (IERR.GE.40) THEN
      WRITE(6,'(A)') ' %%%%XSCALE-err: ill matrix - no solution.'
      CALL XSCPARU(MAXSET,NNSET,KSCAL,BSCAL,HEAP(PARAM),
     &             HEAP(FIX),TGSET,FFK,ERROR,QUPDA,0,
     &             MPERSET)
      ELSE IF (IERR.GE.30) THEN
      WRITE(6,'(A)') ' %%%XSCALE-err: diverging from solution.'
      CALL XSCPARU(MAXSET,NNSET,KSCAL,BSCAL,HEAP(PARAM0),
     &             HEAP(FIX),TGSET,FFK,.FALSE.,QUPDA,0,
     &             MPERSET)
      ELSE IF (IERR.GE.20) THEN
      WRITE(6,'(A)') ' %%XSCALE-err: all parameters are fixed.'
      ELSE IF (IERR.GE.10) THEN
      WRITE(6,'(A)') ' %XSCALE-err: no data selected.'
      CALL XSCPARU(MAXSET,NNSET,KSCAL,BSCAL,HEAP(PARAM),
     &             HEAP(FIX),TGSET,FFK,ERROR,QUPDA,0,
     &             MPERSET)
      END IF
C
      END IF
C
C free heap
      CALL FREHP(H,IREAL8(M1))
      CALL FREHP(A0,IREAL8(M*M1))
      CALL FREHP(AA,IREAL8(M*M1))
      CALL FREHP(PWORK,INTEG4(M1))
      CALL FREHP(FIX,INTEG4(M))
      CALL FREHP(PARAM0,IREAL8(M))
      CALL FREHP(GDERI,IREAL8(M1))
      CALL FREHP(PARAM,IREAL8(M))
      CALL FREHP(ARG,ICPLX8(MAXSET))
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCPARM(MAXSET,NNSET,KSCAL,BSCAL,MARK,PARAM,FIX,
     &                   KINI,BINI,MPERSET,QANISO,QISO,RESTRC,CABGS)
C
C Routine installs PARAM array and checks initially fixed parameters.
C Parameter has been fixed, FIX=1; otherwise non-fixed, FIX=0.
C if KSCAL, BSCAL .eq. MARK means that they need to be refined!!!
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER MAXSET
      DOUBLE PRECISION MARK
      INTEGER NNSET(*)
      DOUBLE PRECISION KSCAL(*), BSCAL(*)
      DOUBLE PRECISION PARAM(*)
      INTEGER FIX(*)
      DOUBLE PRECISION KINI, BINI
      INTEGER MPERSET
      LOGICAL QANISO, QISO
      CHARACTER*4 RESTRC
      DOUBLE PRECISION CABGS(3)
C local
      INTEGER I, I1, J, K
C parameter
      DOUBLE PRECISION ZERO, CSMALL
      PARAMETER (ZERO=0.0D0, CSMALL=0.001D0)
C begin
C
C install K scales
      DO I=1,MAXSET
      K=NNSET(I)
      J=(I-1)*MPERSET+1
      IF (KSCAL(K).EQ.MARK) THEN
      PARAM(J)=KINI
      FIX(J)=0
      ELSE
      PARAM(J)=KSCAL(K)
      FIX(J)=1
      END IF
      END DO
C
C install B scales
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      K=(NNSET(I)-1)*(MPERSET-1)
      DO J=1,MPERSET-1
      IF (BSCAL(K+J).EQ.MARK) THEN
      PARAM(I1+J)=BINI
      FIX(I1+J)=0
      ELSE
      PARAM(I1+J)=BSCAL(K+J)
      FIX(I1+J)=1
      END IF
      END DO
C impose symmetry restrictions for QANISO
      IF (QANISO.AND.(RESTRC.EQ.'ALL'.OR.RESTRC.EQ.'OFFD')) THEN
      IF (ABS(CABGS(1)).LT.CSMALL) THEN
      PARAM(I1+6)=ZERO
      FIX(I1+6)=1
      BSCAL(K+6)=ZERO
      END IF
      IF (ABS(CABGS(2)).LT.CSMALL) THEN
      PARAM(I1+5)=ZERO
      FIX(I1+5)=1
      BSCAL(K+5)=ZERO
      END IF
      IF (ABS(CABGS(3)).LT.CSMALL) THEN
      PARAM(I1+4)=ZERO
      FIX(I1+4)=1
      BSCAL(K+4)=ZERO
      END IF
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCPARU(MAXSET,NNSET,KSCAL,BSCAL,PARAM,FIX,
     &                   TGSET,FFK,ERROR,QUPDA,SILENT,MPERSET)
C
C Routine updates parameters - from PARAM to KSCAL and BSCAL
C KSCAL are multiplied by FFK if QUPDA is true
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER MAXSET
      INTEGER NNSET(*)
      DOUBLE PRECISION KSCAL(*), BSCAL(*)
      DOUBLE PRECISION PARAM(*)
      INTEGER FIX(*), TGSET
      DOUBLE PRECISION FFK
      LOGICAL ERROR, QUPDA
      INTEGER SILENT
      INTEGER MPERSET
C local
      INTEGER I, I1, J, K
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
      IF (.NOT.ERROR) THEN
C copy from PARAM
      DO I=1,MAXSET
      K=NNSET(I)
      J=(I-1)*MPERSET+1
      IF (FIX(J).EQ.0) KSCAL(K)=PARAM(J)
      END DO
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      K=(NNSET(I)-1)*(MPERSET-1)
      DO J=1,MPERSET-1
      IF (FIX(I1+J).EQ.0) BSCAL(K+J)=PARAM(I1+J)
      END DO
      END DO
C
      ELSE
C set k to one and B to zero if error
      DO I=1,MAXSET
      K=NNSET(I)
      J=(I-1)*MPERSET+1
      IF (FIX(J).EQ.0) KSCAL(K)=ONE
      END DO
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      K=(NNSET(I)-1)*(MPERSET-1)
      DO J=1,MPERSET-1
      IF (FIX(I1+J).EQ.0) BSCAL(K+J)=ZERO
      END DO
      END DO
      END IF
C
      IF (QUPDA) THEN
      IF (WRNLEV.GE.SILENT) WRITE(6,'(A,A,F8.3)')
     & ' XSCALE: K<i> has been multiplied by',
     & ' the overall scale ($XSCFFK)',FFK
      DO I=1,MAXSET
      IF (I.NE.TGSET) THEN
      I1=NNSET(I)
      KSCAL(I1)=FFK*KSCAL(I1)
      END IF
      END DO
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCFIXC(M,FIX,PARAM,PWORK,MWORK,MWORK1,ERROR)
C
C Checks non-fixed (FIX=0) parameters and finds out
C number of working parameters.
C Indices store in PWORK.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER M
      INTEGER FIX(*)
      DOUBLE PRECISION PARAM(*)
      INTEGER PWORK(*), MWORK, MWORK1
      LOGICAL ERROR
C local
      INTEGER I
C begin
      ERROR=.FALSE.
C
      MWORK=0
      DO I=1,M
      PWORK(I)=0
      END DO
C
      DO I=1,M
      IF (FIX(I).EQ.0) THEN
      MWORK=MWORK+1
      PWORK(MWORK)=I
      END IF
      END DO
      MWORK1=MWORK+1
      PWORK(MWORK1)=M+1
C
      IF (MWORK.EQ.0) THEN
      ERROR=.TRUE.
      END IF
C
      RETURN
      END
C=====================================================================
      SUBROUTINE XSCLSQ(NCYC,DIAG,M,M1,EPS,
     &                  NNSELE,INDX,XRH,XRK,XRL,ARG,
     &                  PARAM,GDERI,PARAM0,FIX,PWORK,AA,A0,H,
     &                  FIXMOD,FFK,XRFFKQ,KRES,FFK1,FFK2,QFFK1,QFFK2,
     &                  MAXSET,TGSET,NNSET,HPSFSET,TYPESF,SET,
     &                  KSMIN,BFMIN,BFMAX,ERROR,IERR,QTLOW,
     &                  QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,SILENT,
     &                  QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,XRTR)
C
C Routine minimizes the target by the Least-Squares method.
C The target is a constructed non-linear fuction of M parameters.
C
C TARG = sum{ W * ( |FTARG| - FFK * |FTOT| )^2 } / sum{ W * FTARG^2 }
C
C where FTARG is the destined ('obseved') data set
C and FTOT=f(h,PARAM) is the data set ('calculated') to be refined,
C which can be the sum of multiple data sets.
C
C Notes:
C ------
C 1)  M is the number of parameters, M1=M+1.
C 2)  PARAM(M) is the parameter array, PARAM0(M) is the working array.
C     GDERI(M1) is the derivative array and FIX(M) is a flag array.
C 3)  AA(M,M1) is the normal matrix,
C     A0(M,M1) and H(M1) are the working arraies.
C 4)  MPWRK is the actual number of parameters to be determined,
C     and PWORK is the index array of non-fixed parameters.
C 5)  NCYC is the number of least-squares cycles.
C 6)  EPS is for convergence control.
C 7)  DIAG is diagonal approximation (default: full matrix)
C     (0:full matrix, 1:diagonal, 2:block of 2x2).
C 8)  A resistant factor D adding to the diagonal elements
C     would slow down divergence and speed up convergence.
C 9)  When ID is negative the routine will use the old matrix.
C 10) A subroutine of computing derivatives GDERI is required.
C 11) Solve a matrix AA(M,M1) and result in GDERI(M1).
C 12) New parameters update by PARAM=PARAM+GDERI.
C
C Error codes:
C ------------
C     IERR=0   success minimization
C     IERR=1   stop by the limited number of cycles
C     IERR=10  no data were selected or data went to zero
C     IERR=20  all parameters were fixed;
C     IERR=30  divergence
C     IERR=40  no least-squares solution - ill matrix
C
C 1982.10.12  Jiansheng Jiang in Beijing
C ======================================
C Adapted and modified for CNS
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INTEGER NCYC, DIAG, M, M1
      DOUBLE PRECISION EPS
      INTEGER NNSELE
      INTEGER INDX(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX ARG(*)
      DOUBLE PRECISION PARAM(*), GDERI(*), PARAM0(*)
      INTEGER FIX(*), PWORK(*)
      DOUBLE PRECISION AA(M,M1), A0(M,M1), H(M1)
      INTEGER FIXMOD
      DOUBLE PRECISION FFK, KRES, FFK1, FFK2
      LOGICAL XRFFKQ, QFFK1, QFFK2, QTLOW
      INTEGER MAXSET, TGSET
      INTEGER NNSET(*), HPSFSET(*), TYPESF(*)
      CHARACTER*6 SET(*)
      DOUBLE PRECISION KSMIN, BFMIN, BFMAX
      LOGICAL ERROR
      INTEGER IERR
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      INTEGER SILENT
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), CABGS(3)
      CHARACTER*4 RESTRC
      DOUBLE PRECISION XRTR(3,3)
C local
      INTEGER I, J, KK, ICYCLE, ISET
      INTEGER ID, MPWRK, MPWRK1
      DOUBLE PRECISION FOSUM, FO2SUM, FOSUM1, FOSUM2
      DOUBLE PRECISION TARG, RVAL, TARG0, RVAL0, TC, RC
      DOUBLE PRECISION TLOW, RLOW, THIG, RHIG
      DOUBLE PRECISION D, DOWN
      LOGICAL CNVG, DIVG
C parameters
      DOUBLE PRECISION ZERO, ONE
      DOUBLE PRECISION QUART, TEN, RSMALL
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      PARAMETER (QUART=0.25D0, TEN=10.0D0)
      PARAMETER (RSMALL=1.0D-6)
C
C begin
C
C initialization
C
      DO I=1,M
      PARAM0(I)=PARAM(I)
      END DO
      DO I=1,M1
      GDERI(I)=ZERO
      END DO
C
      ERROR=.FALSE.
      IERR=0
C
C calculate the initial overall scale factor
      CALL XSCFFKO(FFK,XRFFKQ,FOSUM,FO2SUM,ERROR,
     &             KRES,FFK1,FFK2,QFFK1,QFFK2,FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      IF (ERROR) IERR=10
C
      IF (.NOT.ERROR) THEN
C
C use the overall scale factor (FFK) if the user specified
      IF (.NOT.XRFFKQ) THEN
      FFK1=FFK
      FFK2=FFK
      END IF
C
C calculate the target and R-factor by using initial FFKs
      CALL XSCTARG(FFK,FOSUM,FO2SUM,KRES,FFK1,FFK2,
     &             FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             RVAL,TARG,RLOW,TLOW,RHIG,THIG,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
C
      ICYCLE=0
      ID=0
      D=0.01D0
C
      IF (WRNLEV.GE.SILENT) THEN
      WRITE(6,'(A,A,I3,2X,F8.4,A)') ' XSCLSQ:',
     & ' ---------------------- Cycle ',ICYCLE,D,
     & ' ----------------------'
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'Target=',TARG,'  Scale= ',FFK,'  R-factor= ',RVAL
      IF (KRES.GT.RSMALL) THEN
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'T-low =',TLOW,'  Scale= ',FFK1,'  R-factor= ',RLOW
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'T-high=',THIG,'  Scale= ',FFK2,'  R-factor= ',RHIG
      END IF
      DO I=1,MAXSET
      ISET=NNSET(I)
      J=(I-1)*MPERSET+1
      IF (QANISO) THEN
      WRITE(6,'(A,I2,3A,F8.3,A)')
     & ' XSCLSQ: set',ISET,' ',SET(ISET),'  kscal=',PARAM(J),
     & '  bscal_tensor='
      WRITE(6,'(A,6(A,F6.2))') ' XSCLSQ:',
     & ' B11=',PARAM(J+1),
     & ' B22=',PARAM(J+2),
     & ' B33=',PARAM(J+3),
     & ' B12=',PARAM(J+4),
     & ' B13=',PARAM(J+5),
     & ' B23=',PARAM(J+6)
      ELSE
      WRITE(6,'(A,I2,3A,F8.3,A,F8.2,5X,F8.3,F8.2)')
     & ' XSCLSQ: set',ISET,' ',SET(ISET),'  kscal=',PARAM(J),
     & '  bscal=',PARAM(J+1),GDERI(J),GDERI(J+1)
      END IF
      END DO
      END IF
C
      END IF
C
C end of initialization
C
C
C least-squares cycles entre here
C
      ICYCLE=1
      CNVG=.FALSE.
      DO WHILE (ICYCLE.LE.NCYC.AND..NOT.CNVG.AND..NOT.ERROR)
C
      IF (ID.GE.0) THEN
      IF (KRES.GT.RSMALL) THEN
      IF (QTLOW) THEN
      TARG0=TLOW
      RVAL0=RLOW
      ELSE
      TARG0=TLOW+THIG
      RVAL0=RLOW+RHIG
      END IF
      ELSE
      TARG0=TARG
      RVAL0=RVAL
      END IF
      DO I=1,M
      IF (FIX(I).EQ.0) THEN
      PARAM0(I)=PARAM(I)
      END IF
      END DO
      END IF
C
C check non-fixed parameters
      CALL XSCFIXC(M,FIX,PARAM,PWORK,MPWRK,MPWRK1,ERROR)
      IF (ERROR) IERR=20
C
      IF (.NOT.ERROR) THEN
C
C calculate derivatives and solve the matrix
      CALL XSCLSQ2(DIAG,M,M1,TARG0,ID,D,
     &             NNSELE,INDX,XRH,XRK,XRL,ARG,
     &             PARAM,GDERI,AA,A0,H,FIX,FIXMOD,
     &             PWORK,MPWRK,MPWRK1,FFK,KRES,FFK1,FFK2,
     &             MAXSET,TGSET,HPSFSET,TYPESF,ERROR,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,XRTR)
      IF (ERROR) IERR=40
C
      IF (.NOT.ERROR) THEN
C
      ID=ID+1
C update the parameters
      DO I=1,M
      IF (FIX(I).EQ.0) THEN
      PARAM(I)=PARAM0(I)+GDERI(I)
      END IF
      END DO
C
C impose symmetry restrictions according to crystal system
      IF (QANISO) THEN
      CALL XSCRSTRC(RESTRC,MAXSET,TGSET,MPERSET,PARAM,ABCS,CABGS,QISO)
      END IF
C
C restraints on scale factors
      DO I=1,M,MPERSET
      IF (FIX(I).EQ.0) THEN
      IF (PARAM(I).LT.KSMIN) PARAM(I)=KSMIN
      END IF
      END DO
C
C restraints on B factors
      DO I=1,M,MPERSET
      DO J=1,MPERSET-1
      IF (FIX(I+J).EQ.0) THEN
      IF (PARAM(I+J).LT.BFMIN) PARAM(I+J)=BFMIN
      IF (PARAM(I+J).GT.BFMAX) PARAM(I+J)=BFMAX
      END IF
      END DO
      END DO
C
C only one scale factor if FIXMOD=3
      IF (FIXMOD.EQ.3) THEN
      KK=0
      DO I=1,M,MPERSET
      IF (FIX(I).EQ.0.AND.KK.EQ.0) KK=I
      IF (FIX(I).EQ.0.AND.KK.GT.0) PARAM(I)=PARAM(KK)
      END DO
C only one B factor if FIXMOD=4
      ELSE IF (FIXMOD.EQ.4) THEN
      KK=0
      DO I=1,M,MPERSET
      DO J=1,MPERSET-1
      IF (FIX(I+J).EQ.0.AND.KK.EQ.0) KK=I+J
      IF (FIX(I+J).EQ.0.AND.KK.GT.0) PARAM(I+J)=PARAM(KK)
      END DO
      END DO
      END IF
C
C update the overall scale factor
      CALL XSCFFKO(FFK,XRFFKQ,FOSUM,FO2SUM,ERROR,
     &             KRES,FFK1,FFK2,QFFK1,QFFK2,FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      IF (ERROR) IERR=10
C
      IF (.NOT.ERROR) THEN
C
C we have to use the over scale factor (FFK) if the user specified
      IF (.NOT.XRFFKQ) THEN
      FFK1=FFK
      FFK2=FFK
      END IF
C
C calculate the target and R-factor
      CALL XSCTARG(FFK,FOSUM,FO2SUM,KRES,FFK1,FFK2,
     &             FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             RVAL,TARG,RLOW,TLOW,RHIG,THIG,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      IF (KRES.GT.RSMALL) THEN
      IF (QTLOW) THEN
      TC=TLOW
      RC=RLOW
      ELSE
      TC=TLOW+THIG
      RC=RLOW+RHIG
      END IF
      ELSE
      TC=TARG
      RC=RVAL
      END IF
C
C output information
      IF (WRNLEV.GE.SILENT) THEN
      WRITE(6,'(A,A,I3,2X,F8.4,A)') ' XSCLSQ:',
     & ' ---------------------- Cycle ',ICYCLE,D,
     & ' ----------------------'
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'Target=',TARG,'  Scale= ',FFK,'  R-factor= ',RVAL
      IF (KRES.GT.RSMALL) THEN
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'T-low =',TLOW,'  Scale= ',FFK1,'  R-factor= ',RLOW
      WRITE(6,'(A,A,F8.4,A,F8.3,A,F8.4)') ' XSCLSQ: ',
     &   'T-high=',THIG,'  Scale= ',FFK2,'  R-factor= ',RHIG
      END IF
      DO I=1,MAXSET
      ISET=NNSET(I)
      J=(I-1)*MPERSET+1
      IF (QANISO) THEN
      WRITE(6,'(A,I2,3A,F8.3,A)')
     & ' XSCLSQ: set',ISET,' ',SET(ISET),'  kscal=',PARAM(J),
     & '  bscal_tensor='
      WRITE(6,'(A,6(A,F6.2))') ' XSCLSQ:',
     & ' B11=',PARAM(J+1),
     & ' B22=',PARAM(J+2),
     & ' B33=',PARAM(J+3),
     & ' B12=',PARAM(J+4),
     & ' B13=',PARAM(J+5),
     & ' B23=',PARAM(J+6)
      ELSE
      WRITE(6,'(A,I2,3A,F8.3,A,F8.2,5X,F8.3,F8.2)')
     & ' XSCLSQ: set',ISET,' ',SET(ISET),'  kscal=',PARAM(J),
     & '  bscal=',PARAM(J+1),GDERI(J),GDERI(J+1)
      END IF
      END DO
      END IF
C
C
C criteria for the minimization and convergence
C
C check if any parameter has a large shift
      CNVG=.TRUE.
      I=1
      DO WHILE (CNVG.AND.I.LE.M)
      IF (FIX(I).EQ.0) THEN
      IF (ABS(GDERI(I))/(ABS(PARAM(I))+RSMALL).GE.EPS) CNVG=.FALSE.
      END IF
      I=I+1
      END DO
C
C check whether the minimum is reached
      IF (.NOT.CNVG) THEN
      IF (ICYCLE.GT.1.AND.ABS(TC-TARG0)/TARG0.LT.EPS) CNVG=.TRUE.
      IF (ICYCLE.GT.1.AND.ABS(RC-RVAL0)/RVAL0.LT.EPS) CNVG=.TRUE.
      IF (RC.LT.EPS.OR.TC.LT.RSMALL) CNVG=.TRUE.
      END IF
C
C
      IF (.NOT.CNVG) THEN
C
      DIVG=.FALSE.
      IF (TC.LT.TARG0.OR.RC.LT.RVAL0) THEN
C speed up convergence
      ID=1
      IF (D.GT.RSMALL) D=D/TEN
      ELSE
C apply the resistant factor to slow down divergence
      ID=-1
      D=TEN*D
      END IF
      IF (D.GT.TEN*TEN) THEN
      DIVG=.TRUE.
      D=TEN*TEN
      END IF
C
C
C rescue by minifing the derivatives if divergence
      DOWN=QUART
      DO WHILE (DIVG.AND.DOWN.GT.EPS)
      DO I=1,M
      IF (FIX(I).EQ.0) THEN
      PARAM(I)=PARAM0(I)+DOWN*GDERI(I)
      END IF
      END DO
C update the overall scale factor
      CALL XSCFFKO(FFK,XRFFKQ,FOSUM,FO2SUM,ERROR,
     &             KRES,FFK1,FFK2,QFFK1,QFFK2,FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
C use the overall scale factor (FFK) if the user specified
      IF (.NOT.XRFFKQ) THEN
      FFK1=FFK
      FFK2=FFK
      END IF
C calculate the target and R-factor
      CALL XSCTARG(FFK,FOSUM,FO2SUM,KRES,FFK1,FFK2,
     &             FOSUM1,FOSUM2,
     &             ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             RVAL,TARG,RLOW,TLOW,RHIG,THIG,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      IF (KRES.GT.RSMALL) THEN
      IF (QTLOW) THEN
      TC=TLOW
      RC=RLOW
      ELSE
      TC=TLOW+THIG
      RC=RLOW+RHIG
      END IF
      ELSE
      TC=TARG
      RC=RVAL
      END IF
      IF (TC.LT.TARG0.OR.RC.LT.RVAL0) DIVG=.FALSE.
      DOWN=DOWN*QUART
      END DO
C
C divergence
      IF (DIVG) THEN
      ERROR=.TRUE.
      IERR=30
      END IF
C
C end of criteria
      END IF
C
C end of error conditions
      END IF
      END IF
      END IF
C
C end of least-squares cycles
      ICYCLE=ICYCLE+1
      IF (ICYCLE.GT.NCYC) THEN
      CNVG=.TRUE.
      IERR=1
      END IF
      END DO
C
C
      RETURN
      END
C =================================================================
      SUBROUTINE XSCLSQ2(DIAG,M,M1,TARG,ID,D,
     &                   NNSELE,INDX,XRH,XRK,XRL,ARG,
     &                   PARAM,GDERI,AA,A0,H,FIX,FIXMOD,
     &                   PWORK,MPWRK,MPWRK1,FFK,KRES,FFK1,FFK2,
     &                   MAXSET,TGSET,HPSFSET,TYPESF,ERROR,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,
     &                   XRTR)
C
C Routine gives a solution of the Least-squares target.
C See routine XSCLSQ above.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER DIAG, M, M1, ID
      DOUBLE PRECISION TARG, D
      INTEGER NNSELE
      INTEGER INDX(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MAXSET
      DOUBLE COMPLEX ARG(MAXSET)
      DOUBLE PRECISION PARAM(*), GDERI(m1)
      DOUBLE PRECISION AA(M,M1), A0(M,M1), H(M1)
      INTEGER FIX(*), FIXMOD
      INTEGER PWORK(M1), MPWRK, MPWRK1
      DOUBLE PRECISION FFK, KRES, FFK1, FFK2
      INTEGER TGSET
      INTEGER HPSFSET(*), TYPESF(*)
      LOGICAL ERROR
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), CABGS(3)
      CHARACTER*4 RESTRC
      DOUBLE PRECISION XRTR(3,3)
C local
      INTEGER I, J, K, I1, II, REFLCT, IREF
      DOUBLE PRECISION C
C parameters
      DOUBLE PRECISION ZERO, ONE, RSMALL
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,RSMALL=1.0D-12)
C
C begin
      ERROR=.FALSE.
C
      IF (ID.GE.0) THEN
C
C initialize matrix
      DO I=1,MPWRK
      DO J=1,MPWRK1
      A0(I,J)=ZERO
      END DO
      END DO
C initialize derivatives
      DO J=1,MPWRK1
      GDERI(J)=ZERO
      END DO
C
      IF (TARG.LT.RSMALL) RETURN
C
C calculate derivatives and collect matrix
!$omp parallel do default(shared) reduction(+:a0)
!$omp& private(iref,reflct,arg,gderi)
!$omp& schedule(guided)
      DO IREF=1,NNSELE
      REFLCT=INDX(IREF)
      CALL XSCDERI(M,M1,REFLCT,ARG,PARAM,GDERI,XRH,XRK,XRL,
     &             FIX,FIXMOD,FFK,KRES,FFK1,FFK2,
     &             MAXSET,TGSET,HPSFSET,TYPESF,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,XRTR)
      DO I=1,MPWRK
      DO J=1,MPWRK1
      A0(I,J)=A0(I,J)+GDERI(PWORK(I))*GDERI(PWORK(J))
      END DO
      END DO
      END DO
C
C normalize matrix
C full matrix approach only
      IF (DIAG.EQ.0) THEN
      DO I=1,MPWRK
      IF (A0(I,I).LE.ZERO) THEN
      ERROR=.TRUE.
      ELSE
      H(I)=ONE/SQRT(A0(I,I))
      END IF
      END DO
      H(MPWRK1)=ONE/SQRT(TARG)
      DO I=1,MPWRK
      I1=I+1
      DO J=I1,MPWRK1
      A0(I,J)=A0(I,J)*H(I)*H(J)
      END DO
      END DO
      END IF
C
      END IF
C
      IF (.NOT.ERROR) THEN
C
C transfer matrix
      DO I=1,MPWRK
      DO J=1,MPWRK1
      AA(I,J)=A0(I,J)
      END DO
      END DO
C
C Diagonal approximation
      IF (DIAG.EQ.1) THEN
C
      DO I=1,MPWRK
      H(I)=AA(I,I)+D
      IF (H(I).LE.RSMALL) THEN
      ERROR=.TRUE.
      ELSE
      GDERI(PWORK(I))=AA(I,MPWRK1)/H(I)
      END IF
      END DO
C
C Diagonal approximation by a block of 2
      ELSE IF (DIAG.EQ.2) THEN
C
      DO I=1,MPWRK-1,2
      J=I+1
      H(I)=AA(I,I)*AA(J,J)-AA(I,J)*AA(J,I)+D
      IF (H(I).LE.RSMALL) THEN
      ERROR=.TRUE.
      ELSE
      C=ONE/H(I)
      END IF
      GDERI(PWORK(I))=(AA(J,J)*AA(I,MPWRK1)-AA(I,J)*AA(J,MPWRK1))*C
      GDERI(PWORK(J))=(AA(I,I)*AA(J,MPWRK1)-AA(J,I)*AA(I,MPWRK1))*C
      END DO
      IF (MOD(MPWRK,2).EQ.1) THEN
      GDERI(PWORK(MPWRK))=AA(MPWRK,MPWRK1)/(AA(MPWRK,MPWRK)+D)
      END IF
C
C full matrix approach
      ELSE
C
C pack diagonal of matrix
      DO I=1,MPWRK
      AA(I,I)=ONE+D
      END DO
C
C inverse the matrix
      DO I=1,MPWRK
      IF (AA(I,I).LE.RSMALL) THEN
      ERROR=.TRUE.
      ELSE
      C=ONE/AA(I,I)
      END IF
      I1=I+1
      DO J=I1,MPWRK1
      IF (J.NE.MPWRK1) AA(J,I)=AA(I,J)
      AA(I,J)=C*AA(I,J)
      END DO
      IF (I1.LT.MPWRK1) THEN
      DO K=I1,MPWRK
      C=AA(K,I)
      DO J=K,MPWRK1
      AA(K,J)=AA(K,J)-C*AA(I,J)
      END DO
      END DO
      END IF
      END DO
C
C solve the matrix
      C=ONE/H(MPWRK1)
      DO II=1,MPWRK
      I=MPWRK1-II
      I1=I+1
      IF (I1.LT.MPWRK1) THEN
      DO J=I1,MPWRK
      AA(I,MPWRK1)=AA(I,MPWRK1)-AA(J,MPWRK1)*AA(I,J)
      END DO
      END IF
      GDERI(PWORK(I))=AA(I,MPWRK1)*H(I)*C
      END DO
C
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCFFKO(FFK,XRFFKQ,FOSUM,FO2SUM,ERROR,
     &                   KRES,FFK1,FFK2,QFFK1,QFFK2,FOSUM1,FOSUM2,
     &                   ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &                   MAXSET,TGSET,HPSFSET,TYPESF,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   QANISO,QISO,MPERSET,ABCS,XRTR)
C
C Calculates the overall scale factor between FTARG and FTOT
C and various scaling constants
C For KRES > 0.0 two scale factors in two resolution ranges
C are also calculated.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER MAXSET
      DOUBLE PRECISION FFK, FOSUM, FO2SUM
      LOGICAL ERROR
      LOGICAL XRFFKQ, QFFK1, QFFK2
      DOUBLE PRECISION KRES, FFK1, FFK2, FOSUM1, FOSUM2
      DOUBLE COMPLEX ARG(maxset)
      DOUBLE PRECISION PARAM(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER NNSELE
      INTEGER INDX(*)
      INTEGER TGSET
      INTEGER HPSFSET(*), TYPESF(*)
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), XRTR(3,3)
C local
      INTEGER REFLCT, IREF
      DOUBLE COMPLEX FTARG, FTOT
      DOUBLE PRECISION FC2SUM, FOCSUM
      DOUBLE PRECISION FC2SUM1, FOCSUM1, FC2SUM2, FOCSUM2
      DOUBLE PRECISION SSQ, AFO1, AFO2, AFC1, AFC2, WFC2, WFOC
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      FOSUM=ZERO
      FOSUM1=ZERO
      FOSUM2=ZERO
      FO2SUM=ZERO
      FC2SUM=ZERO
      FOCSUM=ZERO
      FC2SUM1=ZERO
      FOCSUM1=ZERO
      FC2SUM2=ZERO
      FOCSUM2=ZERO
C
!$omp parallel do private(iref,reflct,FTARG,FTOT,afo2,afc2,afo1,afc1,
!$omp& wfc2,wfoc,ssq,arg) default(shared) reduction(+:fosum,fo2sum,
!$omp& fc2sum,focsum,fosum1,fc2sum1,focsum1,fosum2,fc2sum2,focsum2)
!$omp& schedule(guided)
      DO IREF=1,NNSELE
      REFLCT=INDX(IREF)
C
      CALL XSCFCAL(REFLCT,ARG,PARAM,XRH,XRK,XRL,SSQ,
     &             MAXSET,TGSET,HPSFSET,TYPESF,FTARG,FTOT,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      AFO2=DBLE(FTARG)**2+DIMAG(FTARG)**2
      AFC2=DBLE(FTOT)**2+DIMAG(FTOT)**2
      AFO1=SQRT(AFO2)
      AFC1=SQRT(AFC2)
      FOSUM=FOSUM + AFO1
      FO2SUM=FO2SUM +AFO2
      WFC2=AFC2
      WFOC=AFO1*AFC1
      FC2SUM=FC2SUM +WFC2
      FOCSUM=FOCSUM +WFOC
C
      IF (KRES.GT.RSMALL) THEN
      IF (SSQ.LT.KRES) THEN
      FOSUM1=FOSUM1 + AFO1
      FC2SUM1=FC2SUM1 +WFC2
      FOCSUM1=FOCSUM1 +WFOC
      ELSE
      FOSUM2=FOSUM2 + AFO1
      FC2SUM2=FC2SUM2 +WFC2
      FOCSUM2=FOCSUM2 +WFOC
      END IF
      END IF
C
      END DO
C
      IF (FO2SUM.LT.RSMALL) THEN
      WRITE(6,'(A)')
     &  ' %XSCALE-error: sum over "observed" data is zero'
      ERROR=.TRUE.
      ELSE IF (FC2SUM.LT.RSMALL) THEN
      WRITE(6,'(A)')
     &  ' %XSCALE-error: sum over "calculated" data is zero'
      ERROR=.TRUE.
      ELSE
      ERROR=.FALSE.
      IF (XRFFKQ) THEN
      FFK=FOCSUM/FC2SUM
      END IF
      END IF
C
      IF (KRES.GT.RSMALL) THEN
      IF (FC2SUM1.LT.RSMALL) THEN
      ERROR=.TRUE.
      FFK1=ZERO
      WRITE(6,'(A,A)') ' %XSCALE-warning: sum over',
     &' the low resolution range is zero'
      ELSE
      ERROR=.FALSE.
      IF (QFFK1) THEN
      FFK1=FOCSUM1/FC2SUM1
      END IF
      END IF
C
      IF (FC2SUM2.LT.RSMALL) THEN
      ERROR=.TRUE.
      FFK2=ZERO
      WRITE(6,'(A,A)') ' %XSCALE-warning: sum over',
     &' the high resolution range is zero'
      ELSE
      ERROR=.FALSE.
      IF (QFFK2) THEN
      FFK2=FOCSUM2/FC2SUM2
      END IF
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCTARG(FFK,FOSUM,FO2SUM,KRES,FFK1,FFK2,
     &                   FOSUM1,FOSUM2,
     &                   ARG,PARAM,XRH,XRK,XRL,NNSELE,INDX,
     &                   MAXSET,TGSET,HPSFSET,TYPESF,
     &                   RVAL,TARG,RLOW,TLOW,RHIG,THIG,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   QANISO,QISO,MPERSET,ABCS,XRTR)
C
C Computes R values and least-squares targets.
C
C Target = sum{ W * ( |FTARG| - FFK * |FTOT| )^2 } / sum{ W * FTARG^2 }
C
C      R = sum{ ( |FTARG| - FFK * |FTOT| ) } / sum{ |FTARG| }
C
C where FTOT is the resultant over a number of SF data sets.
C
C For KRES > 0.0 two R values and two targets in two resolution
C ranges are also computed.
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER MAXSET
      DOUBLE PRECISION FFK, FOSUM, FO2SUM
      DOUBLE PRECISION KRES, FFK1, FFK2, FOSUM1, FOSUM2
      DOUBLE COMPLEX ARG(maxset)
      DOUBLE PRECISION PARAM(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER NNSELE
      INTEGER INDX(*)
      INTEGER TGSET
      INTEGER HPSFSET(*), TYPESF(*)
      DOUBLE PRECISION RVAL, TARG, RLOW, TLOW, RHIG, THIG
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), XRTR(3,3)
C local
      INTEGER REFLCT, IREF
      DOUBLE COMPLEX FTARG, FTOT
      DOUBLE PRECISION SSQ, AFTARG, AFTOT, DF0, DF2
C parameters
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
      TARG=ZERO
      RVAL=ZERO
      RLOW=ZERO
      TLOW=ZERO
      RHIG=ZERO
      THIG=ZERO
C
!$omp parallel do default(shared)
!$omp& private(iref,reflct,arg,ssq,ftarg,ftot,aftarg,aftot,df0,df2)
!$omp& reduction(+:rval,targ,rlow,tlow,rhig,thig)
!$omp& schedule(guided)
      DO IREF=1,NNSELE
      REFLCT=INDX(IREF)
C
      CALL XSCFCAL(REFLCT,ARG,PARAM,XRH,XRK,XRL,SSQ,
     &             MAXSET,TGSET,HPSFSET,TYPESF,FTARG,FTOT,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      AFTARG=SQRT(DBLE(FTARG)**2+DIMAG(FTARG)**2)
      AFTOT=SQRT(DBLE(FTOT)**2+DIMAG(FTOT)**2)
      DF0=ABS(AFTARG-FFK*AFTOT)
      DF2=DF0*DF0
      RVAL=RVAL+DF0
      TARG=TARG+DF2
C
      IF (KRES.GT.RSMALL) THEN
      IF (SSQ.LT.KRES) THEN
      DF0=ABS(AFTARG-FFK1*AFTOT)
      DF2=DF0*DF0
      RLOW=RLOW+DF0
      TLOW=TLOW+DF2
      ELSE
      DF0=ABS(AFTARG-FFK2*AFTOT)
      DF2=DF0*DF0
      RHIG=RHIG+DF0
      THIG=THIG+DF2
      END IF
      END IF
C
      END DO
C
      RVAL=RVAL/FOSUM
      TARG=TARG/FO2SUM
      IF (KRES.GT.RSMALL) THEN
      RLOW=RLOW/FOSUM1
      TLOW=TLOW/FO2SUM
      RHIG=RHIG/FOSUM2
      THIG=THIG/FO2SUM
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCFCAL(REFLCT,ARG,PARAM,XRH,XRK,XRL,SSQ,
     &                   MAXSET,TGSET,HPSFSET,TYPESF,FTARG,FTOT,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   QANISO,QISO,MPERSET,ABCS,XRTR)
C
C Routine calculates n-sets SF at given parameters
C FTARG is the target data set ('observed')
C FTOT is the resultant of number of SF data sets (MAXSET-1)
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER REFLCT
      DOUBLE PRECISION PARAM(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER MAXSET, TGSET
      INTEGER HPSFSET(*), TYPESF(*)
      DOUBLE COMPLEX ARG(*)
      DOUBLE COMPLEX FTARG, FTOT
      DOUBLE PRECISION SSQ
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), XRTR(3,3)
C local
      INTEGER HH, KK, LL, I, I1
      DOUBLE PRECISION FSCA, S1, S2, S3
C parameters
      DOUBLE PRECISION ZERO, TWO, QUART
      PARAMETER (ZERO=0.0D0,TWO=2.0D0,QUART=0.25D0)
C
C begin
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      CALL XRSSQ(HH,KK,LL,SSQ,XRTR)
      IF (QVSTACK) THEN
      CALL XSCSTACK(REFLCT,ARG,MAXSET,VLEVEL,VMAX,VSTACK,N,XEDNI)
      ELSE
      CALL XSCDRAW(REFLCT,ARG,MAXSET,HPSFSET,TYPESF)
      END IF
      FTOT=DCMPLX(ZERO,ZERO)
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      IF (QANISO) THEN
      S1=ABCS(1)*HH
      S2=ABCS(2)*KK
      S3=ABCS(3)*LL
      FSCA=PARAM(I1)*EXP(-(
     & PARAM(I1+1)*S1*S1+PARAM(I1+2)*S2*S2+PARAM(I1+3)*S3*S3+
     & TWO*(PARAM(I1+4)*S1*S2+PARAM(I1+5)*S1*S3+
     & PARAM(I1+6)*S2*S3))*QUART)
      ELSE
      FSCA=PARAM(I1)*EXP(-PARAM(I1+1)*SSQ*QUART)
      END IF
      IF (I.EQ.TGSET) THEN
      FTARG=-FSCA*ARG(I)
      ELSE
      FTOT=FTOT+FSCA*ARG(I)
      END IF
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCDERI(M,M1,REFLCT,ARG,PARAM,GDERI,XRH,XRK,XRL,
     &                   FIX,FIXMOD,FFK,KRES,FFK1,FFK2,
     &                   MAXSET,TGSET,HPSFSET,TYPESF,
     &                   QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &                   QANISO,QISO,MPERSET,ABCS,CABGS,RESTRC,
     &                   XRTR)
C
C Routine calculates derivatives of n-sets SF data
C with respect to parameters
C Derivatives store in GDERI
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER M, M1, REFLCT
      DOUBLE COMPLEX ARG(*)
      DOUBLE PRECISION PARAM(*), GDERI(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER FIX(*), FIXMOD
      DOUBLE PRECISION FFK, KRES, FFK1, FFK2
      INTEGER MAXSET, TGSET
      INTEGER HPSFSET(*), TYPESF(*)
      LOGICAL QVSTACK
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
      LOGICAL QANISO, QISO
      INTEGER MPERSET
      DOUBLE PRECISION ABCS(3), CABGS(3)
      CHARACTER*4 RESTRC
      DOUBLE PRECISION XRTR(3,3)
C local
      INTEGER I, I1, I2, J
      DOUBLE COMPLEX FTARG, FTOT
      DOUBLE PRECISION SSQ, EXPB, AFTARG, AFTOT, S1, S2, S3
      DOUBLE PRECISION TEMP, TEMPK
      DOUBLE PRECISION GRADK, GRADB
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, QUART
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, QUART=0.25D0)
C begin
C
C initialization
      DO I=1,M1
      GDERI(I)=ZERO
      END DO
C
C calculate F
      CALL XSCFCAL(REFLCT,ARG,PARAM,XRH,XRK,XRL,SSQ,
     &             MAXSET,TGSET,HPSFSET,TYPESF,FTARG,FTOT,
     &             QVSTACK,VLEVEL,VMAX,VSTACK,N,XEDNI,
     &             QANISO,QISO,MPERSET,ABCS,XRTR)
      AFTARG=SQRT(DBLE(FTARG)**2+DIMAG(FTARG)**2)
      AFTOT=SQRT(DBLE(FTOT)**2+DIMAG(FTOT)**2)
C
C difference between FTARG and FTOT
      TEMPK=FFK
      IF (KRES.GT.RSMALL) THEN
      IF (SSQ.LT.KRES) THEN
      TEMPK=FFK1
      ELSE
      TEMPK=FFK2
      END IF
      END IF
      GDERI(M1)=AFTARG-TEMPK*AFTOT
C
      IF (AFTOT.LT.RSMALL) THEN
      TEMP=ZERO
      ELSE
      TEMP=ONE/AFTOT
      END IF
C
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      AFTARG=TEMPK*TEMP*
     & (DBLE(FTOT)*DBLE(ARG(I))+DIMAG(FTOT)*DIMAG(ARG(I)))
C
C anisotropic
      IF (QANISO) THEN
      S1=ABCS(1)*XRH(REFLCT)
      S2=ABCS(2)*XRK(REFLCT)
      S3=ABCS(3)*XRL(REFLCT)
      EXPB=EXP(-(
     & PARAM(I1+1)*S1*S1+PARAM(I1+2)*S2*S2+PARAM(I1+3)*S3*S3+
     & TWO*(PARAM(I1+4)*S1*S2+PARAM(I1+5)*S1*S3+
     & PARAM(I1+6)*S2*S3))*QUART)
      GDERI(I1)=EXPB*AFTARG
      GDERI(I1+1)=-S1*S1*QUART*PARAM(I1)*GDERI(I1)
      GDERI(I1+2)=-S2*S2*QUART*PARAM(I1)*GDERI(I1)
      GDERI(I1+3)=-S3*S3*QUART*PARAM(I1)*GDERI(I1)
      GDERI(I1+4)=-TWO*S1*S2*QUART*PARAM(I1)*GDERI(I1)
      GDERI(I1+5)=-TWO*S1*S3*QUART*PARAM(I1)*GDERI(I1)
      GDERI(I1+6)=-TWO*S2*S3*QUART*PARAM(I1)*GDERI(I1)
C
C isotropic
      ELSE
      EXPB=EXP(-PARAM(I1+1)*SSQ*QUART)
      GDERI(I1)=EXPB*AFTARG
C      GDERI(I1+1)=-SSQ*QUART*PARAM(I1)*GDERI(I1)
C to allow "negative" B's - Fixed by JSJ  17-JAN-95
      GDERI(I1+1)=-SSQ*QUART*PARAM(I1)*GDERI(I1)*
     &             SIGN(ONE,PARAM(I1+1))
      END IF
      END DO
C
C impose symmetry restrictions according crystal systems
      IF (QANISO) THEN
      CALL XSCRSTRC(RESTRC,MAXSET,TGSET,MPERSET,GDERI,ABCS,CABGS,QISO)
      END IF
C
C set derivative to zero for fixed parameters
      DO I=1,M
      IF (FIX(I).NE.0) GDERI(I)=ZERO
      END DO
C
C gather derivatives (gradients)
      GRADK=ZERO
      GRADB=ZERO
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      IF (FIX(I1).EQ.0) GRADK=GRADK+GDERI(I1)
      DO J=1,MPERSET-1
      IF (FIX(I1+J).EQ.0) GRADB=GRADB+GDERI(I1+J)
      END DO
      END DO
C
C only one uniform scale factor if FIXMOD=3
      IF (FIXMOD.EQ.3) THEN
      I2=0
      DO I=1,MAXSET
      J=(I-1)*MPERSET+1
      IF (FIX(J).EQ.0.AND.I2.EQ.0) I2=J
      IF (FIX(J).EQ.0.AND.I2.GT.0) GDERI(J)=GRADK
      END DO
C only one uniform B factor if FIXMOD=4
      ELSE IF (FIXMOD.EQ.4) THEN
      I2=0
      DO I=1,MAXSET
      I1=(I-1)*MPERSET+1
      DO J=1,MPERSET-1
      IF (FIX(I1+J).EQ.0.AND.I2.EQ.0) I2=I1+J
      IF (FIX(I1+J).EQ.0.AND.I2.GT.0) GDERI(I1+J)=GRADB
      END DO
      END DO
      END IF
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCDRAW(REFLCT,ARG,MAXSET,HPSFSET,TYPESF)
C
C Routine draws SF property from a number of data sets
C for the same reflection and store in ARG
C All data types (1=DC, 2=DP, 3=IN) convert to double complex
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER REFLCT
      DOUBLE COMPLEX ARG(*)
      INTEGER MAXSET
      INTEGER HPSFSET(*), TYPESF(*)
C local
      INTEGER I
      DOUBLE COMPLEX DCOMP
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      DO I=1,MAXSET
      IF (TYPESF(I).EQ.1) THEN
      CALL XSCHPDC(HEAP(HPSFSET(I)),REFLCT,DCOMP)
      ELSEIF (TYPESF(I).EQ.2) THEN
      CALL XSCHPDP(HEAP(HPSFSET(I)),REFLCT,DCOMP)
      ELSEIF (TYPESF(I).EQ.3) THEN
      CALL XSCHPIN(HEAP(HPSFSET(I)),REFLCT,DCOMP)
      ELSE
      DCOMP=DCMPLX(ZERO,ZERO)
      END IF
      ARG(I)=DCOMP
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCHPDC(SFC,REF,DCOMP)
C
C Takes the element (REF) from the dummy DC array
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      DOUBLE COMPLEX SFC(*)
      DOUBLE COMPLEX DCOMP
      INTEGER REF
C
      DCOMP=SFC(REF)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCHPDP(SFD,REF,DCOMP)
C
C Takes the element (REF) from the dummy DP array
C and converts it to DC
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION SFD(*)
      DOUBLE COMPLEX DCOMP
      INTEGER REF
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      DCOMP=DCMPLX(SFD(REF),ZERO)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCHPIN(SFI,REF,DCOMP)
C
C Takes the element (REF) from the dummy INTEGER array
C and converts it to DC
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER SFI(*)
      DOUBLE COMPLEX DCOMP
      INTEGER REF
      DOUBLE PRECISION TEMP
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      TEMP=SFI(REF)
      DCOMP=DCMPLX(TEMP,ZERO)
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSCSTACK(REFLCT,ARG,MAXSET,VLEVEL,VMAX,VSTACK,N,XEDNI)
C
C Routine draws SF property from VSTACK
C for the same reflection and store in ARG
C All data types are double complex (DC)
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      INTEGER REFLCT
      DOUBLE COMPLEX ARG(*)
      INTEGER MAXSET
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER XEDNI(*)
C local
      INTEGER I, VL
C
      VL=VLEVEL-MAXSET+1
      DO I=1,MAXSET
      ARG(I)=VSTACK(XEDNI(REFLCT),VL)
      VL=VL+1
      END DO
C
      RETURN
      END
C===================================================================
      SUBROUTINE XRABCS(ABCS,CABGS,XRCELL,XRVOL)
C
C Compute reciprocal unit cell dimensions: a* , b*  and c* in ABCS,
C and cos(alpha*), cos(beta*) and cos(gamma*) in CABGS.
C
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      DOUBLE PRECISION ABCS(3), CABGS(3), XRCELL(6), XRVOL
C local
      DOUBLE PRECISION CABG1, CABG2, CABG3, SABG1, SABG2, SABG3
C parameter
      DOUBLE PRECISION ONE, TWO, RAD
      PARAMETER (ONE=1.0D0, TWO=2.0D0, RAD=PI/180.0D0)
C begin
      CABG1=COS(XRCELL(4)*RAD)
      CABG2=COS(XRCELL(5)*RAD)
      CABG3=COS(XRCELL(6)*RAD)
      SABG1=SIN(XRCELL(4)*RAD)
      SABG2=SIN(XRCELL(5)*RAD)
      SABG3=SIN(XRCELL(6)*RAD)
      CABGS(1)=(CABG2*CABG3-CABG1)/(SABG2*SABG3)
      CABGS(2)=(CABG3*CABG1-CABG2)/(SABG3*SABG1)
      CABGS(3)=(CABG1*CABG2-CABG3)/(SABG1*SABG2)
      XRVOL=XRCELL(1)*XRCELL(2)*XRCELL(3)*
     &                   SQRT(ONE+TWO*CABG1*CABG2*CABG3
     &                  -CABG1**2-CABG2**2-CABG3**2)
      ABCS(1)=XRCELL(2)*XRCELL(3)*SABG1/XRVOL
      ABCS(2)=XRCELL(1)*XRCELL(3)*SABG2/XRVOL
      ABCS(3)=XRCELL(1)*XRCELL(2)*SABG3/XRVOL
C
      RETURN
      END
C===================================================================
      SUBROUTINE XSCRSTRC(RESTRC,MAXSET,TGSET,MPERSET,PARAM,ABCS,CABGS,
     &                    QISO)
C
C 1. Impose symmetry restrictions according to the crystal system
C
C Triclinc         none
C Monoclinic       B13=B23=0 when beta=alpha=90
C                  B12=B23=0 when gamma=alpha=90
C                  B12=B13=0 when gamma=beta=90
C Orthorhombic     B12=B13=B23=0
C Tetragonal       B11=B22 and B12=B13=B23=0
C Rhombohedral     B11=B22=B33 and B12=B13=B23
C Hexagonal        B11=B22 and B13=B23=0
C Cubic            B11=B22=B33 and B12=B13=B23=0 (=isotropic)
C
C 2. Impose zero trace restriction if QISO is TRUE.
C
C
C Authors: J.-S. Jiang and Axel T. Brunger
C ========================================
C
      IMPLICIT NONE
C I/O
      CHARACTER*4 RESTRC
      INTEGER MAXSET, TGSET, MPERSET
      DOUBLE PRECISION PARAM(*)
      DOUBLE PRECISION ABCS(3), CABGS(3)
      LOGICAL QISO
C local
      INTEGER I, I1
      DOUBLE PRECISION TEMP
C parameter
      DOUBLE PRECISION ZERO, TWO, THREE, CSMALL
      PARAMETER (ZERO=0.0D0, TWO=2.0D0, THREE=3.0D0, CSMALL=0.001D0)
C begin
C
      DO I=1,MAXSET
      IF (I.NE.TGSET) THEN
      I1=(I-1)*MPERSET+1
C
      IF (RESTRC.EQ.'ALL ') THEN
C
      IF (ABCS(1).EQ.ABCS(2).AND.ABCS(2).EQ.ABCS(3)) THEN
      PARAM(I1+1)=(PARAM(I1+1)+PARAM(I1+2)+PARAM(I1+3))/THREE
      PARAM(I1+2)=PARAM(I1+1)
      PARAM(I1+3)=PARAM(I1+1)
      ELSE IF (ABCS(1).EQ.ABCS(2)) THEN
      PARAM(I1+1)=(PARAM(I1+1)+PARAM(I1+2))/TWO
      PARAM(I1+2)=PARAM(I1+1)
      ELSE IF (ABCS(2).EQ.ABCS(3)) THEN
      PARAM(I1+2)=(PARAM(I1+2)+PARAM(I1+3))/TWO
      PARAM(I1+3)=PARAM(I1+2)
      ELSE IF (ABCS(3).EQ.ABCS(1)) THEN
      PARAM(I1+3)=(PARAM(I1+3)+PARAM(I1+1))/TWO
      PARAM(I1+1)=PARAM(I1+3)
      END IF
C
      END IF
C
      IF (RESTRC.EQ.'ALL '.OR.RESTRC.EQ.'OFFD') THEN
C
      IF (ABS(CABGS(1)).LT.CSMALL) PARAM(I1+6)=ZERO
      IF (ABS(CABGS(2)).LT.CSMALL) PARAM(I1+5)=ZERO
      IF (ABS(CABGS(3)).LT.CSMALL) PARAM(I1+4)=ZERO
C
C identify Rhombohedral for restriction
      IF (ABCS(1).EQ.ABCS(2).AND.ABCS(2).EQ.ABCS(3).AND.
     &    CABGS(1).EQ.CABGS(2).AND.CABGS(2).EQ.CABGS(3).AND.
     &    ABS(CABGS(1)).GT.CSMALL) THEN
      PARAM(I1+4)=(PARAM(I1+4)+PARAM(I1+5)+PARAM(I1+6))/THREE
      PARAM(I1+5)=PARAM(I1+4)
      PARAM(I1+6)=PARAM(I1+4)
      END IF
C
      END IF
C
C
C zero trace restriction.
      IF (.NOT.QISO) THEN
      TEMP=(PARAM(I1+1)+PARAM(I1+2)+PARAM(I1+3))/THREE
      PARAM(I1+1)=PARAM(I1+1)-TEMP
      PARAM(I1+2)=PARAM(I1+2)-TEMP
      PARAM(I1+3)=PARAM(I1+3)-TEMP
      END IF
C
      END IF
      END DO
C
      RETURN
      END
