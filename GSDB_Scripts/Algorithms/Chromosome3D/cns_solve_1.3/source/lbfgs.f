      SUBROUTINE LBFGS
C
C     Limited-memory BFGS minimizer
C
C     Uses code from D.C.Liu and J.Nocedal (nocedal@ece.nwu.edu)
C     "On the limited memory BFGS method for large scale optimization",
C     D. Liu and J. Nocedal, Mathematical Programming B 45 (1989) 503-528.
C
C     Initial implementation by Greg McMullan, HPCF
C     University of Cambridge, 30/3/99
C
C     Modified by Paul Adams, 04/99
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'lbfgs.inc'
      INCLUDE 'funct.inc'
C local
      DOUBLE PRECISION ACC, TOLG, TARGET, DFPRED
      INTEGER NDIM, I, MAXFN
      INTEGER NPRINT, M, IW, MAXFEV
C pointer
      INTEGER XCURR, GCURR, DIAG, IWORK, BEST
C external
      EXTERNAL ENFUNC
C begin
C
C defaults
      MAXFN=10
      TOLG=0.0001D0
      DFPRED=0.001D0
      NPRINT=1
      MAXFEV=10
C
C parsing
      CALL PUSEND('LBFGS>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('LBFGS>')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-minimize-lbfgs')
C
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('NSTEp=',MAXFN)
      ELSE IF (WD(1:4).EQ.'STEP'.OR.WD(1:4).EQ.'DROP') THEN
      CALL NEXTF('DROP=',DFPRED)
      ELSE IF (WD(1:4).EQ.'TOLG') THEN
      CALL NEXTF('TOLGradient=',TOLG)
      ELSE IF (WD(1:4).EQ.'NPRI') THEN
      CALL NEXTI('NPRInt=',NPRINT)
      ELSE IF (WD(1:4).EQ.'MAXF') THEN
      CALL NEXTI('MAXFev=',MAXFEV)
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,1200) MAXFN, TOLG,MAXFEV
1200  FORMAT(/' LBFGS-method :  NSTEp=',I7,
     1        ' TOLGradient=',F12.9,' MAXFevalulations=',I4,/)
      ELSE
      CALL CHKEND('LBFGS>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
C NDIM is the dimension of the problem
C
      NDIM=0
      DO I=1,NATOM
      IF (IMOVE(I).EQ.0) NDIM=NDIM+3
      END DO
      IF (NDIM.EQ.0) THEN
      WRITE(6,'(A)') ' %LBFGS-ERR: all atoms fixed. No action'
      ELSE
C
      WRITE(6,'(A,I6)')
     &  ' LBFGS: number of degrees of freedom=',NDIM
C
      ACC=TOLG**2 *NDIM
C
C dynamic allocation
      XCURR = ALLHP(IREAL8(NDIM))
      GCURR = ALLHP(IREAL8(NDIM))
      DIAG  = ALLHP(IREAL8(NDIM))
      BEST  = ALLHP(IREAL8(NDIM))
C
      M = 6
      IW = NDIM*(2*M+1) + 2*M
      IWORK = ALLHP(IREAL8(IW))
C
      CALL ZXLBFGS(ENFUNC, NDIM, M, HEAP(XCURR), HEAP(GCURR), TARGET,
     &             .FALSE., HEAP(DIAG), NPRINT, ACC, MAXFN, HEAP(IWORK),
     &             HEAP(BEST),MAXFEV)
C
C transfer coordinates from best minimum
      CALL PAR2XYZ(NATOM,HEAP(BEST))
C
C free space
      CALL FREHP(IWORK,IREAL8(IW))
      CALL FREHP(BEST,IREAL8(NDIM))
      CALL FREHP(DIAG,IREAL8(NDIM))
      CALL FREHP(GCURR,IREAL8(NDIM))
      CALL FREHP(XCURR,IREAL8(NDIM))
C
      END IF
C
      RETURN
      END
C
C =================================================================
C
      SUBROUTINE ENFUNC(N,XC,F,G,ITER,QPRINT)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'ener.inc'
C I/O
      INTEGER N, ITER
      DOUBLE PRECISION F, XC(*), G(*)
      LOGICAL QPRINT
C local
      INTEGER IC, I
C
      IF (QPRINT) THEN
         WRITE(6,'(A,I6,A)')
     &        ' ------------------------------ cycle=',ITER,
     &        ' -----------------------------------'
         CALL PRINTE(RENR)
         WRITE(6,'(2A)')
     &        ' --------------------------------------------',
     &        '-----------------------------------'
C
      ELSE
C
         IF (ITER.EQ.0) THEN
C
            IC=1
            DO I=1,NATOM
               IF (IMOVE(I) .EQ. 0 ) THEN
                  XC(IC)   = X(I)
                  XC(IC+1) = Y(I)
                  XC(IC+2) = Z(I)
                  IC = IC + 3
               END IF
            END DO
C
         ELSE
C
            IC=1
            DO I=1,NATOM
               IF (IMOVE(I) .EQ. 0 ) THEN
                  X(I) = XC(IC)
                  Y(I) = XC(IC+1)
                  Z(I) = XC(IC+2)
                  IC = IC + 3
               END IF
            END DO
C
         END IF
C
         DO I=1,NATOM
            IF (IMOVE(I).NE.0) THEN
               DX(I) = 0.0D0
               DY(I) = 0.0D0
               DZ(I) = 0.0D0
            END IF
         END DO
C
         CALL ENERGY
         F = RENR(SSENER)
C
         IC=1
         DO I=1,NATOM
            IF (IMOVE(I).EQ.0) THEN
               G(IC)   = DX(I)
               G(IC+1) = DY(I)
               G(IC+2) = DZ(I)
               IC = IC + 3
            END IF
         END DO
C
      END IF
C
      RETURN
      END
C
C =================================================================
C
      SUBROUTINE PAR2XYZ(N,PARAM)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
C I/O
      INTEGER N
      DOUBLE PRECISION PARAM(*)
C local
      INTEGER IC, I
C
      IC=1
      DO I=1,N
         IF (IMOVE(I) .EQ. 0 ) THEN
            X(I) = PARAM(IC)
            Y(I) = PARAM(IC+1)
            Z(I) = PARAM(IC+2)
            IC = IC + 3
         END IF
      END DO
C
      RETURN
      END
C
C =================================================================
C
      SUBROUTINE LBFGSINI
C
      IMPLICIT NONE
C
      INCLUDE 'lbfgs.inc'
C
      GTOL=9.0D-01
      STPMIN=1.0D-20
      STPMAX=1.0D20
C
      RETURN
      END
C
C =================================================================
C
      SUBROUTINE ZXLBFGS(FUNCT,N,M,XC,G,F,DIAGCO,DIAG,NPRINT,EPS,
     &                   MAXFN,W,BEST,MAXFEV)
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C     This subroutine solves the unconstrained minimization problem
C
C                      min F(x),    x= (x1,x2,...,xN),
C
C      using the limited memory BFGS method. The routine is especially
C      effective on problems involving a large number of variables. In
C      a typical iteration of this method an approximation Hk to the
C      inverse of the Hessian is obtained by applying M BFGS updates to
C      a diagonal matrix Hk0, using information from the previous M steps.
C      The user specifies the number M, which determines the amount of
C      storage required by the routine. The user may also provide the
C      diagonal matrices Hk0 if not satisfied with the default choice.
C      The algorithm is described in "On the limited memory BFGS method
C      for large scale optimization", by D. Liu and J. Nocedal,
C      Mathematical Programming B 45 (1989) 503-528.
C
C Adapted by Paul Adams and Greg McMullan for use in CNS.
C Reverse calling to highest level removed to facilitate use as a
C general minimizer (see how ZXCGR is used for example).
C
C Paul Adams 04/99
C
      IMPLICIT NONE
C
      INCLUDE 'machvar.inc'
      INCLUDE 'lbfgs.inc'
C
C I/O
      INTEGER N, M, NPRINT, MAXFN
      DOUBLE PRECISION XC(N), G(N), DIAG(N), BEST(N)
      DOUBLE PRECISION W(N*(2*M+1)+2*M)
      DOUBLE PRECISION F, EPS
      LOGICAL DIAGCO
      EXTERNAL FUNCT
      INTEGER MAXFEV
C
C local
      DOUBLE PRECISION GNORM,DDOT,STP1,FTOL,YS,YY,SQ,YR,BETA,XNORM
      DOUBLE PRECISION STP,FBEST
      INTEGER NFUN,POINT,ISPT,IYPT,INFO
      INTEGER BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,ITER,IFLAG
      LOGICAL FINISH, QPRINT
C
C parameter
      DOUBLE PRECISION ONE, ZERO
      PARAMETER (ONE=1.0D0, ZERO=0.0D0)
C
C initialize some parameters
C
      CALL LBFGSINI
C
      ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      NFUN=1
      POINT=0
      FINISH=.FALSE.
      QPRINT=.FALSE.
C
C     PARAMETERS FOR LINE SEARCH ROUTINE
C
      FTOL= 1.0D-4
C
      CALL FUNCT(N,XC,F,G,ITER,QPRINT)
C
      FBEST=F
      DO I=1,N
         BEST(I)=XC(I)
      END DO
C
      IF(DIAGCO) THEN
         DO I=1,N
            IF (DIAG(I).LE.ZERO) GO TO 195
         END DO
      ELSE
         DO I=1,N
            DIAG(I)= 1.0D0
         END DO
      ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
      ISPT= N+2*M
      IYPT= ISPT+N*M
      DO I=1,N
         W(ISPT+I)= -G(I)*DIAG(I)
      END DO
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      STP1= ONE/GNORM
C
C    --------------------
C     MAIN ITERATION LOOP
C    --------------------
C
 80   CONTINUE
      ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M) BOUND=M
C
      YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO I=1,N
            DIAG(I)= YS/YY
         END DO
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(DIAGCO) THEN
        DO I=1,N
           IF (DIAG(I).LE.ZERO) GO TO 195
        END DO
      ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO I=1,N
         W(I)= -G(I)
      END DO
C
      CP= POINT
      DO I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
      END DO
C
      DO I=1,N
         W(I)=DIAG(I)*W(I)
      END DO
C
      DO I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
      END DO
C
C     STORE THE NEW SEARCH DIRECTION
C     ------------------------------
C
       DO I=1,N
          W(ISPT+POINT*N+I)= W(I)
       END DO
C
C     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION
C     BY USING THE LINE SEARCH ROUTINE MCSRCH
C     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO I=1,N
         W(I)=G(I)
      END DO
 172  CONTINUE
C
      CALL MCSRCH(N,XC,F,G,W(ISPT+POINT*N+1),STP,FTOL,
     *            FPEPS,MAXFEV,INFO,NFEV,DIAG)
C
      IF (INFO .EQ. -1) THEN
C
         CALL FUNCT(N,XC,F,G,ITER,QPRINT)
C
         IFLAG=1
         GOTO 172
      ELSE
C
         IF(NPRINT.GT.0) THEN
            IF(MOD((ITER),NPRINT).EQ.0) THEN
               QPRINT=.TRUE.
               CALL FUNCT(N,XC,F,G,ITER,QPRINT)
               QPRINT=.FALSE.
            END IF
         END IF
C
      ENDIF
C
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
C
C     COMPUTE THE NEW STEP AND GRADIENT CHANGE
C     -----------------------------------------
C
      NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      END DO
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
C
C     TERMINATION TEST
C     ----------------
C
      GNORM= DSQRT(DDOT(N,G,1,G,1))
      XNORM= DSQRT(DDOT(N,XC,1,XC,1))
      XNORM= DMAX1(1.0D0,XNORM)
      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
C
C is this the best minimum so far?
      IF (F.LT.FBEST) THEN
         DO I=1,N
            BEST(I)=XC(I)
         END DO
      END IF
C
C print information
C
      IF (ITER.GE.MAXFN) THEN
      WRITE(6,'(A)')' LBFGS: normal termination - NSTEP limit reached'
      RETURN
      END IF
C
      IF (FINISH) THEN
      WRITE(6,'(A)')' LBFGS: normal termination - gradient converged'
      RETURN
      ENDIF
C
      GO TO 80
C
C     ------------------------------------------------------------
C     END OF MAIN ITERATION LOOP. ERROR EXITS.
C     ------------------------------------------------------------
C
C line search abandonded
C
 190  IFLAG=-1
      WRITE(6,'(A)')' LBFGS: Line search terminated'
      IF (INFO.EQ.0) THEN
      WRITE (6,'(A)')
     &  ' %LBFGS-ERR: Coding error - improper input parameters'
      ELSEIF (INFO.EQ.2) THEN
      WRITE (6,'(A)')
     &  ' LBFGS: interval of uncertainty equals machine precision'
      ELSEIF (INFO.EQ.3) THEN
      WRITE (6,'(A)')
     &  ' LBFGS: maximum number of function evaluations reached'
      ELSEIF (INFO.EQ.4) THEN
      WRITE (6,'(A)')
     &  ' LBFGS: step size has reached lower bound'
      ELSEIF (INFO.EQ.5) THEN
      WRITE (6,'(A)')
     &  ' LBFGS: step size has reached upper bound'
      ELSEIF (INFO.EQ.6) THEN
      WRITE (6,'(A)')
     &  ' LBFGS: rounding errors prevent further progress'
      END IF
      RETURN
C
C diagonal matrix problems
C
 195  IFLAG=-2
      WRITE(6,'(2A)')
     &  ' %LBFGS-ERR: The ',I,'-th diagonal element of the inverse',
     &  ' Hessian approximation is not positive'
      RETURN
C
C input parameters incorrect
C
 196  IFLAG= -3
      WRITE(6,'(A)')
     &  ' %LBFGS-ERR: Coding Error - improper input parameters'
C
      RETURN
      END
C    ------------------------------------------------------------------
C
C     **************************
C     LINE SEARCH ROUTINE MCSRCH
C     **************************
C
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'lbfgs.inc'
C I/O
      INTEGER N,MAXFEV,INFO,NFEV
      DOUBLE PRECISION F,STP,FTOL,XTOL
      DOUBLE PRECISION X(N),G(N),S(N),WA(N)
C
C                     SUBROUTINE MCSRCH
C
C     A slight modification of the subroutine CSRCH of More' and Thuente.
C     The changes are to allow reverse communication, and do not affect
C     the performance of the routine.
C
C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
C
C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
C     MINIMIZER OF THE MODIFIED FUNCTION
C
C          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
C
C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
C     CONTAINS A MINIMIZER OF F(X+STP*S).
C
C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
C     THE SUFFICIENT DECREASE CONDITION
C
C           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
C
C     AND THE CURVATURE CONDITION
C
C           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
C
C     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
C
C     THE SUBROUTINE STATEMENT IS
C
C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
C         X + STP*S.
C
C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
C
C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
C         OF F AT X + STP*S.
C
C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
C         SEARCH DIRECTION.
C
C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
C         STP CONTAINS THE FINAL ESTIMATE.
C
C       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
C         communication implementation GTOL is defined in a COMMON
C         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
C         SATISFIED.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C         IS AT MOST XTOL.
C
C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
C         communication implementatin they are defined in a COMMON
C         statement).
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
C         MAXFEV BY THE END OF AN ITERATION.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C
C         INFO = 0  IMPROPER INPUT PARAMETERS.
C
C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
C
C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
C
C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
C                   IS AT MOST XTOL.
C
C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
C
C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
C
C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
C
C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
C                   TOLERANCES MAY BE TOO SMALL.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       WA IS A WORK ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MCSTEP
C
C       FORTRAN-SUPPLIED...ABS,MAX,MIN
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
C     **********
C local
      INTEGER J
      DOUBLE PRECISION DGM,DGXM,DGYM,FM,FXM,FYM
C
C parameter
      DOUBLE PRECISION P5, P66, XTRAPF, ZERO
      PARAMETER (P5=0.5D0,P66=0.66D0,XTRAPF=4.0D0,ZERO=0.0D0)
C
C begin
C
      IF (INFO.EQ.-1) GO TO 45
      INFOC = 1
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FTOL .LT. ZERO .OR.
     *    GTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO
     *    .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
C
C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
C     AND CHECK THAT S IS A DESCENT DIRECTION.
C
      DGINIT = ZERO
      DO J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
      END DO
      IF (DGINIT .GE. ZERO) then
         WRITE (PUNIT,15)
 15      FORMAT(/' The search direction is not a descent direction')
         RETURN
      ENDIF
C
C     INITIALIZE LOCAL VARIABLES.
C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5
      DO J = 1, N
         WA(J) = X(J)
      END DO
C
C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
C     THE INTERVAL OF UNCERTAINTY.
C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
C
C     START OF ITERATION.
C
   30 CONTINUE
C
C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
         END IF
C
C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
C
C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
C        STP BE THE LOWEST POINT OBTAINED SO FAR.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0
     *      .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX
C
C        EVALUATE THE FUNCTION AND GRADIENT AT STP
C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
C        We return to main program to obtain F and G.
C
         DO J = 1, N
            X(J) = WA(J) + STP*S(J)
         END DO
         INFO=-1
         RETURN
C
   45    CONTINUE
C
         INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO J = 1, N
            DG = DG + G(J)*S(J)
         END DO
         FTEST1 = FINIT + STP*DGTEST
C
C        TEST FOR CONVERGENCE.
C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX))
     *      .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND.
     *       F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND.
     *       (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1
C
C        CHECK FOR TERMINATION.
C
         IF (INFO .NE. 0) RETURN
C
C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND.
     *       DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
C
C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
C
C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
C
C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,
     *                  BRACKT,STMIN,STMAX,INFOC)
C
C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
C
C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
C           AND TO COMPUTE THE NEW STEP.
C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,
     *                  BRACKT,STMIN,STMAX,INFOC)
         END IF
C
C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
C        INTERVAL OF UNCERTAINTY.
C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1)
     *         STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
         END IF
C
C        END OF ITERATION.
C
         GO TO 30
C
C     LAST LINE OF SUBROUTINE MCSRCH.
C
      END
C
C ====================================================================
C
      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
     *                  STPMIN,STPMAX,INFO)
C
      IMPLICIT NONE
C I/O
      INTEGER INFO
      DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
C
C     SUBROUTINE MCSTEP
C
C     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
C     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
C     A MINIMIZER OF THE FUNCTION.
C
C     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
C     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
C     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
C     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
C     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
C     WITH ENDPOINTS STX AND STY.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
C                        STPMIN,STPMAX,INFO)
C
C     WHERE
C
C       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
C         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
C         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
C         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
C
C       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
C         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
C         UPDATED APPROPRIATELY.
C
C       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
C         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
C         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
C         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
C
C       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
C         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
C         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
C         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
C
C       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
C         AND UPPER BOUNDS FOR THE STEP.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
C         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
C         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
C         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
C     JORGE J. MORE', DAVID J. THUENTE
C
      DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR.
     *    STP .GE. MAX(STX,STY))) .OR.
     *    DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
C
C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
C
      SGND = DP*(DX/ABS(DX))
C
C     FIRST CASE. A HIGHER FUNCTION VALUE.
C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
C
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
            STPF = STPC + (STPQ - STPC)/2
         END IF
         BRACKT = .TRUE.
C
C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
C
      ELSE IF (SGND .LT. 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = ((GAMMA - DP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
         END IF
         BRACKT = .TRUE.
C
C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
C
      ELSE IF (ABS(DP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DP
         S = MAX(ABS(THETA),ABS(DX),ABS(DP))
C
C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
C        TO INFINITY IN THE DIRECTION OF THE STEP.
C
         GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DP) + THETA
         Q = (GAMMA + (DX - DP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
         END IF
         STPQ = STP + (DP/(DP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
            END IF
         END IF
C
C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
C
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
         END IF
      END IF
C
C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
C
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
         END IF
         STX = STP
         FX = FP
         DX = DP
      END IF
C
C     COMPUTE THE NEW STEP AND SAFEGUARD IT.
C
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
         END IF
      END IF
C
      RETURN
C
C     LAST LINE OF SUBROUTINE MCSTEP.
C
      END
C
C   ----------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     Forms the dot product of two vectors.
C     Uses unrolled loops for increments equal to one.
C     Jack Dongarra, LINPACK, 3/11/78.
C
      IMPLICIT NONE
C I/O
      INTEGER N
      INTEGER INCX, INCY
      DOUBLE PRECISION DX(*), DY(*)
C local
      DOUBLE PRECISION DTEMP
      INTEGER I, IX, IY, M, MP1
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C code for unequal increments or equal increments
C not equal to 1
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      DDOT = DTEMP
      RETURN
C
C code for both increments equal to 1
C
C clean-up loop
C
   20 M = MOD(N,5)
      IF ( M .EQ. 0 ) GO TO 40
      DO I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
      END DO
      IF ( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO I = MP1,N,5
         DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     &                   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) +
     &                   DX(I + 4)*DY(I + 4)
      END DO
   60 DDOT = DTEMP
      RETURN
      END
C
C =======================================================================
C
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C Constant times a vector plus a vector.
C Uses unrolled loops for increments equal to one.
C Jack Dongarra, LINPACK, 3/11/78.
C
      IMPLICIT NONE
C I/O
      INTEGER N, INCX, INCY
      DOUBLE PRECISION DX(*), DY(*), DA
C local
      INTEGER I,IX,IY,M,MP1
C
      IF (N.LE.0) RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C code for unequal increments or equal increments
C not equal to 1
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO I = 1,N
         DY(IY) = DY(IY) + DA*DX(IX)
         IX = IX + INCX
         IY = IY + INCY
      END DO
C
      RETURN
C
C code for both increments equal to 1
C
C clean-up loop
C
   20 M = MOD(N,4)
      IF ( M .EQ. 0 ) GO TO 40
      DO I = 1,M
         DY(I) = DY(I) + DA*DX(I)
      END DO
      IF ( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO I = MP1,N,4
         DY(I)     = DY(I)     + DA*DX(I)
         DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
         DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
         DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
      END DO
C
      RETURN
      END
C
C ======================================================================
C
