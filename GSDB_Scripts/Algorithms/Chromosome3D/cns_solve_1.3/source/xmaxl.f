! OpenMP KD 2005-2006
      SUBROUTINE XDOMLF(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the maximum likelihood value
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C  mlf(Fobs, sigmaF, Fcalc, D, sigma_delta) from script level
C
C Fobs: real
C sigmaF: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
C arguments: Fobs, sigmaF, Fcalc, D, sigma_delta
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDOMLF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOMLF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
c OpenMP by Kay Diederichs
C
C Copyright 1996  University of Alberta and Yale University
C
C For non-centric reflections
C
C  <Fo> = sqrt(Pi)*sigd*sqrt(e)/2 * 1F1(-1/2,1,-(DFc/(sigd*sqrt(e))**2)
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C For centric reflections,
C
C  <Fo> = sqrt(2/Pi)*sigd*sqrt(e) *
C          1F1(-1/2,1/2,-(DFc)**2/(2*(sigd*sqrt(e))**2))
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C The target is the minus log likelihood:
C
C   -LL = ln(sqrt(variance)) + (<Fo>-Fo)^2/(2*variance)
C       = 1/2 * [ln(variance) + (<Fo>-Fo)^2/variance]
C
      IMPLICIT NONE
C
      INCLUDE 'checof.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION FOBS, SIGMA, FCALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION EFO, VARML, DELTA
      DOUBLE PRECISION ARGSQ, ARGSQINV, SIGDE, TEMP
      DOUBLE PRECISION SQPI2, SQ2PI
C functions
      DOUBLE PRECISION CHEVAL
C parameter
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, THREE, FOUR
      DOUBLE PRECISION EIGHT, SIXTEEN, TWENTY
      DOUBLE PRECISION TWOFIVE, FOURNINE, TWOSEVEN
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (THREE=3.0D0, FOUR=4.0D0, EIGHT=8.0D0)
      PARAMETER (SIXTEEN=16.0D0, TWENTY=20.0D0, TWOFIVE=25.0D0)
      PARAMETER (FOURNINE=49.0D0, TWOSEVEN=27.0D0)
C
C begin
      SQPI2=HALF*SQRT(PI)
      SQ2PI=SQRT(TWO/PI)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
!$omp parallel do 
!$omp& private(i,fobs,sigma,fcalc,d,sigmad,epsilon,sigde,argsq,argsqinv,
!$omp&         efo,delta,varml,temp)
!$omp& shared(n,vstack,vlevel,type,qherm,xrnsym,mult,index,nchc,wchc,
!$omp&        chco,sqpi2,sq2pi)
!$omp& default(none)
      DO I=1,N
      FOBS=DBLE(VSTACK(I,VLEVEL))
      SIGMA=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C trap SIGMAD equal zero
      IF (SIGDE.GT.ZERO) THEN
C
C is this an acentric reflection
         IF (TYPE(INDEX(I)).GE.1) THEN
            ARGSQ = (D*FCALC/SIGDE)**2
            IF (ARGSQ.GT.ZERO) THEN
               ARGSQINV = ONE/ARGSQ
            END IF
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use the series approximation.
C The error in this approximation is less than 1 in 10^7 for
C DFc/sigd*sqrt(e) > 3.4
            IF (ARGSQ.LT.(-WCHC(1))) THEN
               EFO = SQPI2*SIGDE*CHEVAL(NCHC(1),WCHC(1),
     &                                  CHCO(1,1),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               DELTA = D*FCALC*        (ARGSQINV/FOUR   )*
     &                 (ONE +          (ARGSQINV/EIGHT  )*
     &                 (ONE +    THREE*(ARGSQINV/FOUR   )*
     &                 (ONE +  TWOFIVE*(ARGSQINV/SIXTEEN)*
     &                 (ONE + FOURNINE*(ARGSQINV/TWENTY )*
     &                 (ONE + TWOSEVEN*(ARGSQINV/EIGHT  ))))))
               EFO = D*FCALC+DELTA
            END IF
C
C is this a centric reflection
         ELSE
            ARGSQ = ((D*FCALC/SIGDE)**2)/TWO
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use DFc. For DFc/sigd*sqrt(e) > 4.72 <Fo> ~ DFc within
C 1 part in 10^7.
            IF (ARGSQ.LT.(-WCHC(3))) THEN
               EFO = SQ2PI*SIGDE*CHEVAL(NCHC(3),WCHC(3),
     &                                  CHCO(1,3),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               EFO = D*FCALC
               DELTA = ZERO
            END IF
         END IF
      ELSE
C
C if sigma_delta is equal to zero
C for both centric and acentric
         EFO = D*FCALC
         DELTA = ZERO
      END IF
C
C calculate variance
      VARML = -DELTA*(D*FCALC+EFO)+SIGDE**2+SIGMA**2
C
C our target is the log likelihood
      TEMP=HALF*(LOG(VARML)+(EFO-FOBS)**2/VARML)
C
      VSTACK(I,VLEVEL) = DCMPLX(TEMP,ZERO)
C
      END DO
C
      RETURN
      END
C
C ==================================================================
C
      SUBROUTINE XDODMLF(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the maximum likelihood derivatives
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C   dmlf(Fobs, sigmaF, Fcalc, D, sigma_delta) from script level
C
C arguments Fobs, sigmaF, Fcalc, D, sigma_delta
C Fobs: real
C sigmaF: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDODMLF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C ===================================================================
C
      SUBROUTINE XDODMLF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
c OpenMP by Kay Diederichs
C
C Copyright 1996  University of Alberta and Yale University
C
C this routine calculates the derivates, but must also calculate
C <Fo> and variance, thus duplicating some of the work done in
C routine XDOMLF
C
C For non-centric reflections
C
C   <Fo> = sqrt(Pi)*sigd*sqrt(e)/2 * 1F1(-1/2,1,-(DFc/(sigd*sqrt(e))**2)
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C   d<Fo>/d(Fc)*Fc = D**2/(2*sigd*sqrt(e)) *
C                    sqrt(Pi) * 1F1(1/2,2,-(DFc/(sigd*sqrt(e))**2)
C
C   d(variance)/d(Fc)*Fc = (DFc+<Fo>)(D/Fc-(d<Fo>/d(Fc)*Fc)) +
C                          (DFc-<Fo>)(D/Fc+(d<Fo>/d(Fc)*Fc))
C
C For centric reflections,
C
C   <Fo> = sqrt(2/Pi)*sigd*sqrt(e) *
C          1F1(-1/2,1/2,-(DFc)**2/(2*(sigd*sqrt(e))**2))
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C   d<Fo>/d(Fc)*Fc = sqrt(2/Pi)*(D**2/sigd*sqrt(e)) *
C                    1F1(1/2,3/2,-(DFc)**2/(2*(sigd*sqrt(e))**2))
C
C   d(variance)/d(Fc)*Fc = 2*D**2 - 2*<Fo> * d<Fo>/d(Fc)*Fc
C
C The derivatives are returned:
C
C   1/variance * [(<Fo>-Fo)*(d<Fo>/d(Fc)*Fc) +
C                 ((1-(<Fo>-Fo)^2)/(2*variance)) * (d(variance)/d(Fc)*Fc)]
C
C N.B. The partial derivatives all include a factor of Fc,
C to avoid divide by zero problems they have been divided
C by Fc symbolically.
C
      IMPLICIT NONE
C
      INCLUDE 'checof.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION FOBS, SIGMA, FCALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION EFO, VARML, DELTA, DELTAP
      DOUBLE PRECISION ARGSQ, ARGSQINV, SIGDE, PDER, PDER2
      DOUBLE PRECISION SQPI2, SQ2PI
C functions
      DOUBLE PRECISION CHEVAL
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, FIVE
      DOUBLE PRECISION EIGHT, SIXTEEN, TWENTY, HALF
      DOUBLE PRECISION TWOFIVE, FOURNINE, TWOSEVEN
      DOUBLE PRECISION THREEFIVE, SIXTHREE, THREETHREE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, FIVE=5.0D0, EIGHT=8.0D0, SIXTEEN=16.0D0)
      PARAMETER (TWENTY=20.0D0, TWOFIVE=25.0D0, FOURNINE=49.0D0)
      PARAMETER (TWOSEVEN=27.0D0,THREEFIVE=35.0D0, SIXTHREE=63.0D0)
      PARAMETER (THREETHREE=33.0D0, HALF=0.5D0)
C
C begin
      SQPI2=HALF*SQRT(PI)
      SQ2PI=SQRT(TWO/PI)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
!$omp parallel do 
!$omp& private(i,fobs,sigma,fcalc,d,sigmad,epsilon,sigde,argsq,argsqinv,
!$omp&         efo,delta,deltap,pder,pder2,varml)
!$omp& shared(n,vstack,vlevel,type,qherm,xrnsym,mult,index,nchc,wchc,
!$omp&        chco,sqpi2,sq2pi)
!$omp& default(none)
      DO I=1,N
      FOBS=DBLE(VSTACK(I,VLEVEL))
      SIGMA=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C trap SIGMAD equal zero
      IF (SIGDE.GT.ZERO) THEN
C
C is this an acentric reflection
         IF (TYPE(INDEX(I)).GE.1) THEN
            ARGSQ = (D*FCALC/SIGDE)**2
            IF (ARGSQ.GT.ZERO) THEN
               ARGSQINV = ONE/ARGSQ
            END IF
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use the series approximation.
C The error in this approximation is less than 1 in 10^7 for
C DFc/sigd*sqrt(e) > 3.4
            IF (ARGSQ.LT.(-WCHC(1))) THEN
               EFO = SQPI2*SIGDE*CHEVAL(NCHC(1),WCHC(1),
     &                                  CHCO(1,1),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               DELTA = D*FCALC        *(ARGSQINV/FOUR   )*
     &                 (ONE +          (ARGSQINV/EIGHT  )*
     &                 (ONE +    THREE*(ARGSQINV/FOUR   )*
     &                 (ONE +  TWOFIVE*(ARGSQINV/SIXTEEN)*
     &                 (ONE + FOURNINE*(ARGSQINV/TWENTY )*
     &                 (ONE + TWOSEVEN*(ARGSQINV/EIGHT  ))))))
               EFO = D*FCALC+DELTA
            END IF
C
C calculate derivatives
C again use derivative of the Chebychev polynomial approximation
C to 1F1 for low values of ARGSQ, else use the derivative of the
C series approximation
            IF (ARGSQ.LT.(-WCHC(2))) THEN
               PDER = SQPI2*D*D/SIGDE*
     &                CHEVAL(NCHC(2),WCHC(2),CHCO(1,2),-ARGSQ)
               PDER2 = TWO*(D**2 - EFO*PDER)
            ELSE
               DELTAP = (D/FCALC)*        (-(ARGSQINV/FOUR   )*
     &                  (ONE +      THREE*(  ARGSQINV/EIGHT  )*
     &                  (ONE +       FIVE*(  ARGSQINV/FOUR   )*
     &                  (ONE +  THREEFIVE*(  ARGSQINV/SIXTEEN)*
     &                  (ONE +   SIXTHREE*(  ARGSQINV/TWENTY )*
     &                  (ONE + THREETHREE*(  ARGSQINV/EIGHT  )))))))
               PDER = D/FCALC + DELTAP
               PDER2 = -DELTAP * (D*FCALC+EFO) -
     &                  DELTA * (D/FCALC+PDER)
            END IF
C
C is this a centric reflection
         ELSE
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use DFc. For DFc/sigd*sqrt(e) > 4.72 <Fo> ~ DFc within
C 1 part in 10^7.
            ARGSQ = ((D*FCALC/SIGDE)**2)/TWO
            IF (ARGSQ.LT.(-WCHC(3))) THEN
               EFO = SQ2PI*SIGDE*CHEVAL(NCHC(3),WCHC(3),
     &                                  CHCO(1,3),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               EFO = D*FCALC
               DELTA = ZERO
            END IF
C
C calculate derivatives
C again use derivative of the Chebychev polynomial approximation
C to 1F1 for low values of ARGSQ, else use approximation.
            IF (ARGSQ.LT.(-WCHC(4))) THEN
               PDER = SQ2PI*D*D/SIGDE*
     &                CHEVAL(NCHC(4),WCHC(4),CHCO(1,4),-ARGSQ)
               PDER2 = TWO*(D**2 - EFO*PDER)
            ELSE
               PDER = D/FCALC
               PDER2 = ZERO
            END IF
         END IF
      ELSE
C
C if sigma_delta is equal to zero
C for both centric and acentric
         EFO = D*FCALC
         DELTA = ZERO
         PDER2 = ZERO
         IF (FCALC.GT.ZERO) THEN
            PDER = D/FCALC
         ELSE
            PDER = ZERO
         END IF
      END IF
C
C calculate variance
      VARML = -DELTA*(D*FCALC+EFO)+SIGDE**2+SIGMA**2
C
C value returned is the derivative multiplied by complex Fcalc
      VSTACK(I,VLEVEL)=(ONE/VARML)*((EFO-FOBS)*PDER+
     &     HALF*(ONE-(EFO-FOBS)**2/VARML)*PDER2)*VSTACK(I,VLEVEL+2)
C
      END DO
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE XDOMLFF(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate <Fo>
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C  mlff(Fobs, sigmaF, Fcalc, D, sigma_delta) from script level
C
C arguments Fobs, sigmaF, Fcalc, D, sigma_delta
C
C Fobs: real
C sigmaF: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDOMLFF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C =======================================================================
C
      SUBROUTINE XDOMLFF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
c OpenMP by Kay Diederichs
C
C Copyright 1996  University of Alberta and Yale University
C
C For non-centric reflections
C
C <Fo> = sqrt(Pi)*sigd*sqrt(e)/2 * 1F1(-1/2,1,-(DFc/(sigd*sqrt(e))**2)
C
C For centric reflections,
C
C <Fo> = sqrt(2/Pi)*sigd*sqrt(e) *
C          1F1(-1/2,1/2,-(DFc)**2/(2*(sigd*sqrt(e))**2))
C
C N.B. The partial derivatives all include a factor of Fc,
C to avoid divide by zero problems they have been divided
C by Fc symbolically.
C
      IMPLICIT NONE
C
      INCLUDE 'checof.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION FOBS, SIGMA, FCALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION EFO, DELTA
      DOUBLE PRECISION ARGSQ, ARGSQINV, SIGDE
      DOUBLE PRECISION SQPI2, SQ2PI
C functions
      DOUBLE PRECISION CHEVAL
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR
      DOUBLE PRECISION EIGHT, SIXTEEN, TWENTY
      DOUBLE PRECISION TWOFIVE, FOURNINE, TWOSEVEN, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, EIGHT=8.0D0, SIXTEEN=16.0D0)
      PARAMETER (TWENTY=20.0D0, TWOFIVE=25.0D0, FOURNINE=49.0D0)
      PARAMETER (TWOSEVEN=27.0D0, HALF=0.5D0)
C
C begin
      SQPI2=HALF*SQRT(PI)
      SQ2PI=SQRT(TWO/PI)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
!$omp parallel do 
!$omp& private(i,fobs,sigma,fcalc,d,sigmad,epsilon,sigde,argsq,argsqinv,
!$omp&         efo,delta)
!$omp& shared(n,vstack,vlevel,type,qherm,xrnsym,mult,index,nchc,wchc,
!$omp&        chco,sqpi2,sq2pi)
!$omp& default(none)
      DO I=1,N
      FOBS=DBLE(VSTACK(I,VLEVEL))
      SIGMA=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C trap SIGMAD equal zero
      IF (SIGDE.GT.ZERO) THEN
C
C is this an acentric reflection
         IF (TYPE(INDEX(I)).GE.1) THEN
            ARGSQ = (D*FCALC/SIGDE)**2
            IF (ARGSQ.GT.ZERO) THEN
               ARGSQINV = ONE/ARGSQ
            END IF
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use the series approximation.
C The error in this approximation is less than 1 in 10^7 for
C DFc/sigd*sqrt(e) > 3.4
            IF (ARGSQ.LT.(-WCHC(1))) THEN
               EFO = SQPI2*SIGDE*CHEVAL(NCHC(1),WCHC(1),
     &                                  CHCO(1,1),-ARGSQ)
            ELSE
               DELTA = D*FCALC*        (ARGSQINV/FOUR   )*
     &                 (ONE +          (ARGSQINV/EIGHT  )*
     &                 (ONE +    THREE*(ARGSQINV/FOUR   )*
     &                 (ONE +  TWOFIVE*(ARGSQINV/SIXTEEN)*
     &                 (ONE + FOURNINE*(ARGSQINV/TWENTY )*
     &                 (ONE + TWOSEVEN*(ARGSQINV/EIGHT  ))))))
               EFO = D*FCALC+DELTA
            END IF
C
C is this a centric reflection
         ELSE
            ARGSQ = (D*FCALC/SIGDE)**2/TWO
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use DFc. For DFc/sigd*sqrt(e) > 4.72 <Fo> ~ DFc within
C 1 part in 10^7.
            IF (ARGSQ.LT.(-WCHC(3))) THEN
               EFO = SQ2PI*SIGDE*CHEVAL(NCHC(3),WCHC(3),
     &                                  CHCO(1,3),-ARGSQ)
            ELSE
               EFO = D*FCALC
            END IF
         END IF
      ELSE
C
C if sigma_delta is equal to zero
C for both centric and acentric
         EFO = D*FCALC
      END IF
C
C return the value of <Fo>
      VSTACK(I,VLEVEL) = DCMPLX(EFO,ZERO)
C
      END DO
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE XDOMLVF(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the maximum likelihood variance
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C  mlvf(Fobs, sigmaF, Fcalc, D, sigma_delta) from script level
C
C arguments: Fobs, sigmaF, Fcalc, D, sigma_delta
C
C Fobs: real
C sigmaF: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDOMLVF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE XDOMLVF2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
c OpenMP by Kay Diederichs
C
C Copyright 1996  University of Alberta and Yale University
C
C For non-centric reflections
C
C  <Fo> = sqrt(Pi)*sigd*sqrt(e)/2 * 1F1(-1/2,1,-(DFc/(sigd*sqrt(e))**2)
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C For centric reflections,
C
C  <Fo> = sqrt(2/Pi)*sigd*sqrt(e) *
C          1F1(-1/2,1/2,-(DFc)**2/(2*(sigd*sqrt(e))**2))
C
C   variance = (DFc-<Fo>)(DFc+<Fo>) + (sigd*sqrt(e))**2 + sigf**2
C
C The target is the minus log likelihood:
C
C   -LL = ln(sqrt(variance)) + (<Fo>-Fo)^2/(2*variance)
C       = 1/2 * [ln(variance) + (<Fo>-Fo)^2/variance]
C
C N.B. The partial derivatives all include a factor of Fc,
C to avoid divide by zero problems they have been divided
C by Fc symbolically.
C
      IMPLICIT NONE
C
      INCLUDE 'checof.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION FOBS, SIGMA, FCALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION EFO, VARML, DELTA
      DOUBLE PRECISION ARGSQ, ARGSQINV, SIGDE
      DOUBLE PRECISION SQPI2, SQ2PI
C functions
      DOUBLE PRECISION CHEVAL
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR
      DOUBLE PRECISION EIGHT, SIXTEEN, TWENTY
      DOUBLE PRECISION TWOFIVE, FOURNINE, TWOSEVEN, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, EIGHT=8.0D0, SIXTEEN=16.0D0)
      PARAMETER (TWENTY=20.0D0, TWOFIVE=25.0D0, FOURNINE=49.0D0)
      PARAMETER (TWOSEVEN=27.0D0, HALF=0.5D0)
C
C begin
      SQPI2=HALF*SQRT(PI)
      SQ2PI=SQRT(TWO/PI)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
!$omp parallel do 
!$omp& private(i,fobs,sigma,fcalc,d,sigmad,epsilon,sigde,argsq,argsqinv,
!$omp&         efo,delta,varml)
!$omp& shared(n,vstack,vlevel,type,qherm,xrnsym,mult,index,nchc,wchc,
!$omp&        chco,sqpi2,sq2pi)
!$omp& default(none)
      DO I=1,N
      FOBS=DBLE(VSTACK(I,VLEVEL))
      SIGMA=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C trap SIGMAD equal zero
      IF (SIGDE.GT.ZERO) THEN
C
C is this an acentric reflection
         IF (TYPE(INDEX(I)).GE.1) THEN
            ARGSQ = (D*FCALC/SIGDE)**2
            IF (ARGSQ.GT.ZERO) THEN
               ARGSQINV = ONE/ARGSQ
            END IF
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use the series approximation.
C The error in this approximation is less than 1 in 10^7 for
C DFc/sigd*sqrt(e) > 3.4
            IF (ARGSQ.LT.(-WCHC(1))) THEN
               EFO = SQPI2*SIGDE*CHEVAL(NCHC(1),WCHC(1),
     &                                  CHCO(1,1),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               DELTA = D*FCALC*        (ARGSQINV/FOUR   )*
     &                 (ONE +          (ARGSQINV/EIGHT  )*
     &                 (ONE +    THREE*(ARGSQINV/FOUR   )*
     &                 (ONE +  TWOFIVE*(ARGSQINV/SIXTEEN)*
     &                 (ONE + FOURNINE*(ARGSQINV/TWENTY )*
     &                 (ONE + TWOSEVEN*(ARGSQINV/EIGHT  ))))))
               EFO = D*FCALC+DELTA
            END IF
C
C is this a centric reflection
         ELSE
            ARGSQ = ((D*FCALC/SIGDE)**2)/TWO
C
C calculate <Fo>
C use the Chebyshev polynomial approximation to 1F1 for low values of
C ARGSQ, else use DFc. For DFc/sigd*sqrt(e) > 4.72 <Fo> ~ DFc within
C 1 part in 10^7.
            IF (ARGSQ.LT.(-WCHC(3))) THEN
               EFO = SQ2PI*SIGDE*CHEVAL(NCHC(3),WCHC(3),
     &                                  CHCO(1,3),-ARGSQ)
               DELTA = EFO-D*FCALC
            ELSE
               EFO = D*FCALC
               DELTA = ZERO
            END IF
         END IF
      ELSE
C
C if sigma_delta is equal to zero
C for both centric and acentric
         EFO = D*FCALC
         DELTA = ZERO
      END IF
C
C calculate variance
      VARML = -DELTA*(D*FCALC+EFO)+SIGDE**2+SIGMA**2
C
      VSTACK(I,VLEVEL) = DCMPLX(VARML,ZERO)
C
      END DO
C
      RETURN
      END
C
C ===================================================================
C
      DOUBLE PRECISION FUNCTION CHEVAL(N,W,C,Z)
C
C Calculate the Chebychev polynomial approximation to 1F1 given
C the argument Z, coefficients C, range W, and number of terms N.
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C
      DOUBLE PRECISION C(*)
      DOUBLE PRECISION W,Z
      INTEGER N
C
C local
      DOUBLE PRECISION X, NUM1, NUM2, NUM3, XFAC
      INTEGER I
C
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C begin
      X = Z/W
      IF (X.GT.ONE.OR.X.LT.ZERO) THEN
         CALL WRNDIE(-5,'CHEVAL','argument out of range')
      END IF
C
      XFAC = TWO*(TWO*X-ONE)
C
      NUM1=C(N)
      NUM2=XFAC*NUM1+C(N-1)
      DO I=(N-2),1,-1
         NUM3=XFAC*NUM2-NUM1+C(I)
         NUM1=NUM2
         NUM2=NUM3
      END DO
      CHEVAL=NUM2-XFAC*NUM1/TWO
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE CHEINI
C
C Initialize the Chebyshev coefficients to approximate the
C hypergeometric functions that will be used
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C
      INCLUDE 'checof.inc'
C
C local
      DOUBLE PRECISION QUART, HALF, THREQ, ONE, TWO, THREETWO
      PARAMETER (QUART=0.25D0, HALF=0.5D0, THREQ=0.75D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, THREETWO=1.5D0)
C
C begin
C
C number of terms in approximation
      NCHC(1)=36
      NCHC(2)=37
      NCHC(3)=33
      NCHC(4)=35
      NCHC(5)=36
      NCHC(6)=36
      NCHC(7)=36
      NCHC(8)=36
C
C weight (scale) for approximation
      WCHC(1)=-45.0D0
      WCHC(2)=-45.0D0
      WCHC(3)=-30.0D0
      WCHC(4)=-35.0D0
      WCHC(5)=-35.0D0
      WCHC(6)=-35.0D0
      WCHC(7)=-35.0D0
      WCHC(8)=-35.0D0
C
C coefficients for acentric maximum likelihood
C 1F1(-1/2,1,z) over 0>z>-45
      CALL CHECOF(-HALF,      ONE,WCHC(1),MXCHVL,CHCO(1,1))
C
C coefficients for centric maximum likelihood
C 1F1(1/2,2,z) over 0>z>-45
      CALL CHECOF( HALF,      TWO,WCHC(2),MXCHVL,CHCO(1,2))
C
C coefficients for acentric maximum likelihood derivatives
C 1F1(-1/2,1/2,z) over 0>z>-30
      CALL CHECOF(-HALF,     HALF,WCHC(3),MXCHVL,CHCO(1,3))
C
C coefficients for centric maximum likelihood derivatives
C 1F1(1/2,3/2,z) over 0>z>-35
      CALL CHECOF( HALF, THREETWO,WCHC(4),MXCHVL,CHCO(1,4))
C
C coefficients for acentric maximum likelihood intensity based
C 1F1(1/4,1/2,z) over 0>z>-35
      CALL CHECOF( QUART,    HALF,WCHC(5),MXCHVL,CHCO(1,5))
C
C coefficients for centric maximum likelihood intensity based
C 1F1(-1/4,1/2,z) over 0>z>-35
      CALL CHECOF(-QUART,    HALF,WCHC(6),MXCHVL,CHCO(1,6))
C
C coefficients for acentric maximum likelihood derivatives intensity based
C 1F1(3/4,3/2,z) over 0>z>-35
      CALL CHECOF( THREQ,THREETWO,WCHC(7),MXCHVL,CHCO(1,7))
C
C coefficients for centric maximum likelihood derivatives intensity based
C 1F1(1/4,3/2,z) over 0>z>-35
      CALL CHECOF( QUART,THREETWO,WCHC(8),MXCHVL,CHCO(1,8))
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE CHECOF(AP,CP,W,N,C)
C
C Computes the coefficients in the Chebyshev expansion of 1F1(a;c;z)
C
C Based on Chapter 5 of Y.L.Luke, Algorithms for the Computation of
C Mathematical Functions, Academic Press, New York, 1977
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C
      DOUBLE PRECISION AP, CP, W
      DOUBLE PRECISION C(*)
      INTEGER N
C
C local
      INTEGER I, COUNT
      DOUBLE PRECISION NUM1, NUM2, NUM3, Z, X1, X2, XA, X3, X3A
      DOUBLE PRECISION RHO, P
C
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, TINY
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (FOUR=4.0D0, TINY=1.0D-20)
C
C begin
      NUM1=ZERO
      NUM2=ZERO
      NUM3=TINY
      Z=FOUR/W
      COUNT=N
      C(COUNT)=TINY
      X1=N-ONE
      X2=N
      XA=X1+AP
      X3A=X1+THREE-AP
      DO I=1,N-1
         X1=X1-ONE
         X2=X2-ONE
         XA=XA-ONE
         X3=ONE/(X1+TWO)
         X3A=X3A-ONE
         COUNT=COUNT-1
         C(COUNT)=X2*((Z*(X1+CP)      -X3A*X3)*NUM3+
     &                (Z*(X1+THREE-CP) +XA/X2)*NUM2+
     &                 X3A*X3*NUM1)/XA
         NUM1=NUM2
         NUM2=NUM3
         NUM3=C(COUNT)
      END DO
      C(1)=C(1)/TWO
C
      RHO=C(1)
      P=ONE
      DO I=2,N
         RHO=RHO-P*C(I)
         P=-P
      END DO
C
      DO I=1,N
         C(I)=C(I)/RHO
      END DO
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOMLI(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the maximum likelihood value based on intensity
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C  mli(Iobs, sigmaI, Fcalc, D, sigma_delta) from script level
C
C Iobs: real
C sigmaI: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
C arguments: Iobs, sigmaI, Fcalc, D, sigma_delta
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDOMLI2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C =====================================================================
C
      SUBROUTINE XDOMLI2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C Calculate log likelihood for intensity based maximum likelihood target
C See Pannu and Read,
C     Improved Structure Refinement through Maximum Likelihood,
C     Acta Cryst A52, 659-669 (1996)
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION IOBS, SIGI, FCALC, ICALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION ARG, VAL, SIGDE, VAR, LIKELIHOOD, SUM, SIM
      DOUBLE PRECISION SQTWO
      LOGICAL QAPPROX
C functions
      DOUBLE PRECISION LOGI0
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, HALF=0.5D0)
      DOUBLE PRECISION MAXARG, MINVAL, MAXVAL
      PARAMETER (MAXARG = 250.0D0, MINVAL = 10.0D0, MAXVAL = 250.0D0)
C
C begin
      SQTWO=SQRT(TWO)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
      DO I=1,N
      IOBS=DBLE(VSTACK(I,VLEVEL))
      SIGI=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate Icalc
      ICALC=FCALC*FCALC
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C is this an acentric reflection
      IF (TYPE(INDEX(I)).GE.1) THEN
         VAR = SIGDE*SIGDE
         IF ((VAR.GT.ZERO).AND.(SIGI.GT.ZERO)) THEN
            ARG = SIGI/VAR - IOBS/SIGI
            VAL = (D*D*SIGI)/(VAR*VAR)
         ELSE
            ARG = ZERO
            VAL = ONE
         ENDIF
C
         QAPPROX=.FALSE.
         IF (.NOT.(((ABS(ARG).GT.MAXARG).AND.((VAL*ICALC).GT.MINVAL))
     &       .OR.((VAL*ICALC).GT.MAXVAL).OR.(SIGDE.EQ.ZERO)
     &       .OR.(SIGI.EQ.ZERO))) THEN
            IF (ARG.GT.ZERO) THEN
               VAL = VAL/SQTWO
               CALL APARG(ARG,VAL,ICALC,SUM,SIM)
            ELSE
               CALL ANARG(ARG,VAL,ICALC,SUM,SIM)
            END IF
            IF (SUM.GT.ZERO) THEN
               LIKELIHOOD = D*D*ICALC/VAR - LOG(SUM/VAR) +
     &                      IOBS*IOBS/(TWO*SIGI*SIGI)
               IF (ARG.LE.ZERO) THEN
                  LIKELIHOOD = LIKELIHOOD - (ARG*ARG)/TWO
               END IF
            ELSE
               QAPPROX=.TRUE.
            END IF
         ELSE
            QAPPROX=.TRUE.
         END IF
C
         IF (QAPPROX) THEN
            IF (((ARG.LT.ZERO).AND.(SIGDE.GT.ZERO)).OR.
     &          (SIGI.EQ.ZERO)) THEN
               LIKELIHOOD = LOG(VAR/SQRT(TWO*PI)) +
     &                      (IOBS + D*D*ICALC)/VAR -
     &                      LOGI0(TWO*D*SQRT(IOBS*ICALC)/VAR)
            ELSE
               VAR = VAR*(TWO*D*D*ICALC + VAR) + SIGI*SIGI
               LIKELIHOOD = HALF*( LOG(VAR) +
     &                      (IOBS - D*D*ICALC - SIGDE*SIGDE)**2/VAR)
            END IF
         END IF
C
C is this a centric reflection
      ELSE
         IF ((SIGDE.GT.ZERO).AND.(SIGI.GT.ZERO)) THEN
            ARG = SIGI/(TWO*SIGDE*SIGDE) - IOBS/SIGI
            VAL = (D*D*SIGI)/(TWO*SIGDE**4)
         ELSE
            ARG = ZERO
            VAL = ONE
         ENDIF
C
         QAPPROX=.FALSE.
         IF (.NOT.(((ABS(ARG).GT.MAXARG).AND.((VAL*ICALC).GT.MINVAL))
     &       .OR.((VAL*ICALC).GT.MAXVAL).OR.(SIGDE.EQ.ZERO)
     &       .OR.(SIGI.EQ.ZERO))) THEN
            IF (ARG.GT.ZERO) THEN
               VAL = VAL/SQTWO
               CALL CPARG(ARG,VAL,ICALC,SUM,SIM)
            ELSE
               CALL CNARG(ARG,VAL,ICALC,SUM,SIM)
            END IF
            IF (SUM.GT.ZERO) THEN
               VAR = TWO*SIGDE*SIGDE
               LIKELIHOOD = D*D*ICALC/VAR - LOG(SUM/SQRT(SIGI*VAR))+
     &                      IOBS*IOBS/(TWO*SIGI*SIGI)
               IF (ARG.LE.ZERO) THEN
                  LIKELIHOOD = LIKELIHOOD - (ARG*ARG)/TWO
               END IF
            ELSE
               QAPPROX=.TRUE.
            END IF
         ELSE
            QAPPROX=.TRUE.
         END IF
C
         IF (QAPPROX) THEN
            VAR = SIGDE*SIGDE
            IF (((ARG.LT.ZERO).AND.(SIGDE.GT.ZERO)).OR.
     &          (SIGI.EQ.ZERO)) THEN
               LIKELIHOOD = HALF*LOG(VAR*IOBS) +
     &                      (IOBS + D*D*ICALC)/(TWO*VAR) -
     &                      D*SQRT(IOBS*ICALC)/VAR -
     &                      LOG((ONE + EXP(-TWO*D*
     &                          SQRT(IOBS*ICALC)/VAR))/TWO)
            ELSE
               VAR = TWO*VAR
               VAR = VAR*(TWO*D*D*ICALC + VAR) + SIGI*SIGI
               LIKELIHOOD = HALF*( LOG(VAR) +
     &                     (IOBS - D*D*ICALC - TWO*SIGDE*SIGDE)**2/VAR)
            END IF
         END IF
      END IF
C
C
C value returned is the log likelihood
      VSTACK(I,VLEVEL)= DCMPLX(LIKELIHOOD,ZERO)
C
      END DO
C
      RETURN
      END
C
C ==================================================================
C
      SUBROUTINE XDODMLI(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                  INDEX,HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the maximum likelihood derivatives based on intensity
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C to be used as <sf-obj> =
C   dmli(Iobs, sigmaI, Fcalc, D, sigma_delta) from script level
C
C arguments Iobs, sigmaI, Fcalc, D, sigma_delta
C Iobs: real
C sigmaI: real
C Fcalc: complex
C D: real
C sigma_delta: real
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C
      CALL XDODMLI2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &             INDEX,HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM)
C
      RETURN
      END
C
C ===================================================================
C
      SUBROUTINE XDODMLI2(VLEVEL,VMAX,VSTACK,N,ARGS,
     &                   INDEX,TYPE,MULT,XRNSYM,QHERM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C Calculate derivatives for intensity based maximum likelihood target
C See Pannu and Read,
C     Improved Structure Refinement through Maximum Likelihood,
C     Acta Cryst A52, 659-669 (1996)
C
      IMPLICIT NONE
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      LOGICAL QHERM
C local
      INTEGER I
      DOUBLE PRECISION IOBS, SIGI, FCALC, ICALC, D, SIGMAD, EPSILON
      DOUBLE PRECISION ARG, VAL, SIGDE, VAR, DERIC, SUM, SIM
      DOUBLE PRECISION SQTWO, DIFF
      LOGICAL QAPPROX
C function
      DOUBLE PRECISION SIMILAR
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, HALF=0.5D0)
      DOUBLE PRECISION MAXARG, MINVAL, MAXVAL
      PARAMETER (MAXARG = 250.0D0, MINVAL = 10.0D0, MAXVAL = 250.0D0)
C
C begin
      SQTWO=SQRT(TWO)
C
C loop over all selected reflections
      VLEVEL=VLEVEL-4
      DO I=1,N
      IOBS=DBLE(VSTACK(I,VLEVEL))
      SIGI=DBLE(VSTACK(I,VLEVEL+1))
      FCALC=ABS(VSTACK(I,VLEVEL+2))
      D=DBLE(VSTACK(I,VLEVEL+3))
      SIGMAD=DBLE(VSTACK(I,VLEVEL+4))
C
C calculate Icalc
      ICALC=FCALC*FCALC
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C acentric
         EPSILON=2*XRNSYM/MULT(INDEX(I))
      ELSE
C centric
         EPSILON=XRNSYM/MULT(INDEX(I))
      END IF
C
C combine sigmad and epsilon
      SIGDE=SIGMAD*SQRT(EPSILON)
C
C is this an acentric reflection
      IF (TYPE(INDEX(I)).GE.1) THEN
         VAR = SIGDE*SIGDE
         IF ((VAR.GT.ZERO).AND.(SIGI.GT.ZERO)) THEN
            ARG = SIGI/VAR - IOBS/SIGI
            VAL = (D*D*SIGI)/(VAR*VAR)
         ELSE
            ARG = ZERO
            VAL = ONE
         ENDIF
C
         QAPPROX=.FALSE.
         IF (.NOT.(((ABS(ARG).GT.MAXARG).AND.((VAL*ICALC).GT.MINVAL))
     &       .OR.((VAL*ICALC).GT.MAXVAL).OR.(SIGDE.EQ.ZERO)
     &       .OR.(SIGI.EQ.ZERO))) THEN
            IF (ARG.GT.ZERO) THEN
               VAL = VAL/SQTWO
               CALL APARG(ARG,VAL,ICALC,SUM,SIM)
            ELSE
               CALL ANARG(ARG,VAL,ICALC,SUM,SIM)
            END IF
            IF (SUM.GT.ZERO) THEN
               DERIC = D*D/VAR - SIM/SUM
            ELSE
               QAPPROX=.TRUE.
            END IF
         ELSE
            QAPPROX=.TRUE.
         END IF
C
         IF (QAPPROX) THEN
            IF (((ARG.LT.ZERO).AND.(SIGDE.GT.ZERO)).OR.
     &          (SIGI.EQ.ZERO)) THEN
               IF (ICALC.NE.ZERO) THEN
                  DERIC = D/VAR*(D - SQRT(IOBS/ICALC)*
     &                    SIMILAR(TWO*D*SQRT(IOBS*ICALC)/VAR))
               ELSE
                  DERIC = D*D/VAR
               END IF
            ELSE
               DIFF = IOBS - D*D*ICALC - VAR
               VAR = VAR*(TWO*D*D*ICALC + VAR) + SIGI*SIGI
               DERIC = D*D*(SIGDE*SIGDE - DIFF*(ONE +
     &                 DIFF*SIGDE*SIGDE/VAR))/VAR
            END IF
         END IF
C
C is this a centric reflection
      ELSE
         IF ((SIGDE.GT.ZERO).AND.(SIGI.GT.ZERO)) THEN
            ARG = SIGI/(TWO*SIGDE*SIGDE) - IOBS/SIGI
            VAL = (D*D*SIGI)/(TWO*SIGDE**4)
         ELSE
            ARG = ZERO
            VAL = ONE
         ENDIF
C
         QAPPROX=.FALSE.
         IF (.NOT.(((ABS(ARG).GT.MAXARG).AND.((VAL*ICALC).GT.MINVAL))
     &       .OR.((VAL*ICALC).GT.MAXVAL).OR.(SIGDE.EQ.ZERO)
     &       .OR.(SIGI.EQ.ZERO))) THEN
            IF (ARG.GT.ZERO) THEN
               VAL = VAL/SQTWO
               CALL CPARG(ARG,VAL,ICALC,SUM,SIM)
            ELSE
               CALL CNARG(ARG,VAL,ICALC,SUM,SIM)
            END IF
            IF (SUM.GT.ZERO) THEN
               VAR = TWO*SIGDE*SIGDE
               DERIC = D*D/VAR - SIM/SUM
            ELSE
               QAPPROX=.TRUE.
            END IF
         ELSE
            QAPPROX=.TRUE.
         END IF
C
         IF (QAPPROX) THEN
            VAR = SIGDE*SIGDE
            IF (((ARG.LT.ZERO).AND.(SIGDE.GT.ZERO)).OR.
     &          (SIGI.EQ.ZERO)) THEN
               IF (ICALC.NE.ZERO) THEN
                  DERIC = D/(TWO*VAR)*(D - SQRT(IOBS/ICALC)*
     &                    TANH(D*SQRT(IOBS*ICALC)/VAR))
               ELSE
                  DERIC = D*D/(TWO*VAR)
               END IF
            ELSE
               VAR = TWO*VAR
               DIFF = IOBS - D*D*ICALC - VAR
               VAR = VAR*(TWO*D*D*ICALC + VAR) + SIGI*SIGI
               DERIC = D*D*(SIGDE*SIGDE - DIFF*(ONE +
     &                 DIFF*SIGDE*SIGDE/VAR))/VAR
            END IF
         END IF
      END IF
C
C value returned is the derivative multiplied by complex Fcalc
      VSTACK(I,VLEVEL)= TWO*DERIC*VSTACK(I,VLEVEL+2)
C
      END DO
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE APARG(ARG,VAL,ICALC,SUM,SIM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C calculate the theoretically infinite sum needed in computing
C the negative log likelihood and its derivative wrt Icalc
C for acentric reflections when argument is positive
C
      IMPLICIT NONE
C
      INCLUDE 'xmaxl.inc'
C IO
      DOUBLE PRECISION ARG, VAL, ICALC, SUM, SIM
C local
      INTEGER I
      DOUBLE PRECISION DCYLN(0:MXCYLN)
      DOUBLE PRECISION TER, SIMTOT, SUBTOT, BOTH, OLDSUB
      DOUBLE PRECISION ROOTTWO, ARGSQ
      LOGICAL DONE
C parameter
      INTEGER IMAX
      PARAMETER (IMAX=250)
      DOUBLE PRECISION EPS, ZERO, HALF, ONE, TWO, FOUR
      PARAMETER (EPS=1.0D-25, ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, FOUR=4.0D0)
C functions
      DOUBLE PRECISION CHU, CHUS
C begin
      SIM=ZERO
      TER=ONE
      SUBTOT=ZERO
      ROOTTWO=SQRT(TWO)
C
      IF (ARG.LE.FOUR) THEN
         DCYLN(MXCYLN)=CHU(DBLE(MXCYLN)+HALF,HALF,ARG)/ROOTTWO
         DCYLN(MXCYLN-1)=CHU(DBLE(MXCYLN)-HALF,HALF,ARG)/ROOTTWO
         DO I=(MXCYLN-2),0,-1
            DCYLN(I)=DCYLN(I+1)*ARG/ROOTTWO+(HALF*I+ONE)*DCYLN(I+2)
         END DO
         SUM=DCYLN(0)
         IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
            SIM=DCYLN(1)*VAL
            RETURN
         END IF
         I=1
         DONE=.FALSE.
         DO WHILE (I.LE.MXCYLN.AND.(.NOT.DONE))
            OLDSUB=SUBTOT
            TER=TER*VAL
            SIMTOT=TER*DCYLN(I)
            TER=TER*ICALC/I
            SUBTOT=TER*DCYLN(I)
            IF (SUBTOT.GT.BIG) THEN
               SUM=-ONE
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &          I.GT.10) DONE=.TRUE.
            IF (.NOT.DONE) THEN
               SUM=SUM+SUBTOT
               SIM=SIM+SIMTOT
            END IF
            I=I+1
         END DO
         IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
      ELSE
         ARGSQ=ARG*ARG*HALF
         SUM=CHUS(HALF,HALF,ARGSQ)
         IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
            SUM=SUM/ROOTTWO
            SIM=VAL*CHUS(ONE,HALF,ARGSQ)/ROOTTWO
            RETURN
         END IF
         I=1
         DONE=.FALSE.
         DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
            OLDSUB=SUBTOT
            TER=TER*VAL
            BOTH=CHUS((DBLE(I)*HALF+HALF),HALF,ARGSQ)
            SIMTOT=TER*BOTH
            TER=TER*ICALC/I
            SUBTOT=TER*BOTH
            IF (SUBTOT.GT.BIG.OR.SUBTOT.LT.ZERO) THEN
               SUM=-ONE
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &          I.GT.10) THEN
               SUM=SUM/ROOTTWO
               SIM=SIM/ROOTTWO
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE) THEN
               SUM=SUM+SUBTOT
               SIM=SIM+SIMTOT
            END IF
            I=I+1
         END DO
         IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
         IF (.NOT.DONE) THEN
            SUM=SUM/ROOTTWO
            SIM=SIM/ROOTTWO
         END IF
      END IF
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE ANARG(ARG,VAL,ICALC,SUM,SIM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C calculate the theoretically infinite sum needed in computing
C the negative log likelihood and its derivative wrt Icalc
C for acentric reflections when argument is negative
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
      INCLUDE 'xmaxl.inc'
C IO
      DOUBLE PRECISION ARG, VAL, ICALC, SUM, SIM
C local
      INTEGER I
      DOUBLE PRECISION TER, SIMTOT, SUBTOT, OLDSUB
      DOUBLE PRECISION SPIF1, SPIF2, SPIF3, SQPITWO
      LOGICAL DONE
C parameter
      INTEGER MXTERM
      PARAMETER (MXTERM=625)
      DOUBLE PRECISION EPS, ZERO, HALF, ONE, TWO, FOUR
      PARAMETER (EPS=1.0D-25, ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0,FOUR=4.0D0)
C functions
      DOUBLE PRECISION DERFC
      EXTERNAL DERFC
C begin
      SQPITWO=SQRT(PI/TWO)
      SPIF2=SQPITWO*DERFC((ARG/SQRT(TWO)))
      SPIF3=EXP(-ARG*ARG*HALF)-ARG*SPIF2
      SUM=SPIF2
      SIM=ZERO
C
      IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
         SIM=SPIF3*VAL
         RETURN
      END IF
C
      TER=ONE
      SUBTOT=ZERO
C
      I=1
      DONE=.FALSE.
      DO WHILE (I.LE.MXTERM.AND.(.NOT.DONE))
         OLDSUB=SUBTOT
         TER=TER*VAL
         SIMTOT=TER*SPIF3
         TER=TER*ICALC/I
         SUBTOT=TER*SPIF3
         IF (SUBTOT.GT.BIG) THEN
            SUM=-ONE
            DONE=.TRUE.
         END IF
         IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &       I.GT.10) DONE=.TRUE.
         IF (.NOT.DONE) THEN
            SUM=SUM+SUBTOT
            SIM=SIM+SIMTOT
            SPIF1=SPIF2
            SPIF2=SPIF3
            SPIF3=(SPIF1-ARG*SPIF2)/(I+ONE)
         END IF
         I=I+1
      END DO
      IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE CPARG(ARG,VAL,ICALC,SUM,SIM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C calculate the theoretically infinite sum needed in computing
C the negative log likelihood and its derivative wrt Icalc
C for centric reflections when argument is positive
C
      IMPLICIT NONE
C
      INCLUDE 'xmaxl.inc'
C IO
      DOUBLE PRECISION ARG, VAL, ICALC, SUM, SIM
C local
      INTEGER I
      DOUBLE PRECISION DCYLN(0:MXCYLN)
      DOUBLE PRECISION TER, SIMTOT, SUBTOT, BOTH, OLDSUB
      DOUBLE PRECISION ROOTTWO, RRTWO, ARGSQ
      LOGICAL DONE
C parameter
      INTEGER IMAX
      PARAMETER (IMAX=250)
      DOUBLE PRECISION EPS, ZERO, HALF, ONE, TWO, FOUR
      PARAMETER (EPS=1.0D-25, ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0,FOUR=4.0D0)
C functions
      DOUBLE PRECISION CHU, CHUS
C begin
      ROOTTWO=SQRT(TWO)
      RRTWO=TWO**0.25D0
      SIM=ZERO
      TER=ONE
      SUBTOT=ZERO
C
      IF (ARG.LE.FOUR) THEN
         DCYLN(MXCYLN)=CHU(DBLE(MXCYLN),HALF,ARG)/RRTWO
         DCYLN(MXCYLN-1)=CHU(DBLE(MXCYLN-1),HALF,ARG)/RRTWO
         DO I=(MXCYLN-2),0,-1
            DCYLN(I)=(HALF*I+HALF+HALF/TWO)*DCYLN(I+2)+
     &               ARG*DCYLN(I+1)/ROOTTWO
         END DO
         SUM=DCYLN(0)
         IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
            SIM=DCYLN(1)*VAL*HALF
            RETURN
         END IF
         I=1
         DONE=.FALSE.
         DO WHILE (I.LE.MXCYLN.AND.(.NOT.DONE))
            OLDSUB=SUBTOT
            TER=TER*VAL/(TWO*I)
            SIMTOT=TER*DCYLN(I)*I
            TER=TER*ICALC
            SUBTOT=TER*DCYLN(I)
            IF (SUBTOT.GT.BIG) THEN
               SUM=-ONE
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &          I.GT.10) DONE=.TRUE.
            IF (.NOT.DONE) THEN
               SUM=SUM+SUBTOT
               SIM=SIM+SIMTOT
            END IF
            I=I+1
         END DO
         IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
      ELSE
         ARGSQ=ARG*ARG*HALF
         SUM=CHUS((HALF/TWO),HALF,ARGSQ)
         IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
            SUM=SUM/RRTWO
            SIM=HALF*VAL*CHUS((HALF+HALF/TWO),HALF,ARGSQ)/RRTWO
            RETURN
         END IF
         I=1
         DONE=.FALSE.
         DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
            OLDSUB=SUBTOT
            TER=TER*VAL/(TWO*I)
            BOTH=CHUS((DBLE(I)*HALF+HALF/TWO),HALF,ARGSQ)
            SIMTOT=TER*BOTH*I
            TER=TER*ICALC
            SUBTOT=TER*BOTH
            IF (SUBTOT.GT.BIG.OR.SUBTOT.LT.ZERO) THEN
               SUM=-ONE
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &          I.GT.10) THEN
               SUM=SUM/RRTWO
               SIM=SIM/RRTWO
               DONE=.TRUE.
            END IF
            IF (.NOT.DONE) THEN
               SUM=SUM+SUBTOT
               SIM=SIM+SIMTOT
            END IF
            I=I+1
         END DO
         IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
         IF (.NOT.DONE) THEN
            SUM=SUM/RRTWO
            SIM=SIM/RRTWO
         END IF
      END IF
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE CNARG(ARG,VAL,ICALC,SUM,SIM)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C calculate the theoretically infinite sum needed in computing
C the negative log likelihood and its derivative wrt Icalc
C for centric reflections when argument is negative
C
      IMPLICIT NONE
C
      INCLUDE 'consta.inc'
      INCLUDE 'checof.inc'
      INCLUDE 'xmaxl.inc'
C IO
      DOUBLE PRECISION ARG, VAL, ICALC, SUM, SIM
C local
      INTEGER I
      DOUBLE PRECISION ARGSQ, SPIF1, SPIF2, SPIF3
      DOUBLE PRECISION ROOTPI, ROOTTWO, RRTWO, TWOTQ
      DOUBLE PRECISION TER, SUBTOT, OLDSUB, SIMTOT
      LOGICAL DONE
C parameter
      INTEGER MXTERM
      PARAMETER (MXTERM=521)
      DOUBLE PRECISION EPS, ZERO, QUART, HALF, THREQ, ONE
      DOUBLE PRECISION ONEHALF, TWO, FOUR
      PARAMETER (EPS=1.0D-25, ZERO=0.0D0, QUART=0.25D0, HALF=0.5D0)
      PARAMETER (THREQ=0.75D0, ONE=1.0D0, ONEHALF=1.5D0)
      PARAMETER (TWO=2.0D0, FOUR=4.0D0)
      DOUBLE PRECISION GAMMA1, GAMMA2, GAMMA3
      PARAMETER (GAMMA1=1.225416702465177D0)
      PARAMETER (GAMMA2=3.625609908221908D0)
      PARAMETER (GAMMA3=0.906402477055477D0)
C functions
      DOUBLE PRECISION MFUNC, CHEVAL
C begin
      ROOTPI=SQRT(PI)
      ROOTTWO=SQRT(TWO)
      RRTWO=TWO**0.25D0
      TWOTQ=TWO**0.75D0
      ARGSQ=-ARG*ARG*HALF
      IF (ARGSQ.LE.WCHC(5)) THEN
         SPIF2=ROOTPI/RRTWO*
     &        (MFUNC(QUART,HALF,ARGSQ)/GAMMA1 -
     &         ROOTTWO*ARG*
     &         MFUNC(THREQ,ONEHALF,ARGSQ)/GAMMA2)
         SPIF3=ROOTPI/TWOTQ*
     &        (MFUNC(-QUART,HALF,ARGSQ)/GAMMA3 -
     &         ROOTTWO*ARG*
     &         MFUNC(QUART,ONEHALF,ARGSQ)/GAMMA1)
      ELSE
         SPIF2=ROOTPI/RRTWO*
     &        (CHEVAL(NCHC(5),WCHC(5),CHCO(1,5),ARGSQ)/GAMMA1 -
     &         ROOTTWO*ARG*
     &         CHEVAL(NCHC(7),WCHC(7),CHCO(1,7),ARGSQ)/GAMMA2)
         SPIF3=ROOTPI/TWOTQ*
     &        (CHEVAL(NCHC(6),WCHC(6),CHCO(1,6),ARGSQ)/GAMMA3 -
     &         ROOTTWO*ARG*
     &         CHEVAL(NCHC(8),WCHC(8),CHCO(1,8),ARGSQ)/GAMMA1)
      END IF
C
      SUM=SPIF2
      SIM=ZERO
      TER=ONE
C
      IF (VAL.EQ.ZERO.OR.ICALC.EQ.ZERO) THEN
         SIM=HALF*SPIF3*VAL
         RETURN
      END IF
C
      SUBTOT=ZERO
C
      I=1
      DONE=.FALSE.
      DO WHILE (I.LE.MXTERM.AND.(.NOT.DONE))
         OLDSUB=SUBTOT
         TER=TER*VAL/(TWO*I)
         SIMTOT=TER*SPIF3*I
         TER=TER*ICALC
         SUBTOT=TER*SPIF3
         IF (SUBTOT.GT.BIG) THEN
            SUM=-ONE
            DONE=.TRUE.
         END IF
         IF (.NOT.DONE.AND.SUBTOT.LT.EPS.AND.OLDSUB.GE.SUBTOT.AND.
     &       I.GT.10) DONE=.TRUE.
         IF (.NOT.DONE) THEN
            SUM=SUM+SUBTOT
            SIM=SIM+SIMTOT
            SPIF1=SPIF2
            SPIF2=SPIF3
            SPIF3=(SPIF1-ARG*SPIF2)/(I+HALF)
         END IF
         I=I+1
      END DO
      IF (.NOT.DONE.AND.SUBTOT.GT.OLDSUB) SUM=-ONE
C
      RETURN
      END
C
C ====================================================================
C
      DOUBLE PRECISION FUNCTION MDGAMM(X)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C Calculate the gamma function only for those small positive
C arguments used in the MLI target
C
      IMPLICIT NONE
C IO
      DOUBLE PRECISION X
C parameter
      INTEGER IMAX
      PARAMETER (IMAX=100)
      DOUBLE PRECISION ZERO, QUART, HALF, THREQ, ONE
      PARAMETER (ZERO=0.0D0, QUART=0.25D0, HALF=0.5D0, THREQ=0.75D0)
      PARAMETER (ONE=1.0D0)
C local
      INTEGER I
      DOUBLE PRECISION TOTAL, XTEMP
      LOGICAL DONE
C
      XTEMP=X
      TOTAL=ONE
      DONE=.TRUE.
      IF (X.GT.ONE) THEN
         I=1
         DONE=.FALSE.
         DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
            XTEMP=X-I
            TOTAL=XTEMP*TOTAL
            IF (XTEMP.LE.ONE) DONE=.TRUE.
            I=I+1
         END DO
      END IF
C
      IF (DONE) THEN
         IF (XTEMP.EQ.QUART) THEN
            MDGAMM=3.625609908221908D0*TOTAL
         ELSEIF (XTEMP.EQ.HALF) THEN
            MDGAMM=1.772453850905516D0*TOTAL
         ELSEIF (XTEMP.EQ.THREQ) THEN
            MDGAMM=1.225416702465177D0*TOTAL
         ELSEIF (XTEMP.EQ.ONE) THEN
            MDGAMM=TOTAL
         ELSE
            CALL WRNDIE(-5,'MDGAMM','argument out of range')
         END IF
      ELSE
         CALL WRNDIE(-5,'MDGAMM','argument out of range')
      END IF
C
      RETURN
      END
C
C ====================================================================
C
      DOUBLE PRECISION FUNCTION MFUNC(A,B,X)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C IO
      DOUBLE PRECISION A, B, X
C parameter
      INTEGER IMAX
      PARAMETER (IMAX=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.0D-15)
      DOUBLE PRECISION ZERO, QUART, HALF, THREQ, ONE
      PARAMETER (ZERO=0.0D0, QUART=0.25D0, HALF=0.5D0, THREQ=0.75D0)
      PARAMETER (ONE=1.0D0)
C local
      INTEGER I
      DOUBLE PRECISION TERM, OLDTERM, NEWTERM, Y, SUM
      LOGICAL DONE
C functions
      DOUBLE PRECISION MDGAMM
C
      OLDTERM=2.5D5
      Y=ABS(X)
      SUM=ONE
      TERM=ONE
C
      I=1
      DONE=.FALSE.
      DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
         TERM=((A+I-ONE)*(I+A-B)/(Y*I))*TERM
         NEWTERM=ABS(TERM)
         IF (NEWTERM.GT.OLDTERM) DONE=.TRUE.
         IF (.NOT.DONE) THEN
            OLDTERM=NEWTERM
            SUM=SUM+TERM
         END IF
         IF (ABS(TERM).LT.EPS) DONE=.TRUE.
         I=I+1
      END DO
C
      MFUNC=MDGAMM(B)/MDGAMM(B-A)*Y**(-A)*SUM
C
      RETURN
      END
C
C ====================================================================
C
      DOUBLE PRECISION FUNCTION DERFC(ARG)
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996 Yale University
C
C Evaluate the complementary error function
C
C Modified from code of: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
      IMPLICIT NONE
C
      INCLUDE 'xmaxl.inc'
C IO
      DOUBLE PRECISION ARG
C parameter
      DOUBLE PRECISION THRESH
      PARAMETER (THRESH=0.46875D0)
      DOUBLE PRECISION ZERO, HALF, ONE, TWO, FOUR, SQRPI, SIXTEN
      PARAMETER (ZERO=0.0D0, HALF=0.5D0, ONE=1.0D0, TWO=2.0D0)
      PARAMETER (FOUR=4.0D0, SQRPI=5.6418958354775628695D-1)
      PARAMETER (SIXTEN=16.0D0)
C local
      DOUBLE PRECISION Y, YSQ, X, XDEN, XNUM, DEL, RESULT
C begin
      X = ARG
      Y = ABS(ARG)
      IF (Y.LE.THRESH) THEN
         YSQ=ZERO
         IF (Y.GT.XSMALL) YSQ = Y*Y
         XNUM = 1.85777706184603153D-1*YSQ
         XDEN = YSQ
C
         XNUM = (((XNUM + 3.16112374387056560D00)*YSQ
     &                  + 1.13864154151050156D02)*YSQ
     &                  + 3.77485237685302021D02)*YSQ
C
         XDEN = (((XDEN + 2.36012909523441209D01)*YSQ
     &                  + 2.44024637934444173D02)*YSQ
     &                  + 1.28261652607737228D03)*YSQ
C
         RESULT = ONE -  X * (XNUM + 3.20937758913846947D03)/
     &                       (XDEN + 2.84423683343917062D03)
C
      ELSEIF (Y.LE.FOUR) THEN
         XNUM = 2.15311535474403846D-8*Y
         XDEN = Y
C
         XNUM = (((((((XNUM + 5.64188496988670089D-1)*Y
     &                      + 8.88314979438837594D00)*Y
     &                      + 6.61191906371416295D01)*Y
     &                      + 2.98635138197400131D02)*Y
     &                      + 8.81952221241769090D02)*Y
     &                      + 1.71204761263407058D03)*Y
     &                      + 2.05107837782607147D03)*Y
C
         XDEN = (((((((XDEN + 1.57449261107098347D01)*Y
     &                      + 1.17693950891312499D02)*Y
     &                      + 5.37181101862009858D02)*Y
     &                      + 1.62138957456669019D03)*Y
     &                      + 3.29079923573345963D03)*Y
     &                      + 4.36261909014324716D03)*Y
     &                      + 3.43936767414372164D03)*Y
C
         RESULT = (XNUM + 1.23033935479799725D03)/
     &            (XDEN + 1.23033935480374942D03)
C
         YSQ = INT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ)*EXP(-DEL)*RESULT
C
      ELSEIF (Y.LT.XBIG) THEN
         YSQ = ONE / (Y * Y)
         XNUM = 1.63153871373020978D-2*YSQ
         XDEN = YSQ
C
         XNUM = ((((XNUM + 3.05326634961232344D-1)*YSQ
     &                   + 3.60344899949804439D-1)*YSQ
     &                   + 1.25781726111229246D-1)*YSQ
     &                   + 1.60837851487422766D-2)*YSQ
C
         XDEN = ((((XDEN + 2.56852019228982242D00)*YSQ
     &                   + 1.87295284992346047D00)*YSQ
     &                   + 5.27905102951428412D-1)*YSQ
     &                   + 6.05183413124413191D-2)*YSQ
C
         RESULT = YSQ *(XNUM + 6.58749161529837803D-4)/
     &                 (XDEN + 2.33520497626869185D-3)
C
         RESULT = (SQRPI-RESULT)/Y
         YSQ = INT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ)*EXP(-DEL)*RESULT
      ELSE
         RESULT = ZERO
      END IF
C
      IF (Y.GT.THRESH.AND.X.LT.ZERO) RESULT=TWO-RESULT
C
      DERFC = RESULT
C
      RETURN
      END
C
C ====================================================================
C
      DOUBLE PRECISION FUNCTION CHU(A,B,X)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C IO
      DOUBLE PRECISION A, B, X
C local
      INTEGER I
      DOUBLE PRECISION HALFX, GAM(9), P, V, TERM, VAL, SUBTOT
      DOUBLE PRECISION NUMBER, INVNUM, LOGNUM, SUM
C parameter
      DOUBLE PRECISION ZERO, QUART, HALF, ONE, TWO, THREE, FOUR
      DOUBLE PRECISION FIVE, SEVEN, EIGHT, NINE, THRTEN, TWOFIVE
      DOUBLE PRECISION NINEFIVE, FOURTWO, TWOONE, TWOTWO, THREQ
      DOUBLE PRECISION NINETY
      PARAMETER (ZERO=0.0D0, QUART=0.25D0, HALF=0.5D0, ONE=1.0D0)
      PARAMETER (TWO=2.0D0, THREE=3.0D0, FOUR=4.0D0, FIVE=5.0D0)
      PARAMETER (SEVEN=7.0D0, EIGHT=8.0D0, NINE=9.0D0)
      PARAMETER (THRTEN=13.0D0, TWOFIVE=2.5D0, NINEFIVE=9.5D0)
      PARAMETER (FOURTWO=42.0D0, TWOONE=21.0D0, TWOTWO=22.0D0)
      PARAMETER (THREQ=0.75D0, NINETY=90.0D0)
C
      HALFX = HALF*X
C
      GAM(2) = -HALFX*HALFX
      GAM(1) = GAM(2)*HALFX*TWO/THREE
      GAM(4) = GAM(2)*GAM(2)*TWO
      GAM(3) = HALFX*(-ONE+GAM(4)/FIVE)
      GAM(5) = GAM(1)*(GAM(4)*THREE/SEVEN-EIGHT)
      GAM(6) = GAM(2)*(EIGHT*GAM(4)/THREE-NINE)
      GAM(7) = HALFX*(NINEFIVE-GAM(4)*(THRTEN-GAM(4)*TWOFIVE/NINE))
      GAM(8) = GAM(4)*(GAM(4)*FOUR-FOURTWO)
      GAM(9) = GAM(1)*(287.5D0+GAM(4)*(TWOONE*GAM(4)/TWOTWO-NINETY))
C
      P=SQRT(A)
      V=ZERO
      TERM=ONE
      VAL=HALF/P
      DO I=1,9
         TERM=TERM*VAL
         SUBTOT=TERM*GAM(I)
         V=V+SUBTOT
      END DO
C
      NUMBER=THREQ+HALF*A
      INVNUM=ONE/NUMBER
      LOGNUM=LOG(NUMBER)
      SUM=ONE +
     &    INVNUM*(0.083333333333333333D0 +
     &    INVNUM*(0.003472222222222222D0 -
     &    INVNUM*(0.002681327160493827D0 -
     &    INVNUM*(0.000229472093621399D0 +
     &    INVNUM*(0.000784039221720666D0 +
     &    INVNUM*(0.000069728137583659D0 -
     &    INVNUM*(0.000592166437353694D0 -
     &    INVNUM*(0.000051717909082606D0))))))))
C
      CHU=EXP(-P*X+V+X*X*QUART-NUMBER*(LOGNUM-ONE)+
     &    HALF*LOGNUM)/(SQRT(TWO)*SUM)
C
      RETURN
      END
C
C ====================================================================
C
      DOUBLE PRECISION FUNCTION CHUS(A,B,X)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C IO
      DOUBLE PRECISION A, B, X
C local
      INTEGER N, I, KMAX
      DOUBLE PRECISION A1, C, P, Q, R, U(0:50)
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      DO I=0,50
         U(I)=ZERO
      END DO
C
      KMAX=ZERO
      IF (A.LT.ZERO) THEN
         N=INT(A)
      ELSE
         N=INT(ZERO)
      END IF
C
      Q=A-N
      A1=Q
      U(0)=ONE
C
      IF (N.LT.0.AND.A.EQ.B) THEN
         IF (A.GT.ZERO) CALL UABX(A1,A1,X,KMAX,U)
         P=U(0)
         R=P-Q
         DO I=1,-N,1
            R=X*R
            Q=(A1-I)*P
            P=R-Q
         END DO
      ELSE
         IF (A1.GT.ZERO) CALL UABX(A1,B,X,KMAX,U)
         C=ONE+A1-B+X
         A1=A1-ONE
         P=U(0)
         DO I=1,N,1
            R=(C-I)*P-X*Q
            Q=(A1-I)*(Q-P)
            P=R
         END DO
      END IF
C
      IF (P.GT.ZERO) THEN
         CHUS=P
      ELSE
         CHUS=ZERO
      END IF
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE UABX(A,B,X,KMAX,U)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C IO
      INTEGER KMAX
      DOUBLE PRECISION A, B, X, U(0:50)
C local
      INTEGER R, N, I, K
      DOUBLE PRECISION AR, BR, CR, C, ER, MR, P0, P1, P2, V
      DOUBLE PRECISION M0, M1, U1, U2, U3, W, UPRIME
      LOGICAL LARGEX, DONE
C parameter
      INTEGER IMAX
      PARAMETER (IMAX=501)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.0D-15)
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
C functions
      DOUBLE PRECISION MDGAMM
C begin
      N = INT(A)
      IF (A.EQ.DBLE(N)) N=N-1
      A = A-N
      KMAX = KMAX+N
      LARGEX = ((X.GT.6.5D0).AND.(A.NE.B))
C
      IF (LARGEX) THEN
         MR=ONE
      ELSE
         MR=ZERO
         IF(A.EQ.B) THEN
            M0 = A
            M1 = ONE
         ELSE
            M0 = ZERO
            M1 = ONE
            V  = ONE
            R=1
            DONE=.FALSE.
            DO WHILE (R.LE.IMAX.AND.(.NOT.DONE))
               V = V*V/R
               M0 = M0+V
               V = V*(A+R)/(B+R)
               M1 = M1 + V
               IF (V.LT.(M1*EPS)) DONE=.TRUE.
               R=R+1
            END DO
            V  = EXP(-X)*MDGAMM(A+ONE)/MDGAMM(B+ONE)
            M0 = V*(B+A*M0)
            M1 = V*M1
         ENDIF
      ENDIF
      C = A - B
      CR = TWO + C
      BR = X + A + CR
      P0 = ZERO
      V  = ONE
      P1 = ONE
      ER = ONE
      R  = 0
      AR = A + R
      I=0
      DONE=.FALSE.
      DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
        CALL RECURSION(AR,BR,CR,C,ER,MR,P0,P1,P2,V,R,LARGEX)
        IF (R.GT.KMAX) DONE=.TRUE.
        I=I+1
      END DO
      W = P0*P1/ER
      AR = A + R
      I=0
      DONE=.FALSE.
      DO WHILE (I.LE.IMAX.AND.(.NOT.DONE))
         IF (V*(W/P0+MR*(TWO+A/R)).LT.EPS) DONE=.TRUE.
         IF (.NOT.DONE) THEN
            CALL RECURSION(AR,BR,CR,C,ER,MR,P0,P1,P2,V,R,LARGEX)
            AR = AR + ONE
         END IF
         I=I+1
      END DO
      C = C + ONE
      V = X + C
      U2 = ONE
      W = ZERO
      U3 = -TWO*R/(X+SQRT(X*(X+FOUR*R)))
      I = R - 1
      DO R=I,1,-1
         IF (LARGEX) THEN
            W = W+MR*U2
            MR = MR*(R+ONE)/(C+R)
         ENDIF
         U1 = (-X*U3+(V+R)*U2)/(A+R)
         U3 = U3 - U2
         U2 = U1
         IF ((R.GE.N).AND.(R.LE.KMAX)) U(R-N) = U2
         IF (R.EQ.KMAX) UPRIME = U3
      END DO
      U1 = -X*U3+V*U2
      U3 = U3 - U2
      V  = A
      K  = N - 1
      IF (KMAX.EQ.0) UPRIME = U3
      KMAX = KMAX - N
      IF (LARGEX) THEN
         W = X**(-A)/(A*(W+C*U2)+U1)
      ELSE
         W = X**(-B)/(U1*M1-U3*M0)
      ENDIF
      DO R=0,K
         V = V/(A+R)
      END DO
      IF (N.EQ.0) THEN
         K = 1
         U(0)=W*U1
      ELSE
         K = 0
      ENDIF
      W = V*W
      UPRIME = W*UPRIME
      DO R=K,KMAX
         U(R)=W*U(R)
      END DO
C
      RETURN
      END
C
C ====================================================================
C
      SUBROUTINE RECURSION(AR,BR,CR,C,ER,MR,P0,P1,P2,V,R,LARGEX)
C
C N.S. Pannu and R.J. Read
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C IO
      INTEGER R
      DOUBLE PRECISION AR, BR, CR, C, ER, MR, P0, P1, P2, V
      LOGICAL LARGEX
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
      P2=(BR*P1-AR*P0)/CR
      ER=ER*AR/CR
      R=R+1
      IF (LARGEX) MR=MR*(ONE+C/R)
      V=ER/P2
      BR=BR+TWO
      CR=CR+ONE
      P0=P1
      P1=P2
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDOMLHL(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                   HPTYPE,HPMULT,XRNSYM,QHERM)
C
C
C Author: A. T. Brunger
C
C Copyright 1996  Yale University
C
C to be used as <sf-obj> =
C  mlhl(Fobs, Sigma, Pcalc, Fcalc, PA, PB, PC, PD, D, sigma_delta) from script level
C
C modification: inclusion of Sigma and cases without prior phase info, 12/14/08
C
C Fobs: real
C Sigma: real
C Pcalc: real
C Fcalc: complex
C PA: real
C PB: real        {* Hendrickson-Lattman  *}
C PC: real        {* coefficients         *}
C PD: real
C D: real
C sigma_delta: real
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C get integration step size (phistep)
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS(1))
      END IF
C
      NPHI=NINT(S360/PHISTEP)
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOMLHL2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     &              HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &              HEAP(SINP),HEAP(COSP),HEAP(SIN2P),HEAP(COS2P))
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDOMLHL2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     &                    INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                    SINP,COSP,SIN2P,COS2P)
C
C N.S. Pannu, R.J. Read, and A.T. Brunger
c OpenMP by Kay Diederichs
C
C modification: inclusion of Sigma and cases without prior phase info, 12/14/08
C
C Copyright 1997  University of Alberta and Yale University
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL QHERM
C local
      INTEGER I, IPHI
      DOUBLE PRECISION FOBS, SIGF, FCALC, PCALC
      DOUBLE PRECISION D, SIGMAD, EPSILON, SIGDE
      DOUBLE PRECISION PA, PB, PC, PD
      DOUBLE PRECISION ANOT, BNOT, SVAL
      DOUBLE PRECISION LKHOO, VAR, ARG, NABSARG
      DOUBLE PRECISION PHI, MAXVL
      LOGICAL CEN360
C extrinsic functions
      DOUBLE PRECISION LOGI0, BESEI0
C parameter
      DOUBLE PRECISION ZERO, S180, ONE, TWO, R8SMAL
      PARAMETER (ZERO=0.0D0, S180=180.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0, R8SMAL=1.0D-08)
C begin
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI = (IPHI-1)*PHISTEP
      SINP(IPHI) = SIN(PHI*PI/S180)
      COSP(IPHI) = COS(PHI*PI/S180)
      SIN2P(IPHI) = SIN(TWO*PHI*PI/S180)
      COS2P(IPHI) = COS(TWO*PHI*PI/S180)
      END DO
C
      VLEVEL = VLEVEL-9
!$omp parallel do default(none)
!$omp& private(i,fobs,pcalc,fcalc,pa,pb,pc,pd,d,sigmad,cen360,epsilon)
!$omp& private(sigf)
!$omp& private(sigde,arg,anot,bnot,sval,lkhoo,maxvl,iphi,var,nabsarg)
!$omp& shared(n,vstack,vlevel,type,index,xrnsym,mult,nphi,qherm,cosp)
!$omp& shared(sinp,cos2p,sin2p,phistep)
      DO I=1,N
      FOBS = DBLE(VSTACK(I,VLEVEL))
      SIGF = DBLE(VSTACK(I,VLEVEL+1))
      PCALC = DBLE(VSTACK(I,VLEVEL+2))
      FCALC = ABS(VSTACK(I,VLEVEL+3))
      PA = DBLE(VSTACK(I,VLEVEL+4))
      PB = DBLE(VSTACK(I,VLEVEL+5))
      PC = DBLE(VSTACK(I,VLEVEL+6))
      PD = DBLE(VSTACK(I,VLEVEL+7))
      D = DBLE(VSTACK(I,VLEVEL+8))
      SIGMAD = DBLE(VSTACK(I,VLEVEL+9))
C
C trap SIGMAD equal zero
      IF (SIGMAD.GE.RSMALL) THEN
C
C turn PCALC from degrees to radians
      PCALC = PCALC*PI/S180
C
C check for centric reflections with unrestricted phase
      CEN360 = .FALSE.
      IF ((TYPE(INDEX(I)).LT.1)) THEN
         IF (MOD(PCALC,(PI/TWO)).GE.R8SMAL) CEN360 = .TRUE.
      END IF
C
C calculate epsilon
      IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C
C acentric
      EPSILON = TWO*XRNSYM/MULT(INDEX(I))
      ELSE
C
C centric
      EPSILON = XRNSYM/MULT(INDEX(I))
      ENDIF
C
C combine sigmad and epsilon
      SIGDE = SIGMAD*SQRT(EPSILON)
C
      IF ((TYPE(INDEX(I)).GE.1).OR.CEN360) THEN
C
C acentric reflection
      ARG = TWO*FOBS*D*FCALC/(SIGDE*SIGDE+SIGF*SIGF)
      ANOT = ARG*DCOS(PCALC) + PA
      BNOT = ARG*DSIN(PCALC) + PB
C
      IF ((ABS(PC).LT.RSMALL).AND.(ABS(PD).LT.RSMALL)) THEN
C
C calculate likelihood analytically
C check if there is no prior phase information
         IF ((ABS(PA).LT.RSMALL).AND.(ABS(PB).LT.RSMALL)) THEN
            LKHOO = FOBS - D*FCALC
            LKHOO = LKHOO*LKHOO/(SIGDE*SIGDE+SIGF*SIGF) 
     &              - DLOG(BESEI0(ARG))
         ELSE
            SVAL = DSQRT(ANOT*ANOT + BNOT*BNOT)
            LKHOO = (FOBS*FOBS + D*D*FCALC*FCALC)/
     &        (SIGDE*SIGDE+SIGF*SIGF) - LOGI0(SVAL)
         ENDIF
      ELSE
C
C calculate likelihood numerically
      MAXVL = ZERO
      DO IPHI=1,NPHI
      MAXVL = MAX(MAXVL,
     &        ANOT*COSP(IPHI)+
     &        BNOT*SINP(IPHI)+
     &        PC*COS2P(IPHI)+
     &        PD*SIN2P(IPHI))
      END DO
C
      LKHOO = ZERO
      DO IPHI=1,NPHI
      LKHOO = LKHOO + EXP(-MAXVL+
     &        ANOT*COSP(IPHI)+
     &        BNOT*SINP(IPHI)+
     &        PC*COS2P(IPHI)+
     &        PD*SIN2P(IPHI))
      END DO
      LKHOO = LKHOO*PHISTEP
      LKHOO = (FOBS*FOBS + D*D*FCALC*FCALC)/
     &        (SIGDE*SIGDE+SIGF*SIGF) - LOG(LKHOO) - MAXVL
      ENDIF
      ELSE
C
C centric reflection
      VAR = SIGDE*SIGDE+SIGF*SIGF
      ARG = PA*DCOS(PCALC) + PB*DSIN(PCALC) +
     &      FOBS*D*FCALC/VAR
      NABSARG = -DABS(ARG)
      LKHOO = (FOBS*FOBS + D*D*FCALC*FCALC)/
     &        (TWO*VAR) + NABSARG -
     &        DLOG((ONE + DEXP(TWO*NABSARG))/TWO)
      ENDIF
      ELSE
      LKHOO = ZERO
      ENDIF
C
C our target is the log likelihood
      VSTACK(I,VLEVEL) = DCMPLX(LKHOO,ZERO)
C
      END DO
!$omp end parallel do
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDODMLHL(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     &                    HPTYPE,HPMULT,XRNSYM,QHERM)
C
C Routine to calculate the derivative of the maximum likelihood value that
C includes prior phase information.
C
C Author: A. T. Brunger
C
C modification: inclusion of Sigma and cases without prior phase info, 12/14/08
C
C Copyright 1996  Yale University
C
C to be used as <sf-obj> =
C  dmlhl(Fobs, Sigf, Pcalc, Fcalc, PA, PB, PC, PD, D, sigma_delta) from script level
C
C Fobs: real
C Sigf: real
C Pcalc: real
C Fcalc: complex
C PA: real
C PB: real        {* Hendrickson-Lattman  *}
C PC: real        {* coefficients         *}
C PD: real
C D: real
C sigma_delta: real
C
C arguments: Fobs, Pcalc, Fcalc, PA, PB, PC, PD, D, sigma_delta
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), HPTYPE, HPMULT, XRNSYM
      LOGICAL QHERM
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C get integration step size (phistep)
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS(1))
      END IF
C
      NPHI=NINT(S360/PHISTEP)
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDODMLHL2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     &              HEAP(HPTYPE),HEAP(HPMULT),XRNSYM,QHERM,
     &              HEAP(SINP),HEAP(COSP),HEAP(SIN2P),HEAP(COS2P))
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDODMLHL2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     &                     INDEX,TYPE,MULT,XRNSYM,QHERM,
     &                     SINP,COSP,SIN2P,COS2P)
C
C N.S. Pannu, R.J. Read, and A.T. Brunger
c OpenMP by Kay Diederichs
C
C modification: inclusion of Sigma and cases without prior phase info, 12/14/08
C
C Copyright 1997  University of Alberta and Yale University
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), TYPE(*), MULT(*), XRNSYM
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL QHERM
      INTEGER I, IPHI
      DOUBLE PRECISION FOBS, SIGF, FCALC, PCALC
      DOUBLE PRECISION D, SIGMAD, EPSILON, SIGDE
      DOUBLE PRECISION ACALC, BCALC
      DOUBLE PRECISION PA, PB, PC, PD
      DOUBLE PRECISION ANOT, BNOT, SVAL
      DOUBLE PRECISION VAR, ARG
      DOUBLE PRECISION PHI, MAXVL, TEMP, SIM
      DOUBLE PRECISION DERANOT, DERBNOT, DERIV1, DERIV2
      DOUBLE PRECISION LKHOO, DERFC, DERPC
      LOGICAL CEN360
C extrinsic functions
      DOUBLE PRECISION SIMILAR
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO, R8SMAL
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0, R8SMAL=1.0D-08)
C begin
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI = (IPHI-1)*PHISTEP
      SINP(IPHI) = SIN(PHI*PI/S180)
      COSP(IPHI) = COS(PHI*PI/S180)
      SIN2P(IPHI) = SIN(TWO*PHI*PI/S180)
      COS2P(IPHI) = COS(TWO*PHI*PI/S180)
      END DO
C
      VLEVEL = VLEVEL-9
!$omp parallel do default(none)
!$omp& private(i,fobs,pcalc,acalc,bcalc,fcalc,pa,pb,pc,pd,d,sigmad)
!$omp& private(sigf)
!$omp& private(cen360,epsilon,sigde,arg,anot,bnot,sval,derfc,derpc,sim)
!$omp& private(maxvl,iphi,lkhoo,deranot,derbnot,temp,var,deriv1,deriv2)
!$omp& shared(n,vstack,vlevel,type,index,qherm,xrnsym)
!$omp& shared(mult,nphi,cosp,sinp,cos2p,sin2p,phistep)
      DO I=1,N
      FOBS = DBLE(VSTACK(I,VLEVEL))
      SIGF = DBLE(VSTACK(I,VLEVEL+1))
      PCALC = DBLE(VSTACK(I,VLEVEL+2))
      ACALC = DBLE(VSTACK(I,VLEVEL+3))
      BCALC = DIMAG(VSTACK(I,VLEVEL+3))
      FCALC = ABS(VSTACK(I,VLEVEL+3))
      PA = DBLE(VSTACK(I,VLEVEL+4))
      PB = DBLE(VSTACK(I,VLEVEL+5))
      PC = DBLE(VSTACK(I,VLEVEL+6))
      PD = DBLE(VSTACK(I,VLEVEL+7))
      D = DBLE(VSTACK(I,VLEVEL+8))
      SIGMAD = DBLE(VSTACK(I,VLEVEL+9))
C
C trap SIGMAD equal zero
      IF (SIGMAD.GE.RSMALL) THEN
C
C turn PCALC from degrees to radians
         PCALC = PCALC*PI/S180
C
C check for centric reflections with unrestricted phase
         CEN360 = .FALSE.
         IF ((TYPE(INDEX(I)).LT.1)) THEN
            IF (MOD(PCALC,(PI/TWO)).GT.R8SMAL) CEN360 = .TRUE.
         END IF
C
C calculate epsilon
         IF ((TYPE(INDEX(I)).GE.1).AND.(QHERM)) THEN
C
C acentric
            EPSILON = TWO*XRNSYM/MULT(INDEX(I))
         ELSE
C centric
            EPSILON = XRNSYM/MULT(INDEX(I))
         ENDIF
C
C combine sigmad and epsilon
         SIGDE = SIGMAD*SQRT(EPSILON)
C
         IF ((TYPE(INDEX(I)).GE.1).OR.CEN360) THEN
C
C acentric reflection
            ARG = TWO*FOBS*D/(SIGDE*SIGDE+SIGF*SIGF)
            ANOT = ARG*FCALC*DCOS(PCALC) + PA
            BNOT = ARG*FCALC*DSIN(PCALC) + PB
            IF ((ABS(PC).LT.RSMALL).AND.(ABS(PD).LT.RSMALL)) THEN
C
C check if there is no prior phase information
               IF ((ABS(PA).LT.RSMALL).AND.(ABS(PB).LT.RSMALL)) THEN
                  ARG = ARG*FCALC
                  SIM = SIMILAR(ARG)
                  DERFC = -TWO*D*((FOBS-D*FCALC) + 
     &                 (SIM-ONE)*FOBS)/(SIGDE*SIGDE+SIGF*SIGF) 
                  DERPC = ZERO
               ELSE
C calculate the derivatives analytically
                  SVAL = DSQRT(ANOT*ANOT + BNOT*BNOT)
                  IF (SVAL.LT.RSMALL) THEN
                     DERFC = ZERO
                     DERPC = ZERO
                  ELSE
                     SIM = SIMILAR(SVAL)
                     DERFC = SIM*ARG*(ARG*FCALC +
     &                    PA*DCOS(PCALC) +
     &                    PB*DSIN(PCALC))/SVAL
                     DERFC = TWO*D*D*FCALC/(SIGDE*SIGDE+SIGF*SIGF) - 
     &                    DERFC
                     DERPC = SIM*ARG*FCALC*(
     &                    PA*DSIN(PCALC) -
     &                    PB*DCOS(PCALC))/SVAL
                  END IF
               END IF
            ELSE
C calculate the derivatives numerically
               MAXVL = ZERO
               DO IPHI = 1,NPHI
                  MAXVL = MAX(MAXVL, ANOT*COSP(IPHI) +
     &                 BNOT*SINP(IPHI) +
     &                 PC*COS2P(IPHI) +
     &                 PD*SIN2P(IPHI))
               END DO
C
               LKHOO = ZERO
               DO IPHI=1,NPHI
                  LKHOO = LKHOO + EXP(-MAXVL +
     &                 ANOT*COSP(IPHI) +
     &                 BNOT*SINP(IPHI) +
     &                 PC*COS2P(IPHI) +
     &                 PD*SIN2P(IPHI))
               END DO
               LKHOO = LKHOO*PHISTEP
               LKHOO = -LOG(LKHOO) - MAXVL
C
               DERANOT = ZERO
               DERBNOT = ZERO
               DO IPHI=1,NPHI
                  TEMP = EXP(LKHOO +
     &                 ANOT*COSP(IPHI) +
     &                 BNOT*SINP(IPHI) +
     &                 PC*COS2P(IPHI) +
     &                 PD*SIN2P(IPHI))
                  DERANOT = DERANOT + COSP(IPHI)*TEMP
                  DERBNOT = DERBNOT + SINP(IPHI)*TEMP
               END DO
C
               DERANOT = PHISTEP*DERANOT
               DERBNOT = PHISTEP*DERBNOT
               DERFC = ARG*(DERANOT*DCOS(PCALC) +
     &              DERBNOT*DSIN(PCALC))
               DERFC = TWO*D*D*FCALC/(SIGDE*SIGDE+SIGF*SIGF) -
     &              DERFC
               DERPC = ARG*(DERANOT*DSIN(PCALC) -
     &              DERBNOT*DCOS(PCALC))*FCALC
            ENDIF
C
         ELSE
C
C centric reflection
C
            VAR = SIGDE*SIGDE+SIGF*SIGF
            ARG = PA*DCOS(PCALC) + PB*DSIN(PCALC) +
     &           FOBS*D*FCALC/VAR
            DERFC = D*D*FCALC/VAR -
     &           DTANH(ARG)*FOBS*D/VAR
            DERPC = TWO*DTANH(ARG)*(PA*DSIN(PCALC) -
     &           PB*DCOS(PCALC))
         ENDIF
      ELSE
         DERFC = ZERO
         DERPC = ZERO
      ENDIF
C
C our target is the derivative of the log likelihood
      IF (FCALC.LT.RSMALL) THEN
      VSTACK(I,VLEVEL) = DCMPLX(ZERO,ZERO)
      ELSE
      DERIV1 = DERFC*ACALC - DERPC*BCALC/FCALC
      DERIV2 = DERFC*BCALC + DERPC*ACALC/FCALC
      VSTACK(I,VLEVEL) = DCMPLX(DERIV1,DERIV2)/FCALC
      END IF
C
      END DO
!$omp end parallel do
C
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION SIMILAR(Y)
C
C  N.S. Pannu and R.J. Read
C
C  Copyright 1997  University of Alberta
C
C  Function returns BesselI[1,Y]/BesselI[0,Y]
C
      IMPLICIT NONE
C
C constants
      DOUBLE PRECISION EPSILON, LOWERLIM, ZERO, ONE, TWO, FOUR
      INTEGER MAXTERMS
      PARAMETER (EPSILON = 1.0D-15)
      PARAMETER (LOWERLIM = 20.0D0)
      PARAMETER (ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      PARAMETER (MAXTERMS = 150)
C I\O variables
      DOUBLE PRECISION Y
C local variables
      DOUBLE PRECISION X
      DOUBLE PRECISION DPN, TOT0, TOT1, SUBTOT0, SUBTOT1
      INTEGER N
C
C begin
C
C I_1(x)/I_0(x) is an odd function, so just evaluate for positive
C values and fix for negative values at the end.
C
      X = DABS(Y)
      TOT0 = ONE
      SUBTOT0 = ONE
      TOT1 = ONE
      SUBTOT1 = ONE
      IF (X.LT.LOWERLIM) THEN
C
C use Taylor series expansion
      N=1
      DO WHILE ((N.LE.MAXTERMS).AND.(SUBTOT0.GE.EPSILON))
      DPN = DBLE(N)
      SUBTOT0 = X*X*SUBTOT0/(FOUR*DPN*DPN)
      SUBTOT1 = X*X*SUBTOT1/(FOUR*DPN*(DPN+ONE))
      TOT0 = TOT0 + SUBTOT0
      TOT1 = TOT1 + SUBTOT1
      N=N+1
      END DO
      TOT0 = TOT1*X/(TWO*TOT0)
      ELSE
C use asymptotic expansion
      N=1
      DO WHILE ((N.LE.MAXTERMS).AND.(DABS(SUBTOT0).GE.EPSILON))
      DPN = TWO*DBLE(N)
      SUBTOT0 = (DPN - ONE)*(DPN - ONE) /
     &          (FOUR*X*DPN)*SUBTOT0
      TOT0 = TOT0 + SUBTOT0
      TOT1 = (TWO/(ONE - DPN) - ONE) *
     &        SUBTOT0 + TOT1
      N=N+1
      END DO
      TOT0 = TOT1/TOT0
      ENDIF
C
      IF (Y.LT.ZERO) TOT0 = -TOT0
C
      SIMILAR = TOT0
C
      RETURN
      END
C
C======================================================================
C
