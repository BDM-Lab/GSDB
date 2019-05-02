      SUBROUTINE XDOSIGMA(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRNSYM,QHERM,
     &         XRTR,XRMREF,HPH,HPK,HPL,
     &         HPMULT,HPTYPE,MBINS,XBINLOW,XBINHIGH,BINSHELL)
C
C Routine computes SIGMAA coefficient.
C
C Function call:
C     SIGMA(EOBS,ECALC)
C EOBS and ECALC must be normalized structure factor amplitudes.
C
C EOBS: real
C ECALC: real
C SIGMA: real
C
C
C
C Returns: sigmaa.
C
C Note:    This is a special operation, i.e., N has to include
C          all structure factor elements (see routine XDOSPCL).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRMREF, HPH, HPK, HPL
      INTEGER HPMULT, HPTYPE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
C pointer
      INTEGER CI, CJ, CII, CJJ, CIJ, WSUM, ISHELL
      INTEGER SIGMAA, R, DR, SUMM, SIGOLD, SHELL, ASTL2, SUMW
C begin
C
C
C normalize structure factor
      CI=ALLHP(IREAL8(MBINS))
      CII=ALLHP(IREAL8(MBINS))
      CJ=ALLHP(IREAL8(MBINS))
      CJJ=ALLHP(IREAL8(MBINS))
      CIJ=ALLHP(IREAL8(MBINS))
      WSUM=ALLHP(IREAL8(MBINS))
      ISHELL=ALLHP(INTEG4(N))
      SIGMAA=ALLHP(IREAL8(MBINS))
      R=ALLHP(IREAL8(MBINS))
      DR=ALLHP(IREAL8(MBINS))
      SUMM=ALLHP(IREAL8(MBINS))
      SIGOLD=ALLHP(IREAL8(MBINS))
      SHELL=ALLHP(IREAL8(MBINS+1))
      ASTL2=ALLHP(IREAL8(MBINS))
      SUMW=ALLHP(IREAL8(MBINS))
C
      CALL XDOSIGMA2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRNSYM,QHERM,
     &         XRTR,XRMREF,HEAP(HPH),
     &         HEAP(HPK),HEAP(HPL),HEAP(HPMULT),HEAP(HPTYPE),
     &         MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &         HEAP(CI),HEAP(CII),HEAP(CJ),HEAP(CJJ),
     &         HEAP(CIJ),HEAP(WSUM),HEAP(ISHELL),
     &         HEAP(SIGMAA),HEAP(R),HEAP(DR),HEAP(SUMM),
     &         HEAP(SIGOLD),HEAP(SHELL),HEAP(ASTL2),HEAP(SUMW))
C
      CALL FREHP(SUMW,IREAL8(MBINS))
      CALL FREHP(ASTL2,IREAL8(MBINS))
      CALL FREHP(SHELL,IREAL8(MBINS+1))
      CALL FREHP(SIGOLD,IREAL8(MBINS))
      CALL FREHP(SUMM,IREAL8(MBINS))
      CALL FREHP(DR,IREAL8(MBINS))
      CALL FREHP(R,IREAL8(MBINS))
      CALL FREHP(SIGMAA,IREAL8(MBINS))
      CALL FREHP(ISHELL,INTEG4(N))
      CALL FREHP(WSUM,IREAL8(MBINS))
      CALL FREHP(CJJ,IREAL8(MBINS))
      CALL FREHP(CJ,IREAL8(MBINS))
      CALL FREHP(CIJ,IREAL8(MBINS))
      CALL FREHP(CII,IREAL8(MBINS))
      CALL FREHP(CI,IREAL8(MBINS))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOSIGMA2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRNSYM,QHERM,
     &         XRTR,XRMREF,XRH,XRK,XRL,
     &         MULT,TYPE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &         CI,CII,CJ,CJJ,CIJ,WSUM,
     &         ISHELL,SIGMAA,R,DR,SUMM,SIGOLD,SHELL,ASTL2,SUMW)
C
C
C References: R.J. Read, Improved fourier coefficients for maps
C  using phases from partial structures with errors, Acta Cryst. A42,
C 140--149 (1986).
C
C R.J. Read, Structure-factor probabilities for related structures,
C Acta Cryst A46, 900--912 (1990).
C
C
C Author: R.J. Read and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRMREF, XRH(*), XRK(*), XRL(*)
      INTEGER MULT(*), TYPE(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION CI(*), CII(*), CJ(*), CJJ(*), CIJ(*), WSUM(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION SIGMAA(*)
      DOUBLE PRECISION R(*)
      DOUBLE PRECISION DR(*)
      DOUBLE PRECISION SUMM(*)
      DOUBLE PRECISION SIGOLD(*)
      DOUBLE PRECISION SHELL(*)
      DOUBLE PRECISION ASTL2(*)
      DOUBLE PRECISION SUMW(*)
C local
      INTEGER REFLCT, H, K, L, IND, II
      DOUBLE PRECISION SSQ, EO, EC, EO2, EC2, CSUM, DSUM
      DOUBLE PRECISION CIT, CIIT, CJT, CJJT, CIJT, WSUMT, CCTOT, SIGA
      DOUBLE PRECISION FOML, FOMDX, DELMAX, OVERM, DEL
      INTEGER NSM
      INTEGER JLOW, JHIGH
      DOUBLE PRECISION SLOPE, WT, EPSILON, X
      DOUBLE COMPLEX DCPLAC
      DOUBLE PRECISION DPPLAC
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, HALF, TOL, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (HALF=0.5D0, TOL=1.0D-5, FOUR=4.0D0)
C begin
C
      VLEVEL=VLEVEL-1
C
      IF (N.GT.0) THEN
C
      IF (WRNLEV.GE.10) THEN
      IF (XBINLOW.GT.RSMALL) THEN
      WRITE(6,'(A,I8,A,F10.4,A,F10.4,A)')
     & ' XSIGMA: using ',MBINS,' bins between',
     &  ONE/XBINHIGH,' and ',ONE/XBINLOW,' A resolution.'
      ELSE
      WRITE(6,'(A,I8,A,F10.4,A)')
     & ' XSIGMA: using ',MBINS,' bins between',
     &  ONE/XBINHIGH,' and infinity A resolution.'
      END IF
      END IF
C
C
C initialize averages
      CIT=ZERO
      CJT=ZERO
      CIIT=ZERO
      CJJT=ZERO
      CIJT=ZERO
      WSUMT=ZERO
C
      DO IND=1,MBINS
      WSUM(IND)=ZERO
      CI(IND)=ZERO
      CJ(IND)=ZERO
      CII(IND)=ZERO
      CJJ(IND)=ZERO
      CIJ(IND)=ZERO
      ASTL2(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C compute ASTL2
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
      H=XRH(INDEX(REFLCT))
      K=XRK(INDEX(REFLCT))
      L=XRL(INDEX(REFLCT))
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      ASTL2(ISHELL(REFLCT))=ASTL2(ISHELL(REFLCT))+SSQ
      END IF
      END DO
C
C
C compute the structure factor averages for each bin
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
C
      EO2=VSTACK(REFLCT,VLEVEL)**2
      EC2=VSTACK(REFLCT,VLEVEL+1)**2
      WSUM(ISHELL(REFLCT))=WSUM(ISHELL(REFLCT))+ONE
      CI(ISHELL(REFLCT))=CI(ISHELL(REFLCT))+EO2
      CJ(ISHELL(REFLCT))=CJ(ISHELL(REFLCT))+EC2
      CII(ISHELL(REFLCT))=CII(ISHELL(REFLCT))+EO2**2
      CJJ(ISHELL(REFLCT))=CJJ(ISHELL(REFLCT))+EC2**2
      CIJ(ISHELL(REFLCT))=CIJ(ISHELL(REFLCT))+EO2*EC2
C
      WSUMT=WSUMT+ONE
      CIT=CIT+EO2
      CJT=CJT+EC2
      CIIT=CIIT+EO2**2
      CJJT=CJJT+EC2**2
      CIJT=CIJT+EO2*EC2
C
      END IF
      END DO
C
C compute correlation coefficients and average S^2 for each bin
C
C initial sigmaa is the square-root of the bin-wise correlation
C coefficient between Es.
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL) THEN
      DSUM=(CII(IND)-CI(IND)**2/WSUM(IND))*(CJJ(IND)-
     &       CJ(IND)**2/WSUM(IND))
      CSUM=CIJ(IND) - CI(IND)*CJ(IND)/WSUM(IND)
      IF (DSUM.GT.RSMALL) THEN
      DSUM=SQRT(DSUM)
C
      SIGMAA(IND)=CSUM/DSUM
C
      IF (SIGMAA(IND).LT.RSMALL) THEN
      SIGMAA(IND)=ZERO
      ELSE
      SIGMAA(IND)=SQRT(SIGMAA(IND))
      END IF
C
      ELSE
      SIGMAA(IND)=ZERO
      END IF
      ASTL2(IND)=ASTL2(IND)/WSUM(IND)
      ELSE
      SIGMAA(IND)=ZERO
      ASTL2(IND)=ZERO
      END IF
      SIGOLD(IND)=SIGMAA(IND)
      END DO
C
      IF (WSUMT.GT.RSMALL) THEN
      DSUM=(CIIT-CIT**2/WSUMT)*(CJJT-
     &       CJT**2/WSUMT)
      CSUM=CIJT - CIT*CJT/WSUMT
      IF (DSUM.GT.RSMALL) THEN
      DSUM=SQRT(DSUM)
      CCTOT=CSUM/DSUM
      ELSE
      CCTOT=ZERO
      END IF
      ELSE
      CCTOT=ZERO
      END IF
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(2A,F10.5)') ' XSIGMA: overall correlation between',
     &  ' squared arguments: ',CCTOT
      END IF
C
C ITERATIVELY IMPROVE ESTIMATE OF SIGMAA.
C SIGMAA = <|EO|*|EC|*COS(DA)> = <|EO|*|EC|*M> (SRINIVASAN AND
C CHANDRASEKARAN, INDIAN J. PURE APPL. PHYS. 4:178 (1966)), AND
C M IS A FUNCTION OF SIGMAA.  WITH THE CORRECT SIGMAA, THE MEAN
C VALUE OF (SIGMAA-<|EO|*|EC|*M>) SHOULD BE ZERO.  WITH DATA
C NORMALIZED SO THAT WEIGHTED |E|**2 IS EQUAL TO 1 AND
C WEIGHTING THE CENTRIC(1) AND NON-CENTRIC(2) TERMS, THIS
C CORRESPONDS TO THE MAXIMUM LIKELIHOOD ESTIMATE OF LUNIN AND
C URZHUMTSEV (ACTA CRYST. A40:269 (1984)), ALTHOUGH THEY
C CONSIDERED ONLY NON-CENTRIC REFLECTIONS. USE NEWTON'S
C METHOD TO FIND ZERO.  FORM OF DERIVATIVE IS SLIGHTLY DIFFERENT
C FOR CENTRIC AND NON-CENTRIC BECAUSE OF DIFFERENT EXPRESSIONS
C FOR M (SEE S/R SRIN).  ITERATE UNTIL NO SHIFT IS
C GREATER THAN 0.0001, OR FOR A MAXIMUM OF 10 CYCLES.
C
      II=0
      DELMAX=ONE
      DO WHILE (II.LT.10.AND.DELMAX.GE.0.0001D0)
      II=II+1
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A,I3,A)') ' Iteration number',II,' on sigmaa'
      WRITE(6,'(A)')
     & '                  OLD                 NEW       MEAN'
      WRITE(6,'(A)')
     & '    D LIMITS     SIGMAA     SHIFT    SIGMAA     FOM'
      END IF
C
C initialize R, DR, SUMM
      DO IND=1,MBINS
      R(IND)=ZERO
      DR(IND)=ZERO
      SUMM(IND)=ZERO
      SUMW(IND)=ZERO
      END DO
C
C loop through all selected reflections
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
C
      EO=VSTACK(REFLCT,VLEVEL)
      EC=VSTACK(REFLCT,VLEVEL+1)
C current SIGMAA
      SIGA = SIGMAA(ISHELL(REFLCT))
C
      IF (TYPE(INDEX(REFLCT)).GE.1) THEN
C acentric
      EPSILON=2*XRNSYM/MULT(INDEX(REFLCT))
      WT=TWO
      ELSE
C centric
      EPSILON=XRNSYM/MULT(INDEX(REFLCT))
      WT=ONE
      END IF
C
      SUMW(ISHELL(REFLCT))=SUMW(ISHELL(REFLCT))+WT
C
C compute X
      X = TWO*SIGA*EO*EC/(ONE-SIGA**2)
C
C CALCULATE FIGURE OF MERIT, GIVEN X. FORMULA DIFFERS FOR
C CENTRIC AND NON-CENTRIC.  SRINIVASAN (ACTA CRYST 20:143(1966))
C DIFFERS FROM WOOLFSON (ACTA CRYST 9:804(1956)) AND SIM (ACTA
C CRYST 13:511(1960)) ONLY IN THE EXPRESSION FOR X, WHICH IN
C THE PRESENT CASE INCLUDES CONTRIBUTIONS FROM MISSING STRUCTURE
C AND MODEL ERROR.
C FOR NON-CENTRIC REFLECTIONS, USE THE FUNCTION SIM WHICH IS
C THE RATIO OF THE MODIFIED 1ST AND ZERO ORDER BESSEL FUNCTIONS.
C
      IF (TYPE(INDEX(REFLCT)).GE.1) THEN
C acentric reflections
      CALL XSIM(X,FOML)
      ELSE
C centric reflections
      FOML=TANH(X/TWO)
      END IF
C
      SUMM(ISHELL(REFLCT))=SUMM(ISHELL(REFLCT))+FOML
C
C OMIT CONSTANTS IN RESIDUAL AND DERIVATIVE SUMS, BUT
C REMEMBER TO PUT THEM IN AFTER END OF LOOP.
C
      R(ISHELL(REFLCT))=R(ISHELL(REFLCT))-WT*EO*EC*FOML
C
      IF (TYPE(INDEX(REFLCT)).GE.1) THEN
C
C NON-CENTRIC CASE.
C LIMIT AS X TENDS TO ZERO OF I1(X)/(X*I0(X)), IE. FOM/X
C FOR NON-CENTRIC, IS 1/2.  AVOID POTENTIAL DIVIDE BY ZERO
C AND SUBSTITUTE IF APPROPRIATE. INCLUDE WT=2
C
      IF (X.LT.TOL) FOMDX = HALF
      IF (X.GE.TOL) FOMDX = FOML/X
      DR(ISHELL(REFLCT))=DR(ISHELL(REFLCT))
     &          + FOUR*(EO*EC)**2*(FOML**2 + FOMDX - ONE)
      ELSE
C
C CENTRIC CASE. (WT=1)
C
      DR(ISHELL(REFLCT))= DR(ISHELL(REFLCT))
     &          + (EO*EC)**2*(FOML**2 - ONE)
      END IF
      END IF
      END DO
C
C
      DELMAX=ZERO
      OVERM=ZERO
C
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL) THEN
C
C FIRST ADD IN THE CONSTANT TERMS TO R AND DR
C
      R(IND)=R(IND)+SUMW(IND)*SIGOLD(IND)
      DR(IND)=DR(IND)*(ONE+SIGOLD(IND)**2)/((ONE-SIGOLD(IND)**2)**2)
      DR(IND)=DR(IND)+SUMW(IND)
C
C ADJUST CURRENT ESTIMATE OF SIGMAA; CHECK FOR CONVERGENCE
C
      IF (DR(IND).EQ.ZERO) THEN
      IF (SQRT(BINSHELL(IND+1)).GT.RSMALL) THEN
      WRITE(6,'(2A,F10.4,A,F10.4,A)')
     & ' Shift denominator is zero for ',
     &  'bin between ',ONE/SQRT(BINSHELL(IND)),' and ',
     &  ONE/SQRT(BINSHELL(IND+1)),' A.'
      CALL WRNDIE(-5,'XSIGMAA','Shift denominator is zero.')
      ELSE
      WRITE(6,'(2A,F10.4,A,A)')
     & ' Shift denominator is zero for ',
     &  'bin between ',ONE/SQRT(BINSHELL(IND)),' and ',
     &  ' infinity A.'
      END IF
      ELSE
C
      DEL=-R(IND)/DR(IND)
      IF (ABS(DEL).GT.DELMAX) DELMAX=ABS(DEL)
      SIGMAA(IND)=SIGOLD(IND)+DEL
      IF (SIGMAA(IND).GT.ONE-RSMALL) THEN
C
C SIGMAA IS RESTRICTED TO RANGE 0 TO 1
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A,F10.5,A)')
     &  ' Sigmaa=',SIGMAA(IND),' is >= 1.0: Reset to 0.999.'
      END IF
      SIGMAA(IND)=0.999D0
      ELSE
      FOML=SUMM(IND)/WSUM(IND)
      OVERM=OVERM+SUMM(IND)
C
      IF (WRNLEV.GE.10) THEN
      IF (SQRT(BINSHELL(IND+1)).GT.RSMALL) THEN
      WRITE(6,'(F7.2,A,F5.2,4F10.5)')
     & 1/SQRT(BINSHELL(IND)),'--',1/SQRT(BINSHELL(IND+1)),
     & SIGOLD(IND),DEL,SIGMAA(IND),FOML
      ELSE
      WRITE(6,'(F7.2,A,4F10.5)')
     & 1/SQRT(BINSHELL(IND)),'-- infinity',
     & SIGOLD(IND),DEL,SIGMAA(IND),FOML
      END IF
      END IF
C
      SIGOLD(IND) = SIGMAA(IND)
      END IF
      END IF
      END IF
      END DO
C
      END DO
C
C compute overall figure of merit
      OVERM=OVERM/WSUMT
C
      IF (DELMAX.GE.0.0001D0) THEN
      WRITE(6,'(A,F10.5)')
     & '          Overall mean FOM is',OVERM
      CALL WRNDIE(+1,'XSIGMAA',
     &            'Refinement of SIGMAA has not converged.')
      CALL DECLAR( 'STATUS', 'ST', 'INCOMPLETE', DCPLAC, DPPLAC )
      ELSE
      WRITE(6,'(A)')
     & ' XSIGMAA: Refinement of SIGMAA has converged '
      WRITE(6,'(A,F10.5)')
     & '          Overall mean FOM is',OVERM
      CALL DECLAR( 'STATUS', 'ST', 'COMPLETE', DCPLAC, DPPLAC )
      END IF
C
C
      IF (WRNLEV.GE.10) THEN
C
C write SIGMAA as a function of S^2
      WRITE(6,'(A)') ' Ln(SIGMAA) vs. S^2'
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL.AND.SIGMAA(IND).GT.RSMALL) THEN
      WRITE(6,'(2F10.6)') ASTL2(IND),LOG(SIGMAA(IND))
      ELSEIF (WSUM(IND).GT.RSMALL) THEN
      WRITE(6,'(F10.6,A)') ASTL2(IND),' infinity'
      END IF
      END DO
C
C LN(SIGMAA) VS STHOL**2 SHOULD GIVE STRAIGHT LINE WITH
C SLOPE = -26.3189*(MEAN SQUARE COORDINATE ERROR) AND
C INTERCEPT = 0.5*LN(SIGMAP/SIGMAN).
C THE FIRST AND LAST FEW POINTS USUALLY DEVIATE
C AND ARE IGNORED.
C
      END IF
C
C
C SMOOTH OUT ANY NEAR-ZERO VALUES OF SIGMAA, EXCEPT AT ENDS.
C
      NSM=0
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL.AND.SIGOLD(IND).LE.0.02D0) THEN
      NSM=NSM+1
      IF (WRNLEV.GE.10) THEN
      IF (NSM.LE.1) THEN
      WRITE(6,'(A)')
     & ' XSIGMAA: one or more values of SIGMAA are close to zero.',
     & '          This is probably due to numerical instabilities.',
     & '          Values interpolated from curve should be more valid.'
      END IF
      WRITE(6,'(A,/,A,F7.4,A,I3)')
     & ' XSIGMAA: attempting to interpolate SIGMAA curve to replace',
     & '          value of ',SIGOLD(IND),' for bin',IND
      END IF
      JLOW=IND
      DO WHILE (SIGOLD(JLOW).LE.0.02D0.AND.JLOW.GT.1)
      JLOW=JLOW-1
      END DO
      IF (SIGOLD(JLOW).LE.0.02D0) THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' Not possible to interpolate.'
      END IF
      ELSE
C
C WE'VE FOUND A POINT AT LOWER RESOLUTION.  NOW WE NEED ONE AT
C HIGHER RESOLUTION.
      JHIGH=IND
      DO WHILE (SIGOLD(JHIGH).LE.0.02D0.AND.JHIGH.LT.MBINS)
      JHIGH=JHIGH+1
      END DO
      IF (SIGOLD(JHIGH).LE.0.02D0) THEN
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') ' Not possible to interpolate.'
      END IF
      ELSE
C
C NOW WE CAN INTERPOLATE.
C
      SLOPE=(SIGOLD(JHIGH)-SIGOLD(JLOW))/(ASTL2(JHIGH)-ASTL2(JLOW))
      SIGMAA(IND)=SIGOLD(JLOW)+SLOPE*(ASTL2(IND)-ASTL2(JLOW))
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A,I6,A,I6,A,F7.4)')
     &  ' XSIGMAA: Interpolate between bins',JLOW,' and',JHIGH,
     &  ' to get new SIGMAA value of ',SIGMAA(IND)
      END IF
      END IF
      END IF
      END IF
      END DO
C
C store SIGMAA in stack
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
      VSTACK(REFLCT,VLEVEL)=DCMPLX(SIGMAA(ISHELL(REFLCT)),ZERO)
      END IF
      END DO
C
C print final
      IF (WRNLEV.GE.10) THEN
C
C write SIGMAA as a function of S^2
      WRITE(6,'(A)') ' Final Ln(SIGMAA) vs. S^2'
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL.AND.SIGMAA(IND).GT.RSMALL) THEN
      WRITE(6,'(2F10.6)') ASTL2(IND),LOG(SIGMAA(IND))
      ELSEIF (WSUM(IND).GT.RSMALL) THEN
      WRITE(6,'(F10.6,A)') ASTL2(IND),' infinity'
      END IF
      END DO
      END IF
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XSIM(X,SIM)
C
C     CALCULATE SIM AND SRINIVASAN NON-CENTRIC FIGURE OF MERIT
C     AS I1(X)/I0(X), WHERE I1 AND I0 ARE THE MODIFIED 1ST AND ZERO
C     ORDER BESSEL FUNCTIONS.
C     REFERENCES: SIM, G. A. (1960) ACTA CRYST. 13, 511-512;
C       SRINIVASAN, R. (1966) ACTA CRYST. 20, 143-144;
C       ABRAMOWITZ & STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 378.
C
C Authors: W. Kabsch, R.J. Read, A.T. Brunger
C
      IMPLICIT NONE
C I/O
      DOUBLE PRECISION X, SIM
C local
      DOUBLE PRECISION T
C begin
      T = ABS(X)/3.75D0
      IF (T .LE. 1.0D0) THEN
      T = T*T
      SIM=X*(0.5D0+T*(0.87890594D0+T*(0.51498869D0+T*(0.15084934D0+
     &           T*(0.02658733D0+T*(0.00301532D0+T*0.00032411D0))))))/
     &      (1.0D0+T*(3.5156229D0 +T*(3.0899424D0 +T*(1.2067492D0 +
     &           T*(0.2659732D0 +T*(0.0360768D0 +T*0.0045813D0 ))))))
      ELSE
      T = 1.0D0 /T
      SIM=(0.39894228D0+T*(-0.03988024D0+
     &           T*(-0.00362018D0+T*( 0.00163801D0+
     &           T*(-0.01031555D0+T*( 0.02282967D0+T*(-0.02895312D0+
     &           T*( 0.01787654D0-T*0.00420059D0))))))))*SIGN(1.0D0,X)/
     &    (0.39894228D0+T*( 0.01328592D0+
     &            T*( 0.00225319D0+T*(-0.00157565D0+
     &            T*( 0.00916281D0+T*(-0.02057706D0+T*( 0.02635537D0+
     &            T*(-0.01647633D0+T*0.00392377D0))))))))
      END IF
      RETURN
      END
C
C
C======================================================================
      SUBROUTINE XDOSIGMACV(VLEVEL,VMAX,VSTACK,LSTACK,
     &                      N,INDEX,XRNSYM,QHERM,
     &                      XRTR,XRMREF,HPH,HPK,HPL,
     &                      HPMULT,HPTYPE,MBINS,XBINLOW,XBINHIGH,
     &                      BINSHELL,ARGS)
C
C Returns sigmaa coefficient
C using line smoothing and conjugate gradients minimization.
C
C Suitable for cross-validation where EOBS and ECALC are
C calculated from the selected "test" set.
C
C Function call:
C     SIGACV[SIGMA=<real>](EOBS,ECALC)
C EOBS and ECALC must be normalized structure factor amplitudes
C SIGMA is the sigma for the smoothing restraints.
C
C EOBS: real
C ECALC: real
C
C It is assumed that the appropriate reflections are given to this
C routine - no checking is carried out
C
C Paul Adams and Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xsigmaa.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRMREF, HPH, HPK, HPL
      INTEGER HPMULT, HPTYPE
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE COMPLEX ARGS(*)
C pointer
      INTEGER CI, CJ, CII, CJJ, CIJ
      INTEGER SIGMAA, SUMM, SIGOLD, SHELL, SUMW
C
C parameter
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C get weight if an argument
      IF (ABS(DBLE(ARGS(1))).LT.RSMALL) THEN
         XSUSRW=ZERO
      ELSE
C convert sigma value given to weight (1/sigma^2)
         XSUSRW=ONE/(DBLE(ARGS(1))*DBLE(ARGS(1)))
      END IF
C
      XSREF=N
C
C store value of HPTYPE pointer into xsigmaa.inc datastructure
C for access by the minimizer target function
      HPPTYPE=HPTYPE
C normalize structure factor
      CI=ALLHP(IREAL8(MBINS))
      CII=ALLHP(IREAL8(MBINS))
      CJ=ALLHP(IREAL8(MBINS))
      CJJ=ALLHP(IREAL8(MBINS))
      CIJ=ALLHP(IREAL8(MBINS))
      XSWSUM=ALLHP(IREAL8(MBINS))
      XSSHLL=ALLHP(INTEG4(N))
      XSEOBS=ALLHP(IREAL8(N))
      XSECAL=ALLHP(IREAL8(N))
      XSINDX=ALLHP(INTEG4(N))
      SIGMAA=ALLHP(IREAL8(MBINS))
      SUMM=ALLHP(IREAL8(MBINS))
      SIGOLD=ALLHP(IREAL8(MBINS))
      SHELL=ALLHP(IREAL8(MBINS+1))
      XSASTL2=ALLHP(IREAL8(MBINS))
      SUMW=ALLHP(IREAL8(MBINS))
      XGOOD=ALLHP(ILOGIC(MBINS))
      PSIGMAA=ALLHP(IREAL8(MBINS))
C
      CALL XDOSIGMACV2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRNSYM,QHERM,XRTR,XRMREF,HEAP(HPH),
     &         HEAP(HPK),HEAP(HPL),HEAP(HPMULT),HEAP(HPTYPE),
     &         MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &         HEAP(CI),HEAP(CII),HEAP(CJ),HEAP(CJJ),
     &         HEAP(CIJ),HEAP(XSWSUM),HEAP(XSSHLL),
     &         HEAP(XSEOBS), HEAP(XSECAL), HEAP(XSINDX),
     &         HEAP(SIGMAA),HEAP(SUMM),HEAP(SIGOLD),HEAP(SHELL),
     &         HEAP(XSASTL2),HEAP(SUMW),HEAP(XGOOD),HEAP(PSIGMAA))
C
      CALL FREHP(PSIGMAA,IREAL8(MBINS))
      CALL FREHP(XGOOD,ILOGIC(MBINS))
      CALL FREHP(SUMW,IREAL8(MBINS))
      CALL FREHP(XSASTL2,IREAL8(MBINS))
      CALL FREHP(SHELL,IREAL8(MBINS+1))
      CALL FREHP(SIGOLD,IREAL8(MBINS))
      CALL FREHP(SUMM,IREAL8(MBINS))
      CALL FREHP(SIGMAA,IREAL8(MBINS))
      CALL FREHP(XSINDX,INTEG4(N))
      CALL FREHP(XSECAL,IREAL8(N))
      CALL FREHP(XSEOBS,IREAL8(N))
      CALL FREHP(XSSHLL,INTEG4(N))
      CALL FREHP(XSWSUM,IREAL8(MBINS))
      CALL FREHP(CJJ,IREAL8(MBINS))
      CALL FREHP(CJ,IREAL8(MBINS))
      CALL FREHP(CIJ,IREAL8(MBINS))
      CALL FREHP(CII,IREAL8(MBINS))
      CALL FREHP(CI,IREAL8(MBINS))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDOSIGMACV2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRNSYM,QHERM,XRTR,XRMREF,XRH,XRK,XRL,
     &         MULT,TYPE,MBINS,XBINLOW,XBINHIGH,BINSHELL,
     &         CI,CII,CJ,CJJ,CIJ,WSUM,
     &         ISHELL,EOBS,ECAL,LIND,SIGMAA,SUMM,SIGOLD,SHELL,
     &         ASTL2,SUMW,GOOD,XSIGMAA)
C
C
C References: R.J. Read, Improved fourier coefficients for maps
C  using phases from partial structures with errors, Acta Cryst. A42,
C 140--149 (1986).
C
C R.J. Read, Structure-factor probabilities for related structures,
C Acta Cryst A46, 900--912 (1990).
C
C N.S. Pannu and R.J. Read, Improved structure refinement through
C maximum likelihood. Acta Cryst A52, 659-668 (1996).
C
C Paul Adams, Axel T. Brunger, Navraj S. Pannu, and Randy Read
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xsigmaa.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRNSYM
      LOGICAL QHERM
      DOUBLE PRECISION XRTR(3,3)
      INTEGER XRMREF, XRH(*), XRK(*), XRL(*)
      INTEGER MULT(*), TYPE(*)
      INTEGER MBINS
      DOUBLE PRECISION XBINHIGH, XBINLOW, BINSHELL(*)
      DOUBLE PRECISION CI(*), CII(*), CJ(*), CJJ(*), CIJ(*), WSUM(*)
      INTEGER ISHELL(*)
      DOUBLE PRECISION EOBS(*), ECAL(*)
      DOUBLE PRECISION SIGMAA(*), SUMM(*), SIGOLD(*)
      DOUBLE PRECISION SHELL(*), ASTL2(*), SUMW(*)
      INTEGER LIND(*)
      LOGICAL GOOD(*)
      DOUBLE PRECISION XSIGMAA(*)
C local
      INTEGER REFLCT, H, K, L, IND
      DOUBLE PRECISION SSQ, EO2, EC2, CSUM, DSUM
      DOUBLE PRECISION CIT, CIIT, CJT, CJJT, CIJT, WSUMT, CCTOT
      DOUBLE PRECISION FOML, OVERM
      DOUBLE PRECISION X
      DOUBLE COMPLEX DCPLAC
      DOUBLE PRECISION DPPLAC
      DOUBLE PRECISION MLL
      INTEGER IER, I, II1, II2
      DOUBLE PRECISION W1, W2
C pointer
      INTEGER GRDENT, WORK
C target function
      EXTERNAL XSFUNC
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, HALF, FOUR, EIGHT
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (HALF=0.5D0, FOUR=4.0D0, EIGHT=8.0D0)
C begin
C
      VLEVEL=VLEVEL-1
C
      IF (N.GT.0) THEN
C
      IF (WRNLEV.GE.10) THEN
      IF (XBINLOW.GT.RSMALL) THEN
      WRITE(6,'(A,I8,A,F10.4,A,F10.4,A)')
     & ' XSIGMACV: using ',MBINS,' bins between',
     &  ONE/XBINHIGH,' and ',ONE/XBINLOW,' A resolution.'
      ELSE
      WRITE(6,'(A,I8,A,F10.4,A)')
     & ' XSIGMACV: using ',MBINS,' bins between',
     &  ONE/XBINHIGH,' and infinity A resolution.'
      END IF
      END IF
C
C initialize averages
      CIT=ZERO
      CJT=ZERO
      CIIT=ZERO
      CJJT=ZERO
      CIJT=ZERO
      WSUMT=ZERO
C
      DO IND=1,MBINS
      WSUM(IND)=ZERO
      CI(IND)=ZERO
      CJ(IND)=ZERO
      CII(IND)=ZERO
      CJJ(IND)=ZERO
      CIJ(IND)=ZERO
      ASTL2(IND)=ZERO
      END DO
C
C compute shell bin for each reflection
      CALL XDOBINPP(N,INDEX,XRTR,XRH,XRK,XRL,MBINS,
     &              XBINLOW,XBINHIGH,BINSHELL,ISHELL)
C
C compute ASTL2
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
      H=XRH(INDEX(REFLCT))
      K=XRK(INDEX(REFLCT))
      L=XRL(INDEX(REFLCT))
      SSQ=  (XRTR(1,1)*H + XRTR(2,1)*K + XRTR(3,1)*L)**2
     &    + (XRTR(1,2)*H + XRTR(2,2)*K + XRTR(3,2)*L)**2
     &    + (XRTR(1,3)*H + XRTR(2,3)*K + XRTR(3,3)*L)**2
      ASTL2(ISHELL(REFLCT))=ASTL2(ISHELL(REFLCT))+SSQ
      END IF
      END DO
C
C normalization
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
      EOBS(REFLCT)=VSTACK(REFLCT,VLEVEL)
      ECAL(REFLCT)=VSTACK(REFLCT,VLEVEL+1)
      LIND(REFLCT)=INDEX(REFLCT)
      WSUM(ISHELL(REFLCT))=WSUM(ISHELL(REFLCT))+ONE
      EO2=EOBS(REFLCT)**2
      EC2=ECAL(REFLCT)**2
      CI(ISHELL(REFLCT))=CI(ISHELL(REFLCT))+EO2
      CJ(ISHELL(REFLCT))=CJ(ISHELL(REFLCT))+EC2
      CII(ISHELL(REFLCT))=CII(ISHELL(REFLCT))+EO2**2
      CJJ(ISHELL(REFLCT))=CJJ(ISHELL(REFLCT))+EC2**2
      CIJ(ISHELL(REFLCT))=CIJ(ISHELL(REFLCT))+EO2*EC2
C
      WSUMT=WSUMT+ONE
      CIT=CIT+EO2
      CJT=CJT+EC2
      CIIT=CIIT+EO2**2
      CJJT=CJJT+EC2**2
      CIJT=CIJT+EO2*EC2
C
      END IF
      END DO
C
C compute correlation coefficients and average S^2 for each bin
C
C initial sigmaa is the square-root of the bin-wise correlation
C coefficient between Es.
      DO IND=1,MBINS
      IF (WSUM(IND).GT.THREE) THEN
      DSUM=(CII(IND) - CI(IND)**2/WSUM(IND)) *
     &     (CJJ(IND) - CJ(IND)**2/WSUM(IND))
      CSUM=CIJ(IND) - CI(IND)*CJ(IND)/WSUM(IND)
      IF (DSUM.GT.RSMALL) THEN
      DSUM=SQRT(DSUM)
C
      SIGMAA(IND)=CSUM/DSUM
C
C SIGMAA is set to zero if initially too small
      IF (SIGMAA(IND).LT.RSMALL) THEN
      SIGMAA(IND)=ZERO
      ELSE
      SIGMAA(IND)=SQRT(SIGMAA(IND))
      SIGMAA(IND)=MAX(SIGMAA(IND),ZERO)
      END IF
C
      ELSE
      SIGMAA(IND)=ZERO
      END IF
      ASTL2(IND)=ASTL2(IND)/(FOUR*WSUM(IND))
      ELSE
      SIGMAA(IND)=ZERO
      ASTL2(IND)=(BINSHELL(IND)+BINSHELL(IND+1))/EIGHT
      END IF
      SIGOLD(IND)=SIGMAA(IND)
      GOOD(IND)=(SIGMAA(IND).GT.ZERO)
      END DO
C
      IF (WSUMT.GT.RSMALL) THEN
      DSUM=(CIIT-CIT**2/WSUMT)*(CJJT-
     &       CJT**2/WSUMT)
      CSUM=CIJT - CIT*CJT/WSUMT
      IF (DSUM.GT.RSMALL) THEN
      DSUM=SQRT(DSUM)
      CCTOT=CSUM/DSUM
      ELSE
      CCTOT=ZERO
      END IF
      ELSE
      CCTOT=ZERO
      END IF
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(2A,F10.5)') ' XSIGMAACV: overall correlation between',
     &  ' squared arguments: ',CCTOT
      END IF
C
C Interpolate any bins with sigmaa estimates of zero
      DO I=1,MBINS
      IF (.NOT.GOOD(I)) THEN
C Find two nearest neighbours with real data
C
      CALL NAYBRS(GOOD,MBINS,I,II1,II2)
C
C Determine the linear combination of the neighbours that would
C interpolate (or extrapolate, if necessary) the current value
C as a function of astl2.
C
      W1=(ASTL2(II2)-ASTL2(I))/(ASTL2(II2)-ASTL2(II1))
C
C Don't allow extrapolation to higher resolution to increase
C sigmaa values; in that case, use the sigmaa value from the
C nearest neighbour.
C
      IF ((II1.GT.I).AND.(II2.GT.I).AND.(SIGOLD(II1).GT.SIGOLD(II2)))
     &    W1=ONE
      W2=ONE- W1
      SIGMAA(I)=W1*SIGOLD(II1) + W2*SIGOLD(II2)
      SIGMAA(I)=MAX(ZERO,MIN(0.999D0,SIGMAA(I)))
      ENDIF
      END DO
C
C Refine estimates of SIGMAA using conjugate gradients
C using the powell minimizer, restraint to line applied
C if weight greater than zero
      GRDENT=ALLHP(IREAL8(MBINS))
      WORK=ALLHP(IREAL8(MBINS*6))
C
C calculate weight between log likelihood and line restraint
      XSWGHT=ONE
      CALL XSWCAL(MBINS,SIGMAA,MLL,HEAP(GRDENT))
C
C do finite difference test if wrnlev is very high
      IF (WRNLEV.GE.14) THEN
      CALL XSDBG(SIGMAA,HEAP(GRDENT),HEAP(WORK),0.00001D0,MBINS)
      END IF
C
C do the minimization
C
C
C initialize the lowest target value
      XSTARG=R4BIG
      XSTCYC=0
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(2A)')' --------------------------------',
     & '------------------------------'
      WRITE(6,'(A)')' Starting cross-validated SIGMAA refinement'
      WRITE(6,'(2A)')' --------------------------------',
     & '------------------------------'
      END IF
C
      XSTCYC=0
      CALL ZXCGR(XSFUNC,MBINS,XSTOL,XSNSTP,XSSTSZ,SIGMAA,
     &           HEAP(GRDENT),MLL,HEAP(WORK),IER)
C
C set parameters to the values that produces the lowest target value
      DO IND=1,MBINS
      SIGMAA(IND)=XSIGMAA(IND)
      END DO
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(2A)')' --------------------------------',
     & '------------------------------'
      WRITE(6,'(A,I6,A)')
     & ' XDOSIGMACV: using the parameters (at cycle ', XSTCYC,
     & ') which produce the minimum target. '
      END IF
C
      CALL FREHP(WORK,IREAL8(MBINS*6))
      CALL FREHP(GRDENT,IREAL8(MBINS))
C
C check Error Conditions
      IF (IER.EQ.0) THEN
         WRITE(6,'(A)')' XSIGMAACV: gradient converged'
      ELSE IF (IER.EQ.129) THEN
         WRITE(6,'(A)')' XSIGMAACV: Line search terminated'
      ELSE IF (IER.EQ.130) THEN
         WRITE(6,'(A)')' %XSIGMAACV-ERR: Search Direction Uphill'
      ELSE IF (IER.EQ.131) THEN
         WRITE(6,'(A)')' XSIGMAACV: step limit reached'
      ELSE IF (IER.EQ.132) THEN
         WRITE(6,'(A)')' %XSIGMAACV-ERR: Failure to reduce MLL'
      END IF
C
C say its converged if step limit reached or gradient converged
      IF (IER.EQ.0) THEN
         WRITE(6,'(A)')
     &        ' XSIGMAACV: Refinement of SIGMAA has converged '
         CALL DECLAR( 'STATUS', 'ST', 'COMPLETE', DCPLAC, DPPLAC )
      ELSE IF (IER.EQ.131) THEN
         WRITE(6,'(A)')
     &        ' XSIGMAACV: Refinement of SIGMAA has finished '
         CALL DECLAR( 'STATUS', 'ST', 'INCOMPLETE', DCPLAC, DPPLAC )
      ELSE
         CALL WRNDIE(+1,'XSIGMAACV',
     &        'Refinement of SIGMAA has not converged.')
         CALL DECLAR( 'STATUS', 'ST', 'INCOMPLETE', DCPLAC, DPPLAC )
      END IF
C
      OVERM=ZERO
C
C calculate overall FOM
C loop through all selected reflections
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
C
C compute X
      X = TWO*SIGMAA(ISHELL(REFLCT))*
     &          EOBS(REFLCT)*
     &          ECAL(REFLCT)/
     &   (ONE-SIGMAA(ISHELL(REFLCT))**2)
C
C CALCULATE FIGURE OF MERIT, GIVEN X. FORMULA DIFFERS FOR
C CENTRIC AND NON-CENTRIC.  SRINIVASAN (ACTA CRYST 20:143(1966))
C DIFFERS FROM WOOLFSON (ACTA CRYST 9:804(1956)) AND SIM (ACTA
C CRYST 13:511(1960)) ONLY IN THE EXPRESSION FOR X, WHICH IN
C THE PRESENT CASE INCLUDES CONTRIBUTIONS FROM MISSING STRUCTURE
C AND MODEL ERROR.
C FOR NON-CENTRIC REFLECTIONS, USE THE FUNCTION SIM WHICH IS
C THE RATIO OF THE MODIFIED 1ST AND ZERO ORDER BESSEL FUNCTIONS.
C
      IF (TYPE(INDEX(REFLCT)).GE.1) THEN
C acentric reflections
      CALL XSIM(X,FOML)
      ELSE
C centric reflections
      FOML=TANH(X/TWO)
      END IF
C
      OVERM=OVERM+FOML
C
      END IF
      END DO
C
      IF (WSUMT.GT.RSMALL) THEN
         OVERM=OVERM/WSUMT
      ELSE
         OVERM=ZERO
      END IF
      WRITE(6,'(A,F10.5)')
     & '            Overall mean FOM is',OVERM
C
C store SIGMAA in stack
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
      VSTACK(REFLCT,VLEVEL)=DCMPLX(SIGMAA(ISHELL(REFLCT)),ZERO)
      ELSE
      VSTACK(REFLCT,VLEVEL)=DCMPLX(ZERO,ZERO)
      END IF
      END DO
C
C print final information
      IF (WRNLEV.GE.10) THEN
C
C write SIGMAA as a function of S^2
      IF (XSWGHT.GT.ZERO) THEN
      WRITE(6,'(2A,F12.4)') ' Sigmaa refined:  ',
     &   'restrained minus log likelihood= ',MLL
      ELSE
      WRITE(6,'(2A,F12.4)') ' Sigmaa refined:  ',
     &   'minus log likelihood= ',MLL
      END IF
      WRITE(6,'(A)') ' ========================'
      WRITE(6,'(A)') ' Final Ln(SIGMAA) vs. S^2'
      DO IND=1,MBINS
      IF (WSUM(IND).GT.RSMALL.AND.SIGMAA(IND).GT.RSMALL) THEN
      WRITE(6,'(2F10.6)') ASTL2(IND),LOG(SIGMAA(IND))
      ELSEIF (WSUM(IND).GT.RSMALL) THEN
      WRITE(6,'(F10.6,A)') ASTL2(IND),' infinity'
      END IF
      END DO
      END IF
C
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE NAYBRS(GOOD,NXY,II,II1,II2)
C
C     Randy J. Read and Navraj S. Pannu
C
C     Find two neighbours for interpolation or (if necessary)
C     extrapolation
C
      IMPLICIT NONE
C I/O
      INTEGER NXY,II,II1,II2
      LOGICAL GOOD(NXY)
C local
      INTEGER IM1,IM2,IP1,IP2
      INTEGER I
C
C First find the two closest points with data on the low side.
C
      IM1=0
      IM2=0
      IF (II .NE. NXY) THEN
      DO 50 I=II+1,NXY
      IF (.NOT. GOOD(I)) GO TO 50
      IF (IM1 .EQ. 0) THEN
      IM1=I
      ELSE
      IM2=I
C
C Now we've found two
C
      GO TO 60
      ENDIF
50    CONTINUE
60    CONTINUE
      ENDIF
C
C Now find the two closest points with data on the high side.
C
      IP1=0
      IP2=0
      IF (II .NE. 1) THEN
      DO 80 I=II-1,1,-1
      IF (.NOT. GOOD(I)) GO TO 80
      IF (IP1 .EQ. 0) THEN
      IP1=I
      ELSE
      IP2=I
      GO TO 90
      ENDIF
80    CONTINUE
90    CONTINUE
      ENDIF
C
C Finally choose the points to interpolate or, if necessary, to
C extrapolate
C
      II1=0
      II2=0
      IF (IM1.NE.0) THEN
      II1=IM1
      IF (IP1.NE.0) THEN
      II2=IP1
      ELSE
      IF (IM2.NE.0) II2=IM2
      ENDIF
      ELSE
      IF (IP1.NE.0) THEN
      II1=IP1
      IF (IP2.NE.0) II2=IP2
      ENDIF
      ENDIF
      IF (II1.EQ.0.OR.II2.EQ.0) THEN
      WRITE (6,'(A)')' Problems interpolating or extrapolating data.'
      WRITE (6,'(A)')' Input data had the following distribution:'
      DO I=1,NXY
         IF (GOOD(I)) THEN
            WRITE (6,'(A,I5,A)')'Data point ',I,' could be used'
         ELSE
            WRITE (6,'(A,I5,A)')'Data point ',I,' could not be used'
         END IF
      END DO
      CALL WRNDIE(-5,'NAYBRS',
     &  'Less than two data points for interpolation/extrapolation')
      ENDIF
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XSDBG(PARAM,GRDENT,WORK,STEP,DIM)
C
C perform finite difference test on SIGMAA derivatives wrt
C minus log likelihood
C
C Paul Adams and Axel Brunger
C
      IMPLICIT NONE
C
      INCLUDE 'xsigmaa.inc'
C
      DOUBLE PRECISION PARAM(*), GRDENT(*), WORK(*), STEP
      INTEGER DIM
C local
      INTEGER I
      DOUBLE PRECISION TARGET, NUMGRD
C parameter
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
C
      DO I=1,DIM
         PARAM(I)=PARAM(I)+STEP
         XSCYCLE=-1
         CALL XSFUNC(DIM,PARAM,TARGET,GRDENT)
         NUMGRD=TARGET
         PARAM(I)=PARAM(I)-TWO*STEP
         XSCYCLE=-1
         CALL XSFUNC(DIM,PARAM,TARGET,GRDENT)
         NUMGRD=NUMGRD-TARGET
         NUMGRD=NUMGRD/(TWO*STEP)
         PARAM(I)=PARAM(I)+STEP
         XSCYCLE=-1
         CALL XSFUNC(DIM,PARAM,TARGET,GRDENT)
         WRITE(6,'(A,I6,A,F12.5,A,F12.5,A,F12.5)')
     &     ' XSDBG: Bin # ',I,' analyt.=',GRDENT(I),
     &     ' finite=',NUMGRD
      END DO
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XSWCAL(DIM,PARAM,TARGET,GRDENT)
C
C Calculate the weight that balances the RMS gradients of the
C minus log likelihood and the line restraint
C This weight is then multiplied by the user provided scale factor
C
C Paul Adams and Axel Brunger
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xsigmaa.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
C
      INTEGER DIM
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
C pointers
      INTEGER GRDMLL, GRDRES
C
      GRDMLL=ALLHP(IREAL8(DIM))
      GRDRES=ALLHP(IREAL8(DIM))
C
      XSCYCLE=-1
      XSTARG=R4BIG
C
      CALL XSFUNC2(HEAP(XSEOBS),HEAP(XSECAL),XSREF,HEAP(XSINDX),
     &             HEAP(HPPTYPE),HEAP(XSSHLL),HEAP(XSASTL2),
     &             HEAP(XSWSUM),HEAP(XGOOD),HEAP(GRDMLL),HEAP(GRDRES),
     &             XSWGHT,XSMIN,XSMAX,XSCYCLE,
     &             PARAM,DIM,TARGET,GRDENT,XSTARG,XSTCYC,
     &             HEAP(PSIGMAA))
C
      CALL XSWCAL2(DIM,HEAP(GRDMLL),HEAP(GRDRES),XSWGHT)
C
      IF (XSUSRW.LT.0) THEN
         XSWGHT=XSWGHT*ABS(XSUSRW)
      ELSE
         XSWGHT=XSUSRW
      END IF
      WRITE(6,'(A,F16.8)')
     & ' XSIGMAACV: Overall sigma for line restraint is ',
     & SQRT(1.0D0/XSWGHT)
C
      CALL FREHP(GRDRES,IREAL8(DIM))
      CALL FREHP(GRDMLL,IREAL8(DIM))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XSWCAL2(DIM,GRDMLL,GRDRES,WEIGHT)
C
C see above
C
C Paul Adams and Axel Brunger
C
      IMPLICIT NONE
C
      INTEGER DIM
      DOUBLE PRECISION GRDMLL(*), GRDRES(*), WEIGHT
C local
      INTEGER I
      DOUBLE PRECISION RMSMLL,RMSRES
C parameters
      DOUBLE PRECISION ZERO, ONE, RSMALL
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, RSMALL=1.0D-5)
C
C calculate rms gradients
      RMSMLL=ZERO
      DO I=1,DIM
         RMSMLL=RMSMLL+GRDMLL(I)**2
      END DO
      RMSMLL=SQRT(RMSMLL/DIM)
C
      RMSRES=ZERO
      DO I=1,DIM
         RMSRES=RMSRES+GRDRES(I)**2
      END DO
      RMSRES=SQRT(RMSRES/DIM)
C
      IF (RMSRES.GT.RSMALL) THEN
         WEIGHT=RMSMLL/RMSRES
      ELSE
         WEIGHT=ONE
      END IF
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XSFUNC(DIM,PARAM,TARGET,GRDENT)
C
C Function to calculate the minus log likelihood as a function of SIGMAA
C and the derivatives of the likelifood wrt SIGMAA
C A separate term is used to restrain the SIGMAA values to be linear
C
C Paul Adams and Axel Brunger
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'xsigmaa.inc'
      INCLUDE 'funct.inc'
C
      INTEGER DIM
      DOUBLE PRECISION PARAM(*), TARGET, GRDENT(*)
C pointers
      INTEGER GRDMLL, GRDRES
C
      GRDMLL=ALLHP(IREAL8(DIM))
      GRDRES=ALLHP(IREAL8(DIM))
C
      CALL XSFUNC2(HEAP(XSEOBS),HEAP(XSECAL),XSREF,HEAP(XSINDX),
     &             HEAP(HPPTYPE),HEAP(XSSHLL),HEAP(XSASTL2),
     &             HEAP(XSWSUM),HEAP(XGOOD),HEAP(GRDMLL),HEAP(GRDRES),
     &             XSWGHT,XSMIN,XSMAX,XSCYCLE,
     &             PARAM,DIM,TARGET,GRDENT,
     &             XSTARG,XSTCYC,HEAP(PSIGMAA))
C
      CALL FREHP(GRDRES,IREAL8(DIM))
      CALL FREHP(GRDMLL,IREAL8(DIM))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XSFUNC2(EOBS,ECAL,N,INDEX,TYPE,ISHELL,ASTL2,
     &                   WSUM,GOOD,GRDMLL,GRDRES,WEIGHT,
     &                   XSMIN,XSMAX,XSCYCLE,
     &                   SIGMAA,MBINS,TARGET,GRDENT,
     &                   XSTARG,XSTCYC,XSIGMAA)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C Accumulate minus log likelihood and gradients
C of minus log likelihood w.r.t. SIGMAA values
C and linear restraint penalty and its derivatives w.r.t sigmaa
C
      IMPLICIT NONE
C
      INCLUDE 'timer.inc'
C
      INTEGER MBINS
      DOUBLE PRECISION WSUM(*), ASTL2(*)
      LOGICAL GOOD(*)
      DOUBLE PRECISION SIGMAA(*), GRDENT(*), GRDMLL(*), GRDRES(*)
      DOUBLE PRECISION XSMIN, XSMAX, TARGET
      DOUBLE PRECISION WEIGHT
      INTEGER N, XSCYCLE
      DOUBLE PRECISION EOBS(*), ECAL(*)
      INTEGER INDEX(*), TYPE(*), ISHELL(*)
      DOUBLE PRECISION XSTARG
      INTEGER XSTCYC
      DOUBLE PRECISION XSIGMAA(*)
C local
      INTEGER I,REFLCT,II1,II2
      DOUBLE PRECISION EOSQ, ECSQ, WT, SIGA, SIGASQ, FOM
      DOUBLE PRECISION VAR, X, GRDRMS, MLL, RES, FOMT
      DOUBLE PRECISION WGHT, W1, W2, SIGINT, DEV, DEVSQ
C functions
      DOUBLE PRECISION LOGI0, LOGCH
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, HALF, FOUR, SMALL
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (HALF=0.5D0, FOUR=4.0D0, SMALL=1.0D-5)
C
      DO I=1,MBINS
         SIGMAA(I)=MIN(MAX(SIGMAA(I),XSMIN),XSMAX)
      END DO
C
      MLL=ZERO
      DO I=1,MBINS
         GRDMLL(I)=ZERO
      END DO
C
      FOMT=ZERO
C
      DO REFLCT=1,N
      IF (ISHELL(REFLCT).GT.0) THEN
         EOSQ=EOBS(REFLCT)**2
         ECSQ=ECAL(REFLCT)**2
         IF (TYPE(INDEX(REFLCT)).GE.1) THEN
            WT=TWO
         ELSE
            WT=ONE
         END IF
         SIGA=SIGMAA(ISHELL(REFLCT))
         SIGASQ=SIGA**2
         VAR=ONE-SIGASQ
         X=TWO*SIGA*EOBS(REFLCT)*ECAL(REFLCT)/VAR
         MLL=MLL+(WT/TWO)*(LOG(VAR)+(EOSQ+SIGASQ*ECSQ)/VAR)
C
         IF (TYPE(INDEX(REFLCT)).GE.1) THEN
C acentric reflections
            MLL=MLL-LOGI0(X)
            CALL XSIM(X,FOM)
         ELSE
C centric reflections
            MLL=MLL-LOGCH(X/TWO)
            FOM=TANH(X/TWO)
         END IF
C
         FOMT=FOMT+FOM
C
         GRDMLL(ISHELL(REFLCT))=GRDMLL(ISHELL(REFLCT))+
     &        (WT*(SIGA*(EOSQ+ECSQ-ONE+SIGASQ)-
     &         FOM*EOBS(REFLCT)*ECAL(REFLCT)*
     &        (ONE+SIGASQ))/VAR**2)
      END IF
      END DO
C
C overall FOM
      FOMT=FOMT/N
C
C calculate restraint penalty for line fit
      RES=ZERO
      DO I=1,MBINS
         GRDRES(I)=ZERO
      END DO
C
      WGHT=WEIGHT
      DO I=1,MBINS
C
C  Add in the smoothing terms to ll and gradient
C
C  First find the two nearest neighbours with real data
C
         CALL NAYBRS(GOOD,MBINS,I,II1,II2)
C
C  Determine the linear combination of the neighbours that would
C  interpolate (or extrapolate, if necessary) the current value
C  as a function of astl2.
C
         W1=(ASTL2(II2)-ASTL2(I))/(ASTL2(II2)-ASTL2(II1))
C
C  Don't allow extrapolation to higher resolution to increase
C  sigmaa values; in that case, use the sigmaa value from the
C  nearest neighbour.
C
         IF (     (II1.GT.I).AND.(II2.GT.I)
     &       .AND.(SIGMAA(II1).GT.SIGMAA(II2)))
     &     W1=ONE
         W2=ONE-W1
         SIGINT=W1*SIGMAA(II1) + W2*SIGMAA(II2)
         SIGINT=MAX(XSMIN,MIN(XSMAX,SIGINT))
C
C  If this is an empty shell, keep sigmaa equal to the
C  interpolated value at all times, instead of wasting time
C  refining it.
C
         IF (WSUM(I) .EQ. 0) SIGMAA(I) = SIGINT
C
         DEV=SIGMAA(I)-SIGINT
         DEVSQ=DEV**2
C
C  In the case of extrapolation, both ii1 and ii2 will be on the
C  same side of i.  In this case, multiply the target sigma by
C  two to put less weight on the restraint.
C
         IF ((II1-I)*(II2-I).GT.0.0) WEIGHT=WEIGHT/FOUR
         RES=RES+WEIGHT*DEVSQ
         GRDRES(I)=GRDRES(I)+TWO*WEIGHT*DEV
         GRDRES(II1)=GRDRES(II1)-W1*TWO*WEIGHT*DEV
         GRDRES(II2)=GRDRES(II2)-W2*TWO*WEIGHT*DEV
         WEIGHT=WGHT
      END DO
C
C
C sum the two parts of the function - return in TARGET and GRDENT
      TARGET=MLL+RES
C
      DO I=1,MBINS
         GRDENT(I)=GRDMLL(I)+GRDRES(I)
      END DO
C
C correct for maximum and minimum values
      DO I=1,MBINS
         IF (SIGMAA(I).GE.XSMAX) THEN
            IF (GRDENT(I).LT.ZERO) GRDENT(I)=ZERO
         ELSE IF (SIGMAA(I).LE.XSMIN) THEN
            IF (GRDENT(I).GT.ZERO) GRDENT(I)=ZERO
         END IF
      END DO
C
C calculate rms gradient
      GRDRMS=ZERO
      DO I=1,MBINS
         GRDRMS=GRDRMS+GRDENT(I)**2
      END DO
      GRDRMS=SQRT(GRDRMS/MBINS)
C
      XSCYCLE=XSCYCLE+1
C store parameters for the cycle with the lowest target value
      IF (TARGET.LT.XSTARG) THEN
      XSTARG=TARGET
      XSTCYC=XSCYCLE
      DO I=1,MBINS
      XSIGMAA(I)=SIGMAA(I)
      END DO
      END IF
C
C print information
      IF (WRNLEV.GE.10) THEN
      IF (XSCYCLE.GE.1) THEN
      WRITE(6,'(A,I4,2(A,F12.4),A,F8.5)')
     & ' Cycle=',XSCYCLE,' Target=',TARGET,
     & ' E(grad)=',GRDRMS,' <FOM>=',FOMT
      WRITE(6,'(2(A,F12.4))')
     & '            E(MLL)=',MLL,' E(Line)=',RES
      END IF
      END IF
C
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION LOGI0(X)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
      IMPLICIT NONE
C
      DOUBLE PRECISION X
C functions
      DOUBLE PRECISION BESEI0
C
      LOGI0=LOG(BESEI0(X))+X
C
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION LOGCH(X)
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Copyright 1996  University of Alberta and Yale University
C
C calculates LN(COSH(X))
C this is an even function that quickly converges to
C ABS(X)-LN(2) as the magnitude of x increases
C
      IMPLICIT NONE
C
      DOUBLE PRECISION X
C parameter
      DOUBLE PRECISION TWO, TEN
      PARAMETER (TWO=2.0D0, TEN=10.0D0)
C
      IF (ABS(X).LT.TEN) THEN
         LOGCH=LOG(COSH(X))
      ELSE
         LOGCH=ABS(X)-LOG(TWO)
      END IF
C
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION BESEI0(X)
C
C Calculates approximate values for the modified Bessel function
C of the first kind of order zero multiplied by EXP(-ABS(X)).
C
C R.J. Read and N.S. Pannu
C
C Modified by Paul Adams and Axel Brunger
C
C Derived from:
C Authors: W. J. Cody and L. Stoltz
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
      IMPLICIT NONE
C
      DOUBLE PRECISION X
C local
      DOUBLE PRECISION AX, AX2, SUMN, SUMD
C parameter
      DOUBLE PRECISION ONE, XSMALL, FIFTEEN
      PARAMETER (ONE=1.0D0, FIFTEEN=15.0D0)
      PARAMETER (XSMALL=5.55D-17)
C
      AX=ABS(X)
C
      IF (AX.LT.XSMALL) THEN
         BESEI0=ONE
C
      ELSE IF (AX.LT.FIFTEEN) THEN
         AX2=AX*AX
         SUMN=     -2.2335582639474375249D+15 +
     &        AX2*(-5.5050369673018427753D+14 +
     &        AX2*(-3.2940087627407749166D+13 +
     &        AX2*(-8.4925101247114157499D+11 +
     &        AX2*(-1.1912746104985237192D+10 +
     &        AX2*(-1.0313066708737980747D+08 +
     &        AX2*(-5.9545626019847898221D+05 +
     &        AX2*(-2.4125195876041896775D+03 +
     &        AX2*(-7.0935347449210549190D+00 +
     &        AX2*(-1.5453977791786851041D-02 +
     &        AX2*(-2.5172644670688975051D-05 +
     &        AX2*(-3.0517226450451067446D-08 +
     &        AX2*(-2.6843448573468483278D-11 +
     &        AX2*(-1.5982226675653184646D-14 +
     &        AX2*(-5.2487866627945699800D-18))))))))))))))
         AX2=AX2-FIFTEEN**2
         SUMD=     -9.7087946179594019126D+14 +
     &        AX2*( 3.7604188704092954661D+12 +
     &        AX2*(-6.5626560740833869295D+09 +
     &        AX2*( 6.5158506418655165707D+06 +
     &        AX2*(-3.7277560179962773046D+03 + AX2))))
         BESEI0=(SUMN/SUMD)*(EXP(-AX))
C
      ELSE IF (AX.GE.FIFTEEN) THEN
         AX2=(ONE/AX)-(ONE/FIFTEEN)
         SUMN=     -2.1877128189032726730D-06 +
     &        AX2*( 9.9168777670983678974D-05 +
     &        AX2*(-2.6801520353328635310D-03 +
     &        AX2*(-3.7384991926068969150D-03 +
     &        AX2*( 4.7914889422856814203D-01 +
     &        AX2*(-2.4708469169133954315D+00 +
     &        AX2*( 2.9205384596336793945D+00 +
     &        AX2*(-3.9843750000000000000D-01)))))))
         SUMD=     -5.5194330231005480228D-04 +
     &        AX2*( 3.2547697594819615062D-02 +
     &        AX2*(-1.1151759188741312645D+00 +
     &        AX2*( 1.3982595353892851542D+01 +
     &        AX2*(-6.0228002066743340583D+01 +
     &        AX2*( 8.5539563258012929600D+01 +
     &        AX2*(-3.1446690275135491500D+01 + AX2))))))
         BESEI0=((SUMN/SUMD)+3.984375D-01)/SQRT(AX)
      END IF
C
      RETURN
      END
C
C =========================================================================
C
