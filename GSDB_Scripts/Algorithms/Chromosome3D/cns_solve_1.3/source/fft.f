C======================================================================
C FFT routines ========================================================
C =====================================================================
      SUBROUTINE FFTPRP(BASE,PRIME,AVOID)
C
C Returns maximum factor PRIME allowed in integer
C dimensions of 3-d FFT.  Also returns BASE which specifies
C that BASE should be a factor, e.g. if BASE is two then
C the integer has to be even.
C
C ********************************************************
C Universal Version
C ********************************************************
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'xfft.inc'
      INTEGER BASE, PRIME, AVOID
C begin
      BASE=1
      PRIME=5
      AVOID=2
C
      FTCRBX=1
      RETURN
      END
C
C=======================================================================
C
C=======================================================================
C BEGIN OF FFTPACK
C=======================================================================
C
C      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C      *                                                               *
C      *                           FFTPACK                             *
C      *                                                               *
C      *                                                               *
C      *     A package of Fortran subprograms for calculating          *
C      *                                                               *
C      *     fast Fourier transforms for both complex and real         *
C      *                                                               *
C      *      periodic sequences and certain other symmetric           *
C      *                                                               *
C      *             sequences that are listed below                   *
C      *                                                               *
C      *               (Version 4.1 November 1988)                     *
C      *                                                               *
C      *                             by                                *
C      *                                                               *
C      *                      Paul Swarztrauber                        *
C      *                                                               *
C      *                             of                                *
C      *                                                               *
C      *         The National Center for Atmospheric Research          *
C      *                                                               *
C      *                Boulder, Colorado  (80307)  U.S.A.             *
C      *                                                               *
C      *                   which is sponsored by                       *
C      *                                                               *
C      *              the National Science Foundation                  *
C      *                                                               *
C      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C  LATEST REVISION
C  ---------------
C  November 1988  (Version 4.1)
C
C  PURPOSE
C  -------
C  This package consists of programs which perform fast Fourier
C  transforms for both complex and real periodic sequences and
C  certain other symmetric sequences that are listed below.
C
C  USAGE
C  -----
C  1.   RFFTI     initialize  RFFTF and RFFTB
C  2.   RFFTF     forward transform of a real periodic sequence
C  3.   RFFTB     backward transform of a real coefficient array
C
C  4.   EZFFTI    initialize EZFFTF and EZFFTB
C  5.   EZFFTF    a simplified real periodic forward transform
C  6.   EZFFTB    a simplified real periodic backward transform
C
C  7.   SINTI     initialize SINT
C  8.   SINT      sine transform of a real odd sequence
C
C  9.   COSTI     initialize COST
C  10.  COST      cosine transform of a real even sequence
C
C  11.  SINQI     initialize SINQF and SINQB
C  12.  SINQF     forward sine transform with odd wave numbers
C  13.  SINQB     unnormalized inverse of SINQF
C
C  14.  COSQI     initialize COSQF and COSQB
C  15.  COSQF     forward cosine transform with odd wave numbers
C  16.  COSQB     unnormalized inverse of COSQF
C
C  17.  CFFTI     initialize CFFTF and CFFTB
C  18.  CFFTF     forward transform of a complex periodic sequence
C  19.  CFFTB     unnormalized inverse of CFFTF
C
C  SPECIAL CONDITIONS
C  ------------------
C  Before calling routines EZFFTB and EZFFTF for the first time,
C  or before calling EZFFTB and EZFFTF with a different length,
C  users must initialize by calling routine EZFFTI.
C
C  I/O
C  ---
C  None
C
C  PRECISION
C  ---------
C  None
C
C  REQUIRED LIBRARY FILES
C  ----------------------
C  None
C
C  LANGUAGE
C  --------
C  Fortran
C
C  HISTORY
C  -------
C  Developed at NCAR in Boulder, Colorado by Paul N. Swarztrauber
C  of the Scientific Computing Division.  Released on NCAR's public
C  software libraries in January 1980.
C     September 1973    Version 1
C     April     1976    Version 2
C     January   1978    Version 3
C     December  1979    Version 3.1
C     February  1985    Documentation upgrade
C     May       1985    Increased efficiency
C     November  1988    Version 4.1, Fortran 77 changes
C     June      1989    Add function PIMACH for distribution purposes
C
C  PORTABILITY
C  -----------
C  Fortran 77
C
C  INSTALLATION ON UNIX MACHINES
C  -----------------------------
C  A Makefile has been included in this package to facilitate
C  installation on unix systems.  If you use it, you may want to change
C  the LIBDIR = /usr/local/lib to the directory where you want to install
C  FFTPACK.  This Makefile will archive the testing routine with the
C  library archive file.
C
C  INSTALLATION ON non-UNIX MACHINES
C  ---------------------------------
C  To install the library, compile all of the provided .f files into
C  a relocatable binary library in the manner familiar to your system.
C
C  TESTING
C  -------
C  A subroutine TFFTPK has been supplied as a partial test of FFTPACK.
C  We emphasize the partial nature of this test; it exercises little
C  of the package and is only meant to be a simple check of the library
C  installation.  To use TFFTPK, you will need to write a stub program:
C
C        PROGRAM MAIN
C        CALL TFFTPK(IERROR)
C        STOP
C        END
C
C  Routine TFFTPK should return a 0 and print a success message,
C  or return a 1 and print a failure message otherwise.
C
C  REFERENCES
C  ----------
C
C (1) "Vectorizing the Fast Fourier Transforms", by Paul Swarztrauber,
C     Parallel Computations, G. Rodrigue, ed., Academic Press,
C     New York 1982.
C
C (2) "Fast Fourier Transforms Algorithms for Vector Computers", by
C     Paul Swarztrauber, Parallel Computing, (1984) pp.45-63.
C
C                         - End of FFTPACK document -
C
C=======================================================================
C Single and Double Precision Version of FFTPACK 4.1
C   http://www.scd.ucar.edu/softlib/FFTPACK.html (97-07-05)
C
C Modifications:
C   All files: implicit none
C   All variables declared explicitly
C   Put a @ in front of all subroutine names and function names (s-markup)
C   Check OK: No non-C lines longer than 72 characters (s-col72)
C   Check OK: After removal of @ *.m files equal *.f files (s-check-f-vs-m)
C   ezfftf.m:
C     C split next statement to enforce promotion of 2. to double precision
C     C     CF = 2./FLOAT(N)
C   tfftpk.m: replace specific intrinsic AMAX1 by generic intrinsic MAX
C   Generation of single and double precision files: s-single+double
C   Makefile changed: compile single and double precision files
C   .tree: single and double precision file names
C=======================================================================
C     SUBROUTINE DCFFTB(N,C,WSAVE)
C
C     SUBROUTINE DCFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
C     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , DCFFTB COMPUTES
C     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
C
C     A CALL OF DCFFTF FOLLOWED BY A CALL OF DCFFTB WILL MULTIPLY THE
C     SEQUENCE BY N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE DCFFTB MUST BE
C     INITIALIZED BY CALLING SUBROUTINE DCFFTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C
C     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
C            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
C
C     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C
C     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
C             IN THE PROGRAM THAT CALLS DCFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE DCFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY DCFFTF AND DCFFTB.
C
C     OUTPUT PARAMETERS
C
C     C      FOR J=1,...,N
C
C                C(J)=THE SUM FROM K=1,...,N OF
C
C                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
C
C                            WHERE I=SQRT(-1)
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF SUBROUTINE DCFFTF OR DCFFTB
C
      SUBROUTINE DCFFTB (N,C,WSAVE)
CIMPL DIMENSION       C(*)       ,WSAVE(*)
      implicit none
      integer           n
      double precision  c(*)
      double precision  wsave(*)
C local variables
      integer           iw1
      integer           iw2
C
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C=======================================================================
      SUBROUTINE DCFFTB1 (N,C,CH,WA,IFAC)
CIMPL DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      implicit none
      integer           n
      double precision  c(*)
      double precision  ch(*)
      double precision  wa(*)
      integer           ifac(*)
C local variables
      integer           i
      integer           idl1
      integer           ido
      integer           idot
      integer           ip
      integer           iw
      integer           ix2
      integer           ix3
      integer           ix4
      integer           k1
      integer           l1
      integer           l2
      integer           n2
      integer           na
      integer           nac
      integer           nf
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL DPASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL DPASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL DPASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL DPASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL DPASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL DPASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL DPASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL DPASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL DPASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL DPASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C=======================================================================
C     SUBROUTINE DCFFTI(N,WSAVE)
C
C     SUBROUTINE DCFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH DCFFTF AND DCFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH DCFFTF AND DCFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF DCFFTF OR DCFFTB.
C
      SUBROUTINE DCFFTI (N,WSAVE)
CIMPL DIMENSION       WSAVE(*)
      implicit none
      integer           n
      double precision  wsave(*)
C local variables
      integer           iw1
      integer           iw2
C
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL DCFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C=======================================================================
      SUBROUTINE DCFFTI1 (N,WA,IFAC)
CIMPL DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      implicit none
      integer           n
      double precision  wa(*)
      integer           ifac(*)
C local variables
      double precision  arg
      double precision  argh
      double precision  argld
      double precision  dum
      double precision  fi
      integer           i
      integer           i1
      integer           ib
      integer           ido
      integer           idot
      integer           ii
      integer           ip
      integer           ipm
      integer           j
      integer           k1
      integer           l1
      integer           l2
      integer           ld
      integer           nf
      integer           nl
      integer           nq
      integer           nr
      integer           ntry
      integer           ntryh(4)
      double precision  dpimach
      external          dpimach
      double precision  tpi
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 2.*DPIMACH(DUM)
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD+L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE DPASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
CIMPL DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
CIMPL1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
CIMPL2                CH2(IDL1,IP)
      implicit none
      integer           nac
      integer           ido
      integer           ip
      integer           l1
      integer           idl1
      double precision  cc(ido, ip, l1)
      double precision  c1(ido, l1, ip)
      double precision  c2(idl1, ip)
      double precision  ch(ido, l1, ip)
      double precision  ch2(idl1, ip)
      double precision  wa(*)
C local variables
      integer           i
      integer           idij
      integer           idj
      integer           idl
      integer           idlj
      integer           idot
      integer           idp
      integer           ik
      integer           inc
      integer           ipp2
      integer           ipph
      integer           j
      integer           jc
      integer           k
      integer           l
      integer           lc
      integer           nt
      double precision  wai
      double precision  war
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE DPASSB2 (IDO,L1,CC,CH,WA1)
CIMPL DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
CIMPL1                WA1(1)
      implicit none
      integer           ido
      integer           l1
      double precision  cc(ido, 2, l1)
      double precision  ch(ido, l1, 2)
      double precision  wa1(*)
C local variables
      integer           i
      integer           k
      double precision  ti2
      double precision  tr2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE DPASSB3 (IDO,L1,CC,CH,WA1,WA2)
CIMPL DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
CIMPL1                WA1(*)     ,WA2(*)
      implicit none
      integer           ido
      integer           l1
      double precision  cc(ido, 3, l1)
      double precision  ch(ido, l1, 3)
      double precision  wa1(*)
      double precision  wa2(*)
C local variables
      double precision  ci2
      double precision  ci3
      double precision  cr2
      double precision  cr3
      double precision  di2
      double precision  di3
      double precision  dr2
      double precision  dr3
      integer           i
      integer           k
      double precision  taui
      double precision  taur
      double precision  ti2
      double precision  tr2
      DATA TAUR,TAUI /-.5D0,.866025403784439D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE DPASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
CIMPL DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)
      implicit none
      integer           ido
      integer           l1
      double precision  cc(ido, 4, l1)
      double precision  ch(ido, l1, 4)
      double precision  wa1(*)
      double precision  wa2(*)
      double precision  wa3(*)
C local variables
      double precision  ci2
      double precision  ci3
      double precision  ci4
      double precision  cr2
      double precision  cr3
      double precision  cr4
      integer           i
      integer           k
      double precision  ti1
      double precision  ti2
      double precision  ti3
      double precision  ti4
      double precision  tr1
      double precision  tr2
      double precision  tr3
      double precision  tr4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE DPASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
CIMPL DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      implicit none
      integer           ido
      integer           l1
      double precision  cc(ido, 5, l1)
      double precision  ch(ido, l1, 5)
      double precision  wa1(*)
      double precision  wa2(*)
      double precision  wa3(*)
      double precision  wa4(*)
C local variables
      double precision  ci2
      double precision  ci3
      double precision  ci4
      double precision  ci5
      double precision  cr2
      double precision  cr3
      double precision  cr4
      double precision  cr5
      double precision  di2
      double precision  di3
      double precision  di4
      double precision  di5
      double precision  dr2
      double precision  dr3
      double precision  dr4
      double precision  dr5
      integer           i
      integer           k
      double precision  ti11
      double precision  ti12
      double precision  ti2
      double precision  ti3
      double precision  ti4
      double precision  ti5
      double precision  tr11
      double precision  tr12
      double precision  tr2
      double precision  tr3
      double precision  tr4
      double precision  tr5
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,
     1 .951056516295154D0,
     2-.809016994374947D0,.587785252292473D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      double precision  FUNCTION DPIMACH (DUM)
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      implicit none
      double precision  dum
C local variables
      double precision  four
      double precision  one
      parameter(four = 4.0)
      parameter(one = 1.0)
      DPIMACH = four*ATAN(one)
      RETURN
      END
C=======================================================================
C     SUBROUTINE SCFFTB(N,C,WSAVE)
C
C     SUBROUTINE SCFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
C     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , SCFFTB COMPUTES
C     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
C     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
C
C     A CALL OF SCFFTF FOLLOWED BY A CALL OF SCFFTB WILL MULTIPLY THE
C     SEQUENCE BY N.
C
C     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SCFFTB MUST BE
C     INITIALIZED BY CALLING SUBROUTINE SCFFTI(N,WSAVE).
C
C     INPUT PARAMETERS
C
C
C     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
C            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
C
C     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C
C     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
C             IN THE PROGRAM THAT CALLS SCFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE SCFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY SCFFTF AND SCFFTB.
C
C     OUTPUT PARAMETERS
C
C     C      FOR J=1,...,N
C
C                C(J)=THE SUM FROM K=1,...,N OF
C
C                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
C
C                            WHERE I=SQRT(-1)
C
C     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
C             DESTROYED BETWEEN CALLS OF SUBROUTINE SCFFTF OR SCFFTB
C
      SUBROUTINE SCFFTB (N,C,WSAVE)
CIMPL DIMENSION       C(*)       ,WSAVE(*)
      implicit none
      integer           n
      real              c(*)
      real              wsave(*)
C local variables
      integer           iw1
      integer           iw2
C
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL SCFFTB1 (N,C,WSAVE,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C=======================================================================
      SUBROUTINE SCFFTB1 (N,C,CH,WA,IFAC)
CIMPL DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      implicit none
      integer           n
      real              c(*)
      real              ch(*)
      real              wa(*)
      integer           ifac(*)
C local variables
      integer           i
      integer           idl1
      integer           ido
      integer           idot
      integer           ip
      integer           iw
      integer           ix2
      integer           ix3
      integer           ix4
      integer           k1
      integer           l1
      integer           l2
      integer           n2
      integer           na
      integer           nac
      integer           nf
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         CALL SPASSB4 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL SPASSB4 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL SPASSB2 (IDOT,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL SPASSB2 (IDOT,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         CALL SPASSB3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL SPASSB3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         CALL SPASSB5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL SPASSB5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL SPASSB (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL SPASSB (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      DO 117 I=1,N2
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C=======================================================================
C     SUBROUTINE SCFFTI(N,WSAVE)
C
C     SUBROUTINE SCFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH SCFFTF AND SCFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH SCFFTF AND SCFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF SCFFTF OR SCFFTB.
C
      SUBROUTINE SCFFTI (N,WSAVE)
CIMPL DIMENSION       WSAVE(*)
      implicit none
      integer           n
      real              wsave(*)
C local variables
      integer           iw1
      integer           iw2
C
      IF (N .EQ. 1) RETURN
      IW1 = N+N+1
      IW2 = IW1+N+N
      CALL SCFFTI1 (N,WSAVE(IW1),WSAVE(IW2))
      RETURN
      END
C=======================================================================
      SUBROUTINE SCFFTI1 (N,WA,IFAC)
CIMPL DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      implicit none
      integer           n
      real              wa(*)
      integer           ifac(*)
C local variables
      real              arg
      real              argh
      real              argld
      real              dum
      real              fi
      integer           i
      integer           i1
      integer           ib
      integer           ido
      integer           idot
      integer           ii
      integer           ip
      integer           ipm
      integer           j
      integer           k1
      integer           l1
      integer           l2
      integer           ld
      integer           nf
      integer           nl
      integer           nq
      integer           nr
      integer           ntry
      integer           ntryh(4)
      real              spimach
      external          spimach
      real              tpi
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/3,4,2,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 2.*SPIMACH(DUM)
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO 110 K1=1,NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD+L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO 108 II=4,IDOT,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IF (IP .LE. 5) GO TO 109
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SPASSB (NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
CIMPL DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
CIMPL1                C1(IDO,L1,IP)          ,WA(*)      ,C2(IDL1,IP),
CIMPL2                CH2(IDL1,IP)
      implicit none
      integer           nac
      integer           ido
      integer           ip
      integer           l1
      integer           idl1
      real              cc(ido, ip, l1)
      real              c1(ido, l1, ip)
      real              c2(idl1, ip)
      real              ch(ido, l1, ip)
      real              ch2(idl1, ip)
      real              wa(*)
C local variables
      integer           i
      integer           idij
      integer           idj
      integer           idl
      integer           idlj
      integer           idot
      integer           idp
      integer           ik
      integer           inc
      integer           ipp2
      integer           ipph
      integer           j
      integer           jc
      integer           k
      integer           l
      integer           lc
      integer           nt
      real              wai
      real              war
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO
C
      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SPASSB2 (IDO,L1,CC,CH,WA1)
CIMPL DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
CIMPL1                WA1(1)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 2, l1)
      real              ch(ido, l1, 2)
      real              wa1(*)
C local variables
      integer           i
      integer           k
      real              ti2
      real              tr2
      IF (IDO .GT. 2) GO TO 102
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(1,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(1,2,K)
         CH(2,K,1) = CC(2,1,K)+CC(2,2,K)
         CH(2,K,2) = CC(2,1,K)-CC(2,2,K)
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            CH(I-1,K,1) = CC(I-1,1,K)+CC(I-1,2,K)
            TR2 = CC(I-1,1,K)-CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
            TI2 = CC(I,1,K)-CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2+WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2-WA1(I)*TI2
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SPASSB3 (IDO,L1,CC,CH,WA1,WA2)
CIMPL DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
CIMPL1                WA1(*)     ,WA2(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 3, l1)
      real              ch(ido, l1, 3)
      real              wa1(*)
      real              wa2(*)
C local variables
      real              ci2
      real              ci3
      real              cr2
      real              cr3
      real              di2
      real              di3
      real              dr2
      real              dr3
      integer           i
      integer           k
      real              taui
      real              taur
      real              ti2
      real              tr2
      DATA TAUR,TAUI /-.5,.866025403784439D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TR2 = CC(1,2,K)+CC(1,3,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         TI2 = CC(2,2,K)+CC(2,3,K)
         CI2 = CC(2,1,K)+TAUR*TI2
         CH(2,K,1) = CC(2,1,K)+TI2
         CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
         CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
         CH(2,K,2) = CI2+CR3
         CH(2,K,3) = CI2-CR3
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TR2 = CC(I-1,2,K)+CC(I-1,3,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,2,K)+CC(I,3,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SPASSB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
CIMPL DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 4, l1)
      real              ch(ido, l1, 4)
      real              wa1(*)
      real              wa2(*)
      real              wa3(*)
C local variables
      real              ci2
      real              ci3
      real              ci4
      real              cr2
      real              cr3
      real              cr4
      integer           i
      integer           k
      real              ti1
      real              ti2
      real              ti3
      real              ti4
      real              tr1
      real              tr2
      real              tr3
      real              tr4
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,4,K)-CC(2,2,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,2,K)-CC(1,4,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,4,K)-CC(I,2,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,2,K)-CC(I-1,4,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SPASSB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
CIMPL DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 5, l1)
      real              ch(ido, l1, 5)
      real              wa1(*)
      real              wa2(*)
      real              wa3(*)
      real              wa4(*)
C local variables
      real              ci2
      real              ci3
      real              ci4
      real              ci5
      real              cr2
      real              cr3
      real              cr4
      real              cr5
      real              di2
      real              di3
      real              di4
      real              di5
      real              dr2
      real              dr3
      real              dr4
      real              dr5
      integer           i
      integer           k
      real              ti11
      real              ti12
      real              ti2
      real              ti3
      real              ti4
      real              ti5
      real              tr11
      real              tr12
      real              tr2
      real              tr3
      real              tr4
      real              tr5
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,
     1 .951056516295154D0,
     2-.809016994374947D0,.587785252292473D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END
C=======================================================================
      real              FUNCTION SPIMACH (DUM)
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      implicit none
      real              dum
C local variables
      real              four
      real              one
      parameter(four = 4.0)
      parameter(one = 1.0)
      SPIMACH = four*ATAN(one)
      RETURN
      END
C=======================================================================
      SUBROUTINE SRADB2 (IDO,L1,CC,CH,WA1)
CIMPL DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
CIMPL1                WA1(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 2, l1)
      real              ch(ido, l1, 2)
      real              wa1(*)
C local variables
      integer           i
      integer           ic
      integer           idp2
      integer           k
      real              ti2
      real              tr2
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END
C=======================================================================
      SUBROUTINE SRADB3 (IDO,L1,CC,CH,WA1,WA2)
CIMPL DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
CIMPL1                WA1(*)     ,WA2(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 3, l1)
      real              ch(ido, l1, 3)
      real              wa1(*)
      real              wa2(*)
C local variables
      real              ci2
      real              ci3
      real              cr2
      real              cr3
      real              di2
      real              di3
      real              dr2
      real              dr3
      integer           i
      integer           ic
      integer           idp2
      integer           k
      real              taui
      real              taur
      real              ti2
      real              tr2
      DATA TAUR,TAUI /-.5,.866025403784439D0/
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SRADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
CIMPL DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 4, l1)
      real              ch(ido, l1, 4)
      real              wa1(*)
      real              wa2(*)
      real              wa3(*)
C local variables
      real              ci2
      real              ci3
      real              ci4
      real              cr2
      real              cr3
      real              cr4
      integer           i
      integer           ic
      integer           idp2
      integer           k
      real              sqrt2
      real              ti1
      real              ti2
      real              ti3
      real              ti4
      real              tr1
      real              tr2
      real              tr3
      real              tr4
      DATA SQRT2 /1.414213562373095D0/
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END
C=======================================================================
      SUBROUTINE SRADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
CIMPL DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
CIMPL1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      implicit none
      integer           ido
      integer           l1
      real              cc(ido, 5, l1)
      real              ch(ido, l1, 5)
      real              wa1(*)
      real              wa2(*)
      real              wa3(*)
      real              wa4(*)
C local variables
      real              ci2
      real              ci3
      real              ci4
      real              ci5
      real              cr2
      real              cr3
      real              cr4
      real              cr5
      real              di2
      real              di3
      real              di4
      real              di5
      real              dr2
      real              dr3
      real              dr4
      real              dr5
      integer           i
      integer           ic
      integer           idp2
      integer           k
      real              ti11
      real              ti12
      real              ti2
      real              ti3
      real              ti4
      real              ti5
      real              tr11
      real              tr12
      real              tr2
      real              tr3
      real              tr4
      real              tr5
      DATA TR11,TI11,TR12,TI12 /.309016994374947D0,
     1 .951056516295154D0,
     2-.809016994374947D0,.587785252292473D0/
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SRADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
CIMPL DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
CIMPL1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
CIMPL2                CH2(IDL1,IP)           ,WA(*)
      implicit none
      integer           ido
      integer           ip
      integer           l1
      integer           idl1
      real              cc(ido, ip, l1)
      real              c1(ido, l1, ip)
      real              c2(idl1, ip)
      real              ch(ido, l1, ip)
      real              ch2(idl1, ip)
      real              wa(*)
C local variables
      real              ai1
      real              ai2
      real              ar1
      real              ar1h
      real              ar2
      real              ar2h
      real              arg
      real              dc2
      real              dcp
      real              ds2
      real              dsp
      real              dum
      integer           i
      integer           ic
      integer           idij
      integer           idp2
      integer           ik
      integer           ipp2
      integer           ipph
      integer           is
      integer           j
      integer           j2
      integer           jc
      integer           k
      integer           l
      integer           lc
      integer           nbd
      real              spimach
      external          spimach
      real              tpi
      TPI = 2.0*SPIMACH(DUM)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
C=======================================================================
C     SUBROUTINE SRFFTB(N,R,WSAVE)
C
C     SUBROUTINE SRFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS
C     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS DEFINED
C     BELOW AT OUTPUT PARAMETER R.
C
C     INPUT PARAMETERS
C
C     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
C             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
C             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
C
C     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
C             TO BE TRANSFORMED
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             IN THE PROGRAM THAT CALLS SRFFTB. THE WSAVE ARRAY MUST BE
C             INITIALIZED BY CALLING SUBROUTINE SRFFTI(N,WSAVE) AND A
C             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
C             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
C             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
C             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
C             THE SAME WSAVE ARRAY CAN BE USED BY SRFFTF AND SRFFTB.
C
C
C     OUTPUT PARAMETERS
C
C     R       FOR N EVEN AND FOR I = 1,...,N
C
C                  R(I) = R(1)+(-1)**(I-1)*R(N)
C
C                       PLUS THE SUM FROM K=2 TO K=N/2 OF
C
C                        2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                       -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C             FOR N ODD AND FOR I = 1,...,N
C
C                  R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF
C
C                       2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                      -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C      *****  NOTE
C                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF SRFFTF
C                  FOLLOWED BY A CALL OF SRFFTB WILL MULTIPLY THE INPUT
C                  SEQUENCE BY N.
C
C     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
C             CALLS OF SRFFTB OR SRFFTF.
C
C
      SUBROUTINE SRFFTB (N,R,WSAVE)
CIMPL DIMENSION       R(*)       ,WSAVE(*)
      implicit none
      integer           n
      real              r(*)
      real              wsave(*)
C
      IF (N .EQ. 1) RETURN
      CALL SRFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
C=======================================================================
      SUBROUTINE SRFFTB1 (N,C,CH,WA,IFAC)
CIMPL DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      implicit none
      integer           n
      real              c(*)
      real              ch(*)
      real              wa(*)
      integer           ifac(*)
C local variables
      integer           i
      integer           idl1
      integer           ido
      integer           ip
      integer           iw
      integer           ix2
      integer           ix3
      integer           ix4
      integer           k1
      integer           l1
      integer           l2
      integer           na
      integer           nf
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL SRADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL SRADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL SRADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL SRADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL SRADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL SRADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL SRADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL SRADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL SRADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL SRADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
C=======================================================================
C     SUBROUTINE SRFFTI(N,WSAVE)
C
C     SUBROUTINE SRFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
C     BOTH SRFFTF AND SRFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
C     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
C     STORED IN WSAVE.
C
C     INPUT PARAMETER
C
C     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
C
C     OUTPUT PARAMETER
C
C     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
C             THE SAME WORK ARRAY CAN BE USED FOR BOTH SRFFTF AND SRFFTB
C             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
C             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
C             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF SRFFTF OR SRFFTB.
C
      SUBROUTINE SRFFTI (N,WSAVE)
CIMPL DIMENSION       WSAVE(*)
      implicit none
      integer           n
      real              wsave(*)
C
      IF (N .EQ. 1) RETURN
      CALL SRFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
C=======================================================================
      SUBROUTINE SRFFTI1 (N,WA,IFAC)
CIMPL DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
      implicit none
      integer           n
      real              wa(*)
      integer           ifac(*)
C local variables
      real              arg
      real              argh
      real              argld
      real              dum
      real              fi
      integer           i
      integer           ib
      integer           ido
      integer           ii
      integer           ip
      integer           ipm
      integer           is
      integer           j
      integer           k1
      integer           l1
      integer           l2
      integer           ld
      integer           nf
      integer           nfm1
      integer           nl
      integer           nq
      integer           nr
      integer           ntry
      integer           ntryh(4)
      real              spimach
      external          spimach
      real              tpi
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 2.0*SPIMACH(DUM)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
C=======================================================================
C END OF FFTPACK
C=======================================================================
C
C=======================================================================
C
      SUBROUTINE FFT3C(MX, MY, NX, NY, NZ, A, ERROR)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER         MX, MY, NX, NY, NZ
      DOUBLE COMPLEX  A(MX, MY, NZ)
      LOGICAL         ERROR
C
C local
      INTEGER   MAXSEQ
      PARAMETER(MAXSEQ = 4096)
      DOUBLE COMPLEX  SEQ(MAXSEQ)
C
      INTEGER  IX, IY, IZ
      INTEGER  MWSPX, MWSPY, MWSPZ
      INTEGER   WSPX,  WSPY,  WSPZ
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'FFT3C: Using FFTPACK4.1'
C
      ERROR = .FALSE.
C
      IF (MAX(NY, NZ) .GT. MAXSEQ) THEN
        CALL WRNDIE(-5, 'FFT3C',
     &              'MaxSeq too small --> recompile program')
        ERROR = .TRUE.
        RETURN
      END IF
C
      MWSPX = 4 * NX + 15
       WSPX = ALLHP(IREAL8(MWSPX))
      MWSPY = 4 * NY + 15
       WSPY = ALLHP(IREAL8(MWSPY))
      MWSPZ = 4 * NZ + 15
       WSPZ = ALLHP(IREAL8(MWSPZ))
C
C Initialize FFT
      CALL DCFFTI(NX, HEAP(WSPX))
      CALL DCFFTI(NY, HEAP(WSPY))
      CALL DCFFTI(NZ, HEAP(WSPZ))
C
      DO IX = 1, NX
        DO IY = 1, NY
          DO IZ = 1, NZ
            SEQ(IZ) = A(IX, IY, IZ)
          END DO
C
C Transform along z (slow direction)
          CALL DCFFTB(NZ, SEQ, HEAP(WSPZ))
C
          DO IZ = 1, NZ
            A(IX, IY, IZ) = SEQ(IZ)
          END DO
        END DO
C
        DO IZ = 1, NZ
          DO IY = 1, NY
            SEQ(IY) = A(IX, IY, IZ)
          END DO
C
C Transform along y (medium direction)
          CALL DCFFTB(NY, SEQ, HEAP(WSPY))
C
          DO IY = 1, NY
            A(IX, IY, IZ) = SEQ(IY)
          END DO
        END DO
      END DO
C
      DO IZ = 1, NZ
        DO IY = 1, NY
C Transform along x (fast direction)
          CALL DCFFTB(NX, A(1, IY, IZ), HEAP(WSPX))
        END DO
      END DO
C
      CALL FREHP(WSPZ, IREAL8(MWSPZ))
      CALL FREHP(WSPY, IREAL8(MWSPY))
      CALL FREHP(WSPX, IREAL8(MWSPX))
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SFFT1C(N, A, ERROR, MWSP, WSP)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER  N
      COMPLEX           A(*)
      LOGICAL  ERROR
      INTEGER  MWSP
      INTEGER   WSP
C
C 1D complex-to-complex transform, single precision
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1C: Using FFTPACK4.1'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        IF (MWSP .EQ. 0) THEN
          MWSP = 4 * N + 15
           WSP = ALLHP(IREAL4(MWSP))
          CALL SCFFTI(N, HEAP(WSP))
        END IF
        CALL SCFFTB(N, A, HEAP(WSP))
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, IREAL4(MWSP))
        MWSP = 0
         WSP = 0
      END IF
C
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE SFFT1CR(N, A, ERROR, MWSP, WSP)
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INTEGER  N
      REAL              A(0:*)
      LOGICAL  ERROR
      INTEGER  MWSP
      INTEGER   WSP
C
C 1D complex-to-real transform, single precision
C
C local
      INTEGER  I
C
C begin
      IF (WRNLEV .GT. 10)
     &  WRITE(6, '(1X, A)') 'SFFT1CR: Using FFTPACK4.1'
C
      ERROR = .FALSE.
C
      IF (N .GT. 0) THEN
        DO I = 0, N - 1
          A(I + 1) = A(I + 2)
        END DO
        IF (MWSP .EQ. 0) THEN
          MWSP = 2 * N + 15
           WSP = ALLHP(IREAL4(MWSP))
          CALL SRFFTI(N, HEAP(WSP))
        END IF
        CALL SRFFTB(N, A, HEAP(WSP))
      ELSE IF (MWSP .GT. 0) THEN
        CALL FREHP(WSP, IREAL4(MWSP))
        MWSP = 0
         WSP = 0
      END IF
C
      RETURN
      END
