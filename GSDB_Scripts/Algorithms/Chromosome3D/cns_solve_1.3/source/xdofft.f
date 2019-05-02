      SUBROUTINE XDOFT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & RHOMASK,RRHO,IRHO,
     & XRNREF,XRH,XRK,XRL,FCALC,NAP,NBP,NCP,NAPP,NBPP,NCPP,
     & XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
C Computes Fourier Transform of structure factors.
C Front-end routine for XDOFT2 and XDOFT3.
C
C Input is an asymmetric unit of the structure factor
C stored in FCALC, for REFLCT=1,...,XRNREF (all non-zero FCALC
C elements == selected elements).
C Output is an asymmetric unit of the map, stored in the
C brick A=MAASY,...,NAASY, B=MBASY,...,NBASY, C=MCASY,...,NCASY for
C those elements with RHOMASK(A,B,C).NE.0.
C
C I/O
      IMPLICIT NONE
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NA, NB, NC, MAASY,MBASY,MCASY, NAASY, NBASY, NCASY
      INTEGER  RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      INTEGER NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRVOL, XRCELL(9)
C local
      INTEGER NRHOL
C pointers
      INTEGER CWK, RHOL, ITEMP, CTEMP
C begin
C
C
C check grid to make sure that it is fine enough
      CALL XFGRCHK(NA,NB,NC,XRNREF,.FALSE.,1,XRH,XRK,XRL,FCALC)
C
C call memory-saving version (somewhat more CPU-intensive)
      IF (NA.NE.NAP.OR.NB.NE.NBP.OR.NC.NE.NCP) THEN
      CWK=ALLHP(ICPLX8(NA+2))
      ITEMP=ALLHP(INTEG4(NAASY-MAASY+1))
      CTEMP=ALLHP(ICPLX8(NAASY-MAASY+1))
      NRHOL=NAP*NBP*NCP
      RHOL=ALLHP(ICPLX8(NRHOL))
      CALL XDOFT3(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     RHOMASK,RRHO,IRHO,
     &     XRNREF,XRH,XRK,XRL,FCALC,HEAP(CWK),
     &     HEAP(ITEMP),HEAP(CTEMP),
     &     NAP,NBP,NCP,NAPP,NBPP,NCPP,HEAP(RHOL),
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
      CALL FREHP(ITEMP,INTEG4(NAASY-MAASY+1))
      CALL FREHP(CTEMP,ICPLX8(NAASY-MAASY+1))
      CALL FREHP(RHOL,ICPLX8(NRHOL))
      CALL FREHP(CWK,ICPLX8(NA+2))
C
      ELSE
C
C call fast version without sublattices
      CWK=ALLHP(ICPLX8(NA+2))
      IF (QHERM) THEN
      NRHOL=NA*NB*(NC/2+1)
      ELSE
      NRHOL=NA*NB*NC
      END IF
      RHOL=ALLHP(ICPLX8(NRHOL))
      CALL XDOFT2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     RHOMASK,RRHO,IRHO,
     &     XRNREF,XRH,XRK,XRL,FCALC,HEAP(CWK),
     &     HEAP(RHOL),
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
      CALL FREHP(RHOL,ICPLX8(NRHOL))
      CALL FREHP(CWK,ICPLX8(NA+2))
C
      END IF
      RETURN
      END
C=====================================================================
      SUBROUTINE XDOFT2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & RHOMASK,RRHO,IRHO,
     & XRNREF,XRH,XRK,XRL,FCALC,CWK,RHOL,
     & XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
C Structure factor to map.  Fast version
C without sublattices.
C
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER NA, NB, NC, MAASY,MBASY,MCASY, NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNREF,XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      DOUBLE COMPLEX CWK(*)
      DOUBLE COMPLEX RHOL(NA,NB,*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRVOL, XRCELL(9)
C local
      INTEGER H, K, L, HH, KK, LL, NCC, ISYM, REFLCT
      INTEGER HM, KM, LM, LR
      LOGICAL ERROR
      DOUBLE PRECISION RTH, TEMP, SCALE
      DOUBLE PRECISION SEC1, SEC2, SEC3, SECS, SECS2
      DOUBLE COMPLEX FLOC
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      SCALE=ONE/XRVOL
C
C initialize timers
      IF (TIMER.GT.0) THEN
      SEC1=ZERO
      SEC2=ZERO
      SEC3=ZERO
      END IF
C
      IF (QHERM) THEN
      NCC=NC/2+1
      ELSE
      NCC=NC
      END IF
C
C step1: expand structure factors and copy into local map
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C initialize RHO
      DO L=1,NCC
      DO K=1,NB
      DO H=1,NA
      RHOL(H,K,L)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
C
      RTH=TWO*PI/XRSYTH
C
C loop over all symmetry operators
      DO ISYM=1,XRNSYM
C
C in the following loop we compute the complex factors and H,K,L indices
      DO REFLCT=1,XRNREF
C
C only use non-zero elements of FCALC.
      IF (ABS(FCALC(REFLCT)).GT.RSMALL) THEN
C
C apply transpose of inversed symmetry operator to H, K, L
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      H=XRITSY(ISYM,1,1)*HH+XRITSY(ISYM,1,2)*KK+XRITSY(ISYM,1,3)*LL
      K=XRITSY(ISYM,2,1)*HH+XRITSY(ISYM,2,2)*KK+XRITSY(ISYM,2,3)*LL
      L=XRITSY(ISYM,3,1)*HH+XRITSY(ISYM,3,2)*KK+XRITSY(ISYM,3,3)*LL
C
C compute phase shift of translational part of symmetry operator,
C use TRANSFORMED indices
      TEMP=( XRSYMM(ISYM,1,4)*H
     &      +XRSYMM(ISYM,2,4)*K
     &      +XRSYMM(ISYM,3,4)*L ) *RTH
C
C multiply map by scale in order to have it on an absolute scale
      FLOC=DCMPLX(DCOS(TEMP),DSIN(TEMP))*(FCALC(REFLCT))*SCALE
C
C project H,K,L into periodic lattices (1...NA,...)
      H=MOD(H+ 10000*NA,NA) +1
      K=MOD(K+ 10000*NB,NB) +1
      L=MOD(L+ 10000*NC,NC) +1
C
C if QHERM, map H,K,L into L<=NC/2+1 hemisphere (real RHO symmetry)
      HM=NA+2-H
      IF (H.EQ.1) HM=1
      KM=NB+2-K
      IF (K.EQ.1) KM=1
      LM=NC+2-L
      IF (L.EQ.1) LM=1
      IF (QHERM.AND.(L.GT.NC/2+1
     &     .OR.(L.EQ.LM.AND.K.GT.NB/2+1)
     &     .OR.(L.EQ.LM.AND.K.EQ.KM.AND.H.GT.NA/2+1))) THEN
C take conjugate of FLOC for inverse FFT.
      RHOL(HM,KM,LM)=FLOC
      ELSE
C take conjugate of FLOC for inverse FFT.
      RHOL(H,K,L)=DCONJG(FLOC)
      END IF
C
      END IF
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC1=SEC1+SECS2-SECS
      END IF
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C step2: FFT local map
C ++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
      IF (QHERM) THEN
      CALL FFT3CR2(NA,NB,NA,NB,NC,RHOL,CWK,ERROR)
      ELSE
      CALL FFT3C(NA,NB,NA,NB,NC,RHOL,ERROR)
      END IF
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC2=SEC2+SECS2-SECS
      END IF
C
C step3: copy local map into asymmetric unit of main map
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
C take conjugate for inverse FFT
      IF (.NOT.QHERM) THEN
C
      DO L=MCASY,NCASY
      LL=MOD(L+10000*NC,NC)+1
      DO K=MBASY,NBASY
      KK=MOD(K+10000*NB,NB)+1
      DO H=MAASY,NAASY
      HH=MOD(H+10000*NA,NA)+1
      IF (RHOMASK(H,K,L).NE.0) THEN
      RRHO(H,K,L)=DBLE(RHOL(HH,KK,LL))
      IRHO(H,K,L)=-DIMAG(RHOL(HH,KK,LL))
      END IF
      END DO
      END DO
      END DO
C
      ELSE
      DO L=MCASY,NCASY
      LL=MOD(L+10000*NC,NC)+1
      LR=(LL+1)/2
      LL=MOD(LL,2)
      DO K=MBASY,NBASY
      KK=MOD(K+10000*NB,NB)+1
      DO H=MAASY,NAASY
      HH=MOD(H+10000*NA,NA)+1
      IF (RHOMASK(H,K,L).NE.0) THEN
      IF (LL.EQ.1) THEN
      RRHO(H,K,L)=DBLE(RHOL(HH,KK,LR))
      ELSE
      RRHO(H,K,L)=DIMAG(RHOL(HH,KK,LR))
      END IF
      END IF
      END DO
      END DO
      END DO
C
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC3=SEC3+SECS2-SECS
      END IF
C
      IF (TIMER.GT.0) THEN
      WRITE(6,'(A,F10.4,/,A,F10.4,/,A,F10.4)')
     & ' XDOFT2: CPU-time:  transformation1=',SEC1,
     & ' XDOFT2: CPU-time:  FFT=',SEC2,
     & ' XDOFT2: CPU-time:  transformation2=',SEC3
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOFT3(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & RHOMASK,RRHO,IRHO,
     & XRNREF,XRH,XRK,XRL,FCALC,
     & CWK,ITEMP,CTEMP,NAP,NBP,NCP,NAPP,NBPP,NCPP,RHOL,
     & XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
C Structure factor to map.  Sublattice version.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'timer.inc'
      INCLUDE 'consta.inc'
      INTEGER NA, NB, NC, NAASY, NBASY, NCASY, MAASY, MBASY, MCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      DOUBLE COMPLEX CWK(*)
      INTEGER ITEMP(MAASY:NAASY)
      DOUBLE COMPLEX CTEMP(MAASY:NAASY)
      INTEGER NAP, NBP, NCP, NAPP, NBPP, NCPP
      DOUBLE COMPLEX RHOL(NAP,NBP,*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRVOL, XRCELL(9)
C local
      INTEGER H, K, L, HH, KK, LL, ISYM, REFLCT
      INTEGER NAPPP, NBPPP, NCPPP
      INTEGER AL, BL, CL, A, B, C, HM, KM, LM
      LOGICAL ERROR
      DOUBLE PRECISION RTH, TEMP, RNA, RNB, RNC, SCALE
      DOUBLE PRECISION SEC1, SEC2, SEC3, SECS, SECS2
      DOUBLE COMPLEX FLOC
      DOUBLE COMPLEX PHASEA, PHASEB, PHASEC
      DOUBLE COMPLEX PHASEA1, PHASEB1, PHASEC1
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C begin
C
      SCALE=ONE/XRVOL
C
C initialize timers
      IF (TIMER.GT.0) THEN
      SEC1=ZERO
      SEC2=ZERO
      SEC3=ZERO
      END IF
C
C initialize map
      IF (QHERM) THEN
      DO L=MCASY,NCASY
      DO K=MBASY,NBASY
      DO H=MAASY,NAASY
      IF (RHOMASK(H,K,L).NE.0) THEN
      RRHO(H,K,L)=ZERO
      END IF
      END DO
      END DO
      END DO
      ELSE
      DO L=MCASY,NCASY
      DO K=MBASY,NBASY
      DO H=MAASY,NAASY
      IF (RHOMASK(H,K,L).NE.0) THEN
      RRHO(H,K,L)=ZERO
      IRHO(H,K,L)=ZERO
      END IF
      END DO
      END DO
      END DO
      END IF
C
      RTH=TWO*PI/XRSYTH
C
C
C loop over all sublattices
      DO NAPPP=0,NAPP-1
      DO NBPPP=0,NBPP-1
      DO NCPPP=0,NCPP-1
C
C step 1: expand structure factors to P1 and copy into sublattice map
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C initialize sublattice map
      DO L=1,NCP
      DO K=1,NBP
      DO H=1,NAP
      RHOL(H,K,L)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
C
C loop over all symmetry operators
      DO ISYM=1,XRNSYM
C
C in the following loop we compute the complex factors and
C H,K,L indices
      DO REFLCT=1,XRNREF
C
C only use non-zero elements of FCALC
      IF (ABS(FCALC(REFLCT)).GT.RSMALL) THEN
C apply transpose of inversed symmetry operator to H, K, L
      H=XRH(REFLCT)
      K=XRK(REFLCT)
      L=XRL(REFLCT)
      HH=XRITSY(ISYM,1,1)*H+XRITSY(ISYM,1,2)*K+XRITSY(ISYM,1,3)*L
      KK=XRITSY(ISYM,2,1)*H+XRITSY(ISYM,2,2)*K+XRITSY(ISYM,2,3)*L
      LL=XRITSY(ISYM,3,1)*H+XRITSY(ISYM,3,2)*K+XRITSY(ISYM,3,3)*L
C
C project H,K,L into periodic lattices (1...NA,...)
      H=MOD(HH+ 10000*NA,NA) +1
      K=MOD(KK+ 10000*NB,NB) +1
      L=MOD(LL+ 10000*NC,NC) +1
C
      HM=NA+2-H
      IF (H.EQ.1) HM=1
      KM=NB+2-K
      IF (K.EQ.1) KM=1
      LM=NC+2-L
      IF (L.EQ.1) LM=1
C
C use sublattice here
      IF (MOD(H-NAPPP-1,NAPP).EQ.0.AND.
     &    MOD(K-NBPPP-1,NBPP).EQ.0.AND.
     &    MOD(L-NCPPP-1,NCPP).EQ.0) THEN
C
C take conjugate for inverse FFT.
C multiply map by scale in order to have it on an absolute scale.
C compute phase shift of translational part of symmetry operator,
C use TRANSFORMED indices
      TEMP=( XRSYMM(ISYM,1,4)*HH
     &      +XRSYMM(ISYM,2,4)*KK
     &      +XRSYMM(ISYM,3,4)*LL ) *RTH
      FLOC=DCONJG(DCMPLX(DCOS(TEMP),DSIN(TEMP))*FCALC(REFLCT))*SCALE
C
C use sublattice here
      H=(H-NAPPP-1)/NAPP+1
      K=(K-NBPPP-1)/NBPP+1
      L=(L-NCPPP-1)/NCPP+1
      RHOL(H,K,L)=FLOC
      END IF
C
C use sublattice here
      IF (QHERM) THEN
      IF (MOD(HM-NAPPP-1,NAPP).EQ.0.AND.
     &    MOD(KM-NBPPP-1,NBPP).EQ.0.AND.
     &    MOD(LM-NCPPP-1,NCPP).EQ.0) THEN
C
C take conjugate for inverse FFT.
C multiply map by scale in order to have it on an absolute scale.
C compute phase shift of translational part of symmetry operator,
C use TRANSFORMED indices
      TEMP=( XRSYMM(ISYM,1,4)*HH
     &      +XRSYMM(ISYM,2,4)*KK
     &      +XRSYMM(ISYM,3,4)*LL ) *RTH
      FLOC=DCONJG(DCMPLX(DCOS(TEMP),DSIN(TEMP))*FCALC(REFLCT))*SCALE
C
C use sublattice here
      HM=(HM-NAPPP-1)/NAPP+1
      KM=(KM-NBPPP-1)/NBPP+1
      LM=(LM-NCPPP-1)/NCPP+1
      RHOL(HM,KM,LM)=DCONJG(FLOC)
      END IF
      END IF
C
      END IF
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC1=SEC1+SECS2-SECS
      END IF
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C step2: FFT of sublattice map
C ++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL FFT3C(NAP,NBP,NAP,NBP,NCP,RHOL,ERROR)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC2=SEC2+SECS2-SECS
      END IF
C ++++++++++++++++++++++++++++
C
C step3: transform sublattice map and copy into asymmetric unit
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
      RNA=TWO*PI/NA
      RNB=TWO*PI/NB
      RNC=TWO*PI/NC
      RNA=RNA*NAPPP
      RNB=RNB*NBPPP
      RNC=RNC*NCPPP
C
C precompute stuff for the inner-most loop
      PHASEA1=DCMPLX(DCOS(RNA),DSIN(RNA))
      PHASEA=PHASEA1**MAASY
      DO A=MAASY,NAASY
      ITEMP(A)=MOD(A+10000*NAP,NAP)+1
      CTEMP(A)=PHASEA
      PHASEA=PHASEA*PHASEA1
      END DO
C
      PHASEC1=DCMPLX(DCOS(RNC),DSIN(RNC))
      PHASEC=DCMPLX(ONE,ZERO)*(PHASEC1**MCASY)
C
      DO C=MCASY,NCASY
      PHASEB1=DCMPLX(DCOS(RNB),DSIN(RNB))
      PHASEB=PHASEC*(PHASEB1**MBASY)
      CL=MOD(C+10000*NCP,NCP)+1
C
      DO B=MBASY,NBASY
      BL=MOD(B+10000*NBP,NBP)+1
C
      IF (QHERM) THEN
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      AL=ITEMP(A)
C
C use only the real part
      RRHO(A,B,C)=RRHO(A,B,C)+
     &  DBLE(PHASEB*CTEMP(A)*RHOL(AL,BL,CL))
      END IF
      END DO
      ELSE
      DO A=MAASY,NAASY
      IF (RHOMASK(A,B,C).NE.0) THEN
      AL=ITEMP(A)
      FLOC=PHASEB*CTEMP(A)*RHOL(AL,BL,CL)
C
C take the complex conjugate (inverse FFT!)
      RRHO(A,B,C)=RRHO(A,B,C)+DBLE(FLOC)
      IRHO(A,B,C)=IRHO(A,B,C)-DIMAG(FLOC)
      END IF
      END DO
      END IF
C
      PHASEB=PHASEB*PHASEB1
      END DO
      PHASEC=PHASEC*PHASEC1
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC3=SEC3+SECS2-SECS
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      END DO
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      WRITE(6,'(A,F10.4,/,A,F10.4,/,A,F10.4)')
     & ' XDOFT: CPU-time:  transformation1=',SEC1,
     & ' XDOFT: CPU-time:  FFT=',SEC2,
     & ' XDOFT: CPU-time:  transformation2=',SEC3
      END IF
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOIFT(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & RHOMASK,RRHO,IRHO,
     & XRNREF,XRH,XRK,XRL,FCALC,
     & NAP,NBP,NCP,NAPP,NBPP,NCPP,
     & XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
C Map to structure factor.  Front-end routine for XDOIFT2.
C
C Input is an asymmetric unit of the map, stored in the
C brick A=MAASY,...,NAASY, B=MBASY,...,NBASY, C=MCASY,...,NCASY for
C those elements with RHOMASK(A,B,C).NE.0.
C
C Output is an asymmetric unit of the structure factor
C stored in FCALC, for REFLCT=1,...,XRNREF.
C
C Author: Axel T. Brunger
C =======================
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER NA, NB, NC, MAASY,MBASY,MCASY, NAASY, NBASY, NCASY
      INTEGER RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      INTEGER NAP, NBP, NCP, NAPP, NBPP, NCPP
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRVOL, XRCELL(9)
C local
      INTEGER NRHOL
C pointers
      INTEGER CWK, RHOL, ITEMP, JTEMP
C begin
      CWK=ALLHP(ICPLX8(NA+2))
      IF (QHERM) THEN
      NRHOL=NAP*NBP*(NCP/2+1)
      ELSE
      NRHOL=NAP*NBP*NCP
      END IF
      RHOL=ALLHP(ICPLX8(NRHOL))
      ITEMP=ALLHP(INTEG4(NAP))
      JTEMP=ALLHP(INTEG4(NAP))
      CALL XDOIFT2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     &     RHOMASK,RRHO,IRHO,
     &     XRNREF,XRH,XRK,XRL,FCALC,HEAP(CWK),
     &     HEAP(ITEMP),HEAP(JTEMP),
     &     NAP,NBP,NCP,NAPP,NBPP,NCPP,HEAP(RHOL),
     &     XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
      CALL FREHP(JTEMP,INTEG4(NAP))
      CALL FREHP(ITEMP,INTEG4(NAP))
      CALL FREHP(RHOL,ICPLX8(NRHOL))
      CALL FREHP(CWK,ICPLX8(NA+2))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOIFT2(NA,NB,NC,MAASY,MBASY,MCASY,NAASY,NBASY,NCASY,
     & RHOMASK,RRHO,IRHO,
     & XRNREF,XRH,XRK,XRL,FCALC,
     & CWK,ITEMP,JTEMP,
     & NAP,NBP,NCP,NAPP,NBPP,NCPP,RHOL,
     & XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,XRVOL,XRCELL)
C
C Map to structure factor.
C
C Author: Axel T. Brunger
C =======================
C
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER NA, NB, NC, MAASY, MBASY, MCASY, NAASY, NBASY, NCASY
      INTEGER  RHOMASK(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL RRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      REAL IRHO(MAASY:NAASY,MBASY:NBASY,MCASY:NCASY)
      INTEGER XRNREF, XRH(*), XRK(*), XRL(*)
      DOUBLE COMPLEX FCALC(*)
      DOUBLE COMPLEX CWK(*)
      INTEGER ITEMP(*), JTEMP(*)
      INTEGER NAP, NBP, NCP, NAPP, NBPP, NCPP
      DOUBLE COMPLEX RHOL(NAP,NBP,*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      DOUBLE PRECISION XRVOL, XRCELL(9)
C local
      INTEGER H, K, L, ISYM, REFLCT, HH, KK, LL, HM, KM, LM
      INTEGER NAPPP, NBPPP, NCPPP
      INTEGER A, B, C, AL, BL, CL, NCCP
      INTEGER CRL, CCL
      LOGICAL ERROR
      DOUBLE PRECISION RTH, GRDVOL, TEMP, RNA, RNB, RNC
      DOUBLE PRECISION SEC1, SEC2, SEC3, SECS, SECS2
      DOUBLE COMPLEX RLOC
      INTEGER NAPCCC, ALL
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, FOUR=4.0D0)
C begin
C
C
C initialize timers
      IF (TIMER.GT.0) THEN
      SEC1=ZERO
      SEC2=ZERO
      SEC3=ZERO
      END IF
C
      IF (QHERM) THEN
      NCCP=NCP/2+1
      ELSE
      NCCP=NCP
      END IF
C
C initialize FCALC
      DO REFLCT=1,XRNREF
      FCALC(REFLCT)=DCMPLX(ZERO,ZERO)
      END DO
C
      GRDVOL=NA*NB*NC/XRVOL
C
      RTH=TWO*PI/XRSYTH
C
C loop over all sublattices
      DO NAPPP=0,NAPP-1
      DO NBPPP=0,NBPP-1
      DO NCPPP=0,NCPP-1
C
C step1: project ASU of map into sublattice map
C +++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
      RNA=TWO*PI/NA
      RNB=TWO*PI/NB
      RNC=TWO*PI/NC
      RNA=RNA*NAPPP
      RNB=RNB*NBPPP
      RNC=RNC*NCPPP
C
      DO CL=1,NCCP
      DO BL=1,NBP
      DO AL=1,NAP
      RHOL(AL,BL,CL)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
C
C precompute things for the innermost loop
      NAPCCC=0
      DO AL=1,NAP
      A=(AL-1)*NAPP+NAPPP
      A=MOD(A+ 10000*NA -MAASY,NA) +MAASY
      IF (A.LE.NAASY.AND.A.GE.MAASY) THEN
      NAPCCC=NAPCCC+1
      ITEMP(NAPCCC)=AL
      JTEMP(NAPCCC)=A
      END IF
      END DO
C
C*$* OPTIMIZE(0)
      DO CL=1,NCP
      C=(CL-1)*NCPP+NCPPP
      C=MOD(C+ 10000*NC -MCASY,NC) +MCASY
      IF (C.LE.NCASY.AND.C.GE.MCASY) THEN
      IF (QHERM) THEN
      CRL=(CL+1)/2
      CCL=MOD(CL,2)
      ELSE
      CRL=CL
      END IF
C
      DO BL=1,NBP
      B=(BL-1)*NBPP+NBPPP
      B=MOD(B+ 10000*NB -MBASY,NB) +MBASY
      IF (B.LE.NBASY.AND.B.GE.MBASY) THEN
C
      IF (QHERM.AND.CCL.EQ.1) THEN
      DO ALL=1,NAPCCC
      AL=ITEMP(ALL)
      A=JTEMP(ALL)
      IF (RHOMASK(A,B,C).NE.0) THEN
      TEMP=RRHO(A,B,C)/RHOMASK(A,B,C)
      RHOL(AL,BL,CRL)=DCMPLX(TEMP,DIMAG(RHOL(AL,BL,CRL)))
      END IF
      END DO
C
      ELSEIF (QHERM.AND.CCL.EQ.0) THEN
      DO ALL=1,NAPCCC
      AL=ITEMP(ALL)
      A=JTEMP(ALL)
      IF (RHOMASK(A,B,C).NE.0) THEN
      TEMP=RRHO(A,B,C)/RHOMASK(A,B,C)
      RHOL(AL,BL,CRL)=DCMPLX(DBLE(RHOL(AL,BL,CRL)),TEMP)
      END IF
      END DO
C
      ELSE
C
      DO ALL=1,NAPCCC
      AL=ITEMP(ALL)
      A=JTEMP(ALL)
      IF (RHOMASK(A,B,C).NE.0) THEN
      RHOL(AL,BL,CL)=DCMPLX(RRHO(A,B,C),IRHO(A,B,C))/RHOMASK(A,B,C)
      END IF
      END DO
C
      END IF
      END IF
      END DO
C
      END IF
      END DO
C*$* OPTIMIZE(5)
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC1=SEC1+SECS2-SECS
      END IF
C ++++++++++++++++++++++++++++++++++++++
C
C step2: take FFT of sublattice map
C +++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
      IF (QHERM) THEN
      CALL FFT3RC2(NAP,NBP,NAP,NBP,NCP,RHOL,CWK,ERROR)
      ELSE
      CALL FFT3C(NAP,NBP,NAP,NBP,NCP,RHOL,ERROR)
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC2=SEC2+SECS2-SECS
      END IF
C +++++++++++++++++++++++++++++++++
C
C step3: reduce sublattice map and copy into structure factors
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (TIMER.GT.0) CALL VCPU(SECS)
C
C loop over all symmetry operators
      DO ISYM=1,XRNSYM
C
C in the following loop we compute the complex factors
C and H,K,L indices
      DO REFLCT=1,XRNREF
C
C apply transpose of symmetry operator to H, K, L
      HH=XRH(REFLCT)
      KK=XRK(REFLCT)
      LL=XRL(REFLCT)
      H=XRSYMM(ISYM,1,1)*HH + XRSYMM(ISYM,2,1)*KK + XRSYMM(ISYM,3,1)*LL
      K=XRSYMM(ISYM,1,2)*HH + XRSYMM(ISYM,2,2)*KK + XRSYMM(ISYM,3,2)*LL
      L=XRSYMM(ISYM,1,3)*HH + XRSYMM(ISYM,2,3)*KK + XRSYMM(ISYM,3,3)*LL
C
C compute phase shift of translational part of symmetry operator
C and apply sublattice phase shift
      TEMP= (XRSYMM(ISYM,1,4)*HH
     &     + XRSYMM(ISYM,2,4)*KK
     &     + XRSYMM(ISYM,3,4)*LL ) *RTH
     &     + H*RNA + K*RNB + L*RNC
C
C project H,K,L into periodic sublattices (1...NAP,...)
      H=MOD(H+ 10000*NAP,NAP) +1
      K=MOD(K+ 10000*NBP,NBP) +1
      L=MOD(L+ 10000*NCP,NCP) +1
C
C if QHERM, map H,K,L into L<=NC/2+1 hemisphere (real RHO symmetry)
      HM=NAP+2-H
      IF (H.EQ.1) HM=1
      KM=NBP+2-K
      IF (K.EQ.1) KM=1
      LM=NCP+2-L
      IF (L.EQ.1) LM=1
      IF (QHERM.AND.(L.GT.NCP/2+1
     &     .OR.(L.EQ.LM.AND.K.GT.NBP/2+1)
     &     .OR.(L.EQ.LM.AND.K.EQ.KM.AND.H.GT.NAP/2+1))) THEN
      RLOC=DCONJG(RHOL(HM,KM,LM))
      ELSE
      RLOC=RHOL(H,K,L)
      END IF
C
      FCALC(REFLCT)=FCALC(REFLCT)+
     &              DCMPLX(DCOS(TEMP),DSIN(TEMP))*RLOC/GRDVOL
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      SEC3=SEC3+SECS2-SECS
      END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      END DO
      END DO
      END DO
C
      IF (TIMER.GT.0) THEN
      WRITE(6,'(A,F10.4,/,A,F10.4,/,A,F10.4)')
     & ' XDOIFT: CPU-time:  transformation1=',SEC1,
     & ' XDOIFT: CPU-time:  FFT=',SEC2,
     & ' XDOIFT: CPU-time:  transformation2=',SEC3
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE FFT3CR2(MX,MY,NX,NY,NZ,A,CWK,ERROR)
C
C SAVE!!!!!! and put in main version (dervatives and map calc!!!)
C complex array version.
C
C Computes 3-dimensional fast Fourier transform
C
C    X(x+1,y+1,z+1) = SUM {h=0...NX-1,k=0...NY-1,l=0...NZ-1}
C          A(h+1,k+1,l+1) *EXP(2*PI*i*(x*h/NX+y*k/NY+z*l/NZ))
C
C of a hermitian symmetric complex valued sequence A(x,y,z),
C {x=1,...,NX, y=1,...,NY, z=1,...,NZ/2+1}, one hemisphere,
C the hermitian symmetry implies A(x,y,NZ+2-z) = conjg[A(x,y,z)])
C
C The result is a real valued sequence a(h,k,l), which is stored
C in the complex array A as
C
C    A(h,k,l) = a(h,k,{2l-1}) + i*a(h,k,{2l})
C
C where h=1,...,NX, k=1,...,NY, l=1,...,NZ/2.
C
C Special case Z=-Z (i.e., z=NZ+2-z):
C    A(x,y,z) is expected to be defined for x=1,...,NX, y=1...NY/2+1
C Special case Z=-Z and Y=-Y:
C    A(x,y,z) is expected to be defined for x=1,...,NX/2+1
C
C CWK is a complex work array of dimension NX
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER MX, MY, NX, NY, NZ
      DOUBLE COMPLEX A(MX,MY,NZ/2+1)
      DOUBLE COMPLEX CWK(NX)
      LOGICAL ERROR
C local
      INTEGER X, Y, Z, XM, YM, ZM, NXX, NYY
      DOUBLE COMPLEX PHASE1, PHASE, ODD, EVEN
C parameters
      DOUBLE PRECISION HALF, ONE, ZERO, TPI
      PARAMETER (HALF=0.5D0, ONE=1.0D0, ZERO=0.0D0, TPI=2.0D0*PI)
C begin
      ERROR=.FALSE.
C
C the data in RHO have to be symmetrized for Z=-Z.
      DO Z=1,NZ/2+1
      ZM=MOD(NZ+1-Z,NZ)+1
      IF (Z.EQ.ZM) THEN
      DO Y=1,NY/2+1
      YM=MOD(NY+1-Y,NY)+1
      IF (Y.EQ.YM) THEN
      NXX=NX/2+1
      ELSE
      NXX=NX
      END IF
C
      DO X=2,NXX
      XM=NX-X+2
      A(XM,YM,ZM)=DCONJG(A(X,Y,Z))
      END DO
      A(1,YM,ZM)=DCONJG(A(1,Y,Z))
C
      END DO
      END IF
      END DO
C
      PHASE1=DCMPLX(DCOS(TPI/NZ),DSIN(TPI/NZ))
      PHASE=DCMPLX(ONE,ZERO)
      DO Z=1,NZ/4+1
      ZM=NZ/2+2-Z
      IF (ZM.EQ.Z) THEN
      NYY=NY/2+1
      ELSE
      NYY=NY
      END IF
      DO Y=1,NYY
      YM=MOD(NY+1-Y,NY)+1
      IF (YM.EQ.Y.AND.ZM.EQ.Z) THEN
      NXX=NX/2+1
      ELSE
      NXX=NX
      END IF
      CWK(1)=A(1,YM,ZM)
      DO X=2,NXX
      CWK(X)=A(NX-X+2,YM,ZM)
      END DO
      DO X=1,NXX
      EVEN=A(X,Y,Z)+DCONJG(CWK(X))
      ODD=DCONJG(A(X,Y,Z))-CWK(X)
      ODD=DCMPLX(DIMAG(ODD),DBLE(ODD))*PHASE
      A(X,Y,Z)=EVEN+ODD
      CWK(X)=DCONJG(EVEN-ODD)
      END DO
      A(1,YM,ZM)=CWK(1)
      DO X=2,NXX
      A(NX-X+2,YM,ZM)=CWK(X)
      END DO
      END DO
      PHASE=PHASE*PHASE1
      END DO
C
C carry out Fourier transform to obtain real data set represented
C as complex array A
      CALL FFT3C(MX,MY,NX,NY,NZ/2,A,ERROR)
C
      RETURN
      END
C======================================================================
      SUBROUTINE FFT3RC2(MX,MY,NX,NY,NZ,A,CWK,ERROR)
C
C Computes 3-dimensional fast Fourier transform real -> complex.
C complex array version.
C
C    X(x+1,y+1,z+1) = SUM {h=0...NX-1,k=0...NY-1,l=0...NZ-1}
C          A(h+1,k+1,l+1) *EXP(2*PI*i*(x*h/NX+y*k/NY+z*l/NZ))
C
C of a real valued sequence a(x,y,z), represented as the
C complex array A with
C
C     A(x,y,z) = a(x,y,{2z-1}) + i*a(x,y,{2z})
C
C where x=1,...,NX, y=1,...,NY, z=1,...,NZ/2.  The result X is
C returned in array A(x,y,z), where  {x=1,...,NX, y=1,...,NY,
C z=1,...,NZ/2+1}, i.e. only one hemisphere of data is returned.
C The remaining coefficients can by obtained by using
C
C     X(x,y,NZ+2-z) = conjg[X(x,y,z)]
C
C CWK is a work array of dimension NX
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER MX, MY, NX, NY, NZ
      DOUBLE COMPLEX A(MX,MY,NZ/2+1)
      DOUBLE COMPLEX CWK(NX)
      LOGICAL ERROR
C local
      INTEGER X, Y, Z, YM, ZM, NXX, NYY
      DOUBLE COMPLEX PHASE1, PHASE, ODD, EVEN
C parameters
      DOUBLE PRECISION HALF, ONE, ZERO, TPI
      PARAMETER (HALF=0.5D0, ONE=1.0D0, ZERO=0.0D0, TPI=2.0D0*PI)
C begin
      ERROR=.FALSE.
C
C carry out Fourier transform on real data set which is represented as
C complex array A
      CALL FFT3C(MX,MY,NX,NY,NZ/2,A,ERROR)
      IF (.NOT.ERROR) THEN
C
C decompose Fourier transform A
C the inner-most loops are vectorized and over the fast index X
      PHASE1=DCMPLX(DCOS(TPI/NZ),DSIN(TPI/NZ))
      PHASE=DCMPLX(ONE,ZERO)
      DO Y=1,NY
      DO X=1,NX
      A(X,Y,NZ/2+1)=A(X,Y,1)
      END DO
      END DO
C
      DO Z=1,NZ/4+1
      ZM=NZ/2+2-Z
      IF (ZM.EQ.Z) THEN
      NYY=NY/2+1
      ELSE
      NYY=NY
      END IF
      DO Y=1,NYY
      YM=MOD(NY+1-Y,NY)+1
      IF (YM.EQ.Y.AND.ZM.EQ.Z) THEN
      NXX=NX/2+1
      ELSE
      NXX=NX
      END IF
      CWK(1)=A(1,YM,ZM)
      DO X=2,NXX
      CWK(X)=A(NX-X+2,YM,ZM)
      END DO
      DO X=1,NXX
      EVEN=A(X,Y,Z)+DCONJG(CWK(X))
      ODD=-DCONJG(A(X,Y,Z))+CWK(X)
      ODD=DCMPLX(DIMAG(ODD),DBLE(ODD))*PHASE
      A(X,Y,Z)=HALF*(EVEN+ODD)
      CWK(X)=DCONJG(HALF*(EVEN-ODD))
      END DO
      A(1,YM,ZM)=CWK(1)
      DO X=2,NXX
      A(NX-X+2,YM,ZM)=CWK(X)
      END DO
      END DO
      PHASE=PHASE*PHASE1
      END DO
      END IF
C
      RETURN
      END
