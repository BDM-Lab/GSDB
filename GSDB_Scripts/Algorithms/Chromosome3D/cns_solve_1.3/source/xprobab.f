C======================================================================
      SUBROUTINE XDOGETFOM(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_fom phase probability function.
C     GET_FOM[PHIStep=<real>,CEN360=<logical>](PA,PB,PC,PD).
C User has to supply PA,PB,PC,PD.  No error checking is done.
C
C Returns Integ ( exp( i phi )  exp( pa cos(phi)   + pb sin(phi)
C                          + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C           /
C          Integ (   exp(pa cos(phi)   + pb sin(phi)
C                      + pc cos(2 phi) + pd sin(2 phi) ) )  d phi
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
      LOGICAL CEN360
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS(1))
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C get cen360 logical flag
      IF (DBLE(ARGS(2)).GT.RSMALL) THEN
      CEN360=.TRUE.
      ELSE
      CEN360=.FALSE.
      END IF
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOGETFO2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),CEN360)
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETFO2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,CEN360)
C
C get_fom phase probability function.  GET_FOM(PA,PB,PC,PD).
C User has to supply PA,PB,PC,PD.  No error checking is done.
C
C Returns Int ( exp( i phi )  pa cos(phi) + pb sin(phi)
C                   + pc cos(2 phi) + pd sin(2 phi) ) d phi
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL CEN360
C local
      INTEGER I, TEMP, IISYM, IPHI
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHAS, PHI, WORK
      DOUBLE COMPLEX CTEMP
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C
      VLEVEL=VLEVEL-3
!$omp parallel do default(shared)
!$omp& private(i,pa,pb,pc,pd,work,iphi,pk,ctemp,iisym,phas,temp)
      DO I=1,N
      PA=DBLE(VSTACK(I,VLEVEL))
      PB=DBLE(VSTACK(I,VLEVEL+1))
      PC=DBLE(VSTACK(I,VLEVEL+2))
      PD=DBLE(VSTACK(I,VLEVEL+3))
C
      IF (TYPE(INDEX(I)).GT.0.OR.CEN360) THEN
C
C acentric reflection or 360 degree-range anomalous centric!!
C
C get normalization constant
      WORK=ZERO
      DO IPHI=1,NPHI
      WORK=MAX(WORK,PA*COSP(IPHI)+PB*SINP(IPHI)
     & +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
      END DO
C
      PK=ZERO
      DO IPHI=1,NPHI
      PK=PK+EXP(-WORK+PA*COSP(IPHI)+PB*SINP(IPHI)
     & +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
      END DO
      PK=PK*PHISTEP
      PK=-LOG(PK)-WORK
C
      CTEMP=DCMPLX(ZERO,ZERO)
      DO IPHI=1,NPHI
C
      CTEMP=CTEMP+DCMPLX(COSP(IPHI),SINP(IPHI))
     &      *EXP(
     & PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     & +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
      VSTACK(I,VLEVEL)=CTEMP*PHISTEP
C
      ELSE
C
C centric reflection !!
C
C get the centric phase shift
      TEMP=ABS(TYPE(INDEX(I)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the centric phase shift and store in PHAS
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(I))
     & +XRSYMM(IISYM,2,4)*XRK(INDEX(I))
     & +XRSYMM(IISYM,3,4)*XRL(INDEX(I)) ) ) / RTH
C
C get the normalization constant
      WORK=MAX(+ PA * COS(PHAS) + PB * SIN(PHAS),
     & - PA * COS(PHAS) - PB * SIN(PHAS))
C
      PK=  EXP(+ PA * COS(PHAS) + PB * SIN(PHAS) -WORK)
     & + EXP(- PA * COS(PHAS) - PB * SIN(PHAS) -WORK)
C
      PK=-LOG(PK)-WORK
C
C store figure of merit in VSTACK(I,VLEVEL)
      VSTACK(I,VLEVEL)=
     & DCMPLX(DCOS(PHAS),DSIN(PHAS)) *(
     & EXP(PK + PA * COS(PHAS) + PB * SIN(PHAS))
     & -EXP(PK - PA * COS(PHAS) - PB * SIN(PHAS)))
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETNRM(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_norm phase probability function.
C   GET_NORM[PHIStep=<real>,CEN360=<logical>](PA,PB,PC,PD).
C User has to supply PA,PB,PC,PD.  No error checking is done.
C
C Returns Log(Integ ( exp( pa cos(phi)   + pb sin(phi)
C                    + pc cos(2 phi) + pd sin(2 phi) ) ) d phi )
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
      LOGICAL CEN360
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS(1))
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C get cen360 logical flag
      IF (DBLE(ARGS(2)).GT.RSMALL) THEN
      CEN360=.TRUE.
      ELSE
      CEN360=.FALSE.
      END IF
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOGETNR2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),CEN360)
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETNR2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,CEN360)
C
C See routine XDOGETNRM above.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL CEN360
C local
      INTEGER I, TEMP, IISYM, IPHI
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHAS, PHI, WORK
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C
      VLEVEL=VLEVEL-3
!$omp parallel do default(shared)
!$omp& private(i,pa,pb,pc,pd,work,iphi,pk,temp,iisym,phas)
      DO I=1,N
      PA=DBLE(VSTACK(I,VLEVEL))
      PB=DBLE(VSTACK(I,VLEVEL+1))
      PC=DBLE(VSTACK(I,VLEVEL+2))
      PD=DBLE(VSTACK(I,VLEVEL+3))
C
      IF (TYPE(INDEX(I)).GT.0.OR.CEN360) THEN
C
C acentric reflection or 360 degree-range anomalous centric!!
C
C get normalization constant
      WORK=ZERO
      DO IPHI=1,NPHI
      WORK=MAX(WORK,PA*COSP(IPHI)+PB*SINP(IPHI)
     & +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
      END DO
C
      PK=ZERO
      DO IPHI=1,NPHI
      PK=PK+EXP(-WORK+PA*COSP(IPHI)+PB*SINP(IPHI)
     & +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
      END DO
C
      PK=PK*PHISTEP
      PK=-LOG(PK)-WORK
C
      VSTACK(I,VLEVEL)=DCMPLX(PK,ZERO)
C
      ELSE
C
C centric reflection !!
C
C get the centric phase shift
      TEMP=ABS(TYPE(INDEX(I)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the centric phase shift and store in PHAS
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(I))
     & +XRSYMM(IISYM,2,4)*XRK(INDEX(I))
     & +XRSYMM(IISYM,3,4)*XRL(INDEX(I)) ) ) / RTH
C
C get the normalization constant
      WORK=MAX(+ PA * COS(PHAS) + PB * SIN(PHAS),
     & - PA * COS(PHAS) - PB * SIN(PHAS))
C
      PK=  EXP(+ PA * COS(PHAS) + PB * SIN(PHAS) -WORK)
     & + EXP(- PA * COS(PHAS) - PB * SIN(PHAS) -WORK)
C
      PK=-LOG(PK)-WORK
C
C store pk in VSTACK(I,VLEVEL)
      VSTACK(I,VLEVEL)=DCMPLX(PK,ZERO)
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_ml lack-of-closure maximum likelihood function.
C GET_ML[PHIStep=<real>,CEN360=<logical>]
C              (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ ( (abs(fh+combine(fp,_phi))-fph)^2
C                  exp(pk + pa cos(phi)   + pb sin(phi)
C                         + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
      LOGICAL CEN360
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
      NPHI=NINT(S360/PHISTEP)
C
C get cen360 logical flag
      IF (DBLE(ARGS(2)).GT.RSMALL) THEN
      CEN360=.TRUE.
      ELSE
      CEN360=.FALSE.
      END IF
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOGETML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),CEN360)
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,CEN360)
C
C get_ml lack-of-closure maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL CEN360
C local
      INTEGER I, TEMP, IISYM, IPHI
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHAS, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH
      DOUBLE PRECISION FP, FPH
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,ctemp,iphi)
!$omp& private(temp,iisym,phas)
      DO I=1,N
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
      IF (TYPE(INDEX(I)).GT.0.OR.CEN360) THEN
C
C acentric reflection or 360 degree-range anomalous centric!!
C
      CTEMP=DCMPLX(ZERO,ZERO)
      DO IPHI=1,NPHI
C
      CTEMP=CTEMP
     &  +(ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))-FPH)**2
     &  *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &         +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
      VSTACK(I,VLEVEL)=WEIGHT*CTEMP*PHISTEP
C
      ELSE
C
C centric reflection !!
C
C get the centric phase shift
      TEMP=ABS(TYPE(INDEX(I)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the centric phase shift and store in PHAS
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(I))
     &           +XRSYMM(IISYM,2,4)*XRK(INDEX(I))
     &           +XRSYMM(IISYM,3,4)*XRL(INDEX(I)) ) ) / RTH
C
C store GET_ML in VSTACK(I,VLEVEL)
      VSTACK(I,VLEVEL)=WEIGHT*
     & ( (ABS(FH+FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)**2
     & *EXP(PK + PA * DCOS(PHAS) + PB * DSIN(PHAS))
     & + (ABS(FH-FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)**2
     & *EXP(PK - PA * DCOS(PHAS) - PB * DSIN(PHAS)))
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_dml derivative of lack-of-closure maximum likelihood function.
C GET_DML[PHIStep=<real>,CEN360=<logical>]
C                    (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ ( (2*(abs(fh+combine(fp,_phi))-fph)
C                       *(fh+combine(fp,_phi))
C                       /abs(fh+combine(fp,_phi))
C                  exp(pk + pa cos(phi)   + pb sin(phi)
C                         + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
      LOGICAL CEN360
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
      NPHI=NINT(S360/PHISTEP)
C
C get cen360 logical flag
      IF (DBLE(ARGS(2)).GT.RSMALL) THEN
      CEN360=.TRUE.
      ELSE
      CEN360=.FALSE.
      END IF
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOGETDML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),CEN360)
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,CEN360)
C
C get_ml derivative of lack-of-closure maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL CEN360
C local
      INTEGER I, TEMP, IISYM, IPHI
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHAS, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH, CCTEMP
      DOUBLE PRECISION FP, FPH
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,ctemp,iphi)
!$omp& private(cctemp,temp,iisym,phas)
      DO I=1,N
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
      IF (TYPE(INDEX(I)).GT.0.OR.CEN360) THEN
C
C acentric reflection or 360 degree-range anomalous centric!!
C
      CTEMP=DCMPLX(ZERO,ZERO)
      DO IPHI=1,NPHI
C
      CCTEMP=(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))
C
      CTEMP=CTEMP
     & +(  (ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)) )-FPH)   )
     & *CCTEMP/MAX(RSMALL,ABS(CCTEMP))
     & *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &        +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
      VSTACK(I,VLEVEL)=TWO*WEIGHT*CTEMP*PHISTEP
C
      ELSE
C
C centric reflection !!
C
C get the centric phase shift
      TEMP=ABS(TYPE(INDEX(I)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the centric phase shift and store in PHAS
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(I))
     & +XRSYMM(IISYM,2,4)*XRK(INDEX(I))
     & +XRSYMM(IISYM,3,4)*XRL(INDEX(I)) ) ) / RTH
C
      CCTEMP=(FH+FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)))
      VSTACK(I,VLEVEL)=
     &    (ABS(FH+FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)
     &    *CCTEMP/MAX(RSMALL,ABS(CCTEMP))
     &    *EXP(PK + PA * DCOS(PHAS) + PB * DSIN(PHAS))
      CCTEMP=(FH-FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)))
      VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL)+
     &    (ABS(FH-FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)
     &    *CCTEMP/MAX(RSMALL,ABS(CCTEMP))
     &    *EXP(PK - PA * DCOS(PHAS) - PB * DSIN(PHAS))
      VSTACK(I,VLEVEL)=WEIGHT*TWO*VSTACK(I,VLEVEL)
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDSML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_dsml derivative of lack-of-closure maximum likelihood function.
C GET_DSML[PHIStep=<real>,CEN360=<logical>]
C               (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ ( (2*(abs(fh+combine(fp,_phi))-fph)
C                  *exp(pk + pa cos(phi)   + pb sin(phi)
C                         + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI
      DOUBLE PRECISION PHISTEP
      LOGICAL CEN360
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS(1)).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS(1))
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C get cen360 logical flag
      IF (DBLE(ARGS(2)).GT.RSMALL) THEN
      CEN360=.TRUE.
      ELSE
      CEN360=.FALSE.
      END IF
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
C
      CALL XDOGETDSML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),CEN360)
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDSML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,CEN360)
C
C get_ml derivative of lack-of-closure maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS(*)
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      LOGICAL CEN360
C local
      INTEGER I, TEMP, IISYM, IPHI
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHAS, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH
      DOUBLE PRECISION FP, FPH
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,ctemp,iphi)
!$omp& private(temp,iisym,phas)
      DO I=1,N
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
      IF (TYPE(INDEX(I)).GT.0.OR.CEN360) THEN
C
C acentric reflection or 360 degree-range anomalous centric!!
C
      CTEMP=DCMPLX(ZERO,ZERO)
      DO IPHI=1,NPHI
C
      CTEMP=CTEMP
     & +(  (ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)) )-FPH)   )
     & *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &        +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
      VSTACK(I,VLEVEL)=TWO*WEIGHT*CTEMP*PHISTEP
C
      ELSE
C
C centric reflection !!
C
C get the centric phase shift
      TEMP=ABS(TYPE(INDEX(I)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the centric phase shift and store in PHAS
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(I))
     & +XRSYMM(IISYM,2,4)*XRK(INDEX(I))
     & +XRSYMM(IISYM,3,4)*XRL(INDEX(I)) ) ) / RTH
C
      VSTACK(I,VLEVEL)=
     &    (ABS(FH+FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)
     &    *EXP(PK + PA * DCOS(PHAS) + PB * DSIN(PHAS))
      VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL)+
     &    (ABS(FH-FP*DCMPLX(DCOS(PHAS),DSIN(PHAS)) )-FPH)
     &    *EXP(PK - PA * DCOS(PHAS) - PB * DSIN(PHAS))
      VSTACK(I,VLEVEL)=WEIGHT*TWO*VSTACK(I,VLEVEL)
      END IF
C
      END DO
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETAML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_aml fano maximum likelihood function.
C GET_AML[PHIStep=<real>](FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ (
C    (abs(fh+combine(fp,_phi))
C      -abs(friedel(fh)+combine(abs(friedel(fp)),-_phi))
C   -fph+friedel(fph))^2
C      *exp(pk + pa cos(phi)   + pb sin(phi)
C              + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C for acentric reflections.
C
C 0 for centric reflections.
C
C Note: this is a special function, i.e., it requires
C all reflections to be present in the list.
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI, I
      DOUBLE PRECISION PHISTEP
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P, WORK
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C for non-anomalous data, the function is zero
      IF (QHERM) THEN
      VLEVEL=VLEVEL-8
      DO I=1,N
      VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
      END DO
      ELSE
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS)
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
      WORK=ALLHP(INTEG4(XRNREF))
C
      CALL XDOGETAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),HEAP(WORK))
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(WORK,INTEG4(XRNREF))
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,WORK)
C
C get_aml fano maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      INTEGER WORK(*)
C local
      INTEGER I, TEMP, IISYM, IPHI, R2, R22, H, K, L
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH, FDSHIFT, FDFH
      DOUBLE PRECISION FP, FPH, PHAS, FDFP, FDFPH
      LOGICAL FOUND
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C precompute indices and phase shifts for Friedel mates
      DO I=1,XRNREF
      WORK(I)=0
      END DO
      DO I=1,N
      WORK(INDEX(I))=I
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,found,temp)
!$omp& private(iisym,r22,r2,h,k,l,phas,fdshift,fdfp,fdfh)
!$omp& private(fdfph,ctemp,iphi)
      DO I=1,N
C
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
C set the function to zero for centric reflections
C and reflections with no Friedel mate
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=ZERO
C
      IF (TYPE(INDEX(I)).GT.0) THEN
C
C acentric reflection !!
C
C get Friedel mates for FP, FH, and FPH
      FOUND=ABS(TYPE(INDEX(I))).GT.1
      IF (FOUND) THEN
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(INDEX(I)))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R22=TEMP/XRNSYM
      R2=WORK(R22)
      FOUND=R2.GT.0
C
      IF (FOUND) THEN
C
C compute phase shift
      H=XRH(R22)
      K=XRK(R22)
      L=XRL(R22)
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*H
     &              +XRSYMM(IISYM,2,4)*K
     &              +XRSYMM(IISYM,3,4)*L ) / RTH )
      FDSHIFT=ONE/DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
      FDFP=VSTACK(R2,VLEVEL)
      FDFH=VSTACK(R2,VLEVEL+1)*FDSHIFT
      FDFPH=VSTACK(R2,VLEVEL+2)
C
      CTEMP=DCMPLX(ZERO,ZERO)
C
C compute maximum likelihood Fano
      DO IPHI=1,NPHI
C
      CTEMP=CTEMP
     &     +(ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))
     &     -ABS(FDFH+FDFP*DCMPLX(COSP(IPHI),-SINP(IPHI)))
     &     -FPH+FDFPH)**2
     &  *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &         +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
C
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=WEIGHT*CTEMP*PHISTEP
      END IF
      END IF
      END IF
C
      END DO
C
C copy everything from the temporary store to VLEVEL
      DO I=1,N
      VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+8)
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDAML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_daml derivative of fano maximum likelihood function.
C GET_DAML[PHIStep=<real>](FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ (
C    (abs(fh+combine(fp,_phi))
C      -abs(friedel(fh)+combine(friedel(fp),-_phi))
C   -fph+friedel(fph))^2
C      *exp(pk + pa cos(phi)   + pb sin(phi)
C              + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C for acentric reflections.
C
C 0 for centric reflections.
C
C Note: this is a special function, i.e., it requires
C all reflections to be present in the list.
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI, I
      DOUBLE PRECISION PHISTEP
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P, WORK
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C for non-anomalous data, the function is zero
      IF (QHERM) THEN
      VLEVEL=VLEVEL-8
      DO I=1,N
      VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
      END DO
      ELSE
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS)
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
      WORK=ALLHP(INTEG4(XRNREF))
C
      CALL XDOGETDAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),HEAP(WORK))
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(WORK,INTEG4(XRNREF))
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,WORK)
C
C get_aml derivative of fano maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      INTEGER WORK(*)
C local
      INTEGER I, TEMP, IISYM, IPHI, R2, R22, H, K, L
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH, FDSHIFT, FDFH, CCTEMP
      DOUBLE PRECISION FP, FPH, PHAS, FDFP, FDFPH
      LOGICAL FOUND
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C precompute indices and phase shifts for Friedel mates
      DO I=1,XRNREF
      WORK(I)=0
      END DO
      DO I=1,N
      WORK(INDEX(I))=I
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,found,temp)
!$omp& private(iisym,r22,r2,h,k,l,phas,fdshift,fdfp,fdfh)
!$omp& private(fdfph,ctemp,iphi,cctemp)
      DO I=1,N
C
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
C set the function to zero for centric reflections
C and reflections with no Friedel mate
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=ZERO
C
      IF (TYPE(INDEX(I)).GT.0) THEN
C
C acentric reflection !!
C
C get Friedel mates for FP, FH, and FPH
      FOUND=ABS(TYPE(INDEX(I))).GT.1
      IF (FOUND) THEN
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(INDEX(I)))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R22=TEMP/XRNSYM
      R2=WORK(R22)
      FOUND=R2.GT.0
C
      IF (FOUND) THEN
C
C compute phase shift
      H=XRH(R22)
      K=XRK(R22)
      L=XRL(R22)
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*H
     &              +XRSYMM(IISYM,2,4)*K
     &              +XRSYMM(IISYM,3,4)*L ) / RTH )
      FDSHIFT=ONE/DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
      FDFP=VSTACK(R2,VLEVEL)
      FDFH=VSTACK(R2,VLEVEL+1)*FDSHIFT
      FDFPH=VSTACK(R2,VLEVEL+2)
C
      CTEMP=DCMPLX(ZERO,ZERO)
C
C compute maximum likelihood Fano
      DO IPHI=1,NPHI
C
      CCTEMP=(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))
      CTEMP=CTEMP
     &     +(ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))
     &     -ABS(FDFH+FDFP*DCMPLX(COSP(IPHI),-SINP(IPHI)))
     &     -FPH+FDFPH)*CCTEMP/MAX(RSMALL,ABS(CCTEMP))
     &  *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &         +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
C
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=TWO*WEIGHT*CTEMP*PHISTEP
      END IF
      END IF
      END IF
C
      END DO
C
C copy everything from the temporary store to VLEVEL
      DO I=1,N
      VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+8)
      END DO
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDSAML(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX,
     & XRNREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY)
C
C get_daml derivative of fano maximum likelihood function.
C GET_DSAML[PHIStep=<real>](FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD).
C User has to supply all arguments.  No error checking is done.
C
C Returns
C  weight *
C  Integ (
C    (abs(fh+combine(fp,_phi))
C      -abs(friedel(fh)+combine(friedel(fp),-_phi))
C   -fph+friedel(fph))^2
C      *exp(pk + pa cos(phi)   + pb sin(phi)
C              + pc cos(2 phi) + pd sin(2 phi) ) ) d phi
C for acentric reflections.
C
C 0 for centric reflections.
C
C Note: this is a special function, i.e., it requires
C all reflections to be present in the list.
C
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS
      INTEGER INDEX(*), XRNREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      INTEGER NPHI, I
      DOUBLE PRECISION PHISTEP
C pointer
      INTEGER SINP, COSP, SIN2P, COS2P, WORK
C parameter
      DOUBLE PRECISION ZERO, S360, ONE
      PARAMETER (ZERO=0.0D0, S360=360D0, ONE=1.0D0)
C begin
C
C for non-anomalous data, the function is zero
      IF (QHERM) THEN
      VLEVEL=VLEVEL-8
      DO I=1,N
      VSTACK(I,VLEVEL)=DCMPLX(ZERO,ZERO)
      END DO
      ELSE
C
C
C get integration step size (phistep)
      IF (DBLE(ARGS).LT.RSMALL) THEN
      PHISTEP=ONE
      ELSE
      PHISTEP=DBLE(ARGS)
      END IF
      NPHI=NINT(S360/PHISTEP)
C
C allocate space for precomputed SIN, COS, etc arrays
      SINP=ALLHP(IREAL8(NPHI))
      COSP=ALLHP(IREAL8(NPHI))
      SIN2P=ALLHP(IREAL8(NPHI))
      COS2P=ALLHP(IREAL8(NPHI))
      WORK=ALLHP(INTEG4(XRNREF))
C
      CALL XDOGETDSAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,INDEX,
     & XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     & QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,HEAP(SINP),HEAP(COSP),HEAP(SIN2P),
     & HEAP(COS2P),HEAP(WORK))
C
C de-allocate space for precomputed SIN, COS, etc arrays
      CALL FREHP(WORK,INTEG4(XRNREF))
      CALL FREHP(SINP,IREAL8(NPHI))
      CALL FREHP(COSP,IREAL8(NPHI))
      CALL FREHP(SIN2P,IREAL8(NPHI))
      CALL FREHP(COS2P,IREAL8(NPHI))
C
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOGETDSAML2(VLEVEL,VMAX,VSTACK,N,ARGS,PHISTEP,NPHI,
     & INDEX,
     & XRNREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     & XRSYMM,XRITSY,SINP,COSP,SIN2P,COS2P,WORK)
C
C get_aml derivative of fano maximum likelihood function.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      DOUBLE COMPLEX ARGS
      DOUBLE PRECISION PHISTEP
      INTEGER NPHI
      INTEGER INDEX(*), XRNREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      DOUBLE PRECISION SINP(*), COSP(*), SIN2P(*), COS2P(*)
      INTEGER WORK(*)
C local
      INTEGER I, TEMP, IISYM, IPHI, R2, R22, H, K, L
      DOUBLE PRECISION PA, PB, PC, PD, PK, RTH, PHI, WEIGHT
      DOUBLE COMPLEX CTEMP, FH, FDSHIFT, FDFH
      DOUBLE PRECISION FP, FPH, PHAS, FDFP, FDFPH
      LOGICAL FOUND
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720, ONE, TWO
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
      PARAMETER (ONE=1.0D0, TWO=2.0D0)
C begin
      RTH=XRSYTH
C
C precompute SIN, COS, etc.
      DO IPHI=1,NPHI
      PHI=(IPHI-1)*PHISTEP
      SINP(IPHI)=SIN(PHI*PI/S180)
      COSP(IPHI)=COS(PHI*PI/S180)
      SIN2P(IPHI)=SIN(TWO*PHI*PI/S180)
      COS2P(IPHI)=COS(TWO*PHI*PI/S180)
      END DO
C
C precompute indices and phase shifts for Friedel mates
      DO I=1,XRNREF
      WORK(I)=0
      END DO
      DO I=1,N
      WORK(INDEX(I))=I
      END DO
C
C (FP,FH,FPH,WEIGHT,PK,PA,PB,PC,PD)
C
      VLEVEL=VLEVEL-8
!$omp parallel do default(shared)
!$omp& private(i,fp,fh,fph,weight,pk,pa,pb,pc,pd,found,temp)
!$omp& private(iisym,r22,r2,h,k,l,phas,fdshift,fdfp,fdfh)
!$omp& private(fdfph,ctemp,iphi)
      DO I=1,N
C
      FP=VSTACK(I,VLEVEL)
      FH=VSTACK(I,VLEVEL+1)
      FPH=VSTACK(I,VLEVEL+2)
      WEIGHT=VSTACK(I,VLEVEL+3)
      PK=DBLE(VSTACK(I,VLEVEL+4))
      PA=DBLE(VSTACK(I,VLEVEL+5))
      PB=DBLE(VSTACK(I,VLEVEL+6))
      PC=DBLE(VSTACK(I,VLEVEL+7))
      PD=DBLE(VSTACK(I,VLEVEL+8))
C
C set the function to zero for centric reflections
C and reflections with no Friedel mate
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=ZERO
C
      IF (TYPE(INDEX(I)).GT.0) THEN
C
C acentric reflection !!
C
C get Friedel mates for FP, FH, and FPH
      FOUND=ABS(TYPE(INDEX(I))).GT.1
      IF (FOUND) THEN
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(INDEX(I)))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R22=TEMP/XRNSYM
      R2=WORK(R22)
      FOUND=R2.GT.0
C
      IF (FOUND) THEN
C
C compute phase shift
      H=XRH(R22)
      K=XRK(R22)
      L=XRL(R22)
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*H
     &              +XRSYMM(IISYM,2,4)*K
     &              +XRSYMM(IISYM,3,4)*L ) / RTH )
      FDSHIFT=ONE/DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
      FDFP=VSTACK(R2,VLEVEL)
      FDFH=VSTACK(R2,VLEVEL+1)*FDSHIFT
      FDFPH=VSTACK(R2,VLEVEL+2)
C
      CTEMP=DCMPLX(ZERO,ZERO)
C
C compute maximum likelihood Fano
      DO IPHI=1,NPHI
C
      CTEMP=CTEMP
     &     +(ABS(FH+FP*DCMPLX(COSP(IPHI),SINP(IPHI)))
     &     -ABS(FDFH+FDFP*DCMPLX(COSP(IPHI),-SINP(IPHI)))
     &     -FPH+FDFPH)
     &  *EXP(PK+PA*COSP(IPHI)+PB*SINP(IPHI)
     &         +PC*COS2P(IPHI)+PD*SIN2P(IPHI))
C
      END DO
C
C temporarily store it in VLEVEL+8
      VSTACK(I,VLEVEL+8)=TWO*WEIGHT*CTEMP*PHISTEP
      END IF
      END IF
      END IF
C
      END DO
C
C copy everything from the temporary store to VLEVEL
      DO I=1,N
      VSTACK(I,VLEVEL)=VSTACK(I,VLEVEL+8)
      END DO
C
      RETURN
      END
