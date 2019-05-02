C======================================================================
      SUBROUTINE XDOFRIED(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRMREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY,RPNTYP,RPNDOM,MODE)
C
C
C Stores value of Friedel mate in array (mode=1) or returns
C a logical flag array with the subset of reflections for
C which selected Friedel mates exist.
C
C Note:    This is a special operation, i.e., N has to include
C          all structure factor elements (see routine XDOSPCL).
C
C Note:    This routine assumes that the reflection data
C          are stored as an asymmetric unit (using XRASYM).
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*), XRMREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      CHARACTER*2 RPNTYP, RPNDOM
      INTEGER MODE
C pointer
      INTEGER STORE, WORK
C begin
      WORK=ALLHP(INTEG4(XRMREF))
      STORE=ALLHP(ICPLX8(N))
      CALL XDOFRIE2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &      XRMREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     &      QHERM,XRNSYM,XRMSYM,XRSYTH,
     &      XRSYMM,XRITSY,RPNTYP,RPNDOM,HEAP(STORE),HEAP(WORK),MODE)
      CALL FREHP(WORK,INTEG4(XRMREF))
      CALL FREHP(STORE,ICPLX8(N))
      RETURN
      END
C======================================================================
      SUBROUTINE XDOFRIE2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRMREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY,RPNTYP,RPNDOM,STORE,WORK,MODE)
C
C See XDOFRIE2 above.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRMREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      CHARACTER*2 RPNTYP, RPNDOM
      DOUBLE COMPLEX STORE(*)
      INTEGER WORK(*)
      INTEGER MODE
C local
      INTEGER REFLCT, R2, R22, H, K, L, IISYM, TEMP
      LOGICAL FOUND, QSHIFT
      DOUBLE PRECISION PHAS, RTH
C parameter
      DOUBLE PRECISION ZERO, TWO
      PARAMETER (ZERO=0.0D0, TWO=2.0D0)
C begin
C
C
C
C apply phase shift only to operands of type complex.
      QSHIFT=RPNTYP.EQ.'DC'
C
      RTH=XRSYTH
C
      IF (QHERM) THEN
C
C non-anomalous data:
      IF (MODE.EQ.1) THEN
      DO REFLCT=1,N
      VSTACK(REFLCT,VLEVEL)=DCONJG(VSTACK(REFLCT,VLEVEL))
      END DO
      END IF
C
C Remark: no operation required for logical (MODE=2).
C If hermitian symmetry is present all Friedel mates
C exist by definition.
C
      ELSE
C
C MODE=1, anomalous data
      IF (MODE.EQ.1) THEN
      DO REFLCT=1,XRMREF
      WORK(REFLCT)=0
      END DO
      DO REFLCT=1,N
      WORK(INDEX(REFLCT))=REFLCT
      END DO
      DO REFLCT=1,N
C
      FOUND=ABS(TYPE(INDEX(REFLCT))).GT.1
      IF (FOUND) THEN
C
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(INDEX(REFLCT)))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R22=TEMP/XRNSYM
      R2=WORK(R22)
      FOUND=R2.GT.0
C
      IF (FOUND) THEN
C
C store value of Friedel mate
      STORE(REFLCT)=VSTACK(R2,VLEVEL)
C
C apply phase shift if required
      IF (IISYM.NE.1.AND.QSHIFT) THEN
      H=XRH(R22)
      K=XRK(R22)
      L=XRL(R22)
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*H
     &              +XRSYMM(IISYM,2,4)*K
     &              +XRSYMM(IISYM,3,4)*L ) / RTH )
      STORE(REFLCT)=STORE(REFLCT)/DCMPLX(DCOS(PHAS),DSIN(PHAS))
      END IF
      END IF
C
      END IF
C
      IF (.NOT.FOUND.AND.MODE.EQ.1) THEN
      IF (WRNLEV.GE.5) THEN
      WRITE(6,'(A,3I6,A)')
     &  ' XDOFRIED-warning: Friedel mate for reflection ',
     &  XRH(INDEX(REFLCT)),XRK(INDEX(REFLCT)),XRL(INDEX(REFLCT)),
     &  ' not found in selected set.'
      END IF
      STORE(REFLCT)=DCMPLX(ZERO,ZERO)
      END IF
C
      END DO
C
C now we have to copy things.
      DO REFLCT=1,N
      VSTACK(REFLCT,VLEVEL)=STORE(REFLCT)
      END DO
C
      ELSE
C MODE=2, anomalous data
C
      DO REFLCT=1,XRMREF
      WORK(REFLCT)=0
      END DO
      DO REFLCT=1,N
      WORK(INDEX(REFLCT))=REFLCT
      END DO
C
      DO REFLCT=1,N
C
C first test: Friedel mate exists of selected element
      FOUND=ABS(TYPE(INDEX(REFLCT))).GT.1.AND.LSTACK(REFLCT,VLEVEL)
      IF (FOUND) THEN
C
C decode information about symmetry operator and reflection index
      TEMP=ABS(TYPE(INDEX(REFLCT)))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R22=TEMP/XRNSYM
      R2=WORK(R22)
C second test: Friedel mate is in selected set
      FOUND=R2.GT.0.AND.LSTACK(R2,VLEVEL)
      END IF
C
      LSTACK(REFLCT,VLEVEL)=FOUND
      END DO
C
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE XTYPMAP(PERM,TYPE,WORK,N,XRNSYM)
C
C Routine maps the type array elements after the reflection
C list has been sorted.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER PERM(*), TYPE(*), WORK(*), N, XRNSYM
C local
      INTEGER REFLCT, IISYM, R2, TEMP
C begin
      DO REFLCT=1,N
      IF (1.GT.PERM(REFLCT).OR.PERM(REFLCT).GT.N) THEN
      WRITE(6,'(A,I6,I6)') ' XTYPMAP: internal error no. 1.',REFLCT,
     &  PERM(REFLCT)
      END IF
      WORK(PERM(REFLCT))=REFLCT
      END DO
C
      DO REFLCT=1,N
      IF (ABS(TYPE(REFLCT)).GT.1) THEN
      TEMP=ABS(TYPE(REFLCT))-1
      IISYM=MOD(TEMP,XRNSYM)+1
      R2=TEMP/XRNSYM
      R2=WORK(R2)
C
      IF (TYPE(REFLCT).GT.0) THEN
      TYPE(REFLCT)=(R2*XRNSYM+IISYM)
      ELSE
      TYPE(REFLCT)=-(R2*XRNSYM+IISYM)
      END IF
      END IF
      END DO
C
C
      RETURN
      END
C======================================================================
      SUBROUTINE XDOCPHAS(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRMREF,HPH,HPK,HPL,HPTYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY)
C
C Computes phase for centric reflections (up to +-180.).
C Set to zero for all acentric reflections.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*), XRMREF, HPH, HPK, HPL, HPTYPE
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C begin
      CALL XDOCPHA2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRMREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),HEAP(HPTYPE),
     &         QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY)
      RETURN
      END
C======================================================================
      SUBROUTINE XDOCPHA2(VLEVEL,VMAX,VSTACK,LSTACK,N,INDEX,
     &         XRMREF,XRH,XRK,XRL,TYPE,QHERM,XRNSYM,XRMSYM,XRSYTH,
     &         XRSYMM,XRITSY)
C
C See routine XDOCPHAS above.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*)
      LOGICAL LSTACK(N,*)
      INTEGER INDEX(*), XRMREF, XRH(*), XRK(*), XRL(*), TYPE(*)
      LOGICAL QHERM
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
C local
      DOUBLE PRECISION RTH, PHAS
      INTEGER TEMP, IISYM, REFLCT
C parameter
      DOUBLE PRECISION ZERO, S180, S360, S720
      PARAMETER (ZERO=0.0D0, S180=180.D0, S360=360.D0, S720=720.D0)
C begin
C
      RTH=XRSYTH
C
      VLEVEL=VLEVEL+1
C
      DO REFLCT=1,N
C
      IF (TYPE(INDEX(REFLCT)).GT.0) THEN
C
C acentric reflection
      VSTACK(REFLCT,VLEVEL)=DCMPLX(ZERO,ZERO)
C
      ELSE
C
C centric reflection
      TEMP=ABS(TYPE(INDEX(REFLCT)))-1
C
C get the symmetry operator number
      IISYM=MOD(TEMP,XRNSYM)+1
C
C get the phase shift
      PHAS=PI*( ( XRSYMM(IISYM,1,4)*XRH(INDEX(REFLCT))
     &           +XRSYMM(IISYM,2,4)*XRK(INDEX(REFLCT))
     &           +XRSYMM(IISYM,3,4)*XRL(INDEX(REFLCT)) ) ) / RTH
      PHAS=MOD(PHAS*S180/PI+S720,S180)
      VSTACK(REFLCT,VLEVEL)=DCMPLX(PHAS,ZERO)
      END IF
C
      END DO
C
C
      RETURN
      END
C
