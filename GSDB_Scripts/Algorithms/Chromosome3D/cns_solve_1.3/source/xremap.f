C ===================================================================
      SUBROUTINE XDORMAP(VLEVEL,VMAX,VSTACK,N,ARGS,INDEX)
C
C Routine remaps reflections by application of a matrix to the indices
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C I/O
      INTEGER VLEVEL, VMAX, N
      DOUBLE COMPLEX VSTACK(N,*), ARGS(*)
      INTEGER INDEX(*)
C local
      INTEGER I, J, K, MDIM, MATRIX(3,3)
      INTEGER HMAX, KMAX, LMAX, HMIN, KMIN, LMIN
      DOUBLE PRECISION DET
C parameter
      DOUBLE PRECISION SMALL, ONE
      PARAMETER (SMALL=1.0D-4,ONE=1.0D0)
C pointer
      INTEGER REMAP
C begin
C
C initialize matrix
      DO I=1,3
         DO J=1,3
            MATRIX(I,J)=0
         END DO
      END DO
C
C construct the transformation matrix
      DO K=1,3
         I=INT(DBLE(ARGS(K)))
         J=INT(DIMAG(ARGS(K)))
         IF (ABS(I).GE.1.AND.ABS(I).LE.3) MATRIX(K,ABS(I))=SIGN(1,I)
         IF (ABS(J).GE.1.AND.ABS(J).LE.3) MATRIX(K,ABS(J))=SIGN(1,J)
      END DO
C
C check matrix determinant
      DET=(MATRIX(1,1)*MATRIX(2,2) -
     &     MATRIX(1,2)*MATRIX(2,1))*MATRIX(3,3)
     &   +(MATRIX(2,1)*MATRIX(3,2) -
     &     MATRIX(2,2)*MATRIX(3,1))*MATRIX(1,3)
     &   +(MATRIX(3,1)*MATRIX(1,2) -
     &     MATRIX(3,2)*MATRIX(1,1))*MATRIX(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      WRITE(6,'(A,F8.5)') ' %XDOREMAP-ERR: invalid determinant =',DET
      CALL WRNDIE(-5,'XDOREMAP',
     &     'unsuitable remapping matrix - not a rotation.')
      END IF
C
C find extent of indices
      CALL XRRR3(XRNREF,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &           HMAX,KMAX,LMAX,HMIN,KMIN,LMIN)
C
      IF (HMAX.EQ.0) HMAX=1
      IF (KMAX.EQ.0) KMAX=1
      IF (LMAX.EQ.0) LMAX=1
      MDIM=(HMAX-HMIN+1)*(KMAX-KMIN+1)*(LMAX-LMIN+1)
C
C allocate space for the remapping array
      REMAP=ALLHP(ICPLX8(MDIM))
C
C do the remapping
      CALL XDORMAP2(N,VLEVEL,VSTACK,MATRIX,
     &              INDEX,HEAP(HPH),HEAP(HPK),HEAP(HPL),
     &              XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &              HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,HEAP(REMAP))
C
C deallocate remapping array
      CALL FREHP(REMAP,ICPLX8(MDIM))
C
      RETURN
      END
C
C======================================================================
C
      SUBROUTINE XDORMAP2(N,VLEVEL,VSTACK,MATRIX,
     &                    INDEX,XRH,XRK,XRL,
     &                    XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY,QHERM,
     &                    HMAX,KMAX,LMAX,HMIN,KMIN,LMIN,REMAP)
C
C see routine XREMAP above
C
C Author: Paul Adams
C ==================
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C
      INTEGER N, VLEVEL
      DOUBLE COMPLEX VSTACK(N,*)
      INTEGER INDEX(*)
      INTEGER XRH(*), XRK(*), XRL(*)
      INTEGER XRNSYM, XRMSYM, XRSYTH
      INTEGER XRSYMM(XRMSYM,3,4), XRITSY(XRMSYM,3,3)
      LOGICAL QHERM
      INTEGER HMIN, HMAX, KMIN, KMAX, LMIN, LMAX
      DOUBLE COMPLEX REMAP(HMIN:HMAX,KMIN:KMAX,LMIN:LMAX)
      INTEGER MATRIX(3,3)
C local
      INTEGER REFLCT, HH, KK, LL, HH2, KK2, LL2, IFRIED, IISYM
      DOUBLE PRECISION RTH, PHAS
      DOUBLE COMPLEX SHIFT, CTEMP
      LOGICAL SYSAB
C parameter
      DOUBLE PRECISION ZERO, ONE, TWO, S180
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,S180=180.D0)
C begin
C
      RTH=XRSYTH
C
C initialize the remapping array with the current values
C
      DO LL=LMIN,LMAX
      DO KK=KMIN,KMAX
      DO HH=HMIN,HMAX
      REMAP(HH,KK,LL)=DCMPLX(ZERO,ZERO)
      END DO
      END DO
      END DO
C
C remap the indices
      DO REFLCT=1,N
C
C put initial value in temporary space
      CTEMP=VSTACK(REFLCT,VLEVEL)
C
C new index after application of operator
      HH=MATRIX(1,1)*XRH(INDEX(REFLCT))+
     &   MATRIX(1,2)*XRK(INDEX(REFLCT))+
     &   MATRIX(1,3)*XRL(INDEX(REFLCT))
      KK=MATRIX(2,1)*XRH(INDEX(REFLCT))+
     &   MATRIX(2,2)*XRK(INDEX(REFLCT))+
     &   MATRIX(2,3)*XRL(INDEX(REFLCT))
      LL=MATRIX(3,1)*XRH(INDEX(REFLCT))+
     &   MATRIX(3,2)*XRK(INDEX(REFLCT))+
     &   MATRIX(3,3)*XRL(INDEX(REFLCT))
C
C map reflections into the asymmetric unit
      CALL XRASYM(HH,KK,LL,HH2,KK2,LL2,
     &            IISYM,IFRIED,XRNSYM,XRMSYM,XRITSY,QHERM)
C
C check if this is a systematic absence
      CALL XRSYSAB(HH2,KK2,LL2,SYSAB,
     &             XRNSYM,XRMSYM,XRSYTH,XRSYMM,XRITSY)
C
C it this is not an absent reflection then apply shifts etc
      IF (.NOT.SYSAB) THEN
C
C deal with phase shift here
      IF (IISYM.NE.1.OR.IFRIED.EQ.-1) THEN
C
      PHAS=TWO*PI*(( XRSYMM(IISYM,1,4)*HH2
     &              +XRSYMM(IISYM,2,4)*KK2
     &              +XRSYMM(IISYM,3,4)*LL2 ) / RTH )
C
      IF (IFRIED.EQ.-1) THEN
         PHAS=-PHAS
      END IF
C
      SHIFT=DCMPLX(DCOS(PHAS),DSIN(PHAS))
C
      CTEMP=CTEMP*SHIFT
      IF (IFRIED.EQ.-1) CTEMP=DCONJG(CTEMP)
C
      END IF
C
      END IF
C
C fill remapping array with value after transformation
      IF ((HH2.GE.HMIN).AND.(HH2.LE.HMAX).AND.
     &    (KK2.GE.KMIN).AND.(KK2.LE.KMAX).AND.
     &    (LL2.GE.LMIN).AND.(LL2.LE.LMAX)) THEN
      REMAP(HH2,KK2,LL2)=CTEMP
      END IF
C
      END DO
C
C fill VSTACK with the remapped data points
      DO REFLCT=1,N
         HH=XRH(INDEX(REFLCT))
         KK=XRK(INDEX(REFLCT))
         LL=XRL(INDEX(REFLCT))
         VSTACK(REFLCT,VLEVEL)=REMAP(HH,KK,LL)
      END DO
C
      RETURN
      END
C
C ===================================================================
C
