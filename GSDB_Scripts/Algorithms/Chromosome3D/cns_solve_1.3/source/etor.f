      SUBROUTINE EDIHE(EP,EPV,IP,JP,KP,LP,INTERE,NPHI,ACPC,ACPD,ACPB,
     &                 QACPC,QACPD,QACPB,
     &                 MODE,CPO,CPE,ANALYS,WGHT,VWGHT,IAC,
     &                 SEGID,RESID,TYPE)
C
C Calculates dihedral energy.
C Front-end for actual energy routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'update.inc'
      INCLUDE 'param.inc'
      DOUBLE PRECISION EP,EPV
      DOUBLE PRECISION WGHT,VWGHT
      INTEGER IP(*), JP(*), KP(*), LP(*), NPHI, INTERE(*)
      DOUBLE PRECISION ACPC(*)
      INTEGER ACPD(*)
      DOUBLE PRECISION ACPB(*)
      LOGICAL QACPC(*), QACPD(*), QACPB(*)
      CHARACTER*(*) MODE
      DOUBLE PRECISION CPO(*)
      INTEGER CPE(*)
      CHARACTER*4 ANALYS, IAC(*), SEGID(*), RESID(*), TYPE(*)
C
C fill the unkown ACPC, ACPB, ACPD values with type-based parameters
      IF (UPDIHE) THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      UPDIHE=.FALSE.
      END IF
C
      CALL ETOR(EP,EPV,IP,JP,KP,LP,INTERE,NPHI,ACPC,ACPD,ACPB,MODE,
     &                 CPO,CPE,ANALYS,WGHT,VWGHT)
      RETURN
      END
C======================================================================
      SUBROUTINE EIMPR(EM,EMV,IM,JM,KM,LM,INTERE,NIMPHI,ACIC,ACID,ACIB,
     &            QACIC,QACID,QACIB,MODE,CPO,CPE,ANALYS,WGHT,VWGHT,IAC,
     &            SEGID,RESID,TYPE)
C
C Calculates improper energy.
C Front-end for actual energy routine.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'update.inc'
      INCLUDE 'param.inc'
      DOUBLE PRECISION EM,EMV
      DOUBLE PRECISION WGHT,VWGHT
      INTEGER IM(*), JM(*), KM(*), LM(*), NIMPHI, INTERE(*)
      DOUBLE PRECISION ACIC(*)
      INTEGER ACID(*)
      DOUBLE PRECISION ACIB(*)
      LOGICAL QACIC(*), QACID(*), QACIB(*)
      CHARACTER*(*) MODE
      DOUBLE PRECISION CPO(*)
      INTEGER CPE(*)
      CHARACTER*4 ANALYS, IAC(*), SEGID(*), RESID(*), TYPE(*)
C
C fill the unkown ACPC, ACPB, ACPD values with type-based parameters
      IF (UPIMPR) THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,
     &            QACIC,QACIB,QACID,IAC,SEGID,RESID,TYPE)
      UPIMPR=.FALSE.
      END IF
C
      CALL ETOR(EM,EMV,IM,JM,KM,LM,INTERE,NIMPHI,ACIC,ACID,ACIB,MODE,
     &                 CPO,CPE,ANALYS,WGHT,VWGHT)
      RETURN
      END
C======================================================================
      SUBROUTINE ETOR(EP,EPV,IP,JP,KP,LP,INTERE,NPHI,CPC,CPD,CPB,MODE,
     &                CPO,CPE,ANALYS,WGHT,VWGHT)
C
C Computes torsion angles and first derivatives in a vectorizable way.
C ====================================================================
C
C Functional form:
C
C IF CPD(IPHI) not equal zero, we use a sinusoidal function:
C
C       c * (1.0+cos ( n*phi + delta ))
C       ===============================
C        where phi is defined as the dihedral angle (ri, rj, rk, rl).
C
C ELSE we use a harmonic function, where we distinguish between two
C   cases:
C   reduced form:  const * ( phi - phi0 ) **2
C   general form:  const * well(phi-phi0,offset) **exponent
C   where
C                     well(a,b) =    a-b   a > b
C                     well(a,b) =     0    -b < a < b
C                     well(a,b) =    a+b   a < -b
C
C   where phi is defined as the dihedral angle (ri, rj, rk, rl).
C
C    The algorithm uses two alternative ways to compute the derivatives
C with regard to ri,j,k,l.  In this way artificial singularities of
C the numerical expressions are avoided ( at phi = 0.0 and pi/2 ).
C The switching between the two expression is done in a vectorizable
C fashion.
C    This algorithm may be used for any functional form of
C F ( phi ).  One just has to replace the code as indicated below.
C
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'pick.inc'
      DOUBLE PRECISION EP,EPV
      DOUBLE PRECISION WGHT,VWGHT
      INTEGER IP(*), JP(*), KP(*), LP(*), NPHI, INTERE(*)
      DOUBLE PRECISION CPC(*)
      INTEGER CPD(*)
      DOUBLE PRECISION CPB(*)
      CHARACTER*(*) MODE
      DOUBLE PRECISION CPO(*)
      INTEGER CPE(*)
      CHARACTER*4 ANALYS
C local
C difference vectors
      DOUBLE PRECISION RIJX, RIJY, RIJZ, RJKX, RJKY, RJKZ
      DOUBLE PRECISION RKLX, RKLY, RKLZ
C >>MN modification begin
      DOUBLE PRECISION RIJ, RJK, RKL
C <<MN modification end
C coordinate basis vectors (cross products)
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION RAR, RBR, RCR
C derivatives
      DOUBLE PRECISION DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ
      DOUBLE PRECISION DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ
C
C phi, sin ( phi ) reciprocals, etc.
      DOUBLE PRECISION CP, SP, PHI, RECCP, RECSP
C switch
      INTEGER SWITCH
      DOUBLE PRECISION EPS
C >>MN modification begin
      DOUBLE PRECISION EPS2, EPS3, WARNG
C <<MN modification end
C energy and force of function
      DOUBLE PRECISION E, DF, DI, TEMP
C indices
      INTEGER IPHI, ILOC, IL, II
      LOGICAL GENRAL
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, MCONST, RADI
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
      PARAMETER (RADI=180.0D0/PI)
C===>parameter MCONST truncates the reciprocals and makes sure
C===>we don't devide by zero
      PARAMETER (MCONST=0.0001D0)
C===>parameter EPS determines the absolute value of SIN(PHI) at
C===>which we want to switch to the alternative computation for the
C===>derivative
      PARAMETER (EPS=0.1D0)
C >>MN modification begin
C EPS2 and EPS3 determine when to switch off the dihedral energy
C because a bond or angle are close to zero or 180
      PARAMETER (EPS2=0.99D0, EPS3=0.01D0)
C <<MN modification end
C local vector arrays
C===>VECDIM determines the optimal length of the vector loops
C===>on the CRAY this is 64
      INTEGER VECDIM
      PARAMETER (VECDIM=64)
      INTEGER LIND(VECDIM)
      DOUBLE PRECISION XI(VECDIM), XJ(VECDIM), XK(VECDIM), XL(VECDIM)
      DOUBLE PRECISION YI(VECDIM), YJ(VECDIM), YK(VECDIM), YL(VECDIM)
      DOUBLE PRECISION ZI(VECDIM), ZJ(VECDIM), ZK(VECDIM), ZL(VECDIM)
      DOUBLE PRECISION LCPC(VECDIM), LCPB(VECDIM), LCPO(VECDIM)
      INTEGER LCPD(VECDIM), LCPE(VECDIM)
      DOUBLE PRECISION DPRIJX(VECDIM), DPRIJY(VECDIM), DPRIJZ(VECDIM)
      DOUBLE PRECISION DPRJKX(VECDIM), DPRJKY(VECDIM), DPRJKZ(VECDIM)
      DOUBLE PRECISION DPRKLX(VECDIM), DPRKLY(VECDIM), DPRKLZ(VECDIM)
C
      LOGICAL QANAL
C
C begin
C
C initialization
      EP=0.0D0
C
      QANAL=.FALSE.
      IF (ANALYS.EQ.'ANAL') QANAL=.TRUE.
C
C ====================================================================
C sinusoidal form of the dihedral angle potential
C ====================================================================
C
C loop over all active torsions
C make vector loops of length VECDIM
      IPHI=0
      DO WHILE (IPHI.LT.NPHI)
      ILOC=0
      DO WHILE (ILOC.LT.VECDIM.AND.IPHI.LT.NPHI)
      IPHI=IPHI+1
      IF (CPD(IPHI).NE.0) THEN
      IF (QANAL.OR.ABS(INTERE(IP(IPHI))+INTERE(JP(IPHI))
     &    +INTERE(KP(IPHI))+INTERE(LP(IPHI))).LE.+3) THEN
      ILOC=ILOC+1
C
C Gather the data into temporary arrays suitable for vectorization.
C The data that have to be gathered are: index of current torsion
C (for scattering the forces later), coordinates of RI, RJ, RK, RL
C and the parameters CPD, CPC, and CPB for the function F ( PHI )
      LIND(ILOC)=IPHI
      END IF
      END IF
      END DO
      DO IL=1,ILOC
      XI(IL)=X(IP(LIND(IL)))
      XJ(IL)=X(JP(LIND(IL)))
      XK(IL)=X(KP(LIND(IL)))
      XL(IL)=X(LP(LIND(IL)))
      YI(IL)=Y(IP(LIND(IL)))
      YJ(IL)=Y(JP(LIND(IL)))
      YK(IL)=Y(KP(LIND(IL)))
      YL(IL)=Y(LP(LIND(IL)))
      ZI(IL)=Z(IP(LIND(IL)))
      ZJ(IL)=Z(JP(LIND(IL)))
      ZK(IL)=Z(KP(LIND(IL)))
      ZL(IL)=Z(LP(LIND(IL)))
      LCPC(IL)=CPC(LIND(IL))
      LCPD(IL)=CPD(LIND(IL))
      LCPB(IL)=CPB(LIND(IL))
      END DO
C
C =================
C-BEGIN-VECTOR-LOOP
C =================
      DO IL=1,ILOC
C
C first compute differences RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
      RIJX=XI(IL)-XJ(IL)
      RIJY=YI(IL)-YJ(IL)
      RIJZ=ZI(IL)-ZJ(IL)
      RJKX=XJ(IL)-XK(IL)
      RJKY=YJ(IL)-YK(IL)
      RJKZ=ZJ(IL)-ZK(IL)
      RKLX=XK(IL)-XL(IL)
      RKLY=YK(IL)-YL(IL)
      RKLZ=ZK(IL)-ZL(IL)
C
C now compute A= RIJ X RJK, B= RJK X RKL, C= RJK X ( RIJ X RJK )
      AX=RIJY*RJKZ-RIJZ*RJKY
      AY=RIJZ*RJKX-RIJX*RJKZ
      AZ=RIJX*RJKY-RIJY*RJKX
      BX=RJKY*RKLZ-RKLY*RJKZ
      BY=RJKZ*RKLX-RKLZ*RJKX
      BZ=RJKX*RKLY-RKLX*RJKY
      CX=RJKY*AZ-RJKZ*AY
      CY=RJKZ*AX-RJKX*AZ
      CZ=RJKX*AY-RJKY*AX
C
C compute the norm of A, B, C and set to MCONST if too small
C than MCONST
      RAR=ONE/SQRT(MAX(MCONST,AX*AX+AY*AY+AZ*AZ))
      RBR=ONE/SQRT(MAX(MCONST,BX*BX+BY*BY+BZ*BZ))
      RCR=ONE/SQRT(MAX(MCONST,CX*CX+CY*CY+CZ*CZ))
C
C normalize A, B and C
      AX=AX*RAR
      AY=AY*RAR
      AZ=AZ*RAR
      BX=BX*RBR
      BY=BY*RBR
      BZ=BZ*RBR
      CX=CX*RCR
      CY=CY*RCR
      CZ=CZ*RCR
C
C compute CP= cos ( phi )  and SP= sin ( phi )
      CP=AX*BX+AY*BY+AZ*BZ
      SP=CX*BX+CY*BY+CZ*BZ
C
C compute PHI ( make sure CP within boundaries and get sign from SP )
      PHI=-SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
C
C=====================================================================
C==>compute function(phi) and derivative
C==>replace this part if you want to use another functional
C==>form
C==>"proper" function:  const * cos( n phi + delta )
      E=LCPC(IL)*(ONE+COS(LCPD(IL)*PHI+LCPB(IL)))
      DF=-LCPD(IL)*LCPC(IL)*SIN(LCPD(IL)*PHI+LCPB(IL))*WGHT
C=====================================================================
C
C accumulate energy
      EP=EP+E
C
C compute heavyside function:
C
C               /   1   |SP| > EPS
C  SWITCH(SP) = |
C               \   0   |SP| < EPS
C this function allows us to switch between two alternatives to compute
C d phi/ d r
C
      SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
      RECSP=DF*SWITCH*SIGN(ONE,SP)/MAX(ABS(SP),MCONST)
      RECCP=DF*(SWITCH+1)*SIGN(ONE,CP)/MAX(ABS(CP),MCONST)
C
C compute:
C  d cos( PHI )  d cos ( PHI )
C  ------------, -------------
C     d A             d B
C
      DCPAX=-RAR*(BX-CP*AX)
      DCPAY=-RAR*(BY-CP*AY)
      DCPAZ=-RAR*(BZ-CP*AZ)
      DCPBX=-RBR*(AX-CP*BX)
      DCPBY=-RBR*(AY-CP*BY)
      DCPBZ=-RBR*(AZ-CP*BZ)
C
C compute:
C  d sin( PHI )  d sin ( PHI )
C  ------------, -------------
C     d C             d B
      DSPCX=-RCR*(BX-SP*CX)
      DSPCY=-RCR*(BY-SP*CY)
      DSPCZ=-RCR*(BZ-SP*CZ)
      DSPBX=-RBR*(CX-SP*BX)
      DSPBY=-RBR*(CY-SP*BY)
      DSPBZ=-RBR*(CZ-SP*BZ)
C
C compute:
C    d (phi)
C DF -------,  etc. (get rid of singularity by using two alternatives)
C    d rij
      DPRIJX(IL)=
     & RECSP*(RJKY*DCPAZ-DCPAY*RJKZ)+
     & RECCP*((RJKY**2+RJKZ**2)*DSPCX-RJKX*RJKY*DSPCY-RJKX*RJKZ*DSPCZ)
      DPRIJY(IL)=
     & RECSP*(RJKZ*DCPAX-DCPAZ*RJKX)+
     & RECCP*((RJKZ**2+RJKX**2)*DSPCY-RJKY*RJKZ*DSPCZ-RJKY*RJKX*DSPCX)
      DPRIJZ(IL)=
     & RECSP*(RJKX*DCPAY-DCPAX*RJKY)+
     & RECCP*((RJKX**2+RJKY**2)*DSPCZ-RJKZ*RJKX*DSPCX-RJKZ*RJKY*DSPCY)
C
      DPRKLX(IL)=
     & RECSP*(DCPBY*RJKZ-RJKY*DCPBZ)+
     & RECCP*(DSPBY*RJKZ-RJKY*DSPBZ)
      DPRKLY(IL)=
     & RECSP*(DCPBZ*RJKX-RJKZ*DCPBX)+
     & RECCP*(DSPBZ*RJKX-RJKZ*DSPBX)
      DPRKLZ(IL)=
     & RECSP*(DCPBX*RJKY-RJKX*DCPBY)+
     & RECCP*(DSPBX*RJKY-RJKX*DSPBY)
C
      DPRJKX(IL)=
     & RECSP*(DCPAY*RIJZ-DCPAZ*RIJY+DCPBZ*RKLY-DCPBY*RKLZ)+
     & RECCP*(-(RJKY*RIJY+RJKZ*RIJZ)*DSPCX
     &        +(TWO*RJKX*RIJY-RIJX*RJKY)*DSPCY
     &        +(TWO*RJKX*RIJZ-RIJX*RJKZ)*DSPCZ
     &        +DSPBZ*RKLY-DSPBY*RKLZ)
      DPRJKY(IL)=
     & RECSP*(DCPAZ*RIJX-DCPAX*RIJZ+DCPBX*RKLZ-DCPBZ*RKLX)+
     & RECCP*(-(RJKZ*RIJZ+RJKX*RIJX)*DSPCY
     &        +(TWO*RJKY*RIJZ-RIJY*RJKZ)*DSPCZ
     &        +(TWO*RJKY*RIJX-RIJY*RJKX)*DSPCX
     &        +DSPBX*RKLZ-DSPBZ*RKLX)
      DPRJKZ(IL)=
     & RECSP*(DCPAX*RIJY-DCPAY*RIJX+DCPBY*RKLX-DCPBX*RKLY)+
     & RECCP*(-(RJKX*RIJX+RJKY*RIJY)*DSPCZ
     &        +(TWO*RJKZ*RIJX-RIJZ*RJKX)*DSPCX
     &        +(TWO*RJKZ*RIJY-RIJZ*RJKY)*DSPCY
     &        +DSPBY*RKLX-DSPBX*RKLY)
      END DO
C ===============
C-END-VECTOR-LOOP
C ===============
C
C now rescatter forces !!!! ( non-vectorizable )
      DO IL=1,ILOC
      II=LIND(IL)
      DX(IP(II))=DX(IP(II))+DPRIJX(IL)
      DY(IP(II))=DY(IP(II))+DPRIJY(IL)
      DZ(IP(II))=DZ(IP(II))+DPRIJZ(IL)
      DX(JP(II))=DX(JP(II))+DPRJKX(IL)-DPRIJX(IL)
      DY(JP(II))=DY(JP(II))+DPRJKY(IL)-DPRIJY(IL)
      DZ(JP(II))=DZ(JP(II))+DPRJKZ(IL)-DPRIJZ(IL)
      DX(KP(II))=DX(KP(II))+DPRKLX(IL)-DPRJKX(IL)
      DY(KP(II))=DY(KP(II))+DPRKLY(IL)-DPRJKY(IL)
      DZ(KP(II))=DZ(KP(II))+DPRKLZ(IL)-DPRJKZ(IL)
      DX(LP(II))=DX(LP(II))           -DPRKLX(IL)
      DY(LP(II))=DY(LP(II))           -DPRKLY(IL)
      DZ(LP(II))=DZ(LP(II))           -DPRKLZ(IL)
      END DO
C
      END DO
C
C
C ====================================================================
C harmonic form of the dihedral angle potential
C ====================================================================
C
C the mode is for the harmonic form of the potential
      GENRAL=MODE.EQ.'GENERAL'
C
C loop over all active torsions
C make vector loops of length VECDIM
      IPHI=0
      DO WHILE (IPHI.LT.NPHI)
      ILOC=0
      DO WHILE (ILOC.LT.VECDIM.AND.IPHI.LT.NPHI)
      IPHI=IPHI+1
      IF (CPD(IPHI).EQ.0) THEN
      IF (QANAL.OR.ABS(INTERE(IP(IPHI))+INTERE(JP(IPHI))
     &    +INTERE(KP(IPHI))+INTERE(LP(IPHI))).LE.+3) THEN
      ILOC=ILOC+1
C
C Gather the data into temporary arrays suitable for vectorization.
C The data that have to be gathered are: index of current torsion
C (for scattering the forces later), coordinates of RI, RJ, RK, RL
C and the parameters CPD, CPC, and CPB for the function F ( PHI )
      LIND(ILOC)=IPHI
      END IF
      END IF
      END DO
      DO IL=1,ILOC
      XI(IL)=X(IP(LIND(IL)))
      XJ(IL)=X(JP(LIND(IL)))
      XK(IL)=X(KP(LIND(IL)))
      XL(IL)=X(LP(LIND(IL)))
      YI(IL)=Y(IP(LIND(IL)))
      YJ(IL)=Y(JP(LIND(IL)))
      YK(IL)=Y(KP(LIND(IL)))
      YL(IL)=Y(LP(LIND(IL)))
      ZI(IL)=Z(IP(LIND(IL)))
      ZJ(IL)=Z(JP(LIND(IL)))
      ZK(IL)=Z(KP(LIND(IL)))
      ZL(IL)=Z(LP(LIND(IL)))
      LCPC(IL)=CPC(LIND(IL))
      LCPD(IL)=CPD(LIND(IL))
      LCPB(IL)=CPB(LIND(IL))
      END DO
C
C decide whether to use the general form of the potential
      IF (GENRAL) THEN
      DO IL=1,ILOC
      LCPO(IL)=CPO(LIND(IL))
      LCPE(IL)=CPE(LIND(IL))
      END DO
      ELSE
      DO IL=1,ILOC
      LCPO(IL)=ZERO
      LCPE(IL)=TWO
      END DO
      END IF
C
C =================
C-BEGIN-VECTOR-LOOP
C =================
      DO IL=1,ILOC
C
C first compute differences RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
      RIJX=XI(IL)-XJ(IL)
      RIJY=YI(IL)-YJ(IL)
      RIJZ=ZI(IL)-ZJ(IL)
      RJKX=XJ(IL)-XK(IL)
      RJKY=YJ(IL)-YK(IL)
      RJKZ=ZJ(IL)-ZK(IL)
      RKLX=XK(IL)-XL(IL)
      RKLY=YK(IL)-YL(IL)
      RKLZ=ZK(IL)-ZL(IL)
C >>MN modification begin
C check if the dihedral angle is well defined:
C compute norm of vectors, set warning flag if too small
C compute truncated CP = cos( phi ), set the warning flag if there are
C linear angles
      WARNG=ZERO
      CP=RIJX**2+RIJY**2+RIJZ**2
      WARNG=WARNG+MAX(MCONST-ABS(CP),ZERO)
      CP=RJKX**2+RJKY**2+RJKZ**2
      WARNG=WARNG+MAX(MCONST-ABS(CP),ZERO)
      CP=RKLX**2+RKLY**2+RKLZ**2
      WARNG=WARNG+MAX(MCONST-ABS(CP),ZERO)
C
      RIJ=ONE/SQRT(MAX(MCONST,RIJX**2+RIJY**2+RIJZ**2))
      RJK=ONE/SQRT(MAX(MCONST,RJKX**2+RJKY**2+RJKZ**2))
      RKL=ONE/SQRT(MAX(MCONST,RKLX**2+RKLY**2+RKLZ**2))
      CP=(RIJX*RJKX+RIJY*RJKY+RIJZ*RJKZ)*RIJ*RJK
      WARNG=WARNG+MAX(ABS(CP)-EPS2,ZERO)
      CP=(RJKX*RKLX+RJKY*RKLY+RJKZ*RKLZ)*RJK*RKL
      WARNG=WARNG+MAX(ABS(CP)-EPS2,ZERO)
C <<MN modification end
C
C now compute A= RIJ X RJK, B= RJK X RKL, C= RJK X ( RIJ X RJK )
      AX=RIJY*RJKZ-RIJZ*RJKY
      AY=RIJZ*RJKX-RIJX*RJKZ
      AZ=RIJX*RJKY-RIJY*RJKX
      BX=RJKY*RKLZ-RKLY*RJKZ
      BY=RJKZ*RKLX-RKLZ*RJKX
      BZ=RJKX*RKLY-RKLX*RJKY
      CX=RJKY*AZ-RJKZ*AY
      CY=RJKZ*AX-RJKX*AZ
      CZ=RJKX*AY-RJKY*AX
C
C compute the norm of A, B, C and set to MCONST if two small
C than MCONST
      RAR=ONE/SQRT(MAX(MCONST,AX*AX+AY*AY+AZ*AZ))
      RBR=ONE/SQRT(MAX(MCONST,BX*BX+BY*BY+BZ*BZ))
      RCR=ONE/SQRT(MAX(MCONST,CX*CX+CY*CY+CZ*CZ))
C
C normalize A, B and C
      AX=AX*RAR
      AY=AY*RAR
      AZ=AZ*RAR
      BX=BX*RBR
      BY=BY*RBR
      BZ=BZ*RBR
      CX=CX*RCR
      CY=CY*RCR
      CZ=CZ*RCR
C
C compute CP= cos ( phi )  and SP= sin ( phi )
      CP=AX*BX+AY*BY+AZ*BZ
      SP=CX*BX+CY*BY+CZ*BZ
C
C compute PHI ( make sure CP within boundaries and get sign from SP )
      PHI=-SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
C
C=====================================================================
C==>compute function(phi) and derivative
C==>replace this part if you want to use another functional
C==>form
C "improper" function:
C   reduced form:  const * ( phi - phi0 ) **2
C   general form:  const * well(phi-phi0,offset) **exponent
C   where
C                     well(a,b) =    a-b   a > b
C                     well(a,b) =     0    -b < a < b
C                     well(a,b) =    a+b   a < -b
C
C
C take the difference between phi and phi0
      DI=PHI-LCPB(IL)
C >>MN modification begin
      IF (WARNG.GT.RSMALL) DI=ZERO
C <<MN modification end
C the difference phi-phi0 is always mapped between -pi and +pi
      DI=MOD(DI+THREE*PI,TWO*PI)-PI
C include the "offset" LCPO
      DI=MAX(ZERO,DI-LCPO(IL))+MIN(ZERO,DI+LCPO(IL))
C compute energy and first derivative
      DF=LCPC(IL)*DI**(LCPE(IL)-1)
      E=DF*DI
      DF=LCPE(IL)*DF*WGHT
C=====================================================================
C
C accumulate energy
      EP=EP+E
C
C compute heavyside function:
C
C               /   1   |SP| > EPS
C  SWITCH(SP) = |
C               \   0   |SP| < EPS
C this function allows us to switch between two alternatives to compute
C d phi/ d r
C
      SWITCH=-MIN(1,INT(ABS(SP)-EPS+1.0D0))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
      RECSP=DF*SWITCH*SIGN(ONE,SP)/MAX(ABS(SP),MCONST)
      RECCP=DF*(SWITCH+1)*SIGN(ONE,CP)/MAX(ABS(CP),MCONST)
C
C compute:
C  d cos( PHI )  d cos ( PHI )
C  ------------, -------------
C     d A             d B
C
      DCPAX=-RAR*(BX-CP*AX)
      DCPAY=-RAR*(BY-CP*AY)
      DCPAZ=-RAR*(BZ-CP*AZ)
      DCPBX=-RBR*(AX-CP*BX)
      DCPBY=-RBR*(AY-CP*BY)
      DCPBZ=-RBR*(AZ-CP*BZ)
C
C compute:
C  d sin( PHI )  d sin ( PHI )
C  ------------, -------------
C     d C             d B
      DSPCX=-RCR*(BX-SP*CX)
      DSPCY=-RCR*(BY-SP*CY)
      DSPCZ=-RCR*(BZ-SP*CZ)
      DSPBX=-RBR*(CX-SP*BX)
      DSPBY=-RBR*(CY-SP*BY)
      DSPBZ=-RBR*(CZ-SP*BZ)
C
C compute:
C    d (phi)
C DF -------,  etc. (get rid of singularity by using two alternatives)
C    d rij
      DPRIJX(IL)=
     & RECSP*(RJKY*DCPAZ-DCPAY*RJKZ)+
     & RECCP*((RJKY**2+RJKZ**2)*DSPCX-RJKX*RJKY*DSPCY-RJKX*RJKZ*DSPCZ)
      DPRIJY(IL)=
     & RECSP*(RJKZ*DCPAX-DCPAZ*RJKX)+
     & RECCP*((RJKZ**2+RJKX**2)*DSPCY-RJKY*RJKZ*DSPCZ-RJKY*RJKX*DSPCX)
      DPRIJZ(IL)=
     & RECSP*(RJKX*DCPAY-DCPAX*RJKY)+
     & RECCP*((RJKX**2+RJKY**2)*DSPCZ-RJKZ*RJKX*DSPCX-RJKZ*RJKY*DSPCY)
C
      DPRKLX(IL)=
     & RECSP*(DCPBY*RJKZ-RJKY*DCPBZ)+
     & RECCP*(DSPBY*RJKZ-RJKY*DSPBZ)
      DPRKLY(IL)=
     & RECSP*(DCPBZ*RJKX-RJKZ*DCPBX)+
     & RECCP*(DSPBZ*RJKX-RJKZ*DSPBX)
      DPRKLZ(IL)=
     & RECSP*(DCPBX*RJKY-RJKX*DCPBY)+
     & RECCP*(DSPBX*RJKY-RJKX*DSPBY)
C
      DPRJKX(IL)=
     & RECSP*(DCPAY*RIJZ-DCPAZ*RIJY+DCPBZ*RKLY-DCPBY*RKLZ)+
     & RECCP*(-(RJKY*RIJY+RJKZ*RIJZ)*DSPCX
     &        +(TWO*RJKX*RIJY-RIJX*RJKY)*DSPCY
     &        +(TWO*RJKX*RIJZ-RIJX*RJKZ)*DSPCZ
     &        +DSPBZ*RKLY-DSPBY*RKLZ)
      DPRJKY(IL)=
     & RECSP*(DCPAZ*RIJX-DCPAX*RIJZ+DCPBX*RKLZ-DCPBZ*RKLX)+
     & RECCP*(-(RJKZ*RIJZ+RJKX*RIJX)*DSPCY
     &        +(TWO*RJKY*RIJZ-RIJY*RJKZ)*DSPCZ
     &        +(TWO*RJKY*RIJX-RIJY*RJKX)*DSPCX
     &        +DSPBX*RKLZ-DSPBZ*RKLX)
      DPRJKZ(IL)=
     & RECSP*(DCPAX*RIJY-DCPAY*RIJX+DCPBY*RKLX-DCPBX*RKLY)+
     & RECCP*(-(RJKX*RIJX+RJKY*RIJY)*DSPCZ
     &        +(TWO*RJKZ*RIJX-RIJZ*RJKX)*DSPCX
     &        +(TWO*RJKZ*RIJY-RIJZ*RJKY)*DSPCY
     &        +DSPBY*RKLX-DSPBX*RKLY)
      END DO
C ===============
C-END-VECTOR-LOOP
C ===============
C
C now rescatter forces !!!! ( non-vectorizable )
      DO IL=1,ILOC
      II=LIND(IL)
      DX(IP(II))=DX(IP(II))+DPRIJX(IL)
      DY(IP(II))=DY(IP(II))+DPRIJY(IL)
      DZ(IP(II))=DZ(IP(II))+DPRIJZ(IL)
      DX(JP(II))=DX(JP(II))+DPRJKX(IL)-DPRIJX(IL)
      DY(JP(II))=DY(JP(II))+DPRJKY(IL)-DPRIJY(IL)
      DZ(JP(II))=DZ(JP(II))+DPRJKZ(IL)-DPRIJZ(IL)
      DX(KP(II))=DX(KP(II))+DPRKLX(IL)-DPRJKX(IL)
      DY(KP(II))=DY(KP(II))+DPRKLY(IL)-DPRJKY(IL)
      DZ(KP(II))=DZ(KP(II))+DPRKLZ(IL)-DPRJKZ(IL)
      DX(LP(II))=DX(LP(II))           -DPRKLX(IL)
      DY(LP(II))=DY(LP(II))           -DPRKLY(IL)
      DZ(LP(II))=DZ(LP(II))           -DPRKLZ(IL)
      END DO
C
      END DO
C
C weight the energies
      EPV=EP*VWGHT
      EP=EP*WGHT
C
C do some analysis if appropriate ( only for the last torsion !!! )
      IF (ANALYS.EQ.'ANAL') THEN
      IF (CPD(1).NE.0) THEN
      PCDATA(PCGEOM)=PHI*RADI
      PCDATA(PCENER)=EP
      PCDATA(PCCONS)=CPC(1)
      PCDATA(PCPERI)=CPD(1)
C===>ANALSYS FOR F(PHI)
      TEMP=(PHI*CPD(1)+CPB(1)-PI)/(2.0D0*PI)
      TEMP=INT(TEMP+0.5D0*SIGN(ONE,TEMP))
      PCDATA(PCEQUI)=(TEMP*2.0D0*PI-CPB(1)+PI)/CPD(1)*RADI
      PCDATA(PCDEVI)=PCDATA(PCEQUI)-PHI*RADI
      ELSE
      PCDATA(PCGEOM)=PHI*RADI
      PCDATA(PCENER)=E
      PCDATA(PCCONS)=CPC(1)
      PCDATA(PCPERI)=CPD(1)
C===>ANALSYS FOR F(PHI)
      PCDATA(PCDEVI)=-DI*RADI
      PCDATA(PCEQUI)=CPB(1)*RADI
      END IF
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,
     &                  QACPC,QACPB,QACPD,IAC,SEGID,RESID,TYPE)
C
C Gets type-based parameters for dihedrals.
C
C Authors: Axel T. Brunger and Thomas Simonson
C ============================================
C
      IMPLICIT NONE
      INCLUDE 'timer.inc'
      INCLUDE 'param.inc'
      INTEGER NPHI, IP(*), JP(*), KP(*), LP(*)
      DOUBLE PRECISION ACPC(*), ACPB(*)
      INTEGER ACPD(*)
      CHARACTER*4 IAC(*), SEGID(*), RESID(*), TYPE(*)
      LOGICAL QACPC(*), QACPD(*), QACPB(*)
C local
      INTEGER IFLAG, I, ICP, OFFSET
      LOGICAL PHIFLG, COND
      INTEGER MARK
      PARAMETER (MARK=0)
C
C begin
      IFLAG=0
      PHIFLG=.TRUE.
      OFFSET=1
      DO I=1,NPHI
C ignore atom-based parameter entries
      IF (QACPC(I).OR.QACPB(I).OR.QACPD(I)) THEN
C
C multiple dihedral condition
      COND=(I.GT.1)
      IF (COND) THEN
      COND=(IP(I-1).EQ.IP(I).AND.JP(I-1).EQ.JP(I).AND.
     1      KP(I-1).EQ.KP(I).AND.LP(I-1).EQ.LP(I))
      END IF
      IF (COND) THEN
      OFFSET=OFFSET+1
      IF (PHIFLG.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(2A)') ' CODDIH: using multiple potential terms for ',
     1     'some dihedrals'
      END IF
      PHIFLG=.FALSE.
      ELSE
      OFFSET=1
      ENDIF
      CALL BINDIH(ICP,IAC(IP(I)),IAC(JP(I)),IAC(KP(I)),IAC(LP(I)),
     &            NNCCP,KCP1,KCP2,KCP3,KCP4,OFFSET,MARK)
      IF (ICP.NE.MARK) THEN
      IF (QACPC(I)) ACPC(I)=CPC(ICP)
      IF (QACPB(I)) ACPB(I)=CPB(ICP)
      IF (QACPD(I)) ACPD(I)=CPD(ICP)
      ELSE
C dihedral wildcard lookup
      CALL BINDIH(ICP,'X   ',IAC(JP(I)),IAC(KP(I)),'X   ',
     &            NNCCP,KCP1,KCP2,KCP3,KCP4,OFFSET,MARK)
      IF (ICP.NE.MARK) THEN
      IF (QACPC(I)) ACPC(I)=CPC(ICP)
      IF (QACPB(I)) ACPB(I)=CPB(ICP)
      IF (QACPD(I)) ACPD(I)=CPD(ICP)
      ELSE
      WRITE(6,'(A)')
     &' %CODDIH-ERR: missing dihedral parameters %%%%%%%%%%%%%%%%%%%%%'
      IF (QACPC(I)) WRITE(6,'(A)') '  dihedral energy constant missing.'
      IF (QACPB(I)) WRITE(6,'(A)') '  target dihedral value missing.'
      IF (QACPD(I)) WRITE(6,'(A)') '  periodicity missing.'
      WRITE(6,'(9A)')
     & '  ATOM1: SEGId="',SEGID(IP(I)),'",  RESId="',RESID(IP(I)),
     & '",  NAME="',TYPE(IP(I)),'",  CHEMical="',IAC(IP(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM2: SEGId="',SEGID(JP(I)),'",  RESId="',RESID(JP(I)),
     & '",  NAME="',TYPE(JP(I)),'",  CHEMical="',IAC(JP(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM3: SEGId="',SEGID(KP(I)),'",  RESId="',RESID(KP(I)),
     & '",  NAME="',TYPE(KP(I)),'",  CHEMical="',IAC(KP(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM4: SEGId="',SEGID(LP(I)),'",  RESId="',RESID(LP(I)),
     & '",  NAME="',TYPE(LP(I)),'",  CHEMical="',IAC(LP(I)),'"'
      WRITE(6,'(A)')
     &' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      END IF
      END IF
      IF (ICP.EQ.MARK) IFLAG=IFLAG+1
      END IF
      END DO
C
      IF (IFLAG.GT.0) THEN
      CALL WRNDIE(-1,'CODDIH','program will be aborted.')
      ELSE IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)')
     &   ' CODDIH: dihedral type-based parameters retrieved'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE BINDIH(BINPAR,K1,K2,K3,K4,
     &                  NNCCP,KCP1,KCP2,KCP3,KCP4,OFFSET,MARK)
C
C makes binary searches to retrieve type-based bond
C parameters.
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INTEGER BINPAR
      CHARACTER*4 K1, K2, K3, K4
      INTEGER NNCCP
      CHARACTER*4 KCP1(*), KCP2(*), KCP3(*), KCP4(*)
      INTEGER OFFSET, MARK
C local
      LOGICAL LTSTEQ
      EXTERNAL LTSTEQ
      CHARACTER*16 A, B
      INTEGER   I, J, N
C begin
      BINPAR=MARK
      IF (NNCCP.GT.0) THEN
C dihedrals ( permutation 1,2,3,4 <-> 4,3,2,1 )
      IF (LTSTEQ(K2//K1,8,K3//K4,8,.TRUE.)) THEN
      A=K1//K2//K3//K4
      ELSE
      A=K4//K3//K2//K1
      END IF
      I=1
      J=NNCCP
      N=(I+J)/2
      B=KCP1(N)//KCP2(N)//KCP3(N)//KCP4(N)
      DO WHILE (B.NE.A.AND..NOT.I.GE.J)
      IF (.NOT.LTSTEQ(B,16,A,16,.TRUE.)) THEN
      J=N-1
      ELSE
      I=N+1
      END IF
      N=MAX(1,(I+J)/2)
      B=KCP1(N)//KCP2(N)//KCP3(N)//KCP4(N)
      END DO
      IF (A.EQ.B) THEN
C
C search for lowest index plus offset in dihedral list
      DO WHILE (N.GT.1.AND.B.EQ.A)
      N=N-1
      B=KCP1(N)//KCP2(N)//KCP3(N)//KCP4(N)
      END DO
      IF (.NOT.(N.EQ.1.AND.B.EQ.A)) N=N+1
      N=OFFSET-1+N
      B=KCP1(N)//KCP2(N)//KCP3(N)//KCP4(N)
      IF (B.EQ.A) THEN
      BINPAR=N
      END IF
      END IF
C
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,
     &                  QACIC,QACIB,QACID,IAC,SEGID,RESID,TYPE)
C
C Gets type-based parameters for impropers.
C
C Author: Axel T. Brunger and Thomas Simonson
C ===========================================
C
      IMPLICIT NONE
      INCLUDE 'timer.inc'
      INCLUDE 'param.inc'
      INTEGER NIMPHI, IM(*), JM(*), KM(*), LM(*)
      DOUBLE PRECISION ACIC(*), ACIB(*)
      INTEGER ACID(*)
      CHARACTER*4 IAC(*),SEGID(*),RESID(*),TYPE(*)
      LOGICAL QACIC(*), QACIB(*), QACID(*)
C local
      INTEGER IFLAG, I, ICI, OFFSET
      LOGICAL PHIFLG, COND
      INTEGER MARK
      PARAMETER (MARK=0)
C
C begin
      IFLAG=0
      PHIFLG=.TRUE.
      OFFSET=1
      DO I=1,NIMPHI
      IF (QACIC(I).OR.QACIB(I).OR.QACID(I)) THEN
C multiple improper condition
      COND=(I.GT.1)
      IF (COND) THEN
      COND=(IM(I-1).EQ.IM(I).AND.JM(I-1).EQ.JM(I).AND.
     &      KM(I-1).EQ.KM(I).AND.LM(I-1).EQ.LM(I))
      END IF
      IF (COND) THEN
      OFFSET=OFFSET+1
      IF (PHIFLG.AND.WRNLEV.GE.10) THEN
      WRITE(6,'(2A)') ' CODIMP: using multiple potential terms for ',
     &     'some impropers'
      END IF
      PHIFLG=.FALSE.
      ELSE
      OFFSET=1
      ENDIF
      CALL BINIMP(ICI,IAC(IM(I)),IAC(JM(I)),IAC(KM(I)),
     &       IAC(LM(I)),NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
      IF (ICI.NE.MARK) THEN
      IF (QACIC(I)) ACIC(I)=CIC(ICI)
      IF (QACIB(I)) ACIB(I)=CIB(ICI)
      IF (QACID(I)) ACID(I)=CID(ICI)
      ELSE
C improper wildcard lookup a-X-X-d
      CALL BINIMP(ICI,IAC(IM(I)),'X   ','X   ',IAC(LM(I)),
     &       NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
      IF (ICI.NE.MARK) THEN
      IF (QACIC(I)) ACIC(I)=CIC(ICI)
      IF (QACIB(I)) ACIB(I)=CIB(ICI)
      IF (QACID(I)) ACID(I)=CID(ICI)
      ELSE
C improper wildcard lookup X-b-c-d
      CALL BINIMP(ICI,'X   ',IAC(JM(I)),IAC(KM(I)),IAC(LM(I)),
     &       NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
      IF (ICI.NE.MARK) THEN
      IF (QACIC(I)) ACIC(I)=CIC(ICI)
      IF (QACIB(I)) ACIB(I)=CIB(ICI)
      IF (QACID(I)) ACID(I)=CID(ICI)
      ELSE
C improper wildcard lookup X-b-c-X
      CALL BINIMP(ICI,'X   ',IAC(JM(I)),IAC(KM(I)),'X   ',
     &       NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
      IF (ICI.NE.MARK) THEN
      IF (QACIC(I)) ACIC(I)=CIC(ICI)
      IF (QACIB(I)) ACIB(I)=CIB(ICI)
      IF (QACID(I)) ACID(I)=CID(ICI)
      ELSE
C improper wildcard lookup X-X-c-d
      CALL BINIMP(ICI,'X   ','X   ',IAC(KM(I)),IAC(LM(I)),
     &       NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
      IF (ICI.NE.MARK) THEN
      IF (QACIC(I)) ACIC(I)=CIC(ICI)
      IF (QACIB(I)) ACIB(I)=CIB(ICI)
      IF (QACID(I)) ACID(I)=CID(ICI)
      ELSE
      WRITE(6,'(A)')
     &' %CODIMP-ERR: missing improper parameters %%%%%%%%%%%%%%%%%%%%%'
      IF (QACIC(I)) WRITE(6,'(A)') '  improper energy constant missing.'
      IF (QACIB(I)) WRITE(6,'(A)') '  target improper value missing.'
      IF (QACID(I)) WRITE(6,'(A)') '  periodicity missing.'
      WRITE(6,'(9A)')
     & '  ATOM1: SEGId="',SEGID(IM(I)),'",  RESId="',RESID(IM(I)),
     & '",  NAME="',TYPE(IM(I)),'",  CHEMical="',IAC(IM(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM2: SEGId="',SEGID(JM(I)),'",  RESId="',RESID(JM(I)),
     & '",  NAME="',TYPE(JM(I)),'",  CHEMical="',IAC(JM(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM3: SEGId="',SEGID(KM(I)),'",  RESId="',RESID(KM(I)),
     & '",  NAME="',TYPE(KM(I)),'",  CHEMical="',IAC(KM(I)),'"'
      WRITE(6,'(9A)')
     & '  ATOM4: SEGId="',SEGID(LM(I)),'",  RESId="',RESID(LM(I)),
     & '",  NAME="',TYPE(LM(I)),'",  CHEMical="',IAC(LM(I)),'"'
      WRITE(6,'(A)')
     &' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      END IF
      END IF
      END IF
      END IF
      END IF
      IF (ICI.EQ.MARK) IFLAG=IFLAG+1
      END IF
      END DO
C
      IF (IFLAG.GT.0) THEN
      CALL WRNDIE(-1,'CODIMP','program will be aborted.')
      ELSE IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)')
     &   ' CODIMP: improper type-based parameters retrieved'
      END IF
C
      RETURN
      END
C======================================================================
      SUBROUTINE BINIMP(BINPAR,K1,K2,K3,K4,
     &                  NCI,KCI1,KCI2,KCI3,KCI4,OFFSET,MARK)
C
C makes binary searches to retrieve type-based bond
C parameters.
C Author: Axel T. Brunger
C =======================
C
      IMPLICIT NONE
C input/output
      INTEGER BINPAR
      CHARACTER*4 K1, K2, K3, K4
      INTEGER NCI
      CHARACTER*4 KCI1(*), KCI2(*), KCI3(*), KCI4(*)
      INTEGER OFFSET, MARK
C local
      LOGICAL LTSTEQ
      EXTERNAL LTSTEQ
      CHARACTER*16 A, B
      INTEGER   I, J, N
C begin
      BINPAR=MARK
      IF (NCI.GT.0) THEN
C impropers
      IF (LTSTEQ(K1//K2,8,K4//K3,8,.TRUE.)) THEN
      A=K1//K2//K3//K4
      ELSE
      A=K4//K3//K2//K1
      END IF
      I=1
      J=NCI
      N=(I+J)/2
      B=KCI1(N)//KCI2(N)//KCI3(N)//KCI4(N)
      DO WHILE (B.NE.A.AND..NOT.I.GE.J)
      IF (.NOT.LTSTEQ(B,16,A,16,.TRUE.)) THEN
      J=N-1
      ELSE
      I=N+1
      END IF
      N=MAX(1,(I+J)/2)
      B=KCI1(N)//KCI2(N)//KCI3(N)//KCI4(N)
      END DO
      IF (A.EQ.B) THEN
C
C search for lowest index plus offset in dihedral list
      DO WHILE (N.GT.1.AND.B.EQ.A)
      N=N-1
      B=KCI1(N)//KCI2(N)//KCI3(N)//KCI4(N)
      END DO
      IF (.NOT.(N.EQ.1.AND.B.EQ.A)) N=N+1
      N=OFFSET-1+N
      B=KCI1(N)//KCI2(N)//KCI3(N)//KCI4(N)
      IF (B.EQ.A) THEN
      BINPAR=N
      END IF
      END IF
C
      END IF
      RETURN
      END
