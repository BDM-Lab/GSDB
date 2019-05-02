C===============
      SUBROUTINE EANGLEDB (EDB, WHICH)
C
C Calls EANGLEDB, which does the actual energy calculation
C
C by John Kuszewski Aug 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION EDB
      CHARACTER*7 WHICH
C local vbls
      INTEGER CLASS
C begin
C
C zero out partial energy
C
      EDB = ZERO
C
C for each class, call the appropriate
C energy term based on its dimensionality
C
      DO CLASS = 1, NANGLEDBCLASSES
           IF (ANGLEDBCLASSTYPE(CLASS).EQ.ANGLE) THEN
                CALL EANGLEDBANG (EDB, HEAP(ANGLEDBIPTR),
     &               HEAP(ANGLEDBJPTR), HEAP(ANGLEDBKPTR),
     &               HEAP(ANGLEDBLPTR), HEAP(ANGLEDBMPTR),
     &               HEAP(ANGLEDBNPTR), HEAP(ANGLEDBPPTR),
     &               HEAP(ANGLEDBQPTR), HEAP(ANGLEDBRPTR),
     &               HEAP(ANGLEDBSPTR), HEAP(ANGLEDBTPTR),
     &               HEAP(EXPANGLEDBPTRS(CLASS)),
     &               CLASS, HEAP(CALCANGLEDBPTR),
     &               PHISTEPS(CLASS), PSISTEPS(CLASS), WHICH)
           ELSE IF (ANGLEDBCLASSTYPE(CLASS).EQ.TORSION) THEN
                CALL EANGLEDBTOR (EDB, HEAP(ANGLEDBIPTR),
     &               HEAP(ANGLEDBJPTR), HEAP(ANGLEDBKPTR),
     &               HEAP(ANGLEDBLPTR), HEAP(ANGLEDBMPTR),
     &               HEAP(ANGLEDBNPTR), HEAP(ANGLEDBPPTR),
     &               HEAP(ANGLEDBQPTR), HEAP(ANGLEDBRPTR),
     &               HEAP(ANGLEDBSPTR), HEAP(ANGLEDBTPTR),
     &               HEAP(ANGLEDBUPTR),
     &               HEAP(EXPANGLEDBPTRS(CLASS)),
     &               CLASS, HEAP(CALCANGLEDBPTR),
     &               PHISTEPS(CLASS), PSISTEPS(CLASS), WHICH)
           END IF
      END DO
      RETURN
      END
C===============
      SUBROUTINE EANGLEDBANG (EDB, ATOMI, ATOMJ, ATOMK, ATOML, ATOMM,
     &     ATOMN, ATOMP, ATOMQ, ATOMR, ATOMS, ATOMT,
     &     EXPECTEDANGLEDB, CLASS, CALCANGLEDB,
     &     CURPHIS, CURPSIS, WHICH)
C
C Calculates bond angle energies where the equilibrium
C value varies with the value of two dihedral angles.
C
C Based on Andy Karplus, "Experimentally observed conformation-
C dependent geometry and hidden strain in proteins"
C Protein Science (1996) 5:1406-1420
C
C energy is of the form
C
C    E = kf (theta - thetaExpected(phi,psi))^2
C
C WHICH is a flag that switches between energy & force and
C violations calcs
C
C by John Kuszewski June 1996
C
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*),
     &     ATOMN(*), ATOMP(*), ATOMQ(*), ATOMR(*), ATOMS(*), ATOMT(*),
     &     CURPHIS, CURPSIS, CLASS
      DOUBLE PRECISION EXPECTEDANGLEDB(CURPHIS, CURPSIS),
     &     CALCANGLEDB(6, *)
      DOUBLE PRECISION EDB
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, START, FINISH, LOWRPHI,
     &        HIRPHI, RPSI, LOWRPSI, HIRPSI, RPHI
      DOUBLE PRECISION ERR, K1
      DOUBLE PRECISION XI, XJ, XK, XL, XM, YI, YJ, YK, YL, YM, ZI
      DOUBLE PRECISION ZJ, ZK, ZL,ZM
      DOUBLE PRECISION XIJ, XJK, XKL, YIJ, YJK, YKL, ZIJ, ZJK
      DOUBLE PRECISION ZKL
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION AX2, AY2, AZ2, BX2, BY2, BZ2, CX2, CY2, CZ2
      DOUBLE PRECISION RAR, RBR, RCR, CP, SP, PHI
      DOUBLE PRECISION RAR2, RBR2, RCR2, CP2, SP2, PSI
      DOUBLE PRECISION E, DF, DF2
      DOUBLE PRECISION HIPHI, LOWPHI, HIPSI, LOWPSI
      DOUBLE PRECISION SWITCH, RECSP, RECCP
      DOUBLE PRECISION SWITCH2, RECSP2, RECCP2
      DOUBLE PRECISION DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ
      DOUBLE PRECISION DCPAX2, DCPAY2, DCPAZ2, DCPBX2, DCPBY2, DCPBZ2
      DOUBLE PRECISION DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ
      DOUBLE PRECISION DSPCX2, DSPCY2, DSPCZ2, DSPBX2, DSPBY2, DSPBZ2
      DOUBLE PRECISION DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ
      DOUBLE PRECISION DPRIJX2, DPRIJY2, DPRIJZ2, DPRKLX2, DPRKLY2
      DOUBLE PRECISION DPRKLZ2, PHICONV
      DOUBLE PRECISION DPRJKX, DPRJKY, DPRJKZ, PSICONV
      DOUBLE PRECISION DPRJKX2, DPRJKY2, DPRJKZ2
      DOUBLE PRECISION XN, XP, XQ, YN, YP, YQ, ZN, ZP, ZQ
      DOUBLE PRECISION XMN, XNP, XPQ, YMN, YNP, YPQ, ZMN, ZNP, ZPQ
      DOUBLE PRECISION CORRECT1, CORRECT2, CORRECT3, correct4, correct5
      DOUBLE PRECISION correct6, CORRECTPHI, CORRECTPSI
      DOUBLE PRECISION XR, YR, ZR, XS, YS, ZS, XT, YT, ZT
      DOUBLE PRECISION RIJX, RIJY, RIJZ
      DOUBLE PRECISION RKJX, RKJY, RKJZ, RIJ, RKJ, CT, ST, THETA
      DOUBLE PRECISION THETAEQUIL
      DOUBLE PRECISION DEDTHETA, DEDTHETAEQUIL, DPRKJX, DPRKJY, DPRKJZ
      DOUBLE PRECISION DTHETAEQUILDPHI, DTHETAEQUILDPSI, RECST
      LOGICAL USEDERIV
C
C begin
C
      PHICONV = CURPHIS / (TWO * PI)
      PSICONV = CURPSIS / (TWO * PI)
      K1 = ANGLEDBFORCES(CLASS)
      ERR = ANGLEDBERRORS(CLASS)
      CORRECT1 = ANGLEDBPHASEC(1, CLASS)
      CORRECT2 = ANGLEDBPHASEC(2, CLASS)
      CORRECT3 = ANGLEDBPHASEC(3, CLASS)
      CORRECT4 = ANGLEDBPHASEC(4, CLASS)
      CORRECT5 = ANGLEDBPHASEC(5, CLASS)
      CORRECT6 = ANGLEDBPHASEC(6, CLASS)
      USEDERIV = ANGLEDBDERIV(CLASS)
      IF (CLASS.EQ.1) THEN
           START = 1
      ELSE
           START = ANGLEDBASSNDX(CLASS-1)+1
      END IF
      FINISH = ANGLEDBASSNDX(CLASS)
C
C following Axel's code in ETOR,
C
      DO COUNT = START, FINISH
C
C part one:  get the equilibrium angle
C from the database
C
C get coords of each atom
C
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
          XN = X(ATOMN(COUNT))
          XP = X(ATOMP(COUNT))
          XQ = X(ATOMQ(COUNT))
          XR = X(ATOMR(COUNT))
          XS = X(ATOMS(COUNT))
          XT = X(ATOMT(COUNT))
C
C
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
          YN = Y(ATOMN(COUNT))
          YP = Y(ATOMP(COUNT))
          YQ = Y(ATOMQ(COUNT))
          YR = Y(ATOMR(COUNT))
          YS = Y(ATOMS(COUNT))
          YT = Y(ATOMT(COUNT))
C
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
          ZN = Z(ATOMN(COUNT))
          ZP = Z(ATOMP(COUNT))
          ZQ = Z(ATOMQ(COUNT))
          ZR = Z(ATOMR(COUNT))
          ZS = Z(ATOMS(COUNT))
          ZT = Z(ATOMT(COUNT))
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
          XIJ = XI - XJ
          XJK = XJ - XK
          XKL = XK - XL
          XMN = XM - XN
          XNP = XN - XP
          XPQ = XP - XQ
C
          YIJ = YI - YJ
          YJK = YJ - YK
          YKL = YK - YL
          YMN = YM - YN
          YNP = YN - YP
          YPQ = YP - YQ
C
          ZIJ = ZI - ZJ
          ZJK = ZJ - ZK
          ZKL = ZK - ZL
          ZMN = ZM - ZN
          ZNP = ZN - ZP
          ZPQ = ZP - ZQ
C
C now calculate A=RIJ*RJK, B = RJK*RKL, C = RJK*(RIJ*RJK)
C
          AX = YIJ*ZJK-ZIJ*YJK
          AY = ZIJ*XJK-XIJ*ZJK
          AZ = XIJ*YJK-YIJ*XJK
          BX = YJK*ZKL-YKL*ZJK
          BY = ZJK*XKL-ZKL*XJK
          BZ = XJK*YKL-XKL*YJK
          CX = YJK*AZ-ZJK*AY
          CY = ZJK*AX-XJK*AZ
          CZ = XJK*AY-YJK*AX
C
          AX2 = YMN*ZNP-ZMN*YNP
          AY2 = ZMN*XNP-XMN*ZNP
          AZ2 = XMN*YNP-YMN*XNP
          BX2 = YNP*ZPQ-YPQ*ZNP
          BY2 = ZNP*XPQ-ZPQ*XNP
          BZ2 = XNP*YPQ-XPQ*YNP
          CX2 = YNP*AZ2-ZNP*AY2
          CY2 = ZNP*AX2-XNP*AZ2
          CZ2 = XNP*AY2-YNP*AX2
C
C calculate the norm of A, B, & C & set to MCONST if it's too small
C
          RAR = ONE/SQRT(MAX(MCONST, AX**2+AY**2+AZ**2))
          RBR = ONE/SQRT(MAX(MCONST, BX**2+BY**2+BZ**2))
          RCR = ONE/SQRT(MAX(MCONST, CX**2+CY**2+CZ**2))
C
          RAR2 = ONE/SQRT(MAX(MCONST, AX2**2+AY2**2+AZ2**2))
          RBR2 = ONE/SQRT(MAX(MCONST, BX2**2+BY2**2+BZ2**2))
          RCR2 = ONE/SQRT(MAX(MCONST, CX2**2+CY2**2+CZ2**2))
C
C normalize A, B, & C
C
          AX = AX*RAR
          AY = AY*RAR
          AZ = AZ*RAR
          BX = BX*RBR
          BY = BY*RBR
          BZ = BZ*RBR
          CX = CX*RCR
          CY = CY*RCR
          CZ = CZ*RCR
C
          AX2 = AX2*RAR2
          AY2 = AY2*RAR2
          AZ2 = AZ2*RAR2
          BX2 = BX2*RBR2
          BY2 = BY2*RBR2
          BZ2 = BZ2*RBR2
          CX2 = CX2*RCR2
          CY2 = CY2*RCR2
          CZ2 = CZ2*RCR2
C
C calculate cos(phi) & sin(phi)
C
          CP = AX*BX+AY*BY+AZ*BZ
          SP = CX*BX+CY*BY+CZ*BZ
C
C calculate cos(psi) & sin(psi)
C
          CP2 = AX2*BX2+AY2*BY2+AZ2*BZ2
          SP2 = CX2*BX2+CY2*BY2+CZ2*BZ2
C
C calculate phi (make sure cos is within bounds and get sign from sin)
C and keep it it radians
C
          PHI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
          PSI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP2))),SP2)
C
C figure out which bin you should go into, given that the
C expectation values are in the range 1..PHISTEPS
C
C added the possibility of different phase corrections
C for different classes of PROCHECK data (eg., phi/psi
C vs. chi1/chi2 plots) 9/1/95
C
          if (phi.lt.CORRECT1) then
               CORRECTPHI = PHI + CORRECT2
          ELSE
               CORRECTPHI = PHI
          END IF
          RPHI = INT ((CORRECTPHI + CORRECT3) * PHICONV) + ONE
          IF (PSI.LT.CORRECT4) THEN
               CORRECTPSI = PSI + CORRECT5
          ELSE
               CORRECTPSI = PSI
          END IF
          RPSI = INT ((CORRECTPSI + CORRECT6) * PSICONV) + ONE
C
C enforce bounds
C
          IF (RPHI.GT.CURPHIS) THEN
               RPHI = RPHI - CURPHIS
          ELSE IF (RPHI.LT.1) THEN
               RPHI = RPHI + CURPHIS
          END IF
          IF (RPSI.GT.CURPSIS) THEN
               RPSI = RPSI - CURPSIS
          ELSE IF (RPSI.LT.1) THEN
               RPSI = RPSI + CURPSIS
          END IF
C
C look up expectation values
C
          THETAEQUIL = EXPECTEDANGLEDB(RPHI,RPSI)
C
C don't calculate energy & forces if you're at a spot without
C an expectation value
C
          IF (THETAEQUIL.GT.NOEXPTEST) THEN
               E = 0
               DF = 0
               DF2 = 0
C
C Calculate energy & forces
C
          ELSE
C
C part two:  get the current angle
C
C this is lifted from eangle2
C
C first compute the differences RIJ=RU-RV and RKJ=RW-RV
               RIJX=XR-XS
               RIJY=YR-YS
               RIJZ=ZR-ZS
               RKJX=XT-XS
               RKJY=YT-YS
               RKJZ=ZT-ZS
C
C compute the norm of RIJ, RKJ and set to MCONST if it is too small.
               RIJ=ONE/SQRT(MAX(MCONST,RIJX**2+RIJY**2+RIJZ**2))
               RKJ=ONE/SQRT(MAX(MCONST,RKJX**2+RKJY**2+RKJZ**2))
C
C normalize RIJ, RKJ
               RIJX=RIJX*RIJ
               RIJY=RIJY*RIJ
               RIJZ=RIJZ*RIJ
               RKJX=RKJX*RKJ
               RKJY=RKJY*RKJ
               RKJZ=RKJZ*RKJ
C
C compute CT = cos( theta )
               CT=RIJX*RKJX+RIJY*RKJY+RIJZ*RKJZ
C
C compute ST = sin( theta )   ( make sure CT within boundaries )
               ST=SQRT(MAX(ZERO,ONE-CT**2))
C
C compute THETA ( make sure  CT within boundaries )
               THETA=ACOS(MIN(ONE,MAX(-ONE,CT)))
C
C part three:  calculate energy and
C              simple parts of the force
C
               E = K1 * (THETA - THETAEQUIL)**2
               EDB=EDB+E
               DEDTHETA = 2 * K1 * (THETA - THETAEQUIL)
               DEDTHETAEQUIL = -2 * K1 * (THETA - THETAEQUIL)
C
C if we're in print mode, copy out the relevant data
C
          IF (WHICH.EQ.'ANALYZE') THEN
               CALCANGLEDB(1,COUNT) = THETAEQUIL
               CALCANGLEDB(2,COUNT) = THETA
               CALCANGLEDB(3,COUNT) = THETA - THETAEQUIL
               CALCANGLEDB(4,COUNT) = PHI
               CALCANGLEDB(5,COUNT) = PSI
               CALCANGLEDB(6,COUNT) = E
          END IF
C
C part four:  calculate d theta / d xi
C             this is also from EANGLE2
C
               DF = DEDTHETA
C
C compute reciprocal of SP multiplied by the derivative of the
C function DF
               RECST=DF/MAX(MCONST,ST)
C
C compute:
C    d (phi)
C DF -------,  etc.
C    d rij
               DPRIJX=RECST*RIJ*(CT*RIJX-RKJX)
               DPRIJY=RECST*RIJ*(CT*RIJY-RKJY)
               DPRIJZ=RECST*RIJ*(CT*RIJZ-RKJZ)
               DPRKJX=RECST*RKJ*(CT*RKJX-RIJX)
               DPRKJY=RECST*RKJ*(CT*RKJY-RIJY)
               DPRKJZ=RECST*RKJ*(CT*RKJZ-RIJZ)
C
C Since overall force eqn is dE/dtheta dtheta/dxi + ...,
C I can update the forces on the three atoms involved
C in calculating theta right now
C
               IF (WHICH.NE.'ANALYZE') THEN
                    DX(ATOMR(COUNT)) = DX(ATOMR(COUNT)) + DPRIJX
                    DY(ATOMR(COUNT)) = DY(ATOMR(COUNT)) + DPRIJY
                    DZ(ATOMR(COUNT)) = DZ(ATOMR(COUNT)) + DPRIJZ
                    DX(ATOMS(COUNT)) = DX(ATOMS(COUNT)) -
     &                   DPRIJX - DPRKJX
                    DY(ATOMS(COUNT)) = DY(ATOMS(COUNT)) -
     &                   DPRIJY - DPRKJY
                    DZ(ATOMS(COUNT)) = DZ(ATOMS(COUNT)) -
     &                   DPRIJZ - DPRKJZ
                    DX(ATOMT(COUNT)) = DX(ATOMT(COUNT)) + DPRKJX
                    DY(ATOMT(COUNT)) = DY(ATOMT(COUNT)) + DPRKJY
                    DZ(ATOMT(COUNT)) = DZ(ATOMT(COUNT)) + DPRKJZ
               END IF
C
C part 5: get d thetaequil / d phi and d thetaequil / d psi
C
C get slope along phi
C
C figure out what the neighboring squares are
C and get their energies from the grid
C
               LOWRPHI = RPHI - 1
               IF (LOWRPHI.LE.0) LOWRPHI = LOWRPHI + CURPHIS
               HIRPHI = RPHI + 1
               IF (HIRPHI.GT.CURPHIS) HIRPHI = HIRPHI-CURPHIS
               HIPHI = EXPECTEDANGLEDB(HIRPHI, RPSI)
               LOWPHI = EXPECTEDANGLEDB(LOWRPHI, RPSI)
C
C deal with the possibility of being on the edge of
C an unknown area
C
               IF ((HIPHI.GT.NOEXPTEST).AND.(LOWPHI.GT.NOEXPTEST)) THEN
                    DTHETAEQUILDPHI = ZERO
               ELSE IF (HIPHI.GT.NOEXPTEST) THEN
                    DTHETAEQUILDPHI = ((THETAEQUIL - LOWPHI) / ONE)
     &                   * PHICONV
               ELSE IF (LOWPHI.GT.NOEXPTEST) THEN
                    DTHETAEQUILDPHI = ((HIPHI - THETAEQUIL) / ONE)
     &                   * PHICONV
               ELSE
                    DTHETAEQUILDPHI = ((HIPHI - LOWPHI) / TWO)
     &                   * PHICONV
               END IF
C
C do the same thing for the derivative along psi
C
               LOWRPSI = RPSI - 1
               IF (LOWRPSI.LE.0) LOWRPSI = LOWRPSI + CURPSIS
               HIRPSI = RPSI + 1
               IF (HIRPSI.GT.CURPSIS) HIRPSI = HIRPSI-CURPSIS
               HIPSI = EXPECTEDANGLEDB(RPHI, HIRPSI)
               LOWPSI = EXPECTEDANGLEDB(RPHI, LOWRPSI)
               IF (ANGLEDBPOTENTIAL.EQ.SQUARE) THEN
                    HIPSI = MAX(HIPSI-ERR, 0.0D0)
                    LOWPSI = MAX(LOWPSI-ERR,0.0D0)
               END IF
               IF ((HIPSI.GT.NOEXPTEST).AND.(LOWPSI.GT.NOEXPTEST)) THEN
                    DTHETAEQUILDPSI = 0
               ELSE IF (HIPSI.GT.NOEXPTEST) THEN
                    DTHETAEQUILDPSI = ((THETAEQUIL - LOWPSI) / ONE)
     &                   * PSICONV
               ELSE IF (LOWPSI.GT.NOEXPTEST) THEN
                    DTHETAEQUILDPSI = ((HIPSI - THETAEQUIL) / ONE)
     &                   * PSICONV
               ELSE
                    DTHETAEQUILDPSI = ((HIPSI - LOWPSI) / TWO)
     &                   * PSICONV
               END IF
          END IF
C
C part six: calculate d phi / dxi, d psi / dxi etc
C and in the process calcualte dE / dxi etc.
C
C compute heavyside function
C
          SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
          SWITCH2=-MIN(1,INT(ABS(SP2)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
C
C set DF to dE/dPhi and DF2 to dE/dPsi
C
          DF = DEDTHETAEQUIL * DTHETAEQUILDPHI
          DF2 = DEDTHETAEQUIL * DTHETAEQUILDPSI
          RECSP=DF*SWITCH*SIGN(ONE,SP)/MAX(ABS(SP),MCONST)
          RECCP=DF*(SWITCH+1)*SIGN(ONE,CP)/MAX(ABS(CP),MCONST)
          RECSP2=DF2*SWITCH2*SIGN(ONE,SP2)/MAX(ABS(SP2),MCONST)
          RECCP2=DF2*(SWITCH2+1)*SIGN(ONE,CP2)/MAX(ABS(CP2),MCONST)
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
          DCPAX2=-RAR2*(BX2-CP2*AX2)
          DCPAY2=-RAR2*(BY2-CP2*AY2)
          DCPAZ2=-RAR2*(BZ2-CP2*AZ2)
          DCPBX2=-RBR2*(AX2-CP2*BX2)
          DCPBY2=-RBR2*(AY2-CP2*BY2)
          DCPBZ2=-RBR2*(AZ2-CP2*BZ2)
C
C compute:
C  d sin( PHI )  d sin ( PHI )
C  ------------, -------------
C     d C             d B
C
          DSPCX=-RCR*(BX-SP*CX)
          DSPCY=-RCR*(BY-SP*CY)
          DSPCZ=-RCR*(BZ-SP*CZ)
          DSPBX=-RBR*(CX-SP*BX)
          DSPBY=-RBR*(CY-SP*BY)
          DSPBZ=-RBR*(CZ-SP*BZ)
C
          DSPCX2=-RCR2*(BX2-SP2*CX2)
          DSPCY2=-RCR2*(BY2-SP2*CY2)
          DSPCZ2=-RCR2*(BZ2-SP2*CZ2)
          DSPBX2=-RBR2*(CX2-SP2*BX2)
          DSPBY2=-RBR2*(CY2-SP2*BY2)
          DSPBZ2=-RBR2*(CZ2-SP2*BZ2)
C
C compute:
C    d (phi)
C DF -------,  etc. (get rid of singularity by using two alternatives)
C    d rij
C
          DPRIJX=
     &     RECSP*(YJK*DCPAZ-DCPAY*ZJK)+
     &     RECCP*((YJK**2+ZJK**2)*DSPCX-XJK*YJK*DSPCY-XJK*ZJK*DSPCZ)
          DPRIJY=
     &     RECSP*(ZJK*DCPAX-DCPAZ*XJK)+
     &     RECCP*((ZJK**2+XJK**2)*DSPCY-YJK*ZJK*DSPCZ-YJK*XJK*DSPCX)
          DPRIJZ=
     &     RECSP*(XJK*DCPAY-DCPAX*YJK)+
     &     RECCP*((XJK**2+YJK**2)*DSPCZ-ZJK*XJK*DSPCX-ZJK*YJK*DSPCY)
C
          DPRKLX=
     &     RECSP*(DCPBY*ZJK-YJK*DCPBZ)+
     &     RECCP*(DSPBY*ZJK-YJK*DSPBZ)
          DPRKLY=
     &     RECSP*(DCPBZ*XJK-ZJK*DCPBX)+
     &     RECCP*(DSPBZ*XJK-ZJK*DSPBX)
          DPRKLZ=
     &     RECSP*(DCPBX*YJK-XJK*DCPBY)+
     &     RECCP*(DSPBX*YJK-XJK*DSPBY)
C
          DPRJKX=
     &     RECSP*(DCPAY*ZIJ-DCPAZ*YIJ+DCPBZ*YKL-DCPBY*ZKL)+
     &     RECCP*(-(YJK*YIJ+ZJK*ZIJ)*DSPCX
     &        +(TWO*XJK*YIJ-XIJ*YJK)*DSPCY
     &        +(TWO*XJK*ZIJ-XIJ*ZJK)*DSPCZ
     &        +DSPBZ*YKL-DSPBY*ZKL)
          DPRJKY=
     &     RECSP*(DCPAZ*XIJ-DCPAX*ZIJ+DCPBX*ZKL-DCPBZ*XKL)+
     &     RECCP*(-(ZJK*ZIJ+XJK*XIJ)*DSPCY
     &        +(TWO*YJK*ZIJ-YIJ*ZJK)*DSPCZ
     &        +(TWO*YJK*XIJ-YIJ*XJK)*DSPCX
     &        +DSPBX*ZKL-DSPBZ*XKL)
          DPRJKZ=
     &     RECSP*(DCPAX*YIJ-DCPAY*XIJ+DCPBY*XKL-DCPBX*YKL)+
     &     RECCP*(-(XJK*XIJ+YJK*YIJ)*DSPCZ
     &        +(TWO*ZJK*XIJ-ZIJ*XJK)*DSPCX
     &        +(TWO*ZJK*YIJ-ZIJ*YJK)*DSPCY
     &        +DSPBY*XKL-DSPBX*YKL)
C
C
          DPRIJX2=
     &     RECSP2*(YKL*DCPAZ2-DCPAY2*ZKL)+
     &     RECCP2*((YKL**2+ZKL**2)*DSPCX2-XKL*YKL*DSPCY2-
     &     XKL*ZKL*DSPCZ2)
          DPRIJY2=
     &     RECSP2*(ZKL*DCPAX2-DCPAZ2*XKL)+
     &     RECCP2*((ZKL**2+XKL**2)*DSPCY2-YKL*ZKL*DSPCZ2-
     &     YKL*XKL*DSPCX2)
          DPRIJZ2=
     &     RECSP2*(XKL*DCPAY2-DCPAX2*YKL)+
     &     RECCP2*((XKL**2+YKL**2)*DSPCZ2-ZKL*XKL*DSPCX2-
     &     ZKL*YKL*DSPCY2)
C
          DPRKLX2=
     &     RECSP2*(DCPBY2*ZKL-YKL*DCPBZ2)+
     &     RECCP2*(DSPBY2*ZKL-YKL*DSPBZ2)
          DPRKLY2=
     &     RECSP2*(DCPBZ2*XKL-ZKL*DCPBX2)+
     &     RECCP2*(DSPBZ2*XKL-ZKL*DSPBX2)
          DPRKLZ2=
     &     RECSP2*(DCPBX2*YKL-XKL*DCPBY2)+
     &     RECCP2*(DSPBX2*YKL-XKL*DSPBY2)
C
          DPRJKX2=
     &     RECSP2*(DCPAY2*ZJK-DCPAZ2*YJK+DCPBZ2*ypq-DCPBY2*zpq)+
     &     RECCP2*(-(YKL*YJK+ZKL*ZJK)*DSPCX2
     &        +(TWO*XKL*YJK-XJK*YKL)*DSPCY2
     &        +(TWO*XKL*ZJK-XJK*ZKL)*DSPCZ2
     &        +DSPBZ2*ypq-DSPBY2*zpq)
          DPRJKY2=
     &     RECSP2*(DCPAZ2*XJK-DCPAX2*ZJK+DCPBX2*zpq-DCPBZ2*xpq)+
     &     RECCP2*(-(ZKL*ZJK+XKL*XJK)*DSPCY2
     &        +(TWO*YKL*ZJK-YJK*ZKL)*DSPCZ2
     &        +(TWO*YKL*XJK-YJK*XKL)*DSPCX2
     &        +DSPBX2*zpq-DSPBZ2*xpq)
          DPRJKZ2=
     &     RECSP2*(DCPAX2*YJK-DCPAY2*XJK+DCPBY2*xpq-DCPBX2*ypq)+
     &     RECCP2*(-(XKL*XJK+YKL*YJK)*DSPCZ2
     &        +(TWO*ZKL*XJK-ZJK*XKL)*DSPCX2
     &        +(TWO*ZKL*YJK-ZJK*YKL)*DSPCY2
     &        +DSPBY2*xpq-DSPBX2*ypq)
C
C now update forces if in energy & force mode
C
          IF ((WHICH.NE.'ANALYZE').AND.(USEDERIV)) THEN
               DX(ATOMI(COUNT))=DX(ATOMI(COUNT))+DPRIJX
               DY(ATOMI(COUNT))=DY(ATOMI(COUNT))+DPRIJY
               DZ(ATOMI(COUNT))=DZ(ATOMI(COUNT))+DPRIJZ
               DX(ATOMJ(COUNT))=DX(ATOMJ(COUNT))+DPRJKX-DPRIJX
               DY(ATOMJ(COUNT))=DY(ATOMJ(COUNT))+DPRJKY-DPRIJY
               DZ(ATOMJ(COUNT))=DZ(ATOMJ(COUNT))+DPRJKZ-DPRIJZ
               DX(ATOMK(COUNT))=DX(ATOMK(COUNT))+DPRKLX-DPRJKX
               DY(ATOMK(COUNT))=DY(ATOMK(COUNT))+DPRKLY-DPRJKY
               DZ(ATOMK(COUNT))=DZ(ATOMK(COUNT))+DPRKLZ-DPRJKZ
               DX(ATOML(COUNT))=DX(ATOML(COUNT))       -DPRKLX
               DY(ATOML(COUNT))=DY(ATOML(COUNT))       -DPRKLY
               DZ(ATOML(COUNT))=DZ(ATOML(COUNT))       -DPRKLZ
C
               DX(ATOMM(COUNT))=DX(ATOMM(COUNT))+DPRIJX2
               DY(ATOMM(COUNT))=DY(ATOMM(COUNT))+DPRIJY2
               DZ(ATOMM(COUNT))=DZ(ATOMM(COUNT))+DPRIJZ2
               DX(ATOMN(COUNT))=DX(ATOMN(COUNT))+DPRJKX2-DPRIJX2
               DY(ATOMN(COUNT))=DY(ATOMN(COUNT))+DPRJKY2-DPRIJY2
               DZ(ATOMN(COUNT))=DZ(ATOMN(COUNT))+DPRJKZ2-DPRIJZ2
               DX(ATOMP(COUNT))=DX(ATOMP(COUNT))+DPRKLX2-DPRJKX2
               DY(ATOMP(COUNT))=DY(ATOMP(COUNT))+DPRKLY2-DPRJKY2
               DZ(ATOMP(COUNT))=DZ(ATOMP(COUNT))+DPRKLZ2-DPRJKZ2
               DX(ATOMQ(COUNT))=DX(ATOMQ(COUNT))        -DPRKLX2
               DY(ATOMQ(COUNT))=DY(ATOMQ(COUNT))        -DPRKLY2
               DZ(ATOMQ(COUNT))=DZ(ATOMQ(COUNT))        -DPRKLZ2
          END IF
      END DO
      RETURN
      END
C===============
      SUBROUTINE EANGLEDBTOR (EDB, ATOMI, ATOMJ, ATOMK, ATOML, ATOMM,
     &     ATOMN, ATOMP, ATOMQ, ATOMR, ATOMS, ATOMT, ATOMU,
     &     EXPECTEDANGLEDB, CLASS, CALCANGLEDB,
     &     CURPHIS, CURPSIS, WHICH)
C
C Calculates torsion angle energies where the equilibrium
C value varies with the value of two dihedral angles.
C
C Based on Andy Karplus, "Experimentally observed conformation-
C dependent geometry and hidden strain in proteins"
C Protein Science (1996) 5:1406-1420
C
C energy is of the form
C
C    E = kf (OMEGA - OMEGAExpected(phi,psi))^2
C
C WHICH is a flag that switches between energy & force and
C violations calcs
C
C by John Kuszewski June 1996
C
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*)
      INTEGER ATOMN(*), ATOMP(*), ATOMQ(*), ATOMR(*), ATOMS(*)
      INTEGER ATOMT(*), ATOMU(*), CURPHIS, CURPSIS, CLASS
      DOUBLE PRECISION EXPECTEDANGLEDB(CURPHIS, CURPSIS)
      DOUBLE PRECISION CALCANGLEDB(6, *)
      DOUBLE PRECISION EDB
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, START, FINISH, LOWRPHI
      INTEGER HIRPHI, RPSI, LOWRPSI, HIRPSI, RPHI
      DOUBLE PRECISION ER, ERR, K1
      DOUBLE PRECISION XI, XJ, XK, XL, XM, YI, YJ, YK, YL, YM, ZI
      DOUBLE PRECISION ZJ, ZK, ZL,ZM
      DOUBLE PRECISION XIJ, XJK, XKL, YIJ, YJK, YKL, ZIJ, ZJK
      DOUBLE PRECISION ZKL
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION AX2, AY2, AZ2, BX2, BY2, BZ2, CX2, CY2, CZ2
      DOUBLE PRECISION RAR, RBR, RCR, CP, SP, PHI
      DOUBLE PRECISION RAR2, RBR2, RCR2, CP2, SP2, PSI
      DOUBLE PRECISION E, DF, DF2
      DOUBLE PRECISION HIPHI, LOWPHI, HIPSI, LOWPSI
      DOUBLE PRECISION SWITCH, RECSP, RECCP
      DOUBLE PRECISION SWITCH2, RECSP2, RECCP2
      DOUBLE PRECISION DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ
      DOUBLE PRECISION DCPAX2, DCPAY2, DCPAZ2, DCPBX2, DCPBY2, DCPBZ2
      DOUBLE PRECISION DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ
      DOUBLE PRECISION DSPCX2, DSPCY2, DSPCZ2, DSPBX2, DSPBY2, DSPBZ2
      DOUBLE PRECISION DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ
      DOUBLE PRECISION DPRIJX2, DPRIJY2, DPRIJZ2, DPRKLX2, DPRKLY2
      DOUBLE PRECISION DPRKLZ2, PHICONV
      DOUBLE PRECISION DPRJKX, DPRJKY, DPRJKZ, PSICONV
      DOUBLE PRECISION DPRJKX2, DPRJKY2, DPRJKZ2
      DOUBLE PRECISION XN, XP, XQ, YN, YP, YQ, ZN, ZP, ZQ
      DOUBLE PRECISION XMN, XNP, XPQ, YMN, YNP, YPQ, ZMN, ZNP, ZPQ
      DOUBLE PRECISION CORRECT1, CORRECT2, CORRECT3, correct4
      DOUBLE PRECISION correct5
      DOUBLE PRECISION correct6, CORRECT7, CORRECT8, CORRECT9
      DOUBLE PRECISION CORRECTPHI, CORRECTPSI
      DOUBLE PRECISION XR, YR, ZR, XS, YS, ZS, XT, YT, ZT
      DOUBLE PRECISION XU, YU, ZU
      DOUBLE PRECISION OMEGAEQUIL
      DOUBLE PRECISION DEDOMEGA, DEDOMEGAEQUIL
      DOUBLE PRECISION DOMEGAEQUILDPHI, DOMEGAEQUILDPSI
C
      DOUBLE PRECISION XRS, XST, XTU, YRS, YST, YTU, ZRS, ZST, ZTU
      DOUBLE PRECISION AX3, AY3, AZ3, BX3, BY3, BZ3, CX3, CY3, CZ3
      DOUBLE PRECISION RAR3, RBR3, RCR3, CP3, SP3, OMEGA, CORRECTOMEGA
      DOUBLE PRECISION SWITCH3
      DOUBLE PRECISION RECSP3, RECCP3, DCPAX3, DCPAY3, DCPAZ3, DCPBX3
      DOUBLE PRECISION DCPBY3
      DOUBLE PRECISION DCPBZ3, DSPCX3, DSPCY3, DSPCZ3, DSPBX3, DSPBY3
      DOUBLE PRECISION DSPBZ3
      DOUBLE PRECISION DPRIJX3, DPRIJY3, DPRIJZ3, DPRKLX3, DPRKLY3
      DOUBLE PRECISION DPRKLZ3
      DOUBLE PRECISION DPRJKX3, DPRJKY3, DPRJKZ3
      LOGICAL USEDERIV
C
C begin
C
      PHICONV = CURPHIS / (TWO * PI)
      PSICONV = CURPSIS / (TWO * PI)
      K1 = ANGLEDBFORCES(CLASS)
      ERR = ANGLEDBERRORS(CLASS)
      CORRECT1 = ANGLEDBPHASEC(1, CLASS)
      CORRECT2 = ANGLEDBPHASEC(2, CLASS)
      CORRECT3 = ANGLEDBPHASEC(3, CLASS)
      CORRECT4 = ANGLEDBPHASEC(4, CLASS)
      CORRECT5 = ANGLEDBPHASEC(5, CLASS)
      CORRECT6 = ANGLEDBPHASEC(6, CLASS)
      CORRECT7 = ANGLEDBPHASEC(7, CLASS)
      CORRECT8 = ANGLEDBPHASEC(8, CLASS)
      CORRECT9 = ANGLEDBPHASEC(9, CLASS)
      USEDERIV = ANGLEDBDERIV(CLASS)
      IF (CLASS.EQ.1) THEN
           START = 1
      ELSE
           START = ANGLEDBASSNDX(CLASS-1)+1
      END IF
      FINISH = ANGLEDBASSNDX(CLASS)
      ER = 0.0D0
C
C following Axel's code in ETOR,
C
      DO COUNT = START, FINISH
C
C part one:  get the equilibrium angle
C from the database
C
C get coords of each atom
C
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
          XN = X(ATOMN(COUNT))
          XP = X(ATOMP(COUNT))
          XQ = X(ATOMQ(COUNT))
          XR = X(ATOMR(COUNT))
          XS = X(ATOMS(COUNT))
          XT = X(ATOMT(COUNT))
          XU = X(ATOMU(COUNT))
C
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
          YN = Y(ATOMN(COUNT))
          YP = Y(ATOMP(COUNT))
          YQ = Y(ATOMQ(COUNT))
          YR = Y(ATOMR(COUNT))
          YS = Y(ATOMS(COUNT))
          YT = Y(ATOMT(COUNT))
          YU = Y(ATOMU(COUNT))
C
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
          ZN = Z(ATOMN(COUNT))
          ZP = Z(ATOMP(COUNT))
          ZQ = Z(ATOMQ(COUNT))
          ZR = Z(ATOMR(COUNT))
          ZS = Z(ATOMS(COUNT))
          ZT = Z(ATOMT(COUNT))
          ZU = Z(ATOMU(COUNT))
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
          XIJ = XI - XJ
          XJK = XJ - XK
          XKL = XK - XL
          XMN = XM - XN
          XNP = XN - XP
          XPQ = XP - XQ
C
          YIJ = YI - YJ
          YJK = YJ - YK
          YKL = YK - YL
          YMN = YM - YN
          YNP = YN - YP
          YPQ = YP - YQ
C
          ZIJ = ZI - ZJ
          ZJK = ZJ - ZK
          ZKL = ZK - ZL
          ZMN = ZM - ZN
          ZNP = ZN - ZP
          ZPQ = ZP - ZQ
C
C now calculate A=RIJ*RJK, B = RJK*RKL, C = RJK*(RIJ*RJK)
C
          AX = YIJ*ZJK-ZIJ*YJK
          AY = ZIJ*XJK-XIJ*ZJK
          AZ = XIJ*YJK-YIJ*XJK
          BX = YJK*ZKL-YKL*ZJK
          BY = ZJK*XKL-ZKL*XJK
          BZ = XJK*YKL-XKL*YJK
          CX = YJK*AZ-ZJK*AY
          CY = ZJK*AX-XJK*AZ
          CZ = XJK*AY-YJK*AX
C
          AX2 = YMN*ZNP-ZMN*YNP
          AY2 = ZMN*XNP-XMN*ZNP
          AZ2 = XMN*YNP-YMN*XNP
          BX2 = YNP*ZPQ-YPQ*ZNP
          BY2 = ZNP*XPQ-ZPQ*XNP
          BZ2 = XNP*YPQ-XPQ*YNP
          CX2 = YNP*AZ2-ZNP*AY2
          CY2 = ZNP*AX2-XNP*AZ2
          CZ2 = XNP*AY2-YNP*AX2
C
C calculate the norm of A, B, & C & set to MCONST if it's too small
C
          RAR = ONE/SQRT(MAX(MCONST, AX**2+AY**2+AZ**2))
          RBR = ONE/SQRT(MAX(MCONST, BX**2+BY**2+BZ**2))
          RCR = ONE/SQRT(MAX(MCONST, CX**2+CY**2+CZ**2))
C
          RAR2 = ONE/SQRT(MAX(MCONST, AX2**2+AY2**2+AZ2**2))
          RBR2 = ONE/SQRT(MAX(MCONST, BX2**2+BY2**2+BZ2**2))
          RCR2 = ONE/SQRT(MAX(MCONST, CX2**2+CY2**2+CZ2**2))
C
C normalize A, B, & C
C
          AX = AX*RAR
          AY = AY*RAR
          AZ = AZ*RAR
          BX = BX*RBR
          BY = BY*RBR
          BZ = BZ*RBR
          CX = CX*RCR
          CY = CY*RCR
          CZ = CZ*RCR
C
          AX2 = AX2*RAR2
          AY2 = AY2*RAR2
          AZ2 = AZ2*RAR2
          BX2 = BX2*RBR2
          BY2 = BY2*RBR2
          BZ2 = BZ2*RBR2
          CX2 = CX2*RCR2
          CY2 = CY2*RCR2
          CZ2 = CZ2*RCR2
C
C calculate cos(phi) & sin(phi)
C
          CP = AX*BX+AY*BY+AZ*BZ
          SP = CX*BX+CY*BY+CZ*BZ
C
C calculate cos(psi) & sin(psi)
C
          CP2 = AX2*BX2+AY2*BY2+AZ2*BZ2
          SP2 = CX2*BX2+CY2*BY2+CZ2*BZ2
C
C calculate phi (make sure cos is within bounds and get sign from sin)
C and keep it it radians
C
          PHI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
          PSI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP2))),SP2)
C
C figure out which bin you should go into, given that the
C expectation values are in the range 1..PHISTEPS
C
C added the possibility of different phase corrections
C for different classes of PROCHECK data (eg., phi/psi
C vs. CHI1/CHI2 plots) 9/1/95
C
          if (phi.lt.CORRECT1) then
               CORRECTPHI = PHI + CORRECT2
          ELSE
               CORRECTPHI = PHI
          END IF
          RPHI = INT ((CORRECTPHI + CORRECT3) * PHICONV) + ONE
          IF (PSI.LT.CORRECT4) THEN
               CORRECTPSI = PSI + CORRECT5
          ELSE
               CORRECTPSI = PSI
          END IF
          RPSI = INT ((CORRECTPSI + CORRECT6) * PSICONV) + ONE
C
C enforce bounds
C
          IF (RPHI.GT.CURPHIS) THEN
               RPHI = RPHI - CURPHIS
          ELSE IF (RPHI.LT.1) THEN
               RPHI = RPHI + CURPHIS
          END IF
          IF (RPSI.GT.CURPSIS) THEN
               RPSI = RPSI - CURPSIS
          ELSE IF (RPSI.LT.1) THEN
               RPSI = RPSI + CURPSIS
          END IF
C
C look up expectation values
C
          OMEGAEQUIL = EXPECTEDANGLEDB(RPHI,RPSI)
C
C part two:  get the current torsion angle
C
C this is copied from ERAMA1D
C
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
          XRS = XR - XS
          XST = XS - XT
          XTU = XT - XU
C
          YRS = YR - YS
          YST = YS - YT
          YTU = YT - YU
C
          ZRS = ZR - ZS
          ZST = ZS - ZT
          ZTU = ZT - ZU
C
C now calculate A=RIJ*RJK, B = RJK*RKL, C = RJK*(RIJ*RJK)
C
          AX3 = YRS*ZST-ZRS*YST
          AY3 = ZRS*XST-XRS*ZST
          AZ3 = XRS*YST-YRS*XST
          BX3 = YST*ZTU-YTU*ZST
          BY3 = ZST*XTU-ZTU*XST
          BZ3 = XST*YTU-XTU*YST
          CX3 = YST*AZ3-ZST*AY3
          CY3 = ZST*AX3-XST*AZ3
          CZ3 = XST*AY3-YST*AX3
C
C calculate the norm of A, B, & C & set to MCONST if it's too small
C
          RAR3 = ONE/SQRT(MAX(MCONST, AX3**2+AY3**2+AZ3**2))
          RBR3 = ONE/SQRT(MAX(MCONST, BX3**2+BY3**2+BZ3**2))
          RCR3 = ONE/SQRT(MAX(MCONST, CX3**2+CY3**2+CZ3**2))
C
C normalize A, B, & C
C
          AX3 = AX3*RAR3
          AY3 = AY3*RAR3
          AZ3 = AZ3*RAR3
          BX3 = BX3*RBR3
          BY3 = BY3*RBR3
          BZ3 = BZ3*RBR3
          CX3 = CX3*RCR3
          CY3 = CY3*RCR3
          CZ3 = CZ3*RCR3
C
C calculate cos(phi) & sin(phi)
C
          CP3 = AX3*BX3+AY3*BY3+AZ3*BZ3
          SP3 = CX3*BX3+CY3*BY3+CZ3*BZ3
C
C calculate OMEGA (make sure cos is within
C bounds and get sign from sin) and keep it it radians
C
          OMEGA = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP3))),SP3)
C
C handle phase corrections as before
C
          if (OMEGA.LT.CORRECT7) then
               CORRECTOMEGA = OMEGA + CORRECT8
          ELSE
               CORRECTOMEGA = OMEGA
          END IF
C
C part three:  calculate energy and
C              simple parts of the force
C
C
C don't calculate energy & forces if you're at a spot without
C an expectation value
C
          IF (OMEGAEQUIL.GT.NOEXPTEST) THEN
               E = 0
               DEDOMEGA = 0
               DEDOMEGAEQUIL = 0
C
C Calculate energy & forces
C
          ELSE
               E = K1 * (CORRECTOMEGA - OMEGAEQUIL)**2
               EDB=EDB+E
               DEDOMEGA = 2 * K1 * (CORRECTOMEGA - OMEGAEQUIL)
               DEDOMEGAEQUIL = -2 * K1 * (CORRECTOMEGA - OMEGAEQUIL)
          END IF
C
C if we're in print mode, copy out the relevant data
C
          IF (WHICH.EQ.'ANALYZE') THEN
               CALCANGLEDB(1,COUNT) = OMEGAEQUIL
               CALCANGLEDB(2,COUNT) = OMEGA
               CALCANGLEDB(3,COUNT) = CORRECTOMEGA-OMEGAEQUIL
               CALCANGLEDB(4,COUNT) = PHI
               CALCANGLEDB(5,COUNT) = PSI
               CALCANGLEDB(6,COUNT) = E
          END IF
C
C accumulate energy
C
          ER = ER + E
C
C part four:  calculate d OMEGA / d xi
C             this is lifted from ERAMA1D
C
C
C compute heavyside function
C
          SWITCH3=-MIN(1,INT(ABS(SP3)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
          DF = DEDOMEGA
          RECSP3=DF*SWITCH3*SIGN(ONE,SP3)/MAX(ABS(SP3),MCONST)
          RECCP3=DF*(SWITCH3+1)*SIGN(ONE,CP3)/MAX(ABS(CP3),MCONST)
C
C compute:
C  d cos( PHI )  d cos ( PHI )
C  ------------, -------------
C     d A             d B
C
          DCPAX3=-RAR3*(BX3-CP3*AX3)
          DCPAY3=-RAR3*(BY3-CP3*AY3)
          DCPAZ3=-RAR3*(BZ3-CP3*AZ3)
          DCPBX3=-RBR3*(AX3-CP3*BX3)
          DCPBY3=-RBR3*(AY3-CP3*BY3)
          DCPBZ3=-RBR3*(AZ3-CP3*BZ3)
C
C compute:
C  d sin( PHI )  d sin ( PHI )
C  ------------, -------------
C     d C             d B
C
          DSPCX3=-RCR3*(BX3-SP3*CX3)
          DSPCY3=-RCR3*(BY3-SP3*CY3)
          DSPCZ3=-RCR3*(BZ3-SP3*CZ3)
          DSPBX3=-RBR3*(CX3-SP3*BX3)
          DSPBY3=-RBR3*(CY3-SP3*BY3)
          DSPBZ3=-RBR3*(CZ3-SP3*BZ3)
C
C compute:
C    d (phi)
C DF -------,  etc. (get rid of singularity by using two alternatives)
C    d rij
C
          DPRIJX3=
     &     RECSP3*(YST*DCPAZ3-DCPAY3*ZST)+
     &     RECCP3*((YST**2+ZST**2)*DSPCX3-
     &         XST*YST*DSPCY3-XST*ZST*DSPCZ3)
          DPRIJY3=
     &     RECSP3*(ZST*DCPAX3-DCPAZ3*XST)+
     &     RECCP3*((ZST**2+XST**2)*DSPCY3-
     &         YST*ZST*DSPCZ3-YST*XST*DSPCX3)
          DPRIJZ3=
     &     RECSP3*(XST*DCPAY3-DCPAX3*YST)+
     &     RECCP3*((XST**2+YST**2)*DSPCZ3-
     &         ZST*XST*DSPCX3-ZST*YST*DSPCY3)
C
          DPRKLX3=
     &     RECSP3*(DCPBY3*ZST-YST*DCPBZ3)+
     &     RECCP3*(DSPBY3*ZST-YST*DSPBZ3)
          DPRKLY3=
     &     RECSP3*(DCPBZ3*XST-ZST*DCPBX3)+
     &     RECCP3*(DSPBZ3*XST-ZST*DSPBX3)
          DPRKLZ3=
     &     RECSP3*(DCPBX3*YST-XST*DCPBY3)+
     &     RECCP3*(DSPBX3*YST-XST*DSPBY3)
C
          DPRJKX3=
     &     RECSP3*(DCPAY3*ZRS-DCPAZ3*YRS+DCPBZ3*YTU-DCPBY3*ZTU)+
     &     RECCP3*(-(YST*YRS+ZST*ZRS)*DSPCX3
     &        +(TWO*XST*YRS-XRS*YST)*DSPCY3
     &        +(TWO*XST*ZRS-XRS*ZST)*DSPCZ3
     &        +DSPBZ3*YTU-DSPBY3*ZTU)
          DPRJKY3=
     &     RECSP3*(DCPAZ3*XRS-DCPAX3*ZRS+DCPBX3*ZTU-DCPBZ3*XTU)+
     &     RECCP3*(-(ZST*ZRS+XST*XRS)*DSPCY3
     &        +(TWO*YST*ZRS-YRS*ZST)*DSPCZ3
     &        +(TWO*YST*XRS-YRS*XST)*DSPCX3
     &        +DSPBX3*ZTU-DSPBZ3*XTU)
          DPRJKZ3=
     &     RECSP3*(DCPAX3*YRS-DCPAY3*XRS+DCPBY3*XTU-DCPBX3*YTU)+
     &     RECCP3*(-(XST*XRS+YST*YRS)*DSPCZ3
     &        +(TWO*ZST*XRS-ZRS*XST)*DSPCX3
     &        +(TWO*ZST*YRS-ZRS*YST)*DSPCY3
     &        +DSPBY3*XTU-DSPBX3*YTU)
C
C Since overall force eqn is dE/dOMEGA dOMEGA/dxi + ...,
C I can update the forces on the three atoms involved
C in calculating OMEGA right now
C
          IF (WHICH.NE.'ANALYZE') THEN
               DX(ATOMR(COUNT))=DX(ATOMR(COUNT))+DPRIJX3
               DY(ATOMR(COUNT))=DY(ATOMR(COUNT))+DPRIJY3
               DZ(ATOMR(COUNT))=DZ(ATOMR(COUNT))+DPRIJZ3
               DX(ATOMS(COUNT))=DX(ATOMS(COUNT))+DPRJKX3-DPRIJX3
               DY(ATOMS(COUNT))=DY(ATOMS(COUNT))+DPRJKY3-DPRIJY3
               DZ(ATOMS(COUNT))=DZ(ATOMS(COUNT))+DPRJKZ3-DPRIJZ3
               DX(ATOMT(COUNT))=DX(ATOMT(COUNT))+DPRKLX3-DPRJKX3
               DY(ATOMT(COUNT))=DY(ATOMT(COUNT))+DPRKLY3-DPRJKY3
               DZ(ATOMT(COUNT))=DZ(ATOMT(COUNT))+DPRKLZ3-DPRJKZ3
               DX(ATOMU(COUNT))=DX(ATOMU(COUNT))        -DPRKLX3
               DY(ATOMU(COUNT))=DY(ATOMU(COUNT))        -DPRKLY3
               DZ(ATOMU(COUNT))=DZ(ATOMU(COUNT))        -DPRKLZ3
          END IF
C
C part 5: get d OMEGAequil / d phi and d OMEGAequil / d psi
C from slopes of expectation grid along phi and psi
C
C get slope along phi
C
C figure out what the neighboring squares are
C and get their energies from the grid
C
          LOWRPHI = RPHI - 1
          IF (LOWRPHI.LE.0) LOWRPHI = LOWRPHI + CURPHIS
          HIRPHI = RPHI + 1
          IF (HIRPHI.GT.CURPHIS) HIRPHI = HIRPHI-CURPHIS
          HIPHI = EXPECTEDANGLEDB(HIRPHI, RPSI)
          LOWPHI = EXPECTEDANGLEDB(LOWRPHI, RPSI)
C
C deal with the possibility of being on the edge of
C an unknown area
C
          IF ((HIPHI.GT.NOEXPTEST).AND.(LOWPHI.GT.NOEXPTEST)) THEN
               DOMEGAEQUILDPHI = ZERO
          ELSE IF (HIPHI.GT.NOEXPTEST) THEN
               DOMEGAEQUILDPHI = ((OMEGAEQUIL - LOWPHI) / ONE)
     &              * PHICONV
          ELSE IF (LOWPHI.GT.NOEXPTEST) THEN
               DOMEGAEQUILDPHI = ((HIPHI - OMEGAEQUIL) / ONE)
     &              * PHICONV
          ELSE
               DOMEGAEQUILDPHI = ((HIPHI - LOWPHI) / TWO)
     &              * PHICONV
          END IF
C
C do the same thing for the derivative along psi
C
          LOWRPSI = RPSI - 1
          IF (LOWRPSI.LE.0) LOWRPSI = LOWRPSI + CURPSIS
          HIRPSI = RPSI + 1
          IF (HIRPSI.GT.CURPSIS) HIRPSI = HIRPSI-CURPSIS
          HIPSI = EXPECTEDANGLEDB(RPHI, HIRPSI)
          LOWPSI = EXPECTEDANGLEDB(RPHI, LOWRPSI)
          IF (ANGLEDBPOTENTIAL.EQ.SQUARE) THEN
               HIPSI = MAX(HIPSI-ERR, 0.0D0)
               LOWPSI = MAX(LOWPSI-ERR,0.0D0)
          END IF
          IF ((HIPSI.GT.NOEXPTEST).AND.(LOWPSI.GT.NOEXPTEST)) THEN
               DOMEGAEQUILDPSI = 0
          ELSE IF (HIPSI.GT.NOEXPTEST) THEN
               DOMEGAEQUILDPSI = ((OMEGAEQUIL - LOWPSI) / ONE)
     &              * PSICONV
          ELSE IF (LOWPSI.GT.NOEXPTEST) THEN
               DOMEGAEQUILDPSI = ((HIPSI - OMEGAEQUIL) / ONE)
     &              * PSICONV
          ELSE
               DOMEGAEQUILDPSI = ((HIPSI - LOWPSI) / TWO)
     &              * PSICONV
          END IF
C
C part six: calculate d phi / dxi, d psi / dxi etc
C and in the process calcualte dE / dxi etc.
C
C compute heavyside function
C
          SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
          SWITCH2=-MIN(1,INT(ABS(SP2)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
C
C set DF to dE/dPhi and DF2 to dE/dPsi
C
          DF = DEDOMEGAEQUIL * DOMEGAEQUILDPHI
          DF2 = DEDOMEGAEQUIL * DOMEGAEQUILDPSI
          RECSP=DF*SWITCH*SIGN(ONE,SP)/MAX(ABS(SP),MCONST)
          RECCP=DF*(SWITCH+1)*SIGN(ONE,CP)/MAX(ABS(CP),MCONST)
          RECSP2=DF2*SWITCH2*SIGN(ONE,SP2)/MAX(ABS(SP2),MCONST)
          RECCP2=DF2*(SWITCH2+1)*SIGN(ONE,CP2)/MAX(ABS(CP2),MCONST)
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
          DCPAX2=-RAR2*(BX2-CP2*AX2)
          DCPAY2=-RAR2*(BY2-CP2*AY2)
          DCPAZ2=-RAR2*(BZ2-CP2*AZ2)
          DCPBX2=-RBR2*(AX2-CP2*BX2)
          DCPBY2=-RBR2*(AY2-CP2*BY2)
          DCPBZ2=-RBR2*(AZ2-CP2*BZ2)
C
C compute:
C  d sin( PHI )  d sin ( PHI )
C  ------------, -------------
C     d C             d B
C
          DSPCX=-RCR*(BX-SP*CX)
          DSPCY=-RCR*(BY-SP*CY)
          DSPCZ=-RCR*(BZ-SP*CZ)
          DSPBX=-RBR*(CX-SP*BX)
          DSPBY=-RBR*(CY-SP*BY)
          DSPBZ=-RBR*(CZ-SP*BZ)
C
          DSPCX2=-RCR2*(BX2-SP2*CX2)
          DSPCY2=-RCR2*(BY2-SP2*CY2)
          DSPCZ2=-RCR2*(BZ2-SP2*CZ2)
          DSPBX2=-RBR2*(CX2-SP2*BX2)
          DSPBY2=-RBR2*(CY2-SP2*BY2)
          DSPBZ2=-RBR2*(CZ2-SP2*BZ2)
C
C compute:
C    d (phi)
C DF -------,  etc. (get rid of singularity by using two alternatives)
C    d rij
C
          DPRIJX=
     &     RECSP*(YJK*DCPAZ-DCPAY*ZJK)+
     &     RECCP*((YJK**2+ZJK**2)*DSPCX-XJK*YJK*DSPCY-XJK*ZJK*DSPCZ)
          DPRIJY=
     &     RECSP*(ZJK*DCPAX-DCPAZ*XJK)+
     &     RECCP*((ZJK**2+XJK**2)*DSPCY-YJK*ZJK*DSPCZ-YJK*XJK*DSPCX)
          DPRIJZ=
     &     RECSP*(XJK*DCPAY-DCPAX*YJK)+
     &     RECCP*((XJK**2+YJK**2)*DSPCZ-ZJK*XJK*DSPCX-ZJK*YJK*DSPCY)
C
          DPRKLX=
     &     RECSP*(DCPBY*ZJK-YJK*DCPBZ)+
     &     RECCP*(DSPBY*ZJK-YJK*DSPBZ)
          DPRKLY=
     &     RECSP*(DCPBZ*XJK-ZJK*DCPBX)+
     &     RECCP*(DSPBZ*XJK-ZJK*DSPBX)
          DPRKLZ=
     &     RECSP*(DCPBX*YJK-XJK*DCPBY)+
     &     RECCP*(DSPBX*YJK-XJK*DSPBY)
C
          DPRJKX=
     &     RECSP*(DCPAY*ZIJ-DCPAZ*YIJ+DCPBZ*YKL-DCPBY*ZKL)+
     &     RECCP*(-(YJK*YIJ+ZJK*ZIJ)*DSPCX
     &        +(TWO*XJK*YIJ-XIJ*YJK)*DSPCY
     &        +(TWO*XJK*ZIJ-XIJ*ZJK)*DSPCZ
     &        +DSPBZ*YKL-DSPBY*ZKL)
          DPRJKY=
     &     RECSP*(DCPAZ*XIJ-DCPAX*ZIJ+DCPBX*ZKL-DCPBZ*XKL)+
     &     RECCP*(-(ZJK*ZIJ+XJK*XIJ)*DSPCY
     &        +(TWO*YJK*ZIJ-YIJ*ZJK)*DSPCZ
     &        +(TWO*YJK*XIJ-YIJ*XJK)*DSPCX
     &        +DSPBX*ZKL-DSPBZ*XKL)
          DPRJKZ=
     &     RECSP*(DCPAX*YIJ-DCPAY*XIJ+DCPBY*XKL-DCPBX*YKL)+
     &     RECCP*(-(XJK*XIJ+YJK*YIJ)*DSPCZ
     &        +(TWO*ZJK*XIJ-ZIJ*XJK)*DSPCX
     &        +(TWO*ZJK*YIJ-ZIJ*YJK)*DSPCY
     &        +DSPBY*XKL-DSPBX*YKL)
C
C
          DPRIJX2=
     &     RECSP2*(YKL*DCPAZ2-DCPAY2*ZKL)+
     &     RECCP2*((YKL**2+ZKL**2)*DSPCX2-XKL*YKL*DSPCY2-
     &     XKL*ZKL*DSPCZ2)
          DPRIJY2=
     &     RECSP2*(ZKL*DCPAX2-DCPAZ2*XKL)+
     &     RECCP2*((ZKL**2+XKL**2)*DSPCY2-YKL*ZKL*DSPCZ2-
     &     YKL*XKL*DSPCX2)
          DPRIJZ2=
     &     RECSP2*(XKL*DCPAY2-DCPAX2*YKL)+
     &     RECCP2*((XKL**2+YKL**2)*DSPCZ2-ZKL*XKL*DSPCX2-
     &     ZKL*YKL*DSPCY2)
C
          DPRKLX2=
     &     RECSP2*(DCPBY2*ZKL-YKL*DCPBZ2)+
     &     RECCP2*(DSPBY2*ZKL-YKL*DSPBZ2)
          DPRKLY2=
     &     RECSP2*(DCPBZ2*XKL-ZKL*DCPBX2)+
     &     RECCP2*(DSPBZ2*XKL-ZKL*DSPBX2)
          DPRKLZ2=
     &     RECSP2*(DCPBX2*YKL-XKL*DCPBY2)+
     &     RECCP2*(DSPBX2*YKL-XKL*DSPBY2)
C
          DPRJKX2=
     &     RECSP2*(DCPAY2*ZJK-DCPAZ2*YJK+DCPBZ2*ypq-DCPBY2*zpq)+
     &     RECCP2*(-(YKL*YJK+ZKL*ZJK)*DSPCX2
     &        +(TWO*XKL*YJK-XJK*YKL)*DSPCY2
     &        +(TWO*XKL*ZJK-XJK*ZKL)*DSPCZ2
     &        +DSPBZ2*ypq-DSPBY2*zpq)
          DPRJKY2=
     &     RECSP2*(DCPAZ2*XJK-DCPAX2*ZJK+DCPBX2*zpq-DCPBZ2*xpq)+
     &     RECCP2*(-(ZKL*ZJK+XKL*XJK)*DSPCY2
     &        +(TWO*YKL*ZJK-YJK*ZKL)*DSPCZ2
     &        +(TWO*YKL*XJK-YJK*XKL)*DSPCX2
     &        +DSPBX2*zpq-DSPBZ2*xpq)
          DPRJKZ2=
     &     RECSP2*(DCPAX2*YJK-DCPAY2*XJK+DCPBY2*xpq-DCPBX2*ypq)+
     &     RECCP2*(-(XKL*XJK+YKL*YJK)*DSPCZ2
     &        +(TWO*ZKL*XJK-ZJK*XKL)*DSPCX2
     &        +(TWO*ZKL*YJK-ZJK*YKL)*DSPCY2
     &        +DSPBY2*xpq-DSPBX2*ypq)
C
C now update forces if in energy & force mode
C
          IF ((WHICH.NE.'ANALYZE').AND.(USEDERIV)) THEN
               DX(ATOMI(COUNT))=DX(ATOMI(COUNT))+DPRIJX
               DY(ATOMI(COUNT))=DY(ATOMI(COUNT))+DPRIJY
               DZ(ATOMI(COUNT))=DZ(ATOMI(COUNT))+DPRIJZ
               DX(ATOMJ(COUNT))=DX(ATOMJ(COUNT))+DPRJKX-DPRIJX
               DY(ATOMJ(COUNT))=DY(ATOMJ(COUNT))+DPRJKY-DPRIJY
               DZ(ATOMJ(COUNT))=DZ(ATOMJ(COUNT))+DPRJKZ-DPRIJZ
               DX(ATOMK(COUNT))=DX(ATOMK(COUNT))+DPRKLX-DPRJKX
               DY(ATOMK(COUNT))=DY(ATOMK(COUNT))+DPRKLY-DPRJKY
               DZ(ATOMK(COUNT))=DZ(ATOMK(COUNT))+DPRKLZ-DPRJKZ
               DX(ATOML(COUNT))=DX(ATOML(COUNT))       -DPRKLX
               DY(ATOML(COUNT))=DY(ATOML(COUNT))       -DPRKLY
               DZ(ATOML(COUNT))=DZ(ATOML(COUNT))       -DPRKLZ
C
               DX(ATOMM(COUNT))=DX(ATOMM(COUNT))+DPRIJX2
               DY(ATOMM(COUNT))=DY(ATOMM(COUNT))+DPRIJY2
               DZ(ATOMM(COUNT))=DZ(ATOMM(COUNT))+DPRIJZ2
               DX(ATOMN(COUNT))=DX(ATOMN(COUNT))+DPRJKX2-DPRIJX2
               DY(ATOMN(COUNT))=DY(ATOMN(COUNT))+DPRJKY2-DPRIJY2
               DZ(ATOMN(COUNT))=DZ(ATOMN(COUNT))+DPRJKZ2-DPRIJZ2
               DX(ATOMP(COUNT))=DX(ATOMP(COUNT))+DPRKLX2-DPRJKX2
               DY(ATOMP(COUNT))=DY(ATOMP(COUNT))+DPRKLY2-DPRJKY2
               DZ(ATOMP(COUNT))=DZ(ATOMP(COUNT))+DPRKLZ2-DPRJKZ2
               DX(ATOMQ(COUNT))=DX(ATOMQ(COUNT))        -DPRKLX2
               DY(ATOMQ(COUNT))=DY(ATOMQ(COUNT))        -DPRKLY2
               DZ(ATOMQ(COUNT))=DZ(ATOMQ(COUNT))        -DPRKLZ2
          END IF
      END DO
      RETURN
      END
C================
      SUBROUTINE READANGLEDB
C
C reads in bond angle database information
C
C by John Kuszewski June 1996
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXANGLEDBS
      INTEGER THETYPE, CURPSIS, CURPHIS
      INTEGER CLASINDEX
      DOUBLE PRECISION K1, CUTOFF
      CHARACTER*4 THENAME
      CHARACTER*20 CLASNAME
C begin
C
      IF (ANGLEDBFLAG) THEN
      ANGLEDBFLAG=.FALSE.
      CALL ANGLEDBDEFAULTS
      CALL ALLOCANGLEDBS(0, MAXANGLEDBS)
      END IF
C
C this is used by READANGLEDB2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
C now read input
C
      CALL PUSEND('ANGLE DATABASE>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('ANGLE DATABASE>')
           CALL MISCOM('ANGLE DATABASE>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-angledatabase')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURANGLEDBCLASS
                CALL NEXTWD('class name =')
                CLASNAME = WD(1:20)
                ANGLEDBMODE = NEW
                DO COUNT = 1, NANGLEDBCLASSES
                     IF (ANGDBCLASSNAMES(COUNT).EQ.CLASNAME) THEN
                          ANGLEDBMODE = UPDATE
                          CURANGLEDBCLASS = COUNT
                     END IF
                END DO
                IF (ANGLEDBMODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXANGLEDBCLASSES) THEN
                         CALL DSPERR('ANGDB','Too many classes.')
                         CALL DSPERR('ANGDB',
     &                     'Increase MAXANGLEDBCLASSES and recompile.')
                         CALL WRNDIE(-5, 'READANGDB',
     &                               'Too many ANGDB classes.')
                     END IF
                     NANGLEDBCLASSES = NANGLEDBCLASSES + 1
                     CURANGLEDBCLASS = NANGLEDBCLASSES
                     ANGDBCLASSNAMES(CURANGLEDBCLASS) = CLASNAME
                     ANGLEDBASSNDX(CURANGLEDBCLASS) = NANGLEDBS
                END IF
C
C If this isn't the first class, close off the old class
C
C                IF (NANGLEDBCLASSES.GT.1) THEN
C                     ANGLEDBASSNDX(OLDCLASS) = NANGLEDBS
C                END IF
C
C set force constant for current class,
C starting a default class if there isnt one defined
C
           ELSE IF (WD(1:4).EQ.'FORC') THEN
                CALL NEXTF('force constant =', K1)
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                WRITE(PUNIT, '(3(A), F8.3)')
     &            'Setting frce const for class ',
     &             ANGDBCLASSNAMES(CURANGLEDBCLASS), ' to ', K1
                ANGLEDBFORCES(CURANGLEDBCLASS) = K1
C
C deriv flags
C starting a default class if there isnt one defined
C
           ELSE IF (WD(1:4).EQ.'DERI') THEN
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                CALL NEXTWD('on or off')
                IF (WD(1:3).EQ.'OFF') THEN
                     ANGLEDBDERIV(CURANGLEDBCLASS) = .FALSE.
                ELSE IF (WD(1:2).EQ.'ON') THEN
                     ANGLEDBDERIV(CURANGLEDBCLASS) = .TRUE.
                ELSE
                     CALL DSPERR('ANGDB', 'Must be on or off.')
                     ANGLEDBDERIV(CURANGLEDBCLASS) = .FALSE.
                END IF
C
C set phase corrections for current class,
C starting a default class if there isnt one defined
C
           ELSE IF (WD(1:4).EQ.'PHAS') THEN
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                CALL NEXTF('phase correction 1 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(1, CURANGLEDBCLASS) = K1
                CALL NEXTF('phase correction 2 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(2, CURANGLEDBCLASS) = K1
                CALL NEXTF('phase correction 3 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(3, CURANGLEDBCLASS) = K1
                CALL NEXTF('phase correction 4 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(4, CURANGLEDBCLASS) = K1
                CALL NEXTF('phase correction 5 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(5, CURANGLEDBCLASS) = K1
                CALL NEXTF('phase correction 6 (deg) =', K1)
                CALL DG2RAD(K1, K1)
                ANGLEDBPHASEC(6, CURANGLEDBCLASS) = K1
                IF (ANGLEDBCLASSTYPE(CURANGLEDBCLASS).EQ.TORSION) THEN
                     CALL NEXTF('phase correction 7 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     ANGLEDBPHASEC(7, CURANGLEDBCLASS) = K1
                     CALL NEXTF('phase correction 8 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     ANGLEDBPHASEC(8, CURANGLEDBCLASS) = K1
                     CALL NEXTF('phase correction 9 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     ANGLEDBPHASEC(9, CURANGLEDBCLASS) = K1
                ELSE
                     ANGLEDBPHASEC(7, CURANGLEDBCLASS) = ZERO
                     ANGLEDBPHASEC(8, CURANGLEDBCLASS) = ZERO
                     ANGLEDBPHASEC(9, CURANGLEDBCLASS) = ZERO
                END IF
                WRITE(PUNIT, '(3(A))')
     &               'Setting phase corrs for class ',
     &               ANGDBCLASSNAMES(CURANGLEDBCLASS),
     &               ' to (in radians)'
                WRITE(PUNIT, '(3(F8.3))')
     &               ANGLEDBPHASEC(1, CURANGLEDBCLASS),
     &               ANGLEDBPHASEC(2, CURANGLEDBCLASS),
     &               ANGLEDBPHASEC(3, CURANGLEDBCLASS)
                WRITE(PUNIT, '(3(F8.3))')
     &               ANGLEDBPHASEC(4, CURANGLEDBCLASS),
     &               ANGLEDBPHASEC(5, CURANGLEDBCLASS),
     &               ANGLEDBPHASEC(6, CURANGLEDBCLASS)
                IF (ANGLEDBCLASSTYPE(CURANGLEDBCLASS).EQ.TORSION) THEN
                     WRITE(PUNIT, '(3(F8.3))')
     &                    ANGLEDBPHASEC(7, CURANGLEDBCLASS),
     &                    ANGLEDBPHASEC(8, CURANGLEDBCLASS),
     &                    ANGLEDBPHASEC(9, CURANGLEDBCLASS)
                END IF
C
C set size of current class's expectation grid,
C starting a default class if there isn't one defined
C
           ELSE IF (WD(1:4).EQ.'SIZE') THEN
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                CALL NEXTA4('type of angle =', THENAME)
                IF (THENAME.EQ.'ANGL') THEN
                     THETYPE = ANGLE
                ELSE IF (THENAME.EQ.'DIHE') THEN
                     THETYPE = TORSION
                ELSE
                     CALL DSPERR('ANGLE DATABASE',
     &                    'Not understood.  Assuming dihedral.')
                     THETYPE = ANGLE
                END IF
                CALL NEXTI('phi step =', CURPHIS)
                CALL NEXTI('psi step =', CURPSIS)
                WRITE (PUNIT, '(3(A), I6, A, I6)')
     &                    'Setting size of two-D class ',
     &                    ANGDBCLASSNAMES(CURANGLEDBCLASS), ' to ',
     &                    CURPHIS, ' by ', CURPSIS
                CALL ALLOCEXPECTEDANGLEDB(CURANGLEDBCLASS,
     &               PHISTEPS(CURANGLEDBCLASS)*
     &               PSISTEPS(CURANGLEDBCLASS),
     &               CURPHIS*CURPSIS)
                ANGLEDBCLASSTYPE(CURANGLEDBCLASS) = THETYPE
                PHISTEPS(CURANGLEDBCLASS) = CURPHIS
                PSISTEPS(CURANGLEDBCLASS) = CURPSIS
C
C set error for current class,
C starting a default class if there isn't one defined
C
           ELSE IF (WD(1:4).EQ.'ERRO') THEN
                CALL NEXTF('error =', K1)
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                WRITE(PUNIT, '(3(A), F8.3)')
     &            'Setting error for class ',
     &             ANGDBCLASSNAMES(CURANGLEDBCLASS), ' to ', K1
                ANGLEDBERRORS(CURANGLEDBCLASS) = K1
C
C reset angle database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
                CALL ANGLEDBDEFAULTS
                CALL ALLOCANGLEDBS(MAXANGLEDBS, MAXANGLEDBS)
C
C set potential type
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                    WRITE(PUNIT, '(A)') 'using square well potential.'
                    ANGLEDBPOTENTIAL = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                    WRITE(PUNIT, '(A)') 'using harmonic potential.'
                    ANGLEDBPOTENTIAL = HARMONIC
                ELSE
                    CALL DSPERR('ANGLEDB',
     &                        'unknown potential. Using square well.')
                    ANGLEDBPOTENTIAL = SQUARE
                END IF
C
C change number of assignment slot
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXANGLEDBS = MAXANGLEDBS
                CALL NEXTI('number of slots =', MAXANGLEDBS)
                CALL ALLOCANGLEDBS(OLDMAXANGLEDBS, MAXANGLEDBS)
C
C debugging
C
           ELSE IF (WD(1:4).EQ.'DEBU') THEN
                WRITE(6,'(A,I5)')'NANGLEDBS = ',NANGLEDBS
                WRITE(6,'(A,I5)')'NANGLEDBCLASSES = ',NANGLEDBCLASSES
                DO COUNT = 1, NANGLEDBCLASSES
                   WRITE(6,'(A,I5,I5)')'CLASS, NDX =', COUNT,
     &                                  ANGLEDBASSNDX(COUNT)
                END DO
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more constraints
C than you have slots for
C
                IF (NANGLEDBS.EQ.MAXANGLEDBS) THEN
                     CALL DSPERR('ANGLEDB','Too many assignments.')
                     CALL DSPERR('ANGLEDB',
     &                    'Re-run with larger NREStraints.')
                     CALL WRNDIE(-5, 'READANGDB', 'Out of space')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURANGLEDBCLASS.EQ.0) THEN
                     NANGLEDBCLASSES = 1
                     CURANGLEDBCLASS = 1
                END IF
                CALL READANGLEDBASSIGN( HEAP(ANGLEDBIPTR),
     &               HEAP(ANGLEDBJPTR), HEAP(ANGLEDBKPTR),
     &               HEAP(ANGLEDBLPTR), HEAP(ANGLEDBMPTR),
     &               HEAP(ANGLEDBNPTR), HEAP(ANGLEDBPPTR),
     &               HEAP(ANGLEDBQPTR), HEAP(ANGLEDBRPTR),
     &               HEAP(ANGLEDBSPTR), HEAP(ANGLEDBTPTR),
     &               heap(ANGLEDBUPTR),
     &               ANGLEDBCLASSTYPE(CURANGLEDBCLASS), HEAP(SPTR))
C
C read in expectation values
C
           ELSE IF (WD(1:4).EQ.'EXPE') THEN
                CALL READEXPECTEDANGLEDB(CURANGLEDBCLASS,
     &               HEAP(EXPANGLEDBPTRS(CURANGLEDBCLASS)),
     &               PHISTEPS(CURANGLEDBCLASS),
     &               PSISTEPS(CURANGLEDBCLASS))
C
C zero out the expectation value arrays
C
           ELSE IF (WD(1:4).EQ.'ZERO') THEN
                CALL DEFEXPECTEDANGLEDB(
     &               HEAP(EXPANGLEDBPTRS(CURANGLEDBCLASS)),
     &               PHISTEPS(CURANGLEDBCLASS)*
     &               PSISTEPS(CURANGLEDBCLASS))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                DO COUNT = 1, MAXANGLEDBCLASSES
                     PRINTTHISCLASS(COUNT) = .FALSE.
                END DO
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                     CALL DSPERR('ANGLEDB',
     &                           'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                          CALL DSPERR('ANGLEDB',
     &                         'cutoff must be positive.')
                          CUTOFF = ABS(CUTOFF)
                     END IF
                     CALL NEXTWD('ALL, TYPE, or CLASs>')
                     IF (WD(1:3).EQ.'ALL') THEN
                          DO COUNT = 1, NANGLEDBCLASSES
                               PRINTTHISCLASS(COUNT) = .TRUE.
                          END DO
                          CALL PRINTANGLEDBS(CUTOFF,
     &                         HEAP(CALCANGLEDBPTR),
     &                         HEAP(ANGLEDBKPTR), HEAP(ANGLEDBPPTR),
     &                         HEAP(ANGLEDBTPTR))
                     ELSE IF (WD(1:4).EQ.'CLAS') THEN
                          CALL NEXTWD('Class name>')
                          CLASNAME = WD(1:20)
                          CLASINDEX = 0
                          DO COUNT = 1, NANGLEDBCLASSES
                               IF (ANGDBCLASSNAMES(COUNT).EQ.
     &                              CLASNAME) THEN
                                    PRINTTHISCLASS(COUNT) = .TRUE.
                                    CLASINDEX = COUNT
                               END IF
                          END DO
                          IF (CLASINDEX.EQ.0) THEN
                               CALL DSPERR('ANGLEDB',
     &                              'unknown class. Using first.')
                               PRINTTHISCLASS(1) = .TRUE.
                          END IF
                          CALL PRINTANGLEDBS(CUTOFF,
     &                         HEAP(CALCANGLEDBPTR),
     &                         HEAP(ANGLEDBKPTR), HEAP(ANGLEDBPPTR),
     &                         HEAP(ANGLEDBTPTR))
                     ELSE IF (WD(1:4).EQ.'TYPE') THEN
                          CALL NEXTWD('ANGLe or TORSion>')
                          IF (WD(1:4).EQ.'ANGL') THEN
                               DO COUNT = 1, NANGLEDBCLASSES
                                    IF (ANGLEDBCLASSTYPE(COUNT).EQ.
     &                                   ANGLE) THEN
                                         PRINTTHISCLASS(COUNT) = .TRUE.
                                    END IF
                               END DO
                          ELSE IF (WD(1:4).EQ.'TORS') THEN
                               DO COUNT = 1, NANGLEDBCLASSES
                                    IF (ANGLEDBCLASSTYPE(COUNT).EQ.
     &                                   TORSION) THEN
                                         PRINTTHISCLASS(COUNT) = .TRUE.
                                    END IF
                               END DO
                          ELSE
                               CALL DSPERR('ANGLEDB',
     &                              'unknown type.  Assuming torsion.')
                               DO COUNT = 1, NANGLEDBCLASSES
                                    IF (ANGLEDBCLASSTYPE(COUNT).EQ.
     &                                   TORSION) THEN
                                         PRINTTHISCLASS(COUNT) = .TRUE.
                                    END IF
                               END DO
                          END IF
                          CALL PRINTANGLEDBS(CUTOFF,
     &                         HEAP(CALCANGLEDBPTR),
     &                         HEAP(ANGLEDBKPTR),
     &                         HEAP(ANGLEDBPPTR),
     &                         HEAP(ANGLEDBTPTR))
                     ELSE
                          CALL DSPERR('ANGLEDB',
     &                         'Expected ALL or CLASs.')
                     END IF
                END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('ANGLE DATABASE>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCEXPECTEDANGLEDB(CLASS, OLDSIZE, NEWSIZE)
C
C Allocates space for angle database expectation arrays
C based on values of PHISTEP and PSISTEP.
C Fills each array with zeros.
C
C by John Kuszewski May 1994
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE, CLASS
C local variables
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(EXPANGLEDBPTRS(CLASS), IREAL8(OLDSIZE))
      END IF
      EXPANGLEDBPTRS(CLASS) = 0
      IF (NEWSIZE.NE.0) THEN
        EXPANGLEDBPTRS(CLASS) = ALLHP(IREAL8(NEWSIZE))
C
C now zero them out
C
        CALL DEFEXPECTEDANGLEDB(HEAP(EXPANGLEDBPTRS(CLASS)), NEWSIZE)
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE ALLOCANGLEDBS (OLDSIZE, NEWSIZE)
C
C resets angle database to hold SIZE constraint entries
C
C by John Kuszewski Nov 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(ANGLEDBIPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBJPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBKPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBLPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBMPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBNPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBPPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBQPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBRPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBSPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBTPTR, INTEG4(OLDSIZE))
           CALL FREHP(ANGLEDBUPTR, INTEG4(OLDSIZE))
           CALL FREHP(CALCANGLEDBPTR, IREAL8(6*OLDSIZE))
      END IF
C
      ANGLEDBIPTR = 0
      ANGLEDBJPTR = 0
      ANGLEDBKPTR = 0
      ANGLEDBLPTR = 0
      ANGLEDBMPTR = 0
      ANGLEDBNPTR = 0
      ANGLEDBPPTR = 0
      ANGLEDBQPTR = 0
      ANGLEDBRPTR = 0
      ANGLEDBSPTR = 0
      ANGLEDBTPTR = 0
      ANGLEDBUPTR = 0
      CALCANGLEDBPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        ANGLEDBIPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBJPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBKPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBLPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBMPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBNPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBPPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBQPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBRPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBSPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBTPTR = ALLHP(INTEG4(NEWSIZE))
        ANGLEDBUPTR = ALLHP(INTEG4(NEWSIZE))
        CALCANGLEDBPTR = ALLHP(IREAL8(6*NEWSIZE))
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE ANGLEDBINIT
C
C initializes angle stuff
C
C by John Kuszewski Nov 1993
C================
      IMPLICIT NONE
C include files
      INCLUDE 'angledb.inc'
C local vbls
C begin
      ANGLEDBFLAG=.TRUE.
      RETURN
      END
C===============
      SUBROUTINE ANGLEDBFREE
C
C frees up angle stuff
C
C by Axel T. Brunger
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C local
      INTEGER COUNT
C begin
      IF (.NOT.ANGLEDBFLAG) THEN
      IF (MAXANGLEDBS.NE.0) THEN
           CALL FREHP(ANGLEDBIPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBJPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBKPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBLPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBMPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBNPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBPPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBQPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBRPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBSPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBTPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(ANGLEDBUPTR, INTEG4(MAXANGLEDBS))
           CALL FREHP(CALCANGLEDBPTR, IREAL8(6*MAXANGLEDBS))
      END IF
      DO COUNT = 1, MAXANGLEDBCLASSES
           CALL FREHP(EXPANGLEDBPTRS(COUNT),
     &     IREAL8(PHISTEPS(COUNT)*PSISTEPS(COUNT)))
      END DO
      END IF
      RETURN
      END
C==============
      SUBROUTINE ANGLEDBDEFAULTS
C
C sets up defaults
C
C by John Kuszewski Nov 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'angledb.inc'
C local variables
      INTEGER COUNT
C begin
      ANGLEDBMODE = NEW
      MAXANGLEDBS = 5000
      NANGLEDBS = 0
      NANGLEDBCLASSES = 0
      CURANGLEDBCLASS = 0
      ANGLEDBPOTENTIAL = HARMONIC
      DO COUNT = 1, MAXANGLEDBCLASSES
           ANGLEDBDERIV(COUNT) = .FALSE.
           PHISTEPS(COUNT) = 36
           PSISTEPS(COUNT) = 36
           ANGLEDBCLASSTYPE(COUNT) = 1
           ANGDBCLASSNAMES(COUNT) = 'DEFAULT'
           ANGLEDBASSNDX(COUNT) = 0
           ANGLEDBFORCES(COUNT) = 1.0
           ANGLEDBPHASEC(1, COUNT) = 0.0
           ANGLEDBPHASEC(2, COUNT) = 0.0
           ANGLEDBPHASEC(3, COUNT) = 0.0
           ANGLEDBPHASEC(4, COUNT) = 0.0
           ANGLEDBPHASEC(5, COUNT) = 0.0
           ANGLEDBPHASEC(6, COUNT) = 0.0
           ANGLEDBPHASEC(7, COUNT) = 0.0
           ANGLEDBPHASEC(8, COUNT) = 0.0
           ANGLEDBPHASEC(9, COUNT) = 0.0
           CALL ALLOCEXPECTEDANGLEDB(COUNT,0,
     &          PHISTEPS(COUNT)*PSISTEPS(COUNT))
      END DO
      RETURN
      END
C==============
      SUBROUTINE DEFEXPECTEDANGLEDB (EXPECTED, SIZE)
C
C fills arrays that hold expectation values with zero
C
C by John Kuszewski June 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'angledb.inc'
C i/o
      DOUBLE PRECISION EXPECTED(*)
      INTEGER SIZE
C local vbls
      INTEGER COUNT
C begin
      DO COUNT = 1, SIZE
           EXPECTED(COUNT) = NOEXPECTATION
      END DO
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDANGLEDB (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS)
C
C actually reads in the expected values for the bond angle
C or dihedral angle
C
C by John Kuszewski May 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'angledb.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS)
C local vbls
      INTEGER PHIPOS, PSIPOS
      DOUBLE PRECISION VAL
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS(CLASS)).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('ANGLEDB',
     &               'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTI('PSIposition =', PSIPOS)
      IF ((PSIPOS.GT.PSISTEPS(CLASS)).OR.(PSIPOS.LT.1)) THEN
           CALL DSPERR('ANGLEDB',
     &               'psi not in range 1..PSISTEPS. Assuming posn 1.')
           PSIPOS = 1
      END IF
      CALL NEXTF('expected angle (degrees) =', VAL)
      CALL DG2RAD(VAL, VAL)
      EXPECTED(PHIPOS,PSIPOS) = VAL
      RETURN
      END
C==============
      SUBROUTINE READANGLEDBASSIGN (ATOMI, ATOMJ, ATOMK,
     &     ATOML, ATOMM, ATOMN, ATOMP, ATOMQ, ATOMR,
     &     ATOMS, ATOMT, ATOMU, THETYPE, SEL)
C
C reads actual angle database assignments into arrays
C
C by John Kuszewski Nov 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*)
      INTEGER ATOMN(*), ATOMP(*), ATOMQ(*), ATOMR(*), ATOMS(*)
      INTEGER ATOMT(*), ATOMU(*), SEL(*), THETYPE
C local variables
      INTEGER NSEL, INSERTPOS, COUNT
      INTEGER I, J, K, L, M, N, P, Q, r, s, t, u
      LOGICAL USEME
C begin
C
      USEME = .TRUE.
C
C get the atoms and observed values
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom H. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      I = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom J. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      J = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom K. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      K = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom L. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      L = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom M. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      M = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom N. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      N = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom P. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      P = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom Q. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      Q = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom r. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      R = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom s. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      S = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('ANGLEDB',
     &     'more than 1 atom in selection for atom t. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      T = SEL(1)
C
      IF (THETYPE.EQ.TORSION) THEN
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('ANGLEDB',
     &         'more than 1 atom in selection for atom u. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           IF (NSEL.LT.1) USEME = .FALSE.
           U = SEL(1)
      ELSE
           U = 0
      END IF
C
C see if we should ignore this one or not
C
      IF (USEME) THEN
C
C if we're in update mode,
C
           IF (ANGLEDBMODE.EQ.UPDATE) THEN
C
C make a space for the new line in the atom* arrays
C
CCCC modification ATB 4/27/08
CCC                DO COUNT = NANGLEDBS+1,
CCC     &               ANGLEDBASSNDX(CURANGLEDBCLASS)+1, -1
                DO COUNT = NANGLEDBS+1,
     &               MAX(2,ANGLEDBASSNDX(CURANGLEDBCLASS)+1), -1
                     ATOMI(COUNT) = ATOMI(COUNT-1)
                     ATOMJ(COUNT) = ATOMJ(COUNT-1)
                     ATOMK(COUNT) = ATOMK(COUNT-1)
                     ATOML(COUNT) = ATOML(COUNT-1)
                     ATOMM(COUNT) = ATOMM(COUNT-1)
                     ATOMN(COUNT) = ATOMN(COUNT-1)
                     ATOMP(COUNT) = ATOMP(COUNT-1)
                     ATOMQ(COUNT) = ATOMQ(COUNT-1)
                     ATOMR(COUNT) = ATOMR(COUNT-1)
                     ATOMS(COUNT) = ATOMS(COUNT-1)
                     ATOMT(COUNT) = ATOMT(COUNT-1)
                     ATOMU(COUNT) = ATOMU(COUNT-1)
                END DO
C
C now update the index for each class that comes after
C the current class, in which the insertion will take place
C
                DO COUNT = CURANGLEDBCLASS, NANGLEDBCLASSES
                     ANGLEDBASSNDX(COUNT) = ANGLEDBASSNDX(COUNT) + 1
                END DO
                INSERTPOS = ANGLEDBASSNDX(CURANGLEDBCLASS)
                NANGLEDBS = NANGLEDBS + 1
           ELSE
                NANGLEDBS = NANGLEDBS + 1
                INSERTPOS = NANGLEDBS
                ANGLEDBASSNDX(CURANGLEDBCLASS) = INSERTPOS
           END IF
C
C now write the data to the arrays
C
           ATOMI(INSERTPOS) = I
           ATOMJ(INSERTPOS) = J
           ATOMK(INSERTPOS) = K
           ATOML(INSERTPOS) = L
           ATOMM(INSERTPOS) = M
           ATOMN(INSERTPOS) = N
           ATOMP(INSERTPOS) = P
           ATOMQ(INSERTPOS) = Q
           ATOMr(INSERTPOS) = R
           ATOMS(INSERTPOS) = S
           ATOMT(INSERTPOS) = T
           ATOMU(INSERTPOS) = U
      ELSE
           CALL DSPERR('ANGLEDB',
     &     'no atoms in one or more selections. Ignoring this entry.')
      END IF
      RETURN
      END
C=================
      SUBROUTINE PRINTANGLEDBS (CUTOFF, CALC, ATOMK, ATOMP,
     &     ATOMT)
C
C prints angle database constraints with energies greater than CUTOFF.
C calculates RMS deviation and puts it into $RMS
C
C by John Kuszewski July 1995
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'angledb.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, CALC(6,*)
      INTEGER ATOMK(*), ATOMP(*), ATOMT(*)
C local variables
      DOUBLE PRECISION CURENER, CURDELTA
      DOUBLE PRECISION RMS, VIOLS, NUMINRMS, RENERGY
      DOUBLE PRECISION CUREQUIL, CURANGLE, CURPHI, CURPSI
      DOUBLE PRECISION DBPREC
      INTEGER COUNT, CLASS, K, P, T, CURTYPE
      DOUBLE COMPLEX DUMMY2
C begin
      RMS = ZERO
      VIOLS = ZERO
      NUMINRMS = ZERO
C
C make sure that the calc-shift array is up to date
C
      CALL EANGLEDB(RENERGY, 'ANALYZE')
      WRITE (PUNIT, '(A)') 'The following angle database entries'
      WRITE (PUNIT, '(A)') 'have energies greater than the cutoff:'
C
C write out first class heading
C
      CLASS = 1
      CURTYPE = ANGLEDBCLASSTYPE(CLASS)
      IF (PRINTTHISCLASS(CLASS)) THEN
           IF (CURTYPE.EQ.TORSION) THEN
                WRITE (PUNIT, '(A, A)')
     &               'torsion class ', ANGDBCLASSNAMES(CLASS)
           ELSE
                WRITE (PUNIT, '(A, A)')
     &               'angle class ', ANGDBCLASSNAMES(CLASS)
           END IF
           WRITE (PUNIT, '(A, A)') '(atom K)(atom P)(atom T)',
     &          '(equil)(cur)(delta)(phi)(psi)(energy)'
      END IF
C
C loop through energy values
C
      DO COUNT = 1,NANGLEDBS
           CUREQUIL = CALC(1,COUNT)
           CURANGLE = CALC(2,COUNT)
           CURDELTA = CALC(3,COUNT)
           CURPHI = CALC(4,COUNT)
           CURPSI = CALC(5,COUNT)
           CURENER = CALC(6,COUNT)
C
C change from rad to deg
C
           CALL RD2DEG(CUREQUIL, CUREQUIL)
           CALL RD2DEG(CURANGLE, CURANGLE)
           CALL RD2DEG(CURPHI, CURPHI)
           CALL RD2DEG(CURPSI, CURPSI)
           CALL RD2DEG(CURDELTA, CURDELTA)
C
C make sure the delta is positive
C
           CURDELTA = ABS(CURDELTA)
C
C is this the start of a new class?
C
          IF (ANGLEDBASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               CURTYPE = ANGLEDBCLASSTYPE(CLASS)
               IF (PRINTTHISCLASS(CLASS)) THEN
                    IF (CURTYPE.EQ.TORSION) THEN
                         WRITE (PUNIT, '(A, A)')
     &                        'torsion class ', ANGDBCLASSNAMES(CLASS)
                    ELSE
                         WRITE (PUNIT, '(A, A)')
     &                        'angle class ', ANGDBCLASSNAMES(CLASS)
                    END IF
                    WRITE (PUNIT, '(A, A)') '(atom K)(atom P)(atom T)',
     &                   '(equil)(cur)(delta)(phi)(psi)(energy)'
               END IF
          END IF
C
C if we should print this class,
C
          IF (PRINTTHISCLASS(CLASS)) THEN
C
C update RMS values if we can make a guess about this entry
C
               RMS = RMS + CURDELTA**2
               NUMINRMS = NUMINRMS + ONE
C
C if the current value is greater than the CUTOFF value
C and is less than the no-expectation-test value,
C
C print it
C
               IF (CURDELTA.GT.CUTOFF) THEN
                    VIOLS = VIOLS + ONE
                    K = ATOMK(COUNT)
                    P = ATOMP(COUNT)
                    T = ATOMT(COUNT)
                    WRITE(PUNIT,'(12A,6(F8.2 ))')
     &                   SEGID(K),RESID(K),RES(K),TYPE(K),
     &                   SEGID(P),RESID(P),RES(P),TYPE(P),
     &                   SEGID(T),RESID(T),RES(T),TYPE(T),
     &                   CUREQUIL, CURANGLE, CURDELTA,
     &                   CURPHI, CURPSI, CURENER
               END IF
          END IF
      END DO
C
C handle summary info
C
      IF (NUMINRMS.NE.ZERO) THEN
           RMS = SQRT(RMS / NUMINRMS)
      ELSE
           RMS = ZERO
      END IF
      WRITE (PUNIT, '(A, F8.3)') 'RMS deviation (deg) = ', RMS
      WRITE (PUNIT, '(A, F6.1)') 'number of viols = ', VIOLS
      DBPREC=NUMINRMS
      CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, DBPREC)
      CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
      CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      RETURN
      END
C
C ============
C
      SUBROUTINE SCRANGDB
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'angledb.inc'
C
C
C reset angle database
      IF (.NOT.ANGLEDBFLAG) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: angle database erased.'
      CALL ANGLEDBFREE
      CALL ANGLEDBINIT
      END IF
      RETURN
      END
C
