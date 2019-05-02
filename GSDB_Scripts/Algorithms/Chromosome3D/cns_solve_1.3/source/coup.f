C===============
      SUBROUTINE ECOUP (EJ, WHICH)
C
C Calls ECOUP2, which does the actual energy calculation
C
C by John Kuszewski Aug 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION EJ
      CHARACTER*7 WHICH
C begin
C
      CALL ECOUP2(EJ, HEAP(COUPIPTR), HEAP(COUPJPTR), HEAP(COUPKPTR),
     &     HEAP(COUPLPTR), HEAP(COUPMPTR), HEAP(COUPNPTR),
     &     HEAP(COUPPPTR), HEAP(COUPQPTR), HEAP(COUPJOBSPTR),
     &     HEAP(COUPJERRPTR), HEAP(CALCCOUPPTR), WHICH, HEAP(COUPCV))
      RETURN
      END
C===============
      SUBROUTINE ECOUP2 (EJ, ATOMI, ATOMJ, ATOMK, ATOML,
     &     ATOMM, ATOMN, ATOMP, ATOMQ, JOBS, JERR, JCALC, WHICH, JCV)
C
C Calculates coupling constant energies
C
C J energies are of the form
C      E = k1*deltaJ**2
C where
C      Kcoup = main force constant,
C      newK = new force constant
C      deltaJ = calculated J - observed J
C
C which is a flag that switches between energy & force and
C calculated J (for violations) calcs
C
C by John Kuszewski July 1993
C cross-validation by Alexandre Bonvin Dec 1995
C
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     ATOMM(*), ATOMN(*), ATOMP(*), ATOMQ(*), JCV(*)
      DOUBLE PRECISION JOBS(2,*), JERR(2,*), JCALC(2,*)
      DOUBLE PRECISION EJ
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, CLASS, POTENTIAL
      DOUBLE PRECISION XI, XJ, XK, XL, XM, XN, XP, XQ,
     &     YI, YJ, YK, YL, YM, YN, YP, YQ, ZI, ZJ, ZK, ZL,
     &     ZM, ZN, ZP, ZQ,
     &     XIJ, XJK, XKL, XMN, XNP, XPQ,
     &     YIJ, YJK, YKL, YMN, YNP, YPQ,
     &     ZIJ, ZJK, ZKL, ZMN, ZNP, ZPQ,
     &     AX, AY, AZ, AX2, AY2, AZ2,
     &     BX, BY, BZ, BX2, BY2, BZ2,
     &     CX, CY, CZ, CX2, CY2, CZ2,
     &     RAR, RBR, RCR, RAR2, RBR2, RCR2,
     &     CP, SP, CP2, SP2, THETAA, THETAB,
     &     o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13,
     &     o14, o15, o16, o17, o18, o19, o20, o21, o22, o23, o24,
     &     o25, o26, o27, o28,
     &     JcalcA, JcalcB, testterm1, testterm2, testterm3,
     &     Eplus, DEplusDthetaA, DEplusDthetaB,
     &     Eminus1, DEminus1DthetaA, DEminus1DthetaB,
     &     Eminus2, DEminus2DthetaA, DEminus2DthetaB,
     &     Eminus3, DEminus3DthetaA, DEminus3DthetaB,
     &     DEDthetaA, DEDthetaB, firstviol, secondviol,
     &     E, DF, DF2, SWITCH, SWITCH2,
     &     RECSP, RECCP, RECSP2, RECCP2,
     &     DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ,
     &     DCPAX2, DCPAY2, DCPAZ2, DCPBX2, DCPBY2, DCPBZ2,
     &     DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ,
     &     DSPCX2, DSPCY2, DSPCZ2, DSPBX2, DSPBY2, DSPBZ2,
     &     DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ,
     &     DPRIJX2, DPRIJY2, DPRIJZ2, DPRKLX2, DPRKLY2, DPRKLZ2,
     &     DPRJKX, DPRJKY, DPRJKZ, DPRJKX2, DPRJKY2, DPRJKZ2,
     &     A, B, C, PHASE
      DOUBLE PRECISION JOBS1, JOBS2, ERRJ, Kcoup, newK
C begin
C
C following Axel's code in ETOR,
C
C zero out partial energy
C
      EJ = ZERO
C
      CLASS = 1
      Kcoup = COUPFORCES(1,CLASS)
      newK  = COUPFORCES(2,CLASS)
      A = COUPAS(CLASS)
      B = COUPBS(CLASS)
      C = COUPCS(CLASS)
      PHASE = COUPPHASES(CLASS)
      POTENTIAL = COUPPOTENTIAL(CLASS)
      DO COUNT = 1, NCOUPS
          IF (COUPASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               Kcoup = COUPFORCES(1,CLASS)
               newK  = COUPFORCES(2,CLASS)
               A = COUPAS(CLASS)
               B = COUPBS(CLASS)
               C = COUPCS(CLASS)
               PHASE = COUPPHASES(CLASS)
               POTENTIAL = COUPPOTENTIAL(CLASS)
          END IF
C
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
          XN = X(ATOMN(COUNT))
          XP = X(ATOMP(COUNT))
          XQ = X(ATOMQ(COUNT))
C
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
          YN = Y(ATOMN(COUNT))
          YP = Y(ATOMP(COUNT))
          YQ = Y(ATOMQ(COUNT))
C
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
          ZN = Z(ATOMN(COUNT))
          ZP = Z(ATOMP(COUNT))
          ZQ = Z(ATOMQ(COUNT))
C
          JOBS1 = JOBS(1,COUNT)
          JOBS2 = JOBS(2,COUNT)
          ERRJ = JERR(1,COUNT)
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
C calculate cos(theta) & sin(theta)
C
          CP = AX*BX+AY*BY+AZ*BZ
          SP = CX*BX+CY*BY+CZ*BZ
C
C calculate cos(theta) & sin(theta)
C
          CP2 = AX2*BX2+AY2*BY2+AZ2*BZ2
          SP2 = CX2*BX2+CY2*BY2+CZ2*BZ2
C
C calculate theta (make sure cos is within bounds and get sign from
C sin) and keep it it radians
C
          thetaA = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
          thetaB = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP2))),SP2)
C
C calculate energy and forces for normal entries
C
C Note that CALCJ = A*COS(PHI+PHASE)**2 + B*COS(PHI+PHASE) + C
C
          IF (POTENTIAL.EQ.SQUARE) THEN
               o1=phase+thetaA
               o2=cos(o1)
               o3=B*o2
               o4=o2**2
               o5=A*o4
               o6=C-JOBS1+o3+o5
               o7=kcoup*(-ERRJ**2+o6**2)
               JcalcA = C+o3+o5
               if (o7.lt.zero) then
                    E = zero
                    DEDthetaA = zero
                    DEDthetaB = zero
               else
                    E = o7
                    DEDthetaA = -2.d0*kcoup*(B+2.d0*A*o2)*
     &                   o6*sin(o1)
                    DEDthetaB = zero
               end if
           ELSE IF (POTENTIAL.EQ.HARMONIC) THEN
               o1=phase+thetaA
               o2=cos(o1)
               o3=o2**2
               o4=A*o3
               o5=phase+thetaA
               o6=cos(o5)
               o7=B*o6
               o8=C-Jobs1+o4+o7
               JcalcA = C+o4+o7
               E = kcoup*o8**2
               DEDthetaA = 2.d0*kcoup*o8*
     &              (-2.d0*A*o2*sin(o1)-B*sin(o5))
               DEDthetaB = zero
          ELSE
C
C Calculate Constantine energy (harmonic only)
C
C All this is from Mma
C
               o1 = phase + thetaA
               o2 = Cos(o1)
               o3 = B*o2
               o4 = o2**2
               o5 = A*o4
               o6 = phase + thetaB
               o7 = Cos(o6)
               o8 = B*o7
               o9 = o7**2
               o10 = A*o9
               o11 = -o10 + o3 + o5 - o8
               o12 = o11**2
               o13 = Sqrt(o12)
               o14 = Jobs1 - Jobs2
               o15 = o14**2
               o16 = Sqrt(o15)
               o17 = 2*C - Jobs1 - Jobs2 + o10 + o3 + o5 + o8
               o18 = Sin(o1)
               o19 = B*o18
               o20 = A*o18*o2
               o21 = -o19 - 2*o20
               o22 = Sin(o6)
               o23 = B*o22
               o24 = A*o22*o7
               o25 = o13 - o16
               o26 = 1/o13
               o27 = o23 + 2*o24
               o28 = -o13 + o16
C
               JcalcA = C + o3 + o5
               JcalcB = C + o10 + o8
               testterm1 = o13
               testterm2 = o16
               testterm3 = o16/2
               Eplus = kcoup*o17**2
               DEplusDthetaA = 2*kcoup*o17*o21
               DEplusDthetaB = 2*kcoup*o17*(-o23 - 2*o24)
               Eminus1 = kcoup*o25**2
               DEminus1DthetaA = 2*kcoup*o11*o21*o25*o26
               DEminus1DthetaB = 2*kcoup*o11*o25*o26*o27
               Eminus2 = newK*kcoup*o28**2
               DEminus2DthetaA = -2*newK*kcoup*o11*o21*o26*o28
               DEminus2DthetaB = -2*newK*kcoup*o11*o26*o27*o28
               Eminus3 = newK*kcoup*(-o12 + o15/2)
               DEminus3DthetaA = -2*newK*kcoup*o11*o21
               DEminus3DthetaB = -2*newK*kcoup*o11*o27
C
C Apply the rules from the Constantine paper
C Note that there is no square well potential
C for this
C
               if (testterm1.gt.testterm2) then
                    E = Eplus + Eminus1
                    DEDthetaA = DEplusDthetaA + DEminus1DthetaA
                    DEDthetaB = DEplusDthetaB + DEminus1DthetaB
               else if (testterm1.ge.testterm3) then
                    E = Eplus + Eminus2
                    DEDthetaA = DEplusDthetaA + DEminus2DthetaA
                    DEDthetaB = DEplusDthetaB + DEminus2DthetaB
               else
                    E = Eplus + Eminus3
                    DEDthetaA = DEplusDthetaA + DEminus3DthetaA
                    DEDthetaB = DEplusDthetaB + DEminus3DthetaB
               end if
          END IF
C
C If we're printing out the couplings, then we need to
C try both possible assignments of observed with the
C the angle and pick the one with the lower total violation.
C
C For non-MULTIPLE couplings, the calculated J goes
C in the first slot
C
          IF (WHICH.EQ.'ANALYZE') then
               IF (POTENTIAL.NE.MULTIPLE) THEN
                    JCALC(1,COUNT) = JcalcA
                    JCALC(2,COUNT) = 0
               else
                    firstviol  = ABS(Jobs1 - JcalcA)**2 +
     &                   ABS(Jobs2 - JcalcB)**2
                    secondviol = ABS(Jobs1 - JcalcB)**2 +
     &                   ABS(Jobs2 - JcalcA)**2
                    if (firstviol.le.secondviol) then
                         Jcalc(1,count) = JcalcA
                         Jcalc(2,count) = JcalcB
                    else
                         Jcalc(1,count) = JcalcB
                         Jcalc(2,count) = JcalcA
                    end if
               end if
          end if
C
C accumulate energy (only for working set)
C
          IF (JCV(COUNT).NE.JICV) THEN
            EJ = EJ + E
          END IF
          DF = DEDthetaA
          DF2 = DEDthetaB
C
C compute heavyside function
C
          SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
          SWITCH2=-MIN(1,INT(ABS(SP2)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
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
          DPRJKZ=
     &     RECSP*(DCPAX*YIJ-DCPAY*XIJ+DCPBY*XKL-DCPBX*YKL)+
     &     RECCP*(-(XJK*XIJ+YJK*YIJ)*DSPCZ
     &        +(TWO*ZJK*XIJ-ZIJ*XJK)*DSPCX
     &        +(TWO*ZJK*YIJ-ZIJ*YJK)*DSPCY
     &        +DSPBY*XKL-DSPBX*YKL)
C
C now update forces if in energy & force mode
C
          IF (WHICH.NE.'ANALYZE') THEN
            IF (JCV(COUNT).NE.JICV) THEN
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
          END IF
      END DO
      RETURN
      END
C================
      SUBROUTINE READCOUP
C
C reads in coupling constant information
C
C by John Kuszewski July 1993
C cross-validation by Alexandre Bonvin Dec 1995
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'ctitla.inc'
C i/o
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXCOUPS, TEMP
      DOUBLE PRECISION K1, K2, CUTOFF, A, B, C, P
      CHARACTER*4 THENAME
C begin
C
C this is used by READCOUP2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
C reset database only if no couplings have been entered
C ie., this is the first time in the script that
C COUPlings has appeared
C
      IF (COUPIPTR.EQ.0) THEN
           CALL ALLOCCOUPS(0, MAXCOUPS)
      END IF
C
C now read input
C
      CALL PUSEND('COUPLINGS>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('COUPLINGS>')
           CALL MISCOM('COUPLINGS>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-coupling')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURCLASS
                CALL NEXTA4('class name =', THENAME)
                MODE = NEW
                DO COUNT = 1, NCLASSES
                     IF (COUPCLASSNAMES(COUNT).EQ.THENAME) THEN
                          MODE = UPDATE
                          CURCLASS = COUNT
                     END IF
                END DO
                IF (MODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXCOUPCLASSES) THEN
                          CALL DSPERR('COUP','Too many classes.')
                          CALL DSPERR('COUP',
     &                      'Increase MAXCOUPCLASSES and recompile.')
                          CALL WRNDIE(-5, 'READCOUP',
     &                                'Too many J-coupling classes.')
                     END IF
                     NCLASSES = NCLASSES + 1
                     CURCLASS = NCLASSES
                     COUPCLASSNAMES(CURCLASS) = THENAME
C
C If this isn't the first class, close off the old class
C
                IF (NCLASSES.GT.1) THEN
                     COUPASSNDX(OLDCLASS) = NCOUPS
                END IF
              END IF
C
C set Karplus curve coefficients for this class
C
           ELSE IF (WD(1:4).EQ.'COEF') THEN
                CALL NEXTF('coefficient A =', A)
                CALL NEXTF('coefficient B =', B)
                CALL NEXTF('coefficient C =', C)
                CALL NEXTF('phase =', P)
C
C start a default class if there isn't one defined
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                WRITE(PUNIT, '(A, A, A, F8.3, F8.3, F8.3)')
     &             'Setting Karplus coefficients for class ',
     &             COUPCLASSNAMES(CURCLASS), ' to ', A, B, C
                WRITE(PUNIT, '(A, A, A, F8.3)')
     &             'And setting Karplus phase (degrees) for class ',
     &             COUPCLASSNAMES(CURCLASS), ' to ', P
C
C switch the value to radians
C
                CALL DG2RAD(P, P)
                COUPPHASES(CURCLASS) = P
                COUPAS(CURCLASS) = A
                COUPBS(CURCLASS) = B
                COUPCS(CURCLASS) = C
C
C set force constant for current class
C
           ELSE IF (WD(1:4).EQ.'FORC') THEN
                CALL NEXTF('force constant =', K1)
                IF (COUPPOTENTIAL(CURCLASS).EQ.MULTIPLE) THEN
                     CALL NEXTF('second force constant =', K2)
                ELSE
                     K2 = ONE
                END IF
C
C start a default class if there isn't one defined
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                IF (COUPPOTENTIAL(CURCLASS).NE.MULTIPLE) THEN
                     WRITE(PUNIT, '(A, A, A, F8.3)')
     &                    'Setting force const for class ',
     &                    COUPCLASSNAMES(CURCLASS), ' to ', K1
                ELSE
                     WRITE(PUNIT, '(A, A, A, F8.3, F8.3)')
     &                    'Setting force consts for class ',
     &                    COUPCLASSNAMES(CURCLASS), ' to ', K1, K2
                END IF
                COUPFORCES(1,CURCLASS) = K1
                COUPFORCES(2,CURCLASS) = K2
C
C reset couplings database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
C
                CALL COUPHP
                CALL COUPINIT
                CALL ALLOCCOUPS(0, MAXCOUPS)
C
C
C set potential type for current class
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                     WRITE(PUNIT, '(A)')
     &                'using square well potential in current class.'
                     COUPPOTENTIAL(CURCLASS) = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                     WRITE(PUNIT, '(A)')
     &                'using harmonic potential in current class.'
                     COUPPOTENTIAL(CURCLASS) = HARMONIC
                ELSE IF (THENAME.EQ.'MULT') THEN
                     WRITE(PUNIT, '(A)')
     &       'using MULTIPLE harmonic potential in current class.'
                     COUPPOTENTIAL(CURCLASS) = MULTIPLE
                ELSE
                     CALL DSPERR('COUP',
     &       'unknown potential. Using square well for current class.')
                     COUPPOTENTIAL(CURCLASS) = SQUARE
                END IF
C
C
C change number of assignment slots
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXCOUPS = MAXCOUPS
                CALL NEXTI('number of slots =', MAXCOUPS)
                CALL ALLOCCOUPS(OLDMAXCOUPS, MAXCOUPS)
C
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more coupling assignments
C than you have slots for
C
                IF (NCOUPS.EQ.MAXCOUPS) THEN
                     CALL DSPERR('COUP','Too many assignments.')
                     CALL DSPERR('COUP',
     &                    'Increase NREStraints and run again.')
                     CALL WRNDIE(-1,'COUP>',
     &                    'exceeded allocation for J-restraints')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                CALL READCOUP2(HEAP(COUPIPTR), HEAP(COUPJPTR),
     &               HEAP(COUPKPTR), HEAP(COUPLPTR),
     &               HEAP(COUPMPTR), HEAP(COUPNPTR),
     &               HEAP(COUPPPTR), HEAP(COUPQPTR),
     &               HEAP(COUPJOBSPTR), HEAP(COUPJERRPTR),
     &               HEAP(SPTR), HEAP(COUPCV))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                     CALL DSPERR('COUPLINGS',
     &                           'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                          CALL DSPERR('COUPLINGS',
     &                                'cutoff must be positive.')
                          CUTOFF = ABS(CUTOFF)
                     END IF
                     CALL NEXTA4('ALL or CLASs>', THENAME)
                     IF (THENAME(1:3).EQ.'ALL') THEN
                          DO COUNT = 1,NCLASSES
                               PRINTCLASS(COUNT) = .TRUE.
                          END DO
                     ELSE IF (THENAME(1:4).EQ.'CLAS') THEN
                          CALL NEXTA4('class name =', THENAME)
                          DO COUNT = 1,NCLASSES
                               IF (COUPCLASSNAMES(COUNT).EQ.
     &                              THENAME) THEN
                                    PRINTCLASS(COUNT) = .TRUE.
                               ELSE
                                    PRINTCLASS(COUNT) = .FALSE.
                               END IF
                          END DO
                     ELSE
                          CALL DSPERR('COUPLINGS',
     &                         'not understood.  Printing all classes')
                          DO COUNT = 1,NCLASSES
                               PRINTCLASS(COUNT) = .TRUE.
                          END DO
                     END IF
C
                     CALL PRINTCOUPS(CUTOFF, HEAP(CALCCOUPPTR),
     &                    HEAP(COUPJOBSPTR),
     &                    HEAP(COUPIPTR), HEAP(COUPJPTR),
     &                    HEAP(COUPKPTR), HEAP(COUPLPTR),
     &                    HEAP(COUPMPTR), HEAP(COUPNPTR),
     &                    HEAP(COUPPPTR), HEAP(COUPQPTR),
     &                    HEAP(COUPCV), 0)
                     IF (JICV.GT.0) THEN
                       CALL PRINTCOUPS(CUTOFF, HEAP(CALCCOUPPTR),
     &                      HEAP(COUPJOBSPTR),
     &                      HEAP(COUPIPTR), HEAP(COUPJPTR),
     &                      HEAP(COUPKPTR), HEAP(COUPLPTR),
     &                      HEAP(COUPMPTR), HEAP(COUPNPTR),
     &                      HEAP(COUPPPTR), HEAP(COUPQPTR),
     &                      HEAP(COUPCV), 1)
                     END IF
                END IF
C====================================================================
      ELSE IF (WD(1:2).EQ.'CV') THEN
      CALL NEXTI('CV excluded partition number:',JICV)
C====================================================================
      ELSE IF (WD(1:4).EQ.'PART') THEN
      CALL NEXTI('number of PARTitions:',TEMP)
      TEMP=MAX(0,TEMP)
      CALL JCVS(NCOUPS,HEAP(COUPCV),TEMP)
      IF (TEMP.EQ.0) THEN
      JICV=0
      END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('COUPLINGS>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCCOUPS (OLDSIZE, NEWSIZE)
C
C resets coupling constant arrays to hold SIZE entries
C
C by John Kuszewski July 1993
C cross-validation by Alexandre Bonvin Dec 1995
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'funct.inc'
      INCLUDE 'couplings.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(COUPIPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPJPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPKPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPLPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPMPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPNPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPPPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPQPTR, INTEG4(OLDSIZE))
           CALL FREHP(COUPCV  , INTEG4(OLDSIZE))
           CALL FREHP(COUPJOBSPTR, IREAL8(2*OLDSIZE))
           CALL FREHP(COUPJERRPTR, IREAL8(2*OLDSIZE))
           CALL FREHP(CALCCOUPPTR, IREAL8(2*OLDSIZE))
      END IF
C
      COUPIPTR = 0
      COUPJPTR = 0
      COUPKPTR = 0
      COUPLPTR = 0
      COUPMPTR = 0
      COUPNPTR = 0
      COUPPPTR = 0
      COUPQPTR = 0
      COUPCV   = 0
      COUPJOBSPTR = 0
      COUPJERRPTR = 0
      CALCCOUPPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        COUPIPTR = ALLHP(INTEG4(NEWSIZE))
        COUPJPTR = ALLHP(INTEG4(NEWSIZE))
        COUPKPTR = ALLHP(INTEG4(NEWSIZE))
        COUPLPTR = ALLHP(INTEG4(NEWSIZE))
        COUPMPTR = ALLHP(INTEG4(NEWSIZE))
        COUPNPTR = ALLHP(INTEG4(NEWSIZE))
        COUPPPTR = ALLHP(INTEG4(NEWSIZE))
        COUPQPTR = ALLHP(INTEG4(NEWSIZE))
        COUPCV   = ALLHP(INTEG4(NEWSIZE))
        COUPJOBSPTR = ALLHP(IREAL8(2*NEWSIZE))
        COUPJERRPTR = ALLHP(IREAL8(2*NEWSIZE))
        CALCCOUPPTR = ALLHP(IREAL8(2*NEWSIZE))
      END IF
C
      RETURN
      END
C==============
      SUBROUTINE COUPDEFAULTS
C
C sets up defaults
C
C by John Kuszewski July 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'couplings.inc'
C local variables
      INTEGER COUNT
      DOUBLE PRECISION P
C begin
      MODE = NEW
      MAXCOUPS = 200
      NCOUPS = 0
      NCLASSES = 0
      CURCLASS = 0
      JICV = 0
      COUPIPTR = 0
      COUPJPTR = 0
      COUPKPTR = 0
      COUPLPTR = 0
      COUPMPTR = 0
      COUPNPTR = 0
      COUPPPTR = 0
      COUPQPTR = 0
      COUPCV   = 0
      COUPJOBSPTR = 0
      COUPJERRPTR = 0
      CALCCOUPPTR = 0
      DO COUNT = 1, MAXCOUPCLASSES
           COUPCLASSNAMES(COUNT) = 'DEFAULT'
           COUPASSNDX(COUNT) = 0
           COUPFORCES (1,COUNT) = 50
           COUPFORCES (2,COUNT) = 10
           COUPPOTENTIAL(COUNT) = HARMONIC
C
C these values are for the HnHa coupling
C Wang and Bax, JACS 118, 2483-2494 (1996)
C
           COUPAS(COUNT) = 6.98
           COUPBS(COUNT) = -1.38
           COUPCS(COUNT) = 1.72
C
C get the value in radians
C
           CALL DG2RAD(-60.0D0, P)
           COUPPHASES(COUNT) = P
      END DO
      RETURN
      END
C==============
      SUBROUTINE READCOUP2 (ATOMI, ATOMJ, ATOMK, ATOML,
     &     ATOMM, ATOMN, ATOMP, ATOMQ, JOBS, JERR, SEL, JCV)
C
C reads actual J coupling assignments into arrays
C
C by John Kuszewski July 1993
C cross-validation by Alexandre Bonvin Dec 1995
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     ATOMM(*), ATOMN(*), ATOMP(*), ATOMQ(*), SEL(*)
      INTEGER JCV(*)
      DOUBLE PRECISION JOBS(2,*), JERR(2,*)
C local variables
      INTEGER NSEL, INSERTPOS, COUNT, CURSTOP, OTHERSTOP
      DOUBLE PRECISION JO, JE
C begin
C
C if we're in update mode, make a space for the new line
C
      IF (MODE.EQ.UPDATE) THEN
           DO COUNT = NCOUPS+1, COUPASSNDX(CURCLASS)+1, -1
                ATOMI(COUNT) = ATOMI(COUNT-1)
                ATOMJ(COUNT) = ATOMJ(COUNT-1)
                ATOMK(COUNT) = ATOMK(COUNT-1)
                ATOML(COUNT) = ATOML(COUNT-1)
                ATOMM(COUNT) = ATOMM(COUNT-1)
                ATOMN(COUNT) = ATOMN(COUNT-1)
                ATOMP(COUNT) = ATOMP(COUNT-1)
                ATOMQ(COUNT) = ATOMQ(COUNT-1)
                JOBS(1,COUNT) = JOBS(1,COUNT-1)
                JOBS(2,COUNT) = JOBS(2,COUNT-1)
                JERR(1,COUNT) = JERR(1,COUNT-1)
                JERR(2,COUNT) = JERR(2,COUNT-1)
                JCV(COUNT) = JCV(COUNT-1)
           END DO
           CURSTOP = COUPASSNDX(CURCLASS)
           DO COUNT = 1, NCLASSES
                OTHERSTOP = COUPASSNDX(COUNT)
                IF (OTHERSTOP.GT.CURSTOP) THEN
                     COUPASSNDX(COUNT) = OTHERSTOP + 1
                END IF
           END DO
           COUPASSNDX(CURCLASS) = CURSTOP + 1
           INSERTPOS = CURSTOP
           NCOUPS = NCOUPS + 1
      ELSE
           NCOUPS = NCOUPS + 1
           INSERTPOS = NCOUPS
           COUPASSNDX(CURCLASS) = INSERTPOS
      END IF
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom i. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMI(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom j. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMJ(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom k. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMK(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom l. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOML(INSERTPOS) = SEL(1)
C
      IF (COUPPOTENTIAL(CURCLASS).EQ.MULTIPLE) THEN
C
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom m. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           ATOMM(INSERTPOS) = SEL(1)
C
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom n. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           ATOMN(INSERTPOS) = SEL(1)
C
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom p. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           ATOMP(INSERTPOS) = SEL(1)
C
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('COUP',
     &     'more than 1 atom in selection for atom q. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           ATOMQ(INSERTPOS) = SEL(1)
C
      ELSE
           ATOMM(INSERTPOS) = ATOMI(INSERTPOS)
           ATOMN(INSERTPOS) = ATOMJ(INSERTPOS)
           ATOMP(INSERTPOS) = ATOMK(INSERTPOS)
           ATOMQ(INSERTPOS) = ATOML(INSERTPOS)
      END IF
C
      CALL NEXTF('observed J =', JO)
      CALL NEXTF('error in J =', JE)
      JOBS(1,INSERTPOS) = JO
      JERR(1,INSERTPOS) = JE
      IF (COUPPOTENTIAL(CURCLASS).EQ.MULTIPLE) THEN
           CALL NEXTF('observed J =', JO)
           CALL NEXTF('error in J =', JE)
           JOBS(2,INSERTPOS) = JO
           JERR(2,INSERTPOS) = JE
      ELSE
           JOBS(2,INSERTPOS) = JOBS(1,INSERTPOS)
           JERR(2,INSERTPOS) = JERR(1,INSERTPOS)
      END IF
      JCV(INSERTPOS)  = -1
      RETURN
      END
C===============
      SUBROUTINE COUPINIT
C
C initializes J-coupling stuff
C
C by John Kuszewski Aug 1993
C updated for MULTIPLE couplings Feb 1996
C================
      IMPLICIT NONE
C include files
      INCLUDE 'couplings.inc'
C begin
      CALL COUPDEFAULTS
      RETURN
      END
C===============
      SUBROUTINE COUPHP
C
C deallocates J-coupling stuff
C
C by Gregory Warren 6/20/95
C================
      IMPLICIT NONE
C include files
      INCLUDE 'couplings.inc'
C begin
      IF (COUPIPTR.NE.0) THEN
      CALL ALLOCCOUPS(MAXCOUPS,0)
      END IF
      RETURN
      END
C=================
      SUBROUTINE PRINTCOUPS (CUTOFF, JCALC, JOBS,
     &     ATOMI, ATOMJ, ATOMK, ATOML, ATOMM, ATOMN, ATOMP, ATOMQ,
     &     JCV, ITEST)
C
C prints couplings with deltaJ greater than cutoff
C calculates RMS deviation and puts it into $RESULT
C
C by John Kuszewski Aug 1993
C cross-validation by Alexandre Bonvin Dec 1995
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'couplings.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, JCALC(2,*), JOBS(2,*)
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*),
     &     ATOMN(*), ATOMP(*), ATOMQ(*)
      INTEGER JCV(*), ITEST
C local variables
      DOUBLE PRECISION CALCJ, CALCJ1, CALCJ2, OBSJ, OBSJ1, OBSJ2,
     &     DELTAJ, DELTAJ1, DELTAJ2, JENERGY, RMS, VIOLS, DBPREC
      INTEGER COUNT, CLASS, I, J, K, L, M, N, P, Q, CURDEGENERACY,
     &     NASSIGNSINCLUDED
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS
C begin
      RMS = ZERO
      VIOLS = ZERO
      NASSIGNSINCLUDED = ZERO
C
C make sure that the calcJ array is up to date
C
      CALL ECOUP(JENERGY, 'ANALYZE')
      IF (JICV.GT.0) THEN
        IF (ITEST.EQ.0) THEN
           WRITE(PUNIT,'(A)')
     & ' $$$$$$$$$$$$$$$$$$ working set $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        ELSE
           WRITE(PUNIT,'(A,I5,A)')
     & ' $$$$$$$$$$$$$$$$$$$$ test set (TEST=',JICV,
     & ')  $$$$$$$$$$$$$$$$$$$$$$'
        END IF
      END IF
      WRITE (PUNIT, '(A)') 'The following couplings have delta J'
      WRITE (PUNIT, '(A)') 'greater than the cutoff:'
      WRITE (PUNIT, '(A)') '  (calculated J) (observed J) (delta J)'
C
C write out first class heading
C
      CLASS = 1
      CURDEGENERACY = COUPPOTENTIAL(CLASS)
      PRINTTHISCLASS = PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
           WRITE (PUNIT, '(A, A)') 'class ', COUPCLASSNAMES(1)
      END IF
C
C for every coupling entry,
C
      DO COUNT = 1, NCOUPS
C
C is this the start of a new class?
C
        IF (COUPASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               PRINTTHISCLASS = PRINTCLASS(CLASS)
               CURDEGENERACY = COUPPOTENTIAL(CLASS)
               IF (PRINTTHISCLASS) THEN
                WRITE (PUNIT, '(A, A)') 'class ', COUPCLASSNAMES(CLASS)
               END IF
        END IF
C
C check if in test set or not (cross-validation)
C
        IF ((ITEST.EQ.0.AND.JCV(COUNT).NE.JICV).OR.
     &      (ITEST.EQ.1.AND.JCV(COUNT).EQ.JICV)) THEN
C
C if this ASSIgnment is in a class that will be printed,
C
          IF (PRINTTHISCLASS) THEN
C
C if this is a normal coupling,
C
               IF (CURDEGENERACY.NE.MULTIPLE) THEN
C
C update RMS delta J
C
                    CALCJ = JCALC(1,COUNT)
                    OBSJ = JOBS(1,COUNT)
                    DELTAJ = CALCJ-OBSJ
                    RMS = RMS + DELTAJ**2
                    NASSIGNSINCLUDED = NASSIGNSINCLUDED + 1
C
C print out delta Js greater than cutoff
C and update number of violations
C
                    IF (ABS(DELTAJ).GT.CUTOFF) THEN
                         I = ATOMI(COUNT)
                         J = ATOMJ(COUNT)
                         K = ATOMK(COUNT)
                         L = ATOML(COUNT)
                         WRITE(PUNIT,'(9X,16(1X,A), 3(F8.3))')
     &                        SEGID(I),RESID(I),RES(I),TYPE(I),
     &                        SEGID(J),RESID(J),RES(J),TYPE(J),
     &                        SEGID(K),RESID(K),RES(K),TYPE(K),
     &                        SEGID(L),RESID(L),RES(L),TYPE(L),
     &                        CALCJ, OBSJ, DELTAJ
                         VIOLS = VIOLS + ONE
                    END IF
C
C if this is a MULTIPLE coupling,
C
               ELSE
C
C update RMS delta J
C
                    CAlCJ1 = JCALC(1,COUNT)
                    OBSJ1 = JOBS(1,COUNT)
                    DELTAJ1 = CALCJ1 - OBSJ1
                    RMS = RMS + DELTAJ1**2
                    CALCJ2 = JCALC(2,COUNT)
                    OBSJ2 = JOBS(2,COUNT)
                    DELTAJ2 = CALCJ2 - OBSJ2
                    RMS = RMS + DELTAJ2**2
                    NASSIGNSINCLUDED = NASSIGNSINCLUDED + 2
C
C print out entries where either deltaJ is
C greater than the cutoff
C
                    IF (((ABS(DELTAJ1).GT.CUTOFF).OR.
     &                   (ABS(DELTAJ2).GT.CUTOFF))
     &                   .AND.PRINTTHISCLASS) THEN
                         I = ATOMI(COUNT)
                         J = ATOMJ(COUNT)
                         K = ATOMK(COUNT)
                         L = ATOML(COUNT)
                         M = ATOMM(COUNT)
                         N = ATOMN(COUNT)
                         P = ATOMP(COUNT)
                         Q = ATOMQ(COUNT)
                         WRITE(PUNIT,'(9X,16(1X,A), 3(F8.3))')
     &                        SEGID(I),RESID(I),RES(I),TYPE(I),
     &                        SEGID(J),RESID(J),RES(J),TYPE(J),
     &                        SEGID(K),RESID(K),RES(K),TYPE(K),
     &                        SEGID(L),RESID(L),RES(L),TYPE(L),
     &                        CALCJ1, OBSJ1, DELTAJ1
                         WRITE(PUNIT,'(9X,16(1X,A), 3(F8.3))')
     &                        SEGID(M),RESID(M),RES(M),TYPE(M),
     &                        SEGID(N),RESID(N),RES(N),TYPE(N),
     &                        SEGID(P),RESID(P),RES(P),TYPE(P),
     &                        SEGID(Q),RESID(Q),RES(Q),TYPE(Q),
     &                        CALCJ2, OBSJ2, DELTAJ2
                         VIOLS = VIOLS + ONE
                    END IF
               END IF
          END IF
        END IF
      END DO
C
      IF (NASSIGNSINCLUDED.GT.ZERO) THEN
           RMS = SQRT(RMS / NASSIGNSINCLUDED)
      ELSE
           RMS = ZERO
      END IF
      WRITE(PUNIT,'(A,F8.3,A,F5.2,A,F6.0,A,I6,A)')
     & '  RMS diff. =',RMS,
     & ', #(violat.>',CUTOFF,')=',VIOLS,
     & ' of ',NASSIGNSINCLUDED,' J-couplings'
      IF (ITEST.EQ.1) THEN
        CALL DECLAR('TEST_RMS', 'DP', ' ', DUMMY2, RMS)
        CALL DECLAR('TEST_VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      ELSE
        DBPREC = NASSIGNSINCLUDED
        CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, DBPREC)
        CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
        CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, RMS)
        CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      END IF
      RETURN
      END
C
C====================================================================
      SUBROUTINE JCVS(JNUM,JCV,PART)
C
C Routine partitions J-coupling data into PART sets.
C JCV will contain integer numbers between 1 and PART.
C
C Author: Axel T. Brunger
C Modifed for J-couplings: Alexandre Bonvin 12/22/95
C
C
C     IMPLICIT NONE
C I/O
      INTEGER JNUM, JCV(*)
      INTEGER PART
C local
      INTEGER I, P, NP, NRETRY, NTRYTOT
      DOUBLE PRECISION RNUM
      NTRYTOT = 0
C begin
  100 CONTINUE
      NRETRY = 0
      IF (PART.GT.0) THEN
      DO I=1,JNUM
      CALL GGUBFS(RNUM)
      JCV(I)=MAX(1,MIN(PART,INT(RNUM*PART)+1))
      END DO
C
      IF (PART .EQ. JNUM) THEN
      DO I=1,JNUM
      JCV(I)=I
      END DO
      END IF
C
      DO P=1,PART
      NP=0
      DO I=1,JNUM
      IF (JCV(I).EQ.P) THEN
      NP=NP+1
      END IF
      END DO
      IF (NP .EQ. 0) THEN
        NRETRY = 1
      ENDIF
      WRITE(6,'(A,I3,A,I5,A)') ' For set ',P,
     & ' there are ',NP,' J-coupling angle restraints.'
      END DO
      ELSE
      WRITE(6,'(A)')
     & ' Data are not partitioned or partitioning removed.'
      DO I=1,JNUM
      JCV(I)=-1
      END DO
      END IF
      NTRYTOT = NTRYTOT + NRETRY
      IF (NRETRY .GT. 0 .AND. NTRYTOT .LE. 10) THEN
      WRITE(6,'(A)')
     & ' Test set with 0 constraints! New trial...'
      GOTO 100
      ELSE IF (NTRYTOT .GT. 10) THEN
      CALL WRNDIE(-1,'JCVS',
     & 'Unable to partition the J-coupling data within ten trials')
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE SCRCUPP
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'couplings.inc'
C
C
C update J-coupling database
      IF (COUPIPTR.NE.0) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: j-coupling database erased.'
      CALL COUPHP
      CALL COUPINIT
      END IF
      RETURN
      END
C
