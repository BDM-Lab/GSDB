C===============
      SUBROUTINE ECSHIFT (EC, WHICH)
C
C Calls ECSHIFT2, which does the actual energy calculation
C
C by John Kuszewski Aug 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION EC
      CHARACTER*7 WHICH
C begin
C
      CALL ECSHIFT2(EC, HEAP(SHIFTIPTR), HEAP(SHIFTJPTR),
     &              HEAP(SHIFTKPTR), HEAP(SHIFTLPTR),
     &              HEAP(SHIFTMPTR), HEAP(CAOBSPTR),
     &              HEAP(CBOBSPTR), HEAP(CACALCPTR), HEAP(CBCALCPTR),
     &              HEAP(EXCAPTR), HEAP(EXCBPTR), HEAP(EXCAERRPTR),
     &              HEAP(EXCBERRPTR), HEAP(RCAPTR), HEAP(RCBPTR),
     &              WHICH)
      RETURN
      END
C===============
      SUBROUTINE ECSHIFT2 (EC, ATOMI, ATOMJ, ATOMK, ATOML, ATOMM,
     &                     CAOBS, CBOBS, CACALC, CBCALC, EXPECTEDCA,
     &                     EXPECTEDCB, EXPECTEDCAERR, EXPECTEDCBERR,
     &                     RCAS, RCBS, WHICH)
C
C Calculates carbon chemical shift energies
C
C chemical shift energies are of the form
C      E = k1*(deltaCa**2+deltaCb**2)
C where
C      k1 = force constant,
C      deltaC = calculated chem shift - observed chem shift
C
C which is a flag that switches between energy & force and
C calculated chem shift (for violations) calcs
C
C by John Kuszewski July 1993
C
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*)
      DOUBLE PRECISION CAOBS(*), CBOBS(*), CACALC(*), CBCALC(*),
     &     RCAS(*), RCBS(*)
      DOUBLE PRECISION EXPECTEDCA(PHISTEPS, PSISTEPS),
     &     EXPECTEDCB(PHISTEPS, PSISTEPS),
     &     EXPECTEDCAERR(PHISTEPS, PSISTEPS),
     &     EXPECTEDCBERR(PHISTEPS, PSISTEPS)
      DOUBLE PRECISION EC
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, CLASS, LOWRPHI,
     &        HIRPHI, RPSI, LOWRPSI, HIRPSI, RPHI
      DOUBLE PRECISION XI, XJ, XK, XL, XM, YI, YJ, YK, YL, YM, ZI,
     &                 ZJ, ZK, ZL,
     &                 ZM, XIJ, XJK, XKL, XLM, YIJ, YJK, YKL, YLM,
     &                 ZIJ, ZJK,
     &                 ZKL, ZLM, AX, AY, AZ, BX, BY, BZ, CX, CY, CZ,
     &                 AX2, AY2, AZ2, BX2, BY2, BZ2, CX2, CY2, CZ2,
     &                 RAR, RBR, RCR, CP, SP, PHI,
     &                 RAR2, RBR2, RCR2, CP2, SP2, PSI,
     &                 EXPCA, EXPCB, DELTACA, DELTACB, E, DF, DF2,
     &                 DCADPHI, DCADPSI, DCBDPHI, DCBDPSI,
     &                 SWITCH, RECSP, RECCP,
     &                 SWITCH2, RECSP2, RECCP2,
     &                 DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ,
     &                 DCPAX2, DCPAY2, DCPAZ2, DCPBX2, DCPBY2, DCPBZ2,
     &                 DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ,
     &                 DSPCX2, DSPCY2, DSPCZ2, DSPBX2, DSPBY2, DSPBZ2,
     &                 DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ,
     &                 DPRIJX2, DPRIJY2, DPRIJZ2, DPRKLX2, DPRKLY2,
     &                 DPRKLZ2, CASIGN, CBSIGN, PHICONV,
     &                 DPRJKX, DPRJKY, DPRJKZ, PSICONV,
     &                 DPRJKX2, DPRJKY2, DPRJKZ2,
     &                 PHISTEPSPERRAD, PSISTEPSPERRAD, HICA, LOWCA,
     &                 DEDPHI, DEDPSI, CAERR, CBERR, HICB, LOWCB
      DOUBLE PRECISION OBSCA, OBSCB, K1
C
C begin
C
      PHICONV = PHISTEPS / (TWO * PI)
      PSICONV = PSISTEPS / (TWO * PI)
C
C following Axel's code in ETOR,
C
C zero out partial energy
C
      EC = ZERO
C
      CLASS = 1
      K1 = SHIFTFORCES(1)
      DO COUNT = 1, NSHIFTS
          IF (SHIFTASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               K1 = SHIFTFORCES(CLASS)
          END IF
C
C get coords of each atom
C
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
C
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
C
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
C
C new random coil implementation
C
          OBSCA = CAOBS(COUNT) - RCAS(ATOMK(COUNT))
          OBSCB = CBOBS(COUNT) - RCBS(ATOMK(COUNT))
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
          XIJ = XI - XJ
          XJK = XJ - XK
          XKL = XK - XL
          XLM = XL - XM
C
          YIJ = YI - YJ
          YJK = YJ - YK
          YKL = YK - YL
          YLM = YL - YM
C
          ZIJ = ZI - ZJ
          ZJK = ZJ - ZK
          ZKL = ZK - ZL
          ZLM = ZL - ZM
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
          AX2 = YJK*ZKL-ZJK*YKL
          AY2 = ZJK*XKL-XJK*ZKL
          AZ2 = XJK*YKL-YJK*XKL
          BX2 = YKL*ZLM-YLM*ZKL
          BY2 = ZKL*XLM-ZLM*XKL
          BZ2 = XKL*YLM-XLM*YKL
          CX2 = YKL*AZ2-ZKL*AY2
          CY2 = ZKL*AX2-XKL*AZ2
          CZ2 = XKL*AY2-YKL*AX2
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
          PHISTEPSPERRAD = PHISTEPS / (TWO * PI)
          PSISTEPSPERRAD = PSISTEPS / (TWO * PI)
          RPHI = NINT((PHI + PI) * PHISTEPSPERRAD)
          RPSI = NINT((PSI + PI) * PSISTEPSPERRAD)
          IF (RPHI.EQ.ZERO) RPHI = PHISTEPS
          IF (RPSI.EQ.ZERO) RPSI = PSISTEPS
C
C look up expectation values
C
          EXPCA = EXPECTEDCA(RPHI,RPSI)
          EXPCB = EXPECTEDCB(RPHI,RPSI)
          CAERR = EXPECTEDCAERR(RPHI,RPSI)
          CBERR = EXPECTEDCBERR(RPHI,RPSI)
          IF (WHICH.EQ.'ANALYZE') THEN
               CACALC(COUNT) = EXPCA
               CBCALC(COUNT) = EXPCB
          END IF
C
C and calculate their deviation from experiment
C
          DELTACA = EXPCA - OBSCA
          DELTACB = EXPCB - OBSCB
C
C don't calculate energy & forces if you're at a spot without
C an expectation value
C
          IF ((EXPCA.GT.NOEXPTEST).OR.(EXPCB.GT.NOEXPTEST)) THEN
               E = 0
               DF = 0
               DF2 = 0
C
C added 5/25/97 ATB
          ELSE IF ((ABS(DELTACA).GT.NOEXPTEST).OR.
     &             (ABS(DELTACB).GT.NOEXPTEST)) THEN
               E = 0
               DF = 0
               DF2 = 0
C
C don't calculate energy & forces if you're in square-well
C mode and within the error
C
          ELSE IF ((ABS(DELTACA).LT.CAERR).AND.(ABS(DELTACB).LT.CBERR)
     &             .AND.(SHIFTPOTENTIAL.EQ.SQUARE)) THEN
               E = 0
               DF = 0
               DF2 = 0
C
C Calculate energy & forces
C
          ELSE
C
C but make sure that you
C don't calculate energy & forces for entries without
C an observed value (eg., gly, his)
C
               IF (ABS(DELTACA).GT.NOEXPTEST) DELTACA = 0
               IF (ABS(DELTACB).GT.NOEXPTEST) DELTACB = 0
C
C subtract off the errors if you're in square well mode
C
               IF (SHIFTPOTENTIAL.EQ.SQUARE) THEN
                    E = K1*(((ABS(DELTACA)-CAERR)**2) +
     &                      ((ABS(DELTACB)-CBERR)**2))
               ELSE
                    E = K1*((DELTACA**2)+(DELTACB**2))
               END IF
C
C get indexes of next higher and lower points along psi,
C get the values, calculate the local slope,
C making sure that you don't use datapoints that don't
C have an expectation value (slopes at edges will therefore
C be somewhat less accurate, but so what...)
C and calculate the derivative wrt psi
C
               LOWRPHI = RPHI - 1
               IF (LOWRPHI.LE.0) LOWRPHI = LOWRPHI + PHISTEPS
               HIRPHI = RPHI + 1
               IF (HIRPHI.GT.PHISTEPS) HIRPHI = HIRPHI - PHISTEPS
               HICA = EXPECTEDCA(HIRPHI, RPSI)
               LOWCA = EXPECTEDCA(LOWRPHI, RPSI)
C
C Changed here, in dCb/dphi, and in dC*/dpsi below
C 7/14/95 JJK
C
               IF ((HICA.GT.NOEXPTEST).AND.(LOWCA.GT.NOEXPTEST)) THEN
                    DCADPHI = ZERO
               ELSE IF (HICA.GT.NOEXPTEST) THEN
                    DCADPHI = ((EXPCA - LOWCA) / ONE) * PHICONV
               ELSE IF (LOWCA.GT.NOEXPTEST) THEN
                    DCADPHI = ((HICA - EXPCA) / ONE) * PHICONV
               ELSE
                    DCADPHI = ((HICA - LOWCA) / TWO) * PHICONV
               END IF
               HICB = EXPECTEDCB(HIRPHI, RPSI)
               LOWCB = EXPECTEDCB(LOWRPHI, RPSI)
               IF ((HICB.GT.NOEXPTEST).AND.(LOWCB.GT.NOEXPTEST)) THEN
                    DCBDPHI = ZERO
               ELSE IF (HICB.GT.NOEXPTEST) THEN
                    DCBDPHI = ((EXPCB - LOWCB) / ONE) * PHICONV
               ELSE IF (LOWCB.GT.NOEXPTEST) THEN
                    DCBDPHI = ((HICB - EXPCB) / ONE) * PHICONV
               ELSE
                    DCBDPHI = ((HICB - LOWCB) / TWO) * PHICONV
               END IF
C
C turn slopes into derivatives
C
               IF (SHIFTPOTENTIAL.EQ.SQUARE) THEN
C
C avoid division-by-zero in flat areas of the
C expectation grid.
C
                    IF (DELTACA.NE.ZERO) THEN
                         CASIGN = DELTACA / ABS(DELTACA)
                    ELSE
                         CASIGN = ONE
                    END IF
                    IF (DELTACB.NE.ZERO) THEN
                         CBSIGN = DELTACB / ABS(DELTACB)
                    ELSE
                         CBSIGN = ONE
                    END IF
                    DEDPHI = TWO*K1*
     &                       (((ABS(DELTACA)-CAERR) *
     &                              CASIGN*DCADPHI) +
     &                        ((ABS(DELTACB)-CBERR) *
     &                              CBSIGN*DCBDPHI))
               ELSE
                    DEDPHI = TWO*K1*(DELTACA * DCADPHI +
     &                               DELTACB * DCBDPHI)
               END IF
C
C now do the same for psi
C
               LOWRPSI = RPSI - 1
               IF (LOWRPSI.LE.0) LOWRPSI = LOWRPSI + PSISTEPS
               HIRPSI = RPSI + 1
               IF (HIRPSI.GT.PSISTEPS) HIRPSI = HIRPSI - PSISTEPS
               HICA = EXPECTEDCA(RPHI, HIRPSI)
               LOWCA = EXPECTEDCA(RPHI, LOWRPSI)
               IF ((HICA.GT.NOEXPTEST).AND.(LOWCA.GT.NOEXPTEST)) THEN
                    DCADPSI = ZERO
               ELSE IF (HICA.GT.NOEXPTEST) THEN
                    DCADPSI = ((EXPCA - LOWCA) / ONE) * PSICONV
               ELSE IF (LOWCA.GT.NOEXPTEST) THEN
                    DCADPSI = ((HICA - EXPCA) / ONE) * PSICONV
               ELSE
                    DCADPSI = ((HICA - LOWCA) / TWO) * PSICONV
               END IF
               HICB = EXPECTEDCB(RPHI, HIRPSI)
               LOWCB = EXPECTEDCB(RPHI, LOWRPSI)
               IF ((HICB.GT.NOEXPTEST).AND.(LOWCB.GT.NOEXPTEST)) THEN
                    DCBDPSI = ZERO
               ELSE IF (HICB.GT.NOEXPTEST) THEN
                    DCBDPSI = ((EXPCB - LOWCB) / ONE) * PSICONV
               ELSE IF (LOWCB.GT.NOEXPTEST) THEN
                    DCBDPSI = ((HICB - EXPCB) / ONE) * PSICONV
               ELSE
                    DCBDPSI = ((HICB - LOWCB) / TWO) * PSICONV
               END IF
               IF (SHIFTPOTENTIAL.EQ.SQUARE) THEN
                    IF (DELTACA.NE.ZERO) THEN
                         CASIGN = DELTACA / ABS(DELTACA)
                    ELSE
                         CASIGN = ONE
                    END IF
                    IF (DELTACB.NE.ZERO) THEN
                         CBSIGN = DELTACB / ABS(DELTACB)
                    ELSE
                         CBSIGN = ONE
                    END IF
                    DEDPSI = TWO*K1*
     &                       (((ABS(DELTACA)-CAERR) *
     &                              CASIGN*DCADPSI) +
     &                        ((ABS(DELTACB)-CBERR) *
     &                              CBSIGN*DCBDPSI))
               ELSE
                    DEDPSI = TWO*K1*(DELTACA * DCADPSI +
     &                               DELTACB * DCBDPSI)
               END IF
C
C now put the values into the vbls to be used later
C
               DF = DEDPHI
               DF2 = DEDPSI
          END IF
C
C accumulate energy
C
          EC = EC + E
C
C compute heavyside function
C
          SWITCH=-MIN(1,INT(ABS(SP)-CARBEPS+ONE))
          SWITCH2=-MIN(1,INT(ABS(SP2)-CARBEPS+ONE))
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
     &     RECSP2*(DCPAY2*ZJK-DCPAZ2*YJK+DCPBZ2*YLM-DCPBY2*ZLM)+
     &     RECCP2*(-(YKL*YJK+ZKL*ZJK)*DSPCX2
     &        +(TWO*XKL*YJK-XJK*YKL)*DSPCY2
     &        +(TWO*XKL*ZJK-XJK*ZKL)*DSPCZ2
     &        +DSPBZ2*YLM-DSPBY2*ZLM)
          DPRJKY2=
     &     RECSP2*(DCPAZ2*XJK-DCPAX2*ZJK+DCPBX2*ZLM-DCPBZ2*XLM)+
     &     RECCP2*(-(ZKL*ZJK+XKL*XJK)*DSPCY2
     &        +(TWO*YKL*ZJK-YJK*ZKL)*DSPCZ2
     &        +(TWO*YKL*XJK-YJK*XKL)*DSPCX2
     &        +DSPBX2*ZLM-DSPBZ2*XLM)
          DPRJKZ2=
     &     RECSP2*(DCPAX2*YJK-DCPAY2*XJK+DCPBY2*XLM-DCPBX2*YLM)+
     &     RECCP2*(-(XKL*XJK+YKL*YJK)*DSPCZ2
     &        +(TWO*ZKL*XJK-ZJK*XKL)*DSPCX2
     &        +(TWO*ZKL*YJK-ZJK*YKL)*DSPCY2
     &        +DSPBY2*XLM-DSPBX2*YLM)
C
C now update forces if in energy & force mode
C
          IF (WHICH.NE.'ANALYZE') THEN
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
               DX(ATOMJ(COUNT))=DX(ATOMJ(COUNT))+DPRIJX2
               DY(ATOMJ(COUNT))=DY(ATOMJ(COUNT))+DPRIJY2
               DZ(ATOMJ(COUNT))=DZ(ATOMJ(COUNT))+DPRIJZ2
               DX(ATOMK(COUNT))=DX(ATOMK(COUNT))+DPRJKX2-DPRIJX2
               DY(ATOMK(COUNT))=DY(ATOMK(COUNT))+DPRJKY2-DPRIJY2
               DZ(ATOMK(COUNT))=DZ(ATOMK(COUNT))+DPRJKZ2-DPRIJZ2
               DX(ATOML(COUNT))=DX(ATOML(COUNT))+DPRKLX2-DPRJKX2
               DY(ATOML(COUNT))=DY(ATOML(COUNT))+DPRKLY2-DPRJKY2
               DZ(ATOML(COUNT))=DZ(ATOML(COUNT))+DPRKLZ2-DPRJKZ2
               DX(ATOMM(COUNT))=DX(ATOMM(COUNT))        -DPRKLX2
               DY(ATOMM(COUNT))=DY(ATOMM(COUNT))        -DPRKLY2
               DZ(ATOMM(COUNT))=DZ(ATOMM(COUNT))        -DPRKLZ2
          END IF
      END DO
      RETURN
      END
C================
      SUBROUTINE READCSHIFT
C
C reads in chemical shift information
C
C by John Kuszewski July 1993
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXSHIFTS
      DOUBLE PRECISION K1, CUTOFF
      CHARACTER*4 THENAME
      DOUBLE COMPLEX DUMMY
      DOUBLE PRECISION DUMMY2
C begin
C
C allocate memory if we haven't already done so
C
      IF (.NOT.CARBRCALLOCATED) THEN
      CALL CSHIFTDEFAULTS
      CALL ALLOCSHIFTS(0, MAXSHIFTS)
      DUMMY2 = NOEXPECTATION
      CALL DECLAR('NOEXPECTATION', 'DP', ' ', DUMMY, DUMMY2)
      CALL ALLOCEXPECTED
      CARBNATOM=NATOM
      RCAPTR = ALLHP(IREAL8(CARBNATOM))
      RCBPTR = ALLHP(IREAL8(CARBNATOM))
      CARBRCALLOCATED = .TRUE.
      END IF
C
C this is used by READSHIFT2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
C now read input
C
      CALL PUSEND('CSHIFTS>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('CSHIFTS>')
           CALL MISCOM('CSHIFTS>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-carbon')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURSHIFTCLASS
                CALL NEXTA4('class name =', THENAME)
                SHIFTMODE = NEW
                DO COUNT = 1, NSHIFTCLASSES
                     IF (SHIFTCLASSNAMES(COUNT).EQ.THENAME) THEN
                          SHIFTMODE = UPDATE
                          CURSHIFTCLASS = COUNT
                     END IF
                END DO
                IF (SHIFTMODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXSHIFTCLASSES) THEN
                         CALL DSPERR('CHEMSHIFTS','Too many classes.')
                         CALL DSPERR('CHEMSHIFTS',
     &                      'Increase MAXSHIFTCLASSES and recompile.')
                         CALL WRNDIE(-5, 'READCHEM',
     &                               'Too many chem shift classes.')
                     END IF
                     NSHIFTCLASSES = NSHIFTCLASSES + 1
                     CURSHIFTCLASS = NSHIFTCLASSES
                     SHIFTCLASSNAMES(CURSHIFTCLASS) = THENAME
C
C If this isn't the first class, close off the old class
C
                IF (NSHIFTCLASSES.GT.1) THEN
                     SHIFTASSNDX(OLDCLASS) = NSHIFTS
                END IF
            END IF
C
C set force constant for current class
C
           ELSE IF (WD(1:4).EQ.'FORC') THEN
                CALL NEXTF('force constant =', K1)
C
C start a default class if there isn't one defined
C
                IF (CURSHIFTCLASS.EQ.0) THEN
                     NSHIFTCLASSES = 1
                     CURSHIFTCLASS = 1
                END IF
                WRITE(PUNIT, '(3(A), F8.3)')
     &            'Setting frce const for class ',
     &             SHIFTCLASSNAMES(CURSHIFTCLASS), ' to ', K1
                SHIFTFORCES(CURSHIFTCLASS) = K1
C
C reset shifts database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
                CALL CARBHP
                CALL CSHIFTINIT
                CALL CSHIFTDEFAULTS
                CALL ALLOCSHIFTS(0, MAXSHIFTS)
                DUMMY2 = NOEXPECTATION
                CALL DECLAR('NOEXPECTATION', 'DP', ' ', DUMMY, DUMMY2)
                CALL ALLOCEXPECTED
                CARBNATOM=NATOM
                RCAPTR = ALLHP(IREAL8(CARBNATOM))
                RCBPTR = ALLHP(IREAL8(CARBNATOM))
                CARBRCALLOCATED = .TRUE.
C
C set potential type
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                    WRITE(PUNIT, '(A)') 'using square well potential.'
                    SHIFTPOTENTIAL = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                    WRITE(PUNIT, '(A)') 'using harmonic potential.'
                    SHIFTPOTENTIAL = HARMONIC
                ELSE
                    CALL DSPERR('CHEMSHIFTS',
     &                        'unknown potential. Using square well.')
                    SHIFTPOTENTIAL = SQUARE
                END IF
C
C change number of assignment slots
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXSHIFTS = MAXSHIFTS
                CALL NEXTI('number of slots =', MAXSHIFTS)
                CALL ALLOCSHIFTS(OLDMAXSHIFTS, MAXSHIFTS)
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more chemical shifts
C than you have slots for
C
                IF (NSHIFTS.EQ.MAXSHIFTS) THEN
                     CALL DSPERR('CSHIFTS','Too many assignments.')
                     CALL DSPERR('CSHIFTS',
     &                    'Increase NREStraints and run again.')
                     CALL WRNDIE(-1,'CSHIFTS>',
     & 'exceeded allocation for carbon chemical shifts-restraints')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURSHIFTCLASS.EQ.0) THEN
                     NSHIFTCLASSES = 1
                     CURSHIFTCLASS = 1
                END IF
                CALL READCSHIFTASSIGN(HEAP(SHIFTIPTR), HEAP(SHIFTJPTR),
     &                                HEAP(SHIFTKPTR), HEAP(SHIFTLPTR),
     &                                HEAP(SHIFTMPTR), HEAP(SPTR),
     &                                HEAP(CAOBSPTR), HEAP(CBOBSPTR))
C
C set number of steps in phi dimension of expectation value arrays
C
           ELSE IF (WD(1:4).EQ.'PHIS') THEN
                CALL NEXTI('phi step =', PHISTEPS)
                CALL ALLOCEXPECTED
C
C set number of steps in psi dimension of expectation value arrays
C
           ELSE IF (WD(1:4).EQ.'PSIS') THEN
                CALL NEXTI('psi step =', PSISTEPS)
                CALL ALLOCEXPECTED
C
C read in expectation values
C
           ELSE IF (WD(1:4).EQ.'EXPE') THEN
                CALL READEXPECTEDSHIFT(HEAP(EXCAPTR),HEAP(EXCBPTR),
     &                             HEAP(EXCAERRPTR),HEAP(EXCBERRPTR))
C
C read in the random coil values
C
           ELSE IF (WD(1:4).EQ.'RCOI') THEN
                CALL READCARBRCSHIFT(HEAP(RCAPTR), HEAP(RCBPTR),
     &                               HEAP(SPTR))
C
C zero out the expectation value arrays
C
           ELSE IF (WD(1:4).EQ.'ZERO') THEN
                CALL DEFEXPECTEDSHIFT(HEAP(EXCAPTR),HEAP(EXCBPTR),
     &                             HEAP(EXCAERRPTR),HEAP(EXCBERRPTR))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                     CALL DSPERR('CSHIFTS',
     &                           'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                          CALL DSPERR('CSHIFTS',
     &                                'cutoff must be positive.')
                          CUTOFF = ABS(CUTOFF)
                     END IF
                     CALL PRINTCSHIFTS(CUTOFF, HEAP(CACALCPTR),
     &                                 HEAP(CBCALCPTR),
     &                                 HEAP(CAOBSPTR), HEAP(CBOBSPTR),
     &                                 HEAP(SHIFTKPTR),
     &                                 HEAP(RCAPTR), HEAP(RCBPTR))
                END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('CSHIFTS>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCEXPECTED
C
C Allocates space for EXPECTEDCA and EXPECTEDCB and their error arrays
C based on values of PHISTEP and PSISTEP.
C Fills each array with zeros.
C
C by John Kuszewski May 1994
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C local variables
      INTEGER NPOINTS
C begin
      NPOINTS = PHISTEPS*PSISTEPS
      IF (EXCAPTR.NE.0) THEN
            CALL FREHP(EXCAPTR, IREAL8(NPOINTS))
            CALL FREHP(EXCBPTR, IREAL8(NPOINTS))
            CALL FREHP(EXCAERRPTR, IREAL8(NPOINTS))
            CALL FREHP(EXCBERRPTR, IREAL8(NPOINTS))
      END IF
      EXCAPTR = 0
      EXCBPTR = 0
      EXCAERRPTR = 0
      EXCBERRPTR = 0
      IF (NPOINTS.NE.0) THEN
        EXCAPTR = ALLHP(IREAL8(NPOINTS))
        EXCBPTR = ALLHP(IREAL8(NPOINTS))
        EXCAERRPTR = ALLHP(IREAL8(NPOINTS))
        EXCBERRPTR = ALLHP(IREAL8(NPOINTS))
C
C now zero them out
C
        CALL DEFEXPECTEDSHIFT(HEAP(EXCAPTR),HEAP(EXCBPTR),
     &                        HEAP(EXCAERRPTR),HEAP(EXCBERRPTR))
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE ALLOCSHIFTS (OLDSIZE, NEWSIZE)
C
C resets chem shift arrays to hold SIZE entries
C
C by John Kuszewski Nov 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(SHIFTIPTR, INTEG4(OLDSIZE))
           CALL FREHP(SHIFTJPTR, INTEG4(OLDSIZE))
           CALL FREHP(SHIFTKPTR, INTEG4(OLDSIZE))
           CALL FREHP(SHIFTLPTR, INTEG4(OLDSIZE))
           CALL FREHP(SHIFTMPTR, INTEG4(OLDSIZE))
           CALL FREHP(CAOBSPTR, IREAL8(OLDSIZE))
           CALL FREHP(CBOBSPTR, IREAL8(OLDSIZE))
           CALL FREHP(CACALCPTR, IREAL8(OLDSIZE))
           CALL FREHP(CBCALCPTR, IREAL8(OLDSIZE))
      END IF
C
      SHIFTIPTR = 0
      SHIFTJPTR = 0
      SHIFTKPTR = 0
      SHIFTLPTR = 0
      SHIFTMPTR = 0
      CAOBSPTR = 0
      CBOBSPTR = 0
      CACALCPTR = 0
      CBCALCPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        SHIFTIPTR = ALLHP(INTEG4(NEWSIZE))
        SHIFTJPTR = ALLHP(INTEG4(NEWSIZE))
        SHIFTKPTR = ALLHP(INTEG4(NEWSIZE))
        SHIFTLPTR = ALLHP(INTEG4(NEWSIZE))
        SHIFTMPTR = ALLHP(INTEG4(NEWSIZE))
        CAOBSPTR = ALLHP(IREAL8(NEWSIZE))
        CBOBSPTR = ALLHP(IREAL8(NEWSIZE))
        CACALCPTR = ALLHP(IREAL8(NEWSIZE))
        CBCALCPTR = ALLHP(IREAL8(NEWSIZE))
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE CSHIFTINIT
C
C initializes chem shift stuff
C
C by John Kuszewski Nov 1993
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
C begin
      CARBRCALLOCATED=.FALSE.
      RETURN
      END
C==============
      SUBROUTINE CARBHP
C
C deallocates carbon shifts stuff
C
C by Alexandre Bonvin 4/8/96
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C local
      INTEGER NPOINTS
C begin
      IF (CARBRCALLOCATED) THEN
C
      CALL FREHP(RCAPTR,IREAL8(CARBNATOM))
      CALL FREHP(RCBPTR,IREAL8(CARBNATOM))
      CARBRCALLOCATED = .FALSE.
      NPOINTS = PHISTEPS*PSISTEPS
      IF (EXCAPTR.NE.0) THEN
      CALL FREHP(EXCAPTR, IREAL8(NPOINTS))
      CALL FREHP(EXCBPTR, IREAL8(NPOINTS))
      CALL FREHP(EXCAERRPTR, IREAL8(NPOINTS))
      CALL FREHP(EXCBERRPTR, IREAL8(NPOINTS))
      END IF
      CALL ALLOCSHIFTS(MAXSHIFTS, 0)
C
      END IF
      RETURN
      END
C==============
      SUBROUTINE CSHIFTDEFAULTS
C
C sets up defaults
C
C by John Kuszewski Nov 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
C local variables
      INTEGER COUNT
C begin
      SHIFTMODE = NEW
      MAXSHIFTS = 500
      NSHIFTS = 0
      NSHIFTCLASSES = 0
      CURSHIFTCLASS = 0
      PHISTEPS = 180
      PSISTEPS = 180
      SHIFTPOTENTIAL = HARMONIC
      CARBRCALLOCATED = .FALSE.
      DO COUNT = 1, MAXSHIFTCLASSES
           SHIFTCLASSNAMES(COUNT) = 'DEFAULT'
           SHIFTASSNDX(COUNT) = 0
           SHIFTFORCES(COUNT) = 50
      END DO
      EXCAPTR = 0
      EXCBPTR = 0
      EXCAERRPTR = 0
      EXCBERRPTR = 0
      RETURN
      END
C==============
      SUBROUTINE DEFEXPECTEDSHIFT (EXPECTEDCA, EXPECTEDCB,
     &                             EXPECTEDCAERR, EXPECTEDCBERR)
C
C fills arrays that hold expectation values with a big number
C Can't use zero because we need to distinguish between
C an expectation of zero and no expectation.
C
C by John Kuszewski June 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
C i/o
      DOUBLE PRECISION EXPECTEDCA(PHISTEPS, PSISTEPS),
     &     EXPECTEDCB(PHISTEPS, PSISTEPS),
     &     EXPECTEDCAERR(PHISTEPS, PSISTEPS),
     &     EXPECTEDCBERR(PHISTEPS, PSISTEPS)
C local vbls
      INTEGER COUNT1, COUNT2
C begin
      DO COUNT1 = 1, PHISTEPS
      DO COUNT2 = 1, PSISTEPS
           EXPECTEDCA(COUNT1,COUNT2) = NOEXPECTATION
           EXPECTEDCB(COUNT1,COUNT2) = NOEXPECTATION
           EXPECTEDCAERR(COUNT1,COUNT2) = NOEXPECTATION
           EXPECTEDCBERR(COUNT1,COUNT2) = NOEXPECTATION
      END DO
      END DO
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDSHIFT (EXPECTEDCA, EXPECTEDCB,
     &                              EXPECTEDCAERR, EXPECTEDCBERR)
C
C actually reads in the expected values for Ca and Cb
C
C by John Kuszewski May 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
C i/o
      DOUBLE PRECISION EXPECTEDCA(PHISTEPS, PSISTEPS),
     &     EXPECTEDCB(PHISTEPS, PSISTEPS),
     &     EXPECTEDCAERR(PHISTEPS, PSISTEPS),
     &     EXPECTEDCBERR(PHISTEPS, PSISTEPS)
C local vbls
      INTEGER PHIPOS, PSIPOS
      DOUBLE PRECISION CAVAL, CBVAL, CAERR, CBERR
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('CHEMSHIFTS',
     &               'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTI('PSIposition =', PSIPOS)
      IF ((PSIPOS.GT.PSISTEPS).OR.(PSIPOS.LT.1)) THEN
           CALL DSPERR('CHEMSHIFTS',
     &               'psi not in range 1..PSISTEPS. Assuming posn 1.')
           PSIPOS = 1
      END IF
      CALL NEXTF('expected CA value =', CAVAL)
      CALL NEXTF('expected CA error =', CAERR)
      CALL NEXTF('expected CB value =', CBVAL)
      CALL NEXTF('expected CB error =', CBERR)
      EXPECTEDCA(PHIPOS,PSIPOS) = CAVAL
      EXPECTEDCAERR(PHIPOS,PSIPOS) = CAERR
      EXPECTEDCB(PHIPOS,PSIPOS) = CBVAL
      EXPECTEDCBERR(PHIPOS,PSIPOS) = CBERR
      RETURN
      END
C==============
      SUBROUTINE READCARBRCSHIFT (RCAS, RCBS, SEL)
C
C Reads the random coil entries
C
C by John Kuszewski Apr 1995
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cshifts.inc'
C i/o
      DOUBLE PRECISION RCAS(*), RCBS(*)
      INTEGER SEL(*)
C local vbls
      INTEGER NSEL, COUNT
      DOUBLE PRECISION CA, CB
C begin
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      CALL NEXTF('13CA rcoil shift =', CA)
      CALL NEXTF('13CB rcoil shift =', CB)
      DO COUNT = 1, NSEL
           RCAS(SEL(COUNT)) = CA
           RCBS(SEL(COUNT)) = CB
      END DO
      RETURN
      END
C==============
      SUBROUTINE READCSHIFTASSIGN (ATOMI, ATOMJ, ATOMK, ATOML, ATOMM,
     &                            SEL, CAOBS, CBOBS)
C
C reads actual chem shift assignments into arrays
C
C by John Kuszewski Nov 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), ATOMM(*),
     &        SEL(*)
      DOUBLE PRECISION CAOBS(*), CBOBS(*)
C local variables
      INTEGER NSEL, INSERTPOS, COUNT, CURSTOP, OTHERSTOP
      DOUBLE PRECISION CA, CB
      INTEGER I, J, K, L, M
      LOGICAL USEME
C begin
C
      USEME = .TRUE.
C
C get the atoms and observed values
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('CARBSHIFTS',
     &     'more than 1 atom in selection for atom Ci-1. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      I = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('CARBSHIFTS',
     &     'more than 1 atom in selection for atom Ni. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      J = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('CARBSHIFTS',
     &     'more than 1 atom in selection for atom Cai. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      K = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('CARBSHIFTS',
     &     'more than 1 atom in selection for atom Ci. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      L = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('CARBSHIFTS',
     &     'more than 1 atom in selection for atom Ni+1. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      IF (NSEL.LT.1) USEME = .FALSE.
      M = SEL(1)
C
      CALL NEXTF('observed Ca shift =', CA)
      CALL NEXTF('observed Cb shift =', CB)
C
C see if we should ignore this one or not
C
      IF (USEME) THEN
C
C if we're in update mode, make a space for the new line
C
           IF (SHIFTMODE.EQ.UPDATE) THEN
                DO COUNT = NSHIFTS+1, SHIFTASSNDX(CURSHIFTCLASS)+1, -1
                     ATOMI(COUNT) = ATOMI(COUNT-1)
                     ATOMJ(COUNT) = ATOMJ(COUNT-1)
                     ATOMK(COUNT) = ATOMK(COUNT-1)
                     ATOML(COUNT) = ATOML(COUNT-1)
                     ATOMM(COUNT) = ATOMM(COUNT-1)
                     CAOBS(COUNT) = CAOBS(COUNT-1)
                     CBOBS(COUNT) = CBOBS(COUNT-1)
                END DO
                CURSTOP = SHIFTASSNDX(CURSHIFTCLASS)
                DO COUNT = 1, NSHIFTCLASSES
                     OTHERSTOP = SHIFTASSNDX(COUNT)
                     IF (OTHERSTOP.GT.CURSTOP) THEN
                          SHIFTASSNDX(COUNT) = OTHERSTOP + 1
                     END IF
                END DO
                SHIFTASSNDX(CURSHIFTCLASS) = CURSTOP + 1
                INSERTPOS = CURSTOP
                NSHIFTS = NSHIFTS + 1
           ELSE
                NSHIFTS = NSHIFTS + 1
                INSERTPOS = NSHIFTS
                SHIFTASSNDX(CURSHIFTCLASS) = INSERTPOS
           END IF
C
C now write the data to the arrays
C
           ATOMI(INSERTPOS) = I
           ATOMJ(INSERTPOS) = J
           ATOMK(INSERTPOS) = K
           ATOML(INSERTPOS) = L
           ATOMM(INSERTPOS) = M
           CAOBS(INSERTPOS) = CA
           CBOBS(INSERTPOS) = CB
      ELSE
           CALL DSPERR('CARBSHIFTS',
     &     'no atoms in one or more selections. Ignoring this entry.')
      END IF
      RETURN
      END
C=================
      SUBROUTINE PRINTCSHIFTS (CUTOFF, CACALC, CBCALC, CAOBS,
     &                        CBOBS, ATOMK, RCAS, RCBS)
C
C prints chemical shifts with difference between expected and
C actual Ca or Cb shifts greater than cutoff.
C calculates RMS deviation and puts it into $RMSCA AND $RMSCB
C
C by John Kuszewski Nov 1993
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'cshifts.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, CACALC(*), CBCALC(*), CAOBS(*),
     &                 CBOBS(*), RCAS(*), RCBS(*)
      INTEGER ATOMK(*)
C local variables
      DOUBLE PRECISION CALCCA, CALCCB, OBSCA, OBSCB, NUMSHIFTSINCARMS,
     &                 DELTACA, DELTACB, SENERGY, RMSCA, RMSCB, VIOLS,
     &                 NUMSHIFTSINCBRMS, RCCA, RCCB, ABSDELTACA,
     &                 ABSDELTACB, DA, DB, SHIFTCOUNT
      INTEGER COUNT, CLASS, K
      DOUBLE COMPLEX DUMMY2
C parameters
      DOUBLE PRECISION UNKNOWN
      PARAMETER (UNKNOWN=-99.99D0)
C begin
      RMSCA = ZERO
      RMSCB = ZERO
      VIOLS = ZERO
      NUMSHIFTSINCARMS = ZERO
      NUMSHIFTSINCBRMS = ZERO
      SHIFTCOUNT = ZERO
C
C make sure that the calc-shift array is up to date
C
      CALL ECSHIFT(SENERGY, 'ANALYZE')
      WRITE (PUNIT, '(A)') 'The following Ca or Cb delta secondary'
      WRITE (PUNIT, '(2A)') 'shifts are greater than the cutoff ',
     &'(-99.99 = no expectation):'
      WRITE (PUNIT, '(A)')
     &'(Ca) (obsCa) (rcoilCa) (obsdeltaCa) (calcdeltaCa) (deltadeltaCa)'
      WRITE (PUNIT, '(A)')
     &'     (obsCb) (rcoilCb) (obsdeltaCb) (calcdeltaCb) (deltadeltaCb)'
C
C write out first class heading
C
      CLASS = 1
      WRITE (PUNIT, '(A, A)') 'class ', SHIFTCLASSNAMES(1)
C
C for every shift entry,
C
      DO COUNT = 1, NSHIFTS
C
C is this the start of a new class?
C
          IF (SHIFTASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               WRITE (PUNIT, '(A, A)') 'class ',
     &                                 SHIFTCLASSNAMES(CLASS)
          END IF
C
C calc delta shifts
          K = ATOMK(COUNT)
          CALCCA = CACALC(COUNT)
          CALCCB = CBCALC(COUNT)
          RCCA = RCAS(K)
          RCCB = RCBS(K)
          OBSCA = CAOBS(COUNT)
          OBSCB = CBOBS(COUNT)
          DELTACA = CALCCA-(OBSCA-RCCA)
          DELTACB = CALCCB-(OBSCB-RCCB)
          ABSDELTACA = ABS(DELTACA)
          ABSDELTACB = ABS(DELTACB)
          DA=(OBSCA-RCCA)
          DB=(OBSCB-RCCB)
C
C update RMS values if we can make a guess about this shift
C
          IF (ABSDELTACA.LT.NOEXPTEST) THEN
               RMSCA = RMSCA + DELTACA**2
               NUMSHIFTSINCARMS = NUMSHIFTSINCARMS + ONE
          END IF
          IF (ABSDELTACB.LT.NOEXPTEST) THEN
               RMSCB = RMSCB + DELTACB**2
               NUMSHIFTSINCBRMS = NUMSHIFTSINCBRMS + ONE
          END IF
C
          IF (CALCCA.GT.NOEXPTEST) THEN
          CALCCA=UNKNOWN
          DELTACA=UNKNOWN
          END IF
          IF (CALCCB.GT.NOEXPTEST) THEN
          CALCCB=UNKNOWN
          DELTACB=UNKNOWN
          END IF
          IF (RCCA.GT.NOEXPTEST) THEN
          RCCA=UNKNOWN
          DELTACA=UNKNOWN
          DA=UNKNOWN
          END IF
          IF (RCCB.GT.NOEXPTEST) THEN
          RCCB=UNKNOWN
          DELTACB=UNKNOWN
          DB=UNKNOWN
          END IF
          IF (OBSCA.GT.NOEXPTEST) THEN
          OBSCA=UNKNOWN
          DELTACA=UNKNOWN
          DA=UNKNOWN
          END IF
          IF (OBSCB.GT.NOEXPTEST) THEN
          OBSCB=UNKNOWN
          DELTACB=UNKNOWN
          DB=UNKNOWN
          END IF
C
C count the number of shift restraints
C Greg Warren 03/06/98
          IF (OBSCA.GT.ZERO) THEN
          SHIFTCOUNT = SHIFTCOUNT + ONE
          END IF
          IF (OBSCB.GT.ZERO) THEN
          SHIFTCOUNT = SHIFTCOUNT + ONE
          END IF
C print out delta shifts greater than cutoff
C and update number of violations
C
          IF ((ABSDELTACA.GT.CUTOFF).OR.(ABSDELTACB.GT.CUTOFF)) THEN
               WRITE(PUNIT,'(4(A), 3F6.2,2F6.2,3F6.2,2F6.2)')
     &               SEGID(K),RESID(K),RES(K),TYPE(K),
     &     OBSCA, RCCA, DA, CALCCA, DELTACA,
     &     OBSCB, RCCB, DB, CALCCB, DELTACB
               VIOLS = VIOLS + ONE
          END IF
      END DO
      IF (NUMSHIFTSINCARMS.NE.ZERO) THEN
           RMSCA = SQRT(RMSCA / NUMSHIFTSINCARMS)
      ELSE
           RMSCA = ZERO
      END IF
      IF (NUMSHIFTSINCBRMS.NE.ZERO) THEN
           RMSCB = SQRT(RMSCB / NUMSHIFTSINCBRMS)
      ELSE
           RMSCB = ZERO
      END IF
C
      WRITE (PUNIT, '(A, F15.3)') 'RMS Ca error = ', RMSCA
      WRITE (PUNIT, '(A, F15.3)') 'RMS Cb error = ', RMSCB
      WRITE (PUNIT, '(A, F5.2, A, F6.0, A, F6.0, A)')
     & '#violations>',CUTOFF,' = ', VIOLs,
     & ' of ', (NUMSHIFTSINCARMS+NUMSHIFTSINCBRMS) ,
     & ' Ca+Cb chemical shifts'
      CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, SHIFTCOUNT)
      CALL DECLAR('RMSCA', 'DP', ' ', DUMMY2, RMSCA)
      CALL DECLAR('RMSCB', 'DP', ' ', DUMMY2, RMSCB)
      CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      RETURN
      END
C
      SUBROUTINE SCRCSFT
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'cshifts.inc'
C
C reset carbon shift database
      IF (CARBRCALLOCATED) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: Carbon chemical shift database erased.'
      CALL CARBHP
      CALL CSHIFTINIT
      END IF
C
      RETURN
      END
C
