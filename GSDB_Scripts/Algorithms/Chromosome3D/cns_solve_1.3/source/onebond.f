C===============
      SUBROUTINE EONEBOND (EJ, WHICH)
C
C Calls EONEBOND2, which does the actual energy calculation
C
C by John Kuszewski May 1995
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'onebond.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION EJ
      CHARACTER*7 WHICH
C
      CALL EONEBOND2(EJ, HEAP(ONEIPTR), HEAP(ONEJPTR),
     &     HEAP(ONEKPTR), HEAP(ONELPTR), HEAP(ONEMPTR),
     &     HEAP(ONENPTR), HEAP(ONEPPTR), HEAP(ONEQPTR),
     &     HEAP(ONEJOBSPTR), HEAP(ONEJERRPTR), HEAP(CALCONEPTR),
     &     WHICH)
      RETURN
      END
C===============
      SUBROUTINE EONEBOND2 (EJ, ATOMI, ATOMJ, ATOMK, ATOML,
     &     ATOMM, ATOMN, ATOMP, ATOMQ, JOBS, JERR, JCALC, WHICH)
C
C Calculates one bond coupling constant energies
C
C J energies are of the form
C      E = k1*deltaJ**2
C where
C      k1 = force constant,
C      deltaJ = calculated J - observed J
C
C which is a flag that switches between energy & force and
C calculated J (for violations) calcs
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
      INCLUDE 'onebond.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     ATOMM(*), ATOMN(*), ATOMP(*), ATOMQ(*)
      DOUBLE PRECISION JOBS(*), JERR(*), JCALC(*)
      DOUBLE PRECISION EJ
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, CLASS
      DOUBLE PRECISION XI, XJ, XK, XL, YI, YJ, YK, YL, ZI, ZJ, ZK, ZL,
     &     XIJ, XJK, XKL, YIJ, YJK, YKL, ZIJ, ZJK, ZKL,
     &     AX, AY, AZ, BX, BY, BZ, CX, CY, CZ,
     &     RAR, RBR, RCR, CP, SP, PHI,
     &     CALCJ, DELTAJ, E, DF,
     &     SWITCH, RECSP, RECCP,
     &     DCPAX, DCPAY, DCPAZ, DCPBX, DCPBY, DCPBZ,
     &     DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY, DSPBZ,
     &     DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ,
     &     DPRJKX, DPRJKY, DPRJKZ, A, B, C, D, P1, P2, P3
      DOUBLE PRECISION XM, XN, XP, XQ, YM, YN, YP, YQ, ZM, ZN, ZP, ZQ,
     &     XIJ2, XJK2, XKL2, YIJ2, YJK2,
     &     YKL2, ZIJ2, ZJK2, ZKL2,
     &     AX2, AY2, AZ2, BX2, BY2, BZ2, CX2, CY2, CZ2,
     &     RAR2, RBR2, RCR2, CP2, SP2, PSI,
     &     DF2,
     &     SWITCH2, RECSP2, RECCP2,
     &     DCPAX2, DCPAY2, DCPAZ2, DCPBX2, DCPBY2, DCPBZ2,
     &     DSPCX2, DSPCY2, DSPCZ2, DSPBX2, DSPBY2, DSPBZ2,
     &     DPRIJX2, DPRIJY2, DPRIJZ2, DPRKLX2, DPRKLY2, DPRKLZ2,
     &     DPRJKX2, DPRJKY2, DPRJKZ2
      DOUBLE PRECISION OBSJ, ERRJ, K1
C begin
C
C following Axel's code in ETOR,
C
C zero out partial energy
C
      EJ = ZERO
C
      CLASS = 1
      K1 = ONEFORCES(1)
      A = ONEAS(1)
      B = ONEBS(1)
      C = ONECS(1)
      D = ONEDS(1)
      P1 = ONEPHASE1S(1)
      P2 = ONEPHASE2S(1)
      P3 = ONEPHASE3S(1)
      DO COUNT = 1, NONES
          IF (ONEASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               K1 = ONEFORCES(CLASS)
               A = ONEAS(CLASS)
               B = ONEBS(CLASS)
               C = ONECS(CLASS)
               D = ONEDS(CLASS)
               P1 = ONEPHASE1S(CLASS)
               P2 = ONEPHASE2S(CLASS)
               P3 = ONEPHASE3S(CLASS)
          END IF
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
          XN = X(ATOMN(COUNT))
          XP = X(ATOMP(COUNT))
          XQ = X(ATOMQ(COUNT))
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
          YN = Y(ATOMN(COUNT))
          YP = Y(ATOMP(COUNT))
          YQ = Y(ATOMQ(COUNT))
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
          ZN = Z(ATOMN(COUNT))
          ZP = Z(ATOMP(COUNT))
          ZQ = Z(ATOMQ(COUNT))
C
          OBSJ = JOBS(COUNT)
          ERRJ = JERR(COUNT)
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
          XIJ = XI - XJ
          XJK = XJ - XK
          XKL = XK - XL
          YIJ = YI - YJ
          YJK = YJ - YK
          YKL = YK - YL
          ZIJ = ZI - ZJ
          ZJK = ZJ - ZK
          ZKL = ZK - ZL
C
          XIJ2 = XM - XN
          XJK2 = XN - XP
          XKL2 = XP - XQ
          YIJ2 = YM - YN
          YJK2 = YN - YP
          YKL2 = YP - YQ
          ZIJ2 = ZM - ZN
          ZJK2 = ZN - ZP
          ZKL2 = ZP - ZQ
C
C now calculate A=RIJ*RJK, B = RJK*RKL, C = RJK*(RIJ*RJK)
C
          AX = YIJ*ZJK-ZIJ*YJK
          AY = ZIJ*XJK-XIJ*ZJK
          AZ = XIJ*YJK-YIJ*XJK
          BX = YJK*ZKL-YKL*ZJK
          BY = ZJK*XKL-ZKL*XJK
          BZ = XJK*YKL-XKL*YJK
C
          CX = YJK*AZ-ZJK*AY
          CY = ZJK*AX-XJK*AZ
          CZ = XJK*AY-YJK*AX
C
          AX2 = YIJ2*ZJK2-ZIJ2*YJK2
          AY2 = ZIJ2*XJK2-XIJ2*ZJK2
          AZ2 = XIJ2*YJK2-YIJ2*XJK2
          BX2 = YJK2*ZKL2-YKL2*ZJK2
          BY2 = ZJK2*XKL2-ZKL2*XJK2
          BZ2 = XJK2*YKL2-XKL2*YJK2
          CX2 = YJK2*AZ2-ZJK2*AY2
          CY2 = ZJK2*AX2-XJK2*AZ2
          CZ2 = XJK2*AY2-YJK2*AX2
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
          CP2 = AX2*BX2+AY2*BY2+AZ2*BZ2
          SP2 = CX2*BX2+CY2*BY2+CZ2*BZ2
C
C calculate phi (make sure cos is within bounds and get sign from sin)
C and keep it it radians
C
          PHI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
          PSI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP2))),SP2)
C
C calculate energy as a function of delta J
C
          CALCJ =  A+B*SIN(PSI+P1)+C*COS(2*(PSI+P2))+D*COS(2*(PHI+P3))
          IF (WHICH.EQ.'ANALYZE') JCALC(COUNT) = CALCJ
          DELTAJ = CALCJ-OBSJ
          IF ((ABS(DELTAJ).LT.ERRJ).AND.(POTENTIAL.EQ.SQUARE)) THEN
               E = 0
               DF = 0
          ELSE
               E = K1*(DELTAJ**2)
               DF = -4*D*K1*DELTAJ*SIN(2*(PHI + P3))
               DF2 = 2*K1*DELTAJ*(B*COS(PSI+P1)-2*C*SIN(2*(PSI+P2)))
          END IF
C
C accumulate energy
C
          EJ = EJ + E
C
C compute heavyside function
C
          SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
C
          SWITCH2=-MIN(1,INT(ABS(SP2)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
          RECSP=DF*SWITCH*SIGN(ONE,SP)/MAX(ABS(SP),MCONST)
          RECCP=DF*(SWITCH+1)*SIGN(ONE,CP)/MAX(ABS(CP),MCONST)
C
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
     &     RECSP2*(YJK2*DCPAZ2-DCPAY2*ZJK2)+
     &     RECCP2*((YJK2**2+ZJK2**2)*DSPCX2-XJK2*YJK2*DSPCY2-XJK2
     &         *ZJK2*DSPCZ2)
          DPRIJY2=
     &     RECSP2*(ZJK2*DCPAX2-DCPAZ2*XJK2)+
     &     RECCP2*((ZJK2**2+XJK2**2)*DSPCY2-YJK2*ZJK2*DSPCZ2-YJK2
     &         *XJK2*DSPCX2)
          DPRIJZ2=
     &     RECSP2*(XJK2*DCPAY2-DCPAX2*YJK2)+
     &     RECCP2*((XJK2**2+YJK2**2)*DSPCZ2-ZJK2*XJK2*DSPCX2-ZJK2
     &         *YJK2*DSPCY2)
C
          DPRKLX2=
     &     RECSP2*(DCPBY2*ZJK2-YJK2*DCPBZ2)+
     &     RECCP2*(DSPBY2*ZJK2-YJK2*DSPBZ2)
          DPRKLY2=
     &     RECSP2*(DCPBZ2*XJK2-ZJK2*DCPBX2)+
     &     RECCP2*(DSPBZ2*XJK2-ZJK2*DSPBX2)
          DPRKLZ2=
     &     RECSP2*(DCPBX2*YJK2-XJK2*DCPBY2)+
     &     RECCP2*(DSPBX2*YJK2-XJK2*DSPBY2)
C
          DPRJKX2=
     &     RECSP2*(DCPAY2*ZIJ2-DCPAZ2*YIJ2+DCPBZ2*YKL2-DCPBY2*ZKL2)+
     &     RECCP2*(-(YJK2*YIJ2+ZJK2*ZIJ2)*DSPCX2
     &        +(TWO*XJK2*YIJ2-XIJ2*YJK2)*DSPCY2
     &        +(TWO*XJK2*ZIJ2-XIJ2*ZJK2)*DSPCZ2
     &        +DSPBZ2*YKL2-DSPBY2*ZKL2)
          DPRJKY2=
     &     RECSP2*(DCPAZ2*XIJ2-DCPAX2*ZIJ2+DCPBX2*ZKL2-DCPBZ2*XKL2)+
     &     RECCP2*(-(ZJK2*ZIJ2+XJK2*XIJ2)*DSPCY2
     &        +(TWO*YJK2*ZIJ2-YIJ2*ZJK2)*DSPCZ2
     &        +(TWO*YJK2*XIJ2-YIJ2*XJK2)*DSPCX2
     &        +DSPBX2*ZKL2-DSPBZ2*XKL2)
          DPRJKZ2=
     &     RECSP2*(DCPAX2*YIJ2-DCPAY2*XIJ2+DCPBY2*XKL2-DCPBX2*YKL2)+
     &     RECCP2*(-(XJK2*XIJ2+YJK2*YIJ2)*DSPCZ2
     &        +(TWO*ZJK2*XIJ2-ZIJ2*XJK2)*DSPCX2
     &        +(TWO*ZJK2*YIJ2-ZIJ2*YJK2)*DSPCY2
     &        +DSPBY2*XKL2-DSPBX2*YKL2)
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
      SUBROUTINE READONEBOND
C
C reads in one bond coupling constant information
C
C by John Kuszewski May 1995
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'onebond.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C i/o
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXONES
      DOUBLE PRECISION K1, CUTOFF, A, B, C, D, P1, P2, P3
      CHARACTER*4 THENAME
C begin
C
C this is used by READONE2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
C reset database only if no couplings have been entered
C ie., this is the first time in the script that
C COUPlings has appeared
C
      IF (ONEIPTR.EQ.0) THEN
           CALL ALLOCONES(0, MAXONES)
      END IF
C
C now read input
C
      CALL PUSEND('1-BOND COUPLINGS>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('1-BOND COUPLINGS>')
           CALL MISCOM('1-BOND COUPLINGS>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-onebond')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURCLASS
                CALL NEXTA4('class name =', THENAME)
                MODE = NEW
                DO COUNT = 1, NCLASSES
                     IF (ONECLASSNAMES(COUNT).EQ.THENAME) THEN
                          MODE = UPDATE
                          CURCLASS = COUNT
                     END IF
                END DO
                IF (MODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXONECLASSES) THEN
                         CALL DSPERR('1-bond-COUP','Too many classes.')
                         CALL DSPERR('1-bond-COUP',
     &                      'Increase MAXONECLASSES and recompile.')
                         CALL WRNDIE(-5, 'READONE',
     &                                'Too many 1J-coupling classes.')
                     END IF
                     NCLASSES = NCLASSES + 1
                     CURCLASS = NCLASSES
                     ONECLASSNAMES(CURCLASS) = THENAME
C
C If this isn't the first class, close off the old class
C
                IF (NCLASSES.GT.1) THEN
                     ONEASSNDX(OLDCLASS) = NONES
                END IF
           END IF
C
C set Karplus-like-curve coefficients for this class
C
           ELSE IF (WD(1:4).EQ.'COEF') THEN
                CALL NEXTF('coefficient A =', A)
                CALL NEXTF('coefficient B =', B)
                CALL NEXTF('coefficient C =', C)
                CALL NEXTF('coefficient D =', D)
                CALL NEXTF('phase of the first term =', P1)
                CALL NEXTF('phase of the second term =', P2)
                CALL NEXTF('phase of the third term =', P3)
C
C     start a default class if there isn't one defined
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                WRITE(PUNIT, '(A, A, A, F8.3, F8.3, F8.3, F8.3)')
     &             'Setting curve coefficients for class ',
     &             ONECLASSNAMES(CURCLASS), ' to ', A, B, C, D
                WRITE(PUNIT, '(A, A, A, F8.3, F8.3, F8.3)')
     &             'And setting curve phases (degrees) for class ',
     &             ONECLASSNAMES(CURCLASS), ' to ', P1, P2, P3
                ONEAS(CURCLASS) = A
                ONEBS(CURCLASS) = B
                ONECS(CURCLASS) = C
                ONEDS(CURCLASS) = D
C
C     switch the value to radians
C
                CALL DG2RAD(P1, P1)
                CALL DG2RAD(P2, P2)
                CALL DG2RAD(P3, P3)
                ONEPHASE1S(CURCLASS) = P1
                ONEPHASE2S(CURCLASS) = P2
                ONEPHASE3S(CURCLASS) = P3
C
C set force constant for current class
C
           ELSE IF (WD(1:4).EQ.'FORC') THEN
                CALL NEXTF('force constant =', K1)
C
C start a default class if there isn't one defined
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                WRITE(PUNIT, '(A, A, A, F8.3)')
     &             'Setting force const for class ',
     &             ONECLASSNAMES(CURCLASS), ' to ', K1
                ONEFORCES(CURCLASS) = K1
C
C reset couplings database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
                CALL ONEHP
                CALL ONEINIT
                CALL ALLOCONES(0, MAXONES)
C
C set potential type
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                     WRITE(PUNIT, '(A)')
     &                     'using square well potential.'
                     POTENTIAL = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                     WRITE(PUNIT, '(A)') 'using harmonic potential.'
                     POTENTIAL = HARMONIC
                ELSE
                     CALL DSPERR('ONEBOND',
     &                    'unknown potential. Using square.')
                     POTENTIAL = SQUARE
                END IF
C
C
C change number of assignment slots
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXONES = MAXONES
                CALL NEXTI('number of slots =', MAXONES)
                CALL ALLOCONES(OLDMAXONES, MAXONES)
C
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more coupling assignments
C than you have slots for
C
                IF (NONES.EQ.MAXONES) THEN
                     CALL DSPERR('ONEBOND','Too many assignments.')
                     CALL DSPERR('ONEBOND',
     &                        'Increase NREStraints and run again.')
                     CALL WRNDIE(-1,'ONEBOND>',
     &               'exceeded allocation for onebond J-restraints')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                CALL READONE2(HEAP(ONEIPTR), HEAP(ONEJPTR),
     &               HEAP(ONEKPTR), HEAP(ONELPTR), HEAP(ONEMPTR),
     &               HEAP(ONENPTR), HEAP(ONEPPTR), HEAP(ONEQPTR),
     &               HEAP(ONEJOBSPTR), HEAP(ONEJERRPTR), HEAP(SPTR))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                     CALL DSPERR('1-BOND-COUP',
     &                       'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                          CALL DSPERR('1-BOND-COUP',
     &                         'cutoff must be positive.')
                          CUTOFF = ABS(CUTOFF)
                     END IF
                     CALL PRINTONES(CUTOFF, HEAP(CALCONEPTR),
     &                    HEAP(ONEJOBSPTR), HEAP(ONEIPTR),
     &                    HEAP(ONEJPTR), HEAP(ONEKPTR), HEAP(ONELPTR),
     &                    HEAP(ONEMPTR), HEAP(ONENPTR),
     &                    HEAP(ONEPPTR), HEAP(ONEQPTR))
                END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('1-BOND COUPLINGS>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCONES (OLDSIZE, NEWSIZE)
C
C resets coupling constant arrays to hold SIZE entries
C
C by John Kuszewski July 1993
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'funct.inc'
      INCLUDE 'onebond.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(ONEIPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEJPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEKPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONELPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEMPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONENPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEPPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEQPTR, INTEG4(OLDSIZE))
           CALL FREHP(ONEJOBSPTR, IREAL8(OLDSIZE))
           CALL FREHP(ONEJERRPTR, IREAL8(OLDSIZE))
           CALL FREHP(CALCONEPTR, IREAL8(OLDSIZE))
      END IF
C
      ONEIPTR = 0
      ONEJPTR = 0
      ONEKPTR = 0
      ONELPTR = 0
      ONEMPTR = 0
      ONENPTR = 0
      ONEPPTR = 0
      ONEQPTR = 0
      ONEJOBSPTR = 0
      ONEJERRPTR = 0
      CALCONEPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        ONEIPTR = ALLHP(INTEG4(NEWSIZE))
        ONEJPTR = ALLHP(INTEG4(NEWSIZE))
        ONEKPTR = ALLHP(INTEG4(NEWSIZE))
        ONELPTR = ALLHP(INTEG4(NEWSIZE))
        ONEMPTR = ALLHP(INTEG4(NEWSIZE))
        ONENPTR = ALLHP(INTEG4(NEWSIZE))
        ONEPPTR = ALLHP(INTEG4(NEWSIZE))
        ONEQPTR = ALLHP(INTEG4(NEWSIZE))
        ONEJOBSPTR = ALLHP(IREAL8(NEWSIZE))
        ONEJERRPTR = ALLHP(IREAL8(NEWSIZE))
        CALCONEPTR = ALLHP(IREAL8(NEWSIZE))
      END IF
C
      RETURN
      END
C==============
      SUBROUTINE ONEDEFAULTS
C
C sets up defaults
C
C by John Kuszewski July 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'onebond.inc'
C local variables
      INTEGER COUNT
      DOUBLE PRECISION P
C begin
      MODE = NEW
      MAXONES = 200
      NONES = 0
      NCLASSES = 0
      CURCLASS = 0
      POTENTIAL = SQUARE
      ONEIPTR = 0
      ONEJPTR = 0
      ONEKPTR = 0
      ONELPTR = 0
      ONEMPTR = 0
      ONENPTR = 0
      ONEPPTR = 0
      ONEQPTR = 0
      ONEJOBSPTR = 0
      ONEJERRPTR = 0
      CALCONEPTR = 0
      DO COUNT = 1, MAXONECLASSES
           ONECLASSNAMES(COUNT) = 'DEFAULT'
           ONEASSNDX(COUNT) = 0
           ONEFORCES (COUNT) = 50.0D0
C
C these values are for the 1JHaCa coupling
C Vuister, Delaglio, and Bax, JACS 114, 9674 (1992)
C
           ONEAS(COUNT) = 140.3D0
           ONEBS(COUNT) = 1.4D0
           ONECS(COUNT) = -4.1D0
           ONEDS(COUNT) = 2.0D0
C
C switch the value to radians
C
           CALL DG2RAD(138.0D0, P)
           ONEPHASE1S(COUNT) = P
           CALL DG2RAD(138.0D0, P)
           ONEPHASE2S(COUNT) = P
           CALL DG2RAD(30.0D0, P)
           ONEPHASE3S(COUNT) = P
      END DO
      RETURN
      END
C==============
      SUBROUTINE READONE2 (ATOMI, ATOMJ, ATOMK, ATOML,
     &               ATOMM, ATOMN, ATOMP, ATOMQ, JOBS, JERR, SEL)
C
C reads actual J coupling assignments into arrays
C
C by John Kuszewski July 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'onebond.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     ATOMM(*), ATOMN(*), ATOMP(*), ATOMQ(*), SEL(*)
      DOUBLE PRECISION JOBS(*), JERR(*)
C local variables
      INTEGER NSEL, INSERTPOS, COUNT, CURSTOP, OTHERSTOP
      DOUBLE PRECISION JO, JE
C begin
C
C if we're in update mode, make a space for the new line
C
      IF (MODE.EQ.UPDATE) THEN
           DO COUNT = NONES+1, ONEASSNDX(CURCLASS)+1, -1
                ATOMI(COUNT) = ATOMI(COUNT-1)
                ATOMJ(COUNT) = ATOMJ(COUNT-1)
                ATOMK(COUNT) = ATOMK(COUNT-1)
                ATOML(COUNT) = ATOML(COUNT-1)
                ATOMM(COUNT) = ATOMM(COUNT-1)
                ATOMN(COUNT) = ATOMN(COUNT-1)
                ATOMP(COUNT) = ATOMP(COUNT-1)
                ATOMQ(COUNT) = ATOMQ(COUNT-1)
           END DO
           CURSTOP = ONEASSNDX(CURCLASS)
           DO COUNT = 1, NCLASSES
                OTHERSTOP = ONEASSNDX(COUNT)
                IF (OTHERSTOP.GT.CURSTOP) THEN
                     ONEASSNDX(COUNT) = OTHERSTOP + 1
                END IF
           END DO
           ONEASSNDX(CURCLASS) = CURSTOP + 1
           INSERTPOS = CURSTOP
           NONES = NONES + 1
      ELSE
           NONES = NONES + 1
           INSERTPOS = NONES
           ONEASSNDX(CURCLASS) = INSERTPOS
      END IF
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom i. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMI(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom j. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMJ(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom k. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMK(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom l. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOML(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom m. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMM(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom n. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMN(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom p. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMP(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
           CALL DSPERR('1-bond couplings',
     &     'more than 1 atom in selection for atom q. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMQ(INSERTPOS) = SEL(1)
C
      CALL NEXTF('observed J =', JO)
      CALL NEXTF('error in J =', JE)
      JOBS(INSERTPOS) = JO
      JERR(INSERTPOS) = JE
      RETURN
      END
C===============
      SUBROUTINE ONEINIT
C
C initializes 1J-coupling stuff
C
C by John Kuszewski Aub 1993
C================
      IMPLICIT NONE
C include files
      INCLUDE 'onebond.inc'
C begin
      CALL ONEDEFAULTS
      RETURN
      END
C===============
      SUBROUTINE ONEHP
C
C deallocate 1J-coupling stuff
C
C by Alexandre Bonvin Apr 1996
C================
      IMPLICIT NONE
C include files
      INCLUDE 'onebond.inc'
C begin
      IF (ONEIPTR.NE.0) THEN
      CALL ALLOCONES(MAXONES, 0)
      END IF
      RETURN
      END
C=================
      SUBROUTINE PRINTONES (CUTOFF, JCALC, JOBS,
     &               ATOMI, ATOMJ, ATOMK, ATOML,
     &               ATOMM, ATOMN, ATOMP, ATOMQ)
C
C prints couplings with deltaJ greater than cutoff
C calculates RMS deviation and puts it into $RESULT
C
C by John Kuszewski Aug 1993
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'onebond.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, JCALC(*), JOBS(*), DBPREC
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     ATOMM(*), ATOMN(*), ATOMP(*), ATOMQ(*)
C local variables
      DOUBLE PRECISION CALCJ, OBSJ, DELTAJ, JENERGY, RMS, VIOLS
      INTEGER COUNT, CLASS, I, J, K, L, M, N, P, Q
      DOUBLE COMPLEX DUMMY2
C begin
      RMS = ZERO
      VIOLS = ZERO
C
C make sure that the calcJ array is up to date
C
      CALL EONEBOND(JENERGY, 'ANALYZE')
      WRITE (PUNIT, '(A)') 'The following one-bond couplings have'
      WRITE (PUNIT, '(A)') 'delta J greater than the cutoff:'
      WRITE (PUNIT, '(A)') '  (calculated J) (observed J) (delta J)'
C
C write out first class heading
C
      CLASS = 1
      WRITE (PUNIT, '(A, A)') 'class ', ONECLASSNAMES(1)
C
C for every coupling entry,
C
      DO COUNT = 1, NONES
C
C is this the start of a new class?
C
          IF (ONEASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               WRITE (PUNIT, '(A, A)') 'class ', ONECLASSNAMES(CLASS)
          END IF
C
C always update RMS delta J
C
          CALCJ = JCALC(COUNT)
          OBSJ = JOBS(COUNT)
          DELTAJ = CALCJ-OBSJ
          RMS = RMS + DELTAJ**2
C
C print out delta Js greater than cutoff
C and update number of violations
C
          IF (ABS(DELTAJ).GT.CUTOFF) THEN
               I = ATOMI(COUNT)
               J = ATOMJ(COUNT)
               K = ATOMK(COUNT)
               L = ATOML(COUNT)
               M = ATOMM(COUNT)
               N = ATOMN(COUNT)
               P = ATOMP(COUNT)
               Q = ATOMQ(COUNT)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(I),RESID(I),RES(I),TYPE(I)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(J),RESID(J),RES(J),TYPE(J)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(K),RESID(K),RES(K),TYPE(K)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(L),RESID(L),RES(L),TYPE(L)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(M),RESID(M),RES(M),TYPE(M)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(N),RESID(N),RES(N),TYPE(N)
               WRITE(PUNIT,'(4(A))')
     &              SEGID(P),RESID(P),RES(P),TYPE(P)
               WRITE(PUNIT,'(4(A), 3(F8.3))')
     &              SEGID(Q),RESID(Q),RES(Q),TYPE(Q),
     &              CALCJ, OBSJ, DELTAJ
               VIOLS = VIOLS + ONE
          END IF
      END DO
      IF (NONES.GT.ZERO) THEN
           RMS = SQRT(RMS / NONES)
      ELSE
           RMS = ZERO
      END IF
      WRITE(PUNIT,'(A,F8.3,A,F5.2,A,F6.0,A,I6,A)')
     & '  RMS diff. =',RMS,
     & ', #(violat.>',CUTOFF,')=',VIOLS,
     & ' of ',NONES,' J-couplings'
      DBPREC = NONES
      CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, DBPREC)
      CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
      CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, RMS)
      CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      RETURN
      END
C
      SUBROUTINE SCRONEB
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'onebond.inc'
C
C
C update one-bond coupling database
      IF (ONEIPTR.NE.0) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: one-bond coupling database erased.'
      CALL ONEHP
      CALL ONEINIT
      END IF
      RETURN
      END
C
