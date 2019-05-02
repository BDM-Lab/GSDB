C====================================================
      SUBROUTINE ESANI (ESA, WHICH)
C
C Calls ESANI2, which does the actual energy calculation
C
C by Nico Tjandra and Marius Clore June 1997
C====================================================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'susc_anis.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION ESA
      CHARACTER*7 WHICH
C
      IF (SANIIPTR.NE.0) THEN
      CALL ESANI2(ESA, HEAP(SANIIPTR), HEAP(SANIJPTR), HEAP(SANIKPTR),
     &        HEAP(SANILPTR), HEAP(SANIMPTR), HEAP(SANINPTR),
     &        HEAP(SANIOBSPTR), HEAP(SANIERRPTR),
     &        HEAP(CALCSANIPTR), WHICH)
      END IF
      RETURN
      END
C
C=====================================================
      SUBROUTINE ESANI2 (ESA, ATOMI, ATOMJ, ATOMK, ATOML,
     & ATOMM, ATOMN, SANIOBS, SANIERR, SANICALC, WHICH)
C
C Calculates anisotropy energies
C
C Anis energies are of the form
C      E = k1*deltaAnis**2
C where
C      k1 = force constant,
C      deltaAnis = calculated Anis - observed Anis
C Alos note: that one defines the Anisotropy functions as:
C dip=coef0+coef1*(3*cos(theta)^2-1)+(3/2)*coef2*(sin(theta)^2*cos(2*phi))
C
C which is a flag that switches between energy & force and
C calculated Anis (for violations) calcs
C
C by Nico Tjandra and Marius Clore June 1997
C
C=========================================================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'susc_anis.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*)
      INTEGER ATOMM(*), ATOMN(*)
      DOUBLE PRECISION SANIOBS(*), SANIERR(*), SANICALC(*)
      DOUBLE PRECISION ESA
      CHARACTER*7 WHICH
C local variables
      INTEGER COUNT, CLASS
      DOUBLE PRECISION XI, XJ, XK, XL, XM, XN, YI, YJ, YK, YL, YM, YN,
     &                 ZI, ZJ, ZK, ZL, ZM, ZN,
     &                 CALCSANI, DELTASANI, E, DF1, DF2,
     &                 PHI, DPXI, DPYI, DPZI, DPXJ, DPYJ, DPZJ,
     &                 DPXK, DPYK, DPZK, DPXL, DPYL, DPZL,
     &                 DPXM, DPYM, DPZM, DPXN, DPYN, DPZN,
     &                 THETA, DTXI, DTYI, DTZI, DTXJ, DTYJ, DTZJ,
     &                 DTXK, DTYK, DTZK, DTXL, DTYL, DTZL,
     &                 DTXM, DTYM, DTZM, DTXN, DTYN, DTZN,
     &                 o1, o2, o3, o4 ,o5 ,o6 ,o7 ,o8, o9, o10,
     &                 o11, o12, o13, o14, o15, o16, o17, o18, o19,
     &                 o20, o21, o22, o23, o24, o25, o26, o27, o28,
     &                 o29, o30, o31, o32, o33, o34, o35, o36, o37,
     &                 o38, o39, o40, o41, o42, o43, o44, o45, o46,
     &                 o47, o48, o49, o50, o51, o52, o53, o54, o55,
     &                 o56, o57, o58, o59, o60, o61, o62, o63, o64,
     &                 o65, o66, o67, o68, o69, o70, o71, o72, o73,
     &                 o74, o75, o76, o77, o78, o79, o80, o81, o82
      DOUBLE PRECISION OBSSANI, ERRSANI, K1, COEF0,
     &                 COEF1, COEF2, A, B
C
C begin
C
C following Axel's code in ETOR,
C
C zero out partial energy
C
      ESA = ZERO
C
      CLASS = 1
      K1 = SANIFORCES(1)
      COEF0 = SANICOEF0(1)
      COEF1 = SANICOEF1(1)
      COEF2 = SANICOEF2(1)
C
          COUNT=0
C
      DO WHILE(COUNT.LT.NSANI)
          COUNT = COUNT+1
C
          IF (SANIASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
             IF (SANIASSNDX(CLASS).EQ.SANIASSNDX(CLASS-1)) THEN
               COUNT = COUNT - 1
             ELSE
               K1 = SANIFORCES(CLASS)
               COEF0 = SANICOEF0(CLASS)
               COEF1 = SANICOEF1(CLASS)
               COEF2 = SANICOEF2(CLASS)
            END IF
          END IF
C
          XI = X(ATOMI(COUNT))
          XJ = X(ATOMJ(COUNT))
          XK = X(ATOMK(COUNT))
          XL = X(ATOML(COUNT))
          XM = X(ATOMM(COUNT))
          XN = X(ATOMN(COUNT))
          YI = Y(ATOMI(COUNT))
          YJ = Y(ATOMJ(COUNT))
          YK = Y(ATOMK(COUNT))
          YL = Y(ATOML(COUNT))
          YM = Y(ATOMM(COUNT))
          YN = Y(ATOMN(COUNT))
          ZI = Z(ATOMI(COUNT))
          ZJ = Z(ATOMJ(COUNT))
          ZK = Z(ATOMK(COUNT))
          ZL = Z(ATOML(COUNT))
          ZM = Z(ATOMM(COUNT))
          ZN = Z(ATOMN(COUNT))
C
          OBSSANI = SANIOBS(COUNT)
          ERRSANI = SANIERR(COUNT)
C
C Calculate phi and its derivatives with respect to x,y,z of atoms:
C i, j, k, l, m, and n. This was done with Mathematica !!!
C
C
        o1=-XI+XJ
        o2=o1**2
        o3=-YI+YJ
        o4=o3**2
        o5=-ZI+ZJ
        o6=o5**2
        o7=o2+o4+o6
        o8=sqrt(o7)
        o9=1/o8
        o10=-XM+XN
        o11=o1*o10
        o12=-YM+YN
        o13=o12*o3
        o14=-ZM+ZN
        o15=o14*o5
        o16=o11+o13+o15
        o17=o10**2
        o18=o12**2
        o19=o14**2
        o20=o17+o18+o19
        o21=sqrt(o20)
        o22=1/o21
        o23=1/o7
        o24=o16**2
        o25=1/o20
        o26=o23*o24*o25
        o27=sqrt(1.d0-o26)
        o28=1/o27
        o29=XM-XN
        o30=o7**(3.d0/2.d0)
        o31=1/o30
        o32=o1*o16*o22*o31
        o33=YM-YN
        o34=o16*o22*o3*o31
        o35=ZM-ZN
        o36=o16*o22*o31*o5
        o37=o20**(3.d0/2.d0)
        o38=1/o37
        o39=o10*o16*o38*o9
        o40=o12*o16*o38*o9
        o41=o14*o16*o38*o9
        o42=-XI+XK
        o43=o42**2
        o44=-YI+YK
        o45=o44**2
        o46=-ZI+ZK
        o47=o46**2
        o48=o43+o45+o47
        o49=sqrt(o48)
        o50=-XI+XL
        o51=o50**2
        o52=-YI+YL
        o53=o52**2
        o54=-ZI+ZL
        o55=o54**2
        o56=o51+o53+o55
        o57=sqrt(o56)
        o58=1/o57
        o59=o10*o42
        o60=o12*o44
        o61=o14*o46
        o62=o59+o60+o61
        o63=1/o62
        o64=o10*o50
        o65=o12*o52
        o66=o14*o54
        o67=o64+o65+o66
        o68=o62**2
        o69=1/o68
        o70=o56**(3.d0/2.d0)
        o71=1/o70
        o72=o49*o50*o63*o67*o71
        o73=1/o49
        o74=o42*o58*o63*o67*o73
        o75=1/o56
        o76=o67**2
        o77=o48*o69*o75*o76
        o78=1/(1.d0+o77)
        o79=o49*o52*o63*o67*o71
        o80=o44*o58*o63*o67*o73
        o81=o49*o54*o63*o67*o71
        o82=o46*o58*o63*o67*o73
C
        THETA = acos(o16*o22*o9)
C
        DTXI = -(o28*(o32+o22*o29*o9))
        DTYI = -(o28*(o34+o22*o33*o9))
        DTZI = -(o28*(o36+o22*o35*o9))
        DTXJ = -(o28*(-o32+o10*o22*o9))
        DTYJ = -(o28*(-o34+o12*o22*o9))
        DTZJ = -(o28*(-o36+o14*o22*o9))
        DTXK = 0
        DTYK = 0
        DTZK = 0
        DTXL = 0
        DTYL = 0
        DTZL = 0
        DTXM = -(o28*(o39+o22*o9*(XI-XJ)))
        DTYM = -(o28*(o40+o22*o9*(YI-YJ)))
        DTZM = -(o28*(o41+o22*o9*(ZI-ZJ)))
        DTXN = -(o28*(-o39+o1*o22*o9))
        DTYN =  -(o28*(-o40+o22*o3*o9))
        DTZN = -(o28*(-o41+o22*o5*o9))
C
        PHI = atan(o49*o58*o63*o67)
C
        DPXI = (o29*o49*o58*o63-o29*o49*o58*o67*o69+o72-o74)*o78
        DPYI = o78*(o33*o49*o58*o63-o33*o49*o58*o67*o69+o79-o80)
        DPZI = o78*(o35*o49*o58*o63-o35*o49*o58*o67*o69+o81-o82)
        DPXJ = 0
        DPYJ = 0
        DPZJ = 0
        DPXK = (-(o10*o49*o58*o67*o69)+o74)*o78
        DPYK = o78*(-(o12*o49*o58*o67*o69)+o80)
        DPZK = o78*(-(o14*o49*o58*o67*o69)+o82)
        DPXL = (o10*o49*o58*o63-o72)*o78
        DPYL = o78*(o12*o49*o58*o63-o79)
        DPZL = o78*(o14*o49*o58*o63-o81)
        DPXM = o78*(-(o49*o58*o67*o69*(XI-XK))+o49*o58*o63*(XI-XL))
        DPYM = o78*(-(o49*o58*o67*o69*(YI-YK))+o49*o58*o63*(YI-YL))
        DPZM = o78*(-(o49*o58*o67*o69*(ZI-ZK))+o49*o58*o63*(ZI-ZL))
        DPXN = (o49*o50*o58*o63-o42*o49*o58*o67*o69)*o78
        DPYN = (o49*o52*o58*o63-o44*o49*o58*o67*o69)*o78
        DPZN = (o49*o54*o58*o63-o46*o49*o58*o67*o69)*o78
C
C
C
C
C calculate energy as a function of delta anis
C
C
C
C Casting the coefficients to relate to the correct definition
C COEF0 = DFS
C COEF2 = B/A   i.e. rhombicity
C
C************************
C IF COEF1 is defined as the minimum observed value
C
C COEF1 = -A - 1.5*A*COEF2 + DFS
C and
C       A = (COEF0 - COEF1)/(1. + (1.5*COEF2))
C************************
C
C************************
C IF COEF1 is defined as the magnitude of the axial component
C
C COEF1 = A
C
C
        A = COEF1
        B = A * COEF2
C
        CALCSANI = COEF0+A*(3.*COS(THETA)*COS(THETA)-1.)+
     &            B*(3./2.)*(SIN(THETA)*SIN(THETA)*COS(2.*PHI))
          IF (WHICH.EQ.'ANALYZE') SANICALC(COUNT) = CALCSANI
          DELTASANI = CALCSANI-OBSSANI
          IF (     (ABS(DELTASANI).LT.ERRSANI)
     &        .AND.(POTENTIAL.EQ.SQUARE)) THEN
              E = 0
              DF1 = 0
              DF2 = 0
          ELSE
           E = K1*(DELTASANI**2)
           DF1 = 2*K1*DELTASANI*(A*(-6.)*SIN(THETA)*COS(THETA)
     &                +B*3.*SIN(THETA)*COS(THETA)*COS(2.*PHI))
           DF2 = 2*K1*DELTASANI*(B*(-3.)*SIN(THETA)*SIN(THETA)*
     &           SIN(2.*PHI))
C
          END IF
C
C
C accumulate energy
C
          ESA = ESA + E
C
C now update forces if in energy & force mode
C
          IF (WHICH.NE.'ANALYZE') THEN
           DX(ATOMI(COUNT))=DX(ATOMI(COUNT))+DF1*DTXI+DF2*DPXI
           DY(ATOMI(COUNT))=DY(ATOMI(COUNT))+DF1*DTYI+DF2*DPYI
           DZ(ATOMI(COUNT))=DZ(ATOMI(COUNT))+DF1*DTZI+DF2*DPZI
           DX(ATOMJ(COUNT))=DX(ATOMJ(COUNT))+DF1*DTXJ+DF2*DPXJ
           DY(ATOMJ(COUNT))=DY(ATOMJ(COUNT))+DF1*DTYJ+DF2*DPYJ
           DZ(ATOMJ(COUNT))=DZ(ATOMJ(COUNT))+DF1*DTZJ+DF2*DPZJ
           DX(ATOMK(COUNT))=DX(ATOMK(COUNT))+DF1*DTXK+DF2*DPXK
           DY(ATOMK(COUNT))=DY(ATOMK(COUNT))+DF1*DTYK+DF2*DPYK
           DZ(ATOMK(COUNT))=DZ(ATOMK(COUNT))+DF1*DTZK+DF2*DPZK
           DX(ATOML(COUNT))=DX(ATOML(COUNT))+DF1*DTXL+DF2*DPXL
           DY(ATOML(COUNT))=DY(ATOML(COUNT))+DF1*DTYL+DF2*DPYL
           DZ(ATOML(COUNT))=DZ(ATOML(COUNT))+DF1*DTZL+DF2*DPZL
           DX(ATOMM(COUNT))=DX(ATOMM(COUNT))+DF1*DTXM+DF2*DPXM
           DY(ATOMM(COUNT))=DY(ATOMM(COUNT))+DF1*DTYM+DF2*DPYM
           DZ(ATOMM(COUNT))=DZ(ATOMM(COUNT))+DF1*DTZM+DF2*DPZM
           DX(ATOMN(COUNT))=DX(ATOMN(COUNT))+DF1*DTXN+DF2*DPXN
           DY(ATOMN(COUNT))=DY(ATOMN(COUNT))+DF1*DTYN+DF2*DPYN
           DZ(ATOMN(COUNT))=DZ(ATOMN(COUNT))+DF1*DTZN+DF2*DPZN
          END IF
      END DO
      RETURN
      END
C
C================
      SUBROUTINE READSANI
C
C reads in coupling constant information
C
C by Nico Tjandra and Marius Clore June 1997
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'susc_anis.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C i/o
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXSANI
      DOUBLE PRECISION K1, CUTOFF, COEF0, COEF1, COEF2
      CHARACTER*4 THENAME
C begin
C
C this is used by READSANI2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
C reset database only if no couplings have been entered
C ie., this is the first time that SANIsotropy has been called
C
      IF (SANIIPTR.EQ.0) THEN
           CALL SANIDEFAULTS
           CALL ALLOCSANI(0, MAXSANI)
      END IF
C
C now read input
C
      CALL PUSEND('SANI>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('SANI>')
           CALL MISCOM('SANI>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-sanisotropy')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURCLASS
                CALL NEXTA4('class name =', THENAME)
                MODE = NEW
                DO COUNT = 1, NCLASSES
                     IF (SANICLASSNAMES(COUNT).EQ.THENAME) THEN
                          MODE = UPDATE
                          CURCLASS = COUNT
                     END IF
                END DO
                IF (MODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXSANICLASSES) THEN
                          CALL DSPERR('SANI','Too many classes.')
                          CALL DSPERR('SANI',
     &                      'Increase MAXSANICLASSES and recompile.')
                          CALL WRNDIE(-5, 'READSANI',
     &                      'Too many Anisotropy classes.')
                     END IF
                     NCLASSES = NCLASSES + 1
                     CURCLASS = NCLASSES
                     SANICLASSNAMES(CURCLASS) = THENAME
C
C If this isn't the first class, close off the old class
C
                     IF (NCLASSES.GT.1) THEN
                            SANIASSNDX(OLDCLASS) = NSANI
                     END IF
               END IF
C
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
     &            'Setting force const for class ',
     &            SANICLASSNAMES(CURCLASS), ' to ', K1
                SANIFORCES(CURCLASS) = K1
C
C set coefficient constant for current class
C
           ELSE IF (WD(1:4).EQ.'COEF') THEN
                CALL NEXTF('coefficient 0 =', COEF0)
                CALL NEXTF('coefficient 1 =', COEF1)
                CALL NEXTF('coefficient 2 =', COEF2)
C
C start a default class if there isn't one already defined
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                WRITE(PUNIT, '(A,A,A, F8.3, F8.3, F8.3, F8.3)')
     &            'Setting coefficients for class ',
     &            SANICLASSNAMES(CURCLASS), 'to', COEF0, COEF1, COEF2
C
                SANICOEF0(CURCLASS) = COEF0
                SANICOEF1(CURCLASS) = COEF1
                SANICOEF2(CURCLASS) = COEF2
C
C
C reset couplings database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
                CALL SANIHP
                CALL SANIINIT
                CALL ALLOCSANI(0, MAXSANI)
C
C set potential type
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                     WRITE(PUNIT, '(A)') 'using square well potential.'
                     POTENTIAL = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                     WRITE(PUNIT, '(A)') 'using harmonic potential.'
                     POTENTIAL = HARMONIC
                ELSE
                  CALL DSPERR('SANI','unknown potential. Using square.')
                  POTENTIAL = SQUARE
                END IF
C
C
C change number of assignment slots
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXSANI = MAXSANI
                CALL NEXTI('number of slots =', MAXSANI)
                CALL ALLOCSANI(OLDMAXSANI, MAXSANI)
                WRITE(PUNIT, '(A, I8)')
     &            'SANI: Allocating space for', MAXSANI,
     &            'number of constraints '
C
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more coupling assignments
C than you have slots for
C
                IF (NSANI.EQ.MAXSANI) THEN
                     CALL DSPERR('SANI','Too many assignments.')
                     CALL DSPERR('SANI',
     &                    'Increase NREStraints and run again.')
                     CALL WRNDIE(-1,'SANI>',
     &                    'exceeded allocation for SANI restraints')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURCLASS.EQ.0) THEN
                     NCLASSES = 1
                     CURCLASS = 1
                END IF
                CALL READSANI2(HEAP(SANIIPTR), HEAP(SANIJPTR),
     &                  HEAP(SANIKPTR), HEAP(SANILPTR),
     &                  HEAP(SANIMPTR),
     &                  HEAP(SANINPTR), HEAP(SANIOBSPTR),
     &                  HEAP(SANIERRPTR), HEAP(SPTR))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                  CALL DSPERR('SANI',
     &                        'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                       CALL DSPERR('SANI', 'cutoff must be positive.')
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
                               IF (SANICLASSNAMES(COUNT).EQ.
     &                              THENAME) THEN
                                    PRINTCLASS(COUNT) = .TRUE.
                               ELSE
                                    PRINTCLASS(COUNT) = .FALSE.
                               END IF
                          END DO
                     ELSE
                          DO COUNT = 1,NCLASSES
                               PRINTCLASS(COUNT) = .TRUE.
                          END DO
                     END IF
C
           CALL PRINTSANI(CUTOFF, HEAP(CALCSANIPTR), HEAP(SANIOBSPTR),
     &                 HEAP(SANIIPTR), HEAP(SANIJPTR), HEAP(SANIKPTR),
     &                 HEAP(SANILPTR), HEAP(SANIMPTR), HEAP(SANINPTR))
                END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('SANI>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCSANI (OLDSIZE, NEWSIZE)
C
C resets coupling constant arrays to hold SIZE entries
C
C by Nico Tjandra and Marius Clore June 1997
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'funct.inc'
      INCLUDE 'susc_anis.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(SANIIPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANIJPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANIKPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANILPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANIMPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANINPTR, INTEG4(OLDSIZE))
           CALL FREHP(SANIOBSPTR, IREAL8(OLDSIZE))
           CALL FREHP(SANIERRPTR, IREAL8(OLDSIZE))
           CALL FREHP(CALCSANIPTR, IREAL8(OLDSIZE))
      END IF
C
      SANIIPTR = 0
      SANIJPTR = 0
      SANIKPTR = 0
      SANILPTR = 0
      SANIMPTR = 0
      SANINPTR = 0
      SANIOBSPTR = 0
      SANIERRPTR = 0
      CALCSANIPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        SANIIPTR = ALLHP(INTEG4(NEWSIZE))
        SANIJPTR = ALLHP(INTEG4(NEWSIZE))
        SANIKPTR = ALLHP(INTEG4(NEWSIZE))
        SANILPTR = ALLHP(INTEG4(NEWSIZE))
        SANIMPTR = ALLHP(INTEG4(NEWSIZE))
        SANINPTR = ALLHP(INTEG4(NEWSIZE))
        SANIOBSPTR = ALLHP(IREAL8(NEWSIZE))
        SANIERRPTR = ALLHP(IREAL8(NEWSIZE))
        CALCSANIPTR = ALLHP(IREAL8(NEWSIZE))
      END IF
C
      RETURN
      END
C==============
      SUBROUTINE SANIDEFAULTS
C
C sets up defaults
C
C by Nico Tjandra and Marius Clore June 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'susc_anis.inc'
C local variables
      INTEGER COUNT
C begin
      SANIIPTR = 0
      MODE = NEW
      MAXSANI = 200
      NSANI = 0
      NCLASSES = 0
      CURCLASS = 0
      POTENTIAL = SQUARE
      DO COUNT = 1, MAXSANICLASSES
           SANICLASSNAMES(COUNT) = 'DEFAULT'
           SANIASSNDX(COUNT) = 0
           SANIFORCES (COUNT) = 50
           SANICOEF0(COUNT) = 1.0
           SANICOEF1(COUNT) = 1.0
           SANICOEF2(COUNT) = 1.0
      END DO
      RETURN
      END
C==============
      SUBROUTINE READSANI2 (ATOMI, ATOMJ, ATOMK, ATOML,
     &         ATOMM, ATOMN, SANIOBS, SANIERR, SEL)
C
C reads actual J coupling assignments into arrays
C
C by Nico Tjandra and Marius Clore June 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'susc_anis.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*)
      INTEGER ATOMM(*), ATOMN(*),  SEL(*)
      DOUBLE PRECISION SANIOBS(*), SANIERR(*)
C local variables
      INTEGER NSEL, INSERTPOS, COUNT, CURSTOP, OTHERSTOP, NFLAG
      DOUBLE PRECISION SANIO, SANIE
C begin
C
C if we're in update mode, make a space for the new line
C
      NFLAG = 0
C
      IF (MODE.EQ.UPDATE) THEN
           DO COUNT = NSANI+1, SANIASSNDX(CURCLASS)+1, -1
                ATOMI(COUNT) = ATOMI(COUNT-1)
                ATOMJ(COUNT) = ATOMJ(COUNT-1)
                ATOMK(COUNT) = ATOMK(COUNT-1)
                ATOML(COUNT) = ATOML(COUNT-1)
                ATOMM(COUNT) = ATOMM(COUNT-1)
                ATOMN(COUNT) = ATOMN(COUNT-1)
                SANIOBS(COUNT) = SANIOBS(COUNT-1)
                SANIERR(COUNT) = SANIERR(COUNT-1)
           END DO
           CURSTOP = SANIASSNDX(CURCLASS)
           DO COUNT = 1, NCLASSES
                OTHERSTOP = SANIASSNDX(COUNT)
                IF (OTHERSTOP.GT.CURSTOP) THEN
                     SANIASSNDX(COUNT) = OTHERSTOP + 1
                END IF
           END DO
           SANIASSNDX(CURCLASS) = CURSTOP + 1
           INSERTPOS = CURSTOP
           NSANI = NSANI + 1
      ELSE
           NSANI = NSANI + 1
           INSERTPOS = NSANI
           SANIASSNDX(CURCLASS) = INSERTPOS
      END IF
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom i. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMI(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom j. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMJ(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom k. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMK(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom l. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOML(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom l. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMM(INSERTPOS) = SEL(1)
C
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
       CALL DSPERR('SANI',
     &             'more than 1 atom in sel for atom l. Using first')
      END IF
      IF (NSEL.EQ.0) THEN
       CALL DSPERR('SANI','Error with atom selection')
       NFLAG = 1
      ENDIF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMN(INSERTPOS) = SEL(1)
C
      CALL NEXTF('observed Anis =', SANIO)
      CALL NEXTF('error in Anis =', SANIE)
      SANIOBS(INSERTPOS) =SANIO
      SANIERR(INSERTPOS) = SANIE
C
C Check for error atom selection. If there is one than reset the
C counter for restraint
C
      IF (NFLAG.EQ.1) THEN
       NSANI = NSANI-1
      ENDIF
C
      RETURN
      END
C===============
      SUBROUTINE SANIINIT
C
C initializes J-coupling stuff
C
C by nico Tjandra and Marius Clore June 1997
C================
      IMPLICIT NONE
C include files
      INCLUDE 'susc_anis.inc'
C begin
      CALL SANIDEFAULTS
C
      RETURN
      END
C===============
      SUBROUTINE SANIHP
C
C deallocates susceptibility anisotropy stuff
C
C by Gregory Warren 1/9/98
C================
      IMPLICIT NONE
C include files
      INCLUDE 'susc_anis.inc'
C begin
      IF (SANIIPTR.NE.0) THEN
      CALL ALLOCSANI(MAXSANI,0)
      END IF
      RETURN
      END
C=================
      SUBROUTINE PRINTSANI (CUTOFF, SANICALC, SANIOBS,
     &       ATOMI, ATOMJ, ATOMK, ATOML, ATOMM, ATOMN)
C
C prints anisotropy with deltaJ greater than cutoff
C calculates RMS deviation and puts it into $RESULT
C
C by Nico Tjandra and Marius Clore June 1997
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'susc_anis.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, SANICALC(*), SANIOBS(*)
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*)
      INTEGER ATOMM(*), ATOMN(*)
C local variables
      DOUBLE PRECISION CALCSANI, OBSSANI, DELTASANI, SANIENERGY
      INTEGER COUNT, CLASS, I, J, K, L, M, N, NUM
      DOUBLE PRECISION RMS, VIOLS, DBPREC
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTTHISCLASS
C begin
      RMS = ZERO
      VIOLS = ZERO
      NUM = 0
C
C make sure that the calcJ array is up to date
C
      CALL ESANI(SANIENERGY, 'ANALYZE')
      WRITE (PUNIT, '(A)') 'The following anisotropies have'
      WRITE (PUNIT, '(A)') 'delta SANI greater than the cutoff:'
      WRITE (PUNIT, '(A)') '(calc Anis) (obs Anis) (delta Anis)'
C
C write out first class heading
C
      CLASS = 1
      PRINTTHISCLASS = PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
           WRITE (PUNIT, '(A, A)') 'class ', SANICLASSNAMES(1)
      END IF
C
C for every anisotropy entry,
C
      COUNT = 0
      DO WHILE(COUNT.LT.NSANI)
          COUNT=COUNT+1
C
C is this the start of a new class?
C
          IF (SANIASSNDX(CLASS).LT.COUNT) THEN
               CLASS = CLASS + 1
               PRINTTHISCLASS = PRINTCLASS(CLASS)
               IF (SANIASSNDX(CLASS).EQ.SANIASSNDX(CLASS-1)) THEN
                 COUNT = COUNT - 1
               ENDIF
               IF (PRINTTHISCLASS) THEN
                 WRITE (PUNIT, '(A, A)') 'class ', SANICLASSNAMES(CLASS)
               END IF
          END IF
C
C If this assignment is in a class to be printed
C and make sure there is an entry for that class
C
         IF((PRINTTHISCLASS).AND.(SANIASSNDX(CLASS).NE.
     &    SANIASSNDX(CLASS-1))) THEN
C
C always update RMS delta J
C
          CALCSANI = SANICALC(COUNT)
          OBSSANI = SANIOBS(COUNT)
          DELTASANI = (CALCSANI-OBSSANI)
          RMS = RMS + DELTASANI**2
          NUM = NUM + 1
C
C print out delta Js greater than cutoff
C and update number of violations
C
          IF (ABS(DELTASANI).GT.CUTOFF) THEN
               I = ATOMI(COUNT)
               J = ATOMJ(COUNT)
               K = ATOMK(COUNT)
               L = ATOML(COUNT)
               M = ATOMM(COUNT)
               N = ATOMN(COUNT)
               WRITE(PUNIT,'(5X,8(1X,A), 3(F8.3))')
     &               SEGID(M),RESID(M),RES(M),TYPE(M),
     &               SEGID(N),RESID(N),RES(N),TYPE(N),
     &               CALCSANI, OBSSANI, DELTASANI
               VIOLS = VIOLS + ONE
          END IF
       END IF
      END DO

      IF (NUM.GT.0) THEN
          RMS = SQRT(RMS / NUM)
      ELSE
          RMS=0.0
      ENDIF
      DBPREC = NUM
      CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, DBPREC)
      CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, RMS)
      CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      RETURN
      END
C
      SUBROUTINE SCRASANI
C
C Author: Gregory Warren
      IMPLICIT NONE
C I/O
      INCLUDE 'susc_anis.inc'
C reset diffusion anisotropy database
      IF (SANIIPTR.NE.0) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: susceptibility anistropy database erased.'
      CALL SANIHP
      CALL SANIINIT
      END IF
      RETURN
      END
C
