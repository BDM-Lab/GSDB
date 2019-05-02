C===============
      SUBROUTINE EVEAN (EV, WHICH)
C     
C     Calls EVEAN2, which does the actual energy calculation
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C===============
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'vectangl.inc'
      INCLUDE 'heap.inc'
C     i/o
      DOUBLE PRECISION EV
      CHARACTER*7 WHICH
C     begin
C     
C     
      CALL EVEAN2(EV, HEAP(VEANIPTR), HEAP(VEANJPTR), HEAP(VEANKPTR),
     &     HEAP(VEANLPTR),
     &     HEAP(VEANVOBSPTR), HEAP(VEANVERRPTR), HEAP(CALCVEANPTR),
     &     WHICH, HEAP(VEANCV))
      RETURN
      END
C===============
      SUBROUTINE EVEAN2 (EV, ATOMI, ATOMJ, ATOMK, ATOML,
     &     VOBS, VERR, VCALC, WHICH, VCV)
C     
C     Calculates vectorangle constant energies
C     
C     energies are of the form
C     E = kbord*deltaV**2 or E*kcent*(1+cos(deltaV / delta))
C     where
C     kbord = energy constant for angles at border,
C     kcent = energy constant for angles in center,
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C     
C===============
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'vectangl.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'comand.inc'
C     
      INCLUDE 'mtf.inc'
C     i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*), VCV(*)
      DOUBLE PRECISION VOBS(2,*), VERR(2,*), VCALC(2,*)
      DOUBLE PRECISION EV
      CHARACTER*7 WHICH
C     local variables
      INTEGER COUNT, CLASS
      DOUBLE PRECISION Kveancent, Kveanbord, Ealpha, Dalpha
      DOUBLE PRECISION XI, XJ, XK, XL, YI, YJ, YK, YL, ZI, ZJ, ZK, ZL
      DOUBLE PRECISION XIJ, XKL, YIJ, YKL, ZIJ, ZKL, LIJ, LKL
      DOUBLE PRECISION COSALPHA, SINALPHA, RECSIN, ALPHA, VOBSA
      DOUBLE PRECISION VOBSB, VOBSC, VOBSD, DV
      DOUBLE PRECISION DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ
C     begin
C     
C     following Axel's code in ETOR,
C     
C     zero out partial energy
C     
      EV = ZERO
C     
      CLASS = 1
      Kveancent = VEANFORCES(2,CLASS)
      Kveanbord = VEANFORCES(1,CLASS)
      DO COUNT = 1, NVEANS
         IF (VEANASSNDX(CLASS).LT.COUNT) THEN
            CLASS = CLASS + 1
            Kveancent = VEANFORCES(1,CLASS)
            Kveanbord = VEANFORCES(2,CLASS)
         END IF
C     
         XI = X(ATOMI(COUNT))
         XJ = X(ATOMJ(COUNT))
         XK = X(ATOMK(COUNT))
         XL = X(ATOML(COUNT))
C     
         YI = Y(ATOMI(COUNT))
         YJ = Y(ATOMJ(COUNT))
         YK = Y(ATOMK(COUNT))
         YL = Y(ATOML(COUNT))
C     
         ZI = Z(ATOMI(COUNT))
         ZJ = Z(ATOMJ(COUNT))
         ZK = Z(ATOMK(COUNT))
         ZL = Z(ATOML(COUNT))
C     
         VOBSA = (VOBS(1,COUNT) - VERR(1,COUNT)) * PI / 180
         VOBSB = (VOBS(1,COUNT) + VERR(1,COUNT)) * PI / 180
         VOBSC = (VOBS(2,COUNT) - VERR(2,COUNT)) * PI / 180
         VOBSD = (VOBS(2,COUNT) + VERR(2,COUNT)) * PI / 180
         if(VOBSC.lt.VOBSB) then
            VOBSC = PI/2
            VOBSB = PI/2
         end if
C     
C     now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C     
         XIJ = XI - XJ
         XKL = XK - XL
C     
         YIJ = YI - YJ
         YKL = YK - YL
C     
         ZIJ = ZI - ZJ
         ZKL = ZK - ZL
C     
C     calculate the norm of A, B, & C & set to MCONST if it's too small
C     
         LIJ = ONE/SQRT(MAX(MCONST, XIJ**2+YIJ**2+ZIJ**2))
         LKL = ONE/SQRT(MAX(MCONST, XKL**2+YKL**2+ZKL**2))
C     
C     normalize A, B, & C
C     
         XIJ = XIJ * LIJ
         YIJ = YIJ * LIJ
         ZIJ = ZIJ * LIJ
         XKL = XKL * LKL
         YKL = YKL * LKL
         ZKL = ZKL * LKL
C     
C     calculate cos(alpha)
C     
         COSALPHA = XIJ*XKL+YIJ*YKL+ZIJ*ZKL
         SINALPHA = SQRT(MAX(ZERO,ONE-COSALPHA**2))
C     
C     calculate alpha (make sure cos is within bounds and get sign from
C     sin) and keep it it radians
C     
         ALPHA = ACOS(MIN(ONE,MAX(-ONE,COSALPHA)))
C     
C     calculate energy and forces for normal entries
C     
         if (ALPHA.lt.VOBSA) then
            Ealpha = Kveanbord*((ALPHA-VOBSA)**2)
            Dalpha = Kveanbord* (ALPHA-VOBSA) *2
         elseif (ALPHA.lt.VOBSB) then
            Ealpha = 0
            Dalpha = 0
         elseif (ALPHA.lt.VOBSC.AND.VOBSB.lt.VOBSC) then
            Ealpha = 0.5*Kveancent*(1+COS(2*PI*(ALPHA-VOBSB)/
     $           (VOBSC-VOBSB)-PI))
            Dalpha = Kveancent*PI/(VOBSB-VOBSC)*
     $           SIN(2*PI*(ALPHA-VOBSB)/(VOBSC-VOBSB)-PI)
         elseif (ALPHA.lt.VOBSD) then
            Ealpha = 0
            Dalpha = 0
         else
            Ealpha = Kveanbord*((ALPHA-VOBSD)**2)
            Dalpha = Kveanbord* (ALPHA-VOBSD) *2
         end if
C     
C     If we're printing out the couplings, then we need to
C     try both possible assignments of observed with the
C     the angle and pick the one with the lower total violation.
C     
         IF (WHICH.EQ.'ANALYZE') VCALC(1,COUNT) = ALPHA
C     
C     accumulate energy (only for working set)
C     
         EV = Ealpha + EV
         DV = Dalpha
C     
C     compute reciprocal of SP multiplied by the derivative of the
C     function DF
         RECSIN = DV/MAX(MCONST,SINALPHA)
C     
C     compute:
C     d (alpha)
C     DF ---------,  etc.
C     d rij
         DPRIJX = RECSIN*LIJ*(COSALPHA*XIJ-XKL)
         DPRIJY = RECSIN*LIJ*(COSALPHA*YIJ-YKL)
         DPRIJZ = RECSIN*LIJ*(COSALPHA*ZIJ-ZKL)
         DPRKLX = RECSIN*LKL*(COSALPHA*XKL-XIJ)
         DPRKLY = RECSIN*LKL*(COSALPHA*YKL-YIJ)
         DPRKLZ = RECSIN*LKL*(COSALPHA*ZKL-ZIJ)
C     
C     now update forces if in energy & force mode
C     
         IF (WHICH.NE.'ANALYZE') THEN
            IF (VCV(COUNT).NE.VICV) THEN
               DX(ATOMI(COUNT))=DX(ATOMI(COUNT))+DPRIJX
               DY(ATOMI(COUNT))=DY(ATOMI(COUNT))+DPRIJY
               DZ(ATOMI(COUNT))=DZ(ATOMI(COUNT))+DPRIJZ
               DX(ATOMJ(COUNT))=DX(ATOMJ(COUNT))-DPRIJX
               DY(ATOMJ(COUNT))=DY(ATOMJ(COUNT))-DPRIJY
               DZ(ATOMJ(COUNT))=DZ(ATOMJ(COUNT))-DPRIJZ
               DX(ATOMK(COUNT))=DX(ATOMK(COUNT))+DPRKLX
               DY(ATOMK(COUNT))=DY(ATOMK(COUNT))+DPRKLY
               DZ(ATOMK(COUNT))=DZ(ATOMK(COUNT))+DPRKLZ
               DX(ATOML(COUNT))=DX(ATOML(COUNT))-DPRKLX
               DY(ATOML(COUNT))=DY(ATOML(COUNT))-DPRKLY
               DZ(ATOML(COUNT))=DZ(ATOML(COUNT))-DPRKLZ
            END IF
         END IF
      END DO
      RETURN
      END
C================
      SUBROUTINE READVEAN
C     
C     reads in vectorangle constant information
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C================
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'vectangl.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'ctitla.inc'
C     i/o
C     local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXVEANS, TEMP
      DOUBLE PRECISION K1, K2, CUTOFF
      CHARACTER*4 THENAME
C     begin
C     
C     this is used by READVEAN2 to hold the selection
C     
      SPTR=ALLHP(INTEG4(NATOM))
C     
C     reset database only if no vectorangles have been entered
C     ie., this is the first time in the xplor script that
C     vectorangles has appeared
C     
      IF (VEANIPTR.EQ.0) THEN
         CALL VEANDEFAULTS
         CALL ALLOCVEANS(0, MAXVEANS)
      END IF
C     
C     now read input
C     
      CALL PUSEND('VEANORANGLES>')
      DO WHILE (.NOT.DONE)
         CALL NEXTWD('VEANORANGLES>')
         CALL MISCOM('VEANORANGLES>',USED)
         IF (.NOT.USED) THEN
C     
            IF (WD(1:4).EQ.'HELP') THEN
               CALL CNSHELP('cns-vectangle')
C     
C     Get class name.  Determine if it's an already-defined class.
C     Insert a new class if it's not.
C     
            ELSE IF (WD(1:4).EQ.'CLAS') THEN
               OLDCLASS = CURCLASS
               CALL NEXTA4('class name =', THENAME)
               MODE = NEW
               DO COUNT = 1, NCLASSES
                  IF (VEANCLASSNAMES(COUNT).EQ.THENAME) THEN
                     MODE = UPDATE
                     CURCLASS = COUNT
                  END IF
               END DO
               IF (MODE.EQ.NEW) THEN
C     
C     make sure you can't add more than the maximum
C     number of classes
C     
                  IF (OLDCLASS.EQ.MAXVEANCLASSES) THEN
                     CALL DSPERR('VEAN','Too many classes.')
                     CALL DSPERR('VEAN',
     &                    'Increase MAXVEANCLASSES and recompile.')
                     CALL WRNDIE(-5, 'READVEAN',
     &                    'Too many vectorangle classes.')
                  END IF
                  NCLASSES = NCLASSES + 1
                  CURCLASS = NCLASSES
                  VEANCLASSNAMES(CURCLASS) = THENAME
C     
C     If this isn't the first class, close off the old class
C     
                  IF (NCLASSES.GT.1) THEN
                     VEANASSNDX(OLDCLASS) = NVEANS
                  END IF
               END IF
C     
C     set force constant for current class
C     
            ELSE IF (WD(1:4).EQ.'FORC') THEN
               CALL NEXTF('border force constant =', K1)
               CALL NEXTF('center force constant =', K2)
C     
C     start a default class if there isn't one defined
C     
               IF (CURCLASS.EQ.0) THEN
                  NCLASSES = 1
                  CURCLASS = 1
               END IF
               WRITE(DUNIT, '(A, A, A, F8.3, F8.3)')
     &              'Setting force consts for class ',
     &              VEANCLASSNAMES(CURCLASS), ' to ', K1, K2
               VEANFORCES(1,CURCLASS) = K1
               VEANFORCES(2,CURCLASS) = K2
C     
C     reset vectorangle database
C     
            ELSE IF (WD(1:4).EQ.'RESE') THEN
               CALL VEANDEFAULTS
               CALL ALLOCVEANS(MAXVEANS, MAXVEANS)
C     
C     change number of assignment slots
C     
            ELSE IF (WD(1:4).EQ.'NRES') THEN
               OLDMAXVEANS = MAXVEANS
               CALL NEXTI('number of slots =', MAXVEANS)
               CALL ALLOCVEANS(OLDMAXVEANS, MAXVEANS)
C     
C     read in an assignment
C     
            ELSE IF (WD(1:4).EQ.'ASSI') THEN
C     
C     make sure you can't add more vectorangle assignments
C     than you have slots for
C     
               IF (NVEANS.EQ.MAXVEANS) THEN
                  CALL DSPERR('VEAN','Too many assignments.')
                  CALL DSPERR('VEAN',
     &                 'Increase NREStraints and run again.')
                  CALL WRNDIE(-1,'VEAN>',
     &                 'exceeded allocation for vectorangle restr.')
               END IF
C     
C     if there isn't a class specified,
C     start a default class
C     
               IF (CURCLASS.EQ.0) THEN
                  NCLASSES = 1
                  CURCLASS = 1
               END IF
               CALL READVEAN2(HEAP(VEANIPTR), HEAP(VEANJPTR),
     &              HEAP(VEANKPTR), HEAP(VEANLPTR),
     &              HEAP(VEANVOBSPTR), HEAP(VEANVERRPTR),
     &              HEAP(SPTR), HEAP(VEANCV))
C     
C     print violations
C     
            ELSE IF (WD(1:4).EQ.'PRIN') THEN
               CALL NEXTWD('PRINt>')
               IF (WD(1:4).NE.'THRE') THEN
                  CALL DSPERR('VEANORANGLES',
     &                 'print expects THREshold parameter.')
               ELSE
                  CALL NEXTF('THREshold =', CUTOFF)
                  IF (CUTOFF.LT.ZERO) THEN
                     CALL DSPERR('VEANORANGLES',
     &                    'cutoff must be positive.')
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
                        IF (VEANCLASSNAMES(COUNT).EQ.
     &                       THENAME) THEN
                           PRINTCLASS(COUNT) = .TRUE.
                        ELSE
                           PRINTCLASS(COUNT) = .FALSE.
                        END IF
                     END DO
                  ELSE
                     CALL DSPERR('VEANORANGLES',
     &                    'not understood.  Printing all classes')
                     DO COUNT = 1,NCLASSES
                        PRINTCLASS(COUNT) = .TRUE.
                     END DO
                  END IF
C     
                  CALL PRINTVEANS(CUTOFF, HEAP(CALCVEANPTR),
     &                 HEAP(VEANVOBSPTR), HEAP(VEANVERRPTR),
     &                 HEAP(VEANIPTR), HEAP(VEANJPTR),
     &                 HEAP(VEANKPTR), HEAP(VEANLPTR),
     &                 HEAP(VEANCV), 0)
                  IF (VICV.GT.0) THEN
                     CALL PRINTVEANS(CUTOFF, HEAP(CALCVEANPTR),
     &                    HEAP(VEANVOBSPTR), HEAP(VEANVERRPTR),
     &                    HEAP(VEANIPTR), HEAP(VEANJPTR),
     &                    HEAP(VEANKPTR), HEAP(VEANLPTR),
     &                    HEAP(VEANCV), 1)
                  END IF
               END IF
C====================================================================
            ELSE IF (WD(1:2).EQ.'CV') THEN
               CALL NEXTI('CV excluded partition number:',VICV)
C====================================================================
            ELSE IF (WD(1:4).EQ.'PART') THEN
               CALL NEXTI('number of PARTitions:',TEMP)
               TEMP=MAX(0,TEMP)
               CALL VCVS(NVEANS,HEAP(VEANCV),TEMP)
               IF (TEMP.EQ.0) THEN
                  VICV=0
               END IF
C     
C     check for END statement
C     
            ELSE
               CALL CHKEND('VEANANGL>', DONE)
            END IF
         END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCVEANS (OLDSIZE, NEWSIZE)
C     
C     resets vectorangle constant arrays to hold SIZE entries
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C===============
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'vectangl.inc'
C     i/o
      INTEGER OLDSIZE, NEWSIZE
C     begin
      IF (OLDSIZE.NE.0) THEN
         CALL FREHP(VEANIPTR, INTEG4(OLDSIZE))
         CALL FREHP(VEANJPTR, INTEG4(OLDSIZE))
         CALL FREHP(VEANKPTR, INTEG4(OLDSIZE))
         CALL FREHP(VEANLPTR, INTEG4(OLDSIZE))
         CALL FREHP(VEANCV  , INTEG4(OLDSIZE))
         CALL FREHP(VEANVOBSPTR, IREAL8(2*OLDSIZE))
         CALL FREHP(VEANVERRPTR, IREAL8(2*OLDSIZE))
         CALL FREHP(CALCVEANPTR, IREAL8(2*OLDSIZE))
      END IF
C
      IF (NEWSIZE.GT.0) THEN     
      VEANIPTR = ALLHP(INTEG4(NEWSIZE))
      VEANJPTR = ALLHP(INTEG4(NEWSIZE))
      VEANKPTR = ALLHP(INTEG4(NEWSIZE))
      VEANLPTR = ALLHP(INTEG4(NEWSIZE))
      VEANCV   = ALLHP(INTEG4(NEWSIZE))
      VEANVOBSPTR = ALLHP(IREAL8(2*NEWSIZE))
      VEANVERRPTR = ALLHP(IREAL8(2*NEWSIZE))
      CALCVEANPTR = ALLHP(IREAL8(2*NEWSIZE))
      END IF
      RETURN
      END
C==============
      SUBROUTINE VEANDEFAULTS
C
C sets up defaults
C
C by Jens Meiler / Michael Nilges Jan 1999
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'vectangl.inc'
C local variables
      INTEGER COUNT
      DOUBLE PRECISION P
C begin
      MODE = NEW
      MAXVEANS = 1000
      NVEANS = 0
      NCLASSES = 0
      CURCLASS = 0
      VICV = 0
      DO COUNT = 1, MAXVEANCLASSES
           VEANCLASSNAMES(COUNT) = 'DEFAULT'
           VEANASSNDX(COUNT) = 0
           VEANFORCES (1,COUNT) = 200
           VEANFORCES (2,COUNT) =  50
      END DO
      RETURN
      END
C==============
      SUBROUTINE READVEAN2 (ATOMI, ATOMJ, ATOMK, ATOML,
     &     VOBS, VERR, SEL, VCV)
C     
C     reads actual vectorangle assignments into arrays
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C==============
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'vectangl.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
C     i/o
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*),
     &     SEL(*), VCV(*)
      DOUBLE PRECISION VOBS(2,*), VERR(2,*)
C     local variables
      INTEGER NSEL, INSERTPOS, COUNT, CURSTOP, OTHERSTOP
      DOUBLE PRECISION JO, JE
C     begin
C     
C     if we're in update mode, make a space for the new line
C     
      IF (MODE.EQ.UPDATE) THEN
         DO COUNT = NVEANS+1, VEANASSNDX(CURCLASS)+1, -1
            ATOMI(COUNT) = ATOMI(COUNT-1)
            ATOMJ(COUNT) = ATOMJ(COUNT-1)
            ATOMK(COUNT) = ATOMK(COUNT-1)
            ATOML(COUNT) = ATOML(COUNT-1)
            VOBS(1,COUNT) = VOBS(1,COUNT-1)
            VOBS(2,COUNT) = VOBS(2,COUNT-1)
            VERR(1,COUNT) = VERR(1,COUNT-1)
            VERR(2,COUNT) = VERR(2,COUNT-1)
            VCV(COUNT) = VCV(COUNT-1)
         END DO
         CURSTOP = VEANASSNDX(CURCLASS)
         DO COUNT = 1, NCLASSES
            OTHERSTOP = VEANASSNDX(COUNT)
            IF (OTHERSTOP.GT.CURSTOP) THEN
               VEANASSNDX(COUNT) = OTHERSTOP + 1
            END IF
         END DO
         VEANASSNDX(CURCLASS) = CURSTOP + 1
         INSERTPOS = CURSTOP
         NVEANS = NVEANS + 1
      ELSE
         NVEANS = NVEANS + 1
         INSERTPOS = NVEANS
         VEANASSNDX(CURCLASS) = INSERTPOS
      END IF
C     
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
         CALL DSPERR('VEAN',
     &        'more than 1 atom in selection for atom i. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMI(INSERTPOS) = SEL(1)
C     
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
         CALL DSPERR('VEAN',
     &        'more than 1 atom in selection for atom j. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMJ(INSERTPOS) = SEL(1)
C     
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
         CALL DSPERR('VEAN',
     &        'more than 1 atom in selection for atom k. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOMK(INSERTPOS) = SEL(1)
C     
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
         CALL DSPERR('VEAN',
     &        'more than 1 atom in selection for atom l. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      ATOML(INSERTPOS) = SEL(1)
C     
      CALL NEXTF('observed vectorangle =', JO)
      CALL NEXTF('error in vectorangle =', JE)
      VOBS(1,INSERTPOS) = JO
      VERR(1,INSERTPOS) = JE
      CALL NEXTF('observed vectorangle =', JO)
      CALL NEXTF('error in vectorangle =', JE)
      VOBS(2,INSERTPOS) = JO
      VERR(2,INSERTPOS) = JE
      VCV(INSERTPOS)  = -1
      RETURN
      END
C===============
      SUBROUTINE VEANINIT
C     
C     initializes vectorangles
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C==============
      IMPLICIT NONE
C     include files
      INCLUDE 'vectangl.inc'
C     begin
      CALL VEANDEFAULTS
      CALL ALLOCVEANS(0, MAXVEANS)
      RETURN
      END
C===============
      SUBROUTINE VEANHP
C
C deallocates vectorangles
C
C by Jens Meiler / Michael Nilges Jan 1999
C================
      IMPLICIT NONE
C include files
      INCLUDE 'vectangl.inc'
C begin
      CALL ALLOCVEANS(MAXVEANS,0)
      RETURN
      END
C=================
      SUBROUTINE PRINTVEANS (CUTOFF, VCALC, VOBS, VERR,
     &     ATOMI, ATOMJ, ATOMK, ATOML,
     &     VCV, ITEST)
C     
C     prints vectorangles with delta alpha greater than cutoff
C     calculates RMS deviation and puts it into $RESULT
C     
C     by Jens Meiler / Michael Nilges Jan 1999
C=================
      IMPLICIT NONE
C     include files
      INCLUDE 'cns.inc'
      INCLUDE 'vectangl.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'consta.inc'
C     i/o
      DOUBLE PRECISION CUTOFF, VCALC(2,*), VOBS(2,*), VERR(2,*)
      INTEGER ATOMI(*), ATOMJ(*), ATOMK(*), ATOML(*)
      INTEGER VCV(*), ITEST
C     local variables
      DOUBLE PRECISION RMS, VIOLS, VENERGY, DELTA, CALCV, VOBSA
      DOUBLE PRECISION VOBSB, VOBSC, VOBSD
      INTEGER NASSIGNSINCLUDED, COUNT, CLASS, I, J, K, L
      LOGICAL PRINTTHISCLASS
      DOUBLE COMPLEX DUMMY2
C     begin
      RMS = ZERO
      VIOLS = ZERO
      NASSIGNSINCLUDED = ZERO
C     
C     make sure that the calcJ array is up to date
C     
      CALL EVEAN(VENERGY, 'ANALYZE')
      IF (VICV.GT.0) THEN
         IF (ITEST.EQ.0) THEN
            WRITE(PUNIT,'(A)')
     &           ' $$$$$$$$$$$$$$$$$$ working set $$$$$$$$$$$$$$$$$$$$$'
         ELSE
            WRITE(PUNIT,'(A,I5,A)')
     &           ' $$$$$$$$$$$$$$$$$$$$ test set (TEST=',VICV,
     &           ')  $$$$$$$$$$$$$$$$$$$$$$'
         END IF
      END IF
      WRITE (PUNIT, '(A)') 'The following vectorangles have delta'
      WRITE (PUNIT, '(A)') 'greater than the cutoff:'
      WRITE (PUNIT, '(A,A)')
     $     '(calculated)(observed borderA borderB borderC borderD)',
     $     '(delta)'
C     
C     write out first class heading
C     
      CLASS = 1
      PRINTTHISCLASS = PRINTCLASS(CLASS)
      IF (PRINTTHISCLASS) THEN
         WRITE (PUNIT, '(A, A)') 'class ', VEANCLASSNAMES(1)
      END IF
C     
C     for every vectorangle entry,
C     
      DO COUNT = 1, NVEANS
C     
C     is this the start of a new class?
C     
         IF (VEANASSNDX(CLASS).LT.COUNT) THEN
            CLASS = CLASS + 1
            PRINTTHISCLASS = PRINTCLASS(CLASS)
            IF (PRINTTHISCLASS) THEN
               WRITE (PUNIT, '(A, A)') 'class ', VEANCLASSNAMES(CLASS)
            END IF
         END IF
C     
C     check if in test set or not (cross-validation)
C     
         IF ((ITEST.EQ.0.AND.VCV(COUNT).NE.VICV).OR.
     &        (ITEST.EQ.1.AND.VCV(COUNT).EQ.VICV)) THEN
C     IF (1.EQ.1) THEN
C     
C     if this ASSIgnment is in a class that will be printed,
C     
            IF (PRINTTHISCLASS) THEN
C     
C     update RMS delta J
C     
               CALCV = VCALC(1,COUNT) * 180 / PI
               VOBSA = (VOBS(1,COUNT) - VERR(1,COUNT))
               VOBSB = (VOBS(1,COUNT) + VERR(1,COUNT))
               VOBSC = (VOBS(2,COUNT) - VERR(2,COUNT))
               VOBSD = (VOBS(2,COUNT) + VERR(2,COUNT))
               if (CALCV.lt.VOBSA) then
                  DELTA = VOBSA - CALCV
               else if (CALCV.lt.VOBSB) then
                  DELTA = 0
               else if (CALCV.lt.VOBSC.AND.VOBSB.lt.VOBSC) then
                  DELTA = MIN(CALCV - VOBSB, VOBSC - CALCV)
               else if (CALCV.lt.VOBSD) then
                  DELTA = 0
               else
                  DELTA = CALCV - VOBSD
               end if
               RMS = RMS + DELTA**2
               NASSIGNSINCLUDED = NASSIGNSINCLUDED + 1
C     
C     print out delta Js greater than cutoff
C     and update number of violations
C     
               IF (ABS(DELTA).GT.CUTOFF) THEN
                  I = ATOMI(COUNT)
                  J = ATOMJ(COUNT)
                  K = ATOMK(COUNT)
                  L = ATOML(COUNT)
                  WRITE(PUNIT,'(9X,16(1X,A), 6(F8.3))')
     &                 SEGID(I),RESID(I),RES(I),TYPE(I),
     &                 SEGID(J),RESID(J),RES(J),TYPE(J),
     &                 SEGID(K),RESID(K),RES(K),TYPE(K),
     &                 SEGID(L),RESID(L),RES(L),TYPE(L),
     &                 CALCV, VOBSA, VOBSB, VOBSC, VOBSD, DELTA
                  VIOLS = VIOLS + ONE
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
     &     '  RMS diff. =',RMS,
     &     ', #(violat.>',CUTOFF,')=',VIOLS,
     &     ' of ',NASSIGNSINCLUDED,' vectorangles'
      IF (ITEST.EQ.1) THEN
         CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, RMS)
         CALL DECLAR('TEST_RMS', 'DP', ' ', DUMMY2, RMS)
         CALL DECLAR('TEST_VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      ELSE
         CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, RMS)
         CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
         CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      END IF
      RETURN
      END
C     
C====================================================================
      SUBROUTINE VCVS(VNUM,VCV,PART)
C     
C     Routine partitions vectorangle data into PART sets.
C     VCV will contain integer numbers between 1 and PART.
C     
C     Author: Axel T. Brunger
C     Modifed for vectorangles: Jens Meiler / Michael Nilges Jan 1999
C     
C     
C     IMPLICIT NONE
C     I/O
      INTEGER VNUM, VCV(*)
      INTEGER PART
C     local
      INTEGER I, P, NP, NRETRY, NTRYTOT
      DOUBLE PRECISION RNUM
      NRETRY = 0
      NTRYTOT = 0
C     begin
 100  CONTINUE
      IF (PART.GT.0) THEN
         DO I=1,VNUM
            CALL GGUBFS(RNUM)
            VCV(I)=MAX(1,MIN(PART,INT(RNUM*PART)+1))
         END DO
C     
         IF (PART.EQ.VNUM) THEN
            DO I=1,VNUM
               VCV(I)=I
            END DO
         END IF
C     
         DO P=1,PART
            NP=0
            DO I=1,VNUM
               IF (VCV(I).EQ.P) THEN
                  NP=NP+1
               END IF
            END DO
            IF (NP .EQ. 0) THEN
               NRETRY = 1
            ENDIF
            WRITE(6,'(A,I3,A,I5,A)') ' For set ',P,
     &           ' there are ',NP,' vectorangle restraints.'
         END DO
      ELSE
         WRITE(6,'(A)')
     &        ' Data are not partitioned or partitioning removed.'
         DO I=1,VNUM
            VCV(I)=-1
         END DO
      END IF
      NTRYTOT = NTRYTOT + NRETRY
      IF (NTRYTOT .GT. 0 .AND. NTRYTOT .LE. 10) THEN
         WRITE(6,'(A)')
     &        ' Test set with 0 constraints! New trial...'
         GOTO 100
      ELSE IF (NTRYTOT .GT. 10) THEN
         CALL WRNDIE(-1,'VCVS',
     &        'Unable to partition the vectorangle data in ten trials')
      ENDIF
C     
      RETURN
      END
C

