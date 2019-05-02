C===============
      SUBROUTINE ERAMA (ER, WHICH)
C
C Front end for ERAMA2
C
C by John Kuszewski Aug 1996
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'heap.inc'
C i/o
      DOUBLE PRECISION ER
      CHARACTER*7 WHICH
C begin
      CALL ERAMA2(ER, WHICH, HEAP(RAMAATOMPTR))
      RETURN
      END
C===============
      SUBROUTINE ERAMA2 (ER, WHICH, ATOMS)
C
C calls the appropriate energy term for each rama constraint
C
C by John Kuszewski Jan 1997
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*)
      DOUBLE PRECISION ER
      CHARACTER*7 WHICH
C local vbls
      INTEGER COUNT, CLASS
C begin
C
C zero out partial energy
C
      ER = ZERO
C
C for ever rama constraint,
C
      DO COUNT = 1, NRAMAS
C
C determine its class and do the right thing
C
           CLASS = ATOMS(1,COUNT)
           IF (RAMACLASSTYPE(CLASS).EQ.1) THEN
                CALL ERAMA1D(ER, COUNT, HEAP(RAMAATOMPTR),
     &               HEAP(EXPRAMAPTRS(CLASS)),
     &               HEAP(CALCRAMAPTR), PHISTEPS(CLASS),
     &               WHICH)
           ELSE IF (RAMACLASSTYPE(CLASS).EQ.2) THEN
                CALL ERAMA2D(ER, COUNT, HEAP(RAMAATOMPTR),
     &               HEAP(EXPRAMAPTRS(CLASS)),
     &               HEAP(CALCRAMAPTR), PHISTEPS(CLASS),
     &               PSISTEPS(CLASS), WHICH)
           ELSE IF (RAMACLASSTYPE(CLASS).EQ.3) THEN
                CALL ERAMA3D(ER, COUNT, HEAP(RAMAATOMPTR),
     &               HEAP(EXPRAMAPTRS(CLASS)),
     &               HEAP(CALCRAMAPTR), PHISTEPS(CLASS),
     &               PSISTEPS(CLASS), CHISTEPS(CLASS), WHICH)
           ELSE
                CALL ERAMA4D(ER, COUNT, HEAP(RAMAATOMPTR),
     &               HEAP(EXPRAMAPTRS(CLASS)),
     &               HEAP(CALCRAMAPTR), PHISTEPS(CLASS),
     &               PSISTEPS(CLASS), CHISTEPS(CLASS),
     &               THTSTEPS(CLASS), WHICH)
           END IF
      END DO
      RETURN
      END
C===============
      SUBROUTINE ERAMA4D (ER, COUNT, ATOMS, EXPECTEDRAMA,
     &     CALCRAMA, CURPHIS, CURPSIS, CURCHIS, CURTHTS, WHICH)
C
C Calculates 4D torsion angle database energies
C
C Torsion angle database energies are of the form
C      E = k1*expectedRama
C where
C      k1 = force constant,
C
C WHICH is a flag that switches between energy & force and
C calculated chem shift (for violations) calcs
C
C by John Kuszewski Nov 1996 Mar 1997
C
C===============
      IMPLICIT NONE
C function declarations
      DOUBLE PRECISION TORPOSITION
      EXTERNAL TORPOSITION
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*), CURPHIS, CURPSIS, CURCHIS, CURTHTS,
     &     COUNT
      DOUBLE PRECISION EXPECTEDRAMA(CURPHIS,CURPSIS,CURCHIS,CURTHTS)
      DOUBLE PRECISION CALCRAMA(5,*), ER
      CHARACTER*7 WHICH
C for convenience
      INTEGER HI, LO
      PARAMETER (HI = 2, LO = 1)
C local variables
      INTEGER CLASS, I, J, K, L, M, N, P, Q, R, S, T, U, D, F, G, H,
     &     LOWPHI, HIPHI, LOWPSI, HIPSI, LOWCHI, HICHI,
     &     LOWTHT, HITHT,
     &     COUNT1, COUNT2, COUNT3, COUNT4
      DOUBLE PRECISION PHIWIDTH, PSIWIDTH, CHIWIDTH, THTWIDTH,
     &     K1, ERR,
     &     CORRECT1, CORRECT2, CORRECT3, CORRECT4, CORRECT5, CORRECT6,
     &     CORRECT7, CORRECT8, CORRECT9, CORRECT10, CORRECT11,
     &     CORRECT12,
     &     PHI, PSI, CHI, THT, CURPOS, NEXTPOS,
     &     PHIFRACT, PSIFRACT, CHIFRACT, THTFRACT,
     &     ECURRENT, DEDPHI, DEDPSI, DEDCHI, DEDTHT
      DOUBLE PRECISION NEIGHBORES(2,2,2,2)
      LOGICAL OKPOINT
C
C begin
C
      CLASS = ATOMS(1,COUNT)
C
      K1 = RAMAFORCES(CLASS) * RAMASCALE
      ERR = RAMAERRORS(CLASS)
      CORRECT1 = RAMAPHASEC(1, CLASS)
      CORRECT2 = RAMAPHASEC(2, CLASS)
      CORRECT3 = RAMAPHASEC(3, CLASS)
      CORRECT4 = RAMAPHASEC(4, CLASS)
      CORRECT5 = RAMAPHASEC(5, CLASS)
      CORRECT6 = RAMAPHASEC(6, CLASS)
      CORRECT7 = RAMAPHASEC(7, CLASS)
      CORRECT8 = RAMAPHASEC(8, CLASS)
      CORRECT9 = RAMAPHASEC(9, CLASS)
      CORRECT10 = RAMAPHASEC(10, CLASS)
      CORRECT11 = RAMAPHASEC(11, CLASS)
      CORRECT12 = RAMAPHASEC(12, CLASS)
C
C calculate bin widths
C
      IF (CORRECT2.LT.(PI-CORRECT1)) THEN
           PHIWIDTH = (2*PI - CORRECT2) / CURPHIS
      ELSE IF (CORRECT2.GT.(PI+CORRECT1)) THEN
           PHIWIDTH = CORRECT2 / CURPHIS
      ELSE
           PHIWIDTH = (PI+CORRECT1) / CURPHIS
      END IF
C
      IF (CORRECT5.LT.(PI-CORRECT4)) THEN
           PSIWIDTH = (2*PI - CORRECT5) / CURPSIS
      ELSE IF (CORRECT5.GT.(PI+CORRECT4)) THEN
           PSIWIDTH = CORRECT5 / CURPSIS
      ELSE
           PSIWIDTH = (PI+CORRECT4) / CURPSIS
      END IF
C
      IF (CORRECT8.LT.(PI-CORRECT7)) THEN
           CHIWIDTH = (2*PI - CORRECT8) / CURCHIS
      ELSE IF (CORRECT8.GT.(PI+CORRECT7)) THEN
           CHIWIDTH = CORRECT8 / CURCHIS
      ELSE
           CHIWIDTH = (PI+CORRECT7) / CURCHIS
      END IF
C
      IF (CORRECT11.LT.(PI-CORRECT10)) THEN
           THTWIDTH = (2*PI - CORRECT11) / CURTHTS
      ELSE IF (CORRECT11.GT.(PI+CORRECT10)) THEN
           THTWIDTH = CORRECT11 / CURTHTS
      ELSE
           THTWIDTH = (PI+CORRECT10) / CURTHTS
      END IF
C
C get coords of each atom
C
      I = ATOMS(2, COUNT)
      J = ATOMS(3, COUNT)
      K = ATOMS(4, COUNT)
      L = ATOMS(5, COUNT)
      M = ATOMS(6, COUNT)
      N = ATOMS(7, COUNT)
      P = ATOMS(8, COUNT)
      Q = ATOMS(9, COUNT)
      R = ATOMS(10, COUNT)
      S = ATOMS(11, COUNT)
      T = ATOMS(12, COUNT)
      U = ATOMS(13, COUNT)
      D = ATOMS(14, COUNT)
      F = ATOMS(15, COUNT)
      G = ATOMS(16, COUNT)
      H = ATOMS(17, COUNT)
C
C determine what the torsion angles are (not rounded!)
C
      PHI = TORPOSITION(I, J, K, L, CORRECT1, CORRECT2, CORRECT3)
      PSI = TORPOSITION(M, N, P, Q, CORRECT4, CORRECT5, CORRECT6)
      CHI = TORPOSITION(R, S, T, U, CORRECT7, CORRECT8, CORRECT9)
      THT = TORPOSITION(D, F, G, H, CORRECT10, CORRECT11, CORRECT12)
C
C get the neighboring points
C and the fractional coordinates within each square
C
C Note that the center of bin i is assumed to be given by
C (i-1) binWidth, where binWidth is the number of radians/bin.
C This is so that bin 1 will be centered about 0 deg.
C
      DO COUNT1 = 1, CURPHIS
           CURPOS = (COUNT1 - 1)*PHIWIDTH
           NEXTPOS = (COUNT1)*PHIWIDTH
           IF ((CURPOS.LE.PHI).AND.(NEXTPOS.GT.PHI)) THEN
                LOWPHI = COUNT1
                HIPHI = COUNT1+1
                PHIFRACT = (PHI - CURPOS) / PHIWIDTH
           END IF
      END DO
      IF (HIPHI.GT.CURPHIS) HIPHI = HIPHI - CURPHIS
C
      DO COUNT1 = 1, CURPSIS
           CURPOS = (COUNT1 - 1)*PSIWIDTH
           NEXTPOS = (COUNT1)*PSIWIDTH
           IF ((CURPOS.LE.PSI).AND.(NEXTPOS.GT.PSI)) THEN
                LOWPSI = COUNT1
                HIPSI = COUNT1+1
                PSIFRACT = (PSI - CURPOS) / PSIWIDTH
           END IF
      END DO
      IF (HIPSI.GT.CURPSIS) HIPSI = HIPSI - CURPSIS
C
      DO COUNT1 = 1, CURCHIS
           CURPOS = (COUNT1 - 1)*CHIWIDTH
           NEXTPOS = (COUNT1)*CHIWIDTH
           IF ((CURPOS.LE.CHI).AND.(NEXTPOS.GT.CHI)) THEN
                LOWCHI = COUNT1
                HICHI = COUNT1+1
                CHIFRACT = (CHI - CURPOS) / CHIWIDTH
           END IF
      END DO
      IF (HICHI.GT.CURCHIS) HICHI = HICHI - CURCHIS
C
      DO COUNT1 = 1, CURTHTS
           CURPOS = (COUNT1 - 1)*THTWIDTH
           NEXTPOS = (COUNT1)*THTWIDTH
           IF ((CURPOS.LE.THT).AND.(NEXTPOS.GT.THT)) THEN
                LOWTHT = COUNT1
                HITHT = COUNT1+1
                THTFRACT = (THT - CURPOS) / THTWIDTH
           END IF
      END DO
      IF (HITHT.GT.CURTHTS) HITHT = HITHT - CURTHTS
C
C look up energies of neighbors
C
      NEIGHBORES(HI,HI,HI,HI) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,HICHI,HITHT)
      NEIGHBORES(LO,HI,HI,HI) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,HICHI,HITHT)
      NEIGHBORES(HI,LO,HI,HI) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,HICHI,HITHT)
      NEIGHBORES(LO,LO,HI,HI) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,HICHI,HITHT)
      NEIGHBORES(HI,HI,LO,HI) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,LOWCHI,HITHT)
      NEIGHBORES(LO,HI,LO,HI) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,LOWCHI,HITHT)
      NEIGHBORES(HI,LO,LO,HI) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,LOWCHI,HITHT)
      NEIGHBORES(LO,LO,LO,HI) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,LOWCHI,HITHT)
      NEIGHBORES(HI,HI,HI,LO) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,HICHI,LOWTHT)
      NEIGHBORES(LO,HI,HI,LO) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,HICHI,LOWTHT)
      NEIGHBORES(HI,LO,HI,LO) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,HICHI,LOWTHT)
      NEIGHBORES(LO,LO,HI,LO) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,HICHI,LOWTHT)
      NEIGHBORES(HI,HI,LO,LO) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,LOWCHI,LOWTHT)
      NEIGHBORES(LO,HI,LO,LO) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,LOWCHI,LOWTHT)
      NEIGHBORES(HI,LO,LO,LO) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,LOWCHI,LOWTHT)
      NEIGHBORES(LO,LO,LO,LO) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,LOWCHI,LOWTHT)
C
C are they all defined?
C
      OKPOINT = .TRUE.
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                DO COUNT3 = LO,HI
                     DO COUNT4 = LO,HI
                          IF (NEIGHBORES(COUNT1,COUNT2,COUNT3,COUNT4)
     &                         .GT.NOEXPTEST) THEN
                               OKPOINT = .FALSE.
                          END IF
                     END DO
                END DO
           END DO
      END DO
C
C multiply neighbor energies by force constant
C
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                DO COUNT3 = LO,HI
                     DO COUNT4 = LO,HI
                          NEIGHBORES(COUNT1,COUNT2,COUNT3,COUNT4) =
     &                    NEIGHBORES(COUNT1,COUNT2,COUNT3,COUNT4) * K1
                     END DO
                END DO
           END DO
      END DO
C
C interpolate energy
C
      ECURRENT =
     &     NEIGHBORES(HI,HI,HI,HI) *
     &      PHIFRACT * PSIFRACT * CHIFRACT * THTFRACT +
     &     NEIGHBORES(LO,HI,HI,HI) *
     &      (1-PHIFRACT) * PSIFRACT * CHIFRACT * THTFRACT +
     &     NEIGHBORES(HI,LO,HI,HI) *
     &      PHIFRACT * (1-PSIFRACT) * CHIFRACT * THTFRACT +
     &     NEIGHBORES(LO,LO,HI,HI) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * CHIFRACT * THTFRACT +
     &     NEIGHBORES(HI,HI,LO,HI) *
     &      PHIFRACT * PSIFRACT * (1-CHIFRACT) * THTFRACT +
     &     NEIGHBORES(LO,HI,LO,HI) *
     &      (1-PHIFRACT) * PSIFRACT * (1-CHIFRACT) * THTFRACT +
     &     NEIGHBORES(HI,LO,LO,HI) *
     &      PHIFRACT * (1-PSIFRACT) * (1-CHIFRACT) * THTFRACT +
     &     NEIGHBORES(LO,LO,LO,HI) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * (1-CHIFRACT) * THTFRACT +
     &     NEIGHBORES(HI,HI,HI,LO) *
     &      PHIFRACT * PSIFRACT * CHIFRACT * (1-THTFRACT) +
     &     NEIGHBORES(LO,HI,HI,LO) *
     &      (1-PHIFRACT) * PSIFRACT * CHIFRACT * (1-THTFRACT) +
     &     NEIGHBORES(HI,LO,HI,LO) *
     &      PHIFRACT * (1-PSIFRACT) * CHIFRACT * (1-THTFRACT) +
     &     NEIGHBORES(LO,LO,HI,LO) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * CHIFRACT * (1-THTFRACT) +
     &     NEIGHBORES(HI,HI,LO,LO) *
     &      PHIFRACT * PSIFRACT * (1-CHIFRACT) * (1-THTFRACT) +
     &     NEIGHBORES(LO,HI,LO,LO) *
     &      (1-PHIFRACT) * PSIFRACT * (1-CHIFRACT) * (1-THTFRACT) +
     &     NEIGHBORES(HI,LO,LO,LO) *
     &      PHIFRACT * (1-PSIFRACT) * (1-CHIFRACT) * (1-THTFRACT) +
     &     NEIGHBORES(LO,LO,LO,LO) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * (1-CHIFRACT) * (1-THTFRACT)
C
C interpolate derivs
C
      DEDPHI =
     &     ((NEIGHBORES(HI,HI,HI,HI) - NEIGHBORES(LO,HI,HI,HI)) /
     &     PHIWIDTH) * PSIFRACT * CHIFRACT * THTFRACT +
     &     ((NEIGHBORES(HI,LO,HI,HI) - NEIGHBORES(LO,LO,HI,HI)) /
     &     PHIWIDTH) * (1-PSIFRACT) * CHIFRACT * THTFRACT +
     &     ((NEIGHBORES(HI,HI,LO,HI) - NEIGHBORES(LO,HI,LO,HI)) /
     &     PHIWIDTH) * PSIFRACT * (1-CHIFRACT) * THTFRACT +
     &     ((NEIGHBORES(HI,LO,LO,HI) - NEIGHBORES(LO,LO,LO,HI)) /
     &     PHIWIDTH) * (1-PSIFRACT) * (1-CHIFRACT) * THTFRACT +
     &     ((NEIGHBORES(HI,HI,HI,LO) - NEIGHBORES(LO,HI,HI,LO)) /
     &     PHIWIDTH) * PSIFRACT * CHIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(HI,LO,HI,LO) - NEIGHBORES(LO,LO,HI,LO)) /
     &     PHIWIDTH) * (1-PSIFRACT) * CHIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(HI,HI,LO,LO) - NEIGHBORES(LO,HI,LO,LO)) /
     &     PHIWIDTH) * PSIFRACT * (1-CHIFRACT) * (1-THTFRACT) +
     &     ((NEIGHBORES(HI,LO,LO,LO) - NEIGHBORES(LO,LO,LO,LO)) /
     &     PHIWIDTH) * (1-PSIFRACT) * (1-CHIFRACT) * (1-THTFRACT)
C
      DEDPSI =
     &     ((NEIGHBORES(HI,HI,HI,HI) - NEIGHBORES(HI,LO,HI,HI)) /
     &     PSIWIDTH) * PHIFRACT * CHIFRACT * THTFRACT +
     &     ((NEIGHBORES(LO,HI,HI,HI) - NEIGHBORES(LO,LO,HI,HI)) /
     &     PSIWIDTH) * (1-PHIFRACT) * CHIFRACT * THTFRACT +
     &     ((NEIGHBORES(HI,HI,LO,HI) - NEIGHBORES(HI,LO,LO,HI)) /
     &     PSIWIDTH) * PHIFRACT * (1-CHIFRACT) * THTFRACT +
     &     ((NEIGHBORES(LO,HI,LO,HI) - NEIGHBORES(LO,LO,LO,HI)) /
     &     PSIWIDTH) * (1-PHIFRACT) * (1-CHIFRACT) * THTFRACT +
     &     ((NEIGHBORES(HI,HI,HI,LO) - NEIGHBORES(HI,LO,HI,LO)) /
     &     PSIWIDTH) * PHIFRACT * CHIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(LO,HI,HI,LO) - NEIGHBORES(LO,LO,HI,LO)) /
     &     PSIWIDTH) * (1-PHIFRACT) * CHIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(HI,HI,LO,LO) - NEIGHBORES(HI,LO,LO,LO)) /
     &     PSIWIDTH) * PHIFRACT * (1-CHIFRACT) * (1-THTFRACT) +
     &     ((NEIGHBORES(LO,HI,LO,LO) - NEIGHBORES(LO,LO,LO,LO)) /
     &     PSIWIDTH) * (1-PHIFRACT) * (1-CHIFRACT) * (1-THTFRACT)
C
      DEDCHI =
     &     ((NEIGHBORES(HI,HI,HI,HI) - NEIGHBORES(HI,HI,LO,HI)) /
     &     CHIWIDTH) * PHIFRACT * PSIFRACT * THTFRACT +
     &     ((NEIGHBORES(LO,HI,HI,HI) - NEIGHBORES(LO,HI,LO,HI)) /
     &     CHIWIDTH) * (1-PHIFRACT) * PSIFRACT * THTFRACT +
     &     ((NEIGHBORES(HI,LO,HI,HI) - NEIGHBORES(HI,LO,LO,HI)) /
     &     CHIWIDTH) * PHIFRACT * (1-PSIFRACT) * THTFRACT +
     &     ((NEIGHBORES(LO,LO,HI,HI) - NEIGHBORES(LO,LO,LO,HI)) /
     &     CHIWIDTH) * (1-PHIFRACT) * (1-PSIFRACT) * THTFRACT +
     &     ((NEIGHBORES(HI,HI,HI,LO) - NEIGHBORES(HI,HI,LO,LO)) /
     &     CHIWIDTH) * PHIFRACT * PSIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(LO,HI,HI,LO) - NEIGHBORES(LO,HI,LO,LO)) /
     &     CHIWIDTH) * (1-PHIFRACT) * PSIFRACT * (1-THTFRACT) +
     &     ((NEIGHBORES(HI,LO,HI,LO) - NEIGHBORES(HI,LO,LO,LO)) /
     &     CHIWIDTH) * PHIFRACT * (1-PSIFRACT) * (1-THTFRACT) +
     &     ((NEIGHBORES(LO,LO,HI,LO) - NEIGHBORES(LO,LO,LO,LO)) /
     &     CHIWIDTH) * (1-PHIFRACT) * (1-PSIFRACT) * (1-THTFRACT)
C
      DEDTHT =
     &     ((NEIGHBORES(HI,HI,HI,HI) - NEIGHBORES(HI,HI,HI,LO)) /
     &     THTWIDTH) * PHIFRACT * PSIFRACT * CHIFRACT +
     &     ((NEIGHBORES(HI,HI,LO,HI) - NEIGHBORES(HI,HI,LO,LO)) /
     &     THTWIDTH) * PHIFRACT * PSIFRACT * (1-CHIFRACT) +
     &     ((NEIGHBORES(HI,LO,HI,HI) - NEIGHBORES(HI,LO,HI,LO)) /
     &     THTWIDTH) * PHIFRACT * (1-PSIFRACT) * CHIFRACT +
     &     ((NEIGHBORES(HI,LO,LO,HI) - NEIGHBORES(HI,LO,LO,LO)) /
     &     THTWIDTH) * PHIFRACT * (1-PSIFRACT) * (1-CHIFRACT) +
     &     ((NEIGHBORES(LO,HI,HI,HI) - NEIGHBORES(LO,HI,HI,LO)) /
     &     THTWIDTH) * (1-PHIFRACT) * PSIFRACT * CHIFRACT +
     &     ((NEIGHBORES(LO,HI,LO,HI) - NEIGHBORES(LO,HI,LO,LO)) /
     &     THTWIDTH) * (1-PHIFRACT) * PSIFRACT * (1-CHIFRACT) +
     &     ((NEIGHBORES(LO,LO,HI,HI) - NEIGHBORES(LO,LO,HI,LO)) /
     &     THTWIDTH) * (1-PHIFRACT) * (1-PSIFRACT) * CHIFRACT +
     &     ((NEIGHBORES(LO,LO,LO,HI) - NEIGHBORES(LO,LO,LO,LO)) /
     &     THTWIDTH) * (1-PHIFRACT) * (1-PSIFRACT) * (1-CHIFRACT)
C
C if we're in print mode, handle it.
C
      IF (WHICH.EQ.'ANALYZE') THEN
           CALCRAMA(1,COUNT) = ECURRENT
           CALCRAMA(2,COUNT) = PHI
           CALCRAMA(3,COUNT) = PSI
           CALCRAMA(4,COUNT) = CHI
           CALCRAMA(5,COUNT) = THT
      ELSE
C
C don't update energy & forces
C if you're at a spot without an expectation value
C
           IF (OKPOINT) THEN
                ER = ER + ECURRENT
C
C don't update forces if you're in square well mode and
C within error
C
                IF (.NOT.((ECURRENT.LT.ERR).AND.
     &               (RAMAPOTENTIAL.EQ.SQUARE))) THEN
C
C calculate the forces & update the main deriv arrays
C
                     CALL CALCTORSIONDERIVS(I, J, K, L, DEDPHI)
                     CALL CALCTORSIONDERIVS(M, N, P, Q, DEDPSI)
                     CALL CALCTORSIONDERIVS(R, S, T, U, DEDCHI)
                     CALL CALCTORSIONDERIVS(D, F, G, H, DEDTHT)
                END IF
           END IF
      END IF
      RETURN
      END
C===============
      SUBROUTINE ERAMA3D (ER, COUNT, ATOMS, EXPECTEDRAMA,
     &     CALCRAMA, CURPHIS, CURPSIS, CURCHIS, WHICH)
C
C Calculates 3D torsion angle database energies
C
C Torsion angle database energies are of the form
C      E = k1*expectedRama
C where
C      k1 = force constant,
C
C WHICH is a flag that switches between energy & force and
C calculated chem shift (for violations) calcs
C
C by John Kuszewski May 1995 Mar 1997
C
C===============
      IMPLICIT NONE
C function declarations
      DOUBLE PRECISION TORPOSITION
      EXTERNAL TORPOSITION
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*), CURPHIS, CURPSIS, CURCHIS, COUNT
      DOUBLE PRECISION EXPECTEDRAMA(CURPHIS, CURPSIS, CURCHIS)
      DOUBLE PRECISION CALCRAMA(5,*), ER
      CHARACTER*7 WHICH
C for convenience
      INTEGER HI, LO
      PARAMETER (HI = 2, LO = 1)
C local variables
      INTEGER CLASS, I, J, K, L, M, N, P, Q, R, S, T, U,
     &     LOWPHI, HIPHI, LOWPSI, HIPSI, LOWCHI, HICHI,
     &     COUNT1, COUNT2, COUNT3
      DOUBLE PRECISION PHIWIDTH, PSIWIDTH, CHIWIDTH,
     &     K1, ERR,
     &     CORRECT1, CORRECT2, CORRECT3, CORRECT4, CORRECT5, CORRECT6,
     &     CORRECT7, CORRECT8, CORRECT9,
     &     PHI, PSI, CHI, CURPOS, NEXTPOS,
     &     PHIFRACT, PSIFRACT, CHIFRACT,
     &     ECURRENT, DEDPHI, DEDPSI, DEDCHI
      DOUBLE PRECISION NEIGHBORES(2,2,2)
      LOGICAL OKPOINT
C
C begin
C
      CLASS = ATOMS(1,COUNT)
C
      K1 = RAMAFORCES(CLASS) * RAMASCALE
      ERR = RAMAERRORS(CLASS)
      CORRECT1 = RAMAPHASEC(1, CLASS)
      CORRECT2 = RAMAPHASEC(2, CLASS)
      CORRECT3 = RAMAPHASEC(3, CLASS)
      CORRECT4 = RAMAPHASEC(4, CLASS)
      CORRECT5 = RAMAPHASEC(5, CLASS)
      CORRECT6 = RAMAPHASEC(6, CLASS)
      CORRECT7 = RAMAPHASEC(7, CLASS)
      CORRECT8 = RAMAPHASEC(8, CLASS)
      CORRECT9 = RAMAPHASEC(9, CLASS)
C
C calculate bin widths
C
      IF (CORRECT2.LT.(PI-CORRECT1)) THEN
           PHIWIDTH = (2*PI - CORRECT2) / CURPHIS
      ELSE IF (CORRECT2.GT.(PI+CORRECT1)) THEN
           PHIWIDTH = CORRECT2 / CURPHIS
      ELSE
           PHIWIDTH = (PI+CORRECT1) / CURPHIS
      END IF
C
      IF (CORRECT5.LT.(PI-CORRECT4)) THEN
           PSIWIDTH = (2*PI - CORRECT5) / CURPSIS
      ELSE IF (CORRECT5.GT.(PI+CORRECT4)) THEN
           PSIWIDTH = CORRECT5 / CURPSIS
      ELSE
           PSIWIDTH = (PI+CORRECT4) / CURPSIS
      END IF
C
      IF (CORRECT8.LT.(PI-CORRECT7)) THEN
           CHIWIDTH = (2*PI - CORRECT8) / CURCHIS
      ELSE IF (CORRECT8.GT.(PI+CORRECT7)) THEN
           CHIWIDTH = CORRECT8 / CURCHIS
      ELSE
           CHIWIDTH = (PI+CORRECT7) / CURCHIS
      END IF
C
C get coords of each atom
C
      I = ATOMS(2, COUNT)
      J = ATOMS(3, COUNT)
      K = ATOMS(4, COUNT)
      L = ATOMS(5, COUNT)
      M = ATOMS(6, COUNT)
      N = ATOMS(7, COUNT)
      P = ATOMS(8, COUNT)
      Q = ATOMS(9, COUNT)
      R = ATOMS(10, COUNT)
      S = ATOMS(11, COUNT)
      T = ATOMS(12, COUNT)
      U = ATOMS(13, COUNT)
C
C determine what the torsion angles are (not rounded!)
C
      PHI = TORPOSITION(I, J, K, L, CORRECT1, CORRECT2, CORRECT3)
      PSI = TORPOSITION(M, N, P, Q, CORRECT4, CORRECT5, CORRECT6)
      CHI = TORPOSITION(R, S, T, U, CORRECT7, CORRECT8, CORRECT9)
C
C get the neighboring points
C and the fractional coordinates within each square
C
C Note that the center of bin i is assumed to be given by
C (i-1) binWidth, where binWidth is the number of radians/bin.
C This is so that bin 1 will be centered about 0 deg.
C
      DO COUNT1 = 1, CURPHIS
           CURPOS = (COUNT1 - 1)*PHIWIDTH
           NEXTPOS = (COUNT1)*PHIWIDTH
           IF ((CURPOS.LE.PHI).AND.(NEXTPOS.GT.PHI)) THEN
                LOWPHI = COUNT1
                HIPHI = COUNT1+1
                PHIFRACT = (PHI - CURPOS) / PHIWIDTH
           END IF
      END DO
      IF (HIPHI.GT.CURPHIS) HIPHI = HIPHI - CURPHIS
C
      DO COUNT1 = 1, CURPSIS
           CURPOS = (COUNT1 - 1)*PSIWIDTH
           NEXTPOS = (COUNT1)*PSIWIDTH
           IF ((CURPOS.LE.PSI).AND.(NEXTPOS.GT.PSI)) THEN
                LOWPSI = COUNT1
                HIPSI = COUNT1+1
                PSIFRACT = (PSI - CURPOS) / PSIWIDTH
           END IF
      END DO
      IF (HIPSI.GT.CURPSIS) HIPSI = HIPSI - CURPSIS
C
      DO COUNT1 = 1, CURCHIS
           CURPOS = (COUNT1 - 1)*CHIWIDTH
           NEXTPOS = (COUNT1)*CHIWIDTH
           IF ((CURPOS.LE.CHI).AND.(NEXTPOS.GT.CHI)) THEN
                LOWCHI = COUNT1
                HICHI = COUNT1+1
                CHIFRACT = (CHI - CURPOS) / CHIWIDTH
           END IF
      END DO
      IF (HICHI.GT.CURCHIS) HICHI = HICHI - CURCHIS
C
C look up energies of neighbors
C
      NEIGHBORES(HI,HI,HI) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,HICHI)
      NEIGHBORES(LO,HI,HI) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,HICHI)
      NEIGHBORES(HI,LO,HI) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,HICHI)
      NEIGHBORES(LO,LO,HI) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,HICHI)
      NEIGHBORES(HI,HI,LO) =
     &     EXPECTEDRAMA(HIPHI,HIPSI,LOWCHI)
      NEIGHBORES(LO,HI,LO) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI,LOWCHI)
      NEIGHBORES(HI,LO,LO) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI,LOWCHI)
      NEIGHBORES(LO,LO,LO) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI,LOWCHI)
C
C are they all defined?
C
      OKPOINT = .TRUE.
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                DO COUNT3 = LO,HI
                     IF (NEIGHBORES(COUNT1,COUNT2,COUNT3)
     &                    .GT.NOEXPTEST) THEN
                          OKPOINT = .FALSE.
                     END IF
                END DO
           END DO
      END DO
C
C multiply neighbor energies by force constant
C
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                DO COUNT3 = LO,HI
                     NEIGHBORES(COUNT1,COUNT2,COUNT3) =
     &                    NEIGHBORES(COUNT1,COUNT2,COUNT3) * K1
                END DO
           END DO
      END DO
C
C interpolate energy
C
      ECURRENT =
     &     NEIGHBORES(HI,HI,HI) *
     &      PHIFRACT * PSIFRACT * CHIFRACT +
     &     NEIGHBORES(LO,HI,HI) *
     &      (1-PHIFRACT) * PSIFRACT * CHIFRACT +
     &     NEIGHBORES(HI,LO,HI) *
     &      PHIFRACT * (1-PSIFRACT) * CHIFRACT +
     &     NEIGHBORES(LO,LO,HI) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * CHIFRACT +
     &     NEIGHBORES(HI,HI,LO) *
     &      PHIFRACT * PSIFRACT * (1-CHIFRACT) +
     &     NEIGHBORES(LO,HI,LO) *
     &      (1-PHIFRACT) * PSIFRACT * (1-CHIFRACT) +
     &     NEIGHBORES(HI,LO,LO) *
     &      PHIFRACT * (1-PSIFRACT) * (1-CHIFRACT) +
     &     NEIGHBORES(LO,LO,LO) *
     &      (1-PHIFRACT) * (1-PSIFRACT) * (1-CHIFRACT)
C
C interpolate derivs
C
      DEDPHI =
     &     ((NEIGHBORES(HI,HI,HI) - NEIGHBORES(LO,HI,HI)) /
     &     PHIWIDTH) * PSIFRACT * CHIFRACT +
     &     ((NEIGHBORES(HI,LO,HI) - NEIGHBORES(LO,LO,HI)) /
     &     PHIWIDTH) * (1-PSIFRACT) * CHIFRACT +
     &     ((NEIGHBORES(HI,HI,LO) - NEIGHBORES(LO,HI,LO)) /
     &     PHIWIDTH) * PSIFRACT * (1-CHIFRACT) +
     &     ((NEIGHBORES(HI,LO,LO) - NEIGHBORES(LO,LO,LO)) /
     &     PHIWIDTH) * (1-PSIFRACT) * (1-CHIFRACT)
C
      DEDPSI =
     &     ((NEIGHBORES(HI,HI,HI) - NEIGHBORES(HI,LO,HI)) /
     &     PSIWIDTH) * PHIFRACT * CHIFRACT +
     &     ((NEIGHBORES(LO,HI,HI) - NEIGHBORES(LO,LO,HI)) /
     &     PSIWIDTH) * (1-PHIFRACT) * CHIFRACT +
     &     ((NEIGHBORES(HI,HI,LO) - NEIGHBORES(HI,LO,LO)) /
     &     PSIWIDTH) * PHIFRACT * (1-CHIFRACT) +
     &     ((NEIGHBORES(LO,HI,LO) - NEIGHBORES(LO,LO,LO)) /
     &     PSIWIDTH) * (1-PHIFRACT) * (1-CHIFRACT)
C
      DEDCHI =
     &     ((NEIGHBORES(HI,HI,HI) - NEIGHBORES(HI,HI,LO)) /
     &     CHIWIDTH) * PHIFRACT * PSIFRACT +
     &     ((NEIGHBORES(LO,HI,HI) - NEIGHBORES(LO,HI,LO)) /
     &     CHIWIDTH) * (1-PHIFRACT) * PSIFRACT +
     &     ((NEIGHBORES(HI,LO,HI) - NEIGHBORES(HI,LO,LO)) /
     &     CHIWIDTH) * PHIFRACT * (1-PSIFRACT) +
     &     ((NEIGHBORES(LO,LO,HI) - NEIGHBORES(LO,LO,LO)) /
     &     CHIWIDTH) * (1-PHIFRACT) * (1-PSIFRACT)
C
C if we're in print mode, handle it.
C
      IF (WHICH.EQ.'ANALYZE') THEN
           CALCRAMA(1,COUNT) = ECURRENT
           CALCRAMA(2,COUNT) = PHI
           CALCRAMA(3,COUNT) = PSI
           CALCRAMA(4,COUNT) = CHI
           CALCRAMA(5,COUNT) = ZERO
      ELSE
C
C don't update energy & forces
C if you're at a spot without an expectation value
C
           IF (OKPOINT) THEN
                ER = ER + ECURRENT
C
C don't update forces if you're in square well mode and
C within error
C
                IF (.NOT.((ECURRENT.LT.ERR).AND.
     &               (RAMAPOTENTIAL.EQ.SQUARE))) THEN
C
C calculate the forces & update the main deriv arrays
C
                     CALL CALCTORSIONDERIVS(I, J, K, L, DEDPHI)
                     CALL CALCTORSIONDERIVS(M, N, P, Q, DEDPSI)
                     CALL CALCTORSIONDERIVS(R, S, T, U, DEDCHI)
                END IF
           END IF
      END IF
      RETURN
      END
C===============
      SUBROUTINE ERAMA2D (ER, COUNT, ATOMS, EXPECTEDRAMA,
     &     CALCRAMA, CURPHIS, CURPSIS, WHICH)
C
C Calculates 2D torsion angle database energies
C
C Torsion angle database energies are of the form
C      E = k1*expectedRama
C where
C      k1 = force constant,
C
C WHICH is a flag that switches between energy & force and
C calculated chem shift (for violations) calcs
C
C by John Kuszewski May 1995 Mar 1997
C
C===============
      IMPLICIT NONE
C function declarations
      DOUBLE PRECISION TORPOSITION
      EXTERNAL TORPOSITION
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*), CURPHIS, CURPSIS, COUNT
      DOUBLE PRECISION EXPECTEDRAMA(CURPHIS, CURPSIS)
      DOUBLE PRECISION CALCRAMA(5,*), ER
      CHARACTER*7 WHICH
C for convenience
      INTEGER HI, LO
      PARAMETER (HI = 2, LO = 1)
C local variables
      INTEGER CLASS, I, J, K, L, M, N, P, Q,
     &     LOWPHI, HIPHI, LOWPSI, HIPSI,
     &     COUNT1, COUNT2
      DOUBLE PRECISION PHIWIDTH, PSIWIDTH,
     &     K1, ERR,
     &     CORRECT1, CORRECT2, CORRECT3, CORRECT4, CORRECT5, CORRECT6,
     &     PHI, PSI, CURPOS, NEXTPOS,
     &     PHIFRACT, PSIFRACT,
     &     ECURRENT, DEDPHI, DEDPSI
      DOUBLE PRECISION NEIGHBORES(2,2)
      LOGICAL OKPOINT
C
C begin
C
      CLASS = ATOMS(1,COUNT)
C
      K1 = RAMAFORCES(CLASS) * RAMASCALE
      ERR = RAMAERRORS(CLASS)
      CORRECT1 = RAMAPHASEC(1, CLASS)
      CORRECT2 = RAMAPHASEC(2, CLASS)
      CORRECT3 = RAMAPHASEC(3, CLASS)
      CORRECT4 = RAMAPHASEC(4, CLASS)
      CORRECT5 = RAMAPHASEC(5, CLASS)
      CORRECT6 = RAMAPHASEC(6, CLASS)
C
C calculate bin widths
C
      IF (CORRECT2.LT.(PI-CORRECT1)) THEN
           PHIWIDTH = (2*PI - CORRECT2) / CURPHIS
      ELSE IF (CORRECT2.GT.(PI+CORRECT1)) THEN
           PHIWIDTH = CORRECT2 / CURPHIS
      ELSE
           PHIWIDTH = (PI+CORRECT1) / CURPHIS
      END IF
C
      IF (CORRECT5.LT.(PI-CORRECT4)) THEN
           PSIWIDTH = (2*PI - CORRECT5) / CURPSIS
      ELSE IF (CORRECT5.GT.(PI+CORRECT4)) THEN
           PSIWIDTH = CORRECT5 / CURPSIS
      ELSE
           PSIWIDTH = (PI+CORRECT4) / CURPSIS
      END IF
C
C get coords of each atom
C
      I = ATOMS(2, COUNT)
      J = ATOMS(3, COUNT)
      K = ATOMS(4, COUNT)
      L = ATOMS(5, COUNT)
      M = ATOMS(6, COUNT)
      N = ATOMS(7, COUNT)
      P = ATOMS(8, COUNT)
      Q = ATOMS(9, COUNT)
C
C determine what the torsion angles are (not rounded!)
C
      PHI = TORPOSITION(I, J, K, L, CORRECT1, CORRECT2, CORRECT3)
      PSI = TORPOSITION(M, N, P, Q, CORRECT4, CORRECT5, CORRECT6)
C
C get the neighboring points
C and the fractional coordinates within each square
C
C Note that the center of bin i is assumed to be given by
C (i-1) binWidth, where binWidth is the number of radians/bin.
C This is so that bin 1 will be centered about 0 deg.
C
      DO COUNT1 = 1, CURPHIS
           CURPOS = (COUNT1 - 1)*PHIWIDTH
           NEXTPOS = (COUNT1)*PHIWIDTH
           IF ((CURPOS.LE.PHI).AND.(NEXTPOS.GT.PHI)) THEN
                LOWPHI = COUNT1
                HIPHI = COUNT1+1
                PHIFRACT = (PHI - CURPOS) / PHIWIDTH
           END IF
      END DO
      IF (HIPHI.GT.CURPHIS) HIPHI = HIPHI - CURPHIS
C
      DO COUNT1 = 1, CURPSIS
           CURPOS = (COUNT1 - 1)*PSIWIDTH
           NEXTPOS = (COUNT1)*PSIWIDTH
           IF ((CURPOS.LE.PSI).AND.(NEXTPOS.GT.PSI)) THEN
                LOWPSI = COUNT1
                HIPSI = COUNT1+1
                PSIFRACT = (PSI - CURPOS) / PSIWIDTH
           END IF
      END DO
      IF (HIPSI.GT.CURPSIS) HIPSI = HIPSI - CURPSIS
C
C look up energies of neighbors
C
      NEIGHBORES(HI,HI) =
     &     EXPECTEDRAMA(HIPHI,HIPSI)
      NEIGHBORES(LO,HI) =
     &     EXPECTEDRAMA(LOWPHI,HIPSI)
      NEIGHBORES(HI,LO) =
     &     EXPECTEDRAMA(HIPHI,LOWPSI)
      NEIGHBORES(LO,LO) =
     &     EXPECTEDRAMA(LOWPHI,LOWPSI)
C
C are they all defined?
C
      OKPOINT = .TRUE.
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                IF (NEIGHBORES(COUNT1,COUNT2)
     &               .GT.NOEXPTEST) THEN
                     OKPOINT = .FALSE.
                END IF
           END DO
      END DO
C
C multiply neighbor energies by force constant
C
      DO COUNT1 = LO,HI
           DO COUNT2 = LO,HI
                NEIGHBORES(COUNT1,COUNT2) =
     &               NEIGHBORES(COUNT1,COUNT2) * K1
           END DO
      END DO
C
C interpolate energy
C
      ECURRENT =
     &     NEIGHBORES(HI,HI) *
     &      PHIFRACT * PSIFRACT +
     &     NEIGHBORES(LO,HI) *
     &      (1-PHIFRACT) * PSIFRACT +
     &     NEIGHBORES(HI,LO) *
     &      PHIFRACT * (1-PSIFRACT) +
     &     NEIGHBORES(LO,LO) *
     &      (1-PHIFRACT) * (1-PSIFRACT)
C
C interpolate derivs
C
      DEDPHI =
     &     ((NEIGHBORES(HI,HI) - NEIGHBORES(LO,HI)) /
     &     PHIWIDTH) * PSIFRACT +
     &     ((NEIGHBORES(HI,LO) - NEIGHBORES(LO,LO)) /
     &     PHIWIDTH) * (1-PSIFRACT)
C
      DEDPSI =
     &     ((NEIGHBORES(HI,HI) - NEIGHBORES(HI,LO)) /
     &     PSIWIDTH) * PHIFRACT +
     &     ((NEIGHBORES(LO,HI) - NEIGHBORES(LO,LO)) /
     &     PSIWIDTH) * (1-PHIFRACT)
C
C if we're in print mode, handle it.
C
      IF (WHICH.EQ.'ANALYZE') THEN
           CALCRAMA(1,COUNT) = ECURRENT
           CALCRAMA(2,COUNT) = PHI
           CALCRAMA(3,COUNT) = PSI
           CALCRAMA(4,COUNT) = ZERO
           CALCRAMA(5,COUNT) = ZERO
      ELSE
C
C don't update energy & forces
C if you're at a spot without an expectation value
C
           IF (OKPOINT) THEN
                ER = ER + ECURRENT
C
C don't update forces if you're in square well mode and
C within error
C
                IF (.NOT.((ECURRENT.LT.ERR).AND.
     &               (RAMAPOTENTIAL.EQ.SQUARE))) THEN
C
C calculate the forces & update the main deriv arrays
C
                     CALL CALCTORSIONDERIVS(I, J, K, L, DEDPHI)
                     CALL CALCTORSIONDERIVS(M, N, P, Q, DEDPSI)
                END IF
           END IF
      END IF
      RETURN
      END
C===============
      SUBROUTINE ERAMA1D (ER, COUNT, ATOMS, EXPECTEDRAMA,
     &     CALCRAMA, CURPHIS, WHICH)
C
C Calculates 1D torsion angle database energies
C
C Torsion angle database energies are of the form
C      E = k1*expectedRama
C where
C      k1 = force constant,
C
C WHICH is a flag that switches between energy & force and
C calculated chem shift (for violations) calcs
C
C NOTE:  change ramapotential to allow separate potentials for
C different classes
C
C by John Kuszewski May 1995 Mar 1997
C
C===============
      IMPLICIT NONE
C function declarations
      DOUBLE PRECISION TORPOSITION
      EXTERNAL TORPOSITION
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*), CURPHIS, COUNT
      DOUBLE PRECISION EXPECTEDRAMA(CURPHIS)
      DOUBLE PRECISION CALCRAMA(5,*), ER
      CHARACTER*7 WHICH
C for convenience
      INTEGER HI, LO
      PARAMETER (HI = 2, LO = 1)
C local variables
      INTEGER CLASS, I, J, K, L,
     &     LOWPHI, HIPHI,
     &     COUNT1
      DOUBLE PRECISION PHIWIDTH,
     &     K1, ERR,
     &     CORRECT1, CORRECT2, CORRECT3,
     &     PHI, CURPOS, NEXTPOS,
     &     PHIFRACT,
     &     ECURRENT, DEDPHI
      DOUBLE PRECISION NEIGHBORES(2)
      LOGICAL OKPOINT
C
C begin
C
      CLASS = ATOMS(1,COUNT)
C
      K1 = RAMAFORCES(CLASS) * RAMASCALE
      ERR = RAMAERRORS(CLASS)
      CORRECT1 = RAMAPHASEC(1, CLASS)
      CORRECT2 = RAMAPHASEC(2, CLASS)
      CORRECT3 = RAMAPHASEC(3, CLASS)
C
C calculate bin widths
C
      IF (CORRECT2.LT.(PI-CORRECT1)) THEN
           PHIWIDTH = (2*PI - CORRECT2) / CURPHIS
      ELSE IF (CORRECT2.GT.(PI+CORRECT1)) THEN
           PHIWIDTH = CORRECT2 / CURPHIS
      ELSE
           PHIWIDTH = (PI+CORRECT1) / CURPHIS
      END IF
C
C get coords of each atom
C
      I = ATOMS(2, COUNT)
      J = ATOMS(3, COUNT)
      K = ATOMS(4, COUNT)
      L = ATOMS(5, COUNT)
C
C determine what the torsion angles are (not rounded!)
C
      PHI = TORPOSITION(I, J, K, L, CORRECT1, CORRECT2, CORRECT3)
C
C get the neighboring points
C and the fractional coordinates within each square
C
C Note that the center of bin i is assumed to be given by
C (i-1) binWidth, where binWidth is the number of radians/bin.
C This is so that bin 1 will be centered about 0 deg.
C
      DO COUNT1 = 1, CURPHIS
           CURPOS = (COUNT1 - 1)*PHIWIDTH
           NEXTPOS = (COUNT1)*PHIWIDTH
           IF ((CURPOS.LE.PHI).AND.(NEXTPOS.GT.PHI)) THEN
                LOWPHI = COUNT1
                HIPHI = COUNT1+1
                PHIFRACT = (PHI - CURPOS) / PHIWIDTH
           END IF
      END DO
      IF (HIPHI.GT.CURPHIS) HIPHI = HIPHI - CURPHIS
C
C look up energies of neighbors
C
      NEIGHBORES(HI) =
     &     EXPECTEDRAMA(HIPHI)
      NEIGHBORES(LO) =
     &     EXPECTEDRAMA(LOWPHI)
C
C are they all defined?
C
      OKPOINT = .TRUE.
      DO COUNT1 = LO,HI
           IF (NEIGHBORES(COUNT1)
     &          .GT.NOEXPTEST) THEN
                OKPOINT = .FALSE.
           END IF
      END DO
C
C multiply neighbor energies by force constant
C
      DO COUNT1 = LO,HI
           NEIGHBORES(COUNT1) =
     &          NEIGHBORES(COUNT1) * K1
      END DO
C
C interpolate energy
C
      ECURRENT =
     &     NEIGHBORES(HI) *
     &      PHIFRACT +
     &     NEIGHBORES(LO) *
     &      (1-PHIFRACT)
C
C interpolate derivs
C
      DEDPHI =
     &     ((NEIGHBORES(HI) - NEIGHBORES(LO)) /
     &     PHIWIDTH)
C
C if we're in print mode, handle it.
C
      IF (WHICH.EQ.'ANALYZE') THEN
           CALCRAMA(1,COUNT) = ECURRENT
           CALCRAMA(2,COUNT) = PHI
           CALCRAMA(3,COUNT) = ZERO
           CALCRAMA(4,COUNT) = ZERO
           CALCRAMA(5,COUNT) = ZERO
      ELSE
C
C don't update energy & forces
C if you're at a spot without an expectation value
C
           IF (OKPOINT) THEN
                ER = ER + ECURRENT
C
C don't update forces if you're in square well mode and
C within error
C
                IF (.NOT.((ECURRENT.LT.ERR).AND.
     &               (RAMAPOTENTIAL.EQ.SQUARE))) THEN
C
C calculate the forces & update the main deriv arrays
C
                     CALL CALCTORSIONDERIVS(I, J, K, L, DEDPHI)
                END IF
           END IF
      END IF
      RETURN
      END
C===============
      DOUBLE PRECISION FUNCTION TORPOSITION (
     &     ATOMI, ATOMJ, ATOMK, ATOML,
     &     CORRECT1, CORRECT2, CORRECT3)
C
C Given four atoms' coordinates, calculates the value of the
C torsion angle connecting them.
C
C Implemented in order to clean up this code.
C
C by John Kuszewski Oct 1996  Mar 1997
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'consta.inc'
C i/o
      INTEGER ATOMI, ATOMJ, ATOMK, ATOML
      DOUBLE PRECISION CORRECT1, CORRECT2, CORRECT3
C local vbls
      DOUBLE PRECISION XI, XJ, XK, XL, YI, YJ, YK, YL,
     &     ZI, ZJ, ZK, ZL, XIJ, XJK, XKL, YIJ, YJK, YKL,
     &     ZIJ, ZJK, ZKL, AX, AY, AZ, BX, BY, BZ, CX, CY, CZ,
     &     RAR, RBR, RCR, CP, SP, PHI, CORRECTPHI
C begin
C
C get coords of each atom
C
      XI = X(ATOMI)
      XJ = X(ATOMJ)
      XK = X(ATOMK)
      XL = X(ATOML)
C
      YI = Y(ATOMI)
      YJ = Y(ATOMJ)
      YK = Y(ATOMK)
      YL = Y(ATOML)
C
      ZI = Z(ATOMI)
      ZJ = Z(ATOMJ)
      ZK = Z(ATOMK)
      ZL = Z(ATOML)
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
      XIJ = XI - XJ
      XJK = XJ - XK
      XKL = XK - XL
C
      YIJ = YI - YJ
      YJK = YJ - YK
      YKL = YK - YL
C
      ZIJ = ZI - ZJ
      ZJK = ZJ - ZK
      ZKL = ZK - ZL
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
C calculate the norm of A, B, & C & set to MCONST if it's too small
C
      RAR = ONE/SQRT(MAX(MCONST, AX**2+AY**2+AZ**2))
      RBR = ONE/SQRT(MAX(MCONST, BX**2+BY**2+BZ**2))
      RCR = ONE/SQRT(MAX(MCONST, CX**2+CY**2+CZ**2))
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
C calculate cos(phi) & sin(phi)
C
      CP = AX*BX+AY*BY+AZ*BZ
      SP = CX*BX+CY*BY+CZ*BZ
C
C calculate phi (make sure cos is within bounds and get sign from sin)
C and keep it it radians
C
      PHI = -SIGN(ACOS(MIN(ONE,MAX(-ONE,CP))),SP)
C
C figure out which bin you should go into, given that the
C expectation values are in the range 1..PHISTEPS
C
C
C added the possibility of different phase corrections
C for different classes of PROCHECK data (eg., phi/psi
C vs. chi1/chi2 plots) 9/1/95
C
      IF (PHI.LT.CORRECT1) THEN
           CORRECTPHI = PHI + CORRECT2
      ELSE
           CORRECTPHI = PHI
      END IF
      CORRECTPHI = CORRECTPHI + CORRECT3
C
C enforce bounds
C
      IF (CORRECTPHI.GT.(TWO*PI)) THEN
           CORRECTPHI = CORRECTPHI - (TWO*PI)
      ELSE IF (CORRECTPHI.LT.ZERO) THEN
           CORRECTPHI = CORRECTPHI + (TWO*PI)
      END IF
C
C convert to degrees
C
C      CALL RD2DEG(CORRECTPHI,CORRECTPHI)
C
C return it
C
      TORPOSITION = CORRECTPHI
      RETURN
      END
C===============
      SUBROUTINE CALCTORSIONDERIVS(ATOMI, ATOMJ, ATOMK, ATOML, DF)
C
C Given four atoms' coordinates and d Etorsion / d phi,
C generate the forces on each atom and update the main
C derivative arrays.
C
C by John Kuszewski Oct 1996
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'rama.inc'
C i/o
      INTEGER ATOMI, ATOMJ, ATOMK, ATOML
      DOUBLE PRECISION DF
C local vbls
      DOUBLE PRECISION XI, XJ, XK, XL, YI, YJ, YK, YL,
     &     ZI, ZJ, ZK, ZL, XIJ, XJK, XKL, YIJ, YJK, YKL,
     &     ZIJ, ZJK, ZKL, AX, AY, AZ, BX, BY, BZ, CX, CY, CZ,
     &     RAR, RBR, RCR, CP, SP
      DOUBLE PRECISION SWITCH, RECSP, RECCP, DCPAX, DCPAY, DCPAZ,
     &     DCPBX, DCPBY, DCPBZ, DSPCX, DSPCY, DSPCZ, DSPBX, DSPBY,
     &     DSPBZ, DPRIJX, DPRIJY, DPRIJZ, DPRKLX, DPRKLY, DPRKLZ,
     &     DPRJKX, DPRJKY, DPRJKZ
C begin
C
C get coords of each atom
C
      XI = X(ATOMI)
      XJ = X(ATOMJ)
      XK = X(ATOMK)
      XL = X(ATOML)
C
      YI = Y(ATOMI)
      YJ = Y(ATOMJ)
      YK = Y(ATOMK)
      YL = Y(ATOML)
C
      ZI = Z(ATOMI)
      ZJ = Z(ATOMJ)
      ZK = Z(ATOMK)
      ZL = Z(ATOML)
C
C now calculate diffs RIJ=RI-RJ, RJK=RJ-RK, RKL=RK-RL
C
      XIJ = XI - XJ
      XJK = XJ - XK
      XKL = XK - XL
C
      YIJ = YI - YJ
      YJK = YJ - YK
      YKL = YK - YL
C
      ZIJ = ZI - ZJ
      ZJK = ZJ - ZK
      ZKL = ZK - ZL
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
C calculate the norm of A, B, & C & set to MCONST if it's too small
C
      RAR = ONE/SQRT(MAX(MCONST, AX**2+AY**2+AZ**2))
      RBR = ONE/SQRT(MAX(MCONST, BX**2+BY**2+BZ**2))
      RCR = ONE/SQRT(MAX(MCONST, CX**2+CY**2+CZ**2))
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
C calculate cos(phi) & sin(phi)
C
      CP = AX*BX+AY*BY+AZ*BZ
      SP = CX*BX+CY*BY+CZ*BZ
C
C compute heavyside function
C
      SWITCH=-MIN(1,INT(ABS(SP)-EPS+ONE))
C
C compute switched and truncated reciprocals of CP and SP multiplied
C by the derivative of the function DF
C
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
C
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
     &     +(TWO*XJK*YIJ-XIJ*YJK)*DSPCY
     &     +(TWO*XJK*ZIJ-XIJ*ZJK)*DSPCZ
     &     +DSPBZ*YKL-DSPBY*ZKL)
      DPRJKY=
     &     RECSP*(DCPAZ*XIJ-DCPAX*ZIJ+DCPBX*ZKL-DCPBZ*XKL)+
     &     RECCP*(-(ZJK*ZIJ+XJK*XIJ)*DSPCY
     &     +(TWO*YJK*ZIJ-YIJ*ZJK)*DSPCZ
     &     +(TWO*YJK*XIJ-YIJ*XJK)*DSPCX
     &     +DSPBX*ZKL-DSPBZ*XKL)
      DPRJKZ=
     &     RECSP*(DCPAX*YIJ-DCPAY*XIJ+DCPBY*XKL-DCPBX*YKL)+
     &     RECCP*(-(XJK*XIJ+YJK*YIJ)*DSPCZ
     &     +(TWO*ZJK*XIJ-ZIJ*XJK)*DSPCX
     &     +(TWO*ZJK*YIJ-ZIJ*YJK)*DSPCY
     &     +DSPBY*XKL-DSPBX*YKL)
C
C now update forces
C
      DX(ATOMI)=DX(ATOMI)+DPRIJX
      DY(ATOMI)=DY(ATOMI)+DPRIJY
      DZ(ATOMI)=DZ(ATOMI)+DPRIJZ
      DX(ATOMJ)=DX(ATOMJ)+DPRJKX-DPRIJX
      DY(ATOMJ)=DY(ATOMJ)+DPRJKY-DPRIJY
      DZ(ATOMJ)=DZ(ATOMJ)+DPRJKZ-DPRIJZ
      DX(ATOMK)=DX(ATOMK)+DPRKLX-DPRJKX
      DY(ATOMK)=DY(ATOMK)+DPRKLY-DPRJKY
      DZ(ATOMK)=DZ(ATOMK)+DPRKLZ-DPRJKZ
      DX(ATOML)=DX(ATOML)       -DPRKLX
      DY(ATOML)=DY(ATOML)       -DPRKLY
      DZ(ATOML)=DZ(ATOML)       -DPRKLZ
C
C the end
C
      RETURN
      END
C================
      SUBROUTINE READRAMA
C
C reads in PROCHECK Ramachandran information
C
C by John Kuszewski May 1995
C================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
C local variables
      INTEGER COUNT, SPTR, OLDCLASS, OLDMAXRAMAS,
     &     THETYPE, CURPSIS, CURPHIS,
     &     CURCHIS, CURTHTS, CLASINDEX
      INTEGER CLEN
      DOUBLE PRECISION K1, CUTOFF
      CHARACTER*4 THENAME
      CHARACTER*50 CLASNAME
C begin
C
C this is used by READRAMA2 to hold the selection
C
      SPTR=ALLHP(INTEG4(NATOM))
C
      IF (ALLOCRAMA) THEN
      CALL RAMADEFAULTS
      CALL ALLOCRAMAS(0, MAXRAMAS)
      ALLOCRAMA=.FALSE.
      END IF
C
C now read input
C
      CALL PUSEND('RAMA>')
      DO WHILE (.NOT.DONE)
           CALL NEXTWD('RAMA>')
           CALL MISCOM('RAMA>',USED)
           IF (.NOT.USED) THEN
C
           IF (WD(1:4).EQ.'HELP') THEN
C
              CALL CNSHELP('cns-conformation')
C
C Get class name.  Determine if it's an already-defined class.
C Insert a new class if it's not.
C
           ELSE IF (WD(1:4).EQ.'CLAS') THEN
                OLDCLASS = CURRAMACLASS
                CALL NEXTWD('class name =')
                CLASNAME = WD(1:50)
                RAMAMODE = NEW
                DO COUNT = 1, NRAMACLASSES
                     IF (RAMACLASSNAMES(COUNT).EQ.CLASNAME) THEN
                          RAMAMODE = UPDATE
                          CURRAMACLASS = COUNT
                     END IF
                END DO
                IF (RAMAMODE.EQ.NEW) THEN
C
C make sure you can't add more than the maximum
C number of classes
C
                     IF (OLDCLASS.EQ.MAXRAMACLASSES) THEN
                         CALL DSPERR('RAMA','Too many classes.')
                         CALL DSPERR('RAMA',
     &                      'Increase MAXRAMACLASSES and recompile.')
                         CALL WRNDIE(-5, 'READRAMA',
     &                        'Too many torsion angle DB classes.')
                     END IF
                     NRAMACLASSES = NRAMACLASSES + 1
                     CURRAMACLASS = NRAMACLASSES
                     RAMACLASSNAMES(CURRAMACLASS) = CLASNAME
                END IF
C
C set force constant for current class,
C starting a default class if there isnt one defined
C
           ELSE IF (WD(1:4).EQ.'FORC') THEN
                CALL NEXTF('force constant =', K1)
                IF (CURRAMACLASS.EQ.0) THEN
                     NRAMACLASSES = 1
                     CURRAMACLASS = 1
                END IF
                CLEN=LEN(RAMACLASSNAMES(CURRAMACLASS))
                CALL TRIMM(RAMACLASSNAMES(CURRAMACLASS),CLEN)
                WRITE(PUNIT, '(3(A), F8.3)')
     &            'Setting force const for class ',
     &             RAMACLASSNAMES(CURRAMACLASS)(1:CLEN), ' to ', K1
                RAMAFORCES(CURRAMACLASS) = K1
C
C read in overall scale factor for force constants
C
           ELSE IF (WD(1:4).EQ.'SCAL') THEN
                CALL NEXTF('overall force const scale factor =', K1)
                RAMASCALE = K1
C
C debugging support
C
           ELSE IF (WD(1:4).EQ.'DEBU') THEN
                DO COUNT = 1, NRAMACLASSES
                  CLEN=LEN(RAMACLASSNAMES(COUNT))
                  CALL TRIMM(RAMACLASSNAMES(COUNT),CLEN)
                  WRITE(6,'(A,I5,A)')'CLASS NUM', COUNT,
     &                                RAMACLASSNAMES(COUNT)(1:CLEN)
                END DO
C
C set phase corrections for current class,
C starting a default class if there isnt one defined
C
           ELSE IF (WD(1:4).EQ.'PHAS') THEN
                IF (CURRAMACLASS.EQ.0) THEN
                     NRAMACLASSES = 1
                     CURRAMACLASS = 1
                END IF
                CLEN=LEN(RAMACLASSNAMES(CURRAMACLASS))
                CALL TRIMM(RAMACLASSNAMES(CURRAMACLASS),CLEN)
C
                IF (RAMACLASSTYPE(CURRAMACLASS).EQ.1) THEN
                     CALL NEXTF('phase correction 1 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(1, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 2 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(2, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 3 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(3, CURRAMACLASS) = K1
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting phase corrs for 1D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to (in radians)'
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(1, CURRAMACLASS),
     &                    RAMAPHASEC(2, CURRAMACLASS),
     &                    RAMAPHASEC(3, CURRAMACLASS)
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.2) THEN
                     CALL NEXTF('phase correction 1 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(1, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 2 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(2, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 3 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(3, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 4 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(4, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 5 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(5, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 6 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(6, CURRAMACLASS) = K1
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting phase corrs for 2D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to (in radians)'
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(1, CURRAMACLASS),
     &                    RAMAPHASEC(2, CURRAMACLASS),
     &                    RAMAPHASEC(3, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(4, CURRAMACLASS),
     &                    RAMAPHASEC(5, CURRAMACLASS),
     &                    RAMAPHASEC(6, CURRAMACLASS)
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.3) THEN
                     CALL NEXTF('phase correction 1 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(1, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 2 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(2, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 3 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(3, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 4 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(4, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 5 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(5, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 6 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(6, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 7 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(7, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 8 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(8, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 9 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(9, CURRAMACLASS) = K1
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting phase corrs for 3D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to (in radians)'
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(1, CURRAMACLASS),
     &                    RAMAPHASEC(2, CURRAMACLASS),
     &                    RAMAPHASEC(3, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(4, CURRAMACLASS),
     &                    RAMAPHASEC(5, CURRAMACLASS),
     &                    RAMAPHASEC(6, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(7, CURRAMACLASS),
     &                    RAMAPHASEC(8, CURRAMACLASS),
     &                    RAMAPHASEC(9, CURRAMACLASS)
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.4) THEN
                     CALL NEXTF('phase correction 1 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(1, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 2 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(2, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 3 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(3, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 4 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(4, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 5 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(5, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 6 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(6, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 7 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(7, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 8 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(8, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 9 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(9, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 10 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(10, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 11 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(11, CURRAMACLASS) = K1
                     CALL NEXTF('phase correction 12 (deg) =', K1)
                     CALL DG2RAD(K1, K1)
                     RAMAPHASEC(12, CURRAMACLASS) = K1
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting phase corrs for 4D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to (in radians)'
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(1, CURRAMACLASS),
     &                    RAMAPHASEC(2, CURRAMACLASS),
     &                    RAMAPHASEC(3, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(4, CURRAMACLASS),
     &                    RAMAPHASEC(5, CURRAMACLASS),
     &                    RAMAPHASEC(6, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(7, CURRAMACLASS),
     &                    RAMAPHASEC(8, CURRAMACLASS),
     &                    RAMAPHASEC(9, CURRAMACLASS)
                     WRITE(PUNIT, '(3(F8.3))')
     &                    RAMAPHASEC(10, CURRAMACLASS),
     &                    RAMAPHASEC(11, CURRAMACLASS),
     &                    RAMAPHASEC(12, CURRAMACLASS)
                END IF
C
C set size of current class's expectation grid,
C starting a default class if there isn't one defined
C
           ELSE IF (WD(1:4).EQ.'SIZE') THEN
                IF (CURRAMACLASS.EQ.0) THEN
                     NRAMACLASSES = 1
                     CURRAMACLASS = 1
                END IF
                CLEN=LEN(RAMACLASSNAMES(CURRAMACLASS))
                CALL TRIMM(RAMACLASSNAMES(CURRAMACLASS),CLEN)
                CALL NEXTA4('dimensionality =', THENAME)
                IF (THENAME.EQ.'ONED') THEN
                     THETYPE = 1
                ELSE IF (THENAME.EQ.'TWOD') THEN
                     THETYPE = 2
                ELSE IF (THENAME.EQ.'THRE') THEN
                     THETYPE = 3
                ELSE IF (THENAME.EQ.'FOUR') THEN
                     THETYPE = 4
                ELSE
                     CALL DSPERR('RAMA','Not understood. Assuming 1D')
                     THETYPE = 1
                END IF
                IF (THETYPE.EQ.1) THEN
                     CALL NEXTI('phi step =', CURPHIS)
                     CURPSIS = 1
                     CURCHIS = 1
                     CURTHTS = 1
                     WRITE(PUNIT, '(3(A), I6)')
     &                    'Setting size of one-D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to ',CURPHIS
                ELSE IF (THETYPE.EQ.2) THEN
                     CALL NEXTI('phi step =', CURPHIS)
                     CALL NEXTI('psi step =', CURPSIS)
                     CURCHIS = 1
                     CURTHTS = 1
                     WRITE(PUNIT, '(3(A), I6, A, I6)')
     &                    'Setting size of two-D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN),
     &                    ' to ', CURPHIS, ' by ', CURPSIS
                ELSE IF (THETYPE.EQ.3) THEN
                     CALL NEXTI('phi step =', CURPHIS)
                     CALL NEXTI('psi step =', CURPSIS)
                     CALL NEXTI('chi step =', CURCHIS)
                     CURTHTS = 1
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting size of three-D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN), ' to '
                     WRITE(PUNIT, '(I6, A, I6, A, I6)')
     &                    CURPHIS, ' by ', CURPSIS, ' by ', CURCHIS
                ELSE IF (THETYPE.EQ.4) THEN
                     CALL NEXTI('phi step =', CURPHIS)
                     CALL NEXTI('psi step =', CURPSIS)
                     CALL NEXTI('chi step =', CURCHIS)
                     CALL NEXTI('theta step =', CURTHTS)
                     WRITE(PUNIT, '(3(A))')
     &                    'Setting size of four-D class ',
     &                    RAMACLASSNAMES(CURRAMACLASS)(1:CLEN), ' to '
                     WRITE(PUNIT, '(I6, A, I6, A, I6, A, I6)')
     &                    CURPHIS, ' by ', CURPSIS, ' by ', CURCHIS,
     &                    ' by ', CURTHTS
                END IF
                CALL ALLOCEXPECTEDRAMA(CURRAMACLASS,
     &               PHISTEPS(CURRAMACLASS)*
     &               PSISTEPS(CURRAMACLASS)*
     &               CHISTEPS(CURRAMACLASS)*
     &               THTSTEPS(CURRAMACLASS),
     &               CURPHIS*CURPSIS*CURCHIS*CURTHTS)
                RAMACLASSTYPE(CURRAMACLASS) = THETYPE
                PHISTEPS(CURRAMACLASS) = CURPHIS
                PSISTEPS(CURRAMACLASS) = CURPSIS
                CHISTEPS(CURRAMACLASS) = CURCHIS
                THTSTEPS(CURRAMACLASS) = CURTHTS
C
C set error for current class,
C starting a default class if there isn't one defined
C
           ELSE IF (WD(1:4).EQ.'ERRO') THEN
                CALL NEXTF('error =', K1)
                IF (CURRAMACLASS.EQ.0) THEN
                     NRAMACLASSES = 1
                     CURRAMACLASS = 1
                END IF
                CLEN=LEN(RAMACLASSNAMES(CURRAMACLASS))
                CALL TRIMM(RAMACLASSNAMES(CURRAMACLASS),CLEN)
                WRITE(PUNIT, '(3(A), F8.3)')
     &            'Setting error for class ',
     &             RAMACLASSNAMES(CURRAMACLASS)(1:CLEN), ' to ', K1
                RAMAERRORS(CURRAMACLASS) = K1
C
C reset RAMAs database
C
           ELSE IF (WD(1:4).EQ.'RESE') THEN
                CALL RAMAHP
                CALL RAMAINIT
                CALL RAMADEFAULTS
                CALL ALLOCRAMAS(0, MAXRAMAS)
                ALLOCRAMA=.FALSE.
C
C set potential type
C
           ELSE IF (WD(1:4).EQ.'POTE') THEN
                CALL NEXTA4('potential type =', THENAME)
                IF (THENAME.EQ.'SQUA') THEN
                    WRITE(PUNIT, '(A)') 'using square well potential.'
                    RAMAPOTENTIAL = SQUARE
                ELSE IF (THENAME.EQ.'HARM') THEN
                    WRITE(PUNIT, '(A)') 'using harmonic potential.'
                    RAMAPOTENTIAL = HARMONIC
                ELSE
                    CALL DSPERR('RAMA',
     &                        'unknown potential. Using square well.')
                    RAMAPOTENTIAL = SQUARE
                END IF
C
C change number of assignment slots
C
           ELSE IF (WD(1:4).EQ.'NRES') THEN
                OLDMAXRAMAS = MAXRAMAS
                CALL NEXTI('number of slots =', MAXRAMAS)
                CALL ALLOCRAMAS(OLDMAXRAMAS, MAXRAMAS)
C
C read in an assignment
C
           ELSE IF (WD(1:4).EQ.'ASSI') THEN
C
C make sure you can't add more chemical RAMAs
C than you have slots for
C
                IF (NRAMAS.EQ.MAXRAMAS) THEN
                     CALL DSPERR('CRAMAS','Too many assignments.')
                     CALL WRNDIE(-5, 'READRAMA',
     &                    'Too many assignments.')
                END IF
C
C if there isn't a class specified,
C start a default class
C
                IF (CURRAMACLASS.EQ.0) THEN
                     NRAMACLASSES = 1
                     CURRAMACLASS = 1
                END IF
                CALL READRAMAASSIGN(HEAP(RAMAATOMPTR),
     &               CURRAMACLASS, HEAP(SPTR))
C
C sort the assignment table
C
           ELSE IF (WD(1:4).EQ.'SORT') THEN
                CALL SORTRAMAS (HEAP(RAMAATOMPTR))
C
C read in compressed expectation values
C
           ELSE IF (WD(1:4).EQ.'COMP') THEN
                IF (RAMACLASSTYPE(CURRAMACLASS).EQ.1) THEN
                     CALL READCOMPRESSEDRAMA1D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS))
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.2) THEN
                     CALL READCOMPRESSEDRAMA2D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS))
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.3) THEN
                     CALL READCOMPRESSEDRAMA3D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS),
     &                    CHISTEPS(CURRAMACLASS))
                ELSE
                     CALL READCOMPRESSEDRAMA4D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS),
     &                    CHISTEPS(CURRAMACLASS),
     &                    THTSTEPS(CURRAMACLASS))
                END IF
C
C read in expectation values
C
           ELSE IF (WD(1:4).EQ.'EXPE') THEN
                IF (RAMACLASSTYPE(CURRAMACLASS).EQ.1) THEN
                     CALL READEXPECTEDRAMA1D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)))
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.2) THEN
                     CALL READEXPECTEDRAMA2D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS))
                ELSE IF (RAMACLASSTYPE(CURRAMACLASS).EQ.3) THEN
                     CALL READEXPECTEDRAMA3D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS),
     &                    CHISTEPS(CURRAMACLASS))
                ELSE
                     CALL READEXPECTEDRAMA4D(CURRAMACLASS,
     &                    HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &                    PHISTEPS(CURRAMACLASS),
     &                    PSISTEPS(CURRAMACLASS),
     &                    CHISTEPS(CURRAMACLASS),
     &                    THTSTEPS(CURRAMACLASS))
                END IF
C
C zero out the expectation value arrays
C
           ELSE IF (WD(1:4).EQ.'ZERO') THEN
                CALL DEFEXPECTEDRAMA(HEAP(EXPRAMAPTRS(CURRAMACLASS)),
     &               PHISTEPS(CURRAMACLASS)*PSISTEPS(CURRAMACLASS)*
     &               CHISTEPS(CURRAMACLASS)*THTSTEPS(CURRAMACLASS))
C
C print violations
C
           ELSE IF (WD(1:4).EQ.'PRIN') THEN
                CALL NEXTWD('PRINt>')
                IF (WD(1:4).NE.'THRE') THEN
                     CALL DSPERR('RAMA',
     &                           'print expects THREshold parameter.')
                ELSE
                     CALL NEXTF('THREshold =', CUTOFF)
                     IF (CUTOFF.LT.ZERO) THEN
                          CALL DSPERR('RAMA',
     &                         'cutoff must be positive.')
                          CUTOFF = ABS(CUTOFF)
                     END IF
                     CALL NEXTWD('ALL or CLASs>')
                     IF (WD(1:3).EQ.'ALL') THEN
                          CALL PRINTRAMAS(CUTOFF, HEAP(CALCRAMAPTR),
     &                         HEAP(RAMAATOMPTR), 0)
                     ELSE IF (WD(1:4).EQ.'CLAS') THEN
                          CALL NEXTWD('Class name>')
                          CLASINDEX = 0
                          DO COUNT = 1, NRAMACLASSES
                               IF (RAMACLASSNAMES(COUNT).EQ.CLASNAME)
     &                              CLASINDEX = COUNT
                          END DO
                          IF (CLASINDEX.EQ.0) THEN
                               CALL DSPERR('RAMA',
     &                              'unknown class. Using first.')
                               CLASINDEX = 1
                          END IF
                          CALL PRINTRAMAS(CUTOFF, HEAP(CALCRAMAPTR),
     &                         HEAP(RAMAATOMPTR), CLASINDEX)
                     ELSE
                          CALL DSPERR('RAMA',
     &                         'Expected ALL or CLASs.')
                     END IF
                END IF
C
C check for END statement
C
           ELSE
                CALL CHKEND('RAMA>', DONE)
           END IF
           END IF
      END DO
      DONE = .FALSE.
      CALL FREHP(SPTR,INTEG4(NATOM))
      RETURN
      END
C===============
      SUBROUTINE ALLOCEXPECTEDRAMA(CLASS, OLDSIZE, NEWSIZE)
C
C Allocates space for EXPECTEDRAMA
C based on values of PHISTEP and PSISTEP.
C Fills each array with zeros.
C
C by John Kuszewski May 1995
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE, CLASS
C local variables
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(EXPRAMAPTRS(CLASS), IREAL8(OLDSIZE))
      END IF
      EXPRAMAPTRS(CLASS) = 0
      IF (NEWSIZE.NE.0) THEN
        EXPRAMAPTRS(CLASS) = ALLHP(IREAL8(NEWSIZE))
C
C now zero them out
C
        CALL DEFEXPECTEDRAMA(HEAP(EXPRAMAPTRS(CLASS)), NEWSIZE)
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE ALLOCRAMAS (OLDSIZE, NEWSIZE)
C
C resets chem shift arrays to hold SIZE entries
C
C by John Kuszewski Nov 1995
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER OLDSIZE, NEWSIZE
C begin
      IF (OLDSIZE.NE.0) THEN
           CALL FREHP(RAMAATOMPTR, INTEG4(17*OLDSIZE))
           CALL FREHP(CALCRAMAPTR, IREAL8(5*OLDSIZE))
      END IF
C
      RAMAATOMPTR = 0
      CALCRAMAPTR = 0
C
      IF (NEWSIZE.NE.0) THEN
        RAMAATOMPTR = ALLHP(INTEG4(17*NEWSIZE))
        CALCRAMAPTR = ALLHP(IREAL8(5*NEWSIZE))
      END IF
C
      RETURN
      END
C===============
      SUBROUTINE RAMAINIT
C
C initializes Ramachandran stuff
C
C by John Kuszewski Nov 1995
C================
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C local vbls
C begin
      ALLOCRAMA=.TRUE.
      RETURN
      END
C==============
      SUBROUTINE RAMADEFAULTS
C
C sets up defaults
C
C by John Kuszewski Nov 1995
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C local variables
      INTEGER COUNT
C begin
      RAMAMODE = NEW
      MAXRAMAS = 2
      NRAMAS = 0
      NRAMACLASSES = 0
      CURRAMACLASS = 0
      RAMAPOTENTIAL = HARMONIC
      RAMASCALE = 1.0
      DO COUNT = 1, MAXRAMACLASSES
           PHISTEPS(COUNT) = 1
           PSISTEPS(COUNT) = 1
           CHISTEPS(COUNT) = 1
           THTSTEPS(COUNT) = 1
           RAMACLASSTYPE(COUNT) = 1
           RAMACLASSNAMES(COUNT) = 'DEFAULT'
           RAMAFORCES(COUNT) = 1.0
           RAMAPHASEC(1, COUNT) = 0.0
           RAMAPHASEC(2, COUNT) = 0.0
           RAMAPHASEC(3, COUNT) = 0.0
           RAMAPHASEC(4, COUNT) = 0.0
           RAMAPHASEC(5, COUNT) = 0.0
           RAMAPHASEC(6, COUNT) = 0.0
           RAMAPHASEC(7, COUNT) = 0.0
           RAMAPHASEC(8, COUNT) = 0.0
           RAMAPHASEC(9, COUNT) = 0.0
           RAMAPHASEC(10, COUNT) = 0.0
           RAMAPHASEC(11, COUNT) = 0.0
           RAMAPHASEC(12, COUNT) = 0.0
           CALL ALLOCEXPECTEDRAMA(COUNT,0,
     &          PHISTEPS(COUNT)*PSISTEPS(COUNT)*
     &          CHISTEPS(COUNT)*THTSTEPS(COUNT))
      END DO
      RETURN
      END
C==============
      SUBROUTINE DEFEXPECTEDRAMA (EXPECTED, SIZE)
C
C fills arrays that hold expectation values with zero
C
C by John Kuszewski June 1995
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
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
      SUBROUTINE READCOMPRESSEDRAMA4D (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS, CURCHIS, CURTHTS)
C
C Reads in the potential energy values for 4D torsion angle
C database in compressed format
C
C by John Kuszewski Apr 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CURCHIS, CURTHTS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS, CURCHIS, CURTHTS)
C local vbls
      INTEGER PHI, PSI, CHI, THT
      DOUBLE PRECISION VAL
C begin
      DO PHI = 1, CURPHIS
           DO PSI = 1, CURPSIS
                DO CHI = 1, CURCHIS
                     DO THT = 1, CURTHTS
                          CALL NEXTF('expected value =', VAL)
                          EXPECTED(PHI,PSI,CHI,THT) = VAL
                     END DO
                END DO
           END DO
      END DO
      RETURN
      END
C==============
      SUBROUTINE READCOMPRESSEDRAMA3D (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS, CURCHIS)
C
C Reads in the potential energy values for 3D torsion angle
C database in compressed format
C
C by John Kuszewski Apr 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CURCHIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS, CURCHIS)
C local vbls
      INTEGER PHI, PSI, CHI
      DOUBLE PRECISION VAL
C begin
      DO PHI = 1, CURPHIS
           DO PSI = 1, CURPSIS
                DO CHI = 1, CURCHIS
                     CALL NEXTF('expected value =', VAL)
                     EXPECTED(PHI,PSI,CHI) = VAL
                END DO
           END DO
      END DO
      RETURN
      END
C==============
      SUBROUTINE READCOMPRESSEDRAMA2D (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS)
C
C Reads in the potential energy values for 2D torsion angle
C database in compressed format
C
C by John Kuszewski Apr 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS)
C local vbls
      INTEGER PHI, PSI
      DOUBLE PRECISION VAL
C begin
      DO PHI = 1, CURPHIS
           DO PSI = 1, CURPSIS
                CALL NEXTF('expected value =', VAL)
                EXPECTED(PHI,PSI) = VAL
           END DO
      END DO
      RETURN
      END
C==============
      SUBROUTINE READCOMPRESSEDRAMA1D (CLASS, EXPECTED,
     &     CURPHIS)
C
C Reads in the potential energy values for 1D torsion angle
C database in compressed format
C
C by John Kuszewski Apr 1997
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS)
C local vbls
      INTEGER PHI
      DOUBLE PRECISION VAL
C begin
      DO PHI = 1, CURPHIS
           CALL NEXTF('expected value =', VAL)
           EXPECTED(PHI) = VAL
      END DO
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDRAMA4D (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS, CURCHIS, CURTHTS)
C
C actually reads in the potential energy values for 4D torsion angle
C database
C
C by John Kuszewski Nov 1996
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CURCHIS, CURTHTS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS, CURCHIS, CURTHTS)
C local vbls
      INTEGER PHIPOS, PSIPOS, CHIPOS, THTPOS
      DOUBLE PRECISION VAL
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS(CLASS)).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &          'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTI('PSIposition =', PSIPOS)
      IF ((PSIPOS.GT.PSISTEPS(CLASS)).OR.(PSIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &          'psi not in range 1..PSISTEPS. Assuming posn 1.')
           PSIPOS = 1
      END IF
      CALL NEXTI('CHIposition =', CHIPOS)
      IF ((CHIPOS.GT.CHISTEPS(CLASS)).OR.(CHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &          'chi not in range 1..CHISTEPS. Assuming posn 1.')
           CHIPOS = 1
      END IF
      CALL NEXTI('THETAposition =', THTPOS)
      IF ((THTPOS.GT.THTSTEPS(CLASS)).OR.(THTPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &          'theta not in range 1..THTSTEPS. Assuming posn 1.')
           THTPOS = 1
      END IF
      CALL NEXTF('expected value =', VAL)
      EXPECTED(PHIPOS,PSIPOS,CHIPOS,THTPOS) = VAL
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDRAMA3D (CLASS, EXPECTED,
     &     CURPHIS, CURPSIS, CURCHIS)
C
C actually reads in the potential energy values for 3D torsion angle
C database
C
C by John Kuszewski May 1995
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CURCHIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS, CURCHIS)
C local vbls
      INTEGER PHIPOS, PSIPOS, CHIPOS
      DOUBLE PRECISION VAL
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS(CLASS)).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTI('PSIposition =', PSIPOS)
      IF ((PSIPOS.GT.PSISTEPS(CLASS)).OR.(PSIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'psi not in range 1..PSISTEPS. Assuming posn 1.')
           PSIPOS = 1
      END IF
      CALL NEXTI('CHIposition =', CHIPOS)
      IF ((CHIPOS.GT.CHISTEPS(CLASS)).OR.(CHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'chi not in range 1..CHISTEPS. Assuming posn 1.')
           CHIPOS = 1
      END IF
      CALL NEXTF('expected value =', VAL)
      EXPECTED(PHIPOS,PSIPOS,CHIPOS) = VAL
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDRAMA2D (CLASS, EXPECTED, CURPHIS,
     &     CURPSIS)
C
C actually reads in the potential energy values for 2D torsion angle
C database
C
C by John Kuszewski May 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER CURPHIS, CURPSIS, CLASS
      DOUBLE PRECISION EXPECTED(CURPHIS, CURPSIS)
C local vbls
      INTEGER PHIPOS, PSIPOS
      DOUBLE PRECISION VAL
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS(CLASS)).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTI('PSIposition =', PSIPOS)
      IF ((PSIPOS.GT.PSISTEPS(CLASS)).OR.(PSIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'psi not in range 1..PSISTEPS. Assuming posn 1.')
           PSIPOS = 1
      END IF
      CALL NEXTF('expected value =', VAL)
      EXPECTED(PHIPOS,PSIPOS) = VAL
      RETURN
      END
C==============
      SUBROUTINE READEXPECTEDRAMA1D (CLASS, EXPECTED)
C
C actually reads in the potential energy values for 1D torsion angle
C database
C
C by John Kuszewski May 1994
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      DOUBLE PRECISION EXPECTED(*)
      INTEGER CLASS
C local vbls
      INTEGER PHIPOS
      DOUBLE PRECISION VAL
C begin
      CALL NEXTI('PHIposition =', PHIPOS)
      IF ((PHIPOS.GT.PHISTEPS(CLASS)).OR.(PHIPOS.LT.1)) THEN
           CALL DSPERR('RAMA',
     &               'phi not in range 1..PHISTEPS. Assuming posn 1.')
           PHIPOS = 1
      END IF
      CALL NEXTF('expected value =', VAL)
      EXPECTED(PHIPOS) = VAL
      RETURN
      END
C==============
      SUBROUTINE READRAMAASSIGN (ATOMS, CLASS, SEL)
C
C reads actual torsion angle assignments into array
C
C by John Kuszewski Nov 1993
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'numbers.inc'
C i/o
      INTEGER ATOMS(17,*), SEL(*), CLASS
C local variables
      INTEGER THETYPE, NSEL, COUNT, CURASSGN(17)
      LOGICAL USEME
C begin
C
      THETYPE = RAMACLASSTYPE(CLASS)
      USEME = .TRUE.
      CURASSGN(1) = CLASS
      DO COUNT = 2, 17
           CURASSGN(COUNT) = 0
      END DO
C
C now read the appropriate number of atoms
C
      DO COUNT = 1, (THETYPE*4)
           CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
           IF (NSEL.GT.1) THEN
                CALL DSPERR('RAMA',
     &               'more than 1 atom in selection. Using first')
           END IF
           CALL MAKIND(SEL, NATOM, NSEL)
           IF (NSEL.LT.1) USEME = .FALSE.
           CURASSGN(COUNT+1) = SEL(1)
      END DO
C
C see if we should ignore this one or not
C
      IF (USEME) THEN
C
C add it to the end of the constraint table and record its class
C
           NRAMAS = NRAMAS + 1
           DO COUNT = 1, 17
                ATOMS(COUNT, NRAMAS) = CURASSGN(COUNT)
           END DO
      ELSE
           CALL DSPERR('RAMA',
     &     'no atoms in one or more selections. Ignoring this entry.')
      END IF
      RETURN
      END
C=================
      SUBROUTINE SORTRAMAS (ATOMS)
C
C Given a list of rama constraints in random order, sort them by
C their class number.
C
C Based on SORTTRIOS in trio.s
C
C by John Kuszewski Jan 1997
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C i/o
      INTEGER ATOMS(17,*)
C local vbls
      INTEGER TEMP1(17), TEMP2(17)
      INTEGER NSORTEDRAMAS, CURCLASS, COUNT, COUNT2
C begin
      NSORTEDRAMAS = 0
      DO CURCLASS = 1, NRAMACLASSES
           DO COUNT = NSORTEDRAMAS + 1, NRAMAS
                IF (ATOMS(1,COUNT).EQ.CURCLASS) THEN
                     NSORTEDRAMAS = NSORTEDRAMAS + 1
C
C grab both points' data.  Implemented this way
C to ensure there are no cache problems
C
                     DO COUNT2 = 1, 17
                          TEMP1(COUNT2) = ATOMS(COUNT2,NSORTEDRAMAS)
                     END DO
                     DO COUNT2 = 1, 17
                          TEMP2(COUNT2) = ATOMS(COUNT2,COUNT)
                     END DO
                     DO COUNT2 = 1, 17
                          ATOMS(COUNT2,NSORTEDRAMAS) = TEMP2(COUNT2)
                     END DO
                     DO COUNT2 = 1, 17
                          ATOMS(COUNT2,COUNT) = TEMP1(COUNT2)
                     END DO
                END IF
           END DO
      END DO
      RETURN
      END
C=================
      SUBROUTINE PRINTRAMAS (CUTOFF, CALC, ATOMS, WHICHCLASS)
C
C prints torsion angle database entries with energies greater
C than CUTOFF.  Calculates RMS deviation and puts it into $RMS
C
C by John Kuszewski July 1995
C=================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'rama.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'numbers.inc'
C i/o
      DOUBLE PRECISION CUTOFF, CALC(5,*)
      INTEGER ATOMS(17,*), WHICHCLASS
C local variables
      DOUBLE PRECISION CURCALC, CURANG1, CURANG2, CURANG3, CURANG4,
     &     RMS, VIOLS, NUMINRMS, RENERGY, ETOT
      INTEGER CURCLASS, COUNT, CLASS, K, P, T, G, CURTYPE, CLEN
      LOGICAL PRINTTHISCLASS
      DOUBLE COMPLEX DUMMY2
C begin
      RMS = ZERO
      ETOT = ZERO
      VIOLS = ZERO
      NUMINRMS = ZERO
      CLASS = ZERO
C
C make sure that the calc-shift array is up to date
C
      CALL ERAMA(RENERGY, 'ANALYZE')
C
C write the overall header
C
      WRITE (PUNIT, '(A)') 'The following torsion angle DB entries'
      WRITE (PUNIT, '(A)') 'have energies greater than the cutoff:'
C
C loop through calc-Rama values
C
      DO COUNT = 1,NRAMAS
C
C is this the start of a new class?
C
           CURCLASS = ATOMS(1,COUNT)
           IF (CLASS.NE.CURCLASS) THEN
                CLASS = CURCLASS
                CLEN=LEN(RAMACLASSNAMES(CLASS))
                CALL TRIMM(RAMACLASSNAMES(CLASS),CLEN)
                IF ((WHICHCLASS.EQ.CLASS).OR.(WHICHCLASS.EQ.0)) THEN
                     PRINTTHISCLASS = .TRUE.
                ELSE
                     PRINTTHISCLASS = .FALSE.
                END IF
                CURTYPE = RAMACLASSTYPE(CLASS)
                IF (PRINTTHISCLASS) THEN
                     WRITE (PUNIT, '(A, A)')
     &                    'class ', RAMACLASSNAMES(CLASS)(1:CLEN)
                     IF (CURTYPE.EQ.1) THEN
                          WRITE (PUNIT, '(A)') '(atom K) (energy)'
                     ELSE IF (CURTYPE.EQ.2) THEN
                          WRITE (PUNIT, '(A)')
     &                         '(atom K) (atom P) (energy)'
                     ELSE IF (CURTYPE.EQ.3) THEN
                          WRITE (PUNIT, '(A)')
     &                         '(atom K) (atom P) (atomT) (energy)'
                     ELSE IF (CURTYPE.EQ.4) THEN
                          WRITE (PUNIT, '(A)')
     &                        '(atom K)(atom P)(atomT)(atom G)(energy)'
                     END IF
                END IF
           END IF
C
C get the current energy
C
           CURCALC = CALC(1, COUNT)
           CURANG1 = CALC(2, COUNT)
           CALL RD2DEG(CURANG1, CURANG1)
           CURANG2 = CALC(3, COUNT)
           CALL RD2DEG(CURANG2, CURANG2)
           CURANG3 = CALC(4, COUNT)
           CALL RD2DEG(CURANG3, CURANG3)
           CURANG4 = CALC(5, COUNT)
           CALL RD2DEG(CURANG4, CURANG4)
C
C if we should print this class and if
C we have an energy for this entry,
C
C why did Alexandre try to change this line?
C
           IF ((PRINTTHISCLASS).AND.(CURCALC.LT.NOEXPTEST)) THEN
C
C update RMS values if we can make a guess about this entry
C
                ETOT = ETOT + CURCALC
                RMS = RMS + CURCALC**2
                NUMINRMS = NUMINRMS + ONE
C
C if the current value is greater than the CUTOFF value
C and is less than the no-expectation-test value,
C
C print it
C
                IF (CURCALC.GT.CUTOFF) THEN
                     K = ATOMS(4, COUNT)
                     P = ATOMS(8, COUNT)
                     T = ATOMS(12, COUNT)
                     G = ATOMS(16, COUNT)
                     IF (CURTYPE.EQ.1) THEN
                          WRITE(PUNIT,'(4(A),2(F12.3))')
     &                         SEGID(K),RESID(K),RES(K),TYPE(K),
     &                         CURCALC, CURANG1
                     ELSE IF (CURTYPE.EQ.2) THEN
                          WRITE(PUNIT,'(8(A),3(F12.3))')
     &                         SEGID(K),RESID(K),RES(K),TYPE(K),
     &                         SEGID(P),RESID(P),RES(P),TYPE(P),
     &                         CURCALC, CURANG1, CURANG2
                     ELSE IF (CURTYPE.EQ.3) THEN
                          WRITE(PUNIT,'(12(A),4(F12.3))')
     &                         SEGID(K),RESID(K),RES(K),TYPE(K),
     &                         SEGID(P),RESID(P),RES(P),TYPE(P),
     &                         SEGID(T),RESID(T),RES(T),TYPE(T),
     &                         CURCALC, CURANG1, CURANG2, CURANG3
                     ELSE
                          WRITE(PUNIT,'(16(A),5(F12.3))')
     &                         SEGID(K),RESID(K),RES(K),TYPE(K),
     &                         SEGID(P),RESID(P),RES(P),TYPE(P),
     &                         SEGID(T),RESID(T),RES(T),TYPE(T),
     &                         SEGID(G),RESID(G),RES(G),TYPE(G),
     &                         CURCALC, CURANG1, CURANG2, CURANG3,
     &                         CURANG4

                     END IF
                     VIOLS = VIOLS + ONE
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
      WRITE (PUNIT, '(A, F15.3)') 'RMS energy = ', RMS
      WRITE (PUNIT, '(A, F15.3)') 'total energy = ', ETOT
      CALL DECLAR('RESULT', 'DP', ' ', DUMMY2, ETOT)
      CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
      CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      RETURN
      END
C==============
      SUBROUTINE RAMAHP
C
C deallocate memory
C
C by Alexandre Bonvin Apr 1996
C==============
      IMPLICIT NONE
C include files
      INCLUDE 'rama.inc'
C local variables
      INTEGER COUNT
C begin
      IF (.NOT.ALLOCRAMA) THEN
      CALL ALLOCRAMAS(MAXRAMAS,0)
      DO COUNT = 1, MAXRAMACLASSES
      CALL ALLOCEXPECTEDRAMA(COUNT,
     &          PHISTEPS(COUNT)*PSISTEPS(COUNT)*
     &          CHISTEPS(COUNT)*THTSTEPS(COUNT),0)
      END DO
      END IF
      RETURN
      END
C
      SUBROUTINE SCRRAMA
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'rama.inc'
C
C
C reset ramachandran database
      IF (.NOT.ALLOCRAMA) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: ramachandran database erased.'
      CALL RAMAHP
      CALL RAMAINIT
      END IF
      RETURN
      END
C
