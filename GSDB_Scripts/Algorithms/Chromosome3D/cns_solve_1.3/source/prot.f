C ================
      SUBROUTINE EPROTONSHIFT (EH, WHICH)
C
C Front end for EprotonShift2.
C
C by John Kuszewski Nov 1994
C ================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'hshifts.inc'
C i/o
      DOUBLE PRECISION EH
      CHARACTER*7 WHICH
C local vbls
C begin
C
      CALL EPROTONSHIFT2 (EH,
     &     HEAP(DHCDXPTR), HEAP(DHCDYPTR), HEAP(DHCDZPTR),
     &     HEAP(DHC1DXPTR), HEAP(DHC1DYPTR), HEAP(DHC1DZPTR),
     &     HEAP(DHC2DXPTR), HEAP(DHC2DYPTR), HEAP(DHC2DZPTR),
     &     HEAP(AAPTR), HEAP(ACENPTR), HEAP(AVFPTR), HEAP(ACPTR),
     &     HEAP(EEPTR), HEAP(EHPTR), HEAP(ECPTR),
     &     HEAP(RAPTR), HEAP(REPTR),
     &     HEAP(HRCPTR), HEAP(PROTPTR),
     &     HEAP(DHOBS1), HEAP(DHOBS2), HEAP(DHST1), HEAP(DHST2),
     &     HEAP(DHFI1), HEAP(DHFI2),
     &     WHICH)
      RETURN
      END
C ================
      SUBROUTINE EPROTONSHIFT2 (EH,
     &     DHCALCDX, DHCALCDY, DHCALCDZ,
     &     DHCALC1DX, DHCALC1DY, DHCALC1DZ,
     &     DHCALC2DX, DHCALC2DY, DHCALC2DZ,
     &     ANISATOMS, CENTERS, AVERAGEFACTORS, ANISCONSTANTS,
     &     ELECEXCEPTIONS, ELECHEAVIES, ELECCONSTS,
     &     RINGATOMS, RINGEXCEPTIONS,
     &     HSHIFTRC, PROTONS,
     &     PASSOBS1, PASSOBS2, PASSST1, PASSST2, PASSFI1, PASSFI2,
     &     WHICH)
C
C Calculates proton shift energy and forces. Uses approach outlined in
C Williamson and Asakura, J. Mag. Res. B vol 101, p. 63 (1993),
C and in the program by Williamson.
C
C Energy is of the form
C
C      Eprotonshift = k (Sobs - Srcoil - Smagn - Selec - Sring)^2
C
C for non-degenerate entries.  The degenerate version just has
C this wrapped in the Constantine energy function.
C
C Because each of the last three terms is complicated, they're
C calculated in separate subroutines.  Therefore, the final
C calculation of dervatives has to be put off until all of the
C terms have been evaluated.  Their individual contributions to
C the derivatives are stored in the DHCALCD* arrays.
C
C These derivatives are further modified to allow non-stereo-
C assigned protons to be refined, ala Constantine et al.,
C J. Magn. Res. B 108, 176 (1995)
C
C DELTAS holds the delta shift for each proton (indexed by
C proton atom number)
C
C Also handles the actual printing of calculated and observed
C shifts, in order to avoid having to set up another bunch of
C arrays to pass the intermediate results back to the printing
C routine.
C
C by John Kuszewski Oct 1994
C modified for non-stereo-assigned refinement, Feb 1996
C =================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'cnst.inc'
C i/o
      CHARACTER*7 WHICH
      DOUBLE PRECISION EH
      DOUBLE PRECISION
     &     CENTERS(*), AVERAGEFACTORS(*), ANISCONSTANTS(2,*),
     &     ELECCONSTS(*), HSHIFTRC(*),
     &     DHCALCDX(*), DHCALCDY(*), DHCALCDZ(*),
     &     DHCALC1DX(*), DHCALC1DY(*), DHCALC1DZ(*),
     &     DHCALC2DX(*), DHCALC2DY(*), DHCALC2DZ(*),
     &     PASSOBS1(*), PASSOBS2(*)
      INTEGER ANISATOMS(3,*), ELECHEAVIES(*), RINGATOMS(7,*),
     &     ELECEXCEPTIONS(*), RINGEXCEPTIONS(*), PROTONS(*),
     &     PASSST1(*),PASSST2(*),PASSFI1(*),PASSFI2(*)
C local vbls
      INTEGER HCOUNT, COUNT, H, CLASS, ASSIGNCOUNT, NAVGPROTONS,
     &     START, FINISH, DEG
      DOUBLE PRECISION ANIS, ELEC, RING, DELTASHIFT, TEMP,
     &     OBS, RCOIL, RMS, VIOLS, KPROT, KPROTCONST,
     &     ERR, DELTASHIFTSIGN, DELTASHIFT1, DELTASHIFT2,
     &     TOTANIS, TOTELEC, TOTRING, TOTCOIL, TOTOBS,
     &     TOTANIS1, TOTELEC1, TOTRING1, TOTCOIL1, TOTOBS1,
     &     TOTANIS2, TOTELEC2, TOTRING2, TOTCOIL2, TOTOBS2,
     &     NASSIGNSINCLUDED, DBPREC,
     &     o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14,
     &     htestterm1, htestterm2, htestterm3, heplus, dheplusdhcalc1,
     &     dheplusdhcalc2, heminus1, dheminus1dhcalc1,
     &     dheminus1dhcalc2, heminus2, dheminus2dhcalc1,
     &     dheminus2dhcalc2, heminus3, dheminus3dhcalc1,
     &     dheminus3dhcalc2, hcalc1, hcalc2, hobs1, hobs2, DS
      DOUBLE COMPLEX DUMMY2
      LOGICAL PRINTMODE, PRINTTHISCLASS
C
C begin
C
C if we're in printout mode, print the header
C and the first class heading
C
      IF (WHICH.EQ.'ANALYZE') THEN
           PRINTMODE = .TRUE.
      ELSE
           PRINTMODE = .FALSE.
      END IF
C
      CLASS = 1
      PRINTTHISCLASS = PRINTCLASS(CLASS)
      RMS = ZERO
      VIOLS = ZERO
      NASSIGNSINCLUDED = ZERO
C
      IF (PRINTMODE) THEN
           WRITE (PUNIT, '(A)') 'The following proton delta shifts are'
           WRITE (PUNIT, '(A)') 'greater than the cutoff:'
           WRITE (PUNIT, '(A60)')
     &     '(proton)(delta)(calc)(obs)(rcoil)(magn anis)(elec)(ring)'
           IF (PRINTTHISCLASS) THEN
                WRITE (PUNIT, '(A, A)') 'class ', HSHIFTCLASSNAMES(1)
           END IF
      END IF
C
C zero out holders for energy
C
      EH = ZERO
C
C for every observation entry,
C
      DO ASSIGNCOUNT = 1, NPROTASSIGNS
C
C is this the beginning of a new class?
C
           IF (HCLASSNDX(CLASS).LT.ASSIGNCOUNT) THEN
                CLASS = CLASS + 1
                PRINTTHISCLASS = PRINTCLASS(CLASS)
                IF (PRINTMODE.AND.PRINTTHISCLASS) THEN
                     WRITE (PUNIT, '(A, A)') 'class ',
     &                    HSHIFTCLASSNAMES(CLASS)
                END IF
           END IF
           KPROT = HSHIFTFORCES(1,CLASS)
           KPROTCONST = HSHIFTFORCES(2,CLASS)
           ERR = HSHIFTERRORS(CLASS)
           DEG = HSHIFTDEGENERACY(CLASS)
C
C if this is a normal class, deal with it the normal way
C
           IF (DEG.EQ.1) THEN
C
C zero out derivatives, deltashift, and
C holders for the averaged parts of the shift
C
                DO COUNT = 1,NATOM
                     DHCALCDX(COUNT) = ZERO
                     DHCALCDY(COUNT) = ZERO
                     DHCALCDZ(COUNT) = ZERO
                END DO
                DELTASHIFT = ZERO
                TOTANIS = ZERO
                TOTELEC = ZERO
                TOTRING = ZERO
                TOTCOIL = ZERO
                TOTOBS = ZERO
C
C get the starting and finishing indexes for
C the protons of this assignment
C
                START = PASSST1(ASSIGNCOUNT)
                FINISH = PASSFI1(ASSIGNCOUNT)
C
C get the number of protons to be averaged
C
                NAVGPROTONS = FINISH - START + 1
C
C get the observed shift for this assignment
C
                OBS = PASSOBS1(ASSIGNCOUNT)
C
C for every actual proton within this proton entry,
C
                DO HCOUNT = START, FINISH
C
C get the observed and random coil shifts for this proton
C
                     H = PROTONS(HCOUNT)
                     RCOIL = HSHIFTRC(H)
C
C call the individual terms
C
                     IF (TIMER.GT.2) CALL DSPCPU(' starting anis ')
                     CALL MAGNETANIS (ANISATOMS, H, CENTERS,
     &                    AVERAGEFACTORS, ANISCONSTANTS,
     &                    DHCALCDX, DHCALCDY, DHCALCDZ, ANIS,
     &                    HEAP(EEPTR), HEAP(APPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' starting elec ')
                     CALL ELECSHIFT (ELEC, ELECEXCEPTIONS, DHCALCDX,
     &                    DHCALCDY, DHCALCDZ, H,
     &                    ELECHEAVIES, ELECCONSTS)
                     IF (TIMER.GT.2) CALL DSPCPU(' starting ring ')
                     CALL RINGCURRENT (RINGATOMS, H, RINGEXCEPTIONS,
     &                    RING, DHCALCDX, DHCALCDY, DHCALCDZ,
     &                    HEAP(RCPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' done ring ')
C
C update the average values of the various terms
C
                     TOTANIS = TOTANIS + ANIS
                     TOTELEC = TOTELEC + ELEC
                     TOTRING = TOTRING + RING
                     TOTCOIL = TOTCOIL + RCOIL
                     TOTOBS = TOTOBS + OBS
                END DO
C
C Calculate the average values of the terms and the overall
C delta shift from the averages.
C
C This is equivalent to calculating delta shifts
C for each and averaging them.
C
                TOTANIS = TOTANIS / NAVGPROTONS
                TOTELEC = TOTELEC / NAVGPROTONS
                TOTRING = TOTRING / NAVGPROTONS
                TOTCOIL = TOTCOIL / NAVGPROTONS
                TOTOBS = TOTOBS / NAVGPROTONS
                DELTASHIFT = TOTANIS + TOTELEC + TOTRING +
     &               TOTCOIL - TOTOBS
                IF (DELTASHIFT.EQ.ZERO) THEN
                     DELTASHIFTSIGN = 1
                ELSE
                     DELTASHIFTSIGN = DELTASHIFT / ABS(DELTASHIFT)
                END IF
C
C If we need to print this, handle its stuff
C
                IF (PRINTMODE.AND.PRINTTHISCLASS) THEN
C
C always update the RMS difference and
C the number of assignments included in this RMS
C
                     RMS = RMS + DELTASHIFT**2
                     NASSIGNSINCLUDED = NASSIGNSINCLUDED + ONE
C
C print this one out if it's outside the cutoff
C
                     IF (ABS(DELTASHIFT).GT.PROTPRINTCUTOFF) THEN
C
C loop thru the protons that went into this
C observation, printing out their names
C
                          DO HCOUNT = START,(FINISH-1)
                               H = PROTONS(HCOUNT)
                               WRITE(PUNIT,'(4(A))')
     &                              SEGID(H),RESID(H),RES(H),TYPE(H)
                               IF (PRINTRMSD) RMSD(H) = ABS(DELTASHIFT)
                          END DO
C
C print out the final atom and
C the numbers for this observation
C
                          H = PROTONS(FINISH)
                          IF (PRINTRMSD) RMSD(H) = ABS(DELTASHIFT)
                          WRITE(PUNIT,'(4(A),7(F7.3))')
     &                       SEGID(H),RESID(H),RES(H),TYPE(H),
     &                       DELTASHIFT, (DELTASHIFT+TOTOBS), TOTOBS,
     &                       TOTCOIL, TOTANIS, TOTELEC, TOTRING
                          VIOLS = VIOLS + ONE
                     END IF
C
C If we're not in printout mode,
C calculate the energy and the derivatives
C
                ELSE
C
C handle square well mode
C
                     IF (HSHIFTPOTENTIAL(CLASS).EQ.SQUARE) THEN
                          DELTASHIFT = ABS(DELTASHIFT) - ERR
                          IF (DELTASHIFT.LT.ZERO) DELTASHIFT = ZERO
                          DELTASHIFT = DELTASHIFT * DELTASHIFTSIGN
                     END IF
C
C finally calculate energy and forces
C
                     EH = EH + KPROT * DELTASHIFT**2
                     TEMP = TWO * KPROT * DELTASHIFT
                     DO COUNT = 1,NATOM
                          DX(COUNT) = DX(COUNT) + TEMP *
     &                         DHCALCDX(COUNT)
                          DY(COUNT) = DY(COUNT) + TEMP *
     &                         DHCALCDY(COUNT)
                          DZ(COUNT) = DZ(COUNT) + TEMP *
     &                         DHCALCDZ(COUNT)
                     END DO
                END IF
           ELSE
C
C if this is a degenerate entry, deal with it using
C Constantine's trick
C
C
C zero out derivatives, deltashift, and
C the holders for the averaged parts of the shifts
C
                DO COUNT = 1,NATOM
                     DHCALC1DX(COUNT) = ZERO
                     DHCALC1DY(COUNT) = ZERO
                     DHCALC1DZ(COUNT) = ZERO
                     DHCALC2DX(COUNT) = ZERO
                     DHCALC2DY(COUNT) = ZERO
                     DHCALC2DZ(COUNT) = ZERO
                END DO
                DELTASHIFT = ZERO
                TOTANIS1 = ZERO
                TOTELEC1 = ZERO
                TOTRING1 = ZERO
                TOTCOIL1 = ZERO
                TOTOBS1 = ZERO
                TOTANIS2 = ZERO
                TOTELEC2 = ZERO
                TOTRING2 = ZERO
                TOTCOIL2 = ZERO
                TOTOBS2 = ZERO
C
C get the starting and finishing indexes for
C the protons of the first part of this assignment
C
                START = PASSST1(ASSIGNCOUNT)
                FINISH = PASSFI1(ASSIGNCOUNT)
C
C get the number of protons to be averaged
C
                NAVGPROTONS = FINISH - START + 1
C
C get the observed shift for this assignment
C
                OBS = PASSOBS1(ASSIGNCOUNT)
C
C for every actual proton within the first proton entry,
C
                DO HCOUNT = START, FINISH
C
C get the observed and random coil shifts for this proton
C
                     H = PROTONS(HCOUNT)
                     RCOIL = HSHIFTRC(H)
C
C call the individual terms
C
                     IF (TIMER.GT.2) CALL DSPCPU(' starting anis ')
                     CALL MAGNETANIS (ANISATOMS, H, CENTERS,
     &                    AVERAGEFACTORS, ANISCONSTANTS,
     &                    DHCALC1DX, DHCALC1DY, DHCALC1DZ, ANIS,
     &                    HEAP(EEPTR), HEAP(APPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' starting elec ')
                     CALL ELECSHIFT (ELEC, ELECEXCEPTIONS, DHCALC1DX,
     &                    DHCALC1DY, DHCALC1DZ, H,
     &                    ELECHEAVIES, ELECCONSTS)
                     IF (TIMER.GT.2) CALL DSPCPU(' starting ring ')
                     CALL RINGCURRENT (RINGATOMS, H, RINGEXCEPTIONS,
     &                    RING, DHCALC1DX, DHCALC1DY, DHCALC1DZ,
     &                    HEAP(RCPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' done ring ')
C
C update the average values of the various terms
C
                     TOTANIS1 = TOTANIS1 + ANIS
                     TOTELEC1 = TOTELEC1 + ELEC
                     TOTRING1 = TOTRING1 + RING
                     TOTCOIL1 = TOTCOIL1 + RCOIL
                     TOTOBS1 = TOTOBS1 + OBS
                END DO
C
C Calculate the average values of the terms and the overall
C delta shift from the averages.
C
C This is equivalent to calculating delta shifts
C for each and averaging them.
C
                TOTANIS1 = TOTANIS1 / NAVGPROTONS
                TOTELEC1 = TOTELEC1 / NAVGPROTONS
                TOTRING1 = TOTRING1 / NAVGPROTONS
                TOTCOIL1 = TOTCOIL1 / NAVGPROTONS
                TOTOBS1 = TOTOBS1 / NAVGPROTONS
C
C get the starting and finishing indexes for
C the protons of the second part of this assignment
C
                START = PASSST2(ASSIGNCOUNT)
                FINISH = PASSFI2(ASSIGNCOUNT)
C
C get the number of protons to be averaged
C
                NAVGPROTONS = FINISH - START + 1
C
C get the observed shift for this assignment
C
                OBS = PASSOBS2(ASSIGNCOUNT)
C
C for every actual proton within the second proton entry,
C
                DO HCOUNT = START, FINISH
C
C get the observed and random coil shifts for this proton
C
                     H = PROTONS(HCOUNT)
                     RCOIL = HSHIFTRC(H)
C
C call the individual terms
C
                     IF (TIMER.GT.2) CALL DSPCPU(' starting anis ')
                     CALL MAGNETANIS (ANISATOMS, H, CENTERS,
     &                    AVERAGEFACTORS, ANISCONSTANTS,
     &                    DHCALC2DX, DHCALC2DY, DHCALC2DZ, ANIS,
     &                    HEAP(EEPTR), HEAP(APPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' starting elec ')
                     CALL ELECSHIFT (ELEC, ELECEXCEPTIONS, DHCALC2DX,
     &                    DHCALC2DY, DHCALC2DZ, H,
     &                    ELECHEAVIES, ELECCONSTS)
                     IF (TIMER.GT.2) CALL DSPCPU(' starting ring ')
                     CALL RINGCURRENT (RINGATOMS, H, RINGEXCEPTIONS,
     &                    RING, DHCALC2DX, DHCALC2DY, DHCALC2DZ,
     &                    HEAP(RCPTR))
                     IF (TIMER.GT.2) CALL DSPCPU(' done ring ')
C
C update the average values of the various terms
C
                     TOTANIS2 = TOTANIS2 + ANIS
                     TOTELEC2 = TOTELEC2 + ELEC
                     TOTRING2 = TOTRING2 + RING
                     TOTCOIL2 = TOTCOIL2 + RCOIL
                     TOTOBS2 = TOTOBS2 + OBS
                END DO
C
C Calculate the average values of the terms and the overall
C delta shift from the averages.
C
C This is equivalent to calculating delta shifts
C for each and averaging them.
C
                TOTANIS2 = TOTANIS2 / NAVGPROTONS
                TOTELEC2 = TOTELEC2 / NAVGPROTONS
                TOTRING2 = TOTRING2 / NAVGPROTONS
                TOTCOIL2 = TOTCOIL2 / NAVGPROTONS
                TOTOBS2 = TOTOBS2 / NAVGPROTONS
C
C rename things for convenience
C
                HCALC1 = TOTANIS1 + TOTELEC1 + TOTRING1 + TOTCOIL1
                HCALC2 = TOTANIS2 + TOTELEC2 + TOTRING2 + TOTCOIL2
                HOBS1 = TOTOBS1
                HOBS2 = TOTOBS2
C
C If we need to print this, handle its stuff
C
                IF (PRINTMODE.AND.PRINTTHISCLASS) THEN
C
C This uses the same approach as the degenerate-couplings
C print routine.  Try both possibilities and keep the best one.
C
                     DELTASHIFT1 = ABS(HCALC1 - HOBS1) +
     &                    ABS(HCALC2 - HOBS2)
                     DELTASHIFT2 = ABS(HCALC1 - HOBS2) +
     &                    ABS(HCALC2 - HOBS1)
                     DELTASHIFT = MIN(DELTASHIFT1, DELTASHIFT2)
C
C always update the RMS difference and
C the number of observations included in this RMS
C
C count this as two entries for the purposes of RMS calculation
C so that you can compare the results of degenerate and non-
C degenerate classes
C
                     IF (DELTASHIFT.EQ.DELTASHIFT1) THEN
                          RMS = RMS + (HCALC1 - HOBS1)**2 +
     &                         (HCALC2 - HOBS2)**2
                     ELSE
                          RMS = RMS + (HCALC1 - HOBS2)**2 +
     &                         (HCALC2 - HOBS1)**2
                     END IF
                     NASSIGNSINCLUDED = NASSIGNSINCLUDED + TWO
C
C print this one out if it's outside the cutoff
C
                     IF (ABS(DELTASHIFT).GT.PROTPRINTCUTOFF) THEN
                          VIOLS = VIOLS + ONE
C
C pick the correct first observed shift
C
                          IF (DELTASHIFT.EQ.DELTASHIFT1) THEN
                               DS = HCALC1 - HOBS1
                          ELSE
                               DS = HCALC1 - HOBS2
                          END IF
C
C loop thru the protons that went into the first calculated shift,
C printing out their names
C
                          START = PASSST1(ASSIGNCOUNT)
                          FINISH = PASSFI1(ASSIGNCOUNT)
                          DO HCOUNT = START,(FINISH-1)
                               H = PROTONS(HCOUNT)
                               WRITE(PUNIT,'(4(A))')
     &                              SEGID(H),RESID(H),RES(H),TYPE(H)
                               IF (PRINTRMSD) RMSD(H) = DS
                          END DO
C
C print out the final atom and
C the numbers for the first calculated shift,
C including the proper observed value
C
                          H = PROTONS(FINISH)
                          IF (PRINTRMSD) RMSD(H) = DS
                          IF (DELTASHIFT.EQ.DELTASHIFT1) THEN
                               WRITE(PUNIT,'(4(A),7(F7.3))')
     &                            SEGID(H),RESID(H),RES(H),TYPE(H),
     &                            DS, (TOTOBS1+DS), TOTOBS1, TOTCOIL1,
     &                            TOTANIS1, TOTELEC1, TOTRING1
                          ELSE
                               WRITE(PUNIT,'(4(A),7(F7.3))')
     &                            SEGID(H),RESID(H),RES(H),TYPE(H),
     &                            DS, (TOTOBS2+DS), TOTOBS2, TOTCOIL1,
     &                            TOTANIS1, TOTELEC1, TOTRING1
                          END IF
C
C pick the correct second observed shift
C
                          IF (DELTASHIFT.EQ.DELTASHIFT1) THEN
                               DS = HCALC2 - HOBS2
                          ELSE
                               DS = HCALC2 - HOBS1
                          END IF
C
C loop thru the protons that went into the second calculated shift,
C printing out their names
C
                          START = PASSST2(ASSIGNCOUNT)
                          FINISH = PASSFI2(ASSIGNCOUNT)
                          DO HCOUNT = START,(FINISH-1)
                               H = PROTONS(HCOUNT)
                               WRITE(PUNIT,'(4(A))')
     &                              SEGID(H),RESID(H),RES(H),TYPE(H)
                               IF (PRINTRMSD) RMSD(H) = DS
                          END DO
C
C print out the final atom and
C the numbers for the second calculated shift,
C including the proper observed value
C
                          H = PROTONS(FINISH)
                          IF (PRINTRMSD) RMSD(H) = DS
                          IF (DELTASHIFT.EQ.DELTASHIFT1) THEN
                               WRITE(PUNIT,'(4(A),7(F7.3))')
     &                            SEGID(H),RESID(H),RES(H),TYPE(H),
     &                            DS, (TOTOBS2+DS), TOTOBS2, TOTCOIL2,
     &                            TOTANIS2, TOTELEC2, TOTRING2
                          ELSE
                               WRITE(PUNIT,'(4(A),7(F7.3))')
     &                            SEGID(H),RESID(H),RES(H),TYPE(H),
     &                            DS, (TOTOBS1+DS), TOTOBS1, TOTCOIL2,
     &                            TOTANIS2, TOTELEC2, TOTRING2
                          END IF
                     END IF
                ELSE
C
C otherwise, calculate energy & forces
C
C
C now turn these numbers--Hcalc1, dHcalc1/dxi, Hcalc2, dHcalc2/dxi
C into the Constantine energy and its derivatives wrt coordinates
C
C This is obviously from Mma
C
                     o1 = Hcalc1 - Hcalc2
                     o2 = o1**2
                     o3 = Sqrt(o2)
                     o4 = Hobs1 - Hobs2
                     o5 = o4**2
                     o6 = Sqrt(o5)
                     o7 = Hcalc1 + Hcalc2 - Hobs1 - Hobs2
                     o8 = kprot*o7
                     o9 = o3 - o6
                     o10 = 1/o3
                     o11 = kprot*o1*o10*o9
                     o12 = -o3 + o6
                     o13 = kprot*kprotconst*o1*o10*o12
                     o14 = kprot*kprotconst*o1
                     Htestterm1 = o3
                     Htestterm2 = o6
                     Htestterm3 = o6/2
                     HEplus = kprot*o7**2
                     DHEplusDHcalc1 = 2*o8
                     DHEplusDHcalc2 = 2*o8
                     HEminus1 = kprot*o9**2
                     DHEminus1DHcalc1 = 2*o11
                     DHEminus1DHcalc2 = -2*o11
                     HEminus2 = kprot*kprotconst*o12**2
                     DHEminus2DHcalc1 = -2*o13
                     DHEminus2DHcalc2 = 2*o13
                     HEminus3 = kprot*kprotconst*(-o2 + o5/2)
                     DHEminus3DHcalc1 = -2*o14
                     DHEminus3DHcalc2 = 2*o14
C
C apply the tests to get the Constantine energy
C and produce forces:
C
C E = Eplus + Eminus
C
C dE/dxi = dEplus/dxi + dEminus/dxi
C
C        = dEplus/dHcalc1 dHcalc1/dxi + dEplus/dHcalc2 dHcalc2/dxi +
C
C          dEminus/dHcalc1 dHcalc1/dxi + dEminus/dHcalc2 dHcalc2/dxi
C
C see my calculus book p 812
C
                     IF (HTESTTERM1.GT.HTESTTERM2) THEN
                          EH = EH + HEPLUS + HEMINUS1
                          DO COUNT = 1,NATOM
                          DX(COUNT) = DX(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DX(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DX(COUNT) +
     &                         DHEminus1DHcalc1 * DHCALC1DX(COUNT) +
     &                         DHEminus1DHcalc2 * DHCALC2DX(COUNT)
                          DY(COUNT) = DY(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DY(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DY(COUNT) +
     &                         DHEminus1DHcalc1 * DHCALC1DY(COUNT) +
     &                         DHEminus1DHcalc2 * DHCALC2DY(COUNT)
                          DZ(COUNT) = DZ(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DZ(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DZ(COUNT) +
     &                         DHEminus1DHcalc1 * DHCALC1DZ(COUNT) +
     &                         DHEminus1DHcalc2 * DHCALC2DZ(COUNT)
                          END DO
                     ELSE IF (HTESTTERM1.LT.HTESTTERM3) THEN
                          EH = EH + HEPLUS + HEMINUS3
                          DO COUNT = 1,NATOM
                          DX(COUNT) = DX(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DX(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DX(COUNT) +
     &                         DHEminus3DHcalc1 * DHCALC1DX(COUNT) +
     &                         DHEminus3DHcalc2 * DHCALC2DX(COUNT)
                          DY(COUNT) = DY(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DY(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DY(COUNT) +
     &                         DHEminus3DHcalc1 * DHCALC1DY(COUNT) +
     &                         DHEminus3DHcalc2 * DHCALC2DY(COUNT)
                          DZ(COUNT) = DZ(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DZ(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DZ(COUNT) +
     &                         DHEminus3DHcalc1 * DHCALC1DZ(COUNT) +
     &                         DHEminus3DHcalc2 * DHCALC2DZ(COUNT)
                          END DO
                     ELSE
                          EH = EH + HEPLUS + HEMINUS2
                          DO COUNT = 1,NATOM
                          DX(COUNT) = DX(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DX(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DX(COUNT) +
     &                         DHEminus2DHcalc1 * DHCALC1DX(COUNT) +
     &                         DHEminus2DHcalc2 * DHCALC2DX(COUNT)
                          DY(COUNT) = DY(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DY(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DY(COUNT) +
     &                         DHEminus2DHcalc1 * DHCALC1DY(COUNT) +
     &                         DHEminus2DHcalc2 * DHCALC2DY(COUNT)
                          DZ(COUNT) = DZ(COUNT) +
     &                         DHEplusDHcalc1   * DHCALC1DZ(COUNT) +
     &                         DHEplusDHcalc2   * DHCALC2DZ(COUNT) +
     &                         DHEminus2DHcalc1 * DHCALC1DZ(COUNT) +
     &                         DHEminus2DHcalc2 * DHCALC2DZ(COUNT)
                          END DO
                     END IF
C
C end of if over whether this is a print or energy call
C
                END IF
C
C end of if over whether this is a normal or Constantine entry
C
           END IF
C
C end of the main loop over the observation entries
C
      END DO
C
C If we're in printout mode, print the summary information
C
      IF (PRINTMODE) THEN
           IF (NASSIGNSINCLUDED.GT.0) THEN
                RMS = SQRT(RMS / NASSIGNSINCLUDED)
           ELSE
                RMS = ZERO
           END IF
           WRITE (PUNIT, '(A, F15.3)') 'RMS error = ', RMS
           WRITE (PUNIT, '(A, F9.3)')
     &          'Number of violations: ', VIOLS
           DBPREC = NASSIGNSINCLUDED
           CALL DECLAR('NUMBER', 'DP', ' ', DUMMY2, DBPREC)
           CALL DECLAR('RMS', 'DP', ' ', DUMMY2, RMS)
           CALL DECLAR('VIOLATIONS', 'DP', ' ', DUMMY2, VIOLS)
      END IF
      RETURN
      END
C===============
      SUBROUTINE MAGNETANIS (ANISATOMS, PROTON, CENTERS,
     &     AVERAGEFACTORS, CONSTANTS, DSHIFTDX, DSHIFTDY,
     &     DSHIFTDZ, ANIS, AMIDES, ANISPOSITIONS)
C
C Calculates magnetic anisotropy effects on proton
C chemical shifts.
C
C ANISATOMS holds, in order, the three atoms for each
C     backbone Ca C O
C     asn Cb Cg Od1
C     gln Cg Cd Oe1
C     asp Cb Cg Od1/Od2
C     glu Cg Cd Oe1/Oe2
C     backbone O C N (for the C-N anisotropy)
C
C PROTON is the atom number of the proton whose magnetic
C anisotropy is to be calculated.
C
C AVERAGEFACTORS holds the factor to be multiplied by
C the calculated anisotropy (asp sidechain groups' two
C C=Os effects are averaged together).  Ie., for asp and
C glu sidechain COs, averagefactor = 0.5.  In all other
C cases, it's 1.0
C
C CENTERS holds a list of bond centers for each entry
C in ANISATOMS, expressed as a fraction of the C=O
C bond length.
C
C CONSTANTS holds a list of deltachi1 and deltachi2 for
C each entry in ANISATOMS.  The values should be:
C     all C=Os deltaX1 = -13.0 deltaX2 = -4.0
C     all C-Ns    "      -11.0    "      +1.40
C
C AMIDES holds a list of the amide protons (which matter
C for C-N anisotropy--see below)
C
C ANISPOSITIONS holds a list a integers, indexed by trio
C number, which is 1 if that trio's from a sidechain,
C zero otherwise.
C
C The output is ANIS, holding the magnetic anisotropy-caused
C shift on PROTON, and the DSHIFTD* variables, which hold the
C partial derivatives of the anisotropy wrt each atom coordinate
C that was used.  It's indexed by atom number and defined for
C all atoms.
C
C by John Kuszewski Oct 1994
C===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'coord.inc'
C i/o
      INTEGER ANISATOMS(3,*), PROTON,
     &     AMIDES(*), ANISPOSITIONS(*)
      DOUBLE PRECISION AVERAGEFACTORS(*), CONSTANTS(2, *), CENTERS(*),
     &     DSHIFTDX(*), DSHIFTDY(*), DSHIFTDZ(*), ANIS
C local vbls
      INTEGER ANISCOUNT, CA, CO, O, H
      LOGICAL SAME, ISEXCEP, ISAMIDE
      DOUBLE PRECISION DELTACHI1, DELTACHI2, AVEFACT, DIST
      DOUBLE PRECISION CAX, CAY, CAZ, COX, COY, COZ,
     &     OX, OY, OZ,
     &     HX, HY, HZ, RLEN
C
      DOUBLE PRECISION o1, o2, o3, o4, o5, o6, o7, o8, o9,
     &     o10, o11, o12, o13, o14, o15, o16, o17, o18, o19,
     &     o20, o21, o22, o23, o24, o25, o26, o27, o28, o29,
     &     o30, o31, o32, o33, o34, o35, o36, o37, o38, o39,
     &     o40, o41, o42, o43, o44, o45, o46, o47, o48, o49,
     &     o50, o51, o52, o53, o54, o55, o56, o57, o58, o59,
     &     o60, o61, o62, o63, o64, o65, o66, o67, o68, o69,
     &     o70, o71, o72, o73, o74, o75, o76, o77, o78, o79,
     &     o80, o81, o82, o83, o84, o85, o86, o87, o88, o89,
     &     o90, o91, o92, o93, o94, o95, o96, o97, o98, o99,
     &     o100, o101, o102, o103, o104, o105, o106, o107, o108, o109,
     &     o110, o111, o112, o113, o114, o115, o116, o117, o118, o119,
     &     o120, o121, o122, o123, o124, o125, o126, o127, o128, o129,
     &     o130, o131, o132, o133, o134, o135, o136, o137, o138, o139,
     &     o140, o141, o142, o143, o144, o145, o146, o147, o148, o149,
     &     o150, o151, o152, o153, o154, o155, o156, o157, o158
C
      DOUBLE PRECISION DANISDCAX, DANISDCAY, DANISDCAZ,
     &     DANISDCOX, DANISDCOY, DANISDCOZ, DANISDOX,
     &     DANISDOY, DANISDOZ, DANISDHX, DANISDHY, DANISDHZ
C
C begin
C
C zero out output
C
      ANIS = ZERO
C
C get the proton's coordinates
C
      H = PROTON
      HX = X(H)
      HY = Y(H)
      HZ = Z(H)
C
C for every Ca, C, O trio in the list generated by GenerateAnisList
C for this proton,
C
      DO ANISCOUNT = 1, NANIS
C
C get the atom numbers for the next trio
C
           CA = ANISATOMS(1,ANISCOUNT)
           CO = ANISATOMS(2,ANISCOUNT)
           O  = ANISATOMS(3,ANISCOUNT)
C
C skipping rules are different for CO and CN anisotropies
C
           IF (CENTERS(ANISCOUNT).EQ.CNCENTER) THEN
C
C for CN anis, only skip HNs that are in the same residue
C as the N, which is held in the O variable
C ElecExceptions holds the list of amides
C
                CALL ISANEXCEPTION(H,AMIDES,NAMIDES,ISAMIDE)
                CALL SAMERESIDUE(H,O,SAME)
                IF (ISAMIDE.AND.SAME) THEN
                     ISEXCEP = .TRUE.
                ELSE
                     ISEXCEP = .FALSE.
                END IF
C
C for sidechain CO anisotropies, don't do any skipping
C
           ELSE IF (ANISPOSITIONS(ANISCOUNT).EQ.1) THEN
                ISEXCEP = .FALSE.
C
C For backbone CO anis,
C skip all protons that are in the same residue as CO.
C
           ELSE
                CALL SAMERESIDUE(H,CO,ISEXCEP)
           END IF
C
C Calculate the anis if this isn't a trio to be skipped
C
           IF (.NOT.ISEXCEP) THEN
C
C get the various constants for this anisotropic bond
C
                DELTACHI1 = CONSTANTS(1,ANISCOUNT)
                DELTACHI2 = CONSTANTS(2,ANISCOUNT)
                AVEFACT = AVERAGEFACTORS(ANISCOUNT)
                DIST = CENTERS(ANISCOUNT)
C
C get the trio's coordinates
C
                CAX = X(CA)
                CAY = Y(CA)
                CAZ = Z(CA)
                COX = X(CO)
                COY = Y(CO)
                COZ = Z(CO)
                OX  = X(O)
                OY  = Y(O)
                OZ  = Z(O)
C
C everything is now from Mma---
C
                o1=-cox+ox
                o2=dist*o1
                o3=-cox+hx-o2
                o4=o3**2
                o5=-coy+oy
                o6=dist*o5
                o7=-coy+hy-o6
                o8=o7**2
                o9=-coz+oz
                o10=dist*o9
                o11=-coz+hz-o10
                o12=o11**2
                o13=o12+o4+o8
                rlen = sqrt(o13)
C
C if this proton is within the cutoff,
C
                IF (RLEN.LE.ANISCUTOFF) THEN
C
                     o14=o13**(3.d0/2.d0)
                     o15=1/o14
                     o16=cay-coy
                     o17=o1*o16
                     o18=cax-cox
                     o19=o18*o5
                     o20=o17-o19
                     o21=o20**2
                     o22=caz-coz
                     o23=o1*o22
                     o24=o18*o9
                     o25=-o23+o24
                     o26=o25**2
                     o27=o22*o5
                     o28=o16*o9
                     o29=o27-o28
                     o30=o29**2
                     o31=o21+o26+o30
                     o32=1/o31
                     o33=1/o13
                     o34=o25*o3
                     o35=o29*o7
                     o36=o34-o35
                     o37=o36**2
                     o38=o20*o7
                     o39=o11*o25
                     o40=o38-o39
                     o41=o40**2
                     o42=o20*o3
                     o43=o11*o29
                     o44=-o42+o43
                     o45=o44**2
                     o46=o37+o41+o45
                     o47=o32*o33*o46
                     o48=deltachi1*(1.d0-3.d0*(1.d0-o47))
                     o49=o20*o5
                     o50=o25*o9
                     o51=o49-o50
                     o52=o51**2
                     o53=o1*o25
                     o54=o29*o5
                     o55=o53-o54
                     o56=o55**2
                     o57=o1*o20
                     o58=o29*o9
                     o59=-o57+o58
                     o60=o59**2
                     o61=o52+o56+o60
                     o62=1/o61
                     o63=o11*o51
                     o64=o3*o55
                     o65=o63-o64
                     o66=o65**2
                     o67=o51*o7
                     o68=o3*o59
                     o69=-o67+o68
                     o70=o69**2
                     o71=o55*o7
                     o72=o11*o59
                     o73=o71-o72
                     o74=o73**2
                     o75=o66+o70+o74
                     o76=o33*o62*o75
                     o77=deltachi2*(1.d0-3.d0*(1.d0-o76))
                     o78=o48+o77
                     o79=coy-oy
                     o80=o31**2
                     o81=1/o80
                     o82=o5*o79
                     o83=o9**2
                     o84=o82-o83
                     o85=o61**2
                     o86=1/o85
                     o87=coz-oz
                     o88=o1**2
                     o89=o87*o9
                     o90=-o88+o89
                     o91=cox-ox
                     o92=o1*o91
                     o93=o5**2
                     o94=o92-o93
                     o95=caz-oz
                     o96=(-1.d0+dist)*o25
                     o97=-cay+oy
                     o98=(-1.d0+dist)*o20
                     o99=o13**2
                     o100=1/o99
                     o101=o1*o95
                     o102=o101+o23-o24
                     o103=o5*o97
                     o104=o9*o95
                     o105=o103-o104
                     o106=(-1.d0+dist)*o55
                     o107=o1*o97
                     o108=-o107+o17-o19
                     o109=(-1.d0+dist)*o59
                     o110=o13**(5.d0/2.d0)
                     o111=1/o110
                     o112=-caz+oz
                     o113=(-1.d0+dist)*o29
                     o114=cax-ox
                     o115=o112*o5
                     o116=-o115+o27-o28
                     o117=o114*o5
                     o118=o117-o17+o19
                     o119=o1*o114
                     o120=o112*o9
                     o121=-o119+o120
                     o122=(-1.d0+dist)*o51
                     o123=-cax+ox
                     o124=cay-oy
                     o125=o1*o123
                     o126=o124*o5
                     o127=o125-o126
                     o128=o123*o9
                     o129=-o128-o23+o24
                     o130=o124*o9
                     o131=o130-o27+o28
                     o132=-caz+coz
                     o133=dist*o25
                     o134=dist*o20
                     o135=o1*o132
                     o136=o135-o23+o24
                     o137=o16*o5
                     o138=o132*o9
                     o139=o137-o138
                     o140=dist*o55
                     o141=-2.d0*o17+o19
                     o142=dist*o59
                     o143=dist*o29
                     o144=-cax+cox
                     o145=-2.d0*o27+o28
                     o146=o144*o5
                     o147=o146+o17-o19
                     o148=o1*o144
                     o149=o22*o9
                     o150=-o148+o149
                     o151=dist*o51
                     o152=-cay+coy
                     o153=o1*o18
                     o154=o152*o5
                     o155=o153-o154
                     o156=o23-2.d0*o24
                     o157=o152*o9
                     o158=o157+o27-o28
C
C here's the anisotropy and its derivatives
C
            anis = anis -3.33333333333333d-1*avefact*o15*o78
            danisdcax = -3.33333333333333d-1*avefact*o15*
     &           (-3.d0*deltachi1*(o33*o46*(2.d
     &  0*o50+2.d0*o20*o79)*o81-o32*o33*(-2.d0*o3*o44*o79+2.d0*o3*o36*o
     &  9+2.d0*o40*(o7*o79-o11*o9)))-3.d0*deltachi2*(o33*o75*o86*(-2.d0
     &  *o1*o59*o79+2.d0*o51*o84+2.d0*o1*o55*o9)-o33*o62*(2.d0*o69*(-(o
     &  1*o3*o79)-o7*o84)+2.d0*o65*(o11*o84-o1*o3*o9)+2.d0*o73*(o1*o11*
     &  o79+o1*o7*o9))))
            danisdcay = -3.33333333333333d-1*avefact*o15*
     &           (-3.d0*deltachi1*(o33*o46*o81*
     &  (2.d0*o57+2.d0*o29*o87)-o32*o33*(2.d0*o1*o40*o7-2.d0*o36*o7*o87
     &  +2.d0*o44*(-(o1*o3)+o11*o87)))-3.d0*deltachi2*(o33*o75*o86*(2.d
     &  0*o1*o5*o51-2.d0*o5*o55*o87+2.d0*o59*o90)-o33*o62*(2.d0*o65*(o1
     &  *o11*o5+o3*o5*o87)+2.d0*o73*(-(o5*o7*o87)-o11*o90)+2.d0*o69*(-(
     &  o1*o5*o7)+o3*o90))))
            danisdcaz = -3.33333333333333d-1*avefact*o15*
     &           (-3.d0*deltachi1*(o33*o46*o81*
     &  (2.d0*o54+2.d0*o25*o91)-o32*o33*(2.d0*o11*o44*o5-2.d0*o11*o40*o
     &  91+2.d0*o36*(-(o5*o7)+o3*o91)))-3.d0*deltachi2*(o33*o75*o86*(2.
     &  d0*o5*o59*o9-2.d0*o51*o9*o91+2.d0*o55*o94)-o33*o62*(2.d0*o69*(o
     &  3*o5*o9+o7*o9*o91)+2.d0*o65*(-(o11*o9*o91)-o3*o94)+2.d0*o73*(-(
     &  o11*o5*o9)+o7*o94))))
            danisdcox = avefact*(-1.d0+dist)*o111*o3*
     &           o78-3.33333333333333d-1*avefact*o1
     &  5*(-3.d0*deltachi2*(-(o33*o62*(2.d0*(-o106+o105*o11-o102*o3)*o6
     &  5+2.d0*o69*(o109+o108*o3-o105*o7)+2.d0*(-(o108*o11)+o102*o7)*o7
     &  3))+2.d0*(-1.d0+dist)*o100*o3*o62*o75+o33*(2.d0*o105*o51+2.d0*o
     &  102*o55+2.d0*o108*o59)*o75*o86)-3.d0*deltachi1*(2.d0*(-1.d0+dis
     &  t)*o100*o3*o32*o46+o33*o46*o81*(2.d0*o25*o95+2.d0*o20*o97)-o32*
     &  o33*(2.d0*o36*(o3*o95+o96)+2.d0*o40*(-(o11*o95)+o7*o97)+2.d0*o4
     &  4*(-(o3*o97)-o98))))
            danisdcoy = avefact*(-1.d0+dist)*o111*o7*
     &           o78-3.33333333333333d-1*avefact*o1
     &  5*(-3.d0*deltachi2*(-(o33*o62*(2.d0*(o11*o118-o116*o3)*o65+2.d0
     &  *o69*(-o122+o121*o3-o118*o7)+2.d0*(o106-o11*o121+o116*o7)*o73))
     &  +2.d0*(-1.d0+dist)*o100*o62*o7*o75+o33*(2.d0*o118*o51+2.d0*o116
     &  *o55+2.d0*o121*o59)*o75*o86)-3.d0*deltachi1*(2.d0*(-1.d0+dist)*
     &  o100*o32*o46*o7+(2.d0*o114*o20+2.d0*o112*o29)*o33*o46*o81-o32*o
     &  33*(2.d0*(o11*o112-o114*o3)*o44+2.d0*o36*(-o113-o112*o7)+2.d0*o
     &  40*(o114*o7+o98))))
            danisdcoz = avefact*(-1.d0+dist)*o11*o111*
     &           o78-3.33333333333333d-1*avefact*o
     &  15*(-3.d0*deltachi2*(-(o33*o62*(2.d0*(o122+o11*o129-o127*o3)*o6
     &  5+2.d0*o69*(o131*o3-o129*o7)+2.d0*(-o109-o11*o131+o127*o7)*o73)
     &  )+2.d0*(-1.d0+dist)*o100*o11*o62*o75+o33*(2.d0*o129*o51+2.d0*o1
     &  27*o55+2.d0*o131*o59)*o75*o86)-3.d0*deltachi1*(2.d0*(-1.d0+dist
     &  )*o100*o11*o32*o46+(2.d0*o123*o25+2.d0*o124*o29)*o33*o46*o81-o3
     &  2*o33*(2.d0*(o113+o11*o124)*o44+2.d0*o36*(o123*o3-o124*o7)+2.d0
     &  *o40*(-(o11*o123)-o96))))
            danisdox = -(avefact*dist*o111*o3*o78)-
     &           3.33333333333333d-1*avefact*o15*(-3
     &  .d0*deltachi1*(-2.d0*dist*o100*o3*o32*o46-o32*o33*(2.d0*(-o133+
     &  o132*o3)*o36+2.d0*(o134-o16*o3)*o44+2.d0*o40*(-(o11*o132)+o16*o
     &  7))+(2.d0*o16*o20+2.d0*o132*o25)*o33*o46*o81)-3.d0*deltachi2*(-
     &  (o33*o62*(2.d0*(o11*o139+o140-o136*o3)*o65+2.d0*o69*(-o142+o141
     &  *o3-o139*o7)+2.d0*(-(o11*o141)+o136*o7)*o73))-2.d0*dist*o100*o3
     &  *o62*o75+o33*(2.d0*o139*o51+2.d0*o136*o55+2.d0*o141*o59)*o75*o8
     &  6))
            danisdoy = -(avefact*dist*o111*o7*o78)-
     &           3.33333333333333d-1*avefact*o15*(-3
     &  .d0*deltachi1*(-2.d0*dist*o100*o32*o46*o7-o32*o33*(2.d0*(o11*o2
     &  2-o144*o3)*o44+2.d0*o40*(-o134+o144*o7)+2.d0*o36*(o143-o22*o7))
     &  +(2.d0*o144*o20+2.d0*o22*o29)*o33*o46*o81)-3.d0*deltachi2*(-(o3
     &  3*o62*(2.d0*(o11*o147-o145*o3)*o65+2.d0*o69*(o151+o150*o3-o147*
     &  o7)+2.d0*(-o140-o11*o150+o145*o7)*o73))-2.d0*dist*o100*o62*o7*o
     &  75+o33*(2.d0*o147*o51+2.d0*o145*o55+2.d0*o150*o59)*o75*o86))
            danisdoz = -(avefact*dist*o11*o111*o78)-
     &           3.33333333333333d-1*avefact*o15*(-
     &  3.d0*deltachi1*(-2.d0*dist*o100*o11*o32*o46-o32*o33*(2.d0*(o133
     &  -o11*o18)*o40+2.d0*(-o143+o11*o152)*o44+2.d0*o36*(o18*o3-o152*o
     &  7))+(2.d0*o18*o25+2.d0*o152*o29)*o33*o46*o81)-3.d0*deltachi2*(-
     &  (o33*o62*(2.d0*(-o151+o11*o156-o155*o3)*o65+2.d0*o69*(o158*o3-o
     &  156*o7)+2.d0*(o142-o11*o158+o155*o7)*o73))-2.d0*dist*o100*o11*o
     &  62*o75+o33*(2.d0*o156*o51+2.d0*o155*o55+2.d0*o158*o59)*o75*o86)
     &  )
            danisdhx = -3.33333333333333d-1*avefact*
     &           o15*(-3.d0*deltachi1*(-(o32*o33*(2
     &  .d0*o25*o36+2.d0*(-o17+o19)*o44))+2.d0*o100*o3*o32*o46)-3.d0*de
     &  ltachi2*(-(o33*o62*(2.d0*(-o53+o54)*o65+2.d0*o59*o69))+2.d0*o10
     &  0*o3*o62*o75))+avefact*o111*o3*o78
            danisdhy = -3.33333333333333d-1*avefact*
     &           o15*(-3.d0*deltachi1*(-(o32*o33*(2
     &  .d0*(-o27+o28)*o36+2.d0*o20*o40))+2.d0*o100*o32*o46*o7)-3.d0*de
     &  ltachi2*(-(o33*o62*(2.d0*(-o49+o50)*o69+2.d0*o55*o73))+2.d0*o10
     &  0*o62*o7*o75))+avefact*o111*o7*o78
            danisdhz = -3.33333333333333d-1*avefact*
     &           o15*(-3.d0*deltachi1*(-(o32*o33*(2
     &  .d0*(o23-o24)*o40+2.d0*o29*o44))+2.d0*o100*o11*o32*o46)-3.d0*de
     &  ltachi2*(-(o33*o62*(2.d0*o51*o65+2.d0*(o57-o58)*o73))+2.d0*o100
     &  *o11*o62*o75))+avefact*o11*o111*o78
C
C now update the arrays of derivatives
C
                     DSHIFTDX(CA) = DSHIFTDX(CA) + DANISDCAX
                     DSHIFTDY(CA) = DSHIFTDY(CA) + DANISDCAY
                     DSHIFTDZ(CA) = DSHIFTDZ(CA) + DANISDCAZ
                     DSHIFTDX(CO) = DSHIFTDX(CO) + DANISDCOX
                     DSHIFTDY(CO) = DSHIFTDY(CO) + DANISDCOY
                     DSHIFTDZ(CO) = DSHIFTDZ(CO) + DANISDCOZ
                     DSHIFTDX(O) = DSHIFTDX(O) + DANISDOX
                     DSHIFTDY(O) = DSHIFTDY(O) + DANISDOY
                     DSHIFTDZ(O) = DSHIFTDZ(O) + DANISDOZ
                     DSHIFTDX(H) = DSHIFTDX(H) + DANISDHX
                     DSHIFTDY(H) = DSHIFTDY(H) + DANISDHY
                     DSHIFTDZ(H) = DSHIFTDZ(H) + DANISDHZ
                END IF
           END IF
      END DO
      RETURN
      END
C =================
      SUBROUTINE DIRECTLYBONDED (A, B)
C
C Given an atom A, use the bonding arrays IB and JB
C to find an atom B directly bonded to it.
C
C by John Kuszewski Nov 1994
C =================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER A, B
C local vbls
      INTEGER COUNT
C begin
C
C define a default value in case one isn't found
C
      B = 0
C
C Check through the IB and JB lists for A, and take
C B to be its partner
C
      DO COUNT = 1,NBOND
           IF (IB(COUNT).EQ.A) THEN
                B = JB(COUNT)
           END IF
           IF (JB(COUNT).EQ.A) THEN
                B = IB(COUNT)
           END IF
      END DO
      RETURN
      END
C =================
      SUBROUTINE SAMERESIDUE (A, B, ANSWER)
C
C Are atoms A and B in the same residue?
C
C by John Kuszewski Nov 1994
C =================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C i/o
      INTEGER A, B
      LOGICAL ANSWER
C local vbls
C begin
      IF ((RESID(A).EQ.RESID(B)).AND.
     &     (SEGID(A).EQ.SEGID(B))) THEN
           ANSWER = .TRUE.
      ELSE
           ANSWER = .FALSE.
      END IF
      RETURN
      END
C =================
      SUBROUTINE ELECSHIFT (ELEC, ELECEXCEPTIONS, DSHIFTDX,
     &     DSHIFTDY, DSHIFTDZ, PROTON, ELECHEAVIES, ELECCONSTS)
C
C Calculates electric field effects on proton shifts
C
C PROTON is the atom number of the proton to be calculated.
C elecExceptions is a list of all HNs (which don't get calculated).
C elecHeavys  is a list of backbone CO, O, N, and HNs, in arb order.
C elecConsts  is a list of Qc values for all heavy atoms
C             (indexed by atom number)  The values are
C              carbonyl carbon  1.60
C                     nitrogen -1.70
C                       oxygen -2.30
C                 amide proton  0.70
C
C The output is ELEC, the overall shift on PROTON caused by the
C electric field effects, and the DELECD* arrays, which hold the
C partial derivatives of the electric field effect.  They're
C indexed by atom number and defined over all atoms.
C
C by John Kuszewski Oct 1994
C ==================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
C i/o
      DOUBLE PRECISION ELEC, DSHIFTDX(*), DSHIFTDY(*),
     &     DSHIFTDZ(*), ELECCONSTS(*)
      INTEGER ELECHEAVIES(*), ELECEXCEPTIONS(*), PROTON
C local vbls
      INTEGER H, C, J, HEAVYCOUNT
      DOUBLE PRECISION HX, HY, HZ, CX, CY, CZ,
     &     Q, JX, JY, JZ, HJLEN,
     &     DEDHX, DEDHY, DEDHZ, DEDCX, DEDCY, DEDCZ,
     &     DEDJX, DEDJY, DEDJZ,
     &     O1,   O2,  O3,  O4,  O5,  O6,  O7,  O8,  O9, O10,
     &     O11, O12, O13, O14, O15, O16, O17, O18, O19, O20,
     &     O21, O22, O23, O24, O25, O26, O27, O28, O29, O30,
     &     O31, O32
      LOGICAL SAME, EXCEP
C begin
C
C
C get the proton and make sure it's not an amide.
C
      H = PROTON
      CALL ISANEXCEPTION(H, ELECEXCEPTIONS, NAMIDES, EXCEP)
      IF (.NOT.EXCEP) THEN
C
C get its directly bonded carbon.
C
           CALL DIRECTLYBONDED(H,C)
           HX = X(H)
           HY = Y(H)
           HZ = Z(H)
           CX = X(C)
           CY = Y(C)
           CZ = Z(C)
C
C for every heavy atom,
C
           DO HEAVYCOUNT = 1, NHEAVIES
                J = ELECHEAVIES(HEAVYCOUNT)
C
C skip this heavy atom's influence if it's in the same
C resiude as the proton
C
                CALL SAMERESIDUE(H, J, SAME)
                IF (.NOT.SAME) THEN
C
C get its atom type dependent constant and combine
C it with the overall constant
C
                     Q = ELECCONSTS(HEAVYCOUNT) * SIGMAE
                     JX = X(J)
                     JY = Y(J)
                     JZ = Z(J)
                     O1=HX-JX
                     O2=O1**2
                     O3=HY-JY
                     O4=O3**2
                     O5=HZ-JZ
                     O6=O5**2
                     O7=O2+O4+O6
                     HJLEN = SQRT(O7)
C
C if the distance is within the cutoff,
C
                     IF (HJLEN.LE.ELECCUTOFF) THEN
                          O8=-CX+HX
                          O9=O8**2
                          O10=-CY+HY
                          O11=O10**2
                          O12=-CZ+HZ
                          O13=O12**2
                          O14=O11+O13+O9
                          O15=SQRT(O14)
                          O16=1/O15
                          O17=O1*O8
                          O18=O10*O3
                          O19=O12*O5
                          O20=O17+O18+O19
                          O21=O7**(3.D0/2.D0)
                          O22=1/O21
                          O23=O14**(3.D0/2.D0)
                          O24=1/O23
                          O25=O20*O22*O24*O8*Q
                          O26=O10*O20*O22*O24*Q
                          O27=O12*O20*O22*O24*Q
                          O28=O7**(5.D0/2.D0)
                          O29=1/O28
                          O30=O1*O16*O20*O29*Q
                          O31=O16*O20*O29*O3*Q
                          O32=O16*O20*O29*O5*Q
                          ELEC = O16*O20*O22*Q
                          DEDCX = O25+(-HX+JX)*O16*O22*Q
                          DEDCY = O26+(-HY+JY)*O16*O22*Q
                          DEDCZ = O27+(-HZ+JZ)*O16*O22*Q
                          DEDHX = -O25-3.D0*O30+
     &                         (-CX+2.D0*HX-JX)*O16*O22*Q
                          DEDHY = -O26-3.D0*O31+
     &                         (-CY+2.D0*HY-JY)*O16*O22*Q
                          DEDHZ = -O27-3.D0*O32+
     &                         (-CZ+2.D0*HZ-JZ)*O16*O22*Q
                          DEDJX = 3.D0*O30+(CX-HX)*O16*O22*Q
                          DEDJY = 3.D0*O31+(CY-HY)*O16*O22*Q
                          DEDJZ = 3.D0*O32+(CZ-HZ)*O16*O22*Q
C
C     now update the shift derivative arrays
C
                          DSHIFTDX(H) = DSHIFTDX(H) + DEDHX
                          DSHIFTDY(H) = DSHIFTDY(H) + DEDHY
                          DSHIFTDZ(H) = DSHIFTDZ(H) + DEDHZ
                          DSHIFTDX(C) = DSHIFTDX(C) + DEDCX
                          DSHIFTDY(C) = DSHIFTDY(C) + DEDCY
                          DSHIFTDZ(C) = DSHIFTDZ(C) + DEDCZ
                          DSHIFTDX(J) = DSHIFTDX(J) + DEDJX
                          DSHIFTDY(J) = DSHIFTDY(J) + DEDJY
                          DSHIFTDZ(J) = DSHIFTDZ(J) + DEDJZ
                     END IF
                END IF
           END DO
      END IF
      RETURN
      END
C ================
      SUBROUTINE ISANEXCEPTION (A, EXCEPTIONS, NEXCEPTIONS, ANSWER)
C
C is the given atom among those in the list of exceptions?
C
C by John Kuszewski Nov 1994
C ================
      IMPLICIT NONE
C include files
      INCLUDE 'hshifts.inc'
C i/o
      INTEGER A, EXCEPTIONS(*), NEXCEPTIONS
      LOGICAL ANSWER
C local vbls
      INTEGER COUNT
C begin
      ANSWER = .FALSE.
      DO COUNT = 1, NEXCEPTIONS
           IF (EXCEPTIONS(COUNT).EQ.A) THEN
                ANSWER = .TRUE.
           END IF
      END DO
      RETURN
      END
C ================
      SUBROUTINE RINGCURRENT (RINGATOMS, PROTON, EXCEPTIONS,
     &     RINGSHIFT, DSHIFTDX, DSHIFTDY, DSHIFTDZ,
     &     RINGCONSTS)
C
C Calculates ring current shift.
C
C RINGATOMS has positions for up to six atoms of the ring.
C The seventh position of RINGATOMS holds the number of
C atoms in that ring.  The order of the atoms is
C
C     Phe CG  CD1 CE1 CZ  CE2 CD2
C     Tyr CG  CD1 CE1 CZ  CE2 CD2
C     His CG  CD2 NE2 CE1 ND1
C     Trp CG  CD1 NE1 CE2 CD2
C     Trp CE2 CD2 CE3 CZ3 CH2 CZ2
C
C PROTON is the atom number of the proton whose ring-current-
C caused shift is to be calculated.
C
C EXCEPTIONS is a list of all Ha and Hn protons.
C
C RINGSHIFT is the output--the shift of PROTON expected
C due to all the ring currents.
C
C DRINGD* is indexed by atom number and defined for
C all atoms.
C
C by John Kuszewski Oct 1994
C ================
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'coord.inc'
C i/o
      INTEGER RINGATOMS(7,*), PROTON, EXCEPTIONS(*)
      DOUBLE PRECISION RINGSHIFT, DSHIFTDX(*),
     &     DSHIFTDY(*), DSHIFTDZ(*), RINGCONSTS(*)
C local variables
      LOGICAL SAME, ISEXCEP
      INTEGER RINGCOUNT, A, B, C, NATOMSINRING
      integer H, COUNT
      INTEGER CEN1, CEN2, CEN3
      DOUBLE PRECISION AX, AY, AZ, BX, BY, BZ, CX, CY, CZ
      DOUBLE PRECISION HX, HY, HZ, RINGCONSTI
      DOUBLE PRECISION CEN1X, CEN1Y, CEN1Z, CEN2X, CEN2Y, CEN2Z
      DOUBLE PRECISION CEN3X, CEN3Y, CEN3Z, CENX, CENY, CENZ, HCENLEN
      DOUBLE PRECISION O1, O2, O3, O4, O5, O6, O7, O8, O9,
     &     O10, O11, O12, O13, O14, O15, O16, O17, O18, O19,
     &     O20, O21, O22, O23, O24, O25, O26, O27, O28, O29,
     &     O30, O31, O32, O33, O34, O35, O36, O37, O38, O39,
     &     O40, O41, O42, O43, O44, O45, O46, O47, O48, O49,
     &     O50, O51, O52, O53, O54, O55, O56, O57, O58, O59,
     &     O60, O61, O62, O63, O64, O65, O66, O67, O68, O69,
     &     O70, O71
      DOUBLE PRECISION DSHIFTDAX, DSHIFTDAY, DSHIFTDAZ
      DOUBLE PRECISION DSHIFTDBX, DSHIFTDBY, DSHIFTDBZ
      DOUBLE PRECISION DSHIFTDCX, DSHIFTDCY, DSHIFTDCZ
      DOUBLE PRECISION DSHIFTDHX, DSHIFTDHY, DSHIFTDHZ
      DOUBLE PRECISION SIGN
C
C begin
C
C initialize output to zero
C
      RINGSHIFT = ZERO
C
C Get the coords for the proton of interest
C and determine if it's an exception to the
C no-self-shifts rule
C
      H = PROTON
      HX = X(H)
      HY = Y(H)
      HZ = Z(H)
      CALL ISANEXCEPTION (H, EXCEPTIONS, NRINGEXCEPS, ISEXCEP)
C
C for every ring in the list generated by GenerateRingList
C for this proton,
C
      DO RINGCOUNT = 1, NRINGS
C
C get the atoms that are used to calculate the ring's center
C
           CEN1 = RINGATOMS(1,RINGCOUNT)
           CEN2 = RINGATOMS(3,RINGCOUNT)
           CEN3 = RINGATOMS(5,RINGCOUNT)
           RINGCONSTI = RINGCONSTS(RINGCOUNT)
C
C only bother with this if it's not in the same residue as
C PROTON or if PROTON is an Ha or HN proton, in which case
C the residue doesn't matter
C
           CALL SAMERESIDUE (CEN1, H, SAME)
           IF (.NOT.(SAME .AND. .NOT.ISEXCEP)) THEN
                NATOMSINRING = RINGATOMS(7,RINGCOUNT)
                cen1x = X(cen1)
                cen1Y = Y(cen1)
                cen1Z = Z(cen1)
                cen2X = X(cen2)
                cen2Y = Y(cen2)
                cen2Z = Z(cen2)
                cen3X = X(cen3)
                cen3Y = Y(cen3)
                cen3Z = Z(cen3)
                cenx= (cen1X + cen2X + cen3X) / 3.0D0
                ceny= (cen1Y + cen2Y + cen3Y) / 3.0D0
                cenz= (cen1Z + cen2Z + cen3Z) / 3.0D0
                hcenlen = SQRT((HX - cenx)**2 + (HY - ceny)**2 +
     &               (HZ - cenz)**2)
                IF (HCENLEN.LE.RINGCUTOFF) THEN
C
C loop over pairs of neighboring atoms in the ring,
C
                     DO COUNT = 1,NATOMSINRING
                          A = RINGATOMS(COUNT,RINGCOUNT)
                          IF (COUNT.EQ.NATOMSINRING) THEN
                               B = RINGATOMS(1,RINGCOUNT)
                               C = RINGATOMS(2,RINGCOUNT)
                          ELSE IF (COUNT.EQ.(NATOMSINRING-1)) THEN
                               B = RINGATOMS(COUNT + 1, RINGCOUNT)
                               C = RINGATOMS(1,RINGCOUNT)
                          ELSE
                               B = RINGATOMS(COUNT + 1, RINGCOUNT)
                               C = RINGATOMS(COUNT + 2,RINGCOUNT)
                          END IF
                          AX = X(A)
                          AY = Y(A)
                          AZ = Z(A)
                          BX = X(B)
                          BY = Y(B)
                          BZ = Z(B)
                          CX = X(C)
                          CY = Y(C)
                          CZ = Z(C)
C
                          o1=ay-by
                          o2=-bx+cx
                          o3=o1*o2
                          o4=ax-bx
                          o5=-by+cy
                          o6=o4*o5
                          o7=-o3+o6
                          o8=o7**2
                          o9=az-bz
                          o10=o2*o9
                          o11=-bz+cz
                          o12=o11*o4
                          o13=o10-o12
                          o14=o13**2
                          o15=o5*o9
                          o16=o1*o11
                          o17=-o15+o16
                          o18=o17**2
                          o19=o14+o18+o8
                          o20=sqrt(o19)
                          o21=1/o20
                          o22=bx-hx
                          o23=ay-hy
                          o24=o22*o23
                          o25=ax-hx
                          o26=by-hy
                          o27=o25*o26
                          o28=-o24+o27
                          o29=o28*o7
                          o30=az-hz
                          o31=o22*o30
                          o32=bz-hz
                          o33=o25*o32
                          o34=o31-o33
                          o35=o13*o34
                          o36=o26*o30
                          o37=o23*o32
                          o38=-o36+o37
                          o39=o17*o38
                          o40=o29+o35+o39
                          o41=o25**2
                          o42=o23**2
                          o43=o30**2
                          o44=o41+o42+o43
                          o45=o44**(3.d0/2.d0)
                          o46=1/o45
                          o47=o22**2
                          o48=o26**2
                          o49=o32**2
                          o50=o47+o48+o49
                          o51=o50**(3.d0/2.d0)
                          o52=1/o51
                          o53=o46+o52
                          o54=bz-cz
                          o55=o19**(3.d0/2.d0)
                          o56=1/o55
                          o57=o44**(5.d0/2.d0)
                          o58=1/o57
                          o59=bx-cx
                          o60=by-cy
                          o61=ay-cy
                          o62=-az+cz
                          o63=o50**(5.d0/2.d0)
                          o64=1/o63
                          o65=-ax+cx
                          o66=az-cz
                          o67=ax-cx
                          o68=-ay+cy
                          o69=-ay+by
                          o70=-az+bz
                          o71=-ax+bx
C
       IF ((o21*o40).LT.ZERO) THEN
            SIGN = 1
       ELSE
            SIGN = +1
       END IF
C
       ringshift = ringshift + 5.d-1*o21*o40*o53*
     &      ringconstB*ringconstI*sign
       dshiftdax = -1.5d0*o21*o25*o40*o58*ringconstB*
     &      ringconstI*sign+5.d-1*o21*o53
     &      *((-bz+hz)*o13+o28*o5+o34*o54+o26*o7)*
     &      ringconstB*ringconstI*sig
     &      n-2.5d-1*o40*o53*o56*(2.d0*o13*o54+2.d0*o5*o7)*
     &      ringconstB*ringconstI*sign
       dshiftday = -1.5d0*o21*o23*o40*o58*ringconstB*
     &      ringconstI*sign+5.d-1*o21*o53
     &      *(o17*o32+o11*o38+o28*o59+(-bx+hx)*o7)*
     &      ringconstB*ringconstI*sign-2.5d-1*o40*
     &      o53*o56*(2.d0*o11*o17+
     &      2.d0*o59*o7)*ringconstB*ringconstI*sign
       dshiftdaz = -1.5d0*o21*o30*o40*o58*ringconstB*
     &      ringconstI*sign-2.5d-1*o40*o53*o56*
     &      (2.d0*o13*o2+2.d0*o17*o60)*ringconstB*
     &      ringconstI*sign+5.d-1*o21*o53*((-by+hy)*
     &      o17+o13*o22+o2*o34+o38*o60)*ringconstB*ringconstI*sign
       dshiftdbx = -1.5d0*o21*o22*o40*o64*ringconstB*
     &      ringconstI*sign+5.d-1*o21*o53
     &      *(o13*o30+o28*o61+o34*o62+(-ay+hy)*o7)*
     &      ringconstB*ringconstI*sign-2.5d-1*o40*o53*
     &      o56*(2.d0*o13*o62+2.d0*o61*o7)*ringconstB*ringconstI*sign
       dshiftdby = -1.5d0*o21*o26*o40*o64*ringconstB*
     &      ringconstI*sign+5.d-1*o21*o53
     &      *((-az+hz)*o17+o28*o65+o38*o66+o25*o7)*
     &      ringconstB*ringconstI*sign-2.5d-1*o40*
     &      o53*o56*(2.d0*o17*o66+2.d0*o65*o7)*ringconstB*
     &      ringconstI*sign
       dshiftdbz = -1.5d0*o21*o32*o40*o64*ringconstB*
     &      ringconstI*sign-2.5d-1*o40*o53*o56*(2.d0*o13*
     &      o67+2.d0*o17*o68)*ringconstB*ringconstI*sign+
     &      5.d-1*o21*o53*((-ax+hx)*o13+o17*o23+o34*o67+
     &      o38*o68)*ringconstB*ringconstI*sign
       dshiftdcx = -2.5d-1*o40*o53*o56*(2.d0*o69*o7+2.d0*
     &      o13*o9)*ringconstB*ringconstI*sign+5.d-1*o21*
     &      o53*(o28*o69+o34*o9)*ringconstB*ringconstI*sign
       dshiftdcy = -2.5d-1*o40*o53*o56*(2.d0*o4*o7+2.d0*o17*
     &      o70)*ringconstB*ringconstI*sign+5.d-1*o21*o53*
     &      (o28*o4+o38*o70)*ringconstB*ringconstI*sign
       dshiftdcz = -2.5d-1*o40*o53*o56*(2.d0*o1*o17+2.d0*
     &      o13*o71)*ringconstB*ringconstI*sign+5.d-1*o21*
     &      o53*(o1*o38+o34*o71)*ringconstB*ringconstI*sign
       dshiftdhx = 5.d-1*o21*o40*(3.d0*o25*o58+3.d0*o22*o64)*
     &      ringconstB*ringconstI*sign+5.d-1*o21*o53*(o1*o7+o13*
     &      o70)*ringconstB*ringconstI*sign
       dshiftdhy = 5.d-1*o21*o40*(3.d0*o23*o58+3.d0*o26*o64)*
     &      ringconstB*ringconstI*sign+5.d-1*o21*o53*(o7*o71+
     &      o17*o9)*ringconstB*ringconstI*sign
       dshiftdhz = 5.d-1*o21*o40*(3.d0*o30*o58+3.d0*
     &      o32*o64)*ringconstB*ringconstI*sign+5.d-1*o21*o53*
     &      (o13*o4+o17*o69)*ringconstB*ringconstI*sign
C
C now update the ring current shift derivative arrays
C
                          DSHIFTDX(A) = DSHIFTDX(A) + DSHIFTDAX
                          DSHIFTDY(A) = DSHIFTDY(A) + DSHIFTDAY
                          DSHIFTDZ(A) = DSHIFTDZ(A) + DSHIFTDAZ
                          DSHIFTDX(B) = DSHIFTDX(B) + DSHIFTDBX
                          DSHIFTDY(B) = DSHIFTDY(B) + DSHIFTDBY
                          DSHIFTDZ(B) = DSHIFTDZ(B) + DSHIFTDBZ
                          DSHIFTDX(C) = DSHIFTDX(C) + DSHIFTDCX
                          DSHIFTDY(C) = DSHIFTDY(C) + DSHIFTDCY
                          DSHIFTDZ(C) = DSHIFTDZ(C) + DSHIFTDCZ
                          DSHIFTDX(H) = DSHIFTDX(H) + DSHIFTDHX
                          DSHIFTDY(H) = DSHIFTDY(H) + DSHIFTDHY
                          DSHIFTDZ(H) = DSHIFTDZ(H) + DSHIFTDHZ
                     END DO
                END IF
           END IF
      END DO
      RETURN
      END
C ===============
      SUBROUTINE READHSHIFT
C
C front end for READSHIFT2
C
C by John Kuszewski Nov 1994
C ===============
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'hshifts.inc'
C local vbls
      INTEGER SPTR, S1PTR, S2PTR
C begin
C
C allocate memory if we haven't already done so
C
      IF (.NOT.HSHIFTALLOCATED) THEN
      CALL ALLOCHSHIFTS
      END IF
      SPTR=ALLHP(INTEG4(NATOM))
      S1PTR=ALLHP(INTEG4(NATOM))
      S2PTR=ALLHP(INTEG4(NATOM))
      CALL READHSHIFT2 (HEAP(SPTR), HEAP(S1PTR), HEAP(S2PTR),
     &     HEAP(PROTPTR), HEAP(DHOBS1),HEAP(DHOBS2),
     &     HEAP(DHST1),HEAP(DHST2), HEAP(DHFI1),HEAP(DHFI2),
     &     HEAP(HRCPTR), HEAP(AAPTR), HEAP(ACENPTR), HEAP(ACPTR),
     &     HEAP(AVFPTR), HEAP(EEPTR), HEAP(EHPTR), HEAP(ECPTR),
     &     HEAP(RAPTR), HEAP(REPTR), HEAP(APPTR), HEAP(RCPTR))
      CALL FREHP(SPTR,INTEG4(NATOM))
      CALL FREHP(S1PTR,INTEG4(NATOM))
      CALL FREHP(S2PTR,INTEG4(NATOM))
      RETURN
      END
C ===============
      SUBROUTINE READHSHIFT2 (SEL, SEL1, SEL2, PROTONS,
     &     PASSOBS1, PASSOBS2, PASSST1, PASSST2, PASSFI1, PASSFI2,
     &     HSHIFTRC, ANISATOMS, CENTERS, ANISCONSTANTS, ANISAVEFACTORS,
     &     AMIDES, ELECHEAVIES, ELECCONSTS, RINGATOMS,
     &     RINGEXCEPTIONS, ANISPOSITIONS, RINGCONSTS)
C
C reads in proton chemical shift information
C
C by John Kuszewski Nov 1994
C ===============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'numbers.inc'
      INCLUDE 'coord.inc'
C i/o
      INTEGER SEL(*), SEL1(*), SEL2(*), PROTONS(*), ANISATOMS(3,*),
     &     AMIDES(*), ELECHEAVIES(*), RINGATOMS(7,*),
     &     RINGEXCEPTIONS(*), ANISPOSITIONS(*),
     &     PASSST1(*), PASSST2(*), PASSFI1(*), PASSFI2(*)
      DOUBLE PRECISION HSHIFTRC(*), CENTERS(*),
     &     ANISCONSTANTS(2,*), ANISAVEFACTORS(*), ELECCONSTS(*),
     &     RINGCONSTS(*), PASSOBS1(*), PASSOBS2(*)
C local variables
      INTEGER COUNT, NSEL,
     &     CA, CO, O, NRA, A1, A2, A3, A4, A5, A6,
     &     POS, NSEL1, NSEL2
      DOUBLE PRECISION K1, OBS1, OBS2, RC, CEN,
     &     X1, X2, AVE, CURRINGI
      CHARACTER*4 THENAME
      LOGICAL SELOK
C
C begin
C
C now read input
C
      CALL PUSEND('PROTONSHIFTS>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('PROTONSHIFTS>')
      CALL MISCOM('PROTONSHIFTS>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-proton')
C
C get the random coil shift
C
      ELSE IF (WD(1:4).EQ.'RCOI') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      CALL NEXTF('random coil proton shift =', RC)
      DO COUNT = 1,NSEL
      HSHIFTRC(SEL(COUNT)) = RC
      END DO
C
C get the anisotropy entry
C
      ELSE IF (WD(1:4).EQ.'ANIS') THEN
      SELOK = .TRUE.
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.NE.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'not 1 atom in selection. Ignoring entry')
      SELOK = .FALSE.
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      CA = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.NE.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'not 1 atom in selection. Ignoring entry')
      SELOK = .FALSE.
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      CO = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.NE.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'not 1 atom in selection. Ignoring entry')
      SELOK = .FALSE.
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      O = SEL(1)
      CALL NEXTA4('bond type =', THENAME)
      IF (THENAME(1:2).EQ.'CO') THEN
      CEN = COCENTER
      X1 = CODELTAX1
      X2 = CODELTAX2
      ELSE IF (THENAME(1:2).EQ.'CN') THEN
      CEN = CNCENTER
      X1 = CNDELTAX1
      X2 = CNDELTAX2
      ELSE
      CALL DSPERR('PROTONSHIFTS',
     &'unknown bond type.  Using CO')
      CEN = COCENTER
      X1 = CODELTAX1
      X2 = CODELTAX2
      END IF
      CALL NEXTA4('is carboxyl =', THENAME)
      IF (THENAME(1:4).EQ.'COOH') THEN
      AVE = 0.5
      ELSE
      AVE = 1.0
      END IF
      CALL NEXTA4('position =', THENAME)
      IF (THENAME(1:2).EQ.'BB') THEN
      POS = 0
      ELSE
      POS = 1
      END IF
      IF (SELOK) THEN
      NANIS = NANIS + 1
      ANISATOMS(1,NANIS) = CA
      ANISATOMS(2,NANIS) = CO
      ANISATOMS(3,NANIS) = O
      CENTERS(NANIS) = CEN
      ANISCONSTANTS(1,NANIS) = X1
      ANISCONSTANTS(2,NANIS) = X2
      ANISAVEFACTORS(NANIS) = AVE
      ANISPOSITIONS(NANIS) = POS
      END IF
C
C get the list of amide protons
C
      ELSE IF (WD(1:4).EQ.'AMID') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      DO COUNT = 1,NSEL
      NAMIDES = NAMIDES + 1
      AMIDES(NAMIDES) = SEL(COUNT)
      NHEAVIES = NHEAVIES + 1
      ELECHEAVIES(NHEAVIES) = SEL(COUNT)
      ELECCONSTS(NHEAVIES) = QHN
      END DO
C
C get the list of carbons
C
      ELSE IF (WD(1:4).EQ.'CARB') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      DO COUNT = 1,NSEL
      NHEAVIES = NHEAVIES + 1
      ELECHEAVIES(NHEAVIES) = SEL(COUNT)
      ELECCONSTS(NHEAVIES) = QC
      END DO
C
C get the list of nitrogens
C
      ELSE IF (WD(1:4).EQ.'NITR') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      DO COUNT = 1,NSEL
      NHEAVIES = NHEAVIES + 1
      ELECHEAVIES(NHEAVIES) = SEL(COUNT)
      ELECCONSTS(NHEAVIES) = QN
      END DO
C
C get the list of oxygens
C
      ELSE IF (WD(1:4).EQ.'OXYG') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      DO COUNT = 1,NSEL
      NHEAVIES = NHEAVIES + 1
      ELECHEAVIES(NHEAVIES) = SEL(COUNT)
      ELECCONSTS(NHEAVIES) = QO
      END DO
C
C get the ring current entries
C
      ELSE IF (WD(1:4).EQ.'RING') THEN
      CALL NEXTA4('ring type =', THENAME)
      IF (THENAME(1:3).EQ.'PHE') THEN
      NRA = 6
      CURRINGI = PHERINGI
      ELSE IF (THENAME(1:3).EQ.'TYR') THEN
      NRA = 6
      CURRINGI = TYRRINGI
      ELSE IF (THENAME(1:4).EQ.'TRP6') THEN
      NRA = 6
      CURRINGI = TRP6RINGI
      ELSE IF (THENAME(1:4).EQ.'TRP5') THEN
      NRA = 5
      CURRINGI = TRP5RINGI
      ELSE IF (THENAME(1:3).EQ.'HIS') THEN
      NRA = 5
      CURRINGI = HISRINGI
      ELSE IF (THENAME(1:4).EQ.'ADE6') THEN
      NRA = 6
      CURRINGI = ADE6RINGI
      ELSE IF (THENAME(1:4).EQ.'ADE5') THEN
      NRA = 5
      CURRINGI = ADE5RINGI
      ELSE IF (THENAME(1:4).EQ.'GUA6') THEN
      NRA = 6
      CURRINGI = GUA6RINGI
      ELSE IF (THENAME(1:4).EQ.'GUA5') THEN
      NRA = 5
      CURRINGI = GUA5RINGI
      ELSE IF (THENAME(1:3).EQ.'CYT') THEN
      NRA = 6
      CURRINGI = CYTRINGI
      ELSE IF (THENAME(1:3).EQ.'THY') THEN
      NRA = 6
      CURRINGI = THYRINGI
      ELSE IF (THENAME(1:3).EQ.'URA') THEN
      NRA = 6
      CURRINGI = URARINGI
      ELSE
      CALL DSPERR('PROTONSHIFTS',
     &'unknown ring type.  Using PHE')
      NRA = 6
      CURRINGI = PHERINGI
      END IF
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A1 = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A2 = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A3 = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A4 = SEL(1)
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A5 = SEL(1)
      IF (NRA.EQ.6) THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      IF (NSEL.GT.1) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'more than 1 atom in selection. Using first')
      END IF
      CALL MAKIND(SEL, NATOM, NSEL)
      A6 = SEL(1)
      ELSE
      A6 = 0
      END IF
      NRINGS = NRINGS + 1
      RINGATOMS(1,NRINGS) = A1
      RINGATOMS(2,NRINGS) = A2
      RINGATOMS(3,NRINGS) = A3
      RINGATOMS(4,NRINGS) = A4
      RINGATOMS(5,NRINGS) = A5
      RINGATOMS(6,NRINGS) = A6
      RINGATOMS(7,NRINGS) = NRA
      RINGCONSTS(NRINGS) = CURRINGI
C
C get the list of alpha and amide protons
C
      ELSE IF (WD(1:4).EQ.'ALPH') THEN
      CALL SELCTA(SEL, NSEL, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL, NATOM, NSEL)
      DO COUNT = 1,NSEL
      NRINGEXCEPS = NRINGEXCEPS + 1
      RINGEXCEPTIONS(NRINGEXCEPS) = SEL(COUNT)
      END DO
C
C add an observed chemical shift to the current class
C
      ELSE IF (WD(1:4).EQ.'OBSE') THEN
C
C     make sure you can't add more chemical shifts
C     than you have slots for
C
      IF (NPROTONS.EQ.MAXPROTONS) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'No room to add any more proton shift assignments.')
      CALL DSPERR('PROTONSHIFTS',
     &'Increase MAXPROTONS and recompile.')
      CALL WRNDIE(-5, 'READHSHIFT',
     &'Too many proton shift assignments.')
      END IF
C
C     make sure you don't try to add assigns to a closed-off class
C
      IF (HCLASSCLOSED(CURHSHIFTCLASS)) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'trying to add to a closed-off class')
      CALL WRNDIE(-5, 'READHSHIFT',
     &'adding to a closed-off class')
      END IF
C
C     read in the shift selection (s)
C
      CALL SELCTA(SEL1, NSEL1, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL1, NATOM, NSEL1)
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.2) THEN
      CALL SELCTA(SEL2, NSEL2, X, Y, Z,.TRUE.)
      CALL MAKIND(SEL2, NATOM, NSEL2)
      END IF
C
C     get the observed shift(s)
C
      CALL NEXTF('observed proton shift =', OBS1)
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.2) THEN
      CALL NEXTF('observed proton shift =', OBS2)
      END IF
C
C     don't do anything if there aren't any atoms selected
C
      IF ((NSEL1.LT.1).AND.(NSEL2.LT.1)) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'no atoms selected.  Ignoring entry.')
      ELSE
C
C     if there isn't a class specified,
C     start a default class
C
      IF (CURHSHIFTCLASS.EQ.0) THEN
      NHSHIFTCLASSES = 1
      CURHSHIFTCLASS = 1
      END IF
C
C     add a new observation at the end of the assignment table
C
C     Note that for format is:
C        .OBS1 -> first observed shift
C        .OBS2 -> second observed shift
C        .START1 -> position in PROTONS array of first atom
C                   selected by selection 1
C        .FINISH1 -> position in PROTONS array of last atom
C                    selected by selection 1
C        and same for .start2, .finish2
C
      NPROTASSIGNS = NPROTASSIGNS + 1
      PASSOBS1(NPROTASSIGNS) = OBS1
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.2) THEN
      PASSOBS2(NPROTASSIGNS) = OBS2
      END IF
C
C     add the positions from the first selection to the proton list
C
      PASSST1(NPROTASSIGNS)= NPROTONS+1
      DO COUNT = 1, NSEL1
      NPROTONS = NPROTONS + 1
      PROTONS(NPROTONS) = SEL1(COUNT)
      END DO
      PASSFI1(NPROTASSIGNS)= NPROTONS
C
C     add the positions from the second selection to the proton list
C
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.2) THEN
      PASSST2(NPROTASSIGNS) = NPROTONS+1
      DO COUNT = 1, NSEL2
      NPROTONS = NPROTONS + 1
      PROTONS(NPROTONS) = SEL2(COUNT)
      END DO
      PASSFI2(NPROTASSIGNS) = NPROTONS
      END IF
C
C     add this assignment to the current class
C
      HCLASSNDX(CURHSHIFTCLASS) = NPROTASSIGNS
      END IF
C
C Switch classes
C
      ELSE IF (WD(1:4).EQ.'CLAS') THEN
C
C     mark all extant classes as closed
C
      DO COUNT = 1, NHSHIFTCLASSES
      HCLASSCLOSED(COUNT) = .TRUE.
      END DO
C
C     get name and see if it's already extant
C
      CALL NEXTA4('class name =', THENAME)
      HSHIFTMODE = NEW
      DO COUNT = 1, NHSHIFTCLASSES
      IF (HSHIFTCLASSNAMES(COUNT).EQ.THENAME) THEN
      HSHIFTMODE = UPDATE
      CURHSHIFTCLASS = COUNT
      END IF
      END DO
C
C     if it's a new class, make sure there's room for it
C     and then add it
C
      IF (HSHIFTMODE.EQ.NEW) THEN
      IF (NHSHIFTCLASSES.EQ.MAXHSHIFTCLASSES) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'Too many classes.')
      CALL DSPERR('PROTONSHIFTS',
     &'Increase MAXHSHIFTCLASSES and recompile.')
      CALL WRNDIE(-5, 'READHSHIFT',
     &'Too many proton shift classes.')
      END IF
      NHSHIFTCLASSES = NHSHIFTCLASSES + 1
      CURHSHIFTCLASS = NHSHIFTCLASSES
      HSHIFTCLASSNAMES(CURHSHIFTCLASS) = THENAME
      END IF
C
C set force constant for current class,
C starting a new one if necessary
C
      ELSE IF (WD(1:4).EQ.'FORC') THEN
      IF (CURHSHIFTCLASS.EQ.0) THEN
      NHSHIFTCLASSES = 1
      CURHSHIFTCLASS = 1
      END IF
      CALL NEXTF('normal force constant =', K1)
      HSHIFTFORCES(1,CURHSHIFTCLASS) = K1
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.2) THEN
      CALL NEXTF('Constantine force constant =', K1)
      END IF
      HSHIFTFORCES(2,CURHSHIFTCLASS) = K1
      IF (HSHIFTDEGENERACY(CURHSHIFTCLASS).EQ.1) THEN
      WRITE(PUNIT, '(3(A), F8.3)')
     &'Setting force consts for class ',
     &HSHIFTCLASSNAMES(CURHSHIFTCLASS), ' to ',
     &HSHIFTFORCES(1,CURHSHIFTCLASS)
      ELSE
      WRITE(PUNIT, '(3(A), F8.3, A, F8.3)')
     &'Setting force consts for class ',
     &HSHIFTCLASSNAMES(CURHSHIFTCLASS), ' to ',
     &HSHIFTFORCES(1,CURHSHIFTCLASS), ' and ',
     &HSHIFTFORCES(2,CURHSHIFTCLASS)
      END IF
C
C set square well error for current class
C starting a new one if necessary
C
      ELSE IF (WD(1:4).EQ.'ERRO') THEN
      CALL NEXTF('error (in ppm) =', K1)
      IF (CURHSHIFTCLASS.EQ.0) THEN
      NHSHIFTCLASSES = 1
      CURHSHIFTCLASS = 1
      END IF
      HSHIFTERRORS(CURHSHIFTCLASS) = K1
      WRITE(PUNIT, '(3(A), F8.3)')
     &'Setting error for class ',
     &HSHIFTCLASSNAMES(CURHSHIFTCLASS), ' to ', K1
C
C reset shifts database
C
      ELSE IF (WD(1:4).EQ.'RESE') THEN
      CALL HSHIFTSHP
      CALL HSHIFTINIT
      CALL ALLOCHSHIFTS
C
C set potential type
C
      ELSE IF (WD(1:4).EQ.'POTE') THEN
      CALL NEXTA4('potential type =', THENAME)
      IF (THENAME.EQ.'SQUA') THEN
      WRITE(PUNIT, '(A)') 'using square well potential.'
      HSHIFTPOTENTIAL(CURHSHIFTCLASS) = SQUARE
      HSHIFTDEGENERACY(CURHSHIFTCLASS) = 1
      ELSE IF (THENAME.EQ.'HARM') THEN
      WRITE(PUNIT, '(A)') 'using harmonic potential.'
      HSHIFTPOTENTIAL(CURHSHIFTCLASS) = HARMONIC
      HSHIFTDEGENERACY(CURHSHIFTCLASS) = 1
      ELSE IF (THENAME.EQ.'MULT') THEN
      WRITE(PUNIT, '(A)') 'using mulitple shift harmonic potential.'
      HSHIFTPOTENTIAL(CURHSHIFTCLASS) = HARMONIC
      HSHIFTDEGENERACY(CURHSHIFTCLASS) = 2
      ELSE
      CALL DSPERR('PROTONSHIFTS',
     &'unknown potential. Using square well.')
      HSHIFTPOTENTIAL(CURHSHIFTCLASS) = SQUARE
      HSHIFTDEGENERACY(CURHSHIFTCLASS) = 1
      END IF
C
C print violations
C
      ELSE IF (WD(1:4).EQ.'PRIN') THEN
      CALL NEXTWD('PRINt>')
      IF (WD(1:4).NE.'THRE') THEN
      CALL DSPERR('PROTONSHIFTS',
     &'print expects THREshold parameter.')
      ELSE
      CALL NEXTF('THREshold =', PROTPRINTCUTOFF)
      IF (PROTPRINTCUTOFF.LT.ZERO) THEN
      CALL DSPERR('PROTONSHIFTS',
     &'cutoff must be positive.')
      PROTPRINTCUTOFF = ABS(PROTPRINTCUTOFF)
      END IF
      CALL NEXTA4('ALL or CLASs>', THENAME)
      IF (THENAME(1:3).EQ.'ALL') THEN
      DO COUNT = 1,NHSHIFTCLASSES
      PRINTCLASS(COUNT) = .TRUE.
      END DO
      ELSE IF (THENAME(1:4).EQ.'CLAS') THEN
      CALL NEXTA4('class name =', THENAME)
      DO COUNT = 1,NHSHIFTCLASSES
      IF (HSHIFTCLASSNAMES(COUNT).EQ.
     &THENAME) THEN
      PRINTCLASS(COUNT) = .TRUE.
      ELSE
      PRINTCLASS(COUNT) = .FALSE.
      END IF
      END DO
      ELSE
      CALL DSPERR('PROTONSHIFTS',
     &'not understood.  Printing all classes')
      DO COUNT = 1,NHSHIFTCLASSES
      PRINTCLASS(COUNT) = .TRUE.
      END DO
      END IF
      CALL NEXTA4('print to>', THENAME)
      IF (THENAME(1:4).EQ.'RMSD') THEN
      PRINTRMSD = .TRUE.
      ELSE IF (THENAME(1:4).EQ.'NORM') THEN
      PRINTRMSD = .FALSE.
      ELSE
      CALL DSPERR('PROTONSHIFTS',
     &'not understood.  Not printing to RMSD')
      PRINTRMSD = .FALSE.
      END IF
      CALL PRINTHSHIFTS
      END IF
C
C check for END statement
C
      ELSE
      CALL CHKEND('PROTONSHIFTS>', DONE)
      END IF
      END IF
      END DO
      DONE = .FALSE.
      RETURN
      END
C =============
      SUBROUTINE ALLOCHSHIFTS
C
C Allocates memory for proton shift arrays.
C
C If they're already allocated, this call is
C caused by a PROTON RESET call.
C
C by John Kuszewski Nov 1994
C =============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'mtf.inc'
C i/o
C local vbls
C begin
      IF (HSHIFTALLOCATED) THEN
           CALL FREHP(DHCDXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHCDYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHCDZPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DZPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DZPTR, IREAL8(PROTNATOM))
           CALL FREHP(AAPTR, INTEG4(3*PROTNATOM))
           CALL FREHP(APPTR, INTEG4(PROTNATOM))
           CALL FREHP(ACENPTR, IREAL8(PROTNATOM))
           CALL FREHP(AVFPTR, IREAL8(PROTNATOM))
           CALL FREHP(ACPTR, IREAL8(2*PROTNATOM))
           CALL FREHP(EEPTR, INTEG4(PROTNATOM))
           CALL FREHP(EHPTR, INTEG4(PROTNATOM))
           CALL FREHP(ECPTR, IREAL8(PROTNATOM))
           CALL FREHP(RAPTR, INTEG4(7*PROTNATOM))
           CALL FREHP(REPTR, INTEG4(PROTNATOM))
           CALL FREHP(HRCPTR, IREAL8(PROTNATOM))
           CALL FREHP(PROTPTR, INTEG4(PROTNATOM))
           CALL FREHP(DHOBS1, IREAL8(MAXPROTONS))
           CALL FREHP(DHOBS2, IREAL8(MAXPROTONS))
           CALL FREHP(DHST1, INTEG4(MAXPROTONS))
           CALL FREHP(DHST2, INTEG4(MAXPROTONS))
           CALL FREHP(DHFI1, INTEG4(MAXPROTONS))
           CALL FREHP(DHFI2, INTEG4(MAXPROTONS))
           CALL FREHP(RCPTR, IREAL8(PROTNATOM))
      END IF
C
      DHCDXPTR = 0
      DHCDYPTR = 0
      DHCDZPTR = 0
      DHC1DXPTR = 0
      DHC1DYPTR = 0
      DHC1DZPTR = 0
      DHC2DXPTR = 0
      DHC2DYPTR = 0
      DHC2DZPTR = 0
      AAPTR = 0
      APPTR = 0
      ACENPTR = 0
      AVFPTR = 0
      ACPTR = 0
      EEPTR = 0
      EHPTR = 0
      ECPTR = 0
      RAPTR = 0
      REPTR = 0
      HRCPTR = 0
      PROTPTR = 0
      DHOBS1 = 0
      DHOBS2 = 0
      DHST1 = 0
      DHST2 = 0
      DHFI1 = 0
      DHFI2 = 0
      RCPTR = 0
      HSHIFTALLOCATED = .FALSE.
C
      PROTNATOM=NATOM
      IF (PROTNATOM.NE.0) THEN
        DHCDXPTR = ALLHP(IREAL8(PROTNATOM))
        DHCDYPTR = ALLHP(IREAL8(PROTNATOM))
        DHCDZPTR = ALLHP(IREAL8(PROTNATOM))
        DHC1DXPTR = ALLHP(IREAL8(PROTNATOM))
        DHC1DYPTR = ALLHP(IREAL8(PROTNATOM))
        DHC1DZPTR = ALLHP(IREAL8(PROTNATOM))
        DHC2DXPTR = ALLHP(IREAL8(PROTNATOM))
        DHC2DYPTR = ALLHP(IREAL8(PROTNATOM))
        DHC2DZPTR = ALLHP(IREAL8(PROTNATOM))
        AAPTR = ALLHP(INTEG4(3*PROTNATOM))
        APPTR = ALLHP(INTEG4(PROTNATOM))
        ACENPTR = ALLHP(IREAL8(PROTNATOM))
        AVFPTR = ALLHP(IREAL8(PROTNATOM))
        ACPTR = ALLHP(IREAL8(2*PROTNATOM))
        EEPTR = ALLHP(INTEG4(PROTNATOM))
        EHPTR = ALLHP(INTEG4(PROTNATOM))
        ECPTR = ALLHP(IREAL8(PROTNATOM))
        RAPTR = ALLHP(INTEG4(7*PROTNATOM))
        REPTR = ALLHP(INTEG4(PROTNATOM))
        HRCPTR = ALLHP(IREAL8(PROTNATOM))
        PROTPTR = ALLHP(INTEG4(PROTNATOM))
        DHOBS1 = ALLHP(IREAL8(MAXPROTONS))
        DHOBS2 = ALLHP(IREAL8(MAXPROTONS))
        DHST1 = ALLHP(INTEG4(MAXPROTONS))
        DHST2 = ALLHP(INTEG4(MAXPROTONS))
        DHFI1 = ALLHP(INTEG4(MAXPROTONS))
        DHFI2 = ALLHP(INTEG4(MAXPROTONS))
        RCPTR = ALLHP(IREAL8(PROTNATOM))
C
C mark that I've allocated space
C
        HSHIFTALLOCATED = .TRUE.
      END IF
C
      RETURN
      END
C==============
      SUBROUTINE HSHIFTSHP
C
C deallocates carbon shifts stuff
C
C by Alexandre Bonvin 4/8/96
C =============
      IMPLICIT NONE
C include files
      INCLUDE 'cns.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'hshifts.inc'
      INCLUDE 'mtf.inc'
C i/o
C local vbls
C begin
      IF (HSHIFTALLOCATED) THEN
           CALL FREHP(DHCDXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHCDYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHCDZPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC1DZPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DXPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DYPTR, IREAL8(PROTNATOM))
           CALL FREHP(DHC2DZPTR, IREAL8(PROTNATOM))
           CALL FREHP(AAPTR, INTEG4(3*PROTNATOM))
           CALL FREHP(APPTR, INTEG4(PROTNATOM))
           CALL FREHP(ACENPTR, IREAL8(PROTNATOM))
           CALL FREHP(AVFPTR, IREAL8(PROTNATOM))
           CALL FREHP(ACPTR, IREAL8(2*PROTNATOM))
           CALL FREHP(EEPTR, INTEG4(PROTNATOM))
           CALL FREHP(EHPTR, INTEG4(PROTNATOM))
           CALL FREHP(ECPTR, IREAL8(PROTNATOM))
           CALL FREHP(RAPTR, INTEG4(7*PROTNATOM))
           CALL FREHP(REPTR, INTEG4(PROTNATOM))
           CALL FREHP(HRCPTR, IREAL8(PROTNATOM))
           CALL FREHP(PROTPTR, INTEG4(PROTNATOM))
           CALL FREHP(DHOBS1, IREAL8(MAXPROTONS))
           CALL FREHP(DHOBS2, IREAL8(MAXPROTONS))
           CALL FREHP(DHST1, INTEG4(MAXPROTONS))
           CALL FREHP(DHST2, INTEG4(MAXPROTONS))
           CALL FREHP(DHFI1, INTEG4(MAXPROTONS))
           CALL FREHP(DHFI2, INTEG4(MAXPROTONS))
           CALL FREHP(RCPTR, IREAL8(PROTNATOM))
      END IF
      RETURN
      END
C ==============
      SUBROUTINE HSHIFTINIT
C
C Initialization stuff for proton shifts
C
C by John Kuszewski Nov 1994
C ==============
      IMPLICIT NONE
C include files
      INCLUDE 'hshifts.inc'
C i/o
C local vbls
      INTEGER COUNT
C begin
      HSHIFTMOVECUTOFF = 0.0
      HSHIFTALLOCATED = .FALSE.
      HSHIFTMODE = NEW
      NPROTONS = 0
      NPROTASSIGNS = 0
      NANIS = 0
      NHEAVIES = 0
      NRINGS = 0
      NAMIDES = 0
      NRINGEXCEPS = 0
      NHSHIFTCLASSES = 0
      CURHSHIFTCLASS = 0
      DO COUNT = 1, MAXHSHIFTCLASSES
           HSHIFTCLASSNAMES(COUNT) = 'DEFAULT'
           HCLASSNDX(COUNT) = 0
           HSHIFTFORCES(1, COUNT) = 0.1
           HSHIFTFORCES(2, COUNT) = 0.1
           HSHIFTERRORS(COUNT) = 0.3
           HCLASSCLOSED(COUNT) = .FALSE.
           HSHIFTDEGENERACY(COUNT) = 1
           HSHIFTPOTENTIAL(COUNT) = HARMONIC
      END DO
      RETURN
      END
C ===============
      SUBROUTINE PRINTHSHIFTS
C
C Since the number of things to be printed is
C relatively large, the printing is handled by
C EPROTONSHIFT2.  Thus, just passes the work to
C EPROTONSHIFT.
C
C The cutoff value is stored in hshifts.inc for
C convenience.
C
C by John Kuszewski Dec 1994
C ===============
      IMPLICIT NONE
C include files
C i/o
C local variables
      DOUBLE PRECISION DUMMY
C begin
      CALL EPROTONSHIFT(DUMMY, 'ANALYZE')
      RETURN
      END
C
      SUBROUTINE SCRHSFT
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'hshifts.inc'
C
C reset proton shift database
      IF (HSHIFTALLOCATED) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: proton chemical shift database erased.'
      CALL HSHIFTSHP
      CALL HSHIFTINIT
      END IF
C
      RETURN
      END
C
