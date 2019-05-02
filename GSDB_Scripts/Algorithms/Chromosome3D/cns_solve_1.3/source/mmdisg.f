      SUBROUTINE MMDG
C
C Reference for this implementation of metric matrix distance
C geometry:
C
C  J. Kuszewski, M. Nilges, A. T. Brunger, "Sampling and Efficiency
C   of Metric Matrix Distance Geometry: A Novel Partial Metrization
C   Algorithm",
C   J. Biomolecular NMR, in press, 1992.
C
C The MMDG (Metric Matrix Distance Geometry) routine
C calls the distance geometry routines.  See
C the individual routines for explanations of the
C non-obvious variables.
C
C BACC, TACC, PACC, IACC, & RACC are values which are added
C and subtracted to the bond lengths, angles, rigid groups, etc.
C to allow for flexibility in the metrization.  The value
C for each should be given in the native unit (Angstroms for
C length, degrees for angles) and is converted to the
C internal unit (Angstroms and radians).
C
C by John Kuszewski 9/20/90
C ========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'update.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'param.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'seed.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'mmdisg.inc'
      INCLUDE 'timer.inc'
C local variables
      INTEGER VALPTR, NSLCT, NMETR, BUNIT, LENGTH
      INTEGER NGROUP, NSELE, NSUB, MFLAGS
      LOGICAL DODUMP, DORND, QEMBED, OK, QSTOREB, QWRITEB, QRECALL
      LOGICAL QREAD, NOSMOOTH
      DOUBLE PRECISION BACC, TACC, PACC, IMACC, RACC, OLDSEED
      DOUBLE PRECISION XPECTD
      DOUBLE PRECISION SECS, SECS2
      CHARACTER*4 SREFER, SSHORT
      INTEGER NPAIR
C pointer
      INTEGER RNDPTR, PNTPTR, SUBPTR
      INTEGER LSTPTR, ACCPTR, SLCTPTR
      INTEGER DUBPTRS, TDUB, VECS
      INTEGER IPK, JPK, IBLOEX, INBEX
C parameters and dummies
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION DBREAL
C
C begin
C
C allocate heap space for the arrays
C
C DUB:  distance upper bounds (N-by-N double), also used for
C       picked distances. Convention: for i<j DUB(i,j) is the
C       upper bound, otherwise it is the lower bound.
C DUBS: DUB matrix projected into substructure space
C VAL:  eigenvalues of MM (N-by-1 double)
C SLCT: atoms selected to be considered for use as root during
C        re-tightening phase
C PNT:  Pointer array to show beginnings of rigid groups in LST
C       (N+1-by-1 integer)
C LST:  Array of atoms in rigid groups (mflags-by-1 INTEGER)
C ACC:  Array of accuracy values for each rigid group (N-by-1 double)
C RND:  N-by-1 array of double-precision (random) values used to decide
C       the metrization order during random metrization.
C SUB:  selection array for substructures (ATB)
      MFLAGS=10*NATOM
      PNTPTR=ALLHP(INTEG4(NATOM+1))
      LSTPTR=ALLHP(INTEG4(MFLAGS))
      ACCPTR=ALLHP(IREAL8(NATOM))
      VALPTR=ALLHP(IREAL8(NATOM))
      SLCTPTR=ALLHP(INTEG4(NATOM))
      RNDPTR=ALLHP(IREAL8(NATOM))
      SUBPTR=ALLHP(INTEG4(NATOM))
C
C default reference are parameters
      SREFER='PARA'
C
C these are the default values for accuracy
      BACC = 0.01D0
      TACC = 2.0D0
      PACC = 2.0D0
      IMACC = 2.0D0
      RACC = 0.01D0
C no rigid groups by default
      NGROUP = 0
C random metrization by default
      DORND = .TRUE.
      DODUMP = .FALSE.
C zero atoms in metrization stage by default
      NMETR = 0
      NSLCT = 0
      CALL FILL4(HEAP(SLCTPTR), NATOM, 0)
C all atoms in substructure by default (ATB)
      NSUB=NATOM
      CALL FILL4(HEAP(SUBPTR), NSUB, 1)
      CALL MAKIND(HEAP(SUBPTR), NATOM, NSUB)
C
C use automatic mode for shortest path bound smoothing by default
      SSHORT='AUTO'
C
C by default don't store or write the bounds matrices
      QWRITEB=.FALSE.
      QSTOREB=.FALSE.
C
C by default don't read or recall bounds matrices
      QRECALL=.FALSE.
      QREAD=.FALSE.
C
C parse the command line, taking care of stuff as you go
C
      CALL PUSEND('MMDG>')
      DO WHILE (.NOT. DONE)
      CALL NEXTWD('MMDG>')
      CALL MISCOM('MMDG>', USED)
      IF (.NOT. USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
         CALL CNSHELP('cns-mmdg')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)')
     &' -------------------------------------------------------------',
     &'--------------'
      IF (.NOT.QRECALL.AND..NOT.QREAD) THEN
      IF (SREFER.EQ.'PARA') THEN
      WRITE(6,'(A)') ' REFErence=PARAmeter'
      ELSE
      WRITE(6,'(A)') ' REFErence=COORdinate'
      END IF
      WRITE(6,'(A,A)') ' SHORtest-path-algorithm=',SSHORT
      WRITE(6,'(A,F6.4,A,F6.4,A,F6.4,A,F6.4)')
     &' BACCuracy=',BACC,' TACCuracy=',TACC,' PACCuracy=',PACC,
     &' IACCuracy=',IMACC
      WRITE(6,'(A,I5)')
     & ' Number of rigid groups=',NGROUP
      ELSE IF (QRECALL) THEN
      WRITE(6,'(A)') ' Bounds matrices are recalled from memory.'
      ELSE IF (QREAD) THEN
      WRITE(6,'(A)') ' Bounds matrices will be read from file.'
      END IF
      IF (QWRITEB) THEN
      WRITE(6,'(A)')
     &  ' Smoothed bounds matrices will be written to file.'
      ELSE IF (QSTOREB) THEN
      WRITE(6,'(A)')
     &  ' Smoothed bounds matrices will be stored in memory.'
      ELSE
      IF (NMETR.GT.0) THEN
      IF (DORND) THEN
      WRITE(6,'(A,I6,A)')
     & ' RANDom metrization being used for ',NMETR,' atoms.'
      ELSE
      WRITE(6,'(A,I6,A)')
     & ' ORDEred metrization being used for ',NMETR,' atoms.'
      END IF
      END IF
      WRITE(6,'(A,I6,A)')
     & ' Bound smoothing being used for ',NATOM-NMETR,' atoms.'
      IF (NSUB.LT.NATOM) THEN
      WRITE(6,'(A,I6,A)')
     & ' A substructure consisting of ',NSUB,' atoms is embedded.'
      ELSE
      WRITE(6,'(A)') ' All atoms are embedded.'
      END IF
      END IF
      WRITE(6,'(A,/,A,F6.4,A,I4,/,A,I6)')
     &' DG-energy term:',
     &'     SCALe=',DGSCALE,' EXPOnent=',DGEXPO,
     &'     Number of selected atoms=',NDGSELE
      WRITE(6,'(2A)')
     &' -------------------------------------------------------------',
     &'--------------'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'REFE') THEN
      CALL NEXTA4('REFErence=',SREFER)
      IF (SREFER.NE.'PARA'.AND.SREFER.NE.'COOR') THEN
      SREFER='PARA'
      CALL DSPERR('MMDG','unkown REFErence. Set to PARAmeter.')
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'SHOR') THEN
      CALL NEXTA4('SHORtest-path-algorithm=',SSHORT)
      IF (SSHORT.NE.'AUTO'.AND.SSHORT.NE.'SPAR'.AND.SSHORT.NE.'FULL')
     & THEN
      SSHORT='AUTO'
      CALL DSPERR('MMDG','unkown SHORtest-path-algorithm.Set to AUTO.')
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'VERB') THEN
      CALL NEXTLO('VERBose mode = ', DODUMP)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'BACC') THEN
      CALL NEXTF('bond length accurracy=', BACC)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'TACC') THEN
      CALL NEXTF('angle accurracy=', TACC)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'PACC') THEN
      CALL NEXTF('dihedral angles accurracy=', PACC)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'IACC') THEN
      CALL NEXTF('improper angles accurracy=', IMACC)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'GROU') THEN
      CALL SELCTA(HEAP(SLCTPTR), NSLCT, X, Y, Z, .TRUE.)
      CALL MAKIND(HEAP(SLCTPTR), NATOM, NSLCT)
      CALL NEXTF('accuracy for this rigid group= ', RACC)
      CALL DGGROUP (HEAP(SLCTPTR), MFLAGS, HEAP(LSTPTR),
     &              HEAP(PNTPTR), HEAP(ACCPTR), RACC, NGROUP, NSLCT)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'RAND') THEN
      DORND = .TRUE.
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'ORDE') THEN
      DORND = .FALSE.
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'METR') THEN
      CALL SELCTA(HEAP(SLCTPTR), NSLCT, X, Y, Z, .TRUE.)
      CALL NEXTI('numb. of sel. atoms for re-tight. phase =', NMETR)
      NMETR=MIN(NSLCT,NMETR)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'SUBS') THEN
      CALL SELCTA(HEAP(SUBPTR), NSUB, X, Y, Z,.TRUE.)
      CALL MAKIND(HEAP(SUBPTR), NATOM, NSUB)
C----------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'STOREB') THEN
C
C we allow QREAD here!
      IF (QRECALL) THEN
      CALL DSPERR('MMDG',
     & 'this option is incompatible with RECALLbounds or READBOunds.')
      ELSE
      QSTOREB=.TRUE.
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'RECALL') THEN
      IF (QWRITEB.OR.QSTOREB) THEN
      CALL DSPERR('MMDG',
     & 'this option is incompatible with RECALLbounds or READBOunds.')
      ELSE
      IF (DUBLEN.NE.NATOM**2) THEN
      CALL WRNDIE(-5,'MMDG',
     & 'Dimen. of stored bounds matrix doesnt match numb. of atoms.')
      ELSE
      QRECALL=.TRUE.
      END IF
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'READBO') THEN
C
C we allow QSTOREB here!
      IF (QWRITEB) THEN
      CALL DSPERR('MMDG',
     & 'this option is incompatible with STOREBounds or WRITEBounds.')
      ELSE
      CALL NEXTFI('READBOunds-file=',OFILE)
      CALL ASSFIL(OFILE,BUNIT,'READ','UNFORMATTED',ERROR)
      IF (.NOT.ERROR) THEN
      REWIND(UNIT=BUNIT)
      READ(BUNIT) LENGTH
      IF (LENGTH.NE.NATOM**2) THEN
      CALL WRNDIE(-5,'MMDG',
     & 'Dimension of bounds matrix does not match number of atoms.')
      ELSE
      CALL DGRESZ(LENGTH)
      CALL DGREAD(BUNIT,NATOM,HEAP(DUBPTR),DGQVDW)
      QREAD=.TRUE.
      END IF
      END IF
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:6).EQ.'WRITEB') THEN
      IF (QREAD.OR.QRECALL) THEN
      CALL DSPERR('MMDG',
     & 'this option is incompatible with STOREBounds or WRITEBounds.')
      ELSE
      QWRITEB=.TRUE.
      CALL NEXTFI('WRITEBounds-file=',IFILE)
      CALL ASSFIL(IFILE,BUNIT,'WRITE','UNFORMATTED',ERROR)
      IF (ERROR) BUNIT=0
      END IF
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'SELE') THEN
C
C free-up old space if required
      IF (DGSLEN.GT.0) THEN
      CALL FREHP(DGSELE,INTEG4(DGSLEN))
      END IF
      DGSELE=ALLHP(INTEG4(NATOM))
      DGSLEN=NATOM
      CALL FILL4(HEAP(DGSELE),NATOM,0)
      NDGSELE=0
      CALL SELCTA(HEAP(DGSELE),NDGSELE,X,Y,Z,.TRUE.)
      CALL MAKIND(HEAP(DGSELE),NATOM,NDGSELE)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
      CALL NEXTF('SCALe=',DGSCALE)
C----------------------------------------------------------------------
      ELSE IF (WD(1:4).EQ.'EXPO') THEN
      CALL NEXTI('EXPOnent=',DGEXPO)
C----------------------------------------------------------------------
      ELSE
      CALL CHKEND('DISGEO>', DONE)
      END IF
      END IF
      END DO
      DONE = .FALSE.
C
      IF (NSUB.GT.0) THEN
C
C store the current random number seed.  This is done to avoid
C dependencies in the test cases upon changing around sorting routines.
C For example, the nonbonded exclusion list routine makes use of the
C sorting routines which in turn would change the random number seed.
      OLDSEED=SEED
C
C bypass the following bounds and bounds smoothing section if
C QRECALL is specified
      IF (.NOT.QRECALL.AND..NOT.QREAD) THEN
      IF (TIMER.GT.0) CALL VCPU(SECS)
      CALL DGRESZ(NATOM**2)
C
C translate degrees into radians for accuracy
      CALL DG2RAD(TACC, TACC)
      CALL DG2RAD(PACC, PACC)
      CALL DG2RAD(IMACC, IMACC)
C
C make sure accuracy values are positive
      BACC = ABS(BACC)
      TACC = ABS(TACC)
      PACC = ABS(PACC)
      IMACC = ABS(IMACC)
C
C initialize the upper and lower bounds matrices
      CALL DGSETUL(NATOM, HEAP(DUBPTR))
C
C The sequence of putting distance information into the bounds
C matrices depends to some degree on the order.  This applies
C to angles, dihedrals, impropers, and dihedral angle restraints as the
C definition of the 1-3 and 1-4 distances depends on the previously set
C bonds (and angles).  Van der Waals repulsion is the penultimate step.
C A repulsive lower bound is only applied if that does not violate the
C previously determined upper bound.  The RIGID group information
C overwrites all previous information for the selected atoms. (ATB)
C
C include conformational energy terms information for bonded terms
      IF (QENER(SSBOND).OR.QENER(SSANGL).OR.QENER(SSDIHE).OR.
     & QENER(SSIMPR)) THEN
C
C loop over all active conformational energy PIGs.
      DO NSELE=1,NPIG
      CALL DGSETEN(NATOM,HEAP(DUBPTR),BACC,TACC,PACC,IMACC,SREFER )
      END DO
      END IF
C
C include NOE information
      IF (QENER(SSNOE)) THEN
      CALL DGSETNO(NATOM, HEAP(DUBPTR),
     &        HEAP(HPNDIS), HEAP(HPNLOW), HEAP(HPNHIG), HEAP(HPNORR),
     &        HEAP(HPNILS), HEAP(HPNJLS), HEAP(HPNIPR), HEAP(HPNJPR),
     &        HEAP(HPNCND), HEAP(HPNCV), DODUMP)
      END IF
C
C include dihedral angle restraints information
      IF (QENER(SSCDIH)) THEN
      CALL DGSETDH(NATOM, HEAP(DUBPTR), HEAP(ICIP),
     &   HEAP(ICJP), HEAP(ICKP), HEAP(ICLP), HEAP(ICCPB), HEAP(ICCPO))
      END IF
C
C include rigid group information
C (this overwrites all previous info) (ATB)
      CALL DGSETRG(NATOM, HEAP(DUBPTR), HEAP(PNTPTR),
     &             HEAP(LSTPTR), HEAP(ACCPTR), NGROUP)
C
C set the flag that indicates that VDW interactions are
C included on-the-fly in mmdg calculations
C (used only for DG energy term)
      IF (QENER(SSVDW)) THEN
      DGQVDW=-1
      ELSE
      DGQVDW=0
      END IF
C
      IF (DODUMP) WRITE(6,'(A)') 'finished initialization'
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  initial setup of bounds matrix=',SECS2-SECS
      END IF
C
C if substructures are present or the QSTOREB or QWRITEB
C options are set we have to perform boundsmoothing
      IF (NSUB.NE.NATOM.OR.QSTOREB.OR.QWRITEB) THEN
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      END IF
      CALL DGSMOOT(SSHORT, NATOM, HEAP(DUBPTR))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  complete bound smoothing=',SECS2-SECS
      END IF
      END IF
C
      END IF
C
C
C we only continue if QSTOREB or QWRITEB are not set
      IF (.NOT.(QSTOREB.OR.QWRITEB)) THEN
C
C if the substructure option is used we have to allocate
C space for the reduced matrices (ATB)
      IF (NSUB.NE.NATOM) THEN
      DUBPTRS=ALLHP(IREAL4(NSUB**2))
      ELSE
C
C otherwise we identify the "reduced" and full bounds matrices
      DUBPTRS=DUBPTR
      END IF
C
C if substructures are present we have to project the full bounds
C matrices into the substructure space.  We also project the
C metrization selection into the substructure space. Note, that
C the latter operation is irreversible!!
      IF (NSUB.NE.NATOM) THEN
      CALL DGREDUC(NATOM,NSUB,HEAP(DUBPTR),
     &             HEAP(DUBPTRS),
     &             HEAP(SUBPTR),NSLCT,HEAP(SLCTPTR))
      END IF
C
      NMETR=MIN(NSLCT,NMETR)
C
      IF (NMETR.GT.0.AND.DORND) THEN
      WRITE(6,'(A,I6,A)')
     & ' MMDG: ',NMETR,' atoms will be randomly metrized.'
      ELSE IF (NMETR.GT.0.AND..NOT.DORND) THEN
      WRITE(6,'(A,I6,A)')
     & ' MMDG: ',NMETR,' atoms will be metrized in order.'
      ELSE
      WRITE(6,'(A)') ' MMDG: no atoms will be metrized.'
      END IF
C
C make an index list for the metrization selection
      CALL MAKIND(HEAP(SLCTPTR), NSUB, NSLCT)
C
C
C now the metrization and bound-smoothing
      IF (DODUMP) WRITE (6, '(A)') 'starting metrization'
C
C generate a 1-2 and 1-3 nonbonded (symmetric) exclusion list
      IF (QENER(SSVDW)) THEN
C make sure lookup tables for van der Waals parameters are updated
      IF (UPNBLK) THEN
      CALL NBUPDA
      UPNBLK=.FALSE.
      END IF
C
C first we obtain the required number of exclusions without
C actually storing them
      CALL DGSETNX(NPAIR,0,0,0,0,.TRUE.)
C
C now we allocate the required space and then call the routine
C again to store the exclusions in IBLOEX, INBEX
      IPK=ALLHP(INTEG4(NPAIR))
      JPK=ALLHP(INTEG4(NPAIR))
      IBLOEX=ALLHP(INTEG4(NATOM+1))
      INBEX=ALLHP(INTEG4(NPAIR))
      CALL DGSETNX(NPAIR,HEAP(IPK),HEAP(JPK),HEAP(IBLOEX),HEAP(INBEX),
     &             .FALSE.)
      CALL FREHP(IPK,INTEG4(NPAIR))
      CALL FREHP(JPK,INTEG4(NPAIR))
      ELSE
      IBLOEX=1
      INBEX=1
      END IF
C
C re-instate the old seed for the random number generator
      SEED=OLDSEED
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      END IF
      CALL DGMETRI(DORND, DODUMP, OK, HEAP(SLCTPTR), NSLCT, NMETR,
     &       HEAP(RNDPTR), NSUB, HEAP(SUBPTR), HEAP(DUBPTRS),
     &       HEAP(IBLOEX), HEAP(INBEX))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  metrization=',SECS2-SECS
      END IF
      IF (DODUMP) WRITE (6, '(A)') 'filling in unset distances'
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      END IF
C
C test to see if we have to perform bound-smoothing:
C we don't want to do it if bound-smoothing was already
C performed and no metrization was carried out,
C or complete metrization was performed (MN, 2/13/94).
      NOSMOOTH=
     & ((NSUB.NE.NATOM.OR.QRECALL.OR.QREAD).AND.NMETR.EQ.0)
     & .OR.NMETR.EQ.NATOM
      CALL DGNOMET (DODUMP, OK, HEAP(SLCTPTR), NSLCT,
     &       NSUB, HEAP(SUBPTR), HEAP(DUBPTRS),
     &       HEAP(IBLOEX), HEAP(INBEX), SSHORT, NOSMOOTH)
C
C free-up space for the nonbonded exclusion list
      IF (QENER(SSVDW)) THEN
      CALL FREHP(IBLOEX,INTEG4(NATOM+1))
      CALL FREHP(INBEX,INTEG4(NPAIR))
      END IF
C
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  final bound smoothing=',SECS2-SECS
      END IF
C
      IF (DODUMP) WRITE (6, '(A)') 'embedding coordinates'
C
C create metric matrix and store in triangular matrix TDUB
C (VALPTR is used for working space here).
      TDUB=ALLHP(IREAL8(NSUB*(NSUB+1)/2))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      END IF
      CALL DGMMMM(NSUB,HEAP(DUBPTRS),HEAP(TDUB),HEAP(VALPTR),XPECTD)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  computation of metric matrix=',SECS2-SECS
      END IF
C
C we leave the freeing of space until after embedding
C to be able to add unembedded points (MN)
C
C now we're ready for embedding.
      VECS=ALLHP(IREAL8(NSUB**2))
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS)
      END IF
      CALL DGEMBED(NATOM, HEAP(DUBPTR), HEAP(TDUB), HEAP(VECS),
     &    HEAP(VALPTR), NSUB, HEAP(SUBPTR), DODUMP, QEMBED, XPECTD)
      IF (TIMER.GT.0) THEN
      CALL VCPU(SECS2)
      WRITE(6,'(A,F10.4)')
     & ' MMDG: CPU-time:  embedding=',SECS2-SECS
      END IF
      CALL FREHP(VECS,IREAL8(NSUB**2))
      CALL FREHP(TDUB,IREAL8(NSUB*(NSUB+1)/2))
      IF (QEMBED) THEN
      CALL DECLAR('EMBEDDED', 'LO', 'TRUE', DBCOMP, DBREAL)
      ELSE
      CALL DECLAR('EMBEDDED', 'LO', 'FALSE', DBCOMP, DBREAL)
      END IF
C
C free up space for the bounds matrices (ATB)
      IF (NSUB.NE.NATOM) THEN
      IF (.NOT.QRECALL) CALL DGRESZ(0)
      CALL FREHP(DUBPTRS,IREAL4(NSUB**2))
      ELSE
      IF (.NOT.QRECALL) CALL DGRESZ(0)
      END IF
C
      ELSE IF (QSTOREB) THEN
      WRITE(6,'(A)')
     & ' MMDG: storing smoothed bounds matrix in memory. ',
     &'        No metrization or embedding performed. '
C
      ELSE IF (QWRITEB.AND.BUNIT.GT.0) THEN
      WRITE(6,'(A)')
     & ' MMDG: writing smoothed bounds matrices to binary file. ',
     &'        No metrization or embedding performed. '
      WRITE(BUNIT) DUBLEN
      CALL DGWRIT(BUNIT,NATOM,HEAP(DUBPTR),DGQVDW)
      CALL VCLOSE(BUNIT,'KEEP',ERROR)
      CALL DGRESZ(0)
      END IF
C
      END IF
C
C
C deallocate other heap space and quit
C
      CALL FREHP(SUBPTR, INTEG4(NATOM))
      CALL FREHP(PNTPTR, INTEG4(NATOM+1))
      CALL FREHP(LSTPTR, INTEG4(MFLAGS))
      CALL FREHP(ACCPTR, IREAL8(NATOM))
      CALL FREHP(VALPTR, IREAL8(NATOM))
      CALL FREHP(SLCTPTR, INTEG4(NATOM))
      CALL FREHP(RNDPTR, IREAL8(NATOM))
      RETURN
      END
C ======================
      SUBROUTINE DGSETUL(NSUB, DUB)
C
C Initialize all lower bounds to a small number and
C all upper bounds to a big number
C except when it's a bound from a point to itself
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
C local
      INTEGER COUNT, COUNT2
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      DO COUNT=1, NSUB
      DO COUNT2=COUNT, NSUB
      IF (COUNT .EQ. COUNT2) THEN
      DUB(COUNT,COUNT) = ZERO
      ELSE
      DUB(COUNT,COUNT2) = R4BIG
      DUB(COUNT2,COUNT) = ZERO
      END IF
      END DO
      END DO
      RETURN
      END
C ======================
      SUBROUTINE DGSETRG(NSUB, DUB, POINTR, LST, ACCS, NGROUP)
C
C Put rigid group information into DUB matrix.
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      INTEGER POINTR(*), LST(*), NGROUP
      DOUBLE PRECISION  ACCS(*)
C local
      INTEGER GRPNUM, PH, PJ, J, H, PHH, PJJ
      LOGICAL ERR
      DOUBLE PRECISION RACC, DX, DY, DZ, DIST, HILEN, LOLEN, R
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C Calculate distance constraints from rigid bodies and add them to DUB
C
      ERR=.FALSE.
      DO GRPNUM=1,NGROUP
      RACC = ACCS(GRPNUM)
      RACC = ABS(RACC)
C
C do a triangular double loop over all atoms of the group
C exlude self-terms (ATB)
      DO H=POINTR(GRPNUM), POINTR(GRPNUM+1) - 1
      PH=LST(H)
      IF (.NOT.INITIA(PH,X,Y,Z)) THEN
      WRITE(6,'(9A)')
     &' %DGSETRG-ERR: unknown reference coordinates for atom "',
     &SEGID(PH),'-',RESID(PH),'-',RES(PH),'-',TYPE(PH),'"'
      ERR=.TRUE.
      END IF
      DO J=H+1, POINTR(GRPNUM+1) - 1
      PJ=LST(J)
      DX=X(PH)-X(PJ)
      DY=Y(PH)-Y(PJ)
      DZ=Z(PH)-Z(PJ)
      DIST=SQRT(DX**2 + DY**2 + DZ**2)
C
C overwrite all previous info (ATB)
      HILEN=DIST + RACC
      LOLEN=MAX(ZERO,DIST - RACC)
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      R=DUB(PHH,PJJ)
      HILEN=MIN(R,HILEN)
      R=DUB(PJJ,PHH)
      LOLEN=MAX(R,LOLEN)
      DUB(PHH,PJJ)=HILEN
      DUB(PJJ,PHH)=LOLEN
      END DO
      END DO
      END DO
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'DGSETRG','Unknown main coordinate set')
      END IF
C
      RETURN
      END
C ======================
      SUBROUTINE DGSETNO(NSUB,DUB,DIS,LOW,HI,ORR,ILS,JLS,IPR,JPR,
     &                   NOECND,NOECV,DODUMP)
C
C Put NOE information into DUB matrix.
C Automatic pseudoatom correction (ATB).
C Included NOE offset.
C Ignores ambiguous NOEs (MN 20-Feb-1998)
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION DIS(*), LOW(*), HI(*)
      INTEGER ORR(*), ILS(*), JLS(*), IPR(*), JPR(*)
      INTEGER NOECND(*), NOECV(*)
      LOGICAL DODUMP
C local
      INTEGER COUNT, PH, PJ, ICOUNT, JCOUNT, PHH, PJJ
      INTEGER I, II, J, JJ, LI, LJ, NNORR
      DOUBLE PRECISION HILEN, LOLEN, R, HIL, LOL
      DOUBLE PRECISION XI, YI, ZI, XJ, YJ, ZJ, XIJ, YIJ, ZIJ, RI, RJ
      LOGICAL CORRECT, WRN
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
C
C This reads the NOE data into the DUB array.  As noted in noe.inc,
C the lower bound is the distance (HPNDIS) - the lower "error" (HPNLOW).
C The upper bound is the distance + the upper "error" (HPNHIG)
C minus offset (which is normally zero).
C
      WRN=.FALSE.
      DO COUNT=1, NOENUM
      IF (NOECV(COUNT).NE.NOEICV) THEN
      HIL=DIS(COUNT) + HI(COUNT)-NOEOFF(NOECND(COUNT))
      LOL=DIS(COUNT) - LOW(COUNT)
C
C test if multiple atom selections are defined that
C would require correction of the upper and lower bounds.
      NNORR=ORR(COUNT+1)-ORR(COUNT)
      IF (NNORR.GT.1) THEN
      IF (DODUMP) THEN
      WRITE(6,'(A,I6)')
     & ' ambiguous restraint ignored at NOE #',COUNT
      END IF
      ELSE
      II=IPR(ORR(COUNT)+1)+1-IPR(ORR(COUNT)+1+1)
      JJ=JPR(ORR(COUNT)+1)+1-JPR(ORR(COUNT)+1+1)
      CORRECT=((II.NE.0.OR.JJ.NE.0).AND..NOT.WRN)
C if more than one pair of atoms in restraint
      IF (CORRECT) THEN
      IF (DODUMP) THEN
      WRITE(6,'(A,I6)')
     & ' Pseudoatom correction applied to NOE #',COUNT
      END IF
      LI=0
      XI=ZERO
      YI=ZERO
      ZI=ZERO
C compute geometric center for I-atoms
      DO II=IPR(ORR(COUNT)+1)+1,IPR(ORR(COUNT)+1+1)
      I=ILS(II)
      IF (.NOT.INITIA(I,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETNO-WRN: unknown reference coordinates for atom "',
     &SEGID(I),'-',RESID(I),'-',RES(I),'-',TYPE(I),'"'
      WRITE(6,'(2A)')
     &'               Automatic pseudoatom correction',
     &' wont be applied.'
      WRN=.TRUE.
      END IF
      XI=XI+REFX(I)
      YI=YI+REFY(I)
      ZI=ZI+REFZ(I)
      LI=LI+1
      END DO
      XI=XI/LI
      YI=YI/LI
      ZI=ZI/LI
      LJ=0
      XJ=ZERO
      YJ=ZERO
      ZJ=ZERO
C compute geometric center for J-atoms
      DO JJ=JPR(ORR(COUNT)+1)+1,JPR(ORR(COUNT)+1+1)
      J=JLS(JJ)
      IF (.NOT.INITIA(J,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETNO-WRN: unknown reference coordinates for atom "',
     &SEGID(J),'-',RESID(J),'-',RES(J),'-',TYPE(J),'"'
      WRITE(6,'(2A)')
     &'               Automatic pseudoatom correction',
     &' wont be applied.'
      WRN=.TRUE.
      END IF
      XJ=XJ+REFX(J)
      YJ=YJ+REFY(J)
      ZJ=ZJ+REFZ(J)
      LJ=LJ+1
      END DO
      XJ=XJ/LJ
      YJ=YJ/LJ
      ZJ=ZJ/LJ
      END IF
      DO ICOUNT = IPR(ORR(COUNT)+1)+1, IPR(ORR(COUNT)+1+1)
      PH=ILS(ICOUNT)
C
      IF (CORRECT) THEN
C computer distance of atom ICOUNT to geometric center
      XIJ=REFX(ILS(ICOUNT))-XI
      YIJ=REFY(ILS(ICOUNT))-YI
      ZIJ=REFZ(ILS(ICOUNT))-ZI
      RI=SQRT(XIJ**2+YIJ**2+ZIJ**2)
      ELSE
      RI=ZERO
      END IF
      DO JCOUNT = JPR(ORR(COUNT)+1)+1, JPR(ORR(COUNT)+1+1)
C
      HILEN=HIL
      LOLEN=LOL
      IF (CORRECT) THEN
C computer distance of atom JCOUNT to geometric center
      XIJ=REFX(JLS(JCOUNT))-XJ
      YIJ=REFY(JLS(JCOUNT))-YJ
      ZIJ=REFZ(JLS(JCOUNT))-ZJ
      RJ=SQRT(XIJ**2+YIJ**2+ZIJ**2)
C
C update upper and lower bound
      IF (DODUMP) THEN
      WRITE(6,'(A,/,9A,9A)')
     & ' Pseudoatom correction applied to atom pair ',
     &'       "',SEGID(ILS(ICOUNT)),'-',RESID(ILS(ICOUNT)),
     &'-',RES(ILS(ICOUNT)),'-',TYPE(ILS(ICOUNT)),'"',
     &'       "',SEGID(JLS(JCOUNT)),'-',RESID(JLS(JCOUNT)),
     &'-',RES(JLS(JCOUNT)),'-',TYPE(JLS(JCOUNT)),'"'
      WRITE(6,'(A,2F10.4)')
     &  ' distance range before correction: ',HILEN, LOLEN
      END IF
      HILEN=HILEN+RI+RJ
      LOLEN=MAX(ZERO,LOLEN-RI-RJ)
      IF (DODUMP) THEN
      WRITE(6,'(A,2F10.4)')
     &  ' distance range after correction: ',HILEN, LOLEN
      END IF
      END IF
      PJ=JLS(JCOUNT)
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      R=DUB(PHH,PJJ)
      HILEN=MIN(R, HILEN)
      R=DUB(PJJ,PHH)
      LOLEN=MAX(R, LOLEN)
      DUB(PHH,PJJ) = HILEN
      DUB(PJJ,PHH) = LOLEN
      END DO
      END DO
      END IF
      END IF
      END DO
      RETURN
      END
C ======================
      SUBROUTINE DGSETDH(NSUB,DUB, CIP, CJP, CKP, CLP, CCPB, CCPO)
C
C Put dihedral angle restraints information into DUB matrix
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'cnst.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      INTEGER CIP(*), CJP(*), CKP(*), CLP(*)
      DOUBLE PRECISION CCPB(*), CCPO(*)
C local
      INTEGER COUNT, PH, PJ, PK, PL, PHH, PLL
      DOUBLE PRECISION HIPHI, LOPHI, HLMAX, HLMIN, R
C begin
C
C Translates J coupling-derived constraints on dihedral angles
C into distance constraints.
C CCPB (in cnst.inc) is the equilibrium value for the dihedral energy
C (=value of
C dihedral angle).  CCPO (in cnst.inc) is the offset of the dihedral energy.
C
      DO COUNT=1, CNPHI
      PH=CIP(COUNT)
      PJ=CJP(COUNT)
      PK=CKP(COUNT)
      PL=CLP(COUNT)
      HIPHI=CCPB(COUNT) + CCPO(COUNT)
      LOPHI=CCPB(COUNT) - CCPO(COUNT)
C
      CALL DGDIS14(NSUB,DUB,PH,PJ,PK,PL,HIPHI,LOPHI,HLMAX,HLMIN)
      PHH=MIN(PH,PL)
      PLL=MAX(PH,PL)
      R=DUB(PHH,PLL)
      HLMAX = MIN(R, HLMAX)
      R=DUB(PLL,PHH)
      HLMIN = MAX(R, HLMIN)
C
      DUB(PHH,PLL) = HLMAX
      DUB(PLL,PHH) = HLMIN
      END DO
      RETURN
      END
C ======================
      SUBROUTINE DGSETEN(NSUB, DUB, BACC, TACC, PACC, IMACC, SREFER)
C
C Sets up initial constraints from conformational energy terms: bonds,
C angles, dihedrals, impropers.
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'update.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'funct.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION BACC, TACC, PACC, IMACC
      CHARACTER*4 SREFER
C local variables
      INTEGER COUNT, I, J, PH, PJ, PK, PL, PHH, PJJ, PLL, WII
      INTEGER K, WI, WJ, WK, WL, PKK, WLL
      DOUBLE PRECISION HITHETA, LOWTHETA, HIPHI, LOPHI
      DOUBLE PRECISION HILEN, LOLEN, BLO, BUP
      DOUBLE PRECISION CLO, CUP, A1, A2, A3, A4, AMAX, AMIN
      DOUBLE PRECISION A5, A6, A7, A8
      DOUBLE PRECISION HLMAX, HLMIN, MAXLEN, MINLEN, DIST, R
      LOGICAL COND, ERR
C parameters
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
C
C begin
C
      ERR=.FALSE.
C
      IF (QENER(SSBOND)) THEN
C Read in all bonds from CBB (in param.inc) through ICB
      IF (SREFER.EQ.'PARA') THEN
C
C get type-based bond parameters if required
      IF (UPBOND) THEN
      CALL CODBON(NBOND,IB,JB,ACBC,ACBB,QACBC,QACBB,IAC,
     &            SEGID,RESID,TYPE)
      UPBOND=.FALSE.
      END IF
C
      DO COUNT=1, NBOND
      HILEN=ACBB(COUNT) + BACC
      LOLEN=ACBB(COUNT) - BACC
      MAXLEN=MAX(HILEN, LOLEN)
      MINLEN=MIN(HILEN, LOLEN)
      PH=IB(COUNT)
      PJ=JB(COUNT)
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      R=DUB(PHH,PJJ)
      HILEN=MIN(R, MAXLEN)
      R=DUB(PJJ,PHH)
      LOLEN=MAX(R, MINLEN)
      DUB(PHH,PJJ)=HILEN
      DUB(PJJ,PHH)=LOLEN
      END DO
      ELSE
      DO COUNT=1, NBOND
      PH=IB(COUNT)
      PJ=JB(COUNT)
      IF (.NOT.INITIA(PH,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PH),'-',RESID(PH),'-',RES(PH),'-',TYPE(PH),'"'
      ERR=.TRUE.
      END IF
      IF (.NOT.INITIA(PJ,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PJ),'-',RESID(PJ),'-',RES(PJ),'-',TYPE(PJ),'"'
      ERR=.TRUE.
      END IF
      DIST=SQRT((REFX(PH)-REFX(PJ))**2 +
     &          (REFY(PH)-REFY(PJ))**2 +
     &          (REFZ(PH)-REFZ(PJ))**2)
      HILEN=DIST + BACC
      LOLEN=DIST - BACC
      MAXLEN=MAX(HILEN, LOLEN)
      MINLEN=MIN(HILEN, LOLEN)
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      R=DUB(PHH,PJJ)
      HILEN=MIN(R, MAXLEN)
      R=DUB(PJJ,PHH)
      LOLEN=MAX(R, MINLEN)
      DUB(PHH,PJJ)=HILEN
      DUB(PJJ,PHH)=LOLEN
      END DO
      END IF
      END IF
C
      IF (QENER(SSANGL)) THEN
      IF (SREFER.EQ.'PARA') THEN
C Calculate angle constraints on D(hk), given D(jk) and D(hj) and the
C law of cosines.
C get type-based angle parameters if required
      IF (UPANGL) THEN
      CALL CODANG(NTHETA,IT,JT,KT,ACTC,
     &            ACTB,ACTUC,ACTUB,QACTC,QACTB,QACTUC,QACTUB,IAC,
     &            SEGID,RESID,TYPE)
      UPANGL=.FALSE.
      END IF
C
      DO COUNT=1, NTHETA
      PH=IT(COUNT)
      PJ=JT(COUNT)
      PK=KT(COUNT)
      HITHETA=ACTB(COUNT) + TACC
      LOWTHETA=ACTB(COUNT) - TACC
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      BUP=DUB(PHH,PJJ)
      BLO=DUB(PJJ,PHH)
      PJJ=MIN(PJ,PK)
      PKK=MAX(PJ,PK)
      CUP=DUB(PJJ,PKK)
      CLO=DUB(PKK,PJJ)
C
      A1=SQRT(BUP**2 + CUP**2 - TWO * BUP * CUP * COS(HITHETA))
      A2=SQRT(BLO**2 + CUP**2 - TWO * BLO * CUP * COS(HITHETA))
      A3=SQRT(BUP**2 + CLO**2 - TWO * BUP * CLO * COS(HITHETA))
      A4=SQRT(BLO**2 + CLO**2 - TWO * BLO * CLO * COS(HITHETA))
C include LOWTHETA as well (ATB)
      A5=SQRT(BUP**2 + CUP**2 - TWO * BUP * CUP * COS(LOWTHETA))
      A6=SQRT(BLO**2 + CUP**2 - TWO * BLO * CUP * COS(LOWTHETA))
      A7=SQRT(BUP**2 + CLO**2 - TWO * BUP * CLO * COS(LOWTHETA))
      A8=SQRT(BLO**2 + CLO**2 - TWO * BLO * CLO * COS(LOWTHETA))
C
      AMAX=MAX(A1, A2, A3, A4, A5, A6, A7, A8)
      AMIN=MIN(A1, A2, A3, A4, A5, A6, A7, A8)
      PHH=MIN(PH,PK)
      PKK=MAX(PH,PK)
      R=DUB(PHH,PKK)
      AMAX=MIN(R, AMAX)
      R=DUB(PKK,PHH)
      AMIN=MAX(R, AMIN)
      DUB(PHH,PKK) = AMAX
      DUB(PKK,PHH) = AMIN
      END DO
      ELSE
      DO COUNT=1, NTHETA
      PH=IT(COUNT)
      PK=KT(COUNT)
      IF (.NOT.INITIA(PH,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PH),'-',RESID(PH),'-',RES(PH),'-',TYPE(PH),'"'
      ERR=.TRUE.
      END IF
      IF (.NOT.INITIA(PK,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PK),'-',RESID(PK),'-',RES(PK),'-',TYPE(PK),'"'
      ERR=.TRUE.
      END IF
      DIST=SQRT((REFX(PH)-REFX(PK))**2 +
     &          (REFY(PH)-REFY(PK))**2 +
     &          (REFZ(PH)-REFZ(PK))**2)
      HILEN=DIST + BACC
      LOLEN=DIST - BACC
      AMAX=MAX(HILEN, LOLEN)
      AMIN=MIN(HILEN, LOLEN)
      PHH=MIN(PH,PK)
      PKK=MAX(PH,PK)
      R=DUB(PHH,PKK)
      AMAX=MIN(R, AMAX)
      R=DUB(PKK,PHH)
      AMIN=MAX(R, AMIN)
      DUB(PHH,PKK) = AMAX
      DUB(PKK,PHH) = AMIN
      END DO
      END IF
      END IF
C
      IF (QENER(SSDIHE)) THEN
C first, generate a "generic" dihedral angle list based on the
C connectivity (ATB).  The minimum and maximum values are calculated
C at zero and pi radians.
      DO I=1, NBOND
      DO J=1, NBOND
      IF (J.NE.I) THEN
      COND=.TRUE.
      IF (IB(I).EQ.IB(J)) THEN
      WK=JB(I)
      WJ=IB(I)
      WI=JB(J)
      ELSE IF (IB(I).EQ.JB(J)) THEN
      WK=JB(I)
      WJ=IB(I)
      WI=IB(J)
      ELSE IF (JB(I).EQ.IB(J)) THEN
      WK=IB(I)
      WJ=JB(I)
      WI=JB(J)
      ELSE IF (JB(I).EQ.JB(J)) THEN
      WK=IB(I)
      WJ=JB(I)
      WI=IB(J)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
      DO K=J+1, NBOND
      IF (K.NE.I) THEN
      COND=.TRUE.
      IF (IB(K).EQ.WK) THEN
      WL=JB(K)
      ELSE IF (JB(K).EQ.WK) THEN
      WL=IB(K)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
C found a "generic" dihedral angle
      HIPHI = PI
      LOPHI = ZERO
C
      CALL DGDIS14(NSUB,DUB,WI,WJ,WK,WL,HIPHI,LOPHI,HLMAX,HLMIN)
      WII=MIN(WI,WL)
      WLL=MAX(WI,WL)
      R=DUB(WII,WLL)
      HLMAX = MIN(R, HLMAX)
      R=DUB(WLL,WII)
      HLMIN = MAX(R, HLMIN)
      DUB(WII,WLL) = HLMAX
      DUB(WLL,WII) = HLMIN
      END IF
      END IF
      END DO
      END IF
      END IF
      END DO
      END DO
C
C get type-based dihedral parameters if required
C NOTE: we have to get the parameters for both the PARA and the
C COORdinate option as in the latter case we still need the
C periodicities from the parameter database.
      IF (UPDIHE) THEN
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,QACPC,QACPB,QACPD,
     &            IAC,SEGID,RESID,TYPE)
      UPDIHE=.FALSE.
      END IF
C
      IF (SREFER.EQ.'PARA') THEN
C
C Now translate the dihedral angles in the parameter list into
C distance constraints. If there is only one possible configuration of
C the distance, the upper & lower bounds are calculated at the given
C value of the angle (CPB in param.inc).  If there are more than one,
C they are the minimum and maximum values calculated at zero and pi
C radians.  The multiplicity is CPD (in param.inc).
C
      DO COUNT=1, NPHI
      PH=IP(COUNT)
      PJ=JP(COUNT)
      PK=KP(COUNT)
      PL=LP(COUNT)
C if ACPD is zero, this means it is an improper and the equilibrium value
C is ACPB.  If ACPD is one, this is means it is a COS function with a
C single minimum at -ACPB+PI (ATB).
      IF (ACPD(COUNT) .EQ. 0) THEN
      HIPHI=ACPB(COUNT) + PACC
      LOPHI=ACPB(COUNT) - PACC
      ELSE IF (ACPD(COUNT) .EQ. 1) THEN
      HIPHI=-ACPB(COUNT) + PI + PACC
      LOPHI=-ACPB(COUNT) + PI - PACC
      ELSE
      HIPHI = PI
      LOPHI = ZERO
      END IF
C
      CALL DGDIS14(NSUB,DUB,PH,PJ,PK,PL,HIPHI,LOPHI,HLMAX,HLMIN)
      PHH=MIN(PH,PL)
      PLL=MAX(PH,PL)
      R=DUB(PHH,PLL)
      HLMAX = MIN(R, HLMAX)
      R=DUB(PLL,PHH)
      HLMIN = MAX(R, HLMIN)
      DUB(PHH,PLL) = HLMAX
      DUB(PLL,PHH) = HLMIN
      END DO
      ELSE
C
C for the REFErence=COORdinate option we just include the
C improper-like dihedrals.
      DO COUNT=1, NPHI
      IF (ACPD(COUNT) .EQ. 0.OR.ACPD(COUNT).EQ.1) THEN
      PH=IP(COUNT)
      PL=LP(COUNT)
      IF (.NOT.INITIA(PH,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PH),'-',RESID(PH),'-',RES(PH),'-',TYPE(PH),'"'
      ERR=.TRUE.
      END IF
      IF (.NOT.INITIA(PL,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PL),'-',RESID(PL),'-',RES(PL),'-',TYPE(PL),'"'
      ERR=.TRUE.
      END IF
      DIST=SQRT((REFX(PH)-REFX(PL))**2 +
     &          (REFY(PH)-REFY(PL))**2 +
     &          (REFZ(PH)-REFZ(PL))**2)
      HILEN=DIST + BACC
      LOLEN=DIST - BACC
      AMAX=MAX(HILEN, LOLEN)
      AMIN=MIN(HILEN, LOLEN)
      PHH=MIN(PH,PL)
      PLL=MAX(PH,PL)
      R=DUB(PHH,PLL)
      AMAX=MIN(R, AMAX)
      R=DUB(PLL,PHH)
      AMIN=MAX(R, AMIN)
      DUB(PHH,PLL) = AMAX
      DUB(PLL,PHH) = AMIN
      END IF
      END DO
C
      END IF
      END IF
C
      IF (QENER(SSIMPR)) THEN
C get type-based improper parameters if required
C NOTE: we have to get the parameters for both the PARA and the
C COORdinate option as in the latter case we still need the
C periodicities from the parameter database.
      IF (UPIMPR) THEN
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,
     &            QACIC,QACIB,QACID,IAC,SEGID,RESID,TYPE)
      UPIMPR=.FALSE.
      END IF
      IF (SREFER.EQ.'PARA') THEN
C
C Now translate the improper angles in the parameter list into
C distance constraints. If there is only one possible configuration of
C the distance, the upper & lower bounds are calculated at the given
C value of the angle (ACIB in mtf.inc).  If there are more than one,
C they are the minimum and maximum values calculated at zero and pi
C radians.  The multiplicity is ACID (in mtf.inc).
      DO COUNT=1, NIMPHI
      PH=IM(COUNT)
      PJ=JM(COUNT)
      PK=KM(COUNT)
      PL=LM(COUNT)
C if CID is zero, this means it is an improper and the equilibrium value
C is CIB.  If CID is one, this is means it is a COS function with a
C single minimum at -CIB+PI (ATB).
      IF (ACID(COUNT) .EQ. 0) THEN
      HIPHI=ACIB(COUNT) + IMACC
      LOPHI=ACIB(COUNT) - IMACC
      ELSE IF (ACID(COUNT) .EQ. 1) THEN
      HIPHI=-ACIB(COUNT) + PI + IMACC
      LOPHI=-ACIB(COUNT) + PI - IMACC
      ELSE
      HIPHI = PI
      LOPHI = ZERO
      END IF
C
      CALL DGDIS14(NSUB,DUB,PH,PJ,PK,PL,HIPHI,LOPHI,HLMAX,HLMIN)
      PHH=MIN(PH,PL)
      PLL=MAX(PH,PL)
      R=DUB(PHH,PLL)
      HLMAX = MIN(R, HLMAX)
      R=DUB(PLL,PHH)
      HLMIN = MAX(R, HLMIN)
      DUB(PHH,PLL) = HLMAX
      DUB(PLL,PHH) = HLMIN
      END DO
      ELSE
C
C for the REFErence=PARAmeter option we just include the
C improper-like dihedrals.
      DO COUNT=1, NIMPHI
      IF (ACID(COUNT) .EQ. 0.OR.ACID(COUNT).EQ.1) THEN
      PH=IM(COUNT)
      PL=LM(COUNT)
      IF (.NOT.INITIA(PH,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PH),'-',RESID(PH),'-',RES(PH),'-',TYPE(PH),'"'
      ERR=.TRUE.
      END IF
      IF (.NOT.INITIA(PL,REFX,REFY,REFZ)) THEN
      WRITE(6,'(9A)')
     &' %DGSETEN-ERR: unknown reference coordinates for atom "',
     &SEGID(PL),'-',RESID(PL),'-',RES(PL),'-',TYPE(PL),'"'
      ERR=.TRUE.
      END IF
      DIST=SQRT((REFX(PH)-REFX(PL))**2 +
     &          (REFY(PH)-REFY(PL))**2 +
     &          (REFZ(PH)-REFZ(PL))**2)
      HILEN=DIST + BACC
      LOLEN=DIST - BACC
      AMAX=MAX(HILEN, LOLEN)
      AMIN=MIN(HILEN, LOLEN)
      PHH=MIN(PH,PL)
      PLL=MAX(PH,PL)
      R=DUB(PHH,PLL)
      AMAX=MIN(R, AMAX)
      R=DUB(PLL,PHH)
      AMIN=MAX(R, AMIN)
      DUB(PHH,PLL) = AMAX
      DUB(PLL,PHH) = AMIN
      END IF
      END DO
C
      END IF
      END IF
C End of SETUP subroutine
C
      IF (ERR) THEN
      CALL WRNDIE(-5,'DGSETRG','Unknown reference coordinates')
      END IF
C
      RETURN
      END
C ======================
      SUBROUTINE DGSETNX(NPAIR,IPK,JPK,IBLOEX,INBEX,QSCAN)
C
C Routine generates nonbonded exclusion list, i.e., a list
C that excludes all 1--2 and 1--3 interactions.  Note that
C this list is symmetric, i.e., for a given atom it produces
C ALL 1--2 and 1--3 exclusions.
C
C If QSCAN is true, we determine the number of exclusions
C without actually storing them, otherwise we store the
C exclusions in the IBLOEX, INBEX list.
C
C Author: Axel T. Brunger
C
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'funct.inc'
      INTEGER NPAIR, IPK(*), JPK(*)
      INTEGER IBLOEX(*), INBEX(*)
      LOGICAL QSCAN
      EXTERNAL EXCH5
C local
      INTEGER I, J, WI, WK, IATOM
      LOGICAL COND
C begin
C
      NPAIR=0
      IF (QSCAN) THEN
      NPAIR=NPAIR+2*NBOND
      ELSE
      DO I=1, NBOND
      NPAIR=NPAIR+1
      IPK(NPAIR)=IB(I)
      JPK(NPAIR)=JB(I)
      NPAIR=NPAIR+1
      IPK(NPAIR)=JB(I)
      JPK(NPAIR)=IB(I)
      END DO
      END IF
C
      DO I=1, NBOND
      DO J=1, NBOND
      IF (J.NE.I) THEN
      COND=.TRUE.
      IF (IB(I).EQ.IB(J)) THEN
      WK=JB(I)
      WI=JB(J)
      ELSE IF (IB(I).EQ.JB(J)) THEN
      WK=JB(I)
      WI=IB(J)
      ELSE IF (JB(I).EQ.IB(J)) THEN
      WK=IB(I)
      WI=JB(J)
      ELSE IF (JB(I).EQ.JB(J)) THEN
      WK=IB(I)
      WI=IB(J)
      ELSE
      COND=.FALSE.
      END IF
      IF (COND) THEN
C found a 1--3 interaction
      NPAIR=NPAIR+1
      IF (.NOT.QSCAN) THEN
      IPK(NPAIR)=WI
      JPK(NPAIR)=WK
      END IF
      END IF
      END IF
      END DO
      END DO
C
      IF (.NOT.QSCAN) THEN
C
C Sort the pair list.
      CALL SORT(NPAIR,EXCH5,ORDER5,IPK,JPK,0,0,0,0,0,2)
C
C convert it into a pointer list
      IATOM=1
      DO I=1,NPAIR
      IBLOEX(1)=0
      IF (IPK(I).GT.IATOM) THEN
      DO J=IATOM,IPK(I)-1
      IBLOEX(J+1)=I-1
      END DO
      IATOM=IPK(I)
      END IF
      INBEX(I)=JPK(I)
      END DO
      DO J=IATOM,NATOM
      IBLOEX(J+1)=NPAIR
      END DO
C
      WRITE(6,'(A,I6,A)') ' DGSETNX: found ',NPAIR,
     & ' exclusions (1--2, 1--3) for repulsive lower bounds.'
      END IF
      RETURN
      END
C ======================
      SUBROUTINE DGSETNN(ROOT, NSUB, SUB, DUB, IBLOEX, INBEX)
C Sets up initial constraints from van der Waals repulsion for
C all lower bounds between the ROOT atom and any other atom.
C
C Front end for DGSETN2
C
C Author: Axel T. Brunger
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'mtf.inc'
C passed-in variables
      INTEGER ROOT, NSUB, SUB(*)
      INTEGER IBLOEX(*), INBEX(*)
      REAL DUB(NSUB,*)
C local variables
      INTEGER NSELE
C pointer
      INTEGER EXCPTR
C begin
      EXCPTR=ALLHP(INTEG4(NATOM))
      DO NSELE=1,NPIG
      CALL DGSETN2(ROOT, NSUB, SUB, DUB, HEAP(EXCPTR),
     &     IBLOEX, INBEX, HEAP(IINTER(NSELE)))
      END DO
      CALL FREHP(EXCPTR,INTEG4(NATOM))
      RETURN
      END
C ======================
      SUBROUTINE DGSETN2(ROOT, NSUB, SUB, DUB, NONEXCL, IBLOEX, INBEX,
     &                   INTERE )
C
C Sets up initial constraints from van der Waals repulsion for
C all lower bounds between the ROOT atom and any other atom.
C This routine operates on substructures.
C
C by John Kuszewksi 7/20/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'param.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
C passed-in variables
      INTEGER ROOT, NSUB, SUB(*)
      REAL DUB(NSUB,*)
      INTEGER NONEXCL(*)
      INTEGER IBLOEX(*), INBEX(*), INTERE(*)
C local variables
      INTEGER I, PJ, COUNT, INROOT
      DOUBLE PRECISION SUM, R, SCALE
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
C Read in van der Waals radii from CNBVR (in param.inc), multiply
C by the scaling factor FREPEL (in nbonds.inc).  LOOKUP (in mtf.inc)
C gives the chemical
C type of an atom, which is the index used by CNBVR.
C
C get nonbonded exclusions between ROOT and any other atom.
C -- store in bookkeeping array NONEXCL
      DO I=1,NATOM
      NONEXCL(I) = 0
      END DO
      DO COUNT = IBLOEX(SUB(ROOT))+1, IBLOEX(SUB(ROOT)+1)
      NONEXCL(INBEX(COUNT)) = 1
      END DO
C
CMN      RH=CNBVR(LOOKUP(SUB(ROOT)))
      IF (FREPEL.GT.ZERO) THEN
      SCALE=FREPEL
      ELSE
      SCALE=ONE
      END IF
C
C only put a repulsion as a lower bound if
C  1. it is included in the INTERE interaction array
C  2. it is consistent with the current
C     upper bound and is tighter than the current lower bound
C  3. it is not excluded.
      INROOT=INTERE(SUB(ROOT))
      DO PJ=1, ROOT-1
CMN      RJ=CNBVR(LOOKUP(SUB(PJ)))
CMN      SUM=SCALE * (RH+RJ)
      SUM=SCALE * CNBVR(LOOKUP(SUB(ROOT)),LOOKUP(SUB(PJ)))
      IF    ( (DUB(PJ,ROOT) .GT. SUM)
     &  .AND. (NONEXCL(SUB(PJ)).EQ.0)
     &  .AND. (INROOT+INTERE(SUB(PJ)).LE.+1) ) THEN
      R=DUB(ROOT,PJ)
      DUB(ROOT,PJ)=MAX(R, SUM)
      END IF
      END DO
      DO PJ=ROOT+1, NSUB
CMN      RJ=CNBVR(LOOKUP(SUB(PJ)))
CMN      SUM=SCALE * (RH+RJ)
      SUM=SCALE * CNBVR(LOOKUP(SUB(ROOT)),LOOKUP(SUB(PJ)))
      IF ((DUB(ROOT,PJ) .GT. SUM) .AND. (NONEXCL(SUB(PJ)).EQ.0)) THEN
      R=DUB(PJ,ROOT)
      DUB(PJ,ROOT)=MAX(R, SUM)
      END IF
      END DO
C
      RETURN
      END
C ======================
      SUBROUTINE DGGROUP(FLAGS, MFLAGS, LIST, POINTR,
     &                   ACCS, ACC,NGROUP, NFLAGS)
C
C Sets up arrays for rigid groups.  Similar to subroutine MRIGI2 in MRIGID.S
C LIST is an array of atom numbers involved in rigid groups.
C POINTR is an array of indexes to LIST, where the starting point in LIST
C for rigid group X is POINTER(X).  ACCS is a list of accuracies for the
C various rigid groups.  ACC is the accuracy for the current rigid group.
C
C by John Kuszewski 6/27/90
C ======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
C passed-in variables
      INTEGER FLAGS(*), MFLAGS, LIST(*), POINTR(*)
      INTEGER NGROUP, NFLAGS
      DOUBLE PRECISION ACC, ACCS(*)
C local variables
      INTEGER X
C begin
      IF (NGROUP .EQ. 0) POINTR(1) = 1
      IF (NFLAGS.GT.0) THEN
      IF (NGROUP.GE.NATOM) THEN
      CALL WRNDIE(-5,'DGGROUP','number of rigid groups exceeded.')
      ELSE
      NGROUP=NGROUP+1
      IF (NFLAGS+POINTR(NGROUP)-1.GT.MFLAGS) THEN
      CALL WRNDIE(-5,'DGGROUP',
     &   'total number of selectable atoms exceeded.')
      NGROUP=NGROUP-1
      ELSE
      DO X=1, NFLAGS
      LIST(X + POINTR(NGROUP) - 1) = FLAGS(X)
      END DO
      POINTR(NGROUP+1)=POINTR(NGROUP) + NFLAGS
      ACCS(NGROUP) = ABS(ACC)
      END IF
      END IF
      END IF
      RETURN
      END
C ===================
      SUBROUTINE DGDIS14(NSUB,DUB,
     &                   PH,PJ,PK,PL,HIPHI,LOPHI,MAXVAL,MINVAL)
C
C Calculates the 1,4 distance across a dihedral angle.  See my
C notes for the derivation.  HL is the distance between PH & PL
C in the plane of the paper.
C
C Tries all permutations of upper & lower bounds on distances and angles
C to calculate bounds from dihedral & improper angles
C
C by John Kuszewski 2/15/90 7/17/90
C ====================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
C passed-in variables
      INTEGER PH, PJ, PK, PL, NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION MAXVAL, MINVAL
      DOUBLE PRECISION HIPHI, LOPHI
C local variables
      DOUBLE PRECISION HJUP, HJLO, JKUP, JKLO, KLUP, KLLO, HKUP, HKLO
      DOUBLE PRECISION JLUP, JLLO, HJ, JK, KL, HK, JL, THISD
      DOUBLE PRECISION HX, JX, KY, LY, HL, XY, CALC1,CALC2,THETA1,THETA2
      INTEGER A, B, C, D, E, F, PHH, PJJ, PLL, PKK
      DOUBLE PRECISION PHI
      LOGICAL OK
      DOUBLE PRECISION ZERO, ONE, BIG
C parameter
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, BIG=1000D0)
C begin
      PHH=MIN(PH,PJ)
      PJJ=MAX(PH,PJ)
      HJUP=DUB(PHH,PJJ)
      HJLO=DUB(PJJ,PHH)
      PJJ=MIN(PJ,PK)
      PKK=MAX(PJ,PK)
      JKUP=DUB(PJJ,PKK)
      JKLO=DUB(PKK,PJJ)
      PKK=MIN(PK,PL)
      PLL=MAX(PK,PL)
      KLUP=DUB(PKK,PLL)
      KLLO=DUB(PLL,PKK)
      PHH=MIN(PH,PK)
      PKK=MAX(PH,PK)
      HKUP=DUB(PHH,PKK)
      HKLO=DUB(PKK,PHH)
      PJJ=MIN(PJ,PL)
      PLL=MAX(PJ,PL)
      JLUP=DUB(PJJ,PLL)
      JLLO=DUB(PLL,PJJ)
C
      IF ((HJUP .GE. R4BIG-BIG) .OR. (JKUP .GE. R4BIG-BIG) .OR.
     &   (KLUP .GE. R4BIG-BIG) .OR. (HKUP .GE. R4BIG-BIG) .OR.
     &   (JLUP .GE. R4BIG-BIG) .OR. (HJLO .LE. R4SMAL) .OR.
     &   (JKLO .LE. R4SMAL) .OR. (KLLO .LE. R4SMAL) .OR.
     &   (HKLO .LE. R4SMAL) .OR. (JLLO .LE. R4SMAL)) THEN
      MAXVAL = R4BIG
      MINVAL = ZERO
      ELSE
C
      MAXVAL = ZERO
      MINVAL = R4BIG
      OK = .FALSE.
      DO A = 1,2
      DO B = 1,2
      DO C = 1,2
      DO D = 1,2
      DO E = 1,2
      DO F = 1,2
      HJ = HJUP
      JK = JKUP
      KL = KLUP
      HK = HKUP
      JL = JLUP
      PHI = HIPHI
      IF (A .EQ. 2) HJ = HJLO
      IF (B .EQ. 2) JK = JKLO
      IF (C .EQ. 2) KL = KLLO
      IF (D .EQ. 2) HK = HKLO
      IF (E .EQ. 2) JL = JLLO
      IF (F .EQ. 2) PHI = LOPHI
C
      CALC1=(HK**2 - JK**2 - HJ**2)/(-2 * JK * HJ)
      CALC2=(JL**2 - JK**2 - KL**2)/(-2 * JK * KL)
      IF ((CALC1 .LE. ONE) .AND. (CALC1 .GE. -ONE) .AND.
     &    (CALC2 .LE. ONE) .AND. (CALC2 .GE. -ONE)) THEN
      THETA1=ACOS(CALC1)
      THETA2=ACOS(CALC2)
      JX=HJ * COS(THETA1)
      HX=HJ * SIN(THETA1)
      KY=KL * COS(THETA2)
      LY=KL * SIN(THETA2)
      XY = JK - JX - KY
      HL=HX**2 + LY**2 - 2 * HX * LY * COS(PHI)
      IF (HL .GE. ZERO) THEN
      HL = SQRT(HL)
      THISD=SQRT(HL**2 + XY**2)
      MAXVAL = MAX(MAXVAL, THISD)
      MINVAL = MIN(MINVAL, THISD)
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END IF
      RETURN
      END
C =========================
      SUBROUTINE DGSHRTP (NSUB,DUB,PATHLEN,LBLPATHS,ROOT,OLDROOT)
C
C Dijkstra-based single-source shortest path routine.
C See M. Fredman and R.E. Tarjan, "Fibonacci heaps and
C their uses in improved network optimization algorithms" in
C ACM Journal 34:3, p. 596-615 (1987) for more info
C
C This routine is 3-4 times faster than the Dial algorithm for
C an alanine 20-mer or 60-mer. 7/10/90
C
C This has to ensure that the positive nodes are taken care of first
C for the triangle inequality to work. 7/12/90
C
C Instead of keeping a logical array to show which atoms have been
C eliminated, this implementation keeps another copy of the pathlen
C array (called labelpaths) which has the current estimated pathlength
C for labeled atoms and R4BIG for the scanned atoms.  This allows
C the use of ISMIN to find the nearest non-scanned atom.
C
C by John Kuszewski 7/6/90
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
C function names
      INTEGER DGMNDX
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION PATHLEN(-NSUB:NSUB)
      DOUBLE PRECISION LBLPATHS(-NSUB:NSUB)
      INTEGER ROOT, OLDROOT
C local variables
      INTEGER X, W, V, PHH, PKK
      LOGICAL DONEALL, DONEUP
      DOUBLE PRECISION NEWLEN, PV
C parameters
      DOUBLE PRECISION ZERO, BIG
      PARAMETER (ZERO=0.0D0, BIG=1000D0)
C begin
C
C initialize the LBLPATHS array
      DO X=-NSUB,NSUB
      LBLPATHS(X)=R4BIG
      END DO
C
C if this is the same root as before don't touch the PATHLEN array
      IF (OLDROOT .NE. ROOT) THEN
      IF (OLDROOT .NE. 0) THEN
      PHH=MIN(ROOT,OLDROOT)
      PKK=MAX(ROOT,OLDROOT)
      DO X = 1, NSUB
      PATHLEN(X) = MIN(R4BIG, PATHLEN(X) + DUB(PHH, PKK))
      END DO
      DO X = -NSUB, -1
      PATHLEN(X) = MIN(ZERO, PATHLEN(X) + DUB(PHH, PKK))
      END DO
      ELSE
      DO X = 1, NSUB
      PATHLEN(X) = R4BIG
      END DO
      DO X = -NSUB, -1
      PATHLEN(X) = ZERO
      END DO
      END IF
      END IF
C
      PATHLEN(ROOT) = ZERO
      LBLPATHS(ROOT) = ZERO
      DONEALL = .FALSE.
      DONEUP = .FALSE.
C
      DO WHILE (.NOT. DONEALL)
      IF (.NOT. DONEUP) V = DGMNDX(NSUB, LBLPATHS(1))
      IF (LBLPATHS(V) .GE. R4BIG-BIG) DONEUP = .TRUE.
      IF (DONEUP) V = DGMNDX(NSUB, LBLPATHS(-NSUB)) - NSUB - 1
      IF (LBLPATHS(V) .GE. R4BIG-BIG) THEN
      DONEALL = .TRUE.
      ELSE
      LBLPATHS(V) = R4BIG
      PV=PATHLEN(V)
      IF (V .LT. 0) THEN
C
C upper bounds
      DO W = -NSUB, +V-1
      NEWLEN=PV + DUB(-V,-W)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      DO W = +V+1, -1
      NEWLEN=PV + DUB(-W,-V)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      ELSE
C
C lower bounds
      DO W = -NSUB, -V-1
      NEWLEN=PV - DUB(-W,V)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      DO W = -V+1, -1
      NEWLEN=PV - DUB(V,-W)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      END IF
      IF (.NOT. DONEUP) THEN
C
C upper bounds
      DO W = 1, V-1
      NEWLEN=PV + DUB(W,V)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      DO W = V+1, NSUB
      NEWLEN=PV + DUB(V,W)
      IF (NEWLEN .LT. PATHLEN(W)) THEN
      PATHLEN(W) = NEWLEN
      LBLPATHS(W) = NEWLEN
      END IF
      END DO
      END IF
      END IF
      END DO
C
C update DUB from PATHLEN
C
C lower bounds
      DO X = -NSUB, -ROOT-1
      DUB(-X,ROOT) = -PATHLEN(X)
      END DO
      DO X = -ROOT+1, -1
      DUB(ROOT,-X) = -PATHLEN(X)
      END DO
C
C upper bounds
      DO X = 1, ROOT-1
      DUB(X,ROOT) = PATHLEN(X)
      END DO
      DO X = ROOT+1, NSUB
      DUB(ROOT,X) = PATHLEN(X)
      END DO
      RETURN
      END
C =========================
      SUBROUTINE DGSHRT2 (NSUB,DUB,PATHLEN,LOWLBL,UPLBL,LOWIND,UPIND,
     &                    ROOT,OLDROOT,TREEPNT,TREE,UCOUNT,LCOUNT)
C
C Dijkstra-type shortest path routine using a
C tree-like storage for the (sparse) bounds matrix.
C
C by Axel T. Brunger
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'consta.inc'
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION PATHLEN(-NSUB:NSUB)
      INTEGER LOWLBL(NSUB), UPLBL(NSUB), LOWIND(NSUB), UPIND(NSUB)
      INTEGER ROOT, OLDROOT, TREEPNT(*), TREE(*), UCOUNT, LCOUNT
C local variables
      INTEGER X, W, V, PHH, PKK, WW, VV
      DOUBLE PRECISION NEWLEN, MMIN
      INTEGER NUPLBL, NLOWLBL, MLBL
C parameters
      DOUBLE PRECISION ZERO, BIG
      PARAMETER (ZERO=0.0D0, BIG=1000D0)
C begin
C
C initialize the LOWLBL, UPLBL, LOWIND, UPIND arrays
      DO X=1,NSUB
      LOWIND(X)=0
      UPIND(X)=0
      END DO
      NUPLBL=1
      UPLBL(1)=ROOT
      UPIND(ROOT)=1
      NLOWLBL=1
      LOWLBL(1)=ROOT
      LOWIND(ROOT)=1
C
C if this is the same root as before don't touch the PATHLEN array
      IF (OLDROOT .NE. ROOT) THEN
      IF (OLDROOT .NE. 0) THEN
      PHH=MIN(ROOT,OLDROOT)
      PKK=MAX(ROOT,OLDROOT)
      DO X = 1, NSUB
      PATHLEN(X) = MIN(R4BIG, PATHLEN(X) + DUB(PHH, PKK))
      END DO
      DO X = -NSUB, -1
      PATHLEN(X) = MIN(ZERO, PATHLEN(X) + DUB(PHH, PKK))
      END DO
      ELSE
      DO X = 1, NSUB
      PATHLEN(X) = R4BIG
      END DO
      DO X = -NSUB, -1
      PATHLEN(X) = ZERO
      END DO
      END IF
      END IF
C
      PATHLEN(ROOT) = ZERO
      PATHLEN(-ROOT) = ZERO
C
C first determine the upper bounds
C ================================
      DO WHILE (NUPLBL.GT.0)
C
C determine minimum element in UPLBL array
      MLBL=1
      MMIN=PATHLEN(UPLBL(1))
      UCOUNT=MAX(UCOUNT,NUPLBL)
      DO X=2, NUPLBL
      IF (PATHLEN(UPLBL(X)).LT.MMIN) THEN
      MLBL=X
      MMIN=PATHLEN(UPLBL(X))
      END IF
      END DO
      V=UPLBL(MLBL)
C
C remove this element from UPLBL array
      UPIND(UPLBL(MLBL))=0
      UPLBL(MLBL)=UPLBL(NUPLBL)
      NUPLBL=NUPLBL-1
C
C loop through upper bounds
      DO W = TREEPNT(V)+1, TREEPNT(V+1)
      WW=MIN(TREE(W),V)
      VV=MAX(TREE(W),V)
      NEWLEN=PATHLEN(V) + DUB(WW,VV)
      IF (NEWLEN .LT. PATHLEN(TREE(W))) THEN
      PATHLEN(TREE(W)) = NEWLEN
C
C add this element to the UPLBL list
      IF (UPIND(TREE(W)).EQ.0) THEN
      NUPLBL=NUPLBL+1
      UPLBL(NUPLBL)=TREE(W)
      UPIND(TREE(W))=NUPLBL
      END IF
      END IF
      END DO
      END DO
C
C now determine the lower bounds
C ==============================
      DO WHILE (NLOWLBL.GT.0)
C
C determine minimum element in LOWLBL array
      MLBL=1
      MMIN=PATHLEN(LOWLBL(1))
      LCOUNT=MAX(LCOUNT,NLOWLBL)
      DO X=2, NLOWLBL
      IF (PATHLEN(LOWLBL(X)).LT.MMIN) THEN
      MLBL=X
      MMIN=PATHLEN(LOWLBL(X))
      END IF
      END DO
      V=-LOWLBL(MLBL)
C
C remove this element from UPLBL array
      LOWIND(LOWLBL(MLBL))=0
      LOWLBL(MLBL)=LOWLBL(NLOWLBL)
      NLOWLBL=NLOWLBL-1
C
C loop through upper and lower bounds
      DO W = TREEPNT(-V)+1, TREEPNT(-V+1)
      WW=MIN(TREE(W),-V)
      VV=MAX(TREE(W),-V)
      NEWLEN=MIN(PATHLEN(V) + DUB(WW,VV), PATHLEN(-V)-DUB(VV,WW))
      IF (NEWLEN .LT. PATHLEN(-TREE(W))) THEN
      PATHLEN(-TREE(W)) = NEWLEN
C
C add this element to the LOWLBL list
      IF (LOWIND(TREE(W)).EQ.0) THEN
      NLOWLBL=NLOWLBL+1
      LOWLBL(NLOWLBL)=TREE(W)
      LOWIND(TREE(W))=NLOWLBL
      END IF
      END IF
      END DO
      END DO
C
C update DUB from PATHLEN
C
C lower bounds
      DO X = 1, ROOT-1
      DUB(ROOT,X) = -PATHLEN(-X)
      END DO
      DO X = ROOT+1, NSUB
      DUB(X,ROOT) = -PATHLEN(-X)
      END DO
C
C upper bounds
      DO X = 1, ROOT-1
      DUB(X,ROOT) = PATHLEN(X)
      END DO
      DO X = ROOT+1, NSUB
      DUB(ROOT,X) = PATHLEN(X)
      END DO
      RETURN
      END
C =========================
      SUBROUTINE DGTREE(QSCAN, NTREE, TREEPNT, TREE, NSUB, DUB)
C
C Routine determines tree of non-trivial lower and upper bounds.
C It stores the tree as a list of pointers TREEPNT into the
C array TREE.  For a given bound both directions (i.e., atom i->j
C and j->) are stored.
C
C Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'consta.inc'
      LOGICAL QSCAN
      INTEGER NTREE, TREEPNT(*), TREE(*), NSUB
      REAL DUB(NSUB,*)
C local
      INTEGER I, J
C parameter
      DOUBLE PRECISION BIG
      PARAMETER (BIG=1000.D0)
C begin
      NTREE=0
      IF (QSCAN) THEN
C
C just get the number of non-trivial bounds
      DO I=1,NSUB
      DO J=I+1,NSUB
      IF (DUB(I,J).LT.R4BIG-BIG.OR.DUB(J,I).GT.R4SMAL) THEN
      NTREE=NTREE+2
      END IF
      END DO
      END DO
      ELSE
      TREEPNT(1)=0
      DO I=1,NSUB
      TREEPNT(I+1)=TREEPNT(I)
      DO J=1,I-1
      IF (DUB(J,I).LT.R4BIG-BIG.OR.DUB(I,J).GT.R4SMAL) THEN
      TREEPNT(I+1)=TREEPNT(I+1)+1
      TREE(TREEPNT(I+1))=J
      END IF
      END DO
      DO J=I+1,NSUB
      IF (DUB(I,J).LT.R4BIG-BIG.OR.DUB(J,I).GT.R4SMAL) THEN
      TREEPNT(I+1)=TREEPNT(I+1)+1
      TREE(TREEPNT(I+1))=J
      END IF
      END DO
      END DO
      NTREE=TREEPNT(NSUB+1)
      END IF
      RETURN
      END
C==========================
      SUBROUTINE DGMETRI(DORANDOM, DODUMP, OK, SLCT,NSLCT,NMETR, RANDS,
     &                   NSUB, SUB, DUB, IBLOEX, INBEX)
C
C Picks random values between bounds and sets both bounds to the
C picked value, tightening other bounds
C as it goes.  Works in O(n^4) worst-case time.  GGUBFS in UTIL.S is a
C random number generator that gives a flat distribution of values.
C
C by John Kuszewski 8/24/89 12/29/89
C =========================
      IMPLICIT NONE
C function names
      INTEGER DGMNDX
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ener.inc'
C passed-in variables
      LOGICAL DODUMP, DORANDOM, OK
      INTEGER SLCT(*), NSLCT, NMETR, SUB(*), NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION RANDS(*)
      INTEGER IBLOEX(*), INBEX(*)
C local variables
      INTEGER H, J, PH, NUMPAIR, OLDROOT, X, NEWSLCT, HH, JJ
      DOUBLE PRECISION RANDOM, GAP, UP, LOW, TEMP
C pointer
      INTEGER PATHPTR, LBLPTHPTR
C parameters
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION ZERO, BIG
      PARAMETER (ZERO=0.0D0, BIG=1000D0)
C begin
C
C initialize total gap between upper & lower bounds to 0,
C number of distances set to 0
C set up list of random numbers for random metrization
C
      NMETR = MIN(NSLCT, ABS(NMETR))
      PATHPTR=ALLHP(IREAL8(NSUB*2+1))
      LBLPTHPTR=ALLHP(IREAL8(NSUB*2+1))
      GAP = ZERO
      NUMPAIR = 0
      OLDROOT = 0
      OK = .TRUE.
      IF (DORANDOM) THEN
      DO X = 1, NSLCT
      CALL GGUBFS(RANDOM)
      RANDS(X) = RANDOM
      END DO
      END IF
      IF (DODUMP) WRITE(6, '(2A)')
     &'     atom h         atom j            upperbound(h,j)',
     &'      distance(h,j)        lowerbound(h,j)'
C
C go through NMETR atoms in metrize's selection,
C in random order if you're doing random metrization,
C and mark those you actually use in the metrization
C by adding NSUB to them.
C
      DO PH = 1, NMETR
      IF (DORANDOM) THEN
      X = DGMNDX(NSLCT, RANDS)
      RANDS(X) = R4BIG
      ELSE
      X = PH
      END IF
      H = SLCT(X)
      SLCT(X) = SLCT(X) + NSUB
C
C now include repulsive information "on-the-fly" for all
C lower bounds from H to any other atom.
      IF (QENER(SSVDW)) THEN
      CALL DGSETNN(H, NSUB, SUB, DUB, IBLOEX, INBEX)
      END IF
C
C retighten the bounds matrixes for the new tree before picking
C a distance
C otherwise, the distance could be greater than the real upper
C bound (see notes on ala 2-mer trace)
C
      CALL DGSHRTP(NSUB,DUB,HEAP(PATHPTR),HEAP(LBLPTHPTR),H,
     &             OLDROOT)
      OLDROOT = H
C
C go through each other atom in order
C
      DO J = 1, NSUB
C
C if this distance hasn't been set already,
C update number of distances set,
C add this distances upper-lower gap to the total,
C determine a random value between the upper and lower bounds,
C set the distance, upper bounds, and lower bounds matrixes to that value,
C and retighten the bounds matrixes.
C
C if upper bound is not equal lower bound then randomly pick
C a distance and set both bounds to it (ATB)
      HH=MIN(H,J)
      JJ=MAX(H,J)
      UP=DUB(HH,JJ)
      LOW=DUB(JJ,HH)
      IF (ABS(LOW - UP) .GT. R4SMAL) THEN
      NUMPAIR = NUMPAIR + 1
      CALL GGUBFS(RANDOM)
C
      IF (UP .LT. LOW) THEN
      OK = .FALSE.
      TEMP = UP + ((LOW - UP) * RANDOM)
      WRITE(6, '(6(A4,1X),3(F15.4,1X),A)') SEGID(H), RESID(H),
     &TYPE(H), SEGID(J), RESID(J), TYPE(J), UP,TEMP,LOW,'*'
C
C reset the root
      OLDROOT=0
      ELSE
      TEMP = LOW + ((UP - LOW) * RANDOM)
      IF (DODUMP) THEN
      WRITE(6, '(6(A4,1X),3(F15.4,1X))') SEGID(H), RESID(H),
     &TYPE(H), SEGID(J), RESID(J), TYPE(J), UP,TEMP,LOW
      END IF
      END IF
C
      GAP = GAP + DUB(HH,JJ) - DUB(JJ,HH)
      DUB(HH,JJ)=TEMP
      DUB(JJ,HH)=TEMP
C
      CALL DGSHRTP(NSUB,DUB,HEAP(PATHPTR),HEAP(LBLPTHPTR),
     &             H,OLDROOT)
C
      END IF
      END DO
      END DO
C
C Remove all non-marked atoms from the selection.
C This allows the no-metrize routine to just avoid
C the re-tightened columns in its boundsmoothing.
C
      NEWSLCT = 0
      DO X = 1, NSLCT
      IF (SLCT(X) .GT. NSUB) THEN
      NEWSLCT = NEWSLCT + 1
      SLCT(NEWSLCT) = SLCT(X) - NSUB
      END IF
      END DO
      NSLCT = NEWSLCT
C
C calculate the average gap between the upper and lower
C bounds on distances set in this subroutine
C
      IF (NUMPAIR .NE. 0) GAP = GAP / NUMPAIR
      CALL DECLAR('MET_GAP', 'DP', ' ', DBCOMP, GAP)
      CALL FREHP(PATHPTR, IREAL8(NSUB*2+1))
      CALL FREHP(LBLPTHPTR, IREAL8(NSUB*2+1))
      RETURN
      END
C =========================
      SUBROUTINE DGNOMET (DODUMP, OK, SLCT, NSLCT, NSUB, SUB, DUB,
     &                    IBLOEX, INBEX, SSHORT, NOSMOOTH)
C
C Picks random values between bounds but
C does not retighten the matrix after each distance setting.
C Fills both bounds with the picked value (ATB).
C Works in O(n^3) worst-case time
C
C by John Kuszewski 8/24/89 12/29/89 6/19/91
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'ener.inc'
C passed-in variables
      LOGICAL DODUMP, OK
      INTEGER SLCT(*), SUB(*), NSLCT, NSUB
      REAL DUB(NSUB,*)
      INTEGER IBLOEX(*), INBEX(*)
      CHARACTER*4 SSHORT
      LOGICAL NOSMOOTH
C local variables
      INTEGER H, J, OLDROOT, SLCTCNT, HH, JJ
      DOUBLE PRECISION GAP, RANDOM, UP, LOW, TEMP
      INTEGER NUMPAIR, NTREE, UCOUNT, LCOUNT
      DOUBLE COMPLEX DBCOMP
      LOGICAL QDOIT
C pointer
      INTEGER PATHPTR, LBLPTHPTR, LOWLBL, UPLBL, TREE, TREEPNT
      INTEGER LOWIND, UPIND
C parameters
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C begin
C
C initialize
      GAP = ZERO
      NUMPAIR = 0
      OLDROOT = 0
      SLCTCNT = 1
      IF (DODUMP) WRITE(6, '(2A)')
     &'     atom h         atom j            upperbound(h,j) ',
     &'     distance(h,j)        lowerbound(h,j)'
      IF (.NOT.NOSMOOTH) THEN
C
C bound-smoothing.
C go through all atoms ph in the molecule in order
      QDOIT=.FALSE.
      IF (SSHORT.EQ.'AUTO'.OR.SSHORT.EQ.'SPAR') THEN
      CALL DGTREE(.TRUE., NTREE, 0, 0, NSUB, DUB)
      WRITE(6,'(A,I6)') ' DGNOMET: number of connections:',NTREE
      QDOIT=NTREE.LT.(NSUB*NSUB)/5.OR.SSHORT.EQ.'SPAR'
      END IF
C
      IF (QDOIT) THEN
      WRITE(6,'(A)')
     & ' DGNOMET: using sparse matrix bound smoothing algorithm.'
      TREE=ALLHP(INTEG4(NTREE))
      TREEPNT=ALLHP(INTEG4(NSUB+1))
      PATHPTR=ALLHP(IREAL8(NSUB*2+1))
      LOWLBL=ALLHP(INTEG4(NSUB))
      UPLBL=ALLHP(INTEG4(NSUB))
      LOWIND=ALLHP(INTEG4(NSUB))
      UPIND=ALLHP(INTEG4(NSUB))
C
C store connections in TREE, TREEPNT
      CALL DGTREE(.FALSE., NTREE, HEAP(TREEPNT), HEAP(TREE), NSUB, DUB)
C
C partial matrix bound-smoothing.
      UCOUNT=0
      LCOUNT=0
      DO H = 1, NSUB
C
C tighten the bounds matrixes only for columns that weren't touched
C during metrization (note: the following if test is well-defined
C even for NSLCT=0 as the whole SCLT array has been initialized to 0)
C
      IF (SLCT(SLCTCNT) .NE. H) THEN
      CALL DGSHRT2(NSUB,DUB,HEAP(PATHPTR),HEAP(LOWLBL),HEAP(UPLBL),
     &             HEAP(LOWIND),HEAP(UPIND),
     &             H,OLDROOT,HEAP(TREEPNT),HEAP(TREE),UCOUNT,LCOUNT)
      OLDROOT = H
      ELSE
      SLCTCNT = MIN(SLCTCNT + 1, NSLCT)
      END IF
      END DO
      WRITE(6,'(A,I5,I5)')
     &  ' DGNOMET: max. size of network (upper, lower): ',UCOUNT,LCOUNT
      CALL FREHP(UPIND,INTEG4(NSUB))
      CALL FREHP(LOWIND,INTEG4(NSUB))
      CALL FREHP(UPLBL,INTEG4(NSUB))
      CALL FREHP(LOWLBL,INTEG4(NSUB))
      CALL FREHP(PATHPTR, IREAL8(NSUB*2+1))
      CALL FREHP(TREEPNT,INTEG4(NSUB+1))
      CALL FREHP(TREE,INTEG4(NTREE))
      ELSE
      WRITE(6,'(A)')
     & ' DGNOMET: using full matrix bound smoothing algorithm.'
      PATHPTR=ALLHP(IREAL8(NSUB*2+1))
      LBLPTHPTR=ALLHP(IREAL8(NSUB*2+1))
      OLDROOT = 0
C
C full-matrix bound-smoothing.
      DO H = 1, NSUB
C
C tighten the bounds matrixes only for columns that weren't touched
C during metrization (note: the following if test is well-defined
C even for NSLCT=0 as the whole SCLT array has been initialized to 0)
C
      IF (SLCT(SLCTCNT) .NE. H) THEN
      CALL DGSHRTP(NSUB,DUB,HEAP(PATHPTR),HEAP(LBLPTHPTR),
     &             H,OLDROOT)
      OLDROOT = H
      ELSE
      SLCTCNT = MIN(SLCTCNT + 1, NSLCT)
      END IF
C
      END DO
      CALL FREHP(PATHPTR, IREAL8(NSUB*2+1))
      CALL FREHP(LBLPTHPTR, IREAL8(NSUB*2+1))
      END IF
      END IF
C
C
C randomly pick a distance between lower and upper bounds
C and set both bounds to the picked value (only for distances
C that wern't touched during metrization).
      SLCTCNT = 1
      DO H = 1, NSUB
      IF (SLCT(SLCTCNT) .NE. H) THEN
C
C include repulsive information "on-the-fly" for all
C lower bounds from H to any other atom.
      IF (QENER(SSVDW)) THEN
      CALL DGSETNN(H, NSUB, SUB, DUB, IBLOEX, INBEX)
      END IF
C
C
      DO J = 1, NSUB
      HH=MIN(H,J)
      JJ=MAX(H,J)
      UP=DUB(HH,JJ)
      LOW=DUB(JJ,HH)
      IF (ABS(LOW - UP) .GT. R4SMAL) THEN
      NUMPAIR = NUMPAIR + 1
      CALL GGUBFS(RANDOM)
C
      IF (UP .LT. LOW) THEN
      OK = .FALSE.
      TEMP = UP + ((LOW - UP) * RANDOM)
      WRITE(6, '(6(A4,1X),3(F15.4,1X),A)') SEGID(H), RESID(H),
     &TYPE(H), SEGID(J), RESID(J), TYPE(J), UP,TEMP,LOW,'*'
      ELSE
      TEMP = LOW + ((UP - LOW) * RANDOM)
      IF (DODUMP) THEN
      WRITE(6, '(6(A4,1X),3(F15.4,1X))') SEGID(H), RESID(H),
     &TYPE(H), SEGID(J), RESID(J), TYPE(J), UP,TEMP,LOW
      END IF
      END IF
C
      GAP = GAP + DUB(HH,JJ) - DUB(JJ,HH)
      DUB(HH,JJ)=TEMP
      DUB(JJ,HH)=TEMP
C
      END IF
      END DO
      ELSE
      SLCTCNT = MIN(SLCTCNT + 1, NSLCT)
      END IF
      END DO
C
      IF (.NOT. OK) THEN
      WRITE(6, '(2A)')
     &'* Geometric inconsistency--lower bound greater than upper.',
     &' Attempted to continue by inverting these bounds. '
      CALL WRNDIE(+1,'DISGEO',
     &'There are geometric inconsistencies. Proceed at own risk.')
      END IF
C
      IF (NUMPAIR .NE. 0) GAP = GAP / NUMPAIR
      CALL DECLAR('NON_MET_GAP', 'DP', ' ', DBCOMP, GAP)
      RETURN
      END
C =========================
      SUBROUTINE DGSMOOT(SSHORT, NSUB, DUB)
C
C Performs bounds-smoothing on upper and lower bounds matrices.
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
C passed-in variables
      CHARACTER*4 SSHORT
      INTEGER NSUB
      REAL DUB(NSUB,*)
C local variables
      INTEGER H, OLDROOT, NTREE, UCOUNT, LCOUNT
      LOGICAL QDOIT
C pointer
      INTEGER PATHPTR, LBLPTHPTR, TREE, TREEPNT, LOWLBL, UPLBL
      INTEGER LOWIND, UPIND
C begin
C
      QDOIT=.FALSE.
      IF (SSHORT.EQ.'AUTO'.OR.SSHORT.EQ.'SPAR') THEN
      CALL DGTREE(.TRUE., NTREE, 0, 0, NSUB, DUB)
      WRITE(6,'(A,I6)') ' DGSMOOT: number of connections:',NTREE
      QDOIT=NTREE.LT.(NSUB*NSUB)/5.OR.SSHORT.EQ.'SPAR'
      END IF
C
      IF (QDOIT) THEN
      WRITE(6,'(A)')
     & ' DGSMOOT: using sparse matrix bound smoothing algorithm.'
      TREE=ALLHP(INTEG4(NTREE))
      TREEPNT=ALLHP(INTEG4(NSUB+1))
      PATHPTR=ALLHP(IREAL8(NSUB*2+1))
      LOWLBL=ALLHP(INTEG4(NSUB))
      UPLBL=ALLHP(INTEG4(NSUB))
      LOWIND=ALLHP(INTEG4(NSUB))
      UPIND=ALLHP(INTEG4(NSUB))
C
C store connections in TREE, TREEPNT
      CALL DGTREE(.FALSE., NTREE, HEAP(TREEPNT), HEAP(TREE), NSUB, DUB)
      OLDROOT = 0
C
C bound-smoothing.
C go through all atoms ph in the molecule in order
      UCOUNT=0
      LCOUNT=0
      DO H = 1, NSUB
      CALL DGSHRT2(NSUB,DUB,HEAP(PATHPTR),HEAP(LOWLBL),HEAP(UPLBL),
     &             HEAP(LOWIND),HEAP(UPIND),
     &             H,OLDROOT,HEAP(TREEPNT),HEAP(TREE),UCOUNT,LCOUNT)
      OLDROOT = H
      END DO
      WRITE(6,'(A,I5,I5)')
     &  ' DGSMOOT: max. size of network (upper, lower): ',UCOUNT,LCOUNT
      CALL FREHP(UPIND,INTEG4(NSUB))
      CALL FREHP(LOWIND,INTEG4(NSUB))
      CALL FREHP(UPLBL,INTEG4(NSUB))
      CALL FREHP(LOWLBL,INTEG4(NSUB))
      CALL FREHP(PATHPTR, IREAL8(NSUB*2+1))
      CALL FREHP(TREEPNT,INTEG4(NSUB+1))
      CALL FREHP(TREE,INTEG4(NTREE))
      ELSE
      WRITE(6,'(A)')
     & ' DGSMOOT: using full matrix bound smoothing algorithm.'
      PATHPTR=ALLHP(IREAL8(NSUB*2+1))
      LBLPTHPTR=ALLHP(IREAL8(NSUB*2+1))
      OLDROOT = 0
C
C bound-smoothing.
C go through all atoms ph in the molecule in order
      DO H = 1, NSUB
      CALL DGSHRTP(NSUB,DUB,HEAP(PATHPTR),HEAP(LBLPTHPTR),
     &             H,OLDROOT)
      OLDROOT = H
      END DO
      CALL FREHP(PATHPTR, IREAL8(NSUB*2+1))
      CALL FREHP(LBLPTHPTR, IREAL8(NSUB*2+1))
      END IF
      RETURN
      END
C =========================
      SUBROUTINE DGEMBED(NATOM, DUB, TDUB, VECS, EIGVALS, NSUB, SUB,
     &                   DODUMP, QEMBED, XPCTED)
C
C Embedding routine.
C Calulates 3-D coordinates from a metric
C matrix TDUB stored in triangular form.   Note,
C that this assumes that the eigenvalues list
C returned by GIVEIS is in increasing order. Included substructures
C (ATB).
C adds non-embedded points by projection into subspace (MN).
C
C by John Kuszewski 8/24/89 12/29/89
C =========================
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'consta.inc'
C passed-in variables
      INTEGER NATOM
      INTEGER NSUB
      REAL DUB(NATOM,*)
      DOUBLE PRECISION TDUB(*)
      DOUBLE PRECISION VECS(NSUB,*), EIGVALS(NSUB)
      INTEGER SUB(*)
      LOGICAL DODUMP, QEMBED
      DOUBLE PRECISION XPCTED
C local variables
      INTEGER J, C, IERR
      DOUBLE PRECISION ACTUAL,XC,YC,ZC
      INTEGER LASTSET
      INTEGER LOW,UP
      DOUBLE PRECISION RANDOM,TEMP
C pointer
      INTEGER HPBMAT, HPIND
C parameter and dummies
      DOUBLE COMPLEX DBCOMP
      DOUBLE PRECISION ZERO, HALF
      PARAMETER (ZERO=0.0D0, HALF=0.5D0)
C begin
C
C store diagonal of MM for later
      DO C=1,NATOM
      RMSD(C)=0
      END DO
      DO C=1,NSUB
      RMSD(SUB(C))=TDUB(C*(C-1)/2+C)
      END DO
C
C allocate workspace for the EISPACK routine
      HPBMAT=ALLHP(IREAL8(8*NSUB))
      HPIND =ALLHP(INTEG4(NSUB))
C
C call the EISPACK diagonalization routine for square matrices (ATB)
      CALL GIVEIS(NSUB,NSUB,NSUB,TDUB,HEAP(HPBMAT),HEAP(HPIND),
     &            EIGVALS,VECS,IERR)
      IF (IERR.NE.0) THEN
        WRITE(6,'(A,I3)')
     &' %GIVEIS-ERR: diagonalization failed, IERR= ', IERR
      END IF
C
C free-up temporary space for EISPACK routine
      CALL FREHP(HPIND, INTEG4(NSUB))
      CALL FREHP(HPBMAT,IREAL8(8*NSUB))
C
C check if embedding was successful
      IF ((EIGVALS(NSUB) .LT. 0) .OR. (EIGVALS(NSUB-1) .LT. 0)
     &.OR. (EIGVALS(NSUB-2) .LT. 0).OR.IERR.NE.0) THEN
      QEMBED = .FALSE.
      CALL WRNDIE(+1,'DISGEO',
     &'Embedding failed. Try again with different seed.')
      ELSE
      QEMBED = .TRUE.
C
C to avoid stupid machine-dependencies we make sure that
C the first component of each eigenvector is always positive (ATB).
      IF (VECS(1,NSUB).LT.ZERO) THEN
      DO C=1,NSUB
      VECS(C,NSUB)=-VECS(C,NSUB)
      END DO
      END IF
      IF (VECS(1,NSUB-1).LT.ZERO) THEN
      DO C=1,NSUB
      VECS(C,NSUB-1)=-VECS(C,NSUB-1)
      END DO
      END IF
      IF (VECS(1,NSUB-2).LT.ZERO) THEN
      DO C=1,NSUB
      VECS(C,NSUB-2)=-VECS(C,NSUB-2)
      END DO
      END IF
      DO C=1, NSUB
      X(SUB(C))=SQRT(EIGVALS(NSUB))*VECS(C,NSUB)
      Y(SUB(C))=SQRT(EIGVALS(NSUB-1))*VECS(C,NSUB-1)
      Z(SUB(C))=SQRT(EIGVALS(NSUB-2))*VECS(C,NSUB-2)
      END DO
C
C add non-embedded points by projection (MN)
      IF (NSUB.LT.NATOM) THEN
      LASTSET=1
      DO C=1,NATOM
      IF (C.EQ.SUB(LASTSET)) THEN
      LASTSET=MIN(NSUB,LASTSET+1)
      ELSE
      XC=0
      YC=0
      ZC=0
C
      DO J=1,NSUB
      LOW=DUB(C,SUB(J))
      UP =DUB(SUB(J),C)
      IF (ABS(LOW - UP) .GT. R4SMAL) THEN
      CALL GGUBFS(RANDOM)
      IF (UP .LT. LOW) THEN
      TEMP = UP + ((LOW - UP) * RANDOM)
      ELSE
      TEMP = LOW + ((UP - LOW) * RANDOM)
      END IF
      ELSE
      TEMP = LOW
      END IF
C
      XC=XC+X(SUB(J))*(RMSD(SUB(J))-TEMP**2)
     &  /(2*EIGVALS(NSUB))
      YC=YC+Y(SUB(J))*(RMSD(SUB(J))-TEMP**2)
     &  /(2*EIGVALS(NSUB-1))
      ZC=ZC+Z(SUB(J))*(RMSD(SUB(J))-TEMP**2)
     &  /(2*EIGVALS(NSUB-2))
      END DO
      X(C)=XC
      Y(C)=YC
      Z(C)=ZC
      END IF
      END DO
      END IF
C
C now compute the actual radius of gyration and compute
C a "recommended" scale factor
      ACTUAL = ZERO
      DO J = 1, NSUB
      ACTUAL = ACTUAL + X(SUB(J))**2 + Y(SUB(J))**2 +
     &                  Z(SUB(J))**2
      END DO
      ACTUAL = SQRT(ACTUAL / NSUB)
C
C store the value of the scale factor in the symbol DGSCALE
      CALL DECLAR('DGSCALE', 'DP', ' ', DBCOMP, XPCTED / ACTUAL)
C
      DO C=1,NSUB
      RMSD(SUB(C))=EIGVALS(C)
      END DO
C
      END IF
      END
C =========================
      SUBROUTINE DGMMMM(NSUB, DUB, TDUB, RTEMP, XPCTED)
C
C Computes metric matrix from distance matrix.  Stores result
C in triangular matrix TDUB (ATB).
C
C by John Kuszewski
C =================
C
      IMPLICIT NONE
C input/output
C passed-in variables
      INTEGER NSUB
      REAL DUB(NSUB,*)
      DOUBLE PRECISION TDUB(*), RTEMP(*), XPCTED
C local variables
      DOUBLE PRECISION SUM2, SUM1, METRIC, R
      INTEGER J, K, H
C parameters
      DOUBLE PRECISION ZERO, HALF
      PARAMETER (ZERO=0.0D0, HALF=0.5D0)
C begin
C
C first, compute metric matrix from distance matrix DUB
      SUM2 = ZERO
      DO J = 1, NSUB
      DO K = J, NSUB
      R=DUB(K,J)
      SUM2=SUM2+R**2
      END DO
      END DO
      SUM2 = SUM2 / (NSUB**2)
      DO H=1, NSUB
      SUM1=0
      DO J=1, NSUB
      R=DUB(J,H)
      SUM1=SUM1+R**2
      END DO
C
C store in (H,H) element
      RTEMP(H)=(SUM1 / NSUB) - SUM2
      TDUB(H*(H+1)/2)=RTEMP(H)
      END DO
      DO H=1, NSUB
      DO J=H+1, NSUB
      R=DUB(J,H)
      METRIC = HALF*(RTEMP(H) + RTEMP(J) - R**2)
C
C store in J,H element
      TDUB(J*(J-1)/2+H)=METRIC
      END DO
      END DO
C
C compute the expected radius of gyration based on the
C metric matrix
      XPCTED = ZERO
      DO J = 1, NSUB
      XPCTED = XPCTED + RTEMP(J)
      END DO
      XPCTED = SQRT(XPCTED / NSUB)
      RETURN
      END
C =========================
      SUBROUTINE DGREDUC(NATOM, NSUB, DUB, DUBS,
     &                   SUB, NSELCT,SELCT)
C
C Projects full bounds matrices into substructure space.
C Projects atom selection pointer array SELCT into substructure
C space.
C
C Author: Axel T. Brunger
C =======================
      IMPLICIT NONE
C input/output
C passed-in variables
      INTEGER NATOM, NSUB
      REAL DUB(NATOM,*)
      INTEGER SUB(*), NSELCT, SELCT(*)
      REAL DUBS(NSUB,*)
C local variables
      INTEGER H, J
C begin
      DO H=1, NSUB
      DO J=1, NSUB
      DUBS(J,H)=DUB(SUB(J),SUB(H))
      END DO
      END DO
C
      NSELCT=0
      DO H=1,NSUB
      IF (SELCT(SUB(H)).EQ.1) THEN
      NSELCT=NSELCT+1
      SELCT(H)=1
      ELSE
      SELCT(H)=0
      END IF
      END DO
      RETURN
      END
C =========================
      SUBROUTINE DGINIT
C
C Routine initializes the bounds matrix pointers
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'mmdisg.inc'
C begin
C initialize the bounds matrix
      DUBPTR=0
      DUBLEN=0
C initialize the DG energy term atom selection
      DGSLEN=0
      DGSELE=0
      NDGSELE=0
C default scale factor and exponent for DG energy term
      DGSCALE=1.0D0
      DGEXPO=2
C initialize the flag that indicates that VDW lower bounds
C have been added to the bounds matrix (only used for EMMDG at
C present)
      DGQVDW=0
      RETURN
      END
C =========================
      SUBROUTINE DGRESZ(L)
C
C Routine resizes the bounds matrix allocation:
C positive entry: free existing space if appropriate and allocate
C                 new space.
C zero entry: free existing space if appropriate.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mmdisg.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'heap.inc'
      INTEGER L
C begin
C
C free up space for the old bounds matrix
      IF (DUBPTR.NE.0) THEN
      CALL FREHP(DUBPTR,IREAL4(DUBLEN))
      DUBPTR=0
      DUBLEN=0
      END IF
C free up space for the old DG energy term selection (only
C if L=0, i.e., the clean-up call)
      IF (L.EQ.0.AND.DGSLEN.GT.0) THEN
      CALL FREHP(DGSELE,INTEG4(DGSLEN))
      DGSELE=0
      DGSLEN=0
      NDGSELE=0
      END IF
C
      IF (L.GT.0) THEN
      DUBPTR=ALLHP(IREAL4(L))
      DUBLEN=L
      END IF
C
      RETURN
      END
C =========================
      SUBROUTINE DGREAD(BUNIT,NSUB,DUB,DGQVDW)
C
C routine reads bounds matrices from file
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER BUNIT, NSUB
      REAL DUB(NSUB,NSUB)
      INTEGER DGQVDW
C local
      INTEGER I,J
C begin
      DO I=1,NSUB
      READ(BUNIT) (DUB(I,J),J=1,NSUB)
      END DO
      READ(BUNIT) DGQVDW
      RETURN
      END
C =========================
      SUBROUTINE DGWRIT(BUNIT,NSUB,DUB,DGQVDW)
C
C routine writes bounds matrices to file
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INTEGER BUNIT, NSUB
      REAL DUB(NSUB,NSUB)
      INTEGER DGQVDW
C local
      INTEGER I, J
C begin
      DO I=1,NSUB
      WRITE(BUNIT) (DUB(I,J),J=1,NSUB)
      END DO
      WRITE(BUNIT) DGQVDW
      RETURN
      END
C =========================
      INTEGER FUNCTION DGMNDX (N, X)
C
C Finds the lowest index i such that
C x(i) = min( x(j): j = 1, 2, .., n)
C x may be a one-dimensional array or a row or
C column of a two-dimensional array.
C
C Designed for easy vectorization.
C
C by John Kuszewski 8/11/91
C =========================
      IMPLICIT NONE
C passed-in variables
      INTEGER N
      DOUBLE PRECISION X(*)
C local variables
      INTEGER COUNT
      LOGICAL DONE
      DOUBLE PRECISION CURMIN
C begin
      CURMIN = X(1)
      DO COUNT = 1,N
      CURMIN = MIN(CURMIN,X(COUNT))
      END DO
      COUNT = 1
      DONE = .FALSE.
      DO WHILE (.NOT. DONE)
      IF (CURMIN .EQ. X(COUNT)) THEN
      DGMNDX = COUNT
      DONE = .TRUE.
      ELSE
      COUNT = COUNT + 1
      END IF
      END DO
      RETURN
      END
C =========================
      SUBROUTINE EMMDG(E)
C
C Distance geometry energy routine that allows one to
C restrain the coordinates to the bounds matrix stored
C in memory (STOREBounds option). Front-end routine for
C EMMDG2.
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mmdisg.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'ener.inc'
      DOUBLE PRECISION E
C local
      INTEGER ROOT, MM, NSELE, NPAIR
C pointer
      INTEGER SUBPTR, EXCPTR, IBLOEX, INBEX, IPK, JPK
C begin
      IF (DUBLEN.NE.NATOM**2.OR.DUBPTR.EQ.0) THEN
      CALL WRNDIE(-5,
     &    'EMMDG',
     &    'Bounds matrix not well-defined (cf. STOREBounds in MMDG)')
      ELSEIF (NDGSELE.EQ.0) THEN
      CALL WRNDIE(-5,
     &     'EMMDG',
     &    'No atoms selected for MMDG energy (cf. SELEction in MMDG)')
      ELSE
C
C we have to include the vdw information "on-the-fly" if required.
      IF (DGQVDW.EQ.-1) THEN
      DGQVDW=1
      WRITE(6,'(A)')
     & ' EMMDG-info: vdw bounds are added to stored bounds matrix!'
C make sure lookup tables for van der Waals parameters are updated
      IF (UPNBLK) THEN
      CALL NBUPDA
      UPNBLK=.FALSE.
      END IF
C
C first we obtain the required number of exclusions without
C actually storing them
      CALL DGSETNX(NPAIR,0,0,0,0,.TRUE.)
C
C now we allocate the required space and then call the routine
C again to store the exclusions in IBLOEX, INBEX
      IPK=ALLHP(INTEG4(NPAIR))
      JPK=ALLHP(INTEG4(NPAIR))
      IBLOEX=ALLHP(INTEG4(NATOM+1))
      INBEX=ALLHP(INTEG4(NPAIR))
      CALL DGSETNX(NPAIR,HEAP(IPK),HEAP(JPK),HEAP(IBLOEX),HEAP(INBEX),
     &             .FALSE.)
      CALL FREHP(IPK,INTEG4(NPAIR))
      CALL FREHP(JPK,INTEG4(NPAIR))
      EXCPTR=ALLHP(INTEG4(NATOM))
      SUBPTR=ALLHP(INTEG4(NATOM))
      CALL FILL4(HEAP(SUBPTR),NATOM,1)
      CALL MAKIND(HEAP(SUBPTR),NATOM,MM)
      DO NSELE=1,NPIG
      DO ROOT=1,NATOM
      CALL DGSETN2(ROOT, NATOM, HEAP(SUBPTR),
     &            HEAP(DUBPTR), HEAP(EXCPTR),
     &     HEAP(IBLOEX),HEAP(INBEX), HEAP(IINTER(NSELE)))
      END DO
      END DO
      CALL FREHP(IBLOEX,INTEG4(NATOM+1))
      CALL FREHP(INBEX,INTEG4(NPAIR))
      CALL FREHP(SUBPTR,INTEG4(NATOM))
      CALL FREHP(EXCPTR,INTEG4(NATOM))
      END IF
C
      CALL EMMDG2(E,NATOM,DGSCALE,DGEXPO,NDGSELE,HEAP(DUBPTR),
     &            HEAP(DGSELE))
      END IF
      RETURN
      END
C =========================
      SUBROUTINE EMMDG2(E,NATOM,DGSCALE,DGEXPO,NDGSELE,DUB,
     &                  DGSELE)
C
C See routine EMMDG above.
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'coord.inc'
      DOUBLE PRECISION E
      INTEGER NATOM
      DOUBLE PRECISION DGSCALE
      INTEGER DGEXPO, NDGSELE
      REAL DUB(NATOM,*)
      INTEGER DGSELE(*)
C local
      INTEGER IAT, JAT, IIAT, JJAT
      DOUBLE PRECISION UP, LOW, XIJ, YIJ, ZIJ, R, DEVI, DF
      DOUBLE PRECISION DXTMP, DYTMP, DZTMP
C parameter
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C begin
      E=ZERO
      DO IAT=1,NDGSELE
C
C CRI BCB, lowercase are changes.  this vectorizes the loop which now
C runs up to 4 times faster (or more)
C
      IIAT=DGSELE(IAT)
      DXTMP = 0.0
      DYTMP = 0.0
      DZTMP = 0.0
CDIR$ IVDEP
      DO JAT=1,IAT-1
      JJAT=DGSELE(JAT)
      UP=DUB(JJAT,IIAT)
      LOW=DUB(IIAT,JJAT)
      XIJ=X(IIAT)-X(JJAT)
      YIJ=Y(IIAT)-Y(JJAT)
      ZIJ=Z(IIAT)-Z(JJAT)
      R=SQRT(XIJ**2+YIJ**2+ZIJ**2)
      DEVI=MAX(ZERO,R-UP)+MIN(ZERO,R-LOW)
      E=E+DGSCALE*(DEVI**DGEXPO)
      DF=(DGEXPO*DGSCALE*(DEVI**(DGEXPO-1)))/R
      DX(JJAT)=DX(JJAT)-DF*XIJ
      DY(JJAT)=DY(JJAT)-DF*YIJ
      DZ(JJAT)=DZ(JJAT)-DF*ZIJ
      DXTMP = DXTMP + DF*XIJ
      DYTMP = DYTMP + DF*YIJ
      DZTMP = DZTMP + DF*ZIJ
      END DO
      DX(IIAT)=DX(IIAT)+DXTMP
      DY(IIAT)=DY(IIAT)+DYTMP
      DZ(IIAT)=DZ(IIAT)+DZTMP
      END DO
      RETURN
      END
