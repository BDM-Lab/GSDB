      SUBROUTINE DYNMCS
C
C Molecular dynamics
C
C Authors: Luke Rice and Axel T. Brunger.
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C begin
      CALL NEXTWD('DYNAmics-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-dynamics')
C
      ELSE IF (WD(1:4).EQ.'TORS' ) THEN
      CALL TORMD
      ELSE IF (WD(1:4).EQ.'CART' ) THEN
      CALL DCART
      ELSE IF (WD(1:4).EQ.'MERG') THEN
      CALL TRJMERGE
      ELSE
      CALL DSPERR('DYNAMICS','unknown qualifier')
      END IF
      RETURN
      END
C===============================================================
      SUBROUTINE TORMD
C
C Authors: Luke Rice and Axel T. Brunger
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'dtorsion.inc'
C
C local
      INTEGER HPIBXCL,HPJBXCL, HPXCLBND
      INTEGER HPFIXSEL2,HPFIXSEL3,HPFRESEL2,HPFRESEL3,HPGPSEL
      INTEGER HPFLAGS, HPFLAGS2, HPCRDIND
C begin
C
      IF (QRESET) THEN
C
C If QRESET is true this means that the torsion topology is empty and 
C has be to be generated.
C
C The following two arrays are stored in the torsion heap storage area.
C they will stay around until the torsion topology has been reset
      HPLIST=ALLHP(INTEG4(NATOM))
      HPPOINTR=ALLHP(INTEG4(NATOM+2))
      END IF
C
C the following arrays are temporary arrays for TORMD2
      HPFLAGS=ALLHP(INTEG4(NATOM))
C
      HPFLAGS2=ALLHP(INTEG4(NATOM))
      HPCRDIND=ALLHP(INTEG4(NATOM))
      HPFIXSEL2=ALLHP(INTEG4(NATOM))
      HPFIXSEL3=ALLHP(INTEG4(NATOM))
      HPFRESEL2=ALLHP(INTEG4(NATOM))
      HPFRESEL3=ALLHP(INTEG4(NATOM))
      HPGPSEL=ALLHP(INTEG4(NATOM))
C
      HPIBXCL=ALLHP(INTEG4(NBOND))
      HPJBXCL=ALLHP(INTEG4(NBOND))
      HPXCLBND = ALLHP(INTEG4(NBOND))
C
      CALL TORMD2(HEAP(HPFLAGS),HEAP(HPFLAGS2),HEAP(HPCRDIND),
     &            HEAP(HPLIST),HEAP(HPPOINTR),
     &            HEAP(HPFIXSEL2),HEAP(HPFIXSEL3),HEAP(HPFRESEL2),
     &            HEAP(HPFRESEL3),HEAP(HPGPSEL),HEAP(HPIBXCL),
     &            HEAP(HPJBXCL),HEAP(HPXCLBND))
C
      CALL FREHP(HPXCLBND,INTEG4(NBOND))
      CALL FREHP(HPJBXCL,INTEG4(NBOND))
      CALL FREHP(HPIBXCL,INTEG4(NBOND))
C
      CALL FREHP(HPGPSEL,INTEG4(NATOM))
      CALL FREHP(HPFRESEL3,INTEG4(NATOM))
      CALL FREHP(HPFRESEL2,INTEG4(NATOM))
      CALL FREHP(HPFIXSEL3,INTEG4(NATOM))
      CALL FREHP(HPFIXSEL2,INTEG4(NATOM))
      CALL FREHP(HPCRDIND,INTEG4(NATOM))
      CALL FREHP(HPFLAGS2,INTEG4(NATOM))
C
      CALL FREHP(HPFLAGS,INTEG4(NATOM))
C
C QRESET=true means that a RESET command was given in TORMD2- reset the torsion heap
C space storage, including the LIST and POINTR arrays. 
      IF (QRESET) THEN
      CALL TORHP(0,.TRUE.)
      CALL TORINI
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE TORMD2(FLAGS,FLAGS2,CRDIND,
     &                  LIST,POINTR,FIXSEL2,FIXSEL3,FRESEL2,
     &                  FRESEL3,GPSEL,IBXCL,JBXCL,XCLBND)
C
C Parser routine for rigid body torsion dynamics simulation.
C
C Authors: Luke Rice and Axel T. Brunger
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'ener.inc'
      INCLUDE 'timer.inc'
      INTEGER FLAGS(NATOM), FLAGS2(NATOM), CRDIND(NATOM)
      INTEGER LIST(NATOM), POINTR(NATOM+2)
      INTEGER IBXCL(NBOND), JBXCL(NBOND), XCLBND(NBOND)
      INTEGER FIXSEL2(NATOM), FIXSEL3(NATOM)
      INTEGER FRESEL2(NATOM), FRESEL3(NATOM)
      INTEGER GPSEL(NATOM)
C
C local
      INTEGER HPTREE
      INTEGER NDIHE, HPTJP, HPTKP
      INTEGER NDHXCL, HPJDHXCL, HPKDHXCL
      INTEGER NDHINC, HPJDHINC, HPKDHINC
      INTEGER NFIX, NFRE, NGPSEL
      INTEGER NBNDXCL
      INTEGER HPJATOMS,HPJBONDS,HPCONGRP
      INTEGER HPNBND, HPCONATM, HPNDIH, HPCONDIH
      INTEGER HPTAKEN, HPGROUP, HPGPLIST
      INTEGER HPNUMBT, HPBONDT
      INTEGER NSAVC, IUNCRD, CRDNUM
      INTEGER NSTEP, I, II, NFLAGS, NPRINT, NTRFRQ
      INTEGER NUNIQ, J, CMPERIOD
      INTEGER HPXCMO, HPYCMO, HPZCMO, HPRKCRD, HPRKVEL
      INTEGER NTREES,MXCHAINS,MXCHNLEN
      DOUBLE PRECISION SCALE, OFFSET, DT
      DOUBLE PRECISION ANUM,ZERO, SECS
      CHARACTER*(COMMAX) CFILE
      CHARACTER*10 FORM
      LOGICAL QFORM, QFIRSC
      LOGICAL DEFAUL, QDEBUG, UDEF, DONE2
      LOGICAL ERR
      LOGICAL TCOUPL, VSCALE, QUNIQ, CMREMOVE
      LOGICAL QLOOP
      LOGICAL RESETCOMMAND
C
      PARAMETER (ANUM=9999.0D0,ZERO=0.0D0)
C begin
C
      RESETCOMMAND=.FALSE.
C
      IF (QRESET) THEN
C
C topology will be generated 
         QREDO=.TRUE.
         NGPSEL = 0
         NFIX = 0
         NFRE = 0
C
C reset bond exclusion list (start with all bonds included)
         NBNDXCL = 0
         DO I=1,NBOND
         XCLBND(I)=0
         END DO
C
C
         KDIHMAX = 100.0D0
         MAXCHN=-1
         MAXLEN=-1
         MAXTREE=-1
C make a default selection consisting of all free atoms, this will
C overwritten by GETPNTR
         NGP=1
         POINTR(1)=0
         POINTR(2)=0
         DO I=1,NATOM
            IF (IMOVE(I).EQ.0) THEN
               POINTR(2)=POINTR(2)+1
               LIST(POINTR(2))=I
            END IF
         END DO
      END IF
C
C determine maximum number of bonds to any atom
      DO I=1,NATOM
      FLAGS(I)=0
      END DO
      DO I=1,NBOND
      FLAGS(IB(I))=FLAGS(IB(I))+1
      FLAGS(JB(I))=FLAGS(JB(I))+1
      END DO
      MAXBND=0
      DO I=1,NATOM
      MAXBND=MAX(MAXBND,FLAGS(I))
      END DO
C
      NSTEP=0
C
      TCOUPL = .FALSE.
      VSCALE = .FALSE.
C
      QFORM=.FALSE.
      OFFSET=800.0D0
      SCALE=10000.0D0
      FORM='12Z6'
C
      CFILE=' '
      NSAVC=0
      IUNCRD=0
C
      QFIRSC=.TRUE.
C
      CRDNUM=0
      DO I=1,NATOM
      CRDNUM=CRDNUM+1
      CRDIND(I)=1
      END DO
C
      QDEBUG=.FALSE.
      DEFAUL=.TRUE.
      DT=1.0D-3
      TBATH=298.0D0
C
      NPRINT=1
C
      CMPERIOD=0
      CMREMOVE=.TRUE.
C
      CALL PUSEND('Torsion Dynamics>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('Torsion Dynamics>')
      CALL MISCOM('Torsion Dynamics>',USED)
      IF (.NOT.USED) THEN
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-dynamics-torsion')
C
C=====================================================================
      ELSE IF (WD(1:4).EQ.'?   ') THEN
      WRITE(6,'(2A)')
     &' -------------------------------------',
     &'--------------------------------------'
      WRITE(6,'(A)') ' Torsion Dynamics: '
      WRITE(6,'(A,I9,A,F12.7)')
     &  '   NSTEp=',NSTEP,'   TIMEstep=',DT
      WRITE(6,'(A,L1)')'   Initial COM removal: CMREmove = ',CMREMOVE
      WRITE(6,'(A,I6)')'   Periodic COM removal: CMPEriodic = ',CMPERIOD
      IF (TCOUPL) THEN
      WRITE(6,'(A,F10.4)')
     & '   TCOUpling=true  TEMPerature= ',TBATH
      END IF
      IF (VSCALE) THEN
      WRITE(6,'(A,F10.4,A,I10,A)')
     & '   VSCAling=true  TEMPerature= ',TBATH
      END IF
      WRITE(6,'(A,A)')
     & '   Trajectory  Info:  TRAJectory-file = ', CFILE
      WRITE(6,'(A,I6,A)')
     & '                      is written every ',NSAVC,' steps'
      WRITE(6,'(2A)')
     &' -------------------------------------',
     &'--------------------------------------'
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TOPO') THEN
C
C only parse topology statements if it does not exist already (i.e., QRESET=TRUE)
         IF (QRESET) THEN
         CALL PUSEND('Torsion Topology>')
         DONE2=.FALSE.
         DO WHILE (.NOT. DONE2)
            CALL NEXTWD('Torsion Topology>')
            CALL MISCOM('Torsion Topology>',USED)
            IF (.NOT.USED) THEN
               IF (WD(1:4).EQ.'HELP') THEN
                  CALL CNSHELP('cns-dynamics-torsion-topology')
               ELSE IF (WD(1:4).EQ.'FIX') THEN
                  CALL NEXTWD('FIX-qualifier=')
                  IF (WD(1:4).EQ.'TORS') THEN
                     ERR=.FALSE.
                     NFIX = NFIX + 1
                     DO I=1,4
                        CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
                        IF (NFLAGS.NE.1) THEN
                           WRITE(6,'(A)')
     &  ' %TORSION-ERR: selection has to contain exactly one atom.'
                           ERR=.TRUE.
                        END IF
                        CALL MAKIND(FLAGS,NATOM,NFLAGS)
                        IF (I.EQ.2) THEN
                           FIXSEL2(NFIX) = FLAGS(1)
                        ELSE IF (I.EQ.3) THEN
                           FIXSEL3(NFIX) = FLAGS(1)
                        END IF
                     END DO
                     IF (ERR) NFIX = NFIX - 1
                  ELSE IF (WD(1:4).EQ.'GROU') THEN
                     CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
                     CALL MAKIND(FLAGS,NATOM,NFLAGS)
                     nuniq = 0
                     DO I=1,NFLAGS
                        quniq = .true.
                        do j=1,ngpsel
                           if (flags(i).eq.gpsel(j)) quniq = .false.
                        end do
                        if (quniq) then
                           nuniq = nuniq + 1
                           GPSEL(NGPSEL+nuniq)=FLAGS(I)
                        end if
                     END DO
                     NGPSEL = NGPSEL + nuniq
                  END IF
C=====================================================================
               ELSE IF (WD(1:4).EQ.'FREE') THEN
                  CALL NEXTWD('FREE-qualifier=')
                  IF (WD(1:4).EQ.'TORS') THEN
                     ERR=.FALSE.
                     NFRE = NFRE + 1
                     DO I=1,4
                        CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
                        IF (NFLAGS.NE.1) THEN
                           WRITE(6,'(A)')
     &  ' %TORSION-ERR: selection has to contain exactly one atom.'
                           ERR=.TRUE.
                        END IF
                        CALL MAKIND(FLAGS,NATOM,NFLAGS)
                        IF (I.EQ.2) THEN
                           FRESEL2(NFRE) = FLAGS(1)
                        ELSE IF (I.EQ.3) THEN
                           FRESEL3(NFRE) = FLAGS(1)
                        END IF
                     END DO
                     IF (ERR) NFRE = NFRE - 1
                  ELSE IF (WD(1:4).EQ.'BOND') THEN
                     ERR=.FALSE.
                     CALL SELCTA(FLAGS,NFLAGS,X,Y,Z,.TRUE.)
                     IF (NFLAGS.EQ.0.AND.WRNLEV.GE.10) THEN
                           WRITE(6,'(A)')
     &' %TORSION-info: no atoms in 1st selection for FREE BOND.'
                           ERR=.TRUE.
                     END IF
                     CALL SELCTA(FLAGS2,NFLAGS,X,Y,Z,.TRUE.)
                     IF (NFLAGS.EQ.0.AND.WRNLEV.GE.10) THEN
                           WRITE(6,'(A)')
     &' %TORSION-info: no atoms in 2nd selection for FREE BOND.'
                           ERR=.TRUE.
                     END IF
                     IF (.NOT.ERR) THEN
C
C lookup the freed bonds and mark as excluded
                      DO I=1,NBOND
                       IF (((FLAGS(IB(I)).EQ.1.AND.FLAGS2(JB(I)).EQ.1)
     &                           .OR.
     &                      (FLAGS(JB(I)).EQ.1.AND.FLAGS2(IB(I)).EQ.1))
     &                     .AND.XCLBND(I).EQ.0) THEN
                          NBNDXCL=NBNDXCL+1
                          IF (NBNDXCL.GT.NBOND) THEN
                            CALL WRNDIE(-5,'TORSION',
     @                'NBNDXCL exceeds NBONDS. Fatal coding error')
                          END IF

                          IBXCL(NBNDXCL)=IB(I)
                          JBXCL(NBNDXCL)=JB(I)
                          XCLBND(I)=1
                       END IF
                      END DO
C
C
                     END IF
                  ELSE
                   CALL DSPERR('TORMD2','unknown FREE qualifier')
                  END IF
C=====================================================================
               ELSE IF (WD(1:4).EQ.'KDIH') THEN
                  CALL NEXTF('KDIHmax=',KDIHMAX)
C=====================================================================
C this is usually automatic
               ELSE IF (WD(1:4).EQ.'MAXT') THEN
                  CALL NEXTI('MAXTree=',MAXTREE)
C=====================================================================
C this is usually automatic
               ELSE IF (WD(1:4).EQ.'MAXC') THEN
                  CALL NEXTI('MAXChain=',MAXCHN)
C=====================================================================
C this is usually automatic
               ELSE IF (WD(1:4).EQ.'MAXL') THEN
                  CALL NEXTI('MAXLength=',MAXLEN)
C=====================================================================
               ELSE IF (WD(1:4).EQ.'MAXD') THEN
                  CALL NEXTI('MAXDihe=',I)
                  WRITE(6,'(A)') ' Obsolete option: MAXDihe. Ignored'
C=====================================================================
               ELSE IF (WD(1:4).EQ.'MAXJ') THEN
                  CALL NEXTI('MAXJoint=',i)
                  WRITE(6,'(A)') ' Obsolete option: MAXJnt. Ignored'
C=====================================================================
               ELSE IF (WD(1:4).EQ.'MAXB') THEN
                  CALL NEXTI('MAXBond=',I)
                  WRITE(6,'(A)') ' Obsolete option: MAXBond. Ignored'
C=====================================================================
               ELSE
                  CALL CHKEND('Torsion Topology>',DONE2)
               END IF
            END IF
         END DO
         DONE2=.FALSE.
         ELSE
C
         CALL PUSEND('Torsion Topology>')
         DONE2=.FALSE.
         DO WHILE (.NOT. DONE2)
            CALL NEXTWD('Torsion Topology>')
            IF (WD(1:4).EQ.'HELP') THEN
               CALL CNSHELP('cns-dynamics-torsion-topology')
            ELSE IF (WD(1:4).EQ.'RESE') THEN
               QRESET=.TRUE.
               RESETCOMMAND=.TRUE.
            ELSE IF (WD(1:3).NE.'END') THEN
         WRITE(6,'(A,A)')
     @        ' +++++++++++++++++++++++++++++ Torsion Topology +++',
     @        '++++++++++++++++++++++++++++'
         WRITE(6,'(A)')
     @        '   WARNING: The previous topology is inherited. '
         WRITE(6,'(A)')
     @        '            Topology statements will not be parsed.  '
         WRITE(6,'(A)')
     @     ' +++++++++++++++++++++++++++++++++++++++',
     @     '+++++++++++++++++++++++++++++++++++++++'
            ELSE
               CALL CHKEND('Torsion Topology>',DONE2)
            END IF
         END DO
         DONE2=.FALSE.
         END IF
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TIME'.OR.WD(1:4).EQ.'DT') THEN
      CALL NEXTF('TIMEstep=',DT)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NSTE') THEN
      CALL NEXTI('NSTEp=',NSTEP)
C=====================================================================
      ELSE IF (WD(1:5).EQ.'NSAVC') THEN
      CALL NEXTI('NSAVC=',NSAVC)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTST('FORMat=',FORM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'SCAL') THEN
      CALL NEXTF('SCALe=',SCALE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'OFFS') THEN
      CALL NEXTF('OFFSet=',OFFSET)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TRAJ') THEN
      CALL NEXTFI('TRAJectory-file=',CFILE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CSEL') THEN
      CALL SELCTA(CRDIND,CRDNUM,X,Y,Z,.TRUE.)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'ASCI') THEN
      CALL NEXTLO('ASCIi=',QFORM)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'REAS') THEN
      CALL NEXTLO('REASsign=',QREDO)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CMRE') THEN
      CALL NEXTLO('CMREmove=',CMREMOVE)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'CMPE') THEN
      CALL NEXTI('CMPEriodic=',CMPERIOD)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NTRF') THEN
      CALL NEXTI('NTRFrq=',NTRFRQ)
      IF (NTRFRQ.LE.0) CMREMOVE=.FALSE.
      IF (NTRFRQ.GE.1) CMREMOVE=.TRUE.
C=====================================================================
      ELSE IF (WD(1:4).EQ.'NPRI') THEN
      CALL NEXTI('NPRINt=',NPRINT)
C=====================================================================
      ELSE IF (WD(1:4).EQ.'TCOU') THEN
      CALL NEXTLO('TCOUpling=',TCOUPL)
      IF (TCOUPL) THEN
      WRITE(6,'(A)')
     &   ' TORMD: temperature coupling (TCOUpl) enabled'
      IF (VSCALE) THEN
      VSCALE=.FALSE.
      WRITE(6,'(A)')
     &   ' TORMD: velocity rescaling (VSCAle) disabled'
      END IF
      END IF
      ELSE IF (WD(1:4).EQ.'VSCA') THEN
      CALL NEXTLO('VSCAling=',VSCALE)
      IF (VSCALE) THEN
      WRITE(6,'(A)')
     &   ' TORMD: velocity rescaling (VSCAle) enabled'
      IF (TCOUPL) THEN
      TCOUPL=.FALSE.
      WRITE(6,'(A)')
     &   ' TORMD: temperature coupling (TCOUpl) disabled'
      END IF
      END IF
      ELSE IF (WD(1:4).EQ.'TBAT'.OR.WD(1:4).EQ.'TEMP') THEN
      CALL NEXTF('TEMPerature=',TBATH)
C=====================================================================
      ELSE
      CALL CHKEND('Torsion Dynamics>',DONE)
      END IF
      END IF
      END DO
      DONE=.FALSE.
C
      IF (.NOT.RESETCOMMAND) THEN
C
      IF (POINTR(2).EQ.0) THEN
      CALL WRNDIE(-5,'TORMD2',
     & 'All atoms are fixed.  Molecular dynamics not performed.')
      ELSE 
C
C open trajectory files
      IF (CFILE.NE.' ') THEN
      IF (QFORM) THEN
      CALL ASSFIL(CFILE,IUNCRD,'WRITE','FORMATTED',ERROR)
      ELSE
      CALL ASSFIL(CFILE,IUNCRD,'WRITE','UNFORMATTED',ERROR)
      END IF
      END IF
C
C set up the trajectory index list
      CALL MAKIND(CRDIND,NATOM,CRDNUM)
C
C check that masses are well-defined
      UDEF=.FALSE.
      DO I=1,NATOM
      IF (AMASS(I).GT.ANUM-1) UDEF=.TRUE.
      END DO
      IF (UDEF) THEN
      CALL WRNDIE(-5,'TORMD2',
     & 'there are undefined or very large (> 9998 amu) masses.')
      END IF
C
      IF (QRESET) THEN
C
C generate topology
         HPTJP = ALLHP(INTEG4(MAXP))
         HPTKP = ALLHP(INTEG4(MAXP))
         HPNBND = ALLHP(INTEG4(NATOM))
         HPCONATM = ALLHP(INTEG4(NATOM*MAXBND))
         HPJDHXCL = ALLHP(INTEG4(MAXP))
         HPKDHXCL = ALLHP(INTEG4(MAXP))
         HPJDHINC = ALLHP(INTEG4(MAXP))
         HPKDHINC = ALLHP(INTEG4(MAXP))
         HPNDIH = ALLHP(INTEG4(NATOM))
         HPGPLIST = ALLHP(INTEG4(NATOM*MAXBND))
         HPGROUP = ALLHP(INTEG4(NATOM))
         HPNUMBT = ALLHP(INTEG4(NATOM))
         HPBONDT = ALLHP(INTEG4(NATOM*MAXBND))
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to GETDIH',SECS  
         END IF 
         CALL GETDIH(NDIHE,HEAP(HPTJP),HEAP(HPTKP),
     @               HEAP(HPNUMBT),HEAP(HPBONDT))
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to GETCON',SECS  
         END IF 
         CALL GETCON(HEAP(HPNBND),HEAP(HPCONATM),NBNDXCL,
     @                IBXCL,JBXCL,XCLBND)
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to FILLXCL',SECS  
         END IF 
         CALL FILLXCL(FIXSEL2,FIXSEL3,NFIX,NATOM,HEAP(HPCONATM),
     @                HEAP(HPNBND),FRESEL2,FRESEL3,NFRE,GPSEL,NGPSEL,
     @                HEAP(HPJDHXCL),HEAP(HPKDHXCL),NDHXCL,
     @                HEAP(HPJDHINC),HEAP(HPKDHINC),NDHINC)
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to XCLDIH',SECS  
         END IF 
         CALL XCLDIH(NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @               NDIHE,HEAP(HPTJP),HEAP(HPTKP),
     @               NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @               HEAP(HPNBND),HEAP(HPCONATM))
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to CNTDIH',SECS  
         END IF 
         CALL CNTDIH(NDIHE,HEAP(HPNDIH),HEAP(1),
     @               NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @               NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @               NBNDXCL,IBXCL,JBXCL,HEAP(HPNBND),
     @               HEAP(HPTJP),HEAP(HPTKP),.TRUE.)
         HPCONDIH = ALLHP(INTEG4(NATOM*MAXDIHE))
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to CNTDIH',SECS  
         END IF 
         CALL CNTDIH(NDIHE,HEAP(HPNDIH),HEAP(HPCONDIH),
     @               NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @               NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @               NBNDXCL,IBXCL,JBXCL,HEAP(HPNBND),
     @               HEAP(HPTJP),HEAP(HPTKP),.FALSE.)
         IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to GETPNTR',SECS  
         END IF 
         CALL GETPNTR(LIST,POINTR,HEAP(HPCONATM),
     @               HEAP(HPCONDIH),HEAP(HPNBND),HEAP(HPNDIH),
     @               HEAP(HPGPLIST),HEAP(HPGROUP))
C
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to DTORJNT',SECS  
         END IF 
         CALL DTORJNT(XCLBND,LIST,POINTR)
C
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to DTOREST',SECS  
         END IF 
         CALL DTOREST(XCLBND,IBXCL,JBXCL,NBNDXCL,LIST,POINTR,
     &                NTREES,MXCHAINS,MXCHNLEN,QLOOP)
C
         IF (QLOOP) THEN
C
C closed loops were detected and bonds excluded from topology list
C re-run the setup
          WRITE(6,'(A)')
     &' TORMD2: re-running topology setup due to closed loops.'
C
          CALL FREHP(HPCONDIH,INTEG4(NATOM*MAXDIHE))
C
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to GETDIH',SECS  
          END IF 
          CALL GETDIH(NDIHE,HEAP(HPTJP),HEAP(HPTKP),
     @                HEAP(HPNUMBT),HEAP(HPBONDT))
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to GETCON',SECS  
          END IF 
          CALL GETCON(HEAP(HPNBND),HEAP(HPCONATM),NBNDXCL,
     @                 IBXCL,JBXCL,XCLBND)
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to FILLXCL',SECS  
          END IF 
          CALL FILLXCL(FIXSEL2,FIXSEL3,NFIX,NATOM,HEAP(HPCONATM),
     @                 HEAP(HPNBND),FRESEL2,FRESEL3,NFRE,GPSEL,NGPSEL,
     @                 HEAP(HPJDHXCL),HEAP(HPKDHXCL),NDHXCL,
     @                 HEAP(HPJDHINC),HEAP(HPKDHINC),NDHINC)
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to XCLDIH',SECS  
          END IF 
          CALL XCLDIH(NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @                NDIHE,HEAP(HPTJP),HEAP(HPTKP),
     @                NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @                HEAP(HPNBND),HEAP(HPCONATM))
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to CNTDIH',SECS  
          END IF 
          CALL CNTDIH(NDIHE,HEAP(HPNDIH),HEAP(1),
     @                NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @                NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @                NBNDXCL,IBXCL,JBXCL,HEAP(HPNBND),
     @                HEAP(HPTJP),HEAP(HPTKP),.TRUE.)
          HPCONDIH = ALLHP(INTEG4(NATOM*MAXDIHE))
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to CNTDIH',SECS  
          END IF 
          CALL CNTDIH(NDIHE,HEAP(HPNDIH),HEAP(HPCONDIH),
     @                NDHXCL,HEAP(HPJDHXCL),HEAP(HPKDHXCL),
     @                NDHINC,HEAP(HPJDHINC),HEAP(HPKDHINC),
     @                NBNDXCL,IBXCL,JBXCL,HEAP(HPNBND),
     @                HEAP(HPTJP),HEAP(HPTKP),.FALSE.)
          IF (TIMER.GT.0) THEN
            CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to GETPNTR',SECS  
          END IF 
          CALL GETPNTR(LIST,POINTR,HEAP(HPCONATM),
     @                HEAP(HPCONDIH),HEAP(HPNBND),HEAP(HPNDIH),
     @                HEAP(HPGPLIST),HEAP(HPGROUP))
C
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to DTORJNT',SECS  
          END IF 
          CALL DTORJNT(XCLBND,LIST,POINTR)
C
          IF (TIMER.GT.0) THEN
           CALL VCPU(SECS)
           WRITE(6,'(A,F12.5)')' CPU time: before call to DTOREST',SECS  
          END IF 
C
C determine maximum dimensions for large matrices
          CALL DTOREST(XCLBND,IBXCL,JBXCL,NBNDXCL,LIST,POINTR,
     &                 NTREES,MXCHAINS,MXCHNLEN,QLOOP)
         END IF
C
C set the maximum dimensions to explicitly specified values
         IF (MAXTREE.LE.0) THEN
         MAXTREE=NTREES
         END IF
         IF (MAXLEN.LE.0) THEN
         MAXLEN=MXCHNLEN
         END IF
         IF (MAXCHN.LE.0) THEN
         MAXCHN=MXCHAINS
         END IF
C
         CALL FREHP(HPBONDT,INTEG4(NATOM*MAXBND))
         CALL FREHP(HPNUMBT,INTEG4(NATOM))
         CALL FREHP(HPGROUP,INTEG4(NATOM))
         CALL FREHP(HPGPLIST,INTEG4(NATOM*MAXBND))
         CALL FREHP(HPCONDIH,INTEG4(NATOM*MAXDIHE))
         CALL FREHP(HPNDIH,INTEG4(NATOM))
         CALL FREHP(HPKDHINC,INTEG4(MAXP))
         CALL FREHP(HPJDHINC,INTEG4(MAXP))
         CALL FREHP(HPKDHXCL,INTEG4(MAXP))
         CALL FREHP(HPJDHXCL,INTEG4(MAXP))
         CALL FREHP(HPCONATM,INTEG4(NATOM*MAXBND))
         CALL FREHP(HPNBND,INTEG4(NATOM))
         CALL FREHP(HPTKP,INTEG4(MAXP))
         CALL FREHP(HPTJP,INTEG4(MAXP))
C
C allocate space for large matrices based on MAXTREE, MAXLEN, and MAXCHN
         CALL TORHP(NATOM,.FALSE.)
C
         HPTREE = ALLHP(INTEG4(NGP))
         HPTAKEN = ALLHP(INTEG4(NGP))
         HPJATOMS = ALLHP(INTEG4(NGP*MAXJNT))
         HPJBONDS = ALLHP(INTEG4(NGP*MAXJNT))
         HPCONGRP = ALLHP(INTEG4(NGP*MAXJNT))
         IF (TIMER.GT.0) THEN
          CALL VCPU(SECS)
          WRITE(6,'(A,F12.5)')' CPU time: before call to TOPOLOGY',SECS  
         END IF 
         CALL TOPOLOGY(LIST,POINTR,HEAP(HPNJOINT),
     @         HEAP(HPGRPSEQ),HEAP(HPNIN),HEAP(HPNOUT),HEAP(HPSIZE),
     @         HEAP(HPCHNLEN),HEAP(HPTAKEN),HEAP(HPTREE),
     @         HEAP(HPNCHAINS),HEAP(HPJATOMS),HEAP(HPJBONDS),
     @         HEAP(HPCONGRP),IBXCL,JBXCL,NBNDXCL,XCLBND)
         IF (TIMER.GT.0) THEN
            CALL VCPU(SECS)
            WRITE(6,'(A,F12.5)')' CPU time: after call to TOPOLOGY',SECS  
         END IF 
         CALL FREHP(HPCONGRP,INTEG4(NGP*MAXJNT))
         CALL FREHP(HPJBONDS,INTEG4(NGP*MAXJNT))
         CALL FREHP(HPJATOMS,INTEG4(NGP*MAXJNT))
         CALL FREHP(HPTAKEN, INTEG4(NGP))
         CALL FREHP(HPTREE, INTEG4(NGP))
C
         CALL SETUP(HEAP(HPNCHAINS),
     @     HEAP(HPCHNLEN),
     @     HEAP(HPCOEFFI),HEAP(HPCOEFFO),HEAP(HPNJOINT),
     @     HEAP(HPSIZE),HEAP(HPGRPSEQ),
     @     HEAP(HPNIN),HEAP(HPNOUT),
     @     HEAP(HPQTIP))
C
      END IF
C
      WRITE(6,'(A,A)')
     @        ' -------------------------- Torsion Angle Dynamics ',
     @        '-----------------------------'
      IF (QRESET) THEN
      WRITE(6,'(A)') ' WARNING: Current Velocities Taken '
      WRITE(6,'(A)') '          Bad Initial Velocities May Result '
      ELSE IF (.NOT.QREDO) THEN
      WRITE(6,'(2A)')
     @ ' Info: Velocities and coordinates taken from ',
     & 'last torsion angle dynamics step.'
      ELSE
      WRITE(6,'(2A)')
     @ ' Info: Velocities and coordinates taken from ',
     & 'main coordinates (x,y,z) and velocities (vx, vy, vz).'
      END IF
      WRITE(6,'(A,A)')
     @     ' ---------------------------------------',
     @     '----------------------------------------'
C
      IF (QREDO) THEN
         CALL SETUP2(LIST,POINTR,HEAP(HPNCHAINS),HEAP(HPCHNLEN),
     @     HEAP(HPXCM),
     @     HEAP(HPYCM),HEAP(HPZCM),HEAP(HPCOEFFI),HEAP(HPCOEFFO),
     @     HEAP(HPMM),HEAP(HPYY),HEAP(HPQ0),HEAP(HPQ1),
     @     HEAP(HPQ2),HEAP(HPQ3),HEAP(HPSIZE),HEAP(HPGRPSEQ),HEAP(HPXB),
     @     HEAP(HPYB),HEAP(HPZB),HEAP(HPNIN),HEAP(HPNOUT),HEAP(HPQDOT),
     @     HEAP(HPHHX),HEAP(HPHHY),HEAP(HPHHZ),HEAP(HPSSX),HEAP(HPSSY),
     @     HEAP(HPSSZ),HEAP(HPQTIP))
         QREDO = .FALSE.
      END IF
C
      HPXCMO=ALLHP(IREAL8(MAXTREE))
      HPYCMO=ALLHP(IREAL8(MAXTREE))
      HPZCMO=ALLHP(IREAL8(MAXTREE))
      HPRKCRD=ALLHP(IREAL8(MAXTREE*7))
      HPRKVEL=ALLHP(IREAL8(MAXTREE*7))
      CALL TORMD3(NSTEP,LIST,POINTR,QFIRSC,IUNCRD,DT,CRDNUM,
     @     CRDIND,NSAVC,NPRINT,QFORM,FORM,SCALE,OFFSET,
     @     CMREMOVE,CMPERIOD,
     @     HEAP(HPGRPSEQ),HEAP(HPNIN),HEAP(HPNOUT),HEAP(HPSIZE),
     @     HEAP(HPXCM),HEAP(HPYCM),HEAP(HPZCM),HEAP(HPCOEFFI),
     @     HEAP(HPCOEFFO),HEAP(HPMM),HEAP(HPYY),HEAP(HPQ0),
     @     HEAP(HPQ1),HEAP(HPQ2),HEAP(HPQ3),HEAP(HPXB),HEAP(HPYB),
     @     HEAP(HPZB),HEAP(HPFRICT),HEAP(HPQQ),HEAP(HPB1),HEAP(HPB2),
     @     HEAP(HPDD),HEAP(HPKK),HEAP(HPLL),HEAP(HPYYDOT),HEAP(HPQ0O),
     @     HEAP(HPQ1O),HEAP(HPQ2O),HEAP(HPQ3O),HEAP(HPYYO),
     @     HEAP(HPRKQ),HEAP(HPRKQDOT),HEAP(HPCHNLEN),HEAP(HPQDOT),
     @     HEAP(HPQDDOT),HEAP(HPQDOTO),HEAP(HPQO),HEAP(HPQ),
     @     HEAP(HPHHX),HEAP(HPHHY),HEAP(HPHHZ),HEAP(HPSSX),
     @     HEAP(HPSSY),HEAP(HPSSZ),HEAP(HPNCHAINS),TCOUPL,VSCALE,
     @     HEAP(HPXCMO),
     @     HEAP(HPYCMO),HEAP(HPZCMO),HEAP(HPRKCRD),HEAP(HPRKVEL))
      CALL FREHP(HPRKVEL,IREAL8(MAXTREE*7))
      CALL FREHP(HPRKCRD,IREAL8(MAXTREE*7))
      CALL FREHP(HPZCMO,IREAL8(MAXTREE))
      CALL FREHP(HPYCMO,IREAL8(MAXTREE))
      CALL FREHP(HPXCMO,IREAL8(MAXTREE))
C
      QRESET=.FALSE.
C
      END IF
C
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE GETDIH(NDIHE,TJP,TKP,NUMBT,BONDT)
C
C Go through the bond list and find all dihedrals.
C Make a list based of the 2-3 atoms.
C
C Skip the dihedral if the 2nd or 3rd atom is fixed.
C
C      NDIHE    counts the total number of dihedrals
C      TJP, TKP store the 2nd and 3rd atoms
C
C Author: Luke Rice
C                   ( After MAXINB, AUTODI )
C modified: Paul Adams 10/98
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER NDIHE,TJP(*),TKP(*)
      INTEGER NUMBT(*)
      INTEGER BONDT(MAXBND,*)
C
C local
      INTEGER IDIHE,I,J,WJ,WK
      LOGICAL COND,QACPT
C
C initialize the dihedral list
      NDIHE = 0
      DO I=1,MAXP
      TJP(I) = 0
      TKP(I) = 0
      END DO
C
C initialize the temporary bond list
      DO I=1,NATOM
      NUMBT(I)=0
      DO J=1,MAXBND
      BONDT(J,I)=0
      END DO
      END DO
C
C generate temporary bond list indexed on atom number
      DO I=1,NBOND
      NUMBT(IB(I))=NUMBT(IB(I))+1
      NUMBT(JB(I))=NUMBT(JB(I))+1
      BONDT(NUMBT(IB(I)),IB(I))=JB(I)
      BONDT(NUMBT(JB(I)),JB(I))=IB(I)
      END DO
C
C go through the bond list to get angles
      DO I=1,NBOND
      WJ=IB(I)
      WK=JB(I)
      COND=.FALSE.
      IF ((NUMBT(WJ).GT.0).AND.(NUMBT(WK).GT.0)) THEN
         DO J=1,NUMBT(WJ)
            IF ((BONDT(J,WJ).NE.WK).AND.
     &          (IMOVE(BONDT(J,WJ)).EQ.0)) COND=.TRUE.
         END DO
         IF (COND) THEN
            COND=.FALSE.
            DO J=1,NUMBT(WK)
               IF ((BONDT(J,WK).NE.WJ).AND.
     &             (IMOVE(BONDT(J,WK)).EQ.0)) COND=.TRUE.
            END DO
         END IF
      END IF
C
C Put this in the (reduced) 2-3 list if it's not already there
      IF (COND) THEN
C
      IDIHE = 1
      QACPT = .TRUE.
      DO WHILE (IDIHE.LE.NDIHE.AND.QACPT)
      QACPT = ((WJ.NE.TJP(IDIHE).OR.
     &          WK.NE.TKP(IDIHE)).AND.
     &         (WJ.NE.TKP(IDIHE).OR.
     &          WK.NE.TJP(IDIHE)))
      IDIHE = IDIHE + 1
      END DO
C
      IF (QACPT) THEN
      NDIHE = NDIHE + 1
      TJP(NDIHE) = WJ
      TKP(NDIHE) = WK
      END IF
C
      END IF
C
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE GETCON(NBND,CONATM,NBNDXCL,IBXCL,JBXCL,XCLBND)
C
C get the number of bonds each atom is involved in
C     and the connected atoms
C
C     NBND        stores the number of bonds each atom is involved in
C     CONATM      lists (for each atom) the connected ones
C     IBXCL,JBXCL stores the atoms of the
C     NBNDXCL     bonds to be ignored for the topology
C     XCLBND      identifies excluded bonds
C
C Author: Luke Rice
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER NBND(*),CONATM(NATOM,*)
      INTEGER NBNDXCL,IBXCL(*),JBXCL(*),XCLBND(*)
C
C Local
      INTEGER I,J,II,JJ
      LOGICAL QACPT
C
C Begin
      DO I=1,NATOM
      NBND(I) = 0
      DO J=1,MAXBND
      CONATM(I,J) = 0
      END DO
      END DO
C
C add bonds to fixed atoms to the excluded bonds
      DO I=1,NBOND
      IF ( (IMOVE(IB(I)).EQ.1.OR.IMOVE(JB(I)).EQ.1)
     &     .AND.XCLBND(I).EQ.0) THEN
      XCLBND(I)=1
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'GETCON',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      IBXCL(NBNDXCL) = IB(I)
      JBXCL(NBNDXCL) = JB(I)
      END IF
      END DO
C
      DO I=1,NBOND
C
      IF (XCLBND(I).EQ.0) THEN
      II = IB(I)
      JJ = JB(I)
      NBND(II) = NBND(II)+1
      NBND(JJ) = NBND(JJ)+1
C
      IF (NBND(II).GT.MAXBND) THEN
      CALL WRNDIE(-5,'TORSION:GETCON',
     @     'MAXBND exceeded. Fatal coding error')
      END IF
      IF (NBND(JJ).GT.MAXBND) THEN
      CALL WRNDIE(-5,'TORSION:GETCON',
     @     'MAXBND exceeded. Fatal coding error')
      END IF
C
      CONATM(II,NBND(II)) = JJ
      CONATM(JJ,NBND(JJ)) = II
      END IF
C
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE DTORJNT(XCLBND,LIST,POINTR)
C
C determines the maximum number of joints (MAXJNT) to any group
C
C Author: Axel Brunger 5/19/08
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XCLBND(*),LIST(*),POINTR(*)
      INTEGER NTREES, MXCHAINS, MXCHNLEN
C local
      INTEGER MAXLEV
C pointers
      INTEGER GROUPID, NJOINT
C begin
      GROUPID=ALLHP(INTEG4(NATOM))
      NJOINT=ALLHP(INTEG4(NGP))
      CALL DTORJNT2(XCLBND,LIST,POINTR,
     &    HEAP(GROUPID),HEAP(NJOINT))
      CALL FREHP(NJOINT,INTEG4(NGP))
      CALL FREHP(GROUPID,INTEG4(NATOM))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE DTORJNT2(XCLBND,LIST,POINTR,
     &                   GROUPID,NJOINT)
C
C determines the maximum number of joints (MAXJNT) to any group
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER XCLBND(*),LIST(*),POINTR(*)
      INTEGER GROUPID(NATOM), NJOINT(NGP)
C local
      INTEGER I, POINTER, IBOND, POINTER2, NJOINTS
      LOGICAL COND
C begin
C
      MAXJNT=0
C initialize number of joints. 
      NJOINTS=0
      DO I=1,NGP
      NJOINT(I)=0
      END DO
C
C make array GROUPID(1...natom) that specifies the group index for each atom. 
      DO I=1,NGP
      DO POINTER=POINTR(I)+1,POINTR(I+1)
      GROUPID(LIST(POINTER))=I
      END DO
      END DO
C
C go through all included bonds
      DO IBOND=1,NBOND
      IF (XCLBND(IBOND).EQ.0) THEN
      DO I=1,NGP
C
C test if this bond is between distinct groups
      DO POINTER=POINTR(I)+1,POINTR(I+1)
C
C we try both combinations for IB and JB
      IF (IB(IBOND).EQ.LIST(POINTER)) THEN
      COND=.TRUE.
      DO POINTER2=POINTR(I)+1,POINTR(I+1)
      IF (JB(IBOND).EQ.LIST(POINTER2)) THEN
      COND=.FALSE.
      END IF
      END DO
      IF (COND) THEN
C
C this bond is between two distinct groups
      NJOINT(I)=NJOINT(I)+1
      NJOINT(GROUPID(JB(IBOND)))=NJOINT(GROUPID(JB(IBOND)))+1
      MAXJNT=MAX(MAXJNT,NJOINT(I),NJOINT(GROUPID(JB(IBOND))))
      END IF
      END IF
C
      END DO
      END DO
      END IF
      END DO
      RETURN
      END
C===============================================================
      SUBROUTINE DTOREST(XCLBND,IBXCL,JBXCL,NBNDXCL,LIST,POINTR,
     &                   NTREES,MXCHAINS,MXCHNLEN,QLOOP)
C
C estimate the dimensions needed for the topology generation
C
C Author: Axel Brunger 5/19/08
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INTEGER XCLBND(NBOND),IBXCL(NBOND),JBXCL(NBOND),NBNDXCL
      INTEGER LIST(NATOM),POINTR(NATOM+2)
      INTEGER NTREES, MXCHAINS, MXCHNLEN
      LOGICAL QLOOP
C local
      INTEGER MAXLEV
C pointers
      INTEGER RECLIST, VERTEX, IND, RETADR
      INTEGER GRPCONNCHN, GRPCONNIDX
      INTEGER GROUPID, NJOINT, GRPCONNGRP, GRPCONNBND
C begin
      MAXLEV=NBOND+1
      RECLIST=ALLHP(INTEG4(NGP))
      VERTEX=ALLHP(INTEG4(MAXLEV))
      IND=ALLHP(INTEG4(MAXLEV))
      RETADR=ALLHP(INTEG4(MAXLEV))
      GRPCONNCHN=ALLHP(INTEG4(NBOND))
      GROUPID=ALLHP(INTEG4(NATOM))
      NJOINT=ALLHP(INTEG4(NGP))
      GRPCONNGRP=ALLHP(INTEG4(NGP*MAXJNT))
      GRPCONNIDX=ALLHP(INTEG4(NGP*MAXJNT))
      GRPCONNBND=ALLHP(INTEG4(NBOND))
      CALL DTOREST2(XCLBND,IBXCL,JBXCL,NBNDXCL,LIST,POINTR,
     &    HEAP(RECLIST),MAXLEV,HEAP(VERTEX),HEAP(IND),HEAP(RETADR),
     &    HEAP(GRPCONNCHN),
     &    HEAP(GROUPID),HEAP(NJOINT),HEAP(GRPCONNGRP),
     &    HEAP(GRPCONNIDX),HEAP(GRPCONNBND),
     &    NTREES,MXCHAINS,MXCHNLEN,QLOOP)
      CALL FREHP(GRPCONNBND,INTEG4(NBOND))
      CALL FREHP(GRPCONNIDX,INTEG4(NGP*MAXJNT))
      CALL FREHP(GRPCONNGRP,INTEG4(NGP*MAXJNT))
      CALL FREHP(NJOINT,INTEG4(NGP))
      CALL FREHP(GROUPID,INTEG4(NATOM))
      CALL FREHP(GRPCONNCHN,INTEG4(NBOND))
      CALL FREHP(RETADR,INTEG4(MAXLEV))
      CALL FREHP(IND,INTEG4(MAXLEV))
      CALL FREHP(VERTEX,INTEG4(MAXLEV))
      CALL FREHP(RECLIST,INTEG4(NGP))
C
      RETURN
      END
C=====================================================================
      SUBROUTINE DTOREST2(XCLBND,IBXCL,JBXCL,NBNDXCL,LIST,POINTR,
     &                   RECLIST,MAXLEV,VERTEX,IND,RETADR,
     &                   GRPCONNCHN,
     &                   GROUPID,NJOINT,
     &                   GRPCONNGRP,GRPCONNIDX,GRPCONNBND,
     &                   NTREES,MXCHAINS,MXCHNLEN,QLOOP)
C
C Author: Axel T. Brunger
C
      IMPLICIT NONE
C input/output
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER XCLBND(NBOND),IBXCL(NBOND),JBXCL(NBOND),NBNDXCL
      INTEGER LIST(NATOM),POINTR(NATOM+2)
      INTEGER RECLIST(NGP), MAXLEV, VERTEX(MAXLEV)
      INTEGER IND(MAXLEV), RETADR(MAXLEV)
      INTEGER GRPCONNCHN(NBOND)
      INTEGER GROUPID(NATOM), NJOINT(NGP)
      INTEGER GRPCONNGRP(NGP,MAXJNT), GRPCONNIDX(NGP,MAXJNT)
      INTEGER GRPCONNBND(NBOND)
      INTEGER NTREES, MXCHAINS, MXCHNLEN
      LOGICAL QLOOP
C local
      INTEGER I, II, LENGTH, LEVEL, PASSED, NEXT
      INTEGER POINTER, POINTER2, IBOND
      INTEGER NJOINTS, NCHAINS, NTCHAINS, NLOOPS
      LOGICAL COND
      DOUBLE PRECISION SECS
C parameter
      INTEGER MARK
      PARAMETER (MARK=-1)
C begin
C
C notes for input information:
C NGP is the number of groups
C POINTR, LIST is the linked list that defines the groups. 
C XCLBND(1...,nbond) is a flag indicating if a bond is included
C
C dimensions:
C MAXJNT is the maximum number joints connected to any group
C MAXLEV is the maximum number of levels during recursion.  It
C        cannot be larger than the number of covalent bonds.  
C        Hence it is set to NBOND in the calling routine. 
C 
C definitions:
C
C RECLIST(i) array indicating that a group i has already been visited
C VERTEX(level) index of group that is considered at particular recursion level
C IND(level) index of joints connected to group that is considered at particular recursion level
C RETADR(level) return address for recursion level
C
C NJOINTS is total number of joints
C NJOINT(i) is number of joints connected to group i
C GRPCONNGRP(i,1...n) is the indices of the groups 1...n connected to group i 
C GRPCONNIDX(i,1...n) are the indices of the joints between group i and groups 1...n
C GRPCONNBND(j) is the bond index for joint j
C GRPCONNCHN(j) is the chain index for joint j.  
C               It will remain zero if this joint would create a closed loop
C
C initialize number of joints. 
      NJOINTS=0
      DO I=1,NGP
      NJOINT(I)=0
      END DO
C
C make array GROUPID(1...natom) that specifies the group index for each atom. 
      DO I=1,NGP
      DO POINTER=POINTR(I)+1,POINTR(I+1)
      GROUPID(LIST(POINTER))=I
      END DO
      END DO
C
C go through all included bonds
      DO IBOND=1,NBOND
      IF (XCLBND(IBOND).EQ.0) THEN
      DO I=1,NGP
C
C test if this bond is between distinct groups
      DO POINTER=POINTR(I)+1,POINTR(I+1)
C
C we try both combinations for IB and JB
      IF (IB(IBOND).EQ.LIST(POINTER)) THEN
      COND=.TRUE.
      DO POINTER2=POINTR(I)+1,POINTR(I+1)
      IF (JB(IBOND).EQ.LIST(POINTER2)) THEN
      COND=.FALSE.
      END IF
      END DO
      IF (COND) THEN
C
C this bond is between two distinct groups
      NJOINT(I)=NJOINT(I)+1
      NJOINT(GROUPID(JB(IBOND)))=NJOINT(GROUPID(JB(IBOND)))+1
      NJOINTS=NJOINTS+1
C
      GRPCONNGRP(I,NJOINT(I))=GROUPID(JB(IBOND))
      GRPCONNGRP(GROUPID(JB(IBOND)),NJOINT(GROUPID(JB(IBOND))))=I
      GRPCONNIDX(I,NJOINT(I))=NJOINTS
      GRPCONNIDX(GROUPID(JB(IBOND)),NJOINT(GROUPID(JB(IBOND))))=NJOINTS
      GRPCONNBND(NJOINTS)=IBOND
      END IF
      END IF
C
      END DO
      END DO
      END IF
      END DO
C
C total number of distinct chains (NCHAINS)
C maximum number of distinct chains per tree (MXCHAINS)
C number of independent trees (NTREES)
      NCHAINS=0
      MXCHAINS=0
      NTREES=0
C
      DO I=1,NJOINTS
      GRPCONNCHN(I)=0
      END DO

C
C do recursive search over graph of NGP groups connected by joints 
C
C initalize indicator that tells if a group has already been visited
      DO I=1,NGP
      RECLIST(I)=MARK
      END DO
C
C
      LEVEL=1
C
C * *loop head: FOR {I=1,NGP}  ( loop over all groups)
      I=1
6363  IF (I.GT.NGP) GOTO 6364
C
      IF (RECLIST(I).NE.MARK) GOTO 6365
      NTREES=NTREES+1
      NTCHAINS=0
C
C *invoke PROCEDURE SEARCH (PASSED:=VERTEX)
      PASSED=I
      RETADR(LEVEL)=1
      GOTO 2000
2001  CONTINUE
C *return address
C
C
6365  CONTINUE
      I=I+1
      GOTO 6363
6364  CONTINUE
C
C * *go to exit
      GOTO 2999
C
C---- PROCEDURE SEARCH ( CALL BY VALUE: PASSED->VERTEX ) --------------
C     LOCAL: IND
C *head
2000  CONTINUE
      LEVEL=LEVEL+1
      IF (LEVEL.GT.MAXLEV) CALL WRNDIE(-2,'<DTOREST2>',
     @  'Coding error: MAXLEV exceeded')
      VERTEX(LEVEL)=PASSED
C
C *begin
      RECLIST(VERTEX(LEVEL))=NTREES
C
C * *loop head: FOR { IND(LEVEL)=1,NJOINT(VERTEX(LEVEL) }
C         (loop over all joints connected to group VERTEX(LEVEL))
      IND(LEVEL)=0
7711  CONTINUE
C * *loop iteration
      IND(LEVEL)=IND(LEVEL)+1
      IF (IND(LEVEL).GT.NJOINT(VERTEX(LEVEL))) GOTO 7722
C * *loop begin
C
      NEXT=GRPCONNGRP(VERTEX(LEVEL),IND(LEVEL))
      IF (RECLIST(NEXT).NE.MARK.AND.RECLIST(NEXT).NE.NTREES) THEN
      CALL WRNDIE(0,'<DTOREST2>','check algorithm')
      END IF
C
C go to end of loop if joint NEXT has already been visited (there
C are two situations: (1) we already traversed this joint in the reverse
C direction, or (2) we did not traverse it in the reverse direction
C but the joint NEXT was already visited: this means that we found
C a closed loop ):
      IF (RECLIST(NEXT).NE.MARK) GOTO 3567
C
C otherwise we found a joint that has not been visited.
C
      IF (NJOINT(VERTEX(LEVEL)).EQ.2) THEN
C special case of a group with two joints.   This joint is part of a chain.
      IF (IND(LEVEL).EQ.1) THEN
      IF (GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),2)).NE.0) THEN
C joint 2 has already been visited and a chain assigned.  
C assign the same chain id to the current joint. 
      GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),1))=
     & GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),2))
      ELSE
C start of a new chain
      NCHAINS=NCHAINS+1
      NTCHAINS=NTCHAINS+1
      MXCHAINS=MAX(NTCHAINS,MXCHAINS)
      GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),1))=NCHAINS
      END IF
      ELSE
C by definition joint 1 has already been visited.  Assign its
C chain id to the current joint
      GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),2))=
     & GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),1))
      END IF
C
      ELSE 
C
C start of a new chain
      NCHAINS=NCHAINS+1
      NTCHAINS=NTCHAINS+1
      MXCHAINS=MAX(NTCHAINS,MXCHAINS)
      GRPCONNCHN(GRPCONNIDX(VERTEX(LEVEL),IND(LEVEL)))=NCHAINS
      END IF
C
C * * *invoke procedure SEARCH ( PASSED:=NEXT )
      PASSED=NEXT
      RETADR(LEVEL)=2
      GOTO 2000
2002  CONTINUE
C * * *return label
C
CCCCC     END IF
3567  CONTINUE
C
C * *end loop
      GOTO 7711
7722  CONTINUE
C * *exit loop
C
C *return to address
      LEVEL=LEVEL-1
      IF (LEVEL.LE.0) CALL WRNDIE(-1,'<DTOREST2>','Level underflow')
      IF (RETADR(LEVEL).EQ.2) GOTO 2002
      IF (RETADR(LEVEL).EQ.1) GOTO 2001
      CALL WRNDIE(-2,'<DTOREST2>','Unknown return address, check code')
C--END PROCEDURE SEARCH-----------------------------------------------
C
2999  CONTINUE
C
C determine chain lengths and maximum length of any chain, 
C and determine number of closed loops
      QLOOP=.FALSE.
      MXCHNLEN=0
      NLOOPS=0
      DO I=1,NJOINTS
C
C if GRPCONNCHN(I) is zero this means that this is the last joint of
C a loop that was traversed 
       IF (GRPCONNCHN(I).EQ.0) THEN
         QLOOP=.TRUE.
         NLOOPS=NLOOPS+1
C
C free the bond that is the last link of that loop
      IBOND=GRPCONNBND(I)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'DTOREST2',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=IB(IBOND)
      JBXCL(NBNDXCL)=JB(IBOND)
C
      WRITE(6,'(17A)')
     @' Loop detected: bond "', SEGID(IB(IBOND)),'-',
     @ RESID(IB(IBOND)),'-',RES(IB(IBOND)),
     @ '-',TYPE(IB(IBOND)),'" and "',
     @ SEGID(JB(IBOND)),'-',RESID(JB(IBOND)),
     @ '-',RES(JB(IBOND)),'-',TYPE(JB(IBOND)),'" unfixed.'
C
       END IF
       ELSE
C
C otherwise, this group is part of a chain that does not form a loop.
C We now determine the length of this chain:
         LENGTH=1
         DO II=I+1,NJOINTS
            IF (GRPCONNCHN(I).EQ.GRPCONNCHN(II)) THEN
               LENGTH=LENGTH+1
            END IF
         END DO
         MXCHNLEN=MAX(MXCHNLEN,LENGTH)
C
       END IF
C
      END DO
C
C summary statements
      WRITE(6,'(A,I10)')
     &' DTOREST2: total number of groups ',NGP
      WRITE(6,'(A,I10)')
     &' DTOREST2: total number of joints ',NJOINTS
      WRITE(6,'(A,I10)')
     &' DTOREST2: total number of closed loops ',NLOOPS
      WRITE(6,'(A,I10)')
     &' DTOREST2: total number of distinct bodies ',NTREES
      WRITE(6,'(A,I10)')
     &' DTOREST2: total number of chains ',NCHAINS
      WRITE(6,'(A,I10)')
     &' DTOREST2: maximum number of chains (per distinct body) ',
     & MXCHAINS
      WRITE(6,'(A,I10)')
     &' DTOREST2: maximum unbranched chain length ',MXCHNLEN
C
C adjustments
      IF (.NOT.QLOOP) THEN
C
C increasing MXCHAINS by NLOOPS and 2
      MXCHAINS=MXCHAINS+NLOOPS+2
C
C increasing MXCHNLEN by 2 
      MXCHNLEN=MXCHNLEN+2
C
      WRITE(6,'(A,I10)')
     &' DTOREST2: adj. max. number of chains (per distinct body) ',
     & MXCHAINS
      WRITE(6,'(A,I10)')
     &' DTOREST2: adj. maximum unbranched chain length ',
     & MXCHNLEN
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE FILLXCL(FIXSEL2,FIXSEL3,NFIX,NATOM,CONATM,NBND,
     @     FRESEL2,FRESEL3,NFRE,GPSEL,NGPSEL,
     @     JDHXCL,KDHXCL,NDHXCL,JDHINC,KDHINC,NDHINC)
C
C Fill the exclusion lists for use in topology setup:
C
C      GPSEL            stores the NGPSEL user specified
C      FIXSEL2, FIXSEL3 store the 2nd and 3rd atoms of the
C               NFIX    user specified fixed dihedrals
C      FRESEL2, FRESEL3 store the 2nd and 3rd atoms of the
C               NFRE    user specified fixed dihedrals
C      JDHXCL,KDHXCL    store the 2nd and 3rd atoms of the
C             NDHXCL    explicitly fixed dihedrals
C      JDHINC,KDHINC    store the 2nd and 3rd atoms of the
C             NDHXCL    explicitly freed dihedrals
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'timer.inc'
      INTEGER NATOM
      INTEGER FIXSEL2(*),FIXSEL3(*),FRESEL2(*),FRESEL3(*)
      INTEGER JDHXCL(*),KDHXCL(*),CONATM(NATOM,*)
      INTEGER JDHINC(*),KDHINC(*),GPSEL(*),NBND(*)
      INTEGER NFIX,NFRE,NGPSEL,NDHXCL,NDHINC
C
C Local
      INTEGER I,ISEL,J,K,KSEL,II,JJ
      LOGICAL QBOMB
C
C Initialize
      NDHXCL = 0
      NDHINC = 0
C
C Process group list first
      DO I=1,NGPSEL
      ISEL = GPSEL(I)
      DO J=1,NBND(ISEL)
      DO K=1,NGPSEL
      KSEL = GPSEL(K)
      IF (CONATM(ISEL,J).EQ.KSEL) THEN
      NDHXCL = NDHXCL + 1
      JDHXCL(NDHXCL) = ISEL
      KDHXCL(NDHXCL) = KSEL
      END IF
      END DO
      END DO
      END DO
C
C Now add in the explicit exclusions
      DO I=1,NFIX
      II = FIXSEL2(I)
      JJ = FIXSEL3(I)
      QBOMB = .TRUE.
      DO J=1,NBND(II)
      IF (CONATM(II,J).EQ.JJ) QBOMB=.FALSE.
      END DO
      IF (QBOMB) THEN
      CALL DIE
      ELSE
      NDHXCL = NDHXCL + 1
      JDHXCL(NDHXCL) = II
      KDHXCL(NDHXCL) = JJ
      END IF
      END DO
C
C Process the include list and check connectivity
      DO I=1,NFRE
      II = FRESEL2(I)
      JJ = FRESEL3(I)
      QBOMB = .TRUE.
      DO J=1,NBND(II)
      IF (CONATM(II,J).EQ.JJ) QBOMB=.FALSE.
      END DO
      IF (QBOMB) THEN
      CALL DIE
      ELSE
      NDHINC = NDHINC + 1
      JDHINC(NDHINC) = II
      KDHINC(NDHINC) = JJ
      END IF
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE XCLDIH(NDHXCL,JDHXCL,KDHXCL,NDIHE,TJP,TKP,
     @                  NDHINC,JDHINC,KDHINC,NBND,CONATM)
C
C make a list of exclusions , and add them to the explicit ones -
C   the list is based on dihedrals with ( force constant > kdihmax )
C   and on 1-2-3-4 connected impropers with ( force constant > kdihmax )
C   also allow for user-specified exclusions/inclusions
C
C Author: Luke Rice
C modified: Paul Adams 10/98
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER NDIHE,TJP(*),TKP(*)
      INTEGER NDHXCL,JDHXCL(*),KDHXCL(*)
      INTEGER NDHINC,JDHINC(*),KDHINC(*)
      INTEGER NBND(*),CONATM(NATOM,*)
C Local
      INTEGER I,J,K,IXCL,IINC
      INTEGER JXCL,KXCL,JINC,KINC
      DOUBLE PRECISION ZERO
      LOGICAL QIJ,QJK,QKL,QEXCL,Q1234,QXTRA
      PARAMETER (ZERO = 0.0D0)
C
C first retrieve atom based parameters
      CALL CODDIH(NPHI,IP,JP,KP,LP,ACPC,ACPB,ACPD,
     @            QACPC,QACPB,QACPD,IAC,SEGID,RESID,TYPE)
      CALL CODIMP(NIMPHI,IM,JM,KM,LM,ACIC,ACIB,ACID,
     @            QACIC,QACIB,QACID,IAC,SEGID,RESID,TYPE)
C
C check to see if there are any unknown dihedrals
      DO I=1,NDIHE
         QXTRA = .TRUE.
         DO J=1,NPHI
            IF ( ((TJP(I).EQ.JP(J)).AND.(TKP(I).EQ.KP(J))).OR.
     @           ((TJP(I).EQ.KP(J)).AND.(TKP(I).EQ.JP(J))) ) THEN
               QXTRA = .FALSE.
            END IF
         END DO
         IF (QXTRA) THEN
         DO J=1,NIMPHI
            QIJ = .FALSE.
            QJK = .FALSE.
            QKL = .FALSE.
            DO K = 1,NBND(IM(J))
               IF (CONATM(IM(J),K).EQ.JM(J)) THEN
                  QIJ = .TRUE.
               END IF
            END DO
            DO K = 1,NBND(JM(J))
               IF (CONATM(JM(J),K).EQ.KM(J)) THEN
                  QJK = .TRUE.
               END IF
            END DO
            DO K = 1,NBND(KM(J))
               IF (CONATM(KM(J),K).EQ.LM(J)) THEN
                  QKL = .TRUE.
               END IF
            END DO
            Q1234 = QIJ.AND.QJK.AND.QKL
C
C fixed wrong reference to JP, KP, mrt 4/15/10
            IF (Q1234) THEN
               IF ( ((TJP(I).EQ.JM(J)).AND.(TKP(I).EQ.KM(J))).OR.
     @              ((TJP(I).EQ.KM(J)).AND.(TKP(I).EQ.JM(J))) ) THEN
                  QXTRA = .FALSE.
               END IF
            END IF
         END DO
         END IF
C
         IF (QXTRA) THEN
            IF (ZERO.GE.KDIHMAX) THEN
               QEXCL = .TRUE.
               DO IXCL=1,NDHXCL
                  JXCL=JDHXCL(IXCL)
                  KXCL=KDHXCL(IXCL)
                  IF (((JXCL.EQ.TJP(I)).AND.(KXCL.EQ.TKP(I))).OR.
     @                ((JXCL.EQ.TKP(I)).AND.(KXCL.EQ.TJP(I)))) THEN
                     QEXCL=.FALSE.
                  END IF
               END DO
C
               IF (QEXCL) THEN
               DO IINC=1,NDHINC
                  JINC=JDHINC(IINC)
                  KINC=KDHINC(IINC)
                  IF (((JINC.EQ.TJP(I)).AND.(KINC.EQ.TKP(I))).OR.
     @                ((JINC.EQ.TKP(I)).AND.(KINC.EQ.TJP(I)))) THEN
                     QEXCL=.FALSE.
                  END IF
               END DO
               END IF
C
               IF (QEXCL) THEN
                  NDHXCL = NDHXCL + 1
                  JDHXCL(NDHXCL) = TJP(I)
                  KDHXCL(NDHXCL) = TKP(I)
               END IF
            END IF
         END IF
      END DO
C
C go through the lists, and flag atoms to be excluded
      DO I=1,NPHI
         IF (ACPC(I).GE.KDIHMAX) THEN
            QEXCL=.TRUE.
            DO IXCL=1,NDHXCL
               JXCL=JDHXCL(IXCL)
               KXCL=KDHXCL(IXCL)
               IF (((JXCL.EQ.JP(I)).AND.(KXCL.EQ.KP(I))).OR.
     @             ((JXCL.EQ.KP(I)).AND.(KXCL.EQ.JP(I)))) THEN
                  QEXCL=.FALSE.
               END IF
            END DO
C
            IF (QEXCL) THEN
            DO IINC=1,NDHINC
               JINC=JDHINC(IINC)
               KINC=KDHINC(IINC)
               IF (((JINC.EQ.JP(I)).AND.(KINC.EQ.KP(I))).OR.
     @             ((JINC.EQ.KP(I)).AND.(KINC.EQ.JP(I)))) THEN
                  QEXCL=.FALSE.
               END IF
            END DO
            END IF
C
            IF (QEXCL) THEN
               NDHXCL = NDHXCL + 1
               JDHXCL(NDHXCL) = JP(I)
               KDHXCL(NDHXCL) = KP(I)
            END IF
         END IF
      END DO
C
      DO I=1,NIMPHI
         IF (ACIC(I).GE.KDIHMAX) THEN
            QIJ = .FALSE.
            QJK = .FALSE.
            QKL = .FALSE.
            DO J=1,NBND(IM(I))
               IF (CONATM(IM(I),J).EQ.JM(I)) THEN
                  QIJ = .TRUE.
               END IF
            END DO
            DO J=1,NBND(JM(I))
               IF (CONATM(JM(I),J).EQ.KM(I)) THEN
                  QJK = .TRUE.
               END IF
            END DO
            DO J=1,NBND(KM(I))
               IF (CONATM(KM(I),J).EQ.LM(I)) THEN
                  QKL = .TRUE.
               END IF
            END DO
            Q1234 = QIJ.AND.QJK.AND.QKL
C
            IF (Q1234) THEN
               QEXCL=.TRUE.
               DO IXCL=1,NDHXCL
                  JXCL=JDHXCL(IXCL)
                  KXCL=KDHXCL(IXCL)
                  IF (((JXCL.EQ.JM(I)).AND.(KXCL.EQ.KM(I))).OR.
     @                ((JXCL.EQ.KM(I)).AND.(KXCL.EQ.JM(I)))) THEN
                     QEXCL=.FALSE.
                  END IF
               END DO
C
               IF (QEXCL) THEN
               DO IINC=1,NDHINC
                  JINC=JDHINC(IINC)
                  KINC=KDHINC(IINC)
                  IF (((JINC.EQ.JM(I)).AND.(KINC.EQ.KM(I))).OR.
     @                ((JINC.EQ.KM(I)).AND.(KINC.EQ.JM(I)))) THEN
                     QEXCL=.FALSE.
                  END IF
               END DO
               END IF
C
               IF (QEXCL) THEN
                  NDHXCL = NDHXCL + 1
                  JDHXCL(NDHXCL) = JM(I)
                  KDHXCL(NDHXCL) = KM(I)
               END IF
            END IF
         END IF
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE CNTDIH(NDIHE,NDIH,CONDIH,NDHXCL,JDHXCL,KDHXCL,
     &                  NDHINC,JDHINC,KDHINC,NBNDXCL,IBXCL,JBXCL,
     &                  NBND,TJP,TKP,QDRY)
C
C get the number of dihedrals each atom is involved in
C    don't count the excluded ones
C also collect the atom numbers for the dihedrals
C
C      NDIH   stores the number of dihedrals each atom is involved in
C      CONDIH lists (for each atom) the atoms connected through dihedrals
C      QDRY: if set to TRUE - dry run - determine the maximum 
C            number of dihedrals per atom.  If set to FALSE, fill the
C            properly allocated array.
C
C Author: Luke Rice
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER NDIHE,NDIH(*),TJP(*),TKP(*)
      INTEGER NDHXCL,JDHXCL(*),KDHXCL(*)
      INTEGER NDHINC,JDHINC(*),KDHINC(*)
      INTEGER NBNDXCL,IBXCL(*),JBXCL(*)
      INTEGER CONDIH(NATOM,*),NBND(*)
      LOGICAL QDRY
C
C Local
      INTEGER I,J,k,JJ,KK,IX,JX,KX,JI,KI
      LOGICAL QACPT
C
      DO I=1,NATOM
         NDIH(I) = 0
         IF (.NOT.QDRY) THEN
           DO J=1,MAXDIHE
              CONDIH(I,J) = 0
           END DO
         END IF
      END DO
C
      DO I=1,NDIHE
         JJ = TJP(I)
         KK = TKP(I)
         QACPT = .TRUE.
C
         DO J=1,NDHXCL
            JX = JDHXCL(J)
            KX = KDHXCL(J)
            IF ( (JX.EQ.JJ).AND.(KX.EQ.KK).OR.
     &           (KX.EQ.JJ).AND.(JX.EQ.KK) ) THEN
               QACPT = .FALSE.
            END IF
         END DO
C
         DO J=1,NDHINC
            JI = JDHINC(J)
            KI = KDHINC(J)
            IF ( (JI.EQ.JJ).AND.(KI.EQ.KK).OR.
     &           (KI.EQ.JJ).AND.(JI.EQ.KK) ) THEN
               QACPT = .TRUE.
            END IF
         END DO
C
         DO J=1,NBNDXCL
            IX = IBXCL(J)
            JX = JBXCL(J)
            IF ( (IX.EQ.JJ).AND.(JX.EQ.KK).OR.
     &           (JX.EQ.JJ).AND.(IX.EQ.KK) ) THEN
               QACPT = .FALSE.
            END IF
         END DO
C
         IF (QACPT) THEN
            NDIH(JJ) = NDIH(JJ) + 1
            NDIH(KK) = NDIH(KK) + 1
C
C begin modification ATB 5/28/08
            IF (QDRY) THEN
               MAXDIHE=MAX(MAXDIHE,NDIH(JJ),NDIH(KK))
            ELSE
               CONDIH(JJ,NDIH(JJ)) = KK
               CONDIH(KK,NDIH(KK)) = JJ
            END IF
C
         END IF
C
      END DO
C
C
      IF (.NOT.QDRY) THEN
C
C Some final cleaning up
C (nedded for some of the customization stuff)
        DO I=1,NATOM
           IF ( (NDIH(I).EQ.1).AND.(NBND(I).EQ.1) ) THEN
              NDIH(I) = 0
              CONDIH(I,1) = 0
           END IF
        END DO
C
        DO I=1,NATOM
           DO J=1,NDIH(I)
              JJ=CONDIH(I,J)
              IF (NBND(JJ).EQ.1) THEN
                 NDIH(I) = NDIH(I) - 1
                 DO K=J,NDIH(I)
                    CONDIH(I,K) = CONDIH(I,K+1)
                 END DO
              END IF
           END DO
        END DO
      END IF
C end modifcation
C
      RETURN
      END
C===============================================================
      SUBROUTINE GETPNTR(LIST,POINTR,CONATM,CONDIH,
     @                   NBND,NDIH,GPLIST,GROUP)
C
C Uses the list of dihedrals to get atom groupings
C   The first step is to fill GPLIST, an array which keeps
C   track of the group dependencies (through one bond at a time).
C   Then go through GPLIST to get the full dependencies.
C
C Author: Luke Rice
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER LIST(*),POINTR(*),CONATM(NATOM,*)
      INTEGER GROUP(*),GPLIST(NATOM,*),CONDIH(NATOM,*)
      INTEGER NBND(*),NDIH(*)
C
C Local
      INTEGER I,J,K,ILOOK,ITAKEN,NFIXED
      INTEGER BNDCON,DIHCON,IGROUP,II,TEMPGP
      LOGICAL QLOOK,QINCR
C
C   Go through the connectivity and get dependencies (through one bond)
      DO I = 1, NATOM
         GROUP(I) = 0
         DO J = 1, MAXBND
            GPLIST(I,J) = 0
         END DO
      END DO
C
      DO I = 1, NATOM
C
         IF (IMOVE(I).EQ.0) THEN
            ILOOK = 0
            DO J = 1, NBND(I)
               QLOOK = .TRUE.
               BNDCON = CONATM(I,J)
               DO K = 1, NDIH(I)
                  DIHCON = CONDIH(I,K)
                  IF (BNDCON.EQ.DIHCON) THEN
                     QLOOK = .FALSE.
                  END IF
               END DO
               IF (QLOOK) THEN
                  ILOOK = ILOOK+1
                  GPLIST(I,ILOOK) = BNDCON
               END IF
            END DO
         END IF
      END DO
C
C Figure out groupings
      IGROUP = 0
      DO I=1,NATOM
C
         IF (IMOVE(I).EQ.0) THEN
            ILOOK = 1
            IF (GROUP(I).EQ.0) THEN
               IGROUP = IGROUP + 1
               GROUP(I) = IGROUP
            END IF
            DO WHILE 
     &      ((GPLIST(I,MIN(MAXBND,ILOOK)).NE.0).AND.(ILOOK.LE.MAXBND))
               II = GPLIST(I,ILOOK)
               IF (II.GT.I) THEN
                  IF (GROUP(II).EQ.0) THEN
                     GROUP(II) = GROUP(I)
                  ELSE
                     TEMPGP = GROUP(I)
                     DO J=1,NATOM
                        IF (GROUP(J).EQ.TEMPGP) GROUP(J)=GROUP(II)
                     END DO
                  END IF
               END IF
               ILOOK = ILOOK + 1
            END DO
         END IF
      END DO
C
C Now fill up the pointr and list arrays
      NFIXED = 0
      DO I=1,NATOM
         IF (IMOVE(I).EQ.1) NFIXED = NFIXED + 1
      END DO
C
      IGROUP = 0
      NGP = 0
      ITAKEN = 0
      POINTR(1) = 0
      POINTR(2) = 0
      DO WHILE (ITAKEN.LT.(NATOM-NFIXED))
         II = POINTR(NGP+1)
         QINCR = .FALSE.
         DO I=1,NATOM
C
            IF (IMOVE(I).EQ.0) THEN
               IF (GROUP(I).EQ.(IGROUP+1)) THEN
                  QINCR = .TRUE.
                  ITAKEN = ITAKEN + 1
                  II = II + 1
                  LIST(II) = I
               END IF
            END IF
         END DO
         IGROUP = IGROUP + 1
         IF (QINCR) THEN
            NGP = NGP + 1
            POINTR(NGP+1) = II
         END IF
      END DO
C
      RETURN
      END
C===============================================================
      SUBROUTINE TORINI
C
C Routine initializes HEAP pointers for
C   torsion angle dynamics
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'dtorsion.inc'
C
C set the QRESET flag to true: this will trigger generation of the torsion topology
      QRESET=.TRUE.
C
      NDEGF=0
C
      HPLIST = 0
      HPPOINTR = 0
C
      HPNJOINT = 0
      HPGRPSEQ = 0
      HPNIN = 0
      HPNOUT = 0
      HPSIZE = 0
      HPCHNLEN = 0
      HPNCHAINS = 0
      HPXCM = 0
      HPYCM = 0
      HPZCM = 0
      HPCOEFFI = 0
      HPCOEFFO = 0
      HPQTIP = 0
      HPMM = 0
      HPYY = 0
      HPQ0 = 0
      HPQ1 = 0
      HPQ2 = 0
      HPQ3 = 0
      HPXB = 0
      HPYB = 0
      HPZB = 0
      HPQDOT = 0
      HPHHX = 0
      HPHHY = 0
      HPHHZ = 0
      HPSSX = 0
      HPSSY = 0
      HPSSZ = 0
      HPFRICT = 0
      HPQQ = 0
      HPB1 = 0
      HPB2 = 0
      HPDD = 0
      HPKK = 0
      HPLL = 0
      HPYYDOT = 0
      HPQ0O = 0
      HPQ1O = 0
      HPQ2O = 0
      HPQ3O = 0
      HPYYO = 0
      HPRKQ = 0
      HPRKQDOT = 0
      HPQDDOT = 0
      HPQDOTO = 0
      HPQO = 0
      HPQ = 0
C
      RETURN
      END
C===============================================================
      SUBROUTINE TORHP(NATOM,QLIST)
C
C Routine (re-)allocates or deallocates various heap arrays 
C for torsion-angle dynamics
C
C QLIST: flag indicating if the LIST, POINTR arrays should
C deallocated. The arrays are allocated in subroutine TORMD.
C 
C Author: Luke Rice
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'dtorsion.inc'
      INTEGER NATOM
      LOGICAL QLIST
C
C see if space already allocated
      IF (QLIST) THEN
         IF (HPLIST.NE.0) CALL FREHP(HPLIST,INTEG4(ONATOM))
         IF (HPPOINTR.NE.0) CALL FREHP(HPPOINTR,INTEG4(ONATOM+2))
      END IF
      IF (HPNJOINT.NE.0) CALL FREHP(HPNJOINT,INTEG4(ONGP))
      IF (HPGRPSEQ.NE.0)
     &  CALL FREHP(HPGRPSEQ,INTEG4(OMAXTREE*OMAXCHN*OMAXLEN))
      IF (HPNIN.NE.0) 
     & CALL FREHP(HPNIN,INTEG4(OMAXTREE*OMAXCHN*(OMAXLEN+1)))
      IF (HPNOUT.NE.0) 
     & CALL FREHP(HPNOUT,INTEG4(OMAXTREE*OMAXCHN*OMAXLEN))
      IF (HPSIZE.NE.0) CALL FREHP(HPSIZE,INTEG4(ONGP))
      IF (HPCHNLEN.NE.0) 
     & CALL FREHP(HPCHNLEN,INTEG4(OMAXTREE*OMAXCHN))
      IF (HPNCHAINS.NE.0) CALL FREHP(HPNCHAINS,INTEG4(OMAXTREE))
      IF (HPXCM.NE.0) CALL FREHP(HPXCM,IREAL8(ONGP))
      IF (HPYCM.NE.0) CALL FREHP(HPYCM,IREAL8(ONGP))
      IF (HPZCM.NE.0) CALL FREHP(HPZCM,IREAL8(ONGP))
      IF (HPCOEFFI.NE.0)
     &  CALL FREHP(HPCOEFFI,IREAL8(OMAXTREE*OMAXCHN*OMAXLEN))
      IF (HPCOEFFO.NE.0)
     &  CALL FREHP(HPCOEFFO,IREAL8(OMAXTREE*OMAXCHN*OMAXLEN))
      IF (HPQTIP.NE.0) CALL FREHP(HPQTIP,INTEG4(ONGP))
      IF (HPMM.NE.0) CALL FREHP(HPMM,IREAL8(36*ONGP))
      IF (HPYY.NE.0) CALL FREHP(HPYY,IREAL8(6*ONGP))
      IF (HPQ0.NE.0) CALL FREHP(HPQ0,IREAL8(ONGP))
      IF (HPQ1.NE.0) CALL FREHP(HPQ1,IREAL8(ONGP))
      IF (HPQ2.NE.0) CALL FREHP(HPQ2,IREAL8(ONGP))
      IF (HPQ3.NE.0) CALL FREHP(HPQ3,IREAL8(ONGP))
      IF (HPXB.NE.0) CALL FREHP(HPXB,IREAL8(ONATOM))
      IF (HPYB.NE.0) CALL FREHP(HPYB,IREAL8(ONATOM))
      IF (HPZB.NE.0) CALL FREHP(HPZB,IREAL8(ONATOM))
      IF (HPQDOT.NE.0) CALL FREHP(HPQDOT,IREAL8(ONGP))
      IF (HPHHX.NE.0) CALL FREHP(HPHHX,IREAL8(ONGP))
      IF (HPHHY.NE.0) CALL FREHP(HPHHY,IREAL8(ONGP))
      IF (HPHHZ.NE.0) CALL FREHP(HPHHZ,IREAL8(ONGP))
      IF (HPSSX.NE.0) CALL FREHP(HPSSX,IREAL8(ONGP))
      IF (HPSSY.NE.0) CALL FREHP(HPSSY,IREAL8(ONGP))
      IF (HPSSZ.NE.0) CALL FREHP(HPSSZ,IREAL8(ONGP))
      IF (HPFRICT.NE.0) CALL FREHP(HPFRICT,IREAL8(ONATOM))
      IF (HPQQ.NE.0) CALL FREHP(HPQQ,IREAL8(6*ONGP))
      IF (HPB1.NE.0) CALL FREHP(HPB1,IREAL8(36*ONGP))
      IF (HPB2.NE.0) CALL FREHP(HPB2,IREAL8(6*ONGP))
      IF (HPDD.NE.0) CALL FREHP(HPDD,IREAL8(6*ONGP))
      IF (HPKK.NE.0) CALL FREHP(HPKK,IREAL8(36*ONGP))
      IF (HPLL.NE.0) CALL FREHP(HPLL,IREAL8(6*ONGP))
      IF (HPYYDOT.NE.0) CALL FREHP(HPYYDOT,IREAL8(6*ONGP))
      IF (HPQ0O.NE.0) CALL FREHP(HPQ0O,IREAL8(ONGP))
      IF (HPQ1O.NE.0) CALL FREHP(HPQ1O,IREAL8(ONGP))
      IF (HPQ2O.NE.0) CALL FREHP(HPQ2O,IREAL8(ONGP))
      IF (HPQ3O.NE.0) CALL FREHP(HPQ3O,IREAL8(ONGP))
      IF (HPYYO.NE.0) CALL FREHP(HPYYO,IREAL8(6*ONGP))
      IF (HPRKQ.NE.0) CALL FREHP(HPRKQ,IREAL8(ONGP))
      IF (HPRKQDOT.NE.0) CALL FREHP(HPRKQDOT,IREAL8(ONGP))
      IF (HPQDDOT.NE.0) CALL FREHP(HPQDDOT,IREAL8(ONGP))
      IF (HPQDOTO.NE.0) CALL FREHP(HPQDOTO,IREAL8(ONGP))
      IF (HPQO.NE.0) CALL FREHP(HPQO,IREAL8(ONGP))
      IF (HPQ.NE.0) CALL FREHP(HPQ,IREAL8(ONGP))
C
C Store values of NATOM, NGROUP, MAXTREE, MAXCHN, MAXLEN
      ONATOM = NATOM
      ONGP = NGP
      OMAXTREE = MAXTREE
      OMAXCHN = MAXCHN
      OMAXLEN = MAXLEN
C
C Now allocate new space
      IF (NATOM.GT.0) THEN
         HPNJOINT = ALLHP(INTEG4(NGP))
         HPGRPSEQ = ALLHP(INTEG4(MAXTREE*MAXCHN*MAXLEN))
         HPNIN = ALLHP(INTEG4(MAXTREE*MAXCHN*(MAXLEN+1)))
         HPNOUT = ALLHP(INTEG4(MAXTREE*MAXCHN*MAXLEN))
         HPSIZE = ALLHP(INTEG4(NGP))
         HPCHNLEN = ALLHP(INTEG4(MAXTREE*MAXCHN))
         HPNCHAINS = ALLHP(INTEG4(MAXTREE))
         HPXCM = ALLHP(IREAL8(NGP))
         HPYCM = ALLHP(IREAL8(NGP))
         HPZCM = ALLHP(IREAL8(NGP))
         HPCOEFFI = ALLHP(IREAL8(MAXTREE*MAXCHN*MAXLEN))
         HPCOEFFO = ALLHP(IREAL8(MAXTREE*MAXCHN*MAXLEN))
         HPQTIP = ALLHP(INTEG4(NGP))
         HPMM = ALLHP(IREAL8(36*NGP))
         HPYY = ALLHP(IREAL8(6*NGP))
         HPQ0 = ALLHP(IREAL8(NGP))
         HPQ1 = ALLHP(IREAL8(NGP))
         HPQ2 = ALLHP(IREAL8(NGP))
         HPQ3 = ALLHP(IREAL8(NGP))
         HPXB = ALLHP(IREAL8(NATOM))
         HPYB = ALLHP(IREAL8(NATOM))
         HPZB = ALLHP(IREAL8(NATOM))
         HPQDOT = ALLHP(IREAL8(NGP))
         HPHHX = ALLHP(IREAL8(NGP))
         HPHHY = ALLHP(IREAL8(NGP))
         HPHHZ = ALLHP(IREAL8(NGP))
         HPSSX = ALLHP(IREAL8(NGP))
         HPSSY = ALLHP(IREAL8(NGP))
         HPSSZ = ALLHP(IREAL8(NGP))
         HPFRICT=ALLHP(IREAL8(NATOM))
         HPQQ=ALLHP(IREAL8(6*NGP))
         HPB1=ALLHP(IREAL8(36*NGP))
         HPB2=ALLHP(IREAL8(6*NGP))
         HPDD=ALLHP(IREAL8(6*NGP))
         HPKK=ALLHP(IREAL8(36*NGP))
         HPLL=ALLHP(IREAL8(6*NGP))
         HPYYDOT=ALLHP(IREAL8(6*NGP))
         HPQ0O=ALLHP(IREAL8(NGP))
         HPQ1O=ALLHP(IREAL8(NGP))
         HPQ2O=ALLHP(IREAL8(NGP))
         HPQ3O=ALLHP(IREAL8(NGP))
         HPYYO=ALLHP(IREAL8(6*NGP))
         HPRKQ=ALLHP(IREAL8(NGP))
         HPRKQDOT=ALLHP(IREAL8(NGP))
         HPQDDOT=ALLHP(IREAL8(NGP))
         HPQDOTO=ALLHP(IREAL8(NGP))
         HPQO=ALLHP(IREAL8(NGP))
         HPQ=ALLHP(IREAL8(NGP))
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE TOPOLOGY(LIST,POINTR,NJOINT,GRPSEQ,
     @     NIN,NOUT,SIZE,CHNLEN,TAKEN,TREE,NCHAINS,JATOMS,
     @     JBONDS,CONGRP,IBXCL,JBXCL,NBNDXCL,XCLBND)
C
C Routine sets up the topology for torsion angle dynamics.
C
C     GRPSEQ(I,J,K): stores the ttopology structure as follows:
C            I       is the tree label
C              J     is the chain label for tree I
C                K   is the node number along chain J
C     NCHAINS        counts the number of chains in a given tree
C     CHNLEN         counts the length of the chains in a given tree
C     NJOINT         counts the number of connections for each group
C     NIN, NOUT      lists the atoms (for each chain node) involved in joints
C     JATOMS         lists the atoms (for each group) involved in joints
C     JBONDS         stores the bond numbers (for each group) of the JATOMS
C     CONGRP         lists (for each group) the connected groups
C     SIZE           is equal to (natms - 1) for each group
C     TAKEN          keeps track of what groups have been put in a tree
C     TREE           is an identifier to say which tree a group belongs in
C
C Author: Luke Rice
C
      IMPLICIT NONE
C
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
C I/O
      INTEGER LIST(*),POINTR(*),XCLBND(*)
      INTEGER NIN(MAXTREE,MAXCHN,*),NOUT(MAXTREE,MAXCHN,*)
      INTEGER NJOINT(*),SIZE(*),K,NCHAIN,CHNLEN(MAXTREE,*)
      INTEGER GRPSEQ(MAXTREE,MAXCHN,*),NCHAINS(*),TREE(*)
      INTEGER TAKEN(*),JATOMS(NGP,*),JBONDS(NGP,*)
      INTEGER THISONE,LASTONE,CONGRP(NGP,*),ONJOINT(MAXA)
      INTEGER IBXCL(*),JBXCL(*),NBNDXCL
C
C Local
      INTEGER I,POINTER,POINTER2,IBOND,JOINTS,SKIPJOINT,INBOND
      INTEGER II,BASE,J,ICHAIN,IGP,CHAINS,IJOINT,III
      INTEGER ITEST,TOTCHN,INATM,OUTATM
      INTEGER TREETAG,ITREE,NTAKEN
      INTEGER MXCHAINS, MXCHNLEN
      DOUBLE PRECISION ZERO
      LOGICAL QTOPO2
      LOGICAL COND,COND2,QLOOP,QCLOSED,qincr,qskbnd
      PARAMETER (ZERO=0.0D0)
C
C begin
      QTOPO2=.FALSE.
      QCLOSED = .FALSE.
C
      MXCHAINS=0
      MXCHNLEN=0
C
      IF (QTOPO2) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A,A/)') '           - Regenerating molecular topology',
     @               ' ignoring closed loops -'
      WRITE(6,'(A,A)')
     @     ' ---------------------------------------',
     @     '----------------------------------------'
      WRITE(6,'(A)') ' '
      END IF
C
C
C go through the list of bonds, and accumulate interesting ones
C    (interesting means bonds between groups)
      DO I=1,NGP
      NJOINT(I)=ZERO
      ONJOINT(I)=ZERO
      TAKEN(I) = 0
      END DO
C
      DO IBOND=1,NBOND
      IF (XCLBND(IBOND).EQ.0) THEN
      DO I=1,NGP
      DO POINTER=POINTR(I)+1,POINTR(I+1)
      IF (IB(IBOND).EQ.LIST(POINTER)) THEN
      COND=.TRUE.
      DO POINTER2=POINTR(I)+1,POINTR(I+1)
      IF (JB(IBOND).EQ.LIST(POINTER2)) THEN
      COND=.FALSE.
      END IF
      END DO
      IF (COND) THEN
      NJOINT(I)=NJOINT(I)+1
C
      IF (NJOINT(I).GT.MAXJNT) THEN
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'MAXJNT exceeded. Fatal coding error')
      END IF
C
      JATOMS(I,NJOINT(I))=IB(IBOND)
      JBONDS(I,NJOINT(I))=IBOND
      END IF
      END IF
      IF (JB(IBOND).EQ.LIST(POINTER)) THEN
      COND=.TRUE.
      DO POINTER2=POINTR(I)+1,POINTR(I+1)
      IF (IB(IBOND).EQ.LIST(POINTER2)) THEN
      COND=.FALSE.
      END IF
      END DO
      IF (COND) THEN
      NJOINT(I)=NJOINT(I)+1
C
      IF (NJOINT(I).GT.MAXJNT) THEN
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'MAXJNT exceeded. Fatal coding error')
      END IF
C
      JATOMS(I,NJOINT(I))=JB(IBOND)
      JBONDS(I,NJOINT(I))=IBOND
      END IF
      END IF
      END DO
      ONJOINT(I) = NJOINT(I)
      END DO
      END IF
      END DO
C
C find out what group(s) this one is connected to
C store the list of connected groups in CONGRP
      DO I=1,NGP
      DO J=1,NJOINT(I)
      DO II=1,NGP
      IF (II.NE.I) THEN
      DO POINTER=POINTR(II)+1,POINTR(II+1)
      IF (JB(JBONDS(I,J)).EQ.LIST(POINTER)) THEN
      CONGRP(I,J)=II
      ELSE IF (IB(JBONDS(I,J)).EQ.LIST(POINTER)) THEN
      CONGRP(I,J)=II
      END IF
      END DO
      END IF
      END DO
      END DO
      END DO
C
C Put a tree label on each group
C If two bodies are connected, they belong to the same tree
      NTREE = 0
      DO I=1,NGP
      TREE(I) = 0
      END DO
C
      DO I=1,NGP
      IF (TREE(I).EQ.0) THEN
      NTREE = NTREE + 1
      TREE(I) = NTREE
      END IF
      DO J=1,NJOINT(I)
      IF (TREE(CONGRP(I,J)).EQ.0) THEN
      TREE(CONGRP(I,J)) = TREE(I)
      ELSE
      IF (TREE(CONGRP(I,J)).LT.TREE(I)) THEN
      TREETAG = TREE(CONGRP(I,J))
      ELSE
      TREETAG = TREE(I)
      END IF
      DO K=1,NGP
      IF ( (TREE(K).EQ.TREE(I)).OR.
     @     (TREE(K).EQ.TREE(CONGRP(I,J))) ) THEN
      TREE(K) = TREETAG
      END IF
      END DO
      END IF
      END DO
      END DO
C
C Loop through all groups and build tree(s)
      NTAKEN = 0
      NTREE = 0
      DO WHILE (NTAKEN.LT.NGP)
      NTREE = NTREE + 1
C
      IF (NTREE.GT.MAXTREE) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A)') 
     & ' Please send a bug report to brunger@stanford.edu'
      WRITE(6,'(2A)') 
     & '  Please increase the "maximum number of distinct bodies"',
     & ' (torsion_maxtree) and re-run.'
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
C
C various checks for the base body
      BASE=0
      DO I=1,NGP
C
C first, check if we have a stand-alone body (i.e. zero joints)
      IF ((NJOINT(I).EQ.0).AND.(TAKEN(I).EQ.0)) THEN
C
C stand-alone bodies with at least two atoms are allowed, ATB 3/28/10
      IF ((POINTR(I+1)-POINTR(I)-1).GT.0) THEN
      BASE=I
C
      ELSE
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A)') '   ERROR: A single-group tree does not have an'
      WRITE(6,'(A)') '            acceptable base.'
      WRITE(6,'(A)') '   There are certain limitations with the current'
      WRITE(6,'(A/)')'   implementation of torsion angle dynamics.'
C
      WRITE(6,'(A)') '   The problem is caused by isolated'
      WRITE(6,'(A/)')'   (not covalently bonded) atoms.'
C
      WRITE(6,'(A)') '   We have made every effort to capture these'
      WRITE(6,'(A)') '   situations automatically.  However, in some'
      WRITE(6,'(A/)')'   cases you may still encounter this problem.'
C
      WRITE(6,'(A)') '   In order to fix this problem'
      WRITE(6,'(A)') '   (a) fix (fix atom selection) the offending'
      WRITE(6,'(A)') '       atoms'
      WRITE(6,'(A)') '   (b) if everything fails, use the Cartesian'
      WRITE(6,'(A/)')'       dynamics option.'
C
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
      END IF
      END DO
C
C second, for the groups with joints, this finds a base group with at least 3 
C atoms and with at least one joint
      DO I=1,NGP
      IF ((NJOINT(I).GE.1).AND.(TAKEN(I).EQ.0)) THEN
      IF ((POINTR(I+1)-POINTR(I)-1).GT.1) THEN
      BASE=I
      END IF
      END IF
      END DO
C
C we also allow a base group with two atoms and with at least one joint, ATB 3/26/10
      IF (BASE.EQ.0) THEN
      DO I=1,NGP
      IF ((NJOINT(I).GE.1).AND.(TAKEN(I).EQ.0)) THEN
      IF ((POINTR(I+1)-POINTR(I)-1).GT.0) THEN
      BASE=I
      END IF
      END IF
      END DO
      END IF
C
      IF (BASE.EQ.0) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A)') '   ERROR: There are no suitable base groups.     '
      WRITE(6,'(A)') '   This problem can be caused by isolated '
      WRITE(6,'(A)') '   bonding networks with undefined or weak '
      WRITE(6,'(A)') '   dihedral force constants. '
      WRITE(6,'(A)') '   The atoms that cannot be placed in a tree '
      WRITE(6,'(A)') '   are listed below: '
      DO I=1,NGP
      IF (TAKEN(I).EQ.0) THEN
      DO POINTER=POINTR(I)+1,POINTR(I+1)
      WRITE(6,'(8A)')
     @     '    %atoms "', SEGID(LIST(POINTER)),'"-',
     @                     RESID(LIST(POINTER)),'-',
     @                     RES(LIST(POINTER)),'-',
     @                     TYPE(LIST(POINTER))
      END DO
      END IF
      END DO
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
C
C now get the big one with the most joints
      DO I=1,NGP
      IF (TAKEN(I).EQ.0) THEN
      IF ((NJOINT(I).GE.NJOINT(BASE)).AND.
     @     ((POINTR(I+1)-POINTR(I)-1).GE.
     @     (POINTR(BASE+1)-POINTR(BASE)-1)) ) THEN
      BASE=I
      END IF
      END IF
      END DO
C
C Count the total expected number of chains in this tree
C (to be corrected later for closed loops)
      TOTCHN=0
      TREETAG = TREE(BASE)
      DO I=1,NGP
      IF (TREE(I).EQ.TREETAG) THEN
      IF (I.EQ.BASE) THEN
      TOTCHN=TOTCHN+NJOINT(I)
      ELSE
      IF (NJOINT(I).GE.3) THEN
      TOTCHN=TOTCHN+NJOINT(I)-1
      END IF
      END IF
      END IF
      END DO
C
C Define processing sequence
C   Treat the base first:
C   Every joint of the base contributes a new chain
C   (unless the connected group is already taken).
C   Each chain is built by adding consecutive bodies
C   until the most recently added group is taken or doesn't
C   have two joints (i.e. is a tip).
C   A closed loop is treated by erasing the information about
C   the closing connection.
      NCHAINS(NTREE) = 0
C
      IF (NJOINT(BASE).EQ.0) THEN
      NCHAIN=1
      NCHAINS(NTREE)=NCHAIN
      CHNLEN(NTREE,NCHAIN)=1
      GRPSEQ(NTREE,1,1)=BASE
      END IF
C
      TAKEN(BASE)=1
      NTAKEN = NTAKEN + 1
      ICHAIN = 1
C
C Loop over all joints of the base
      DO WHILE (ICHAIN.LE.NJOINT(BASE))
      GRPSEQ(NTREE,ICHAIN,1)=BASE
      IF (TAKEN(CONGRP(BASE,ICHAIN)).EQ.0) THEN
      NCHAINS(NTREE) = NCHAINS(NTREE) + 1
      NCHAIN = NCHAINS(NTREE)
      GRPSEQ(NTREE,ICHAIN,2)=CONGRP(BASE,ICHAIN)
      TAKEN(GRPSEQ(NTREE,ICHAIN,2)) = 1
      NTAKEN = NTAKEN + 1
      CHNLEN(NTREE,ICHAIN)=2
C
      QLOOP = .FALSE.
C
C Continue adding to the chain until NJOINT.NE.2 or a loop is closed
      DO WHILE
     @     ( (NJOINT(GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN))).EQ.2)
     @     .AND.(.NOT.QLOOP) )
      THISONE = GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN))
      LASTONE = GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN)-1)
      DO I=1,2
      IF (CONGRP(THISONE,I).EQ.LASTONE) THEN
      IF (TAKEN(CONGRP(THISONE,3-I)).EQ.0) THEN
      NTAKEN = NTAKEN + 1
      CHNLEN(NTREE,ICHAIN)=CHNLEN(NTREE,ICHAIN)+1
C
      MXCHNLEN=MAX(MXCHNLEN,CHNLEN(NTREE,ICHAIN))
C
      IF (CHNLEN(NTREE,ICHAIN).GT.MAXLEN) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A)') 
     & ' Please send a bug report to brunger@stanford.edu'
      WRITE(6,'(2A)') 
     & ' Please increase the "maximum unbranched chain length"',
     & ' (torsion_maxlength) and re-run.'
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
      GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN)) = CONGRP(THISONE,3-I)
      TAKEN(GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN))) = 1
      ELSE
C
C A closed loop.  Fix up.
      QLOOP = .TRUE.
      QTOPO2 = .TRUE.
C
      IF (.NOT.QCLOSED) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A,A)') ' Closed bonding networks have been detected.',
     @               ' The bonds between the following'
      WRITE(6,'(A/)') '   atoms will not be constrained:'
      QCLOSED = .TRUE.
      END IF
C
      DO J=1,NJOINT(THISONE)
      DO K=1,NJOINT(LASTONE)
      COND = ( (IB(JBONDS(LASTONE,K)).EQ.JATOMS(THISONE,J)).OR.
     @     (JB(JBONDS(LASTONE,K)).EQ.JATOMS(THISONE,J)) )
      IF (COND) THEN
      INATM = JATOMS(THISONE,J)
      INBOND = JBONDS(LASTONE,K)
      END IF
      END DO
      END DO
C
      DO J=1,NJOINT(THISONE)
      IF (JBONDS(THISONE,J).NE.INBOND) OUTATM = JATOMS(THISONE,J)
      END DO
C
      qskbnd = .false.
      DO J=1,NJOINT(CONGRP(THISONE,3-I))
      if (.not.qskbnd) then
      IF (IB(JBONDS(CONGRP(THISONE,3-I),J)).EQ.OUTATM) THEN
      qskbnd = .true.
C
C add to bond exclusion list
      IBOND=JBONDS(CONGRP(THISONE,3-I),J)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=IB(IBOND)
      JBXCL(NBNDXCL)=JB(IBOND)
C
      WRITE(6,'(17A)')
     @     '    %atoms "', SEGID(OUTATM),'-',RESID(OUTATM),'-',
     @     RES(OUTATM),'-',TYPE(OUTATM),'" and "',
     @     SEGID(JB(IBOND)),'-',
     @     RESID(JB(IBOND)),'-',
     @     RES(JB(IBOND)),'-',
     @     TYPE(JB(IBOND)),'"'
      END IF
C
C
      SKIPJOINT = J
C
      ELSE IF  (JB(JBONDS(CONGRP(THISONE,3-I),J)).EQ.OUTATM) THEN
      qskbnd = .true.
C
C add to bond exclusion list
      IBOND=JBONDS(CONGRP(THISONE,3-I),J)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=OUTATM
      JBXCL(NBNDXCL)=IB(IBOND)
C
      WRITE(6,'(17A)')
     @     '    %atoms "', SEGID(OUTATM),'-',RESID(OUTATM),'-',
     @     RES(OUTATM),'-',TYPE(OUTATM),'" and "',
     @     SEGID(IB(IBOND)),'-',
     @     RESID(IB(IBOND)),'-',
     @     RES(IB(IBOND)),'-',
     @     TYPE(IB(IBOND)),'"'
      END IF
C
C
      SKIPJOINT = J
      END IF
      end if
      END DO
      qskbnd=.false.
      JOINTS = 0
      DO J=1,NJOINT(CONGRP(THISONE,3-I))
      IF (J.NE.SKIPJOINT) THEN
      JOINTS = JOINTS+1
      JATOMS(CONGRP(THISONE,3-I),JOINTS) = JATOMS(CONGRP(THISONE,3-I),J)
      JBONDS(CONGRP(THISONE,3-I),JOINTS) = JBONDS(CONGRP(THISONE,3-I),J)
      CONGRP(CONGRP(THISONE,3-I),JOINTS) = CONGRP(CONGRP(THISONE,3-I),J)
      END IF
      END DO
      NJOINT(CONGRP(THISONE,3-I)) = JOINTS
      IF (ONJOINT(CONGRP(THISONE,3-I)).GE.2) TOTCHN = TOTCHN - 1
C
      NJOINT(THISONE) = 1
      JATOMS(THISONE,1) = INATM
      JBONDS(THISONE,1) = INBOND
      CONGRP(THISONE,1) = LASTONE
      congrp(thisone,2) = 0
      END IF
      END IF
      END DO
      END DO
      ICHAIN = ICHAIN + 1
      ELSE
C
C A closed loop.  Fix up.
      IF (.NOT.QCLOSED) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A,A)') ' Closed bonding networks have been detected.',
     @               ' The bonds between the following'
      WRITE(6,'(A/)') '   atoms will not be constrained:'
      QCLOSED = .TRUE.
      END IF
C
C
      QTOPO2 = .TRUE.
      DO J=1,NJOINT(CONGRP(BASE,ICHAIN))
      IF (JBONDS(BASE,ICHAIN).EQ.JBONDS(CONGRP(BASE,ICHAIN),J)) THEN
C
C add to bond exclusion list
      IBOND=JBONDS(BASE,ICHAIN)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      INATM = IB(IBOND)
      OUTATM = JB(IBOND)
      IBXCL(NBNDXCL)=OUTATM
      JBXCL(NBNDXCL)=INATM
C
      WRITE(6,'(17A)')
     @     '    %atoms "', SEGID(OUTATM),'-',RESID(OUTATM),'-',
     @     RES(OUTATM),'-',TYPE(OUTATM),'" and "', SEGID(INATM),'-',
     @     RESID(INATM),'-', RES(INATM),'-', TYPE(INATM),'"'
      END IF
C
C
      SKIPJOINT = J
      END IF
      END DO
C
C only correct one joint at a time
      JOINTS = 0
      DO J=1,NJOINT(CONGRP(BASE,ICHAIN))
      IF (SKIPJOINT.NE.J) THEN
      JOINTS = JOINTS+1
      JATOMS(CONGRP(BASE,ICHAIN),JOINTS) = JATOMS(CONGRP(BASE,ICHAIN),J)
      JBONDS(CONGRP(BASE,ICHAIN),JOINTS) = JBONDS(CONGRP(BASE,ICHAIN),J)
      CONGRP(CONGRP(BASE,ICHAIN),JOINTS) = CONGRP(CONGRP(BASE,ICHAIN),J)
      END IF
      END DO
      NJOINT(CONGRP(BASE,ICHAIN)) = JOINTS
      IF (ONJOINT(CONGRP(BASE,ICHAIN)).GE.3) TOTCHN = TOTCHN - 1
      JOINTS = 0
      DO J=1,NJOINT(BASE)
      IF (J.NE.ICHAIN) THEN
      JOINTS = JOINTS+1
      JATOMS(BASE,JOINTS) = JATOMS(BASE,J)
      JBONDS(BASE,JOINTS) = JBONDS(BASE,J)
      CONGRP(BASE,JOINTS) = CONGRP(BASE,J)
      END IF
      END DO
      NJOINT(BASE) = JOINTS
      TOTCHN = TOTCHN - 1
      END IF
      END DO
C
C Now proceed outward from the already defined chain tips:
C   New chains are added to existing chain tips subject to the
C   following conditions:
C     + The tip has at least two joints
C     + The tip is not already listed as a base
C     + A chain is not built from the joint connecting "inwards"
C     + The connected group is not already taken
      DO WHILE (NCHAIN.LT.TOTCHN)
C
      CHAINS=NCHAIN
      DO ICHAIN=1,CHAINS
C
C Scan all existing chain tips.
      III=GRPSEQ(NTREE,ICHAIN,CHNLEN(NTREE,ICHAIN))
      IF (NJOINT(III).GE.2) THEN
      COND=.TRUE.
C
C If it's a base of a chain, skip it.
      DO ITEST=1,CHAINS
      IF (III.EQ.GRPSEQ(NTREE,ITEST,1)) THEN
      COND=.FALSE.
      END IF
      END DO
C
      IF (COND) THEN
      IJOINT = 1
C
C Loop over all joints.  Skip the inward connection.
      DO WHILE (IJOINT.LE.NJOINT(III))
      COND2=.TRUE.
      DO J=1,CHAINS
      IF ( (CONGRP(III,IJOINT).EQ.GRPSEQ(NTREE,J,CHNLEN(NTREE,J)-1))
     @     .AND.(III.EQ.GRPSEQ(NTREE,J,CHNLEN(NTREE,J))) ) THEN
      COND2=.FALSE.
      END IF
      END DO
C
      IF (COND2) THEN
      IF (TAKEN(CONGRP(III,IJOINT)).EQ.1) THEN
C
C A closed loop.  Fix up.
      QTOPO2=.TRUE.
      igp = congrp(iii,ijoint)
      IF (.NOT.QCLOSED) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A,A)') ' Closed bonding networks have been detected.',
     @               ' The bonds between the following'
      WRITE(6,'(A/)') '   atoms will not be constrained:'
      QCLOSED = .TRUE.
      END IF
C
C Again, only correct one joint at a time
      J=1
      DO WHILE (J.LE.NJOINT(III))
      qincr = .true.
      DO K=1,NJOINT(igp)
      IF (JBONDS(III,J).EQ.JBONDS(igp,K)) THEN

C
C add to bond exclusion list
      IBOND=JBONDS(III,J)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=IB(IBOND)
      JBXCL(NBNDXCL)=JB(IBOND)
C
      WRITE(6,'(17A)')
     @'    %atoms "', SEGID(IB(IBOND)),'-',
     @RESID(IB(IBOND)),'-',RES(IB(IBOND)),
     @'-',TYPE(IB(IBOND)),'" and "',
     @SEGID(JB(IBOND)),'-',RESID(JB(IBOND)),
     @'-',RES(JB(IBOND)),'-',TYPE(JB(IBOND)),'"'
      END IF
C
C
      SKIPJOINT = K
      qincr = .false.
      END IF
      END DO
C
      if (.not.qincr) then
      JOINTS = 0
      DO K=1,NJOINT(igp)
      IF (K.NE.SKIPJOINT) THEN
      JOINTS = JOINTS+1
      JATOMS(igp,JOINTS) = JATOMS(igp,K)
      JBONDS(igp,JOINTS) = JBONDS(igp,K)
      CONGRP(igp,JOINTS) = CONGRP(igp,K)
      END IF
      END DO
      NJOINT(igp) = JOINTS
      IF (ONJOINT(igp).GE.3) TOTCHN = TOTCHN - 1
      JOINTS = 0
      DO K=1,NJOINT(III)
      IF (K.NE.IJOINT) THEN
      JOINTS = JOINTS + 1
      JATOMS(III,JOINTS) = JATOMS(III,K)
      JBONDS(III,JOINTS) = JBONDS(III,K)
      CONGRP(III,JOINTS) = CONGRP(III,K)
      END IF
      END DO
      NJOINT(III) = JOINTS
      TOTCHN = TOTCHN - 1
      else
      J = J + 1
      end if
      END DO
      ELSE
C
C Not a closed loop.  Start a new chain.
      NCHAIN=NCHAIN+1
C
      MXCHAINS=MAX(MXCHAINS,NCHAIN)
C
      IF (NCHAIN.GT.MAXCHN) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A)') 
     & ' Please send a bug report to brunger@stanford.edu'
      WRITE(6,'(2A)') 
     & ' Please increase the "maximum number of chains"', 
     & ' (torsion_maxchain) and re-run.'
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
      NCHAINS(NTREE) = NCHAINS(NTREE)+1
      GRPSEQ(NTREE,NCHAIN,1)=III
      GRPSEQ(NTREE,NCHAIN,2)=CONGRP(III,IJOINT)
      TAKEN(GRPSEQ(NTREE,NCHAIN,2)) = 1
      NTAKEN = NTAKEN + 1
      CHNLEN(NTREE,NCHAIN)=2
C
      QLOOP = .FALSE.
C
C Same as before.  Continue adding bodies.
      DO WHILE
     @     ((NJOINT(GRPSEQ(NTREE,NCHAIN,CHNLEN(NTREE,NCHAIN))).EQ.2)
     @     .AND.(.NOT.QLOOP))
C
      THISONE = GRPSEQ(NTREE,NCHAIN,CHNLEN(NTREE,NCHAIN))
      LASTONE = GRPSEQ(NTREE,NCHAIN,CHNLEN(NTREE,NCHAIN)-1)
C
      DO I=1,2
      IF (CONGRP(THISONE,I).EQ.LASTONE) THEN
      IF (TAKEN(CONGRP(THISONE,3-I)).EQ.0) THEN
      NTAKEN = NTAKEN + 1
      CHNLEN(NTREE,NCHAIN)=CHNLEN(NTREE,NCHAIN)+1
C
      MXCHNLEN=MAX(MXCHNLEN,CHNLEN(NTREE,NCHAIN))
C
      IF (CHNLEN(NTREE,NCHAIN).GT.MAXLEN) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(2A)') 
     & ' Please increase the "maximum unbranched chain length"',
     & ' (torsion_maxlength) and re-run.'
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
      GRPSEQ(NTREE,NCHAIN,CHNLEN(NTREE,NCHAIN)) = CONGRP(THISONE,3-I)
      TAKEN(GRPSEQ(NTREE,NCHAIN,CHNLEN(NTREE,NCHAIN))) = 1
      ELSE
C
C A closed loop.  Fix up.
      QLOOP = .TRUE.
      QTOPO2 = .TRUE.
C
      IF (.NOT.QCLOSED) THEN
      WRITE(6,'(A,A)')
     @     ' -------------------------- Torsion Topology ',
     @     '-----------------------------------'
      WRITE(6,'(A,A)') ' Closed bonding networks have been detected.',
     @               ' The bonds between the following'
      WRITE(6,'(A/)') '   atoms will not be constrained:'
      QCLOSED = .TRUE.
      END IF
C
      DO J=1,NJOINT(THISONE)
      DO K=1,NJOINT(LASTONE)
      COND = ( (IB(JBONDS(LASTONE,K)).EQ.JATOMS(THISONE,J)).OR.
     @     (JB(JBONDS(LASTONE,K)).EQ.JATOMS(THISONE,J)) )
      IF (COND) THEN
      INATM = JATOMS(THISONE,J)
      INBOND = JBONDS(LASTONE,K)
      END IF
      END DO
      END DO
C
      DO J=1,NJOINT(THISONE)
      IF (JBONDS(THISONE,J).NE.INBOND) OUTATM = JATOMS(THISONE,J)
      END DO
C
      qskbnd = .false.
      DO J=1,NJOINT(CONGRP(THISONE,3-I))
      if (.not.qskbnd) then
      IF (IB(JBONDS(CONGRP(THISONE,3-I),J)).EQ.OUTATM) THEN
C
C add to bond exclusion list
      IBOND=JBONDS(CONGRP(THISONE,3-I),J)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=OUTATM
      JBXCL(NBNDXCL)=JB(IBOND)
C
      WRITE(6,'(17A)')
     @     '    %atoms "', SEGID(OUTATM),'-',RESID(OUTATM),'-',
     @     RES(OUTATM),'-',TYPE(OUTATM),'" and "',
     @     SEGID(JB(IBOND)),'-',
     @     RESID(JB(IBOND)),'-',
     @     RES(JB(IBOND)),'-',
     @     TYPE(JB(IBOND)),'"'
C
      END IF
      qskbnd = .true.
      SKIPJOINT = J
C
      ELSE IF  (JB(JBONDS(CONGRP(THISONE,3-I),J)).EQ.OUTATM) THEN
C
C add to bond exclusion list
      IBOND=JBONDS(CONGRP(THISONE,3-I),J)
      IF (XCLBND(IBOND).EQ.0) THEN
      NBNDXCL=NBNDXCL+1
      IF (NBNDXCL.GT.NBOND) THEN
      CALL WRNDIE(-5,'TOPOLOGY',
     @     'NBNDXCL exceeds NBONDS. Fatal coding error')
      END IF
      XCLBND(IBOND)=1
      IBXCL(NBNDXCL)=OUTATM
      JBXCL(NBNDXCL)=IB(IBOND)
C
      WRITE(6,'(17A)')
     @     '    %atoms "', SEGID(OUTATM),'-',RESID(OUTATM),'-',
     @     RES(OUTATM),'-',TYPE(OUTATM),'" and "',
     @     SEGID(IB(IBOND)),'-',
     @     RESID(IB(IBOND)),'-',
     @     RES(IB(IBOND)),'-',
     @     TYPE(IB(IBOND)),'"'
C
      END IF
C
      qskbnd = .true.
      SKIPJOINT = J
      END IF
      end if
      END DO
      qskbnd = .false.
      JOINTS = 0
      DO J=1,NJOINT(CONGRP(THISONE,3-I))
      IF (J.NE.SKIPJOINT) THEN
      JOINTS = JOINTS+1
      JATOMS(CONGRP(THISONE,3-I),JOINTS) = JATOMS(CONGRP(THISONE,3-I),J)
      JBONDS(CONGRP(THISONE,3-I),JOINTS) = JBONDS(CONGRP(THISONE,3-I),J)
      CONGRP(CONGRP(THISONE,3-I),JOINTS) = CONGRP(CONGRP(THISONE,3-I),J)
      END IF
      END DO
      NJOINT(CONGRP(THISONE,3-I)) = JOINTS
      IF (ONJOINT(CONGRP(THISONE,3-I)).GE.2) TOTCHN = TOTCHN - 1
C
      NJOINT(THISONE) = 1
      JATOMS(THISONE,1) = INATM
      JBONDS(THISONE,1) = INBOND
      CONGRP(THISONE,1) = LASTONE
      congrp(thisone,2) = 0
      END IF
      END IF
      END DO
      END DO
      END IF
      END IF
      IJOINT = IJOINT + 1
      END DO
      END IF
      END IF
      END DO
      END DO
C
      END DO
C
      IF (QCLOSED) THEN
      WRITE(6,'(A,A)')
     @     ' ---------------------------------------',
     @     '----------------------------------------'
      WRITE(6,'(A)') ' '
      END IF
C
C now order the important atoms for each group
C i.e. fill NIN, NOUT, which basically duplicate JATOMS
C in a form more suitable for later use.
      DO ITREE=1,NTREE
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=1,CHNLEN(ITREE,ICHAIN)
      NIN(ITREE,ICHAIN,IGP) = 0
      NOUT(ITREE,ICHAIN,IGP) = 0
      END DO
      END DO
      END DO
C
      DO ITREE=1,NTREE
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
C
      DO J=1,NJOINT(GRPSEQ(ITREE,ICHAIN,IGP))
      DO K=1,NJOINT(GRPSEQ(ITREE,ICHAIN,IGP-1))
C
      COND=((IB(JBONDS(GRPSEQ(ITREE,ICHAIN,IGP-1),K)).EQ.
     @     JATOMS(GRPSEQ(ITREE,ICHAIN,IGP),J)).OR.(
     @     JB(JBONDS(GRPSEQ(ITREE,ICHAIN,IGP-1),K)).EQ.
     @     JATOMS(GRPSEQ(ITREE,ICHAIN,IGP),J)))
C
      IF (COND) THEN
      NIN(ITREE,ICHAIN,IGP)=JATOMS(GRPSEQ(ITREE,ICHAIN,IGP),J)
      IF (NIN(ITREE,ICHAIN,IGP).EQ.
     @     IB(JBONDS(GRPSEQ(ITREE,ICHAIN,IGP-1),K))) THEN
      NOUT(ITREE,ICHAIN,IGP-1) =
     @     JB(JBONDS(GRPSEQ(ITREE,ICHAIN,IGP-1),K))
      ELSE
      NOUT(ITREE,ICHAIN,IGP-1) =
     @     IB(JBONDS(GRPSEQ(ITREE,ICHAIN,IGP-1),K))
      END IF
      END IF
C
      END DO
      END DO
C
      END DO
      END DO
      END DO
C
C     Print out topology information if desired
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') 'TREE TOPOLOGY FOLLOWS:'
      WRITE(6,'(A,I5)') 'The number of trees is ',NTREE
      DO ITREE=1,NTREE
      WRITE(6,'(A,I5,A,I5,A)')
     @     'Tree ',ITREE,' has ',NCHAINS(ITREE),' chains'
      WRITE(6,'(A,I5)') 'The base is group ', GRPSEQ(ITREE,1,1)
      DO ICHAIN=1,NCHAINS(ITREE)
      WRITE(6,'(A,I5)') 'Chain Number ',ICHAIN
      DO IGP=1,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      WRITE(6,'(A,I5,A,I5,A,I5,A,I5)') 'Link ',IGP,' is group ',III,
     @     ' INATM = ',NIN(ITREE,ICHAIN,IGP),
     @     ' OUTATM = ',NOUT(ITREE,ICHAIN,IGP)
      END DO
      END DO
      END DO
      END IF
C
C define an array to keep track of the size of each group
      DO ITREE=1,NTREE
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=1,CHNLEN(ITREE,ICHAIN)
      III=GRPSEQ(ITREE,ICHAIN,IGP)
      SIZE(III)=POINTR(III+1)-POINTR(III)-1
      END DO
      END DO
      END DO
C
C=========================
C TOPOLOGY SETUP COMPLETE
C=========================
C
      WRITE(6,'(A,I10)')
     &' TOPOLOGY: number of distinct bodies ',NTREE
      WRITE(6,'(A,I10)')
     &' TOPOLOGY: maximum number of chains (per distinct body) ',
     & MXCHAINS
      WRITE(6,'(A,I10)')
     &' TOPOLOGY: maximum unbranched chain length ',MXCHNLEN
C
      IF (QTOPO2) THEN
      WRITE(6,'(A)') 
     & ' You have reached this point due to a program error.'
      WRITE(6,'(A)') 
     & ' Reason: loops detected in TOPOLOGY.'
      WRITE(6,'(A)') 
     & ' Please send a bug report to brunger@stanford.edu'
      CALL WRNDIE(-5,'TORSION:TOPOLOGY',
     @     'Fatal Topology Error')
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE SETUP(NCHAINS,
     @     CHNLEN,COEFFI,COEFFO,NJOINT,
     @     SIZE,GRPSEQ,NIN,
     @     NOUT,QTIP)
C
C Routine uses the topology information to set up dynamics stuff.
C In particular, the following are calculated for each group:
C    COEFFI,COEFFO  Difficult to explain.  They are coefficients which multiply
C                   dynamical quantities of the NIN and NOUT atoms.  They are
C                   either 1.0 (what you would expect) or 0.5.  They are necessary
C                   because it is possible to have a group comprised of a single
C                   atom (i.e. no finite extent).  The problem is that the
C                   algorithm requires bodies with finite extent.  The solution
C                   (read: kludge) is to kind of share certain atoms between
C                   groups - give half of their mass and force contribution to
C                   two different groups.
C
C Author: Luke Rice
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
      INTEGER CHNLEN(MAXTREE,*),NCHAINS(*)
      INTEGER GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER SIZE(*),NIN(MAXTREE,MAXCHN,*),NOUT(MAXTREE,MAXCHN,*)
      INTEGER NJOINT(*),QTIP(*)
      DOUBLE PRECISION COEFFI(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION COEFFO(MAXTREE,MAXCHN,*)
C
C Local
      INTEGER III,IIIM,IIIP,ICHAIN,JCHAIN,TEMPCHN,TEMP
      INTEGER NNIN,NNINP,NNOUT,IGP,LENGTH,COUNT,ITREE,erratm
      LOGICAL QSKIPGP,QDONE
      DOUBLE PRECISION ONE,HALF
      PARAMETER (ONE = 1.0D0, HALF = 0.5D0)
C
C
      DO IGP=1,NGP
      QTIP(IGP) = 0
      END DO
C Loop over trees
      DO ITREE=1,NTREE
C
C First take care of base: the easy part
      III = GRPSEQ(ITREE,1,1)
      DO ICHAIN=1,NJOINT(III)
      COEFFO(ITREE,ICHAIN,1) = ONE
      COEFFI(ITREE,ICHAIN,1) = ONE
      END DO
C
C Now do the rest of the tree:
C   There are several possibilities, treated in the following order:
C      + First treat groups that are the second node in a chain with
C        a small base, because COEFFI might need special attention
C      + Then treat groups that are chain tips, again because the
C        coefficients might need special attention
C      + Then treat groups whose outboard body has just one joint.
C        This is again for treatment of the coefficients.
C      + Finally treat the rest.
C
C   The way to figure out coefficients is as follows:
C      + CASE 1: a big group
C         - COEFFO is always 1 (no sharing of the in atom of the next group)
C         - COEFFI is 1 (no sharing of the in atom) if
C            = the inboard body is big
C            = the inboard body is small, but shares with a different group
C              which group shares is decided by the higher chain number
C              This case only occurs if the inboard body is a base
C         - COEFFI is 0.5 (sharing of the in atom) if
C            = the inboard body is small, and needs to share
C              This is usually the case, unless the inboard body is a base
C              and shares with someone else (dictated by chain number)
C      + CASE 2: a small group
C         - COEFFO is 1 (absorb entirely the in atom of the next group) if
C             = the next group is small, and has only one joint
C               This makes sense, since the nexy guy has no one to share with.
C               In this case, the chain is shortened by one link.
C         - COEFFO is 0.5 (share the in atom of the next group) if
C             = the preceeding is not true (i.e. usually)
C         - COEFFI is 1 (no sharing of the in atom of this group) if
C             = the inboard body is big
C             = the inboard body is small, but shares with someone else
C         - COEFFI is 0.5 (share the in atom of this group) if
C            = the inboard body is small, and needs to share
C
      DO ICHAIN=1,NCHAINS(ITREE)
      QSKIPGP = .FALSE.
      LENGTH = CHNLEN(ITREE,ICHAIN)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      QDONE = .FALSE.
      IF (.NOT.QSKIPGP) THEN
C
      III = GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM = GRPSEQ(ITREE,ICHAIN,IGP-1)
      IF (IGP.NE.CHNLEN(ITREE,ICHAIN)) THEN
      IIIP = GRPSEQ(ITREE,ICHAIN,IGP+1)
      NNINP = NIN(ITREE,ICHAIN,IGP+1)
      ELSE
      IIIP = 0
      NNINP = 0
      END IF
      NNIN = NIN(ITREE,ICHAIN,IGP)
      NNOUT = NOUT(ITREE,ICHAIN,IGP-1)
C
      IF (IGP.EQ.2.AND.SIZE(IIIM).EQ.0) THEN
      IF (SIZE(III).GT.0) THEN
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      COUNT = 0
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      COUNT = COUNT + 1
      END IF
      END IF
      END DO
      IF (COUNT.EQ.0) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      TEMP = JCHAIN
      END IF
      END DO
      IF (ICHAIN.EQ.TEMP) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      ELSE IF (COUNT.GE.1) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      TEMP = JCHAIN
      END IF
      END IF
      END DO
      IF (ICHAIN.EQ.TEMP) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      END IF
      ELSE
C          SMALL GROUP
      IF (NJOINT(III).EQ.1.AND.IGP.EQ.CHNLEN(ITREE,ICHAIN)) THEN
      WRITE(6,'(A,A)')
     @        ' ---------------------------- Torsion Setup ',
     @        '------------------------------------'
      WRITE(6,'(A)')
     & '   ERROR: A forbidden topology has been detected.'
      ERRATM = nnin
      WRITE(6,'(9A)')
     @  '   Atom "', SEGID(ERRATM),'-',RESID(ERRATM),'-',
     @  RES(ERRATM),'-',TYPE(ERRATM),'" is involved.'
      WRITE(6,'(A/)') '            Please alter topology and re-run.'
      CALL WRNDIE(-5,'TORSION:SETUP',
     @               'Fatal Setup Error')
      END IF
      IF (IGP.EQ.CHNLEN(ITREE,ICHAIN)) THEN
      COUNT = 0
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      TEMPCHN = JCHAIN
      END IF
      END IF
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      COUNT = COUNT + 1
      END IF
      END IF
      END DO
      IF (ICHAIN.EQ.TEMPCHN) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      IF (COUNT.EQ.0) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      NNINP = NIN(ITREE,JCHAIN,2)
      NIN(ITREE,ICHAIN,IGP+1) = NNINP
      END IF
      END DO
      ELSE IF (COUNT.GE.1) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      NNINP = NIN(ITREE,JCHAIN,2)
      NIN(ITREE,ICHAIN,IGP+1) = NNINP
      END IF
      END IF
      END DO
      END IF
      ELSE
      COUNT = 0
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      COUNT = COUNT + 1
      END IF
      END IF
      END DO
      IF (COUNT.EQ.0) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      TEMPCHN = JCHAIN
      END IF
      END DO
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      IF (ICHAIN.EQ.TEMPCHN) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      ELSE IF (COUNT.GE.1) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (IIIM.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      TEMPCHN = JCHAIN
      END IF
      END IF
      END DO
      IF (NJOINT(IIIP).EQ.1) THEN
      IF (SIZE(IIIP).EQ.0) THEN
      QTIP(III) = 1
      QSKIPGP = .TRUE.
      CHNLEN(ITREE,ICHAIN) = CHNLEN(ITREE,ICHAIN) - 1
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      ELSE
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      END IF
      ELSE
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      END IF
      IF (ICHAIN.EQ.TEMPCHN) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      END IF
      END IF
      END IF
      QDONE = .TRUE.
      END IF
C
      IF (.NOT.QDONE.AND.IGP.EQ.CHNLEN(ITREE,ICHAIN)) THEN
      IF (SIZE(III).EQ.0.AND.NJOINT(III).EQ.1) THEN
      WRITE(6,'(A,A)')
     @        ' ---------------------------- Torsion Setup ',
     @        '------------------------------------'
      erratm = nnin
      WRITE(6,'(A)')
     &'   ERROR: A forbidden topology has been detected.'
      WRITE(6,'(9A)')
     @  '   Atom "', SEGID(ERRATM),'-',RESID(ERRATM),'-',
     @  RES(ERRATM),'-',TYPE(ERRATM),'" is involved.'
      WRITE(6,'(A/)') '            Please alter topology and re-run.'
      CALL WRNDIE(-5,'TORSION:SETUP',
     @               'Fatal Setup Error')
      END IF
      IF (SIZE(III).GT.0) THEN
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      IF (SIZE(IIIM).GT.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      END IF
      ELSE
      IF (SIZE(IIIM).GT.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      END IF
      COUNT = 0
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      COUNT = COUNT + 1
      END IF
      END IF
      END DO
      IF (COUNT.EQ.0) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      TEMP = GRPSEQ(ITREE,JCHAIN,2)
      NNINP = NIN(ITREE,JCHAIN,2)
      NIN(ITREE,ICHAIN,IGP+1) = NNINP
      END IF
      END DO
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      ELSE IF (COUNT.GE.1) THEN
      DO JCHAIN = 1,NCHAINS(ITREE)
      IF (III.EQ.GRPSEQ(ITREE,JCHAIN,1)) THEN
      IF (SIZE(GRPSEQ(ITREE,JCHAIN,2)).EQ.0) THEN
      TEMP = GRPSEQ(ITREE,JCHAIN,2)
      NNINP = NIN(ITREE,JCHAIN,2)
      NIN(ITREE,ICHAIN,IGP+1) = NNINP
      END IF
      END IF
      END DO
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      END IF
      END IF
      QDONE = .TRUE.
      END IF
C
      IF (.NOT.QDONE) THEN
      IF (NJOINT(IIIP).EQ.1) THEN
      IF (SIZE(III).GT.0) THEN
      IF (SIZE(IIIP).EQ.0) THEN
      WRITE(6,'(A,A)')
     @        ' ---------------------------- Torsion Setup ',
     @        '------------------------------------'
      WRITE(6,'(A)')
     &  '   ERROR: A forbidden topology has been detected.'
      erratm = NIN(ITREE,ICHAIN,IGP+1)
      WRITE(6,'(9A)')
     @  '   Atom "', SEGID(ERRATM),'-',RESID(ERRATM),'-',
     @  RES(ERRATM),'-',TYPE(ERRATM),'" is involved.'
      WRITE(6,'(A/)') '            Please alter topology and re-run.'
      CALL WRNDIE(-5,'TORSION:SETUP',
     @               'Fatal Setup Error')
      ELSE
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      END IF
      IF (SIZE(IIIM).EQ.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      ELSE
      IF (SIZE(IIIP).EQ.0) THEN
      QTIP(III) = 1
      QSKIPGP = .TRUE.
      CHNLEN(ITREE,ICHAIN) = CHNLEN(ITREE,ICHAIN) - 1
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      ELSE
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      END IF
      IF (SIZE(IIIM).EQ.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      END IF
      QDONE = .TRUE.
      END IF
C
      IF (.NOT.QDONE.AND.NJOINT(IIIP).GE.2) THEN
      IF (SIZE(III).GT.0) THEN
      COEFFO(ITREE,ICHAIN,IGP) = ONE
      IF (SIZE(IIIM).EQ.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      ELSE
      COEFFO(ITREE,ICHAIN,IGP) = HALF
      IF (SIZE(IIIM).EQ.0) THEN
      COEFFI(ITREE,ICHAIN,IGP) = HALF
      ELSE
      COEFFI(ITREE,ICHAIN,IGP) = ONE
      END IF
      END IF
      QDONE = .TRUE.
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
C
      IF (WRNLEV.GE.10) THEN
      WRITE(6,'(A)') 'More topology information'
      DO ITREE=1,NTREE
      WRITE(6,'(A,I5)') 'Tree number ',ITREE
      DO ICHAIN=1,NCHAINS(ITREE)
      WRITE(6,'(A,I5)') 'Chain number ',ICHAIN
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      WRITE(6,'(A,I5,A,A,F3.1,A,F3.1)') 'Link ',IGP,':',
     @     ' COEFIN = ' ,COEFFI(ITREE,ICHAIN,IGP),
     @     ' COEFOUT = ',COEFFO(ITREE,ICHAIN,IGP)
      END DO
      END DO
      END DO
      END IF
C
      RETURN
      END
C===============================================================
      SUBROUTINE SETUP2(LIST,POINTR,NCHAINS,
     @     CHNLEN,XCM,YCM,ZCM,COEFFI,COEFFO,
     @     MM,YY,Q0,Q1,Q2,Q3,SIZE,GRPSEQ,XB,YB,ZB,NIN,
     @     NOUT,QDOT,HHX,HHY,HHZ,SSX,SSY,SSZ,QTIP)
C
C Routine uses the topology information to set up dynamics stuff.
C In particular, the following are calculated for each group:
C    XCM,YCM,ZCM    Center of mass coordinates
C    Q0,Q1,Q2,Q3    Orientational coordinates (quaternions)
C    MM             6x6 inertia matrix
C    YY             6x1 velocity vector
C    XB,YB,ZB       Coordinates for each atom in the principal axis frame
C    HHX,HHY,HHZ    Components of the bond vector leading into this body
C    SSX,SSY,SSZ    Components of the vector from COM to NIN
C
      IMPLICIT NONE
C
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
      INCLUDE 'dtorsion.inc'
C
      INTEGER LIST(*),POINTR(*),CHNLEN(MAXTREE,*),NCHAINS(*)
      INTEGER GRPSEQ(MAXTREE,MAXCHN,*)
      INTEGER SIZE(*),NIN(MAXTREE,MAXCHN,*),NOUT(MAXTREE,MAXCHN,*)
      INTEGER QTIP(*)
      DOUBLE PRECISION XCM(*),YCM(*),ZCM(*),MM(6,6,NGP),QDOT(*)
      DOUBLE PRECISION Q0(*),Q1(*),Q2(*),Q3(*),COEFFI(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION XB(*),YB(*),ZB(*),YY(6,NGP)
      DOUBLE PRECISION COEFFO(MAXTREE,MAXCHN,*)
      DOUBLE PRECISION HHX(*),HHY(*),HHZ(*),SSX(*),SSY(*),SSZ(*)
C
C Local
      INTEGER III,IIIM,ICHAIN
      INTEGER NNIN,NNINP,NNOUT,IGP,ITREE
      LOGICAL QBASE,QBIG,QTIP2
      DOUBLE PRECISION ONE,HALF,COEFIN,COEFOU
      PARAMETER (ONE = 1.0D0, HALF = 0.5D0)
C
      NDEGF = 0
      DO ITREE=1,NTREE
      III = GRPSEQ(ITREE,1,1)
      QBASE = .TRUE.
      QBIG  = .TRUE.
      QTIP2 = .FALSE.
C
      CALL GRPROP(0,III,LIST,POINTR,XCM,YCM,ZCM,XB,YB,ZB,
     @     MM,NDEGF,Q0,Q1,Q2,Q3,YY,QBIG,QBASE,0,0,0,ONE,ONE,
     @     QTIP2,NGP,QDOT,HHX,HHY,HHZ,SSX,SSY,SSZ)
C
      DO ICHAIN=1,NCHAINS(ITREE)
      DO IGP=2,CHNLEN(ITREE,ICHAIN)
      III = GRPSEQ(ITREE,ICHAIN,IGP)
      IIIM = GRPSEQ(ITREE,ICHAIN,IGP-1)
      NNIN = NIN(ITREE,ICHAIN,IGP)
      NNINP = NIN(ITREE,ICHAIN,IGP+1)
      NNOUT = NOUT(ITREE,ICHAIN,IGP-1)
      COEFIN = COEFFI(ITREE,ICHAIN,IGP)
      COEFOU = COEFFO(ITREE,ICHAIN,IGP)
      QBASE = .FALSE.
      IF (SIZE(III).GT.0) THEN
      QBIG = .TRUE.
      ELSE
      QBIG = .FALSE.
      END IF
      IF (QTIP(III).EQ.0) THEN
      QTIP2 = .FALSE.
      ELSE
      QTIP2 = .TRUE.
      END IF
C
      CALL GRPROP(IIIM,III,LIST,POINTR,XCM,YCM,ZCM,XB,YB,ZB,
     @     MM,NDEGF,Q0,Q1,Q2,Q3,YY,QBIG,QBASE,NNIN,NNINP,NNOUT,
     @     COEFIN,COEFOU,QTIP2,NGP,QDOT,HHX,HHY,HHZ,SSX,SSY,SSZ)
C
      END DO
      END DO
      END DO
C
C
      RETURN
      END
C===============================================================
      SUBROUTINE GRPROP(IIIM,III,LIST,POINTR,XCM,
     @   YCM,ZCM,XB,YB,ZB,MM,NDEGF,Q0,Q1,Q2,Q3,YY,QBIG,
     @   QBASE,NNIN,NNINP,NNOUT,COEFIN,COEFOU,QTIP,NGP,QDOT,
     @   HHX,HHY,HHZ,SSX,SSY,SSZ)
C
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'coord.inc'
      INCLUDE 'deriv.inc'
      INCLUDE 'consta.inc'
      INCLUDE 'timer.inc'
C
      INTEGER III,LIST(*),POINTR(*),NDEGF,NNIN,NNINP
      INTEGER IIIM,NGP,NNOUT
      DOUBLE PRECISION XCM(*),YCM(*),ZCM(*),XB(*),YB(*)
      DOUBLE PRECISION ZB(*),MM(6,6,NGP),Q0(*),Q1(*),Q2(*)
      DOUBLE PRECISION Q3(*),YY(6,NGP),COEFIN,COEFOU,QDOT(*)
      DOUBLE PRECISION HHX(*),HHY(*),HHZ(*),SSX(*),SSY(*),SSZ(*)
      LOGICAL QBIG,QBASE,QTIP
C
C Local
      INTEGER POINTER,K,IERR,L,J, IND(3)
      DOUBLE PRECISION MASS,XREL,YREL,ZREL,A(6),B(3,9),DET
      DOUBLE PRECISION ZERO,ROOT(3),VECT(3,3),SMALL
      DOUBLE PRECISION IXP,IYP,IZP,ROT(3,3),JX,JY,JZ,WX,WY,WZ
      DOUBLE PRECISION VCMX,VCMY,VCMZ,Q02,Q12,Q22,Q32,ONE,QSQR
      DOUBLE PRECISION NORM,HNORM,QUART,A11,A12,rrelp2,dot
      DOUBLE PRECISION A13,A21,A22,A23,A31,A32,A33,JBX,JBY,JBZ
      DOUBLE PRECISION XRELP,YRELP,ZRELP,TWO,random,rrel2
      DOUBLE PRECISION XDUM, YDUM, ZDUM, MDUM
      DOUBLE PRECISION XDIFF, YDIFF, ZDIFF
      LOGICAL COND, Q2BASE
C parameter
      PARAMETER (ONE=1.0D0,QUART=0.25D0,ZERO=0.0D0)
      PARAMETER (SMALL=1.0D-2,TWO=2.0D0)
C begin
      XCM(III)=ZERO
      YCM(III)=ZERO
      ZCM(III)=ZERO
      MASS=ZERO
      JX=ZERO
      JY=ZERO
      JZ=ZERO
      VCMX=ZERO
      VCMY=ZERO
      VCMZ=ZERO
      DO J=1,6
      A(J)=ZERO
      END DO
C
      DO J=1,3
      DO K=1,9
      B(J,K)=ZERO
      END DO
      END DO
C
      IF (QBASE) THEN
C
C base body
C
C if the base body has only two atoms, then we define a dummy
C third atom with a small mass ATB 3/28/10
      IF (POINTR(III+1)-(POINTR(III)+1).EQ.1) THEN
        Q2BASE=.TRUE.
C
C define a vector that is orthogonal to the vector between
C the two atoms (cross product between the difference vector and the unit vector 1,0,0)
        XDIFF=X(LIST(POINTR(III+1)))-X(LIST(POINTR(III)+1))
        YDIFF=Y(LIST(POINTR(III+1)))-Y(LIST(POINTR(III)+1))
        ZDIFF=Z(LIST(POINTR(III+1)))-Z(LIST(POINTR(III)+1))
        XDUM=ZERO
        YDUM=ZDIFF
        ZDUM=-YDIFF
C
C normalize the vector
        NORM=MAX(R4SMAL,SQRT(XDUM**2+YDUM**2+ZDUM**2))
        XDUM=XDUM/NORM
        YDUM=YDUM/NORM
        ZDUM=ZDUM/NORM
C
C add it to the first atom
        XDUM=XDUM+X(LIST(POINTR(III)+1))
        YDUM=YDUM+Y(LIST(POINTR(III)+1))
        ZDUM=ZDUM+Z(LIST(POINTR(III)+1))
C
C set the mass of the dummy atom to a small number
        MDUM=SMALL
      ELSE
       Q2BASE=.FALSE.
      END IF
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        XCM(III)=XCM(III)+X(LIST(POINTER))*AMASS(LIST(POINTER))
        YCM(III)=YCM(III)+Y(LIST(POINTER))*AMASS(LIST(POINTER))
        ZCM(III)=ZCM(III)+Z(LIST(POINTER))*AMASS(LIST(POINTER))
        MASS=MASS+AMASS(LIST(POINTER))
      END DO
C
      IF (Q2BASE) THEN
        XCM(III)=XCM(III)+XDUM*MDUM
        YCM(III)=YCM(III)+YDUM*MDUM
        ZCM(III)=ZCM(III)+ZDUM*MDUM
        MASS=MASS+MDUM
      END IF
C
      IF (MASS.LT.SMALL) THEN
        WRITE(6,'(A)') 'FATAL ERROR: MASS OF A BODY IS ZERO'
      ELSE
        XCM(III)=XCM(III)/MASS
        YCM(III)=YCM(III)/MASS
        ZCM(III)=ZCM(III)/MASS
      END IF
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        XREL=X(LIST(POINTER))-XCM(III)
        YREL=Y(LIST(POINTER))-YCM(III)
        ZREL=Z(LIST(POINTER))-ZCM(III)
        A(1)=A(1)+AMASS(LIST(POINTER))*(YREL**2+ZREL**2)
        A(2)=A(2)-AMASS(LIST(POINTER))*XREL*YREL
        A(3)=A(3)+AMASS(LIST(POINTER))*(XREL**2+ZREL**2)
        A(4)=A(4)-AMASS(LIST(POINTER))*ZREL*XREL
        A(5)=A(5)-AMASS(LIST(POINTER))*YREL*ZREL
        A(6)=A(6)+AMASS(LIST(POINTER))*(XREL**2+YREL**2)
      END DO
C
      IF (Q2BASE) THEN
C add dummy atom to mass tensor
        XREL=XDUM-XCM(III)
        YREL=YDUM-YCM(III)
        ZREL=ZDUM-ZCM(III)
        A(1)=A(1)+MDUM*(YREL**2+ZREL**2)
        A(2)=A(2)-MDUM*XREL*YREL
        A(3)=A(3)+MDUM*(XREL**2+ZREL**2)
        A(4)=A(4)-MDUM*ZREL*XREL
        A(5)=A(5)-MDUM*YREL*ZREL
        A(6)=A(6)+MDUM*(XREL**2+YREL**2)      
      END IF
C
C
      DO K=1,6
      IF (ABS(A(K)).LT.R4SMAL) THEN
        A(K)=ZERO
      END IF
      END DO
C
C get principal moments, and define body coordinates
      CALL GIVEIS(3,3,3,A,B,IND,ROOT,VECT,IERR)
      IF (IERR.NE.0) THEN
        WRITE(6,'(A,I3)')
     @' %GRROP-ERR: diagonalization of A failed, IERR= ', IERR
        CALL WRNDIE(-5,'GRPROP','Diagonalization failed.')
      END IF
      IXP=ROOT(1)
      IYP=ROOT(2)
      IZP=ROOT(3)
      DO K=1,3
      DO L=1,3
        ROT(K,L)=VECT(L,K)
      END DO
      END DO
      DET=(ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3)
     @    +(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3)
     @    +(ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      DO K=1,3
      DO L=1,3
        ROT(K,L)=-VECT(L,K)
      END DO
      END DO
      END IF
C
      JX=ZERO
      JY=ZERO
      JZ=ZERO
      WX=ZERO
      WY=ZERO
      WZ=ZERO
      VCMX=ZERO
      VCMY=ZERO
      VCMZ=ZERO
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        XREL=X(LIST(POINTER))-XCM(III)
        YREL=Y(LIST(POINTER))-YCM(III)
        ZREL=Z(LIST(POINTER))-ZCM(III)
        XB(LIST(POINTER))=XREL*ROT(1,1)+YREL*ROT(1,2)
     @                   +ZREL*ROT(1,3)
        YB(LIST(POINTER))=XREL*ROT(2,1)+YREL*ROT(2,2)
     @                   +ZREL*ROT(2,3)
        ZB(LIST(POINTER))=XREL*ROT(3,1)+YREL*ROT(3,2)
     @                   +ZREL*ROT(3,3)
        JX=JX+AMASS(LIST(POINTER))*(YREL*ZV(LIST(POINTER))
     @                             -ZREL*YV(LIST(POINTER)))
        JY=JY+AMASS(LIST(POINTER))*(ZREL*XV(LIST(POINTER))
     @                             -XREL*ZV(LIST(POINTER)))
        JZ=JZ+AMASS(LIST(POINTER))*(XREL*YV(LIST(POINTER))
     @                             -YREL*XV(LIST(POINTER)))
        VCMX=VCMX+AMASS(LIST(POINTER))*XV(LIST(POINTER))
        VCMY=VCMY+AMASS(LIST(POINTER))*YV(LIST(POINTER))
        VCMZ=VCMZ+AMASS(LIST(POINTER))*ZV(LIST(POINTER))
      END DO
C
C fill generalized mass matrix
      DO J=1,6
         DO K=1,6
            MM(K,J,III)=0.
         END DO
      END DO
      DO J=1,3
         MM(J,J,III)=MASS
      END DO
          MM(4,4,III)=ROT(1,1)*ROT(1,1)*IXP+ROT(2,1)*ROT(2,1)
     @               *IYP+ROT(3,1)*ROT(3,1)*IZP
          MM(5,4,III)=ROT(1,1)*ROT(1,2)*IXP+ROT(2,1)*ROT(2,2)
     @               *IYP+ROT(3,1)*ROT(3,2)*IZP
          MM(6,4,III)=ROT(1,1)*ROT(1,3)*IXP+ROT(2,1)*ROT(2,3)
     @               *IYP+ROT(3,1)*ROT(3,3)*IZP
          MM(4,5,III)=ROT(1,1)*ROT(1,2)*IXP+ROT(2,1)*ROT(2,2)
     @               *IYP+ROT(3,1)*ROT(3,2)*IZP
          MM(5,5,III)=ROT(1,2)*ROT(1,2)*IXP+ROT(2,2)*ROT(2,2)
     @               *IYP+ROT(3,2)*ROT(3,2)*IZP
          MM(6,5,III)=ROT(1,2)*ROT(1,3)*IXP+ROT(2,2)*ROT(2,3)
     @               *IYP+ROT(3,2)*ROT(3,3)*IZP
          MM(4,6,III)=ROT(1,1)*ROT(1,3)*IXP+ROT(2,1)*ROT(2,3)
     @               *IYP+ROT(3,1)*ROT(3,3)*IZP
          MM(5,6,III)=ROT(1,2)*ROT(1,3)*IXP+ROT(2,2)*ROT(2,3)
     @               *IYP+ROT(3,2)*ROT(3,3)*IZP
          MM(6,6,III)=ROT(1,3)*ROT(1,3)*IXP+ROT(2,3)*ROT(2,3)
     @               *IYP+ROT(3,3)*ROT(3,3)*IZP
C
C get number of degrees of freedom
      IF (IXP.LT.SMALL.OR.IYP.LT.
     @       SMALL.OR.IZP.LT.SMALL) THEN
        NDEGF=NDEGF+5
      ELSE
        NDEGF=NDEGF+6
      END IF
C
      Q02=QUART*(ONE+ROT(1,1)+ROT(2,2)+ROT(3,3))
      Q12=QUART*(ONE+ROT(1,1)-ROT(2,2)-ROT(3,3))
      Q22=QUART*(ONE-ROT(1,1)+ROT(2,2)-ROT(3,3))
      Q32=QUART*(ONE-ROT(1,1)-ROT(2,2)+ROT(3,3))
      IF (ABS(Q02).GT.R4SMAL) THEN
        Q0(III)=SQRT(Q02)
        Q1(III)=QUART*(ROT(2,3)-ROT(3,2))/Q0(III)
        Q2(III)=QUART*(ROT(3,1)-ROT(1,3))/Q0(III)
        Q3(III)=QUART*(ROT(1,2)-ROT(2,1))/Q0(III)
      ELSE IF (ABS(Q12).GT.R4SMAL) THEN
        Q1(III)=SQRT(Q12)
        Q0(III)=QUART*(ROT(2,3)-ROT(3,2))/Q1(III)
        Q2(III)=QUART*(ROT(1,2)+ROT(2,1))/Q1(III)
        Q3(III)=QUART*(ROT(1,3)+ROT(3,1))/Q1(III)
      ELSE IF (ABS(Q22).GT.R4SMAL) THEN
        Q2(III)=SQRT(Q22)
        Q0(III)=QUART*(ROT(3,1)-ROT(1,3))/Q2(III)
        Q1(III)=QUART*(ROT(1,2)+ROT(2,1))/Q2(III)
        Q3(III)=QUART*(ROT(2,3)+ROT(3,2))/Q2(III)
      ELSE IF (ABS(Q32).GT.R4SMAL) THEN
        Q3(III)=SQRT(Q32)
        Q0(III)=QUART*(ROT(1,2)-ROT(2,1))/Q3(III)
        Q1(III)=QUART*(ROT(1,3)+ROT(3,1))/Q3(III)
        Q2(III)=QUART*(ROT(2,3)+ROT(3,2))/Q3(III)
      ELSE
        WRITE(6,'(A)') 
     &   '%GRPROP-ERR: matrix ROT is not a rotation matrix'
        CALL WRNDIE(-5,'GRPROP','ROT is not a rotation matrix.')
      END IF
      QSQR=Q0(III)*Q0(III)+Q1(III)*Q1(III)+
     @     Q2(III)*Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
      A11=Q0(III)*Q0(III)+Q1(III)*Q1(III)-
     @    Q2(III)*Q2(III)-Q3(III)*Q3(III)
      A21=TWO*(Q1(III)*Q2(III)-Q0(III)*Q3(III))
      A31=TWO*(Q1(III)*Q3(III)+Q0(III)*Q2(III))
      A12=TWO*(Q1(III)*Q2(III)+Q0(III)*Q3(III))
      A22=Q0(III)*Q0(III)-Q1(III)*Q1(III)+
     @    Q2(III)*Q2(III)-Q3(III)*Q3(III)
      A32=TWO*(Q2(III)*Q3(III)-Q0(III)*Q1(III))
      A13=TWO*(Q1(III)*Q3(III)-Q0(III)*Q2(III))
      A23=TWO*(Q2(III)*Q3(III)+Q0(III)*Q1(III))
      A33=Q0(III)*Q0(III)-Q1(III)*Q1(III)-
     @    Q2(III)*Q2(III)+Q3(III)*Q3(III)
      COND=
     @     ABS(ROT(1,1)-A11).GT.SMALL
     @     .OR.ABS(ROT(1,2)-A12).GT.SMALL
     @     .OR.ABS(ROT(1,3)-A13).GT.SMALL
     @     .OR.ABS(ROT(2,1)-A21).GT.SMALL
     @     .OR.ABS(ROT(2,2)-A22).GT.SMALL
     @     .OR.ABS(ROT(2,3)-A23).GT.SMALL
     @     .OR.ABS(ROT(3,1)-A31).GT.SMALL
     @     .OR.ABS(ROT(3,2)-A32).GT.SMALL
     @     .OR.ABS(ROT(3,3)-A33).GT.SMALL
      IF (COND) THEN
        WRITE(6,*)
     @     ' %GRPOP-ERR: inconsistent QS '
        WRITE(6,'(a,i5)') ' GROUP = ',III
        WRITE(6,*) ' QBASE=',QBASE
        WRITE(6,*) ' ROT=',ROT
        WRITE(6,*) ' A=',A11,A12,A13,A21,A22,A23,A31,A32,A33
        CALL WRNDIE(-5,'GRPROP','Inconsistent QS.')
      END IF
      JBX=A11*JX+A12*JY+A13*JZ
      JBY=A21*JX+A22*JY+A23*JZ
      JBZ=A31*JX+A32*JY+A33*JZ
      IF (IXP.LT.SMALL) THEN
        WX=ZERO
      ELSE
        WX=JBX/IXP
      END IF
      IF (IYP.LT.SMALL) THEN
        WY=ZERO
      ELSE
        WY=JBY/IYP
      END IF
      IF (IZP.LT.SMALL) THEN
        WZ=ZERO
      ELSE
        WZ=JBZ/IZP
      END IF
C
C base com velocities
      YY(1,III)=VCMX/MASS
      YY(2,III)=VCMY/MASS
      YY(3,III)=VCMZ/MASS
      YY(4,III)=A11*WX+A21*WY+A31*WZ
      YY(5,III)=A12*WX+A22*WY+A32*WZ
      YY(6,III)=A13*WX+A23*WY+A33*WZ
C
C Now update initial velocities so the call to TEMPOPT will be ok
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A11
     @          +YB(LIST(POINTER))*A21+ZB(LIST(POINTER))*A31
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A12
     @          +YB(LIST(POINTER))*A22+ZB(LIST(POINTER))*A32
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A13
     @          +YB(LIST(POINTER))*A23+ZB(LIST(POINTER))*A33
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      ELSE IF (.NOT.QBASE) THEN
      IF (QBIG) THEN
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        IF (LIST(POINTER).NE.NNIN) THEN
        XCM(III)=XCM(III)+X(LIST(POINTER))*AMASS(LIST(POINTER))
        YCM(III)=YCM(III)+Y(LIST(POINTER))*AMASS(LIST(POINTER))
        ZCM(III)=ZCM(III)+Z(LIST(POINTER))*AMASS(LIST(POINTER))
        MASS=MASS+AMASS(LIST(POINTER))
        ELSE
        XCM(III)=XCM(III)+X(LIST(POINTER))*
     @                    AMASS(LIST(POINTER))*COEFIN
        YCM(III)=YCM(III)+Y(LIST(POINTER))*
     @                    AMASS(LIST(POINTER))*COEFIN
        ZCM(III)=ZCM(III)+Z(LIST(POINTER))*
     @                    AMASS(LIST(POINTER))*COEFIN
        MASS=MASS+AMASS(LIST(POINTER))*COEFIN
        END IF
      END DO
      IF (MASS.LT.SMALL) THEN
        WRITE(6,'(A)') 'FATAL ERROR: MASS OF A BODY IS ZERO'
      ELSE
        XCM(III)=XCM(III)/MASS
        YCM(III)=YCM(III)/MASS
        ZCM(III)=ZCM(III)/MASS
      END IF
C
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        XREL=X(LIST(POINTER))-XCM(III)
        YREL=Y(LIST(POINTER))-YCM(III)
        ZREL=Z(LIST(POINTER))-ZCM(III)
        IF (LIST(POINTER).NE.NNIN) THEN
        A(1)=A(1)+AMASS(LIST(POINTER))*(YREL**2+ZREL**2)
        A(2)=A(2)-AMASS(LIST(POINTER))*XREL*YREL
        A(3)=A(3)+AMASS(LIST(POINTER))*(XREL**2+ZREL**2)
        A(4)=A(4)-AMASS(LIST(POINTER))*ZREL*XREL
        A(5)=A(5)-AMASS(LIST(POINTER))*YREL*ZREL
        A(6)=A(6)+AMASS(LIST(POINTER))*(XREL**2+YREL**2)
        ELSE
        A(1)=A(1)+(YREL**2+ZREL**2)*
     @                       AMASS(LIST(POINTER))*COEFIN
        A(2)=A(2)-XREL*YREL*AMASS(LIST(POINTER))*COEFIN
        A(3)=A(3)+(XREL**2+ZREL**2)*
     @                       AMASS(LIST(POINTER))*COEFIN
        A(4)=A(4)-ZREL*XREL*AMASS(LIST(POINTER))*COEFIN
        A(5)=A(5)-YREL*ZREL*AMASS(LIST(POINTER))*COEFIN
        A(6)=A(6)+(XREL**2+YREL**2)*
     @                       AMASS(LIST(POINTER))*COEFIN
        END IF
      END DO
C
      DO K=1,6
      IF (ABS(A(K)).LT.R4SMAL) THEN
        A(K)=ZERO
      END IF
      END DO
C
C get principal moments, and define body coordinates
      CALL GIVEIS(3,3,3,A,B,IND,ROOT,VECT,IERR)
      IF (IERR.NE.0) THEN
        WRITE(6,'(A,I3)')
     @' %GRPROP-ERR: diagonalization failed, IERR= ', IERR
        CALL WRNDIE(-5,'GRPROP','Diagonalization failed.')
      END IF
      IXP=ROOT(1)
      IYP=ROOT(2)
      IZP=ROOT(3)
C
C test to see if GIVEIS choked - happens sometimes if IXP=0
      IF (IXP.LT.R4SMAL) THEN
         VECT(1,1) = VECT(2,2)*VECT(3,3)-VECT(3,2)*VECT(2,3)
         VECT(2,1) = VECT(3,2)*VECT(1,3)-VECT(1,2)*VECT(3,3)
         VECT(3,1) = VECT(1,2)*VECT(2,3)-VECT(2,2)*VECT(1,3)
      END IF
C
C make sure vect is unitary
      DO K=1,3
      DO L=1,3
        ROT(K,L)=VECT(L,K)
      END DO
      END DO
      DET=(ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3)
     @    +(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3)
     @    +(ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)
      IF (ABS(DET-ONE).GT.SMALL) THEN
      DO K=1,3
      DO L=1,3
        ROT(K,L)=-VECT(L,K)
      END DO
      END DO
      END IF
C
      JX=ZERO
      JY=ZERO
      JZ=ZERO
      WX=ZERO
      WY=ZERO
      WZ=ZERO
      VCMX=ZERO
      VCMY=ZERO
      VCMZ=ZERO
      DO POINTER=POINTR(III)+1,POINTR(III+1)
        XREL=X(LIST(POINTER))-XCM(III)
        YREL=Y(LIST(POINTER))-YCM(III)
        ZREL=Z(LIST(POINTER))-ZCM(III)
        XB(LIST(POINTER))=XREL*ROT(1,1)+YREL*ROT(1,2)
     @                   +ZREL*ROT(1,3)
        YB(LIST(POINTER))=XREL*ROT(2,1)+YREL*ROT(2,2)
     @                   +ZREL*ROT(2,3)
        ZB(LIST(POINTER))=XREL*ROT(3,1)+YREL*ROT(3,2)
     @                   +ZREL*ROT(3,3)
        IF (LIST(POINTER).NE.NNIN) THEN
        JX=JX+AMASS(LIST(POINTER))*(YREL*ZV(LIST(POINTER))
     @                             -ZREL*YV(LIST(POINTER)))
        JY=JY+AMASS(LIST(POINTER))*(ZREL*XV(LIST(POINTER))
     @                             -XREL*ZV(LIST(POINTER)))
        JZ=JZ+AMASS(LIST(POINTER))*(XREL*YV(LIST(POINTER))
     @                             -YREL*XV(LIST(POINTER)))
        VCMX=VCMX+AMASS(LIST(POINTER))*XV(LIST(POINTER))
        VCMY=VCMY+AMASS(LIST(POINTER))*YV(LIST(POINTER))
        VCMZ=VCMZ+AMASS(LIST(POINTER))*ZV(LIST(POINTER))
        ELSE
        JX=JX+AMASS(LIST(POINTER))*(YREL*ZV(LIST(POINTER))
     @                     -ZREL*YV(LIST(POINTER)))*COEFIN
        JY=JY+AMASS(LIST(POINTER))*(ZREL*XV(LIST(POINTER))
     @                     -XREL*ZV(LIST(POINTER)))*COEFIN
        JZ=JZ+AMASS(LIST(POINTER))*(XREL*YV(LIST(POINTER))
     @                     -YREL*XV(LIST(POINTER)))*COEFIN
        VCMX=VCMX+AMASS(LIST(POINTER))*COEFIN*
     @                                 XV(LIST(POINTER))
        VCMY=VCMY+AMASS(LIST(POINTER))*COEFIN*
     @                                 YV(LIST(POINTER))
        VCMZ=VCMZ+AMASS(LIST(POINTER))*COEFIN*
     @                                 ZV(LIST(POINTER))
        END IF
      END DO
C
      ELSE
C
      XCM(III)=COEFIN*AMASS(NNIN)*X(NNIN)+
     @         COEFOU*AMASS(NNINP)*X(NNINP)
      YCM(III)=COEFIN*AMASS(NNIN)*Y(NNIN)+
     @         COEFOU*AMASS(NNINP)*Y(NNINP)
      ZCM(III)=COEFIN*AMASS(NNIN)*Z(NNIN)+
     @         COEFOU*AMASS(NNINP)*Z(NNINP)
      MASS=COEFIN*AMASS(NNIN)+COEFOU*AMASS(NNINP)
           IF (MASS.LT.SMALL) THEN
            WRITE(6,'(A)') 'FATAL ERROR: MASS OF A BODY IS ZERO'
           ELSE
            XCM(III)=XCM(III)/MASS
            YCM(III)=YCM(III)/MASS
            ZCM(III)=ZCM(III)/MASS
           END IF
         XREL=X(NNIN)-XCM(III)
         YREL=Y(NNIN)-YCM(III)
         ZREL=Z(NNIN)-ZCM(III)
         XRELP=X(NNINP)-XCM(III)
         YRELP=Y(NNINP)-YCM(III)
         ZRELP=Z(NNINP)-ZCM(III)
C
         vect(1,1)=x(nninp)-x(nnin)
         vect(2,1)=y(nninp)-y(nnin)
         vect(3,1)=z(nninp)-z(nnin)
         norm=vect(1,1)**2+vect(2,1)**2+vect(3,1)**2
         norm=sqrt(norm)
         do k=1,3
         vect(k,1)=vect(k,1)/norm
         end do
C
C now get a perpendicular vector
      do k=1,3
      call ggubfs(random)
      vect(k,2)=random
      end do
C
      dot=zero
      do k=1,3
      dot=dot+vect(k,1)*vect(k,2)
      end do
C
      do k=1,3
      vect(k,2)=vect(k,2)-dot*vect(k,1)
      end do
      norm=vect(1,2)**2+vect(2,2)**2+vect(3,2)**2
      norm=sqrt(norm)
      do k=1,3
      vect(k,2)=vect(k,2)/norm
      end do
C
C finally the third vector
         vect(1,3)=vect(2,1)*vect(3,2)-vect(2,2)*vect(3,1)
         vect(2,3)=vect(3,1)*vect(1,2)-vect(3,2)*vect(1,1)
         vect(3,3)=vect(1,1)*vect(2,2)-vect(1,2)*vect(2,1)
C
C normalize
         norm=vect(1,3)**2+vect(2,3)**2+vect(3,3)**2
         norm=sqrt(norm)
         do k=1,3
         vect(k,3)=vect(k,3)/norm
         end do
C
C now get the moment of inertia about y and z
         ixp=zero
         iyp=zero
         izp=zero
         XREL=X(NNIN)-XCM(III)
         YREL=Y(NNIN)-YCM(III)
         ZREL=Z(NNIN)-ZCM(III)
         rrel2=xrel**2+yrel**2+zrel**2
         iyp=iyp+coefin*amass(nnin)*rrel2
         izp=izp+coefin*amass(nnin)*rrel2
         XRELP=X(NNINP)-XCM(III)
         YRELP=Y(NNINP)-YCM(III)
         ZRELP=Z(NNINP)-ZCM(III)
         rrelp2=xrelp**2+yrelp**2+zrelp**2
         iyp=iyp+coefou*amass(nninp)*rrelp2
         izp=izp+coefou*amass(nninp)*rrelp2
C
C double check that vect is unitary
         DO K=1,3
         DO L=1,3
           ROT(K,L)=VECT(L,K)
         END DO
         END DO
         DET=(ROT(1,1)*ROT(2,2)-ROT(1,2)*ROT(2,1))*ROT(3,3)
     @    +(ROT(2,1)*ROT(3,2)-ROT(2,2)*ROT(3,1))*ROT(1,3)
     @    +(ROT(3,1)*ROT(1,2)-ROT(3,2)*ROT(1,1))*ROT(2,3)
         IF (ABS(DET-ONE).GT.SMALL) THEN
         DO K=1,3
         DO L=1,3
           ROT(K,L)=-VECT(L,K)
         END DO
         END DO
         END IF
C
         XB(NNIN)=XREL*ROT(1,1)+YREL*ROT(1,2)
     @                    +ZREL*ROT(1,3)
         YB(NNIN)=XREL*ROT(2,1)+YREL*ROT(2,2)
     @                    +ZREL*ROT(2,3)
         ZB(NNIN)=XREL*ROT(3,1)+YREL*ROT(3,2)
     @                    +ZREL*ROT(3,3)
C
         IF (QTIP) THEN
         XB(NNINP)=XRELP*ROT(1,1)+YRELP*ROT(1,2)
     @                    +ZRELP*ROT(1,3)
         YB(NNINP)=XRELP*ROT(2,1)+YRELP*ROT(2,2)
     @                    +ZRELP*ROT(2,3)
         ZB(NNINP)=XRELP*ROT(3,1)+YRELP*ROT(3,2)
     @                    +ZRELP*ROT(3,3)
         END IF
C
         JX=COEFIN*AMASS(NNIN)*(YREL*ZV(NNIN)-ZREL*YV(NNIN))+
     @    COEFOU*AMASS(NNINP)*(YRELP*ZV(NNINP)-ZRELP*YV(NNINP))
         JY=COEFIN*AMASS(NNIN)*(ZREL*XV(NNIN)-XREL*ZV(NNIN))+
     @    COEFOU*AMASS(NNINP)*(ZRELP*XV(NNINP)-XRELP*ZV(NNINP))
         JZ=COEFIN*AMASS(NNIN)*(XREL*YV(NNIN)-YREL*XV(NNIN))+
     @    COEFOU*AMASS(NNINP)*(XRELP*YV(NNINP)-YRELP*XV(NNINP))
         VCMX=COEFIN*AMASS(NNIN)*XV(NNIN)+
     @        COEFOU*AMASS(NNINP)*XV(NNINP)
         VCMY=COEFIN*AMASS(NNIN)*YV(NNIN)+
     @        COEFOU*AMASS(NNINP)*YV(NNINP)
         VCMZ=COEFIN*AMASS(NNIN)*ZV(NNIN)+
     @        COEFOU*AMASS(NNINP)*ZV(NNINP)
C
      END IF
C
C fill generalized mass matrix
      DO J=1,6
         DO K=1,6
            MM(K,J,III)=0.
         END DO
      END DO
      DO J=1,3
         MM(J,J,III)=MASS
      END DO
          MM(4,4,III)=ROT(1,1)*ROT(1,1)*IXP+ROT(2,1)*ROT(2,1)
     @               *IYP+ROT(3,1)*ROT(3,1)*IZP
          MM(5,4,III)=ROT(1,1)*ROT(1,2)*IXP+ROT(2,1)*ROT(2,2)
     @               *IYP+ROT(3,1)*ROT(3,2)*IZP
          MM(6,4,III)=ROT(1,1)*ROT(1,3)*IXP+ROT(2,1)*ROT(2,3)
     @               *IYP+ROT(3,1)*ROT(3,3)*IZP
          MM(4,5,III)=ROT(1,1)*ROT(1,2)*IXP+ROT(2,1)*ROT(2,2)
     @               *IYP+ROT(3,1)*ROT(3,2)*IZP
          MM(5,5,III)=ROT(1,2)*ROT(1,2)*IXP+ROT(2,2)*ROT(2,2)
     @               *IYP+ROT(3,2)*ROT(3,2)*IZP
          MM(6,5,III)=ROT(1,2)*ROT(1,3)*IXP+ROT(2,2)*ROT(2,3)
     @               *IYP+ROT(3,2)*ROT(3,3)*IZP
          MM(4,6,III)=ROT(1,1)*ROT(1,3)*IXP+ROT(2,1)*ROT(2,3)
     @               *IYP+ROT(3,1)*ROT(3,3)*IZP
          MM(5,6,III)=ROT(1,2)*ROT(1,3)*IXP+ROT(2,2)*ROT(2,3)
     @               *IYP+ROT(3,2)*ROT(3,3)*IZP
          MM(6,6,III)=ROT(1,3)*ROT(1,3)*IXP+ROT(2,3)*ROT(2,3)
     @               *IYP+ROT(3,3)*ROT(3,3)*IZP
C
C get number of degrees of freedom
        NDEGF=NDEGF+1
C
      Q02=QUART*(ONE+ROT(1,1)+ROT(2,2)+ROT(3,3))
      Q12=QUART*(ONE+ROT(1,1)-ROT(2,2)-ROT(3,3))
      Q22=QUART*(ONE-ROT(1,1)+ROT(2,2)-ROT(3,3))
      Q32=QUART*(ONE-ROT(1,1)-ROT(2,2)+ROT(3,3))
      IF (ABS(Q02).GT.R4SMAL) THEN
        Q0(III)=SQRT(Q02)
        Q1(III)=QUART*(ROT(2,3)-ROT(3,2))/Q0(III)
        Q2(III)=QUART*(ROT(3,1)-ROT(1,3))/Q0(III)
        Q3(III)=QUART*(ROT(1,2)-ROT(2,1))/Q0(III)
      ELSE IF (ABS(Q12).GT.R4SMAL) THEN
        Q1(III)=SQRT(Q12)
        Q0(III)=QUART*(ROT(2,3)-ROT(3,2))/Q1(III)
        Q2(III)=QUART*(ROT(1,2)+ROT(2,1))/Q1(III)
        Q3(III)=QUART*(ROT(1,3)+ROT(3,1))/Q1(III)
      ELSE IF (ABS(Q22).GT.R4SMAL) THEN
        Q2(III)=SQRT(Q22)
        Q0(III)=QUART*(ROT(3,1)-ROT(1,3))/Q2(III)
        Q1(III)=QUART*(ROT(1,2)+ROT(2,1))/Q2(III)
        Q3(III)=QUART*(ROT(2,3)+ROT(3,2))/Q2(III)
      ELSE IF (ABS(Q32).GT.R4SMAL) THEN
        Q3(III)=SQRT(Q32)
        Q0(III)=QUART*(ROT(1,2)-ROT(2,1))/Q3(III)
        Q1(III)=QUART*(ROT(1,3)+ROT(3,1))/Q3(III)
        Q2(III)=QUART*(ROT(2,3)+ROT(3,2))/Q3(III)
      ELSE
        WRITE(6,'(A)') 'ERROR: matrix is not a rotation matrix'
      END IF
      QSQR=Q0(III)*Q0(III)+Q1(III)*Q1(III)+
     @     Q2(III)*Q2(III)+Q3(III)*Q3(III)
      NORM=SQRT(QSQR)
      Q0(III)=Q0(III)/NORM
      Q1(III)=Q1(III)/NORM
      Q2(III)=Q2(III)/NORM
      Q3(III)=Q3(III)/NORM
      A11=Q0(III)*Q0(III)+Q1(III)*Q1(III)-
     @    Q2(III)*Q2(III)-Q3(III)*Q3(III)
      A21=TWO*(Q1(III)*Q2(III)-Q0(III)*Q3(III))
      A31=TWO*(Q1(III)*Q3(III)+Q0(III)*Q2(III))
      A12=TWO*(Q1(III)*Q2(III)+Q0(III)*Q3(III))
      A22=Q0(III)*Q0(III)-Q1(III)*Q1(III)+
     @    Q2(III)*Q2(III)-Q3(III)*Q3(III)
      A32=TWO*(Q2(III)*Q3(III)-Q0(III)*Q1(III))
      A13=TWO*(Q1(III)*Q3(III)-Q0(III)*Q2(III))
      A23=TWO*(Q2(III)*Q3(III)+Q0(III)*Q1(III))
      A33=Q0(III)*Q0(III)-Q1(III)*Q1(III)-
     @    Q2(III)*Q2(III)+Q3(III)*Q3(III)
      COND=
     @     ABS(ROT(1,1)-A11).GT.SMALL
     @     .OR.ABS(ROT(1,2)-A12).GT.SMALL
     @     .OR.ABS(ROT(1,3)-A13).GT.SMALL
     @     .OR.ABS(ROT(2,1)-A21).GT.SMALL
     @     .OR.ABS(ROT(2,2)-A22).GT.SMALL
     @     .OR.ABS(ROT(2,3)-A23).GT.SMALL
     @     .OR.ABS(ROT(3,1)-A31).GT.SMALL
     @     .OR.ABS(ROT(3,2)-A32).GT.SMALL
     @     .OR.ABS(ROT(3,3)-A33).GT.SMALL
      IF (COND) THEN
        WRITE(6,*)
     @     ' %GRPROP-ERR: inconsistent qs ',ROT,
     @       A11,A12,A13,A21,A22,A23,A31,A32,A33
      END IF
      JBX=A11*JX+A12*JY+A13*JZ
      JBY=A21*JX+A22*JY+A23*JZ
      JBZ=A31*JX+A32*JY+A33*JZ
      IF (IXP.LT.SMALL) THEN
        WX=ZERO
      ELSE
        WX=JBX/IXP
      END IF
      IF (IYP.LT.SMALL) THEN
        WY=ZERO
      ELSE
        WY=JBY/IYP
      END IF
      IF (IZP.LT.SMALL) THEN
        WZ=ZERO
      ELSE
        WZ=JBZ/IZP
      END IF
C
C Fill the YY vector correctly
      HHX(III)=X(NNIN)-X(NNOUT)
      HHY(III)=Y(NNIN)-Y(NNOUT)
      HHZ(III)=Z(NNIN)-Z(NNOUT)
      HNORM=HHX(III)**2+HHY(III)**2+HHZ(III)**2
      HNORM=SQRT(HNORM)
C
      SSX(III)=X(NNIN)-XCM(III)
      SSY(III)=Y(NNIN)-YCM(III)
      SSZ(III)=Z(NNIN)-ZCM(III)
C
C Use lower half of YY for temporary storage
      YY(4,III)=A11*WX+A21*WY+A31*WZ
      YY(5,III)=A12*WX+A22*WY+A32*WZ
      YY(6,III)=A13*WX+A23*WY+A33*WZ
C
      QDOT(III)=(YY(4,III)-YY(4,IIIM))*HHX(III)/HNORM+
     @          (YY(5,III)-YY(5,IIIM))*HHY(III)/HNORM+
     @          (YY(6,III)-YY(6,IIIM))*HHZ(III)/HNORM
C
      YY(4,III)=YY(4,IIIM)+QDOT(III)*HHX(III)/HNORM
      YY(5,III)=YY(5,IIIM)+QDOT(III)*HHY(III)/HNORM
      YY(6,III)=YY(6,IIIM)+QDOT(III)*HHZ(III)/HNORM
C
C Now fill the translational part
      YY(1,III)=YY(1,IIIM)+(-ZCM(IIIM)+ZCM(III))*YY(5,IIIM)+
     @         (YCM(IIIM)-YCM(III))*YY(6,IIIM)+
     @         ((HHZ(III)*SSY(III)-HHY(III)*SSZ(III))/HNORM)*QDOT(III)
      YY(2,III)=YY(2,IIIM)+(ZCM(IIIM)-ZCM(III))*YY(4,IIIM)+
     @         (-XCM(IIIM)+XCM(III))*YY(6,IIIM)+
     @         ((HHX(III)*SSZ(III)-HHZ(III)*SSX(III))/HNORM)*QDOT(III)
      YY(3,III)=YY(3,IIIM)+(-YCM(IIIM)+YCM(III))*YY(4,IIIM)+
     @         (XCM(IIIM)-XCM(III))*YY(5,IIIM)+
     @         ((HHY(III)*SSX(III)-HHX(III)*SSY(III))/HNORM)*QDOT(III)
C
C Update atom coordinates, velocities
      DO POINTER=POINTR(III)+1,POINTR(III+1)
      X(LIST(POINTER))=XCM(III)+XB(LIST(POINTER))*A11
     @          +YB(LIST(POINTER))*A21+ZB(LIST(POINTER))*A31
      Y(LIST(POINTER))=YCM(III)+XB(LIST(POINTER))*A12
     @          +YB(LIST(POINTER))*A22+ZB(LIST(POINTER))*A32
      Z(LIST(POINTER))=ZCM(III)+XB(LIST(POINTER))*A13
     @          +YB(LIST(POINTER))*A23+ZB(LIST(POINTER))*A33
      XREL=X(LIST(POINTER))-XCM(III)
      YREL=Y(LIST(POINTER))-YCM(III)
      ZREL=Z(LIST(POINTER))-ZCM(III)
      XV(LIST(POINTER))=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(LIST(POINTER))=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(LIST(POINTER))=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END DO
C
      IF (QTIP) THEN
      X(NNINP)=XCM(III)+XB(NNINP)*A11
     @          +YB(NNINP)*A21+ZB(NNINP)*A31
      Y(NNINP)=YCM(III)+XB(NNINP)*A12
     @          +YB(NNINP)*A22+ZB(NNINP)*A32
      Z(NNINP)=ZCM(III)+XB(NNINP)*A13
     @          +YB(NNINP)*A23+ZB(NNINP)*A33
      XREL=X(NNINP)-XCM(III)
      YREL=Y(NNINP)-YCM(III)
      ZREL=Z(NNINP)-ZCM(III)
      XV(NNINP)=YY(1,III)+YY(5,III)*ZREL-YY(6,III)*YREL
      YV(NNINP)=YY(2,III)+YY(6,III)*XREL-YY(4,III)*ZREL
      ZV(NNINP)=YY(3,III)+YY(4,III)*YREL-YY(5,III)*XREL
      END IF
C
      END IF
C
      RETURN
      END
C===============================================================
