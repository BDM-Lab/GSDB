      SUBROUTINE SCRATC
C
C After changes in atom numbers, fixing atoms, this routine
C is called to reset the facilities (fragile atom selections)
C that are not mapped during atom number changes.
C
C Author: Axel Brunger.
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'nbonds.inc'
      INCLUDE 'xcrystal.inc'
      INCLUDE 'xtarget.inc'
      INCLUDE 'mtf.inc'
      INCLUDE 'update.inc'
      INCLUDE 'cnst.inc'
      INCLUDE 'vector.inc'
      INCLUDE 'heap.inc'
      INCLUDE 'funct.inc'
      INCLUDE 'noe.inc'
      INCLUDE 'ncs.inc'
      INCLUDE 'plane.inc'
      INCLUDE 'learn.inc'
      INCLUDE 'trtraj.inc'
      INCLUDE 'mmdisg.inc'
      INCLUDE 'dtorsion.inc'
      INCLUDE 'timer.inc'
C
C local
      INTEGER N, I
      LOGICAL QWRITE
C begin
C
C set bond, angle, dihedral, improper type-based parameter update flags
      UPBOND=.TRUE.
      UPANGL=.TRUE.
      UPDIHE=.TRUE.
      UPIMPR=.TRUE.
C
C reset all nonbonded list flags, i.e., generate exclusion and 1-4
C list, atom pair list.
      DO N=1,NPIG
      UPNBLS(N)=.TRUE.
      UPNBSS(N)=.TRUE.
      UPNBEX(N)=.TRUE.
      ENDDO
C set the nonbonded type-based parameter lookup flag
      UPNBLK=.TRUE.
C
C reset x-ray diffraction atomic scatterers
      IF (XRSN.GT.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: XRAY SCATter database erased.'
      END IF
      XRUPAT=.TRUE.
      XRQCHK=.TRUE.
C
C put the Pair of Interacting Group lists back to their default values
      IF (QPIGRST) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: Pairs of Interacting Groups erased.'
      END IF
      CALL PIGRSET('RESE')
C
C reset dihedral angle restraints list
      IF (CNPHI.GT.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: RESTraints DIHEdral database erased.'
      END IF
      CALL CDIRES
C
C reset NOE restraints list
      IF (NOENUM.GT.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: NOE restraints database erased.'
      END IF
      CALL NOERES
C
C delete the DG bounds matrices if appropriate
      IF (DUBPTR.NE.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: distance geometry bounds matrix erased.'
      END IF
      CALL DGRESZ(0)
C
C reset the parameter learning facility
      IF (PRIND.NE.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: PARAmeter LEARn selection erased.'
      END IF
      CALL PRTFRE
C
C reset the trajectory write facility
      IF (TRIND.NE.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: WRITe TRAJectory selection erased.'
      END IF
      CALL TRTFRE
C
C free and reset NCS group lists
      IF (HPSPNT(1).NE.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: NCS RESTraints database erased.'
      END IF
      CALL NCSFIN
      CALL NCSINR
C
C free PLN group lists
      IF (HPSPNP(1).NE.0.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: RESTraints PLANar database erased.'
      END IF
      CALL PLNFIN
      CALL PLNINR
C
C free all torsion angle lists
      CALL TORHP(0,.TRUE.)
      CALL TORINI
C
C free buffer stores and re-initialize
      QWRITE=.FALSE.
      DO I=1,MAXSTO
      IF (PTRSTO(I).NE.0) THEN
      CALL FREHP(PTRSTO(I),IREAL8(LENSTO(I)))
      QWRITE=.TRUE.
      END IF
      PTRSTO(I)=0
      LENSTO(I)=0
      END DO
      IF (QWRITE.AND.WRNLEV.GE.5) THEN
      WRITE(6,'(A)')
     & ' SCRATC-warning: STORe selections erased.'
      END IF
C
C note: the FOR symbol in ID statement selection is left alone
C to allow for conservative patches (see warning in routine parser.s).
C
C=====================================================================
C #if defined(CNS_SOLVE_COMPILE)
C=====================================================================
      CALL SCRCSFT
      CALL SCRHSFT
      CALL SCRCUPP
      CALL SCRONEB
      CALL SCRRAMA
      CALL SCRANGDB
      CALL SCRADANI
      CALL SCRASANI
C=====================================================================
C #endif
C=====================================================================
      RETURN
      END
C
