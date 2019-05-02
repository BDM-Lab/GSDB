C======================================================================
      SUBROUTINE READD
C
C parse READ qualifiers
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
      INCLUDE 'mtf.inc'
C begin
      CALL NEXTWD('READ-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-read')
C read dynamics file
      ELSE IF (WD(1:4).EQ.'DYNA'.OR.WD(1:4).EQ.'TRAJ') THEN
      CALL READDP
      ELSE
      CALL DSPERR('CNS','unknown qualifier')
      END IF
      RETURN
      END
C======================================================================
      SUBROUTINE WRITT
C
C parses WRITE qualifiers
      IMPLICIT NONE
C I/O
      INCLUDE 'cns.inc'
      INCLUDE 'comand.inc'
C local
      CHARACTER*4 FORMAT
C begin
C
      FORMAT='CNS'
C
      CALL NEXTWD('WRITE-qualifier=')
C
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write')
C
C======================================================================
      ELSE IF (WD(1:4).EQ.'STRU') THEN
C
C defaults
      OFILE='OUTPUT'
C
      CALL PUSEND('WRITe-STRUcture>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe-STRUcture>')
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-structure')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE
      CALL CHKEND('WRITe-STRUcture>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL MTFWRT(OFILE)
C======================================================================
      ELSE IF (WD(1:4).EQ.'PSF ') THEN
C
C defaults
      OFILE='OUTPUT'
C
      CALL PUSEND('WRITe-PSF>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe-PSF>')
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-psf')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE
      CALL CHKEND('WRITe-PSF>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL PSFWRT(OFILE,0)
C======================================================================
      ELSE IF (WD(1:4).EQ.'PARA') THEN
C
C defaults
      OFILE='OUTPUT'
C
      CALL PUSEND('WRITe-PARAmeter>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe-PARAmeter>')
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-parameter')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE IF (WD(1:4).EQ.'FORM') THEN
      CALL NEXTA4('FORMat=',FORMAT)
      IF (FORMAT(1:3).NE.'CNS'.AND.
     &    FORMAT(1:4).NE.'MMCI') THEN
         CALL WRNDIE(1,'WRITE','Unknown format specified.')
         FORMAT='CNS'
      END IF
      ELSE
      CALL CHKEND('WRITe-PARAmeter>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      IF (FORMAT(1:3).EQ.'CNS') THEN
         CALL PARWTR(OFILE,'NORMAL')
      ELSEIF (FORMAT(1:4).EQ.'MMCI') THEN
         CALL PARWTR(OFILE,'VERBOSE')
      END IF
C======================================================================
      ELSE IF (WD(1:4).EQ.'TOPO') THEN
C
C defaults
      OFILE='OUTPUT'
C
      CALL PUSEND('WRITe-TOPOlogy>')
      DO WHILE (.NOT.DONE)
      CALL NEXTWD('WRITe-TOPOlogy>')
      IF (WD(1:4).EQ.'HELP') THEN
C
      CALL CNSHELP('cns-write-topology')
C
      ELSE IF (WD(1:4).EQ.'OUTP') THEN
      CALL NEXTFI('OUTPut=',OFILE)
      ELSE
      CALL CHKEND('WRITe-TOPOlogy>',DONE)
      END IF
      END DO
      DONE=.FALSE.
C
      CALL RTWRIT(OFILE)
C======================================================================
      ELSE IF (WD(1:4).EQ.'COOR') THEN
      CALL CWRITE
C======================================================================
      ELSE IF (WD(1:4).EQ.'TRAJ') THEN
      CALL TRTDYN
C======================================================================
      ELSE
      CALL DSPERR('WRITE','unknown qualifier')
      END IF
      RETURN
      END
C======================================================================
