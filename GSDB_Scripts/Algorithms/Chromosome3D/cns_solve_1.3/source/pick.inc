C
C pick.inc
C ========
C
C Parameters and variables for picking bond, angle, dihedral
C properties.
C
      INTEGER PCGEOM, PCENER, PCEQUI, PCCONS, PCDEVI, PCPERI
      PARAMETER (PCGEOM=1, PCENER=2, PCEQUI=3, PCCONS=4, PCDEVI=5)
      PARAMETER (PCPERI=6)
      DOUBLE PRECISION PCDATA(6)
C
C parameters for nonbonded distance analysis routines
      DOUBLE PRECISION RCUTON, RCUTOF
      CHARACTER*4 SDISP
      INTEGER RUNIT, NISEL, NJSEL, RESDELTA, RESLOW, NDISTANCE
C pointer
      INTEGER IISEL, IJSEL, IMAPI, IMAPJ, ISMAT
C
C double precision common block
      COMMON /PICKC/ PCDATA, RCUTON, RCUTOF
      SAVE /PICKC/
C
C integer common block
      COMMON /PICKI/ RUNIT, NISEL, NJSEL, IISEL, IJSEL, IMAPI,
     &               IMAPJ, ISMAT, RESDELTA, RESLOW, NDISTANCE
C
C character common block
      COMMON /PICKS/ SDISP
