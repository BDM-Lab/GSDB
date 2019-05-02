C  consta.inc
C
C this file contains all physical and mathematical constants
C and conversion factors.
C
C at present the following units are used:
C
C   length: Angstroms
C   time: ps
C   energy: Kcal/mol
C   mass: atomic-mass-unit
C   charge: electron-charge
C
C
      DOUBLE PRECISION RSMALL
      PARAMETER (RSMALL=1.0D-10)
      DOUBLE PRECISION R4SMAL,R4BIG
      PARAMETER (R4SMAL=0.0001D0,R4BIG=1.0D+10)
C
C physical constants in SI units
C ------------------------------
C     Kb = 1.380662 E-23 J/K
C     Na = 6.022045 E23  1/mol
C     e = 1.6021892 E-19 C
C     eps = 8.85418782 E-12 F/m
C
C     1 Kcal = 4184.0 J
C     1 amu = 1.6605655 E-27 Kg
C     1 A = 1.0 E-10 m
C
C reference: CRC Handbook for Chemistry and Physics, 1983/84
C
C
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535898D0)
C
C     TIMFAC is the conversion factor from AKMA time to picoseconds.
C            (TIMFAC = SQRT ( ( 1A )**2 * 1amu * Na  / 1Kcal )
C            this factor is used only intrinsically, all I/O is in ps.
C
      DOUBLE PRECISION TIMFAC
      PARAMETER (TIMFAC=0.04888821D0)
C
C KBOLTZ is Boltzman constant AKMA units (KBOLTZ = N *K  / 1 Kcal)
C                                                   a  b
      DOUBLE PRECISION KBOLTZ
      PARAMETER (KBOLTZ=1.987191D-03)
C
C CCELEC is 1/ (4 pi eps ) in AKMA units, conversion from SI
C units: CCELEC = e*e*Na / (4*pi*eps*1Kcal*1A)
C
      DOUBLE PRECISION CCELEC
      PARAMETER (CCELEC=332.0636D0)
C
C CDEBHU is used in the Debye-Hueckel approximation:
C      DIV GRAD phi = kappa**2 phi
C      kappa**2 = CDEBHU * ionic_strength [M] / ( T [K] eps   )
C                                                          ext
C where CDEBHU is defined as CDEBHU=2E+3 Na e**2 / (eps0 Kb )
C (in SI units, ref: Gordon M.Barrow, Physical Chemistry,
C McGraw Hill (1979) ) and ionic_strength is given in molar units.
C The conversion to AKMA units brings another factor 1.0E-20.
C
      DOUBLE PRECISION CDEBHU
      PARAMETER (CDEBHU=2529.09702D0)
