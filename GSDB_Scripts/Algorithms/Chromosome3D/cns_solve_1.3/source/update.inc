C update.inc
C ==========
C
C flags that indicate to retrieve type-based parameters
C for bonds, angles, dihedrals, and impropers
       LOGICAL UPBOND, UPANGL, UPDIHE, UPIMPR
C
       COMMON /UPCOD/ UPBOND, UPANGL, UPDIHE, UPIMPR
C
       SAVE /UPCOD/
