C Maximum number of PStores
      INTEGER   MAXPSTORES
      PARAMETER(MAXPSTORES=10)
C Dimensions of Flagmap when PStore was filled
      INTEGER   RNPSTORE(3, MAXPSTORES)
C Asymmetric unit when PStore was filled
      INTEGER   NAPSTORE(2, 3, MAXPSTORES)
C Number of peaks in PStores
      INTEGER    NPPSTORE(MAXPSTORES)
C Pointers for PStores
      INTEGER    PTPSTORE(MAXPSTORES)
C
      COMMON /PSTINT/ NPPSTORE, RNPSTORE, NAPSTORE
      COMMON /PSTPTR/ PTPSTORE
C
      SAVE /PSTINT/
      SAVE /PSTPTR/
