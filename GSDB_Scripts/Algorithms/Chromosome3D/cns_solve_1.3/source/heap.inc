C
C heap.inc
C ========
C
C reference pointers for dynamic memory allocation
      INTEGER                    HEAP(0:15)
      CHARACTER*1              C1HEAP(0:15)
      CHARACTER*4              C4HEAP(0:1)
      CHARACTER*80            C80HEAP(0:1)
      CHARACTER*(WORD_SIZE)   CWSHEAP(0:1)
      CHARACTER*(STRING_SIZE) CSSHEAP(0:1)
      EQUIVALENCE(C1HEAP(0), C4HEAP(0), C80HEAP(0),
     &            CWSHEAP(0), CSSHEAP(0))
C
C common blocks
      COMMON /IHEAP/ HEAP
      COMMON /SHEAP/ C1HEAP
