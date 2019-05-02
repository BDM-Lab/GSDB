#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>

void outbuf_(void)
{
  /* routine to call setlinebuf and update disc after each write
   */

  setlinebuf(stdout);
}

typedef struct
  {
    double re;
    double im;
  }
  zomplex;

void zfft3difree_(int *nx, int *ny, int *nz, zomplex *work, int *error)
{
#if defined(_ABIN32) || defined(_ABI64)
  int             *addlw;
  struct utsname  name;

  *error = 1;

  if (uname(&name) < 0) return;
  if (   strcmp(name.release, "6.3") == 0
      || strcmp(name.release, "6.4") == 0)
  {
    addlw = *(((int **)(work + *nx)) - 1);
    if(*addlw != *nx) return;
    free(addlw);

    addlw = *(((int **)(work + *nx + *ny)) - 1);
    if(*addlw != *ny) return;
    free(addlw);

    addlw = *(((int **)(work + *nx + *ny + *nz)) - 1);
    if(*addlw != *nz) return;
    free(addlw);
  }
#endif

  *error = 0;
}
