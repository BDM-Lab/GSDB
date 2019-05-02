#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/param.h>


void outbuf_(void)
{
  /*
   *  routine to call setlinebuf and update disc after each write
   */

    setlinebuf(stdout);
}


void fhostname(char *hostnm, int *maxlen, int *hnlen)
{
  int   i;
  char  *lhn0, *lhn;

  lhn0 = lhn = malloc(MAXHOSTNAMELEN + 1);

  if (gethostname(lhn, MAXHOSTNAMELEN + 1) != 0)
    (void) strcpy(lhn, "unknown");

  *hnlen = *maxlen;

  for (i = 0; i < *maxlen; i++) {
    if (*lhn) {
      hostnm[i] = *lhn++;
      *hnlen = i + 1;
    }
    else
      hostnm[i] = ' ';
  }

  free(lhn0);
}
