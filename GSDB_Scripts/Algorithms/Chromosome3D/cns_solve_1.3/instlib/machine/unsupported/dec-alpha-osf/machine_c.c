#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/param.h>
#include <sys/time.h>

void outbuf_(void)
{
  /* routine to call setlinebuf and update disc after each write
   */

  setlinebuf(stdout);
}

int fsatty_(int *desc)
{
  return isatty(*desc);
}

void fhostname_(char *hostnm, int *maxlen, int *hnlen)
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

int vdate_(char *a, int *amax, int *alen)
{
  /* this routine returns the date */
  /* +++++++++++++++++++++++++++++++++++++++ */
  /* machine dependent NeXT/Unix  C  version */
  /* H.Treutlein                    May 1991 */
  /* +++++++++++++++++++++++++++++++++++++++ */
  
  /* local */

  struct timeval tim;
  struct timezone zone;
  struct tm *LocalTime;
  char month[3];
  
  /* begin */
  if (*amax < 9) {
    printf(" %VDATE-ERR: string to small");
  } else {
    gettimeofday(&tim, &zone);
    LocalTime=localtime(&tim.tv_sec);
    switch (LocalTime->tm_mon){
    case 0: strcpy(month, "Jan");
      break;
    case 1: strcpy(month, "Feb");
      break;
    case 2: strcpy(month, "Mar");
      break;
    case 3: strcpy(month, "Apr");
      break;
    case 4: strcpy(month, "May");
      break;
    case 5: strcpy(month, "Jun");
      break;
    case 6: strcpy(month, "Jul");
      break;
    case 7: strcpy(month, "Aug");
      break;
    case 8: strcpy(month, "Sep");
      break;
    case 9: strcpy(month, "Oct");
      break;
    case 10: strcpy(month, "Nov");
      break;
    case 11: strcpy(month, "Dec");
      break;
    }
    sprintf(a, "%2d-%3s-%2.2d  ",
	    (*LocalTime).tm_mday, month,
	    (*LocalTime).tm_year % 100 );
    *alen = 9;
  }
  
  return 0;
}
