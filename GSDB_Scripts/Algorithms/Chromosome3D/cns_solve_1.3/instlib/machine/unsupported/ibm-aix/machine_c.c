#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/time.h>


void outbuf(void)
{
}


int fsatty(int *desc)
{
	return isatty(*desc);
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


int ibmmv(char *from, char *to)
{
/* Renames file from to to */
/* +++++++++++++++++++++++++++++++++++++++ */
/* machine dependent NeXT/Unix  C  version */
/* H.Treutlein               Sept 11, 1991 */
/* +++++++++++++++++++++++++++++++++++++++ */
    
    if (rename(from, to) != 0) perror("rename err: ");
    return 0;
} /* rename_ */


int vdate(char *a, int *amax, int *alen)
{


/* this routine returns the date */
/* +++++++++++++++++++++++++++++++++++++++ */
/* machine dependent NeXT/Unix  C  version */
/* H.Treutlein                    May 1991 */
/* +++++++++++++++++++++++++++++++++++++++ */

/* input/output */
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
	sprintf(a, "%2d-%3s-%2.2d",
               (*LocalTime).tm_mday, month,
               (*LocalTime).tm_year % 100 ); 
	*alen = 9;
    }
    return 0;
} /* vdate_ */


int vtime(char a[], int *amax, int *alen)
{


/* this routine returns day time */
/* +++++++++++++++++++++++++++++++++++++++ */
/* machine dependent NeXT/Unix  C  version */
/* H.Treutlein                    May 1991 */
/* +++++++++++++++++++++++++++++++++++++++ */

/* input/output */
/* local */
   struct timeval tim;
   struct timezone zone;
   struct tm *LocalTime;
/* begin */
    if (*amax < 8) {
	printf(" %VTIME-ERR: string to small");
    } else {
	*alen = *amax;
	gettimeofday(&tim, &zone);
	LocalTime=localtime(&tim.tv_sec);
	sprintf(a, "%2d:%2d:%2d", (*LocalTime).tm_hour, (*LocalTime).tm_min, (*LocalTime).tm_sec);
	*alen = 8;
    }
    return 0;
} /* vtime_ */
