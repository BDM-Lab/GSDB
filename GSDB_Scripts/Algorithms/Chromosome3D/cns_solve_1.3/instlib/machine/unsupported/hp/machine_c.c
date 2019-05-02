#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/param.h>
#include <sys/times.h>


void outbuf(void)
{
  /* routine to call setlinebuf and update disc after each write
   */

  setvbuf(stdout, NULL, _IOLBF, 0);
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


void getnv(char *name, char *alias, int *length)
{
  /* routine to get value of environment variable name and
   * its length. T. Simonson, Sept. 1992
   */

  char *result;

  result=getenv(name);
  *length=strlen(result);
  if (*length>0) {strcpy(alias,result);};
}


int vdate(char *a, int *amax, int *alen)
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
      sprintf(a, "%2d-%3s-%2.2d",
        (*LocalTime).tm_mday, month,
        (*LocalTime).tm_year % 100 );
      *alen = 9;
  }

  return 0;
}


int vtime(char a[], int *amax, int *alen)
{
  /* this routine returns day time */
  /* +++++++++++++++++++++++++++++++++++++++ */
  /* machine dependent NeXT/Unix  C  version */
  /* H.Treutlein                    May 1991 */
  /* +++++++++++++++++++++++++++++++++++++++ */

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
      sprintf(a, "%2d:%2d:%2d",
        (*LocalTime).tm_hour,
        (*LocalTime).tm_min,
        (*LocalTime).tm_sec);
      *alen = 8;
  }

  return 0;
}


int getlog(char *usernm, int *len)
{
  /* Gets the username associated with the current process and places */
  /* it in USERNM. */
  /* +++++++++++++++++++++++++++++++++++++++ */
  /* machine dependent NeXT/Unix  C  version */
  /* H.Treutlein                June 5, 1991 */
  /* +++++++++++++++++++++++++++++++++++++++ */

  /* local */
  struct passwd *user;

  user =  getpwuid(getuid());
  strcpy(usernm, user->pw_name);
  *len = strlen(user->pw_name);

  return 0;
}


float etime(float *tarray)
{
  struct tms buffer;
  times(&buffer);
  tarray[0] = (double)(buffer.tms_utime)/CLK_TCK;
  tarray[1] = (double)(buffer.tms_stime)/CLK_TCK;
  return (tarray[0] + tarray[1]);
}
