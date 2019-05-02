#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/resource.h>


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


/* date */
/*===========================================================================*/
/* string from ctime looks like:
 * Sun Sep 16 01:03:52 1973\n\0
 * 01234567890123456789012345
 * 0         1          2
 *
 * we return:
 * dd-mmm-yy
 */
void date_( out_str, len )
char *out_str;
int  *len;
{
    long pvec[2];
    char temp[100];

    time(pvec);
    strcpy( temp, ctime(pvec));

    *out_str++ = temp[8];
    *out_str++ = temp[9];
    *out_str++ = '-';
    *out_str++ = temp[4];
    *out_str++ = temp[5];
    *out_str++ = temp[6];
    *out_str++ = '-';
    *out_str++ = temp[22];
    *out_str++ = temp[23];

} /*date*/
void fdate_( out_str, len )
char *out_str;
int  *len;
{
    long pvec[2];
    time(pvec);
    strcpy( out_str, ctime(pvec));
}

int csatty_(fildes)
int *fildes;
{
  return isatty(*fildes);
}

static long clk_tck = 0;

/*  etime    ..  Vinay    ... */
float etime_ (timear)
     float timear[2];
{
  struct tms buffer;
  time_t utime, stime;
  if (! clk_tck) clk_tck = sysconf(_SC_CLK_TCK);
  (void) times(&buffer);
  timear[0] = (float) buffer.tms_utime / (float)clk_tck;
  timear[1] = (float) buffer.tms_stime / (float)clk_tck;
  return (timear[0]+timear[1]);
}


outbuf_()
{
  /*
   *  routine to call setlinebuf and update disc after each write
   */

    setlinebuf(stdout);

}


/*   routine for  Linux rename  .. Vinay   */
int rename_(const char *from, const char *to)
{
 if (rename(from, to) !=0) perror(" rename error: ");
 
 return 0 ;
}
