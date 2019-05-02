#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifndef INTEGER
#define INTEGER int32_t
#endif

void outbuf_(void) {
  (void) setvbuf(stdout, NULL, _IONBF, 0);
}

/* MinGW handles isatty(), but it always returns zero. */
INTEGER csatty_(INTEGER *fildes) {
  return isatty(0);
}

void getsys_(char *sysname, INTEGER *max, INTEGER *length) {
#ifdef __CYGWIN__
  const char *SYSNAME = "Cygwin";
#else
  const char *SYSNAME = "MinGW";
#endif
  *length = strlen(SYSNAME);
  memcpy(sysname,SYSNAME,*length);
}

