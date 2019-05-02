#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <time.h>

int outbuf_(void)
{
  (void) setvbuf(stdout, NULL, _IONBF, 0);
}

int fsatty_(int *dummy)
{
  (void) clock(); /* also initialize clock function */
  if(_isatty(_fileno(stdin))) return 1;
  return 0;
}

float cnsclock_(int *dummy)
{
  return (float) clock() / (float) CLOCKS_PER_SEC;
}
