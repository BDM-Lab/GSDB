#include <stdio.h>
#include <stdlib.h>
#include <io.h>

void __stdcall OUTBUF(void)
{
  (void) setvbuf(stdout, NULL, _IONBF, 0);
}

int __stdcall FSATTY(int *dummy)
{
  if(_isatty(_fileno(stdin))) return 1;
  return 0;
}
