#include <fortran.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int FSATTY(int *desc)
{
  return isatty(*desc);
}
