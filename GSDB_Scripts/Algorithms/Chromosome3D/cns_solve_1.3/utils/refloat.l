%{
void showfloat(const char *text);

#undef  YYLMAX
#define YYLMAX 65536
%}

dreall  ([0-9]+"."[0-9]*)
drealr  ([0-9]*"."[0-9]+)
ereall  ([0-9]+"."[0-9]*[eE][+-]?[0-9]+)
erealr  ([0-9]*"."[0-9]+[eE][+-]?[0-9]+)
float   ({dreall}|{drealr}|{ereall}|{erealr})

%%
[ ]*-?{float}   { showfloat((const char *) yytext); }
[A-Za-z]{float} { printf("%s", yytext); }
\015            ;
.               { putchar(yytext[0]); }
%%


#include <stdlib.h>
#include <string.h>


static const char  *progn = "refloat";


static void usage(void)
{
  fprintf(stderr, "usage: %s [-#] [file]\n", progn);
  exit(1);
}


static int  outprec = 2;


int main(int argc, char *argv[])
{
  int     iarg, n;
  char    *fnin;
  FILE    *fpin;
  char    xtrac;


  fnin = NULL;
  fpin = NULL;
  outprec = 2;

  for (iarg = 1; iarg < argc; iarg++)
  {
    if (argv[iarg][0] == '-')
    {
          n = sscanf(&argv[iarg][1], "%d %c", &outprec, &xtrac);
      if (n != 1 || outprec < 1) usage();
    }
    else
    {
      if (fnin) usage();
      fnin = argv[iarg];
    }
  }

  if (fnin == NULL)
    fpin = stdin;
  else
  {
        fpin = fopen(fnin, "r");
    if (fpin == NULL)
    {
      fprintf(stderr, "%s: Can't open %s\n", progn, fnin);
      exit(1);
    }
  }

  yyin = fpin;
  yylex();

  if (fpin != stdin)
    fclose(fpin);

  return 0;
}


void showfloat(const char *text)
{
  int     l;
  double  val;


  l = strlen(text);

  if (text[0] == ' ') {
    putchar(' ');
    l--;
    text++;
  }

  sscanf(text, "%lf", &val);
  printf("%*.*g", l, outprec, val);
}
