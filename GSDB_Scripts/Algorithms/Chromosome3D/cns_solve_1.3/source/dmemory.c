#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef INTEGER
#define INTEGER int
#endif
#define CHARACTER unsigned char


#if   defined(CNS_ARCH_TYPE_CRAY)

#include <fortran.h>
#define FCHARPRO _fcd FDES
#define FCHARPTR (_fcdtocp(FDES))
#define FCHARLEN (_fcdlen(FDES))
#define cnsialloc_ CNSIALLOC
#define cnscalloc_ CNSCALLOC
#define cns0alloc_ CNS0ALLOC
#define cnsqalloc_ CNSQALLOC
#define cnsqptrsz_ CNSQPTRSZ

#elif defined(CNS_ARCH_TYPE_VMS)

#include <descrip.h>
#define FCHARPRO struct dsc$descriptor_s *FDES
#define FCHARPTR (FDES->dsc$a_pointer)
#define FCHARLEN (FDES->dsc$w_length)
#define cnsialloc_ cnsialloc
#define cnscalloc_ cnscalloc
#define cns0alloc_ cns0alloc
#define cnsqalloc_ cnsqalloc
#define cnsqptrsz_ cnsqptrsz

#elif defined(CNS_ARCH_TYPE_WIN32)

#include <malloc.h>
#define cnsialloc_ __stdcall CNSIALLOC
#define cnscalloc_ __stdcall CNSCALLOC
#define cns0alloc_ __stdcall CNS0ALLOC
#define cnsqalloc_ __stdcall CNSQALLOC
#define cnsqptrsz_ __stdcall CNSQPTRSZ

#elif defined(CNS_ARCH_TYPE_HP) \
   || defined(CNS_ARCH_TYPE_IBM_AIX)

#define cnsialloc_ cnsialloc
#define cnscalloc_ cnscalloc
#define cns0alloc_ cns0alloc
#define cnsqalloc_ cnsqalloc
#define cnsqptrsz_ cnsqptrsz

#endif


#ifndef FCHARPRO
#define FCHARPRO CHARACTER *refptr, unsigned INTEGER len
#define FCHARPTR (refptr)
#define FCHARLEN (len)
#endif


static long   dBytes = 0;
static long CurBytes = 0;
static long MaxBytes = 0;
static long   dOverh = 0;
static long CurOverh = 0;
static long MaxOverh = 0;


/* Normal layout of allocated memory:
   <overhead><newelem elements of size selem><velem><pointer...><overtail>
   ^         ^                               ^      ^
   ptrH      ptrU                            ptrV   ptrP

   For debugging, the number of copies of <pointer...> can be increased
   with the NPTR define.
 */
#ifndef NPTR
#define NPTR 1
#endif


static INTEGER cnsalloc(INTEGER oldelem, INTEGER newelem, int selem,
                        INTEGER *offs,
                        unsigned char *refptr)
{
  /*
    cnsalloc handles malloc(), realloc() and free() requests.

    newelem: number of new elements to allocate
    oldelem: number of elements allocated previously
    selem  : size of one element
    *offs  : ptrU = &refptr[(long)(*offs) * selem];
    refptr : address of FORTRAN reference object

    if      (oldelem < 0)  =>  malloc()
    else if (newelem < 0)  =>  free()
    else                   =>  realloc()

    if oldelem is not less than 0, *offs must be defined.

    Author: R.W. Grosse-Kunstleve
   */

  int            i;
  INTEGER        oelem,  nelem, melem, velem;
  long           obytes, nbytes;
  long           oOverHead, nOverHead;
  long           l0offs, loffs;
  unsigned char  *ptrH, *ptrU, *ptrV, *ptrP;

  const int  SzINTEGER = sizeof (INTEGER);
  const int  SzPOINTER = sizeof (unsigned char *);
  const int  madm = SzINTEGER + NPTR * SzPOINTER;


  if (oldelem < 0 && newelem < 0) return -__LINE__;

  oelem = oldelem; if (oelem < 0) oelem = 0;
  nelem = newelem; if (nelem < 0) nelem = 0;

  obytes = (long) oelem * selem;
  nbytes = (long) nelem * selem;

  dBytes = 0;
  dOverh = 0;

  if (oldelem < 0)
  {
    dBytes = nbytes;
    dOverh = madm;
        ptrH = malloc(nbytes + madm);
    if (ptrH == NULL) return 1;
                   CurBytes += nbytes;
    if (MaxBytes < CurBytes) MaxBytes = CurBytes;
                   CurOverh += madm;
    if (MaxOverh < CurOverh) MaxOverh = CurOverh;
  }
  else
  {
    ptrU = &refptr[(long)(*offs) * selem];
    ptrV = &ptrU[(long) oelem * selem];
    ptrP = &ptrV[SzINTEGER];
    (void) memcpy(&velem, ptrV, SzINTEGER);
    if (abs(velem) != oelem) return 2;
    for (i = 1; i < NPTR; i++)
      if (memcmp(ptrP, &ptrP[i * SzPOINTER], SzPOINTER) != 0)
        return 2;
    (void) memcpy(&ptrH, ptrP, SzPOINTER);
    if (ptrH == NULL) return 2;
        oOverHead = ptrU - ptrH;
    if (oOverHead < 0 || oOverHead > (selem - 1)) return 2;

                   dBytes = nbytes - obytes;
    if (velem < 0) dOverh = -(selem - 1);

    if (newelem < 0)
    {
      dOverh -= madm;
      free(ptrH);
      CurBytes += dBytes;
      CurOverh += dOverh;
      (*offs) = 0;
      return 0;
    }

        ptrH = realloc(ptrH, nbytes + madm);
    if (ptrH == NULL) return 1;
    CurBytes += dBytes;
    CurOverh += dOverh;
    if (MaxBytes < CurBytes) MaxBytes = CurBytes;
  }

  loffs = ptrH - refptr;

  velem = nelem;

  if (loffs % selem == 0)
    loffs /= selem;
  else {
    dOverh += (selem - 1);
        ptrH = realloc(ptrH, nbytes + madm + (selem - 1));
    if (ptrH == NULL) return 1;
                   CurOverh += (selem - 1);
    if (MaxOverh < CurOverh) MaxOverh = CurOverh;
    velem *= -1;
    l0offs = loffs = ptrH - refptr;
    loffs /= selem;
    if (l0offs > 0 && loffs * selem != l0offs) loffs++;
  }

  (*offs) = (INTEGER)  loffs;
  if ((long)(*offs) != loffs) return 3;

  ptrU = &refptr[(long)(*offs) * selem];
  ptrV = &ptrU[(long) nelem * selem];
  ptrP = &ptrV[SzINTEGER];

      nOverHead = ptrU - ptrH;
  if (nOverHead <           0) return -__LINE__;
  if (nOverHead > (selem - 1)) return -__LINE__;

  if (oldelem >= 0 && nOverHead != oOverHead) {
    melem = (oelem < nelem ? oelem : nelem);
    (void) memmove(ptrU, &ptrH[oOverHead], (long) melem * selem);
  }

  (void) memcpy(ptrV, &velem, SzINTEGER);
  for (i = 0; i < NPTR; i++)
    (void) memcpy(&ptrP[i * SzPOINTER], &ptrH, SzPOINTER);

  return 0;
}


INTEGER cnsialloc_(INTEGER *oldelem, INTEGER *newelem, INTEGER *offs,
                   INTEGER *refptr)
{
  return cnsalloc(*oldelem, *newelem, sizeof (INTEGER), offs,
                  (unsigned char *) refptr);
}


INTEGER cnscalloc_(INTEGER *oldelem, INTEGER *newelem, INTEGER *offs,
                   INTEGER *strlen, FCHARPRO)
{
  if (FCHARLEN == 0) return -__LINE__;
  return cnsalloc(*oldelem, *newelem, (*strlen) * sizeof (CHARACTER), offs,
                  (unsigned char *) FCHARPTR);
}


void cns0alloc_(void)
{
    dBytes = 0;
  CurBytes = 0;
  MaxBytes = 0;
    dOverh = 0;
  CurOverh = 0;
  MaxOverh = 0;
}


INTEGER cnsqalloc_(INTEGER *VarID, INTEGER *Digits, INTEGER *mDigits)
{
  long  var;
  int   Sign, nDigits;

  switch (*VarID) {
    case  1: var =   dBytes; break;
    case  2: var = CurBytes; break;
    case  3: var = MaxBytes; break;
    case -1: var =   dOverh; break;
    case -2: var = CurOverh; break;
    case -3: var = MaxOverh; break;
    default: return 0;
  }

  Sign = 1; if (var < 0) { Sign = -1; var = -var; }

  for (nDigits = 0;;) {
    if (nDigits >= *mDigits) return 0;
    Digits[nDigits++] = var % 10;
    var /= 10;
    if (var == 0) break;
  }

  return (INTEGER)(Sign * nDigits);
}


void cnsqptrsz_(INTEGER *size)
{
  *size = (INTEGER)(sizeof (unsigned char *));
}
