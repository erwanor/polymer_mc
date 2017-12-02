/*
  A C-program for MT19937, with initialization improved 2002/1/26.
  Coded by Takuji Nishimura and Makoto Matsumoto.

  Before using, initialize the state by using init_genrand(seed)
  or init_by_array(init_key, key_length).

  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  3. The names of its contributors may not be used to endorse or promote
  products derived from this software without specific prior written
  permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


  Any feedback is very welcome.
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

  The original version of http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c was modified by Takahiro Omi as
  - delete line 47 "#include<stdio.h>"
  - delete line 174 int main(void){...}
  - change N -> MT_N
  - change N -> MT_N
  - change the file name "mt19937ar.c" -> "MT.h"

  The original version of MT.h (www.sat.t.u-tokyo.ac.jp/~omi/code/MT.h) was modified by Koh M. Nakagawa.
  - split MT.h into MT.c + MT.h.
  - I added prototype declarations of functions in MT.h.
  - I added struct MTstate_t.
*/

#include "mt_rand.h"

#include "utils.h"

/* Period parameters */
#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

struct MTstate_t {
  unsigned long mt[MT_N]; /* the array for the state vector  */
  int mti; /* mti==MT_N+1 means mt[MT_N] is not initialized */
};

MTstate* newMTstate(void)
{
  MTstate* mtst = (MTstate*)xmalloc(sizeof(MTstate));
  mtst->mti = MT_N + 1;
  return mtst;
}

void deleteMTstate(MTstate* mtst)
{
  xfree(mtst);
}

void init_genrand(MTstate* mtst,
                  unsigned long s)
{
  mtst->mt[0]= s & 0xffffffffUL;
  for (mtst->mti=1; mtst->mti<MT_N; mtst->mti++) {
    mtst->mt[mtst->mti] =
      (1812433253UL * (mtst->mt[mtst->mti-1] ^ (mtst->mt[mtst->mti-1] >> 30)) + mtst->mti);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mtst->mt[mtst->mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
}

void init_by_array(MTstate* mtst,
                   unsigned long init_key[],
                   int key_length)
{
  int i, j, k;
  init_genrand(mtst, 19650218UL);
  i=1; j=0;
  k = (MT_N>key_length ? MT_N : key_length);
  for (; k; k--) {
    mtst->mt[i] = (mtst->mt[i] ^ ((mtst->mt[i-1] ^ (mtst->mt[i-1] >> 30)) * 1664525UL))
      + init_key[j] + j; /* non linear */
    mtst->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++; j++;
    if (i>=MT_N) { mtst->mt[0] = mtst->mt[MT_N-1]; i=1; }
    if (j>=key_length) j=0;
  }
  for (k=MT_N-1; k; k--) {
    mtst->mt[i] = (mtst->mt[i] ^ ((mtst->mt[i-1] ^ (mtst->mt[i-1] >> 30)) * 1566083941UL))
      - i; /* non linear */
    mtst->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
    i++;
    if (i>=MT_N) { mtst->mt[0] = mtst->mt[MT_N-1]; i=1; }
  }

  mtst->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

unsigned long genrand_int32(MTstate* mtst)
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mtst->mti >= MT_N) { /* generate N words at one time */
    int kk;

    if (mtst->mti == MT_N+1)   /* if init_genrand() has not been called, */
      init_genrand(mtst, 5489UL); /* a default initial seed is used */

    for (kk=0;kk<MT_N-MT_M;kk++) {
      y = (mtst->mt[kk]&UPPER_MASK)|(mtst->mt[kk+1]&LOWER_MASK);
      mtst->mt[kk] = mtst->mt[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<MT_N-1;kk++) {
      y = (mtst->mt[kk]&UPPER_MASK)|(mtst->mt[kk+1]&LOWER_MASK);
      mtst->mt[kk] = mtst->mt[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mtst->mt[MT_N-1]&UPPER_MASK)|(mtst->mt[0]&LOWER_MASK);
    mtst->mt[MT_N-1] = mtst->mt[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    mtst->mti = 0;
  }

  y = mtst->mt[mtst->mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

long genrand_int31(MTstate* mtst)
{
  return (long)(genrand_int32(mtst)>>1);
}

long genrand_int31_range(MTstate* mtst,
                         const long lo,
                         const long hi)
{
  return (long)(genrand_int32(mtst)>>1) % (hi - lo + 1) + lo;
}

double genrand_real1(MTstate* mtst)
{
  return genrand_int32(mtst)*(1.0/4294967295.0);
  /* divided by 2^32-1 */
}

double genrand_real2(MTstate* mtst)
{
  return genrand_int32(mtst)*(1.0/4294967296.0);
  /* divided by 2^32 */
}

double genrand_real3(MTstate* mtst)
{
  return (((double)genrand_int32(mtst)) + 0.5)*(1.0/4294967296.0);
  /* divided by 2^32 */
}

double genrand_res53(MTstate* mtst)
{
  unsigned long a = genrand_int32(mtst) >> 5, b = genrand_int32(mtst) >> 6;
  return(a*67108864.0 + b)*(1.0 / 9007199254740992.0);
}
