#ifndef MT_RAND_H
#define MT_RAND_H

struct MTstate_t;
typedef struct MTstate_t MTstate;

MTstate* newMTstate(void);
void deleteMTstate(MTstate* mtst);

/* initializes mt[MT_N] with a seed */
void init_genrand(MTstate* mtst, unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(MTstate* mtst, unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(MTstate* mtst);

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(MTstate* mtst);

/* generates a random number on [0,1]-real-interval */
double genrand_real1(MTstate* mtst);

/* generates a random number on [0,1)-real-interval */
double genrand_real2(MTstate* mtst);

/* generates a random number on (0,1)-real-interval */
double genrand_real3(MTstate* mtst);

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(MTstate* mtst);
#endif