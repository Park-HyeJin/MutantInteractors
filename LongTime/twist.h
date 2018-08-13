#ifndef _twist_h
#define _twist_h

#ifndef _WIN32
#include <sys/time.h>
#else

#include <time.h>
#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

// Definition of a gettimeofday function

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	// Define a structure to receive the current Windows filetime
	FILETIME ft;

	// Initialize the present time to 0 and the timezone to UTC
	unsigned __int64 tmpres = 0;
	static int tzflag = 0;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		// The GetSystemTimeAsFileTime returns the number of 100 nanosecond 
		// intervals since Jan 1, 1601 in a structure. Copy the high bits to 
		// the 64 bit tmpres, shift it left by 32 then or in the low 32 bits.
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		// Convert to microseconds by dividing by 10
		tmpres /= 10;

		// The Unix epoch starts on Jan 1 1970.  Need to subtract the difference 
		// in seconds from Jan 1 1601.
		tmpres -= DELTA_EPOCH_IN_MICROSECS;

		// Finally change microseconds to seconds and place in the seconds value. 
		// The modulus picks up the microseconds.
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}

	if (NULL != tz)
	{
		if (!tzflag)
		{
			_tzset();
			tzflag++;
		}

		// Adjust for the timezone west of Greenwich
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}

	return 0;
}

#endif // _LINUX

#include <math.h>
#define INLINE inline

#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL
#define UMASK 0x80000000UL
#define LMASK 0x7fffffffUL
#define MIXBITS(u,v) (((u) & UMASK) | ((v) & LMASK))
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v) & 1UL ? MATRIX_A : 0UL))

static unsigned long state[MT_N];
static int twist_left = 1;
static int initf = 0;
static unsigned long *twist_next;

//------------------------------------------------------------------------

int gus () {struct timeval tp; gettimeofday(&tp, 0); return tp.tv_usec;}

//------------------------------------------------------------------------

double dtime () {
	struct timeval tp;

	gettimeofday(&tp, 0);
	return (tp.tv_sec + 1.0e-6 * tp.tv_usec);
}

//------------------------------------------------------------------------

void init_rnd (unsigned long s) {
	int j;

	state[0] = s & 0xffffffffUL;
	for (j = 1; j < MT_N; j++) {
		state[j] = (1812433253UL * (state[j - 1] ^ (state[j - 1] >> 30)) + j); 
		state[j] &= 0xffffffffUL;
	}
	twist_left = initf = 1;
}

//------------------------------------------------------------------------

INLINE static void twist_next_state (void) {
	unsigned long *p = state;
	int j;

	if (initf == 0) init_rnd(dtime());

	twist_left = MT_N;
	twist_next = state;

	for (j = MT_N - MT_M + 1; --j; p++) *p = p[MT_M] ^ TWIST(p[0], p[1]);

	for (j = MT_M; --j; p++) *p = p[MT_M - MT_N] ^ TWIST(p[0], p[1]);

	*p = p[MT_M - MT_N] ^ TWIST(p[0], state[0]);
}

//------------------------------------------------------------------------

INLINE unsigned long rnd32 () {
	unsigned long y;

	if (--twist_left == 0) twist_next_state();
	y = *twist_next++;

	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}

//------------------------------------------------------------------------

INLINE double drnd () {
	return rnd32() * (1.0 / 4294967296.0);
}

//------------------------------------------------------------------------

static double beta = 1.0;
void init_exprnd (unsigned long s, double bbeta) {
	init_rnd(s);
	beta = bbeta;
}

//------------------------------------------------------------------------

INLINE double exprnd () {
	double u;
	do u=drnd(); while(u==0.);
	return -log(u)/beta;
}

//------------------------------------------------------------------------

static double mu = 0.;
static double sig = 1.0;
static double storedval = 0.;
void init_normalrnd (unsigned long s, double mmu, double ssig) {
	init_rnd(s);
	mu = mmu;
	sig = ssig;
	storedval = 0.;
}

//------------------------------------------------------------------------

INLINE double normalrnd () {
	double v1,v2,rsq,fac;
	if(storedval == 0.) {
		do {
			v1 = 2.0*drnd() - 1.0;
			v2 = 2.0*drnd() - 1.0;
			rsq = v1*v1 + v2*v2;
		} while (rsq>=1.0 || rsq==0.0);

		fac = sqrt(-2.0*log(rsq)/rsq);
		storedval = v1*fac;
		return mu + sig*v2*fac;
	}
	else {
		fac = storedval;
		storedval = 0.;
		return mu + sig*fac;
	}
}

//------------------------------------------------------------------------

#endif
