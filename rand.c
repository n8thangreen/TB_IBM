/*----------------------------------------------------------------------------*
RANDOM NUMBERS

Here is a short collection of functions to generate random numbers.  These
functions may be used when standard library functions for random numbers are not
available or when precise characteristics of the numbers are important.  With
these functions, the random sequences should not vary with the machine, with the
operating system, or with the compiler.  The collection comprises seven
functions:

    1. Rand             Generate random number on the unit interval.
    2. RandStart        Start the sequence at a specified point.
    3. RandStartArb     Start the sequence at an arbitrary point.
    4. RandStartNext    Start the sequence where it left off last.
    5. RandStopNext     Record the ending random seed for the next run.
    6. RandEndingSeed   Return the ending random seed.
    7. RandInteger      Generate a random integer (the inner function).
/*

/*----------------------------------------------------------------------------*
RETURN RANDOM NUMBER IN UNIT INTERVAL

The first function in the collection is the most common in scientific
programming.  It merely generates a series of pseudo-random numbers on the
interval 0<=X<1.  If it is the only function called, it will generate the same
sequence of numbers each time the program is run, thus providing repeatability
in numerical experiments.  A simple program that uses it to list ten random
numbers looks like this and prints the list that follows.

    main()
     {
     int n;

     for (n = 0; n < 10; n++)
       printf("%f\n", Rand());
     }

    0.211325
    0.544479
    0.220742
    0.111617
    0.893342
    0.290086
    0.212657
    0.105951
    0.686732
    0.749347

ENTRY: No significant conditions.

EXIT:  'Rand' contains the next random number in the sequence, uniformly
         distributed with '0<=Rand<1'. Only the first 32 bits are significant;
         the remainder are zero (if the machine's arithmetic hardware is
         accurate).
*/

static unsigned long seed;

double Rand()
 {
 unsigned long RandInteger();

 return((double)RandInteger() / 4294967296.);
 }

/*----------------------------------------------------------------------------*
INITIALIZE TO A KNOWN STARTING POINT

By default, the function 'Rand' starts at the beginning of a long deterministic
sequence, continues for 4,294,967,296 steps, then repeats. 'Rand' can be
started at any point within the sequence using the function 'RandStart'.  The
function is passed a ``seed," which is an integer in the range 0 to
4,294,967,295.  The default seed is zero.  With a seed of one, the following
program produces a different sequence.

    main()
     {
     int n;

     RandStart((unsigned long) 1);
     for (n = 0; n < 10; n++)
       printf("%f\n", Rand());
     }

    0.215868
    0.177158
    0.910775
    0.598857
    0.739466
    0.119943
    0.829061
    0.617727
    0.337687
    0.408679

ENTRY: 'k' contains a starting random number seed.

EXIT:  The random sequence is initialized to 'k'.
*/

unsigned long RandStart(unsigned long k)
 {
 seed = k;
 return(seed);
 }

/*----------------------------------------------------------------------------*
INITIALIZE TO AN ARBITRARY STARTING POINT

A series of arbitrary starting points is taken from the time of day. Function
'RandStartArb' sets the random number seed to a value that will change from
time to time, so the sequence will likely start at a different point each time
its program is run (if enough time elapses).  The function returns the
beginning seed, which can be used with function 'RandStart' to repeat the
sequence.

    main()
     {
     unsigned long start;
     int n;

     start = RandStartArb(0);
     for (n = 0; n < 10; n++)
       printf("%f\n", Rand());

     printf("\n");

     RandStart(start);
     for (n = 0; n < 10; n++)
       printf("%f\n", Rand());
     }

    0.302036
    0.718284
    0.199497
    0.233959
    0.267500
    0.025543
    0.451980
    0.901616
    0.583625
    0.535184

    0.302036
    0.718284
    0.199497
    0.233959
    0.267500
    0.025543
    0.451980
    0.901616
    0.583625
    0.535184

ENTRY: 'offset' contains a value to be added to the initial seed derived from
         the time of day. For a series of jobs starting in parallel at
         essentially the same time, this can be a job number to make each seed
         different.

EXIT:  The random number sequence has been set to an arbitrary point.
       'RandStartArb' returns the random number seed, which can be
         used to restart the sequence.
*/

static unsigned long reverse();
unsigned long time();

unsigned long RandStartArb(unsigned long offset)
 {
 static unsigned long base = 1234567;

 base = base*5 + 1;
 seed = base + offset + reverse(time((long*)0));
 return(seed);
 }

/*----------------------------------------------------------------------------*
SAVE AND RESTORE RANDOM NUMBER SEEDS

When a large simulation is divided into many parts, run in sequence, it may be
important that the random numbers remain independent from run to run.  The way
to assure that is to recording the ending random number seed at the end of one
run and reload it at the beginning of the next.  Routines in this section are
provided for that purpose.

Note that it is not safe to start any independent run except the first with a
call to 'RandStartArb'. The clock values obtained by that routine are not
random, but are related to one another, and will be identical if separate runs
occur in rapid succession.

LOAD BEGINNING SEED

This routine looks for a file containing the next seed and, if one is present
in the current directory, starts with that seed.  If none is present, then an
arbitrary seed is taken via 'RandStartArb'.

ENTRY: 's' points to a file name containing the starting seed, or is null if
         the default file name is to be used.

EXIT:  'RandStartNext' contains a result code.
         0 No starting seed was on file. An arbitrary seed has been used instead.
         1 A starting seed was retrieved from the file.
*/

#include <stdio.h>
static char rnd[] = "nextseed.rnd";
static char *file = rnd;

int RandStartNext(char *s)
 {
 unsigned long randseed; FILE *pf;

 if(s) file = s;                              //Retrieve any passed file name.

 pf = fopen(file, "r");                       //If no previous seed exists,
 if(pf==0)                                    //start with an arbitrary one.
 { RandStartArb(0); return(0); }

 fscanf(pf, "%lu", &randseed); fclose(pf);    //Otherwise resume the sequence
 RandStart(randseed); return(1);              //where it left off.
 }

/*
SAVE ENDING SEED

This routine saves the ending random number seed in a file that will be picked
upon a call to 'RandStartNext' the next time the program runs.

ENTRY: 's' points to a file name to receive the ending seed, or is null if
         the default file name is to be used.

EXIT:  The ending seed has been saved.
*/

RandStopNext(char *s)
 {
 unsigned long randseed, RandEndingSeed(); FILE *pf;

 if(s) file = s;                              //Retrieve any passed file name.

 unlink(file);                                //Delete any existing file.

 pf = fopen(file, "w");                       //Record the seed for the next
 if(pf)                                       //time the program runs.
 { fprintf(pf, "%lu\n", RandEndingSeed());
   fclose(pf); }
 }

/*----------------------------------------------------------------------------*
GENERATE INTEGER SEQUENCE

The actual random numbers are generated here as integers of 32 bits each,
created according to a certain nonlinear difference equation (a linear
congruential scheme).  Each successive number 'x(n+1)' in the series is
generated from the preceding number 'x(n)' by the relation
'x(n+1)=a x(n)+c (mod m)'.  The starting point 'x(0)' is called the seed.  The
multiplier 'a' is chosen so that '0<a-sqrt(m)<m', so that 'a==5 (mod 8)', and
so that the arrangement of its bits appears haphazard.  The number 'c' is an
odd integer as close to '(1-1/sqrt(3) m/2' as possible.  The number 'm' is the
length of the random number series to be generated, and is taken here to be
'2**32'. These choices are sufficient to assure that the numbers generated
will reasonably imitate a random series.  See Donald E. Knuth, The Art of
Computer Programming, Volume 2 (Addison Wesley, 1969) for a careful description
of the method.

The integer values can be used directly.  When they are, values less than
32 bits are usually needed, so only a portion of each generated number is used.
The portion taken should consist of a string of consecutive bits from the same
place in each number.  If a string of bits either from the high order or from
the central part of the number is taken, the resulting series will simulate
random numbers selected with replacement.  That is, duplicates and runs of the
same number can be expected.  If a string starting at the low order bit is
taken, the resulting series will simulate random numbers selected without
replacement.  In other words, each possible number will be generated once and
only once (then the series will repeat).

On many machines, unsigned long integers are precisely 32 bits in length.  If
they are longer than 32 bits, the constant 'OVERSIZE' should be defined so
that the operation '& 0xFFFFFFFF' is added to the code.  If they are shorter
than 32 bits, the function will not work.

Note:  In ANSI C, the file 'limits.h' has information on the sizes of integers.

ENTRY: Static 'seed' contains the previous number in the sequence.

EXIT:  'RandInteger' returns a new random number.
       'seed' contains the new random number for use in the next call.
*/

#define OVERSIZE

unsigned long RandInteger()
 {

#ifdef OVERSIZE
 return( seed = (seed*19513957 + 907633385) & 0xffffffff);
#else
 return(seed = seed*19513957 + 907633385);
#endif
 }

/*----------------------------------------------------------------------------*
RETURN ENDING RANDOM NUMBER SEED

The random number sequence can be restarted later if the ending random seed
is obtained from this routine and passed on a subsequent run to routine
'RandStart'.

ENTRY: No significant conditions.

EXIT:  'RandEndingSeed' contains the random number seed, which can be used
         to restart the sequence.
*/

unsigned long RandEndingSeed()
 {
 return(seed);
 }

/*----------------------------------------------------------------------------*
REVERSE ORDER OF BITS

To convert slightly different times to largely different values, this function
reverses the order of bits in an integer.  Values of at least 32 bits are
assumed.

ENTRY: 'k' contains a value to be reversed.

EXIT:  'reverse' contains the value on entry with the low 32 bits in the
         opposite order. That is, bit 0 becomes bit 31, bit 1 becomes bit 30,
         and so forth.
*/

static unsigned long reverse(unsigned long k)
 {
 int n; unsigned long h;

 for (h=n=0; n<32; n++)  { h = (h<<1) + (k&1); k >>= 1; }
 return(h);
 }


/* CLARENCE LEHMAN, JULY 1972 [Originally on the IBM System 360].

This code has been placed in the public domain by its author.  If you distribute
it, please include all the comments.  If you make any changes, please note them
below with your name or initials.

Modifications

 1. Converted to 8080 language, May 1976 [CLL].
 2. Converted to document format, December 1987 [CLL].
 3. Converted to C, April 1988 [CLL].
 4. Expanded for scientific use, June 1991 [CLL].
 5. Ending random seed, March 1994 [CLL].
 6. Description of the method modernized slightly, October 1994 [CLL].
 7. Saving and restoring of ending seed, October 2001 [CLL].
 8. Converted to ANSI C, January 2010 [CLL].
 9. Multiple arbitrary starting seeds, April 2011 [CLL].
*/

