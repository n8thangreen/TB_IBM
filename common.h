/* DATA STRUCTURES COMMON AMONG THE MODULES

This module defines the main data structure in the program, the list of
individuals and their characteristics, also the main space-consuming entity in
the program.

Time is encoded in 8-byte double-precision floating point values ('dec'), which
is simple but a little extravagant. At the cost of a clarity time could be
encoded in 4-byte unsigned integers instead, with a single floating point number
defining the base year and time.
*/

#ifndef TYPEDEF
#define TYPEDEF
typedef double        dec;
typedef short         ints;
typedef unsigned char charu;
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//#define indiv 20000000    //Maximum population size.
#define indiv 75000000  //Maximum population size.
#define INDIV indiv       //Maximum population size.

// FUTURE TIMES:
#define iBirth        0   //Time of initiation of this record
#define iExit         1   //Time for exit from this state
#define iDeath        2   //Time for closure of this record
#define iDisease      3   //Time of progression to disease
#define iTransm       4   //Time to transmit infection to another
#define iMutate       5   //Time of strain type mutation
#define iEmigrate     6   //Time of emigration
#define iRep          7   //Time to report disease case
#define tBirth      t[0]  //Time of initiation of this record
#define tExit       t[1]  //Time for exit from this state
#define tDeath      t[2]  //Time for closure of this record
#define tDisease    t[3]  //Time of progression to disease
#define tTransm     t[4]  //Time to transmit infection to another
#define tMutate     t[5]  //Time of strain type mutation
#define tEmigrate   t[6]  //Time of emigration
#define tRep        t[7]  //Time to report disease case

//#define tImm      t[0]  //Time of immigration to UK
//#define tInfected t[0]  //Time individual was most recently infected
//#define strain    t[0]  //Strain type identification number

struct Indiv           //STRUCTURE OF EACH RECORD                   BYTES
{
  dec  t[8];           //Separate times for individual               64
  charu pending;       //Number of pending event                      1
  charu state;         //Number of present state                      1
  char sex;            //Sex of this individual (0=female, 1=male)    1
  char rob;            //Region of birth (0=Foreign-born, 1=UK-born)  1
  //char inf;          //Region infection acquired (0=abroad, 1=UK)   1
  char ssa;            //0=UK & non-UK other(HIV-), 1=SSA (HIV-) 2=SSA (HIV+) 1
};                     //                                            69 *

extern struct Indiv *A;   //List of individuals.

// STATES:
#define qU         1   //Uninfected
#define qV         2   //Immune
#define qI1        3   //Recent infection
#define qI2        4   //Remote infection
#define qI3        5   //Reinfection
#define qD1        6   //Primary disease, pulmonary
#define qD2        7   //Reactivation disease, pulmonary
#define qD3        8   //Reinfection disease, pulmonary
#define qD4        9   //Primary non-pulmonary disease
#define qD5       10   //Reactivation non-pulmonary disease
#define qD6       11   //Reinfection non-pulmonary disease
#define q0        qU   //Lowest numbered state
#define q1       qD6   //Highest numbered state

// TRANSITIONS:
#define pVaccin    1   //Pending vaccination
#define pTransm    2   //Pending transmission of an infection
#define pRemote    3   //Pending transition to remote infection
#define pDisease   4   //Pending progression to disease
#define pDeath     5   //Pending death
#define pMutate    6   //Pending strain type mutation
#define pEmigrate  7   //Pending emigration
#define pBirth     8   //Pending birth
#define pImmig     9   //Pending immigration
#define pRep      10   //Pending reporting of case

// PROTOTYPES:
dec Rand();                                    //Random number generator
dec Uniform(dec, dec);                         //Uniform random number
dec Expon(dec);                                //Poisson events in time
dec Cauchy(dec, dec);                          //Cauchy distribution
dec Gauss(dec, dec);                           //Gaussian distribution
dec LogNormal(dec, dec);                       //Lognormal distribution
dec RecovDsn(int, dec, dec);                   //Recovery time from disease
dec LifeDsn (int, dec, dec);                   //Lifespan of individual
dec Gompertz(int, dec, dec);                   //Lifespan of individual
char *Tval(dec);                               //Time conversion
unsigned long RandStartArb();                  //Random number initializers
unsigned long RandStart(unsigned long);
dec Val(int, dec, dec[], dec[], int, int);
dec RandF(dec[], dec[], int, dec);
int Loc(dec[], int, int, dec);
dec Tdis(int,int,int,int,dec);
dec GetAge(int, int, int);
dec EmDsn(int, int, dec, dec);
int Error (dec);
int Error1(dec, char*,dec);
int Error2(dec, char*,dec, char*,dec);
int Error3(dec, char*,dec, char*,dec, char*,dec);
int StrainNum(int);

// LOCAL FUNCTIONS:
#define min(a,b) ((a)<(b)? (a):(b))            //Minimum
#define max(a,b) ((a)>(b)? (a):(b))            //Maximum
#define abs(a)   ((a)>=0?  (a):-(a))           //Absolute value
#define round(a) ((a)>=0?  (a)+0.5: (a)-0.5)   //Rounding to integer

// Clarence Lehman, August 2009, modified by Adrienne Keen, 2010

