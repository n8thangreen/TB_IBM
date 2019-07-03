/*============================================================================*
INDIVIDUAL-BASED MODEL FOR TUBERCULOSIS DYNAMICS IN THE UK

This individual-based model (IBM) is used to simulate tuberculosis
dynamics in the UK. Although not currently implemented, the model is designed
to follow individual strain types for every infection to reproduce genotype
clustering patterns seen in disease cases. The first application of this model
is fit to England and Wales case notifications, without considering genetic
typing data, but modelling a large population, about 55 million individuals.
Subsequent applications use data from Scotland and the West Midlands, both areas
with around five million people. One major characteristic of individuals is
their region of birth, UK or Non-UK. In the 'SSAV' version of the model, the
Non-UK-born region of birth is divided into Sub-Saharan African born (SSA-born)
and other non-UK born (ONUK-born).

1. BACKGROUND ALGORITHMS. The core IBM algorithms, including the event
scheduler, event queue, and event dispatcher were written by Clarence Lehman
(CL) in 2009. The method was developed for simulation of HIV dynamics in the US.
Beginning in January 2010, and with intial help from CL, Adrienne Keen (AK) has
adapted this IBM skeleton for modelling TB dynamics in the UK.

2. SUMMARY OF METHOD. This is an event-based simulation; continuous time is
simulated directly. Thus only one event occurs in any second. The time variable
't' is time in years, with arbitrarily high resolution down to small fractions
of a instant and all complexities and inaccuracies associated with multiple
events during a finite time step vanish -- such as undershooting zero when the
sum of the rates times the width of the time step exceeds unity. Continuous time
also allows all activities to occur in a single data array, rather than having
to swap old and new arrays at each time step. Events are assumed to following
probability distributions that vary through time and space.

States of the system never change spontaneously -- all changes are induced by
some other event in the system and usually scheduled in advance. The scheduled
times are determined stochastically from functions whose characteristics may
depend on the state of the individual and the environment at the time. An
individual's age, sex, infection history, or any other considerations can be
incorporated into the functions.

For example, death is scheduled in advance at the time of birth, with the time
chosen randomly from a life-span distribution for babies born in the simulated
year. But the scheduled time of death is not immutable, nor are any other
scheduled times in the system. If the individual is infected, the scheduled time
of death may be cancelled and a time to disease development may be scheduled
instead. At no time does the program visit an individual when it does not need
to, and therein lies its speed.

Many future events may apply to each individual and are saved for that
individual, but only the earliest among each individual's events enters a global
"list of future events."

3. MAIN DATA STRUCTURES. Each individual is assigned a number 1 through 'n' and
recorded in a linear array 'A' of structures 'Indiv'. Each element of 'A'
defines the state of the correspondng individual.

ILLUSTRATIVE DIAGRAM: struct A[indiv+3]
This is the main array of individuals. In this example there are 3 UK-born
and 3 Non-UK-born individuals, with a total maximum population size ('indiv')
of 14. Note, in SSAV version of the model, SSAs are stored as Non-UK born and
it is not possible to tell from their ID number alone whether they are SSA or
other non-UK born.

--What is in the array------Generic array index--------------------------------
0  [Reserved]          (null pointer for list)
1  (Non-UK-born)
2  (Non-UK-born)
3  (Non-UK-born)             immid-1
4  [Empty]                   immid
5  [Empty]
6  [Empty]
7  [Empty]                   maximm
------------------------------------------(Imaginary line b/n non-UK/UK born)
8  (UK-born)                 maximm+1
9  (UK-born)
10 (UK-born)                 ukbid-1
11 [Empty]                   ukbid
12 [Empty]
13 [Empty]
14 [Empty]                   indiv
15 Ext event, birth          indiv+1 or BIRTH
16 Ext event, immigration    indiv+2 or IMM
-----------------------------------------------------------------------------

4. MODEL FITTING
The tb32 version of model accepts four variable disease parameters. These are:
'df', 'd1uk20[M]', 'd2uk20[M]', and 'd3uk20[M]'. See the Data() function for
more information on these. Briefly, 'df' is the factor by which UK-born
disease risks are multiplied to obtain Non-UK born disease risks. The other
three parameters are UK-born disease risks for Primary, Reactivation, and
Reinfection Disease respectively. Risks are for those aged 20 years and over.
Disease risk are fixed for children under ten, allowing for fewer variable
parameters. This program is called as a function from 'fit5i.c' (not as an
independent executable, as done previously) so that it is compatible with
parallele runs using MPI commands.

Example new program call for executable (slightly different syntax is used
for the call inside 'fit5i.c'):
tb32 df=2.5 d1uk20=0.15 d2uk20=0.0004 d3uk20=0.09

5. OTHER

To run TB IBM stand-alone, as opposed to inside the fitting routine ('fit5i'),
simply comment out the '#define main mainiac' and everything else is
automatically handled.

Note that the following program substitutes the term "dec" for the C term
"double". It is short for "decimal", in contrast and parallel with "integer",
saving valuable coding columns at the left of the line and helping data names
line up.

The sequence of random numbers is specified on the command line with phrases
like 'randseq=0', 'randseq=1', 'randseq=-6', 'randseq=239702397623', and the
like. Fixed sequences that are the same each time the program runs occur when
'randseq' is 0 or any positive integer. Each integer gives a different sequence
of random numbers. (Actually just a different starting point in a single long
sequence of random numbers.) Positive integers, or zero, are typically used in
testing, because program results are repeatable.

Arbitrary sequences that are different, with high probability, each time the
program runs occur when 'randseq' is negative. The date and time, measured to
the nearest second, selects the starting sequence, then the negative value
modifies that sequence. Thus if several instances of the program were started on
separate processors at the same time, the first with 'randseq=-1', the second
with 'randseq=-2', and so forth, each instance of the program is guaranteed a
different random number sequence. Unlike the case with non-negative integers,
however, the sequence be different each time the program runs, with very high
probability.

The actual starting seed, incorporating the time of day if requested by a
negative value of 'randseq', is stored in 'rand0' and reported at the end of the
run. That allows a run to be repeated exactly even if it was started with an
arbitrary sequence.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"
#include "fileio.h"

#define PN (q1+1)              //Number of elements in array 'N'.
#define T0 1981                //Start time of model, years.
#define T1 2010                //End time of model, years. The simulation
                               //ends -before- reaching this year.
#define SSAV   1               //Switch model version depending on existence
                               //of separate SubSaharan African group,
                               //0=non-SSA, 1=SSA.
#define SUPER  1               //Notes whether model is run on supercomputer,
                               //0=no, 1=yes (this changes population sizes).
#define DPARAM 1               //Allows model to accept disease progression
                               //parameters (4 in this version), 0=no, 1=yes.
#define BIRTH (indiv+1)        //Index used for scheduling births.
#define IMM   (indiv+2)        //Index for scheduling arrival of immigrants.
#define RT    (T1-T0)          //Running time of model, calendar years.
#define NUK   0                //Array index for non-UK born.
#define UK    1                //Array index for UK-born.
#define HIV   2                //Array index for HIV+.
#define SSA   2                //Array index for SSA-born.
#define M     0                //Array index for males.
#define F     1                //Array index for females.
#define E     0.0000000001;    //Small number added to some event times to
                               //ensure they happen in the future.
#define is1 (iprev1*iimm*0.70) //# strains to draw from for immigrants.
#define is2 (iprev2*iukb*0.70) //Initial # strains to draw from for UK-born.
#define AC 122                 //Age classes for mortality data.
#define LAT 5                  //Years to Remote from recent (re)infection.
#define BY (2010-1870+1)       //Number of birth cohorts for mortality data.

dec N[PN];                     //Current number in each disease state.
dec N2[4][2][3][RT];           //Population sizes in the model at end of year by
                               //age, sex, and rob.
dec N3[4][2][3][RT];           //Population sizes observed in England and Wales,
                               //which are compared to model population sizes
                               //and used to correct case numbers produced by
                               //the model.
dec age1[2],age2[2],agec[2];   //Accumulators for 1st and 2nd moments of age.

dec repc[4][2][3][2][RT];      //Array which holds reported cases for E&W
                               //version of model (simpler version) by age
                               //category, sex, rob, disease site & year.

int deaths;                    //Current number of deaths.
int events;                    //Current number of events dispatched.
int immid;                     //Next available ID number for immigrants.
int ukbid;                     //Next available ID number for UK-born.
int stid;                      //Next available ID for new strain types.

extern dec t;                  //Current time (Managed by 'EventSchedule').
dec pt;                        //Time of previous report.
dec t0 = T0;                   //Beginning time of simulation.
dec t1 = T1;                   //End time of simulation.
int lup;                       //Year of last update to birth & immigration
                               //rates, which are sensitive to calendar year.

unsigned long startsec;        //Starting clock time, seconds of Unix.
unsigned long rand0;           //Starting random number seed.

struct Indiv *A;               //State of each individual, including their
                               //characterisitics, saved event times, etc.

FILE *cc;                      //Create files pointers for output of
FILE *rc;                      //cumulative cases (all) and reported cases.

/*----------------------------------------------------------------------------*
PARAMETERS AND CONTROL VARIABLES
*/

/* Population initialization */
int iimm = 1000000;            //Approximate initial number of immigrants (only
                               //used to set number of intial strain types).
int iukb = 10000000;           //Approximate initial number of UK-born (only
                               //used to set number of initial strain types).
int maximm;                    //Maximum immigrants in pop'n at any time.
dec inf1981[121][3][2][9];     //Cumulative probabilities of the 9 disease
                               //states for pop. initialization (by a,s,rob).
dec n1981[121][2][2];          //Numbers in each age/sex/rob category at
                               //population initialization, 1981.
dec ssa1981[121][2];           //Proportion SSA by age/sex category.
dec iprev1 = 0.15;             //Initial prevalence of infection (immigrants).
dec iprev2 = 0.01;             //" " (UK-born), both used for setting
                               //number of initial strain types.

/* Infection transmission */
dec c[2][2];                   //Effective contacts per year per pulmonary
                               //case (smear+) by sex and region of birth.
dec pcc = 0.50;                //Probability effective contact is close contact
                               //(drawn from within own region of birth,
                               //UK or non-UK).
dec  s2[] = {1,1};             //Relative susceptibility to reinfection (s).
dec smear[121];                //Proportion smear positive by age.

/* Vaccination */
dec  v1[] = {0.71,0.71};       //Efficacy of vaccine (rob).
dec  v2[] = {0.80,0.80};       //Portion vaccinated at designated age (rob).
dec  v3[] = {13,13};           //Average age of vaccination (rob).

/* Disease progression */
dec d1[2][3][121];             //Proportion Recently Infected who progress to
                               //disease over first 5 years of infection
                               //(by sex,rob,age)
dec d3[2][3][121];             //Proportion Reinfected who progress to disease
                               //over first 5 yrs of reinfection (by a,s,r).
dec drr[6];                    //Cumulative, relative risk of disease
                               //progression by year since infection, used
                               //with d1 & d3 (for first 5 years of infection).
dec B1[6];                     //Array of values for finding cumulative risk of
                               //in first five years of infection/reinfection,
                               //used with drr.
dec d2[2][3][AC+2];            //Proportion of those Remotely Infected who
                               //progress to disease, cumulative dsn
                               //by sex,rob (where r =0,1,or2. 2=HIV+ & SSA in
                               //SSA version of model) and age.
dec A2[AC+2];                  //Array of values for finding random time to
                               //disease for Remote Infection.
dec ehiv;                      //Factor by which non-UK born disease risks are
                               //multiplied for HIV+ SSA individuals.
dec df;                        //Factor by which UK-born disease progression
                               //rates are multiplied to get immigrant rates.
dec d1uk10[2];                 //Rates of disease progression for primary (1),
dec d1uk20[2];                 //Reactivation (2), and Reinfection (3) disease
dec d2uk10[2];                 //for those aged 0-10 (10) and 20+ (20) by sex.
dec d2uk20[2];                 //These help construct 'd1', 'd2', and 'd3'.
dec d3uk10[2];
dec d3uk20[2];
dec sdf1[2];                   //Sex factors, ratios of female:male disease
dec sdf2[2];                   //progression risks/rates (by age 0-10,20+).
dec sdf3[2];
dec presp = 0.77;              //Proportion of all TB which is respiratory --
                               //for correcting disease risks in Vynn & Fine to
                               //pulmonary disease risks.
dec p1[121][2][2];             //Portion pulm -primary disease (a,s,rob)
dec p2[121][2][2];             //Portion pulm -reactivation disease (a,s,rob)
dec p3[121][2][2];             //Portion pulm -reinfection disease (a,s,rob)
dec duk1p[2][2];               //Intermediate parameters for pulmonary-only
dec duk2p[2][2];               //rates of disease progression. Used for
dec duk3p[2][2];               //incorporating estimated rates from Vynnycky &
                               //Fine into rates for this model, which are
                               //combined pulmonary/non-pulmonary. Indexed by
                               //age (0-10yrs, 20+yrs) and sex.
dec d1p[121][2][2];            //Intermediate param's for pulmonary-only rates
dec d2p[121][2][2];            //of disease progression. Used to get overall
dec d3p[121][2][2];            //rates using those estimated by Vynn. & Fine
                               //(which are pulmonary rates). Indexed by age,
                               //sex, and rob. Precursors to d1, d2, d3.

/* Disease recovery, indexed by sex */
dec  r3[] = {0.5,0.5};         //Primary disease recovery rate
dec  r4[] = {0.5,0.5};         //Reactivation disease recovery rate
dec  r5[] = {0.5,0.5};         //Reinfection disease recovery rate
dec  r6[] = {0.5,0.5};         //Primary non-pulmonary disease recovery rate
dec  r7[] = {0.5,0.5};         //Reactivation non-pulm. disease recovery rate
dec  r8[] = {0.5,0.5};         //Reinfection non-pulm. disease recovery rate

/* Mortality */
dec A1[AC];                    //Holds ages 0-121 which correspond to the
                               //cumulative probabilities in M1.
dec M1[BY][2][AC];             //Cumulative probabilities of death by a given
                               //birth cohort, sex and age.

dec cft[121][2][RT];           //Case fatality rate due to TB (a,type dis,y)
/*These are old mortality rates used to generate lifetimes with exponential
distribution, left to keep compatibility with testing version of model */
dec  m1 [2][RT];               //Mortality of uninf/vacc/inf ind's (sex, y)
dec  m6 [2][RT];               //Mortality of primary disease (sex,y)
dec  m7 [2][RT];               //Mortality of reactivation disease (sex,y)
dec  m8 [2][RT];               //Mortality of reinfection disease (sex,y)
dec  m9 [2][RT];               //Mortality of primary non-pulm. disease (sex,y)
dec  m10[2][RT];               //Mortality of reactivated non-pulm. disease (s,y)
dec  m11[2][RT];               //Mortality of reinfection non-pulm. disease (s,y)

/*Birth and migration */
dec bcy[RT];                   //Births by calendar year.
dec pmale[RT];                 //Portion of newborns who are male by year.
dec immig[RT];                 //Total (uk+non-uk-born) immigrants by year.
dec pimm[RT];                  //Proportion of immigrants non-UK born by year.
dec ssaim[RT];                 //Proportion of non-UK born immigrants born in
                               //SubSaharan Africa by year.
dec hivp[2][RT];               //HIV Prevalence in SubSaharan African born
                               //immigrants by sex,year.
dec immsex[RT][3];             //Proportion immigrants who are male by yr and
                               //rob:0=non-UK,1=UK, 2=non-UK SSA. Note there
                               //are different input files for SSA and non-
                               //SSA version of the model.
dec immage[RT][2][3][7];       //Cumulative proportion of immigrants (by yr,
                               //sex,rob) in age classes, for use with RandF.
//Ay[7];                       //State variables to accompany 'immage'.
dec immageX[RT][2][3][6];      //Probabilities of 6 age classes (by yr,sex,
                               //rob) from ONS inflow data, ---as in 'immsex'
                               //note there are 2 versions of the input file
                               //for this array. Precursor to 'immage'.
dec infimm[121][3][RT][9];     //Cumulative probabilities immigrants enter
                               //disease states, by age,rob and year.
dec Ax[9];                     //State variables which accompany 'infimm'.
dec ypb, ypi;                  //Years per birth, years per immigrant.
dec em[2][3];                  //Annual emigration rate by sex, rob.

/* Assorted. */
/* The following two recovery rates will not be used in version of model that
defines Remote Infection as 'LAT' years after most recent infection */
dec  r1[] = {0.20,0.20};       //Rate Recent Infection moves to remote (s)
dec  r2[] = {0.20,0.20};       //Rate Reinfection moves to remote (s)

dec md = 0.01;                 //Mutations per year per strain (diseased)
dec mi = 0.1*.01;              //Mutations per year per strain (infected)

dec proprep = 0.75;            //Proportion of cases reported.
//-dec propcp  = 0.70;         //Proportion of reported cases culture +.

dec relativetime = 0;          //Set for relative time reporting.
dec randseq = 0;               //Random number sequence (set with 'randseq=N').
dec tgap    = 0.5;             //Time between reports, years.
dec kernel  = 0;               //Contagion kernel, 0=Panmictic, 1=Cauchy.
dec sigma   = 1;               //Width of contagion kernel, where applicable.
/*
static char *fn[] =            //Output files.
{ "allcases.txt",              // 0
  "repcases.txt",              // 1
  "cases.txt",                 // 2
  "pop.txt",                   // 3
   0 };
*/
struct IO fmt[] =              //Format statements for input/output.
{ /*00*/ { (dec*)bcy,     {-'i',RT} },
  /*01*/ { (dec*)immig,   {-'i',RT} },
  /*02*/ { (dec*)pimm,    {-'i',RT} },
  /*03*/ { (dec*)ssaim,   {-'i',RT} },
  /*04*/ { (dec*)pmale,   {-'i',RT} },
  /*05*/ { (dec*)hivp,    {-'s',2,-'Y',RT}, {-'y',-'S'} },
  /*06*/ { (dec*)infimm,  {-'a',121,-'r',3,-'y',RT,-'q',9}, {-'R',0,SSAV+1,-'Y',-'Q',-'A'} },
  /*--*/ { (dec*)inf1981, {-'a',121,-'r',2,-'q',9},  {-'a',120,0,-'r',UK,UK,  -'q',1,7} },
  /*--*/ { (dec*)inf1981, {-'a',121,-'r',2,-'q',9},  {-'a',120,0,-'r',NUK,NUK,-'q',1,7} },
  /*09*/ { (dec*)ssa1981, {-'a',121,-'s',2}, {-'s',-'A'} },
  /*10*/ { (dec*)n1981,   {-'a',121,-'s',2,-'r',2},  {-'s',-'a',-'R',1,0,1} },
  /*11*/ { (dec*)immsex,  {-'i',RT,-'r',3}, {-'i',-'R',0,SSAV+1} },
  /*12*/ { (dec*)immageX, {-'y',RT,-'s',2,-'r',3,-'a',6}, {-'y',-'r',0,SSAV+1,-'s',-'a'} },
  /*13*/ { (dec*)immage,  {-'i',RT,-'s',2,-'r',3,-'a',7}, {-'i',-'R',0,SSAV+1,-'A',-'s'} },
///*14*/ { (dec*)M1,      {-'b',BY, -'s',2,-'a',AC}, {-'s',-'a',-'B'} },
  /*14*/ { (dec*)M1,      {-'b',BY, -'s',2,-'a',AC}, {-'s',-'b',-'A'} },
  /*15*/ { (dec*)cft,     {-'a',121,-'d',2,-'y',RT} }, //{-'a',-'d',-'Y'} },
  /*16*/ { (dec*)d1,      {-'s',2,-'r',3,-'a',121}, {-'s',-'r',0,1,-'A'} },
  /*17*/ { (dec*)d2,      {-'s',2,-'r',3,-'a',124}, {-'s',-'r',0,2,-'A'} },
  /*18*/ { (dec*)d3,      {-'s',2,-'r',3,-'a',121}, {-'s',-'r',0,1,-'A'} },
  /*19*/ { (dec*)inf1981, {-'a',121,-'s',2,-'r',3,-'q',9}, {-'r',-'s',-'A',120,0,-'Q',1,8} },
  /*20*/ { (dec*)smear,   {-'a',121} },
  /*21*/ { (dec*)N3,      {-'a',4,  -'s',2,-'r',3,-'y',RT} },
  /*22*/ { (dec*)repc,    {-'a',4, -'s',2, -'r',3, -'d',2, -'i', RT} },
         { }
};


/*----------------------------------------------------------------------------*
MAIN INITIALIZATION

This routine should be called each time the program is reused, to clear static
variables for the next run. The function was added when 'tb30i.c' was
made into a function of the fitting routine, to implement parallel,
replicate runs of the TB program. This would not be necessary if the program
were called as independent executable, as before.
*/

MainInit()
{ int i,j,k,l,m;

  for(i=0; i<PN; i++) N[i] = 0;
  for(i=0; i<2;  i++) age1[i] = age2[i] = agec[i] = 0;

  for(i=0; i<4;  i++)
  for(j=0; j<2;  j++)
  for(k=0; k<3;  k++)
  for(l=0; l<2;  l++)
  for(m=0; m<RT; m++)
  { repc[i][j][k][l][m] = 0;
    N2[i][j][k][m]      = 0;
    N3[i][j][k][m]      = 0;}

  deaths = events = immid = ukbid = stid = 0;
  t = pt = 0;
}



/*----------------------------------------------------------------------------*
MAIN PROGRAM

Due to an MPI bug not allowing 'popen' etc to work, the TB program is
defined as a function which returns an array of output (rather than stand-alone
executable) for use with the fitting routine. 'define' statement is used to
control whether TB program is a stand-alone executable or function within
the fitting routine.

*/

//#define main dec *mainiac                    //Make this main not the real main.
#ifdef  main
static int fit5i = 1;                        //Flag set when linked with fitting.
#else
static int fit5i = 0;                        //Flag set when not linked.
#endif
static int fitm  = 0;                        //Flag set when fitting to rates,
                                             //(0=Numbers, 1=Rates).

dec  out[1000]; int outi;                    //Main output array (all).
dec outn[1000]; int outni;                   //Main output array (case numbers).

main(int argc, char *argv[])
{ int i, j, k, l, n, sid;

  startsec = time(NULL);                     //Retrieve the wall-clock time.

  if(fit5i==0) ErrorInit();                  //Trap system failures.
  MainInit();                                //Start the main program.
  EventInit();                               //Start the event queue.
  FinalInit();                               //Start the final reports.
  ReportInit();                              //Start the output reports.

  A = (struct Indiv *)                       //Allocate array of individuals.
      calloc(indiv+3, sizeof(struct Indiv)); //(Not static because of gcc bug
  if(A==0) Error(911.);                      //restricting such arrays to 2GB.)
/*
* cc = fopen(fn[0], "w");                    //Open output files and write
* rc = fopen(fn[1], "w");                    //file headers to them.
* if(cc==0) Error1(510., fn[0],0);
* if(rc==0) Error1(510., fn[1],0);
*/
  Init();                                    //Initialize output files.

  c[M][NUK]=c[F][NUK]=6.0;                   //Fill c, df, ehiv with values for
  c[M][UK] =c[F][UK] =6.0;                   //testing.
  df = 2.0; ehiv = 7.0;

  if(SUPER) maximm = 10000000;               //Adjust 'maximm' depending on
  else      maximm =  5000000;               //whether running on supercomp.

  Data();                                    //Read in appropriate data files
                                             //and store to arrays.
  gparam(argc, argv);                        //Collect parameters for this run
                                             //which have been specified on
                                             //command line.
  Param();                                   //Update variables/distributions
                                             //affected by parameters which
                                             //can change with each model run.

  if(bcy[0]<=0.0001)                         //Calculate years per birth and
  { ypb = RT*100;                            //per immigrant at t=t0 for
    printf("Births are zero!\n"); }          //scheduling them regularly. If
  else ypb = 1./bcy[0];                      //births (immigrants) per year
  if(immig[0]<=0.0001)                       //are zero, make interval very
  { ypi = RT*100;                            //large so they do not ever
    printf("Immigrants are zero!\n"); }      //happen.
  else ypi = 1./immig[0];

                                             //Update time of last update for
  lup = t0;                                  //parameters sensitive to
                                             //calendar year.

  rand0 = abs(randseq);                      //Start the random number sequence
  if(randseq>=0)  RandStart(rand0);          //from a specified or an arbitrary
  else rand0 = RandStartArb(rand0);          //place.

  EventStartTime(t0);                        //Initialize the event queues.

  t = t0;                                    //Set the starting time
//stid = is1+is2;                            //Set first available new strain
                                             //strain type ID.
  InitPop();                                 //Set up initial population.

  Report(argv[0]); pt = t;                   //Report initial conditions.

  BirthG();                                  //Start external event generators
  ImmigrateG();                              //for birth and immigration.

  for(t=t0; t<t1; Dispatch())                //Main loop: process events,
  { if(t-pt<tgap) continue;                  //reporting results periodically.
    pt = t; Report(argv[0]); }

  Report(argv[0]);                           //Get final report.

  Final();                                   //Close processing and return to
  free(A);                                   //caller.

  if(fit5i)                                  //If linked with the fitter, return
  { if(fitm) return out;                     //an array of notification rates or
    else     return outn; }                  //numbers for fitting.

  return 0;                                  //Otherwise return a success code.
}

/*----------------------------------------------------------------------------*
DISPATCH NEXT EVENT

All events pass through this routine. It picks the earliest event in the
queue, sets the time to match that event, and performs the operations called
for by that event. Typically that will result in other events being scheduled,
to be seen in the future as they arrive at the front of the queue.

ENTRY: The system is initialized with all events in the list ready for
         processing.
       't' contains the present time.
       't1' contains the ending time.

EXIT:  The next event has been processed and 'events' incremented, if the
         event's time is less than 't1'.
       't' is advanced to the next event, which may be an unprocessed event at
         time greater than 't1'.
*/

Dispatch()
{ int n; dec tw;

//-printf("About to dispatch an event...\n"); fflush(stdout);
  tw = t;                                    //Remember the previous time.
  n = EventNext(); if(t>t1) return;          //Advance time to the next event.
  tstep(tw, t);                              //Record the size of the time step.
  events += 1;                               //Increment the events counter.

  switch(A[n].pending)                       //Process the event.
  { case pVaccin:   Vaccination(n);  break;  //[vaccination]
    case pTransm:   Transmission(n); break;  //[transmission of an infection]
    case pRemote:   Remote(n);       break;  //[transition to latency]
    case pDisease:  Disease(n);      break;  //[progression to disease]
    case pDeath:    Death(n);        break;  //[death]
    case pMutate:   Mutate(n);       break;  //[strain type mutation]
    case pEmigrate: Emigrate(n);     break;  //[emigration from UK]
    case pBirth:    BirthG();        break;  //[birth generator]
    case pImmig:    ImmigrateG();    break;  //[immigration generator]
    case pRep:      Rep(n);          break;  //[case report]

    default: Error2(921.2,                   //[system error]
      "`A[",n,"].pending=",A[n].pending); }
}



/*----------------------------------------------------------------------------*
BIRTH

This routine is dispatched when an individual is to be born. All newborns are
Uninfected; exit from the Uninfected compartment is by vaccination to the
Immune compartment, by infection to the Recent Infection compartment, and by
death or emigration from the study population.

ENTRY: 'n' indexes an individual being born.
       'b' contains the time of birth. Presently, this is the current time,
         though with some set up, it could be earlier than present (notably,
         'pmale' would have to be indexed differently).
       't' contains the current time.
       'A[n].state' contains the present state of the record (can be any state,
         including 0, which is not a state but marks records not yet assigned).
       'm1' contains the mortality rate for susceptible individuals (if app.).
       'em' contains the emigration rate.
       'v1' contains vaccine efficacy.
       'v2' contains the probability that an individual will be vaccinated.
       'v3' contains the average age of vaccination.
       'VTYPE' is zero if vaccinations are to match ODE conventions.
       No event is scheduled for individual 'n'.

EXIT:  'Birth' contains a status code.
         0 The individual would die before the current time so no birth has been
            recorded and no event scheduled.
         1 Entry 'n' is initialized as a susceptible newborn and its first event
            is scheduled, either vaccination, emigration or death.
       'A[n].state' marks a susceptible individual.
       Counters in 'N' are updated.
*/

#define VTYPE 1                              //Vaccination type.

int Birth(int n, dec b)
{ int y, s, v, e; dec wd, we, wv;

//- printf("Starting birth routine \n");

  if(n<maximm+1) Error1(610.1, "n=",(dec)n); //Check for appropriate 'n', this
  if(n>indiv)    Error1(610.2, "n=",(dec)n); //routine does not allow immigrant
                                             //births or births to those with
                                             //index number >'indiv'.

  //-A[n].bstate = 0;                        //Clear bstate.
  y = (int)t - (int)t0;                      //Retrieve year index for arrays.
  A[n].sex = Rand()<pmale[y]? 0: 1;          //Assign the newborn's sex.
  s = A[n].sex;

  //-printf("Mortality rate: %f\n", m1[s][y]);
  wd = b+LifeDsn(s,t-b,m1[s][y]);            //Schedule a time of death and
  if(wd<t) Error(850.);                      //check for errors.
  //-we = b+Expon(em[s][1]);                 //Schedule time to emigration.
  we = b+EmDsn(1,s,t-b,em[s][UK]);           //Calculate time of emigration.

  A[n].tBirth    = b;                        //Record the time of birth.
  A[n].tDeath    = wd;                       //Record the time of death.
  //-A[n].tEntry = b;                        //Record time of entry into state.
  A[n].tEmigrate = we;                       //Record the time of emigration.
  A[n].tExit     = 0;                        //Clear any other saved event
  A[n].tDisease  = 0;                        //times.
  A[n].tTransm   = 0;
  A[n].tMutate   = 0;
  //-A[n].tInfected = 0;
  A[n].rob       = 1;                        //Set as born in UK.
  NewState(n, qU);                           //Mark as Uninfected.

  v = 0; switch(VTYPE)                       //Select the type of vaccination
  {                                          //scheduling.
//-/*  case 0:
//-      wv = b+Expon(1/v3[UK]);
//-      if(Rand()<(v1[UK]*v2[UK]) && wv<wd && wv<we) v = 1;
//-      break;
//-// This does not produce the same results as the one below!!
//-*/
  case 0:
    wv = b+Expon(v1[UK]*v2[UK]/v3[UK]);      //Generate a time for vaccination
    if(wv<wd && wv<we) v = 1;                //compatible with ODE models (for
    break;                                   //testing).

  case 1:                                    //Generate a vaccination sometime
    wv = b+v3[UK]+Rand();                    //within the specified year if
    if(b<1993 && Rand()<(v1[UK]*v2[UK])      //probabilities allow.
              && wv<wd && wv<we) v = 1;
    break;

  default: Error1(611., "",(dec)VTYPE);      //Improper vaccination type.
  }

  if(v)                                      //If vaccination occurs before
  { A[n].pending = pVaccin;                  //death and emigration, schedule
//-printf("About to schedule vaccination from Birth()\n"); fflush(stdout);
    EventSchedule(n, wv);                    //the vaccination.
    return 1; }

  if(we<wd)                                  //Schedule emigration if that
  { A[n].pending = pEmigrate;                //is the earliest event.
//-printf("About to schedule emigration from Birth()\n"); fflush(stdout);
    EventSchedule(n, we);
    return 1; }

  { A[n].pending = pDeath;                   //Otherwise, schedule death.
//-printf("About to schedule death from Birth()\n"); fflush(stdout);
    EventSchedule(n, wd);
    return 1; }
}



/*----------------------------------------------------------------------------*
IMMIGRATION

This routine brings a new individual into the population. The individual is
assigned all demographic and infection-related attributes according to the
appropriate probability distributions. They are also scheduled for their
earliest event. All information stored for this individual is written over, in
case their index number is being recycled from an individual leaving the study
population through death or emigration.

ENTRY: 'n' contains the index number of the new immigrant. The contents of the
         record is undefined.
       't' contains the current time.
       'indiv' contains the highest index number for any individual.
       'maximm' contains the highest index number for an immigrant.
       't0' contains begining time of the simulation.
       'SSAV' contains model version, 0=non-SSA, 1=SSA
       'ssaim' contains proportion of non-UK born immigrants from SSA by year.
       'immsex[r]' contains the proportion of immigrants who are male,
         r=0, non-UK born; r=1, UK-born; r=2, SSA-born. Note in SSA & non-SSA
         versions of the model, 'non-UK born' will be defined differently.
       'hivp' contains the HIV prevalence, by sex & year, for SSA immigrants.
       'is1' contains the possible strain types for immigrants, 0 to (is1-1).
       'infimm' contains cumulative probabilities immigrants enter disease
         states, by age, rob and calendar year.
       'Ax' contains state variables which accompany 'infimm'.
       'M1' contains the mortality table for non-diseased (real runs only).
       'A1' contains state variables which accompany 'M1'.
       'm1' contains the mortality rate (ODE validation only).
       'em' contains the emigration rate.
       'v1' contains vaccine efficacy.
       'v2' contains the probability that an individual will be vaccinated.
       'v3' contains the average age of vaccination.
       No event is scheduled for individual 'n'.

EXIT:  An event is scheduled for individual 'n'.
       'A[n].state' contains the disease state.
       'A[n].tEntry' contains the time of entry to the new state.
       'A[n].tBirth' contains the time of birth.
       'A[n].tImm' contains the time of immigration.
       'A[n].tDeath' contains the time of death.
       'A[n].tEmigrate' contains the time of emigration.
       'A[n].sex' contains the sex.
*/

int Immigrate(int n)
{ int y,s,rob,rob2,ac,a,st; dec r,age,wd,we,wv,tinf;

//- printf("Starting Immigrate routine...\n"); fflush(stdout);
  if(n>indiv) Error1(610.3, "n=",(dec)n);    //Check for appropriate n.
  if(n<1)     Error1(610.4, "n=",(dec)n);

  //-A[n].bstate = 0;                        //Clear bstate.
// SET UP BASIC, UNINFECTED IMMIGRANT

  NewState(n,qU);                            //Assign to Uninfected state
                                             //to start with.
  //- printf ("time (t) is %f\n",t);
  //-printf("tDeath is %f\ttEmigrate is %f\n",A[n].tDeath,A[n].tEmigrate);
  y = (int)t - (int)t0;                      //Get array index for year.
  //-A[n].tEntry = t;                        //Time of entry into this state.
  //-A[n].tImm   = t;                        //Assign time of immigration.
  if(n<=maximm) A[n].rob = rob = 0;          //Assign rob=0 to all non-UK born.
  else          A[n].rob = rob = 1;          //Assign rob=1 to UK-born.

  s = 0;                                     //Set sex as male to begin.

  if(rob==0 && SSAV==1)                      //If non-UK born and running
  { A[n].ssa = 0;                            //'SSA' version of model, check
    if(Rand()<ssaim[y])                      //to see if SSA. If so, get their
    { A[n].ssa = 1;                          //sex and HIV status.
      if(Rand()>immsex[y][SSA]) s = 1;       //Get sex of ssa.
      if(Rand()<hivp[s][y]) A[n].ssa = 2; }  //Get HIV status of SSA.
    else                                     //If not SSA, leave as other non-
      if(Rand()>immsex[y][0]) s = 1;         //UK and assign sex.
  }                                          //*Might split into SSAV/nonSSAV for clarity!

  else                                       //Assign sex to UK-born and
    if(Rand()>immsex[y][rob]) s = 1;         //non-UK born in non-SSA model.

  A[n].sex = s;                              //Assign sex to record.

  age=GetAge(n,s,rob);                       //Assign age.
  a = (int) age;                             //Save integer age.
/*-printf("GetAge() returned %f when n=%d, s=%s, and rob=%d\n",age,n,s,rob);
* if(age<5)       ac=0;                      //Get age class for 'infimm'.
* else if(age<15) ac=1;                      //(For old way of storing)
* else if(age<30) ac=2;
* else if(age<45) ac=3;
* else if(age<65) ac=4;
* else            ac=5;
* -printf("age is %f and age class is %d\n",age,ac);
*/

  rob2=rob;                                  //Create 'rob2' (0,1 and 2)
  if(SSAV && A[n].ssa) rob2=2;               //since here 'rob' is only 0 or 1).

  A[n].tBirth = t-age;                       //Save time of birth based on age.

  A[n].tDeath = wd                           //Assign time of death and check
              = t+LifeDsn(s,age,m1[s][y]);   //death time is ok.
  //-A[n].tEmigrate = we = t+Expon em[s][rob2]  //Assign time of emigration.
  if(wd<A[n].tBirth+age) Error(612.1);

  A[n].tEmigrate                             //Assign emigration time.
       = we
       = t+EmDsn(rob2,s,age,em[s][rob2]);
//-printf("tDeath = %f\ttEmigrate = %f\n",A[n].tDeath,A[n].tEmigrate);

  if(age<v3[rob] && Rand()<v1[rob]*v2[rob]   //Determine if vaccination should
                 && t<2005-(v3[rob]-age))    //occur and assign vaccination
    wv = t+(v3[rob]-age)+Rand();             //time if so. If it should not
  else wv = t+2*RT+Rand();                   //occur, set to time which will
                                             //not happen in the model.
  if(wv<wd && wv<we)                         //Schedule vaccination if it is
  { A[n].pending = pVaccin;                  //the earliest event.
//-printf("About to schedule vaccination from Immigration()\n"); fflush(stdout);
    EventSchedule(n, wv); }

  else if(wd<we)                             //Schedule death if it is the
  { A[n].tExit = wd;                         //earliest event.
    A[n].pending = pDeath;
//-printf("About to schedule death from Immigration()\n"); fflush(stdout);
    EventSchedule(n, wd); }

  else                                       //Otherwise schedule emigration if
  { A[n].tExit = we;                         //it is the earliest event.
    A[n].pending = pEmigrate;
//-printf("About to schedule emig from Immigration()\n"); fflush(stdout);
    EventSchedule(n, we); }

  A[n].tDisease  = 0;                        //Clear time to disease.
  A[n].tTransm   = 0;                        //Clear time to transmit.
  A[n].tMutate   = 0;                        //Clear time of strain mutation.
  //-A[n].tInfected = 0;                     //Clear time of infection.

  //if(rob2==0 && A[n].ssa<1 && Rand()<0.0)  //Add correction so that a fraction
  //  rob2=1;                                //of immigrants from low-risk
                                             //countries are treated as though
                                             //they are UK-born for disease
                                             //state assignment purposes.

//ASSIGN DISEASE STATE TO IMMIGRANT AND PROCESS ACCORDINGLY

  st = 1+
    (int)RandF(Ax,infimm[a][rob2][y],9,1.0);   //Get random disease state.
//-printf("Disease state is %d\n",st); fflush(stdout);

  if(st==1)                                  //Do nothing if Uninfected.
    return 0;
  else if(st==2)                             //Process Immune.
  { EventCancel(n);
    Vaccination(n);
    return 1; }
  else if(st==3)                             //Process Recently infected.
  { tinf = Rand()*5;                         //Set random time since infection
    Infect(n,tinf,0);                        //within the last five years.
    //-A[n].inf=0;
    return 2; }
  else if(st==4)                             //Process Remotely infected.
  { EventCancel(n);
    NewState(n,qD1);                         //Put in temporary disease state
    Remote(n);                               //to facilitate re-scheduling of
    //-A[n].inf=0;                           //events in Remote() function.
    return 3; }
  else if(st==5)                             //Process Reinfected.
  { NewState(n,qI2);                         //Put in 'Remote infection' so
    tinf = Rand()*5;                         //that Infect() picks this up
    Infect(n,tinf,0);                        //as reinfection, also draw
    //-A[n].inf=0;                           //random infection time.
    return 4; }
  else if(st>5 && st<9)                      //Process Primary,Reactivation,
  { EventCancel(n);                          //Reinfection disease classes.
    NewState(n,st-3);                        //First put in correct infection
    Disease(n);                              //class and then send to Disease().
    //-A[n].inf=0;
    return 5; }
  else
    Error(618.1);
    return 0;                                //(Will never reach this).
}



/*----------------------------------------------------------------------------*
VACCINATION

This routine is dispatched when an individual is scheduled for an effective
vaccination. Ineffective vaccinations (wasted vaccinations) are already
accounted for and never scheduled to arrive at this routine. Effective
vaccinations are assumed to impart lifelong immunity; therefore individuals
will never leave this state, except by dying or emigrating from the
population.

ENTRY: 'n' indexes an individual being born.
       't' contains the current time.
       'A[n].state' contains the present state (always 'qU').
       'A[n].sex' contains the individual's sex.
       'A[n].tEmigrate' contains time of emigration.
       'A[n].tDeath' contains time of death.
       No event is scheduled for individual 'n'.

EXIT:  'A[n].state' contains the new state, qV.
       Counters in 'N' are updated.
*/

int Vaccination(int n)
{
//- printf("Starting vaccination routine...\n"); fflush(stdout);

  NewState(n, qV);                           //Change states.

  if(A[n].tEmigrate<A[n].tDeath)             //Schedule emigration if that is
  { A[n].pending = pEmigrate;                //the earliest event.
//-printf("About to schedule emigration from Vaccination()\n"); fflush(stdout);
    EventSchedule(n, A[n].tEmigrate); }
  else                                       //Otherwise, schedule death.
  { A[n].pending = pDeath;
//-printf("About to schedule death from Vaccination()\n"); fflush(stdout);
    EventSchedule(n, A[n].tDeath); }

  return 0;
}


/*----------------------------------------------------------------------------*
INFECT A SPECIFIED INDIVIDUAL

Individuals receive infections from others via this routine. If the targeted
individual is uninfected (U) or remotely infected (I2), it acquires a new
infection and moves to compartment I1 or I3. That individual is then scheduled
to return to remote infection, develop disease, emigrate, die, or to have its
infection strain mutate, depending on probabilities of each and random chance.
If they are not susceptible to infection, then the transmission event has no
consequences.

Idea:
    Change arguments to Infect() to include place of infection (char 'poi')
    and possibly time of infection? Then Infect() can be better used to infect
    those immigrating and those infected from initialization. Could leave
    default '1' for UK-infected and 't' for infected at current time if that
    is possible to specify some params and leave others as default). Could
    also have two different function prototypes (C++)??.

ENTRY: 'n' indexes the individual to be infected.
       'tinf' contains the time of infection, <=5 years ago.
       'strain' contains the strain ID number of infecting strain.
       't' contains the current time.
       'A[n].state' contains the state of infection target.
       'A[n].tEmigrate' contains the time of emigration of infection target.
       'A[n].tDeath' contains the time of death for infection target.
       'r1' and 'r2' contain the latency rates for 'qI1' and 'qI3'.
       'mi' contains the mutation rate for infection strains.
       An event is still scheduled for individual 'n'.

EXIT:  'Infect' describes the result.
         0 The specified individual could not be infected or reinfected.
         1 Return to remote infection is scheduled.
         2 Disease is scheduled.
         3 Death is scheduled.
         4 Strain mutation is scheduled.
         5 Emigration is scheduled.
       'A[n].state' contains the new state, if applicable.
       Entry 'n' is infected (if 'Infect' is nonzero).
       Counters in 'N' are updated.
*/

int Infect(int n, dec tinf, int strain)
{ int s, rob, a, q; dec d, r, wd, we, wdis, wr, wm;

//-printf("Starting Infect routine...\n"); fflush(stdout);
//-printf("n=%d\tstrain=%d\n",n,strain); fflush(stdout);

  if(n>indiv||n<1)   Error1(610.3,"",n);     //Check for appropriate n.
  if(strain>stid)    Error1(616.0,"",strain);//Check for appropriate strain ID.
  if(tinf>5||tinf<0) Error1(617.0,"",tinf);  //Check for appropriate 'tinf'.
  if(tinf==5) tinf=tinf-E;                   //Correct 'tinf' if equal to 5.

  s   = A[n].sex;                            //Retrieve sex.
  rob = A[n].rob;                            //Retrieve region of birth.
  a   = (int)(t-A[n].tBirth);                //Retrieve integer age.

  switch(A[n].state)                         //Determine the new state and
  {                                          //its associated parameters.
    //-case qI2: d=d3[s][rob][a]; r=r2[s]; q=qI3; break;
    //-case qU:  d=d1[s][rob][a]; r=r1[s]; q=qI1; break;
    case qI2: r=r2[s]; q=qI3; break;
    case qU:  r=r1[s]; q=qI1; break;
    default:  return 0;                      //Avoid uninfectable states.
  }

//-printf("Calling EventCancel from Infect(), n=%d\n", n); fflush(stdout);
//-printf("Info on n: A[n].pending = %d at time %f, current time is %f\n",
//-A[n].pending, A[n].t,t);

  EventCancel(n);                            //Else cancel the pending event and
  NewState(n, q);                            //mark this individual as infected.
  //A[n].tInfected = t-tinf;                 //Save time of infection.
  //A[n].inf        = 1;                     //Save place of infection as UK
                                             //(will be changed outside routine
                                             //if infection is acquired abroad).
  //-A[n].strain = strain;                   //Assign strain type ID number.

  wd = A[n].tDeath;                          //Retrieve time of death.
  we = A[n].tEmigrate;                       //Retrieve time of emigration.
//-wr = t+Expon(r);                          //Calculate recovery (old way).
  wr = t+LAT-tinf;                           //After 'LAT' years individual is
                                             //defined as remotely infected.
//-wdis = t+Expon(d);                        //Calculate disease (old way).
  wdis = t+Tdis(n,a,s,rob,tinf)+E;           //Calculate time to disease.
  if(wdis<=t) Error2(620.0,"t=",t, " wdis=",wdis);
  wm   = t+Expon(mi);                        //Calculate strain mutation time.

  if(wd<we && wd<wr && wd<wdis && wd<wm)     //If death is earliest event,
  { A[n].pending = pDeath;                   //schedule the death and
//-printf("About to schedule death from Infect()\n"); fflush(stdout);
    EventSchedule(n, wd);                    //ignore everything else.
    return 3; }

  if(we<wr && we<wdis && we<wm)              //If emigration is the earliest
  { A[n].pending = pEmigrate;                //event, schedule it and
//-printf("About to schedule emigration from Infect()\n"); fflush(stdout);
    EventSchedule(n, we);                    //ignore everything else.
    return 5; }

  if(wr<wdis && wr<wm)                       //Otherwise, if transition to
  { A[n].pending = pRemote;                  //remote infection would occur
//-printf("About to schedule remote from Infect()\n"); fflush(stdout);
    EventSchedule(n, wr);                    //before disease and mutation,
    A[n].tMutate = wm;                       //schedule latency, save mutation
    return 1; }                              //time, and ignore disease.

  if(wm<wdis)                                //Otherwise, if mutation should
  { A[n].pending = pMutate;                  //occur before disease, schedule
//-printf("About to schedule mutation from Infect()\n"); fflush(stdout);
    EventSchedule(n, wm);                    //mutation and save time to disease
    A[n].tDisease = wdis;                    //onset and time to remote.
    A[n].tExit = wr;
    return 4; }

  { A[n].pending = pDisease;                 //Otherwise, schedule disease and
//-printf("About to schedule disease from Infect()\n"); fflush(stdout);
    EventSchedule(n, wdis);                  //do not save others, as they will
    return 2; }                              //be recalculated at disease onset.
}

/*----------------------------------------------------------------------------*
ENTER COMPARTMENT REMOTE

This routine is dispatched when an infection becomes 'Remote Infection,'
entering compartment I2. Recent Infection (I1), Reinfection (I3), and all
disease compartments (D1-D6) can lead to this state. The Remote Infection state
has four events which can be scheduled -- disease, strain mutation, death or
emigration. Also, remotely-infected individuals can be reinfected but this is
induced by a transmission event dispatched independently, so is not handled
here.

ENTRY: 'n' indexes the individual.
       't' contains the current time.
       'A[n].state' contains the present state (can be 'qI1', 'qI3',
         'qD1'-'qD6').
       'A[n].tMutate' contains the strain mutation time.
       'd2' contains the disease progression rate for 'qI2'.
       'm4' contains the mortality rate for 'qI2'.
       No event is scheduled for individual 'n'.

EXIT:  'Remote' contains a status code:
         2 Disease is scheduled.
         3 Death is scheduled.
         4 Mutation is scheduled.
         5 Emigration is scheduled.
       'A[n].state' represents Remote Infection (qI2).
       'A[n].tDeath' is updated as necessary.
       'A[n].tMutate' is updated as necessary.
       Counters in 'N' are updated.
*/

int Remote(int n)
{ int y, a, s, rob, q; dec age, wdis, wd, we, wm;

//-printf("Starting Remote routine...\n"); fflush(stdout);
  y   = (int)t - (int)t0;                    //Retrieve array index for year.
  age = t-A[n].tBirth;                       //Retrieve age.
  a   = (int) age;                           //Integer age.
  s   = A[n].sex;                            //Retrieve sex.
  rob = A[n].rob;                            //Retrieve region of birth.

  q = A[n].state;                            //Remember the previous state.
  NewState(n, qI2);                          //Mark the individual as remote.

  if(q>=qD1)                                 //Establish a new time for strain
  {  //A[n].tDeath = t+LifeDsn(s,age,m1[s][y]);//mutation if the prior state
     A[n].tMutate = t+Expon(mi); }           //was disease.

//-wdis = t+Expon(d2[s][rob][a]);            //Calculate time to disease (old).
  wdis = t+Tdis(n,a,s,rob,0);                //Calculate time to disease.
  wd = A[n].tDeath;                          //Retrieve time of death.
  we = A[n].tEmigrate;                       //Retrieve time of emigration.
  wm = A[n].tMutate;                         //Retrieve time of mutation.

//-printf("wdis = %f\twd = %f\twe = %f\twm = %f\n",wdis,wd,we,wm);

  if(wd<wdis && wd<wm && wd<we)              //If death would occur before
  { A[n].pending = pDeath;                   //disease, strain mutation, and
//-printf("About to schedule death from Remote()\n"); fflush(stdout);
    EventSchedule(n, wd);                    //emigration, then schedule death.
    return 3; }

  if(wm<wdis && wm<we)                       //Otherwise, if strain mutation
  { A[n].pending = pMutate;                  //occurs before disease and
//-printf("About to schedule mutation from Remote()\n"); fflush(stdout);
    EventSchedule(n, wm);                    //emigration, schedule strain
    A[n].tDisease = wdis;                    //mutation and save disease time.
    return 4; }

  if(we<wdis)                                //Otherwise, if emigration occurs
  { A[n].pending = pEmigrate;                //first then schedule that and
//-printf("About to schedule emigration from Remote()\n"); fflush(stdout);
    EventSchedule(n,we);                     //ignore everything else.
    return 5; }

  { A[n].pending = pDisease;                 //Otherwise, schedule progression
//-printf("About to schedule disease from Remote()\n"); fflush(stdout);
    EventSchedule(n, wdis);                  //to disease and ignore
    return 2; }                              //everything else.
}


/*----------------------------------------------------------------------------*
DISEASE

This routine is dispatched when an infection progresses to active disease.
There are six distinct disease compartments, primary (D1), reactivation (D2),
and reinfection (D3), plus non-pulmonary classes which correspond to each of
these (D4, D5, and D6). Individuals enter disease compartments from three
infection compartments -- I1, I2, and I3 -- which determine the disease
compartment they enter. This routine handles all transitions to disease.

Future events for diseased individuals are scheduled here. All transmission
originates from the disease compartments. In addition to infecting others,
diseased individuals can recover to remote infection, die, emigrate, be
reported or have their infection strain mutate.

ENTRY: 'n' indexes the individual progressing to disease.
       't' contains the current time.
       'A[n].state' contains the state progressing to disease (can be
         'qI1', 'qI2', or 'qI3').
       'A[n].tBirth' contains the time of birth of the individual.
       'A[n].sex' contains the sex of the individual.
       'A[n].rob' contains region of birth of the individual.
       'A[n].tEmigrate' contains the time of emigration.
       'A[n].tDeath' contains the scheduled time of death.
       'cft' contains the case fatality rate (actually a proportion).
       'r3', 'r4', 'r5','r6','r7', and 'r8' contain recovery  rates for 'qD1',
       'qD2','qD3', 'qD4','qD5', and 'qD6' respectively.
       'm6', 'm7', 'm8', 'm9','m10', and 'm11' contain mortality rates for the
         states named above in the ODE-compatible version of the model.
       'p1', 'p2', and 'p3' contain the fraction of disease which becomes
         pulmonary, from sources 'I1', 'I2', and 'I3', respectively.
       'c' contains the average number of new infections produced by this
         individual per year.
       'proprep' contains the proportion of cases reported.
       'md' contains the mutation rate for strains involved in disease.
       No event is scheduled for individual 'n'.

EXIT:  'Disease' contains a status code.
         1 A transmission is scheduled.
         2 Recovery is scheduled.
         3 Death is scheduled.
         4 Strain mutation is scheduled.
         5 Emigration is scheduled.
         6 Case report is scheduled.
       'A[n].state' contains the new state.
       'A[n].tDeath' contains a possibly updated time of death.
       'A[n].tMutate' contains the new time of strain mutation, if applicable.
       'A[n].tTransm' contains the next time of transmission, if applicable.
       'A[n].tExit' contains the time of recovery to remote infection, or if
         would happen after death, the time of death.
       'A[n].tRep' contains the time of disease case report, if applicable.
       Counters in 'N' are updated.
*/

int Disease(int n)
{ int a, s, rob, y, ds, q; dec age, r, m, p, wm, we, wd, wr, wt, e, wrep;

//-printf("Starting  Disease routine...\n"); fflush(stdout);

  age = t-A[n].tBirth;                       //Retrieve age.
  a   = (int)age;                            //Calculate integer age.
  s   = A[n].sex;                            //Retrieve sex.
  rob = A[n].rob;                            //Retrieve region of birth, 0 or 1.
  y   = (int)t - (int)t0;                    //Retrieve array index for year.

  switch(A[n].state)                         //Determine the new state and its
  {                                          //associated parameters.
    case qI1: r=r3[s]; m=m6[s][y]; p=p1[a][s][rob]; q=qD1; break;
    case qI2: r=r4[s]; m=m7[s][y]; p=p2[a][s][rob]; q=qD2; break;
    case qI3: r=r5[s]; m=m8[s][y]; p=p3[a][s][rob]; q=qD3; break;
    default:  Error(922.0); }

  if(Rand()>p)                               //If this should be non-pulmonary
  { switch(A[n].state)                       //disease, update new state and
    {                                        //associated parameters.
      case qI1: r=r6[s]; m= m9[s][y]; q=qD4; break;
      case qI2: r=r7[s]; m=m10[s][y]; q=qD5; break;
      case qI3: r=r8[s]; m=m11[s][y]; q=qD6; break;
      default:  Error(922.0); } }

  NewState(n, q);                            //Mark the individual as diseased.
//-printf("About to put case in cumulative cases....\n"); fflush(stdout);
  Cumul(n,t);                                //Add individual to cumulative cases.                  *Note name....

  wr = A[n].tExit = t+RecovDsn(s,age,r);     //Establish time to remote.
  we = A[n].tEmigrate;                       //Retrieve emigration time.
  wd = A[n].tDeath;                          //Retrieve time of death.
  A[n].tMutate = wm = t+Expon(md);           //Establish new mutation time.

  if(q>=qD4) ds = 0;                         //Set disease site to non-pulmonary
  else       ds = 1;                         //or pulmonary.

  if(Rand()<cft[a][ds][y])                   //If person should die from disease,
  { if(wr<wd && wr<we) e = wr;               //find earliest of death, emigration
    else if(wd<we)     e = wd;               //and recovery to assign death
    else               e = we;               //before these occur, assigning
                                             //death time close to the end of
                                             //disease duration.
    wd = t+0.99*(e-t);

//-printf("Disease death time is %f and recovery time is %f\n",wd,wr);
                                             //Since disease death is before
    A[n].tDeath = wd; }                      //natural death time, replace it.

  if(Rand()<proprep)                         //If case should be reported, find
  { if(wr<wd && wr<we) e = wr;               //reporting time.
                                             //First find earliest of death,
    else if(wd<we)     e = wd;               //emigration and disease recovery
                                             //to get range for reporting time.
    else               e = we;

    A[n].tRep = t+Rand()*(e-t); }            //Randomly assign reporting time.

  else                                       //If case should not be reported,
    A[n].tRep = t+2*RT+Rand();               //assign a reporting time beyond
                                             //running time of model.

  if(A[n].tRep==0) Error1(619., "n=",n);
  wrep = A[n].tRep;                          //Save reporting time.
//-printf("tRep is %f\n", A[n].tRep); fflush(stdout);

  if(wd<wr)      /* Delete if 'Earliest' */  //If death would occur before
     wr = wd;    /* is incorporated.     */  //recovery, give death precedence.

//-printf("About to calculate tRep...\n"); fflush(stdout);
  if(q<qD4 && Rand()<smear[a])               //If this is pulmonary disease and
    wt = t+Expon(c[s][rob]);                 //it is smear positive, set time
  else wt = t+2*RT + Rand();                 //to transmit if smear negative,
  A[n].tTransm = wt;

  if(wt<wr && wt<wm && wt<we && wt<wrep)     //set 'wt' so transmission never
  { A[n].pending = pTransm;                  //happens. If transmission is
//-printf("About to schedule transmission from Disease()\n"); fflush(stdout);
    EventSchedule(n, wt);                    //earliest even, schedule it.
    return 1; }

  if(wrep<wr && wrep<wm && wrep<we)          //If case reporting will occur
  { A[n].pending = pRep;                     //first, schedule it, save the
//-printf("About to schedule case report from Disease() with time: %f\n",wrep);
    EventSchedule(n, wrep);                  //mutation time.
    return 6; }

  if(wr<wd && wr<wm && wr<we)                //If recovery will occur before
  { A[n].pending = pRemote;                  //death, emigration, and mutation,
//-printf("About to schedule remote from Disease()\n"); fflush(stdout);
    EventSchedule(n, wr);                    //schedule recovery.
    return 2; }

  if(wm<wd && wm<we)                         //If mutation will occur before
  { A[n].pending = pMutate;                  //death and emigration, schedule
//-printf("About to schedule mutation from Disease()\n"); fflush(stdout);
    EventSchedule(n, wm);                    //mutation.
    return 4; }

  if(we<wd)                                  //If emigration will occur before
  { A[n].pending = pEmigrate;                //death, schedule emigration and
//-printf("About to schedule emigration from Disease()\n"); fflush(stdout);
    EventSchedule(n,we);                     //ignore all other events.
    return 5; }

  { A[n].pending = pDeath;                   //Otherwise schedule the
//-printf("About to schedule death from Disease()\n"); fflush(stdout);
    EventSchedule(n, wd);                    //individual's death and ignore
    return 3; }                              //all other events.
}



/*----------------------------------------------------------------------------*
TRANSMISSION

Infectious individuals transmit infection via this routine. A target individual
is selected to be infected, either randomly from within the same region of birth
as the infectious individual (only divided into UK-born & non-UK born, whether
or not running SSA version of model) or randomly from the entire population. If
the target individual is susceptible (Uninfected or Remote Infection states),
infection is established. If not, the infection is lost.

The original infectious individual is then scheduled for another transmission,
strain mutation, recovery, case report emigration or death.

ENTRY: 'n' indexes the individual to transmit an infection.
       't' contains the current time.
       'A[n].tExit' contains the time of recovery, or if that time
         is equal to or greater than the time of death, contains time of death.
       'A[n].tDeath' contains the time of death.
       'A[n].tMutate' contains the strain type mutation time.
       'A[n].tEmigrate' contains the emigration time.
       'A[n].rob' contains the region of birth.
       'A[n].sex' contains the sex.
       'A[n].strain' contains the infection strain ID number.
       'pcc' contains the proportion of close contacts.
       Non-UK born individuals are indexed from 1 to (immid-1), a total of
         (immid-1) individuals.
       UK-born individuals are indexed from (maximm+1) to (ukbid-1), a total
         of (ukbid-maximm-1) individuals.
       No event is scheduled for individual 'n'.'

EXIT:  A new individual has been targeted for infection. If the infection takes
         hold, that individual is scheduled for strain mutation, disease,
         remote infection, emigration or death.
       'Transmission' contains a status code.
         1 Another infection is scheduled.
         2 Recovery to remote infection is scheduled.
         3 Death is scheduled.
         4 Strain mutation is scheduled.
         5 Emigration is scheduled.
         6 Case report is scheduled.
       Counters in 'N' are updated.
*/

#define SCHED(X,Y,Z)  { A[n].pending = X; EventSchedule(n,Y); return Z; }

int Transmission(int n)
{ int i, j, low, tot; dec age;
  static int v[] = { iTransm,iDeath,iEmigrate,iExit,iMutate,iRep, -1 };

//- printf("Starting Transmission routine...\n"); fflush(stdout);

  if(Rand()<pcc)                             //If targetting 'close contact'
  { if(A[n].rob)                             //choose random individual from
    { low  = maximm+1;                       //individual's own region of birth.
      tot = (ukbid-1) - low + 1; }
    else
    { low  = 1;
      tot = immid-1 - low + 1; }

    do i=low+(int)(Rand()*tot);              //Find person other than self to
      while(i==n); }                         //infect.

  else                                       //If not a 'close contact', choose
  { do                                       //random person to infect from
    { tot = (immid-1) + (ukbid-maximm-1);    //entire population.
      j = 1 + (int)(Rand()*tot);
      if(j>=immid) i = j+(maximm+1-immid);   //Adjust ID numbers for UK-born.
      else i = j; }
    while (i==n);                            //Avoid infecting self.
  }
//-printf("About to infect chosen person from Trans()\n"); fflush(stdout);
//- Infect(i, A[n].strain);                  //Infect chosen individual.

                                             //(Infect for non-genetic model)
  Infect(i,0,0);                             //Infect chosen individual.

  A[n].tTransm=t+Expon(c[A[n].sex][A[n].rob]); //Establish time to transmit again.

  switch(i=Earliest(A[n].t, v))              //Schedule the earliest event.
  { case iRep:      SCHED(pRep,      A[n].tRep,      6);
    case iTransm:   SCHED(pTransm,   A[n].tTransm,   1);
    case iExit:     SCHED(pRemote,   A[n].tExit,     2);
    case iMutate:   SCHED(pMutate,   A[n].tMutate,   4);
    case iEmigrate: SCHED(pEmigrate, A[n].tEmigrate, 5);
    case iDeath:    SCHED(pDeath,    A[n].tDeath,    3);
    default:        Error1(922., "m=",(dec)i);           }

  return 0;                                //(Will never reach this.)
}



/*----------------------------------------------------------------------------*
MUTATION

This routine is dispatched when the strain type of an infected or diseased
individual mutates. The mutation does not affect any other event.

ENTRY: 'n' indexes an individual whose strain type is to mutate
       't' contains the current time.
       'mi' contains the mutation rate for strains not in an active disease
         case (merely infection).
       'md' contains the mutation rate for strains in a disease case.
       'A[n].state' contains the current state (can be any infection or disease
         state).
       'A[n].strain' contains the strain identification number of the current
         strain of infection or disease.
       'A[n].tDeath' contains the saved time of death.
       'A[n].tEmigrate' contains the saved time of emigration.
       'A[n].tExit' contains the saved time to exit state.
       'A[n].tTransm' contains the saved time to transmit again.
       No event is scheduled for individual 'n'.

EXIT:  The next event for individual 'n' is scheduled.
       'A[n].tMutate' contains the time of next scheduled strain mutation.
       'Mutation' contains a status code.
         1 Recovery to remote infection is scheduled.
         2 Progression to disease is scheduled.
         3 Death is scheduled.
         4 Strain mutation is scheduled.
         5 Emigration is scheduled.
         6 Case report is scheduled.
       Counters in 'N' are updated.
*/

Mutate(int n)
{ dec m, wm, wd, we, wdis, wr, wt, wrep;

//-printf("Starting Mutation routine...\n"); fflush(stdout);

  //-A[n].strain = stid;                     //Assign new, mutant strain type.
  stid += 1;                                 //Update next available strain type
                                             //ID number.
  if(A[n].state<=qI3) m = mi;                //Determine appropriate mutation
  else m = md;                               //rate and calculate time to
  wm = t+Expon(m);                           //mutate again.

  wd   = A[n].tDeath;                        //Get time of death.
  we   = A[n].tEmigrate;                     //Get time of emigration.
  wdis = A[n].tDisease;                      //Get time of disease.
  wr   = A[n].tExit;                         //Get time to remote infection.

  if(A[n].state==qI2)                        //Schedule events for remotely
  {                                          //infected individuals (qI2).
    if(wd<we && wd<wdis && wd<wm)
    { A[n].pending = pDeath;                 //If death would occur before
//-printf("About to schedule death from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wd);                  //emigration, disease and strain
      return 3; }                            //mutation, schedule death.

    if(wm<we && wm<wdis)                     //Otherwise, if strain mutation
    { A[n].pending = pMutate;                //occurs before disease and
//-printf("About to schedule mutation from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wm);                  //emigration, schedule mutation.
      return 4; }

    if(wdis<we)                              //Otherwise, if disease occurs
    { A[n].pending = pDisease;               //before emigration, schedule
//-printf("About to schedule disease from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wdis);                //progression to disease.
      return 2; }

    { A[n].pending = pEmigrate;              //Otherwise, schedule emigration.
//-printf("About to schedule emigration from Mutate()\n"); fflush(stdout);
      EventSchedule(n, we);
      return 5; }
  }

  if(A[n].state<=qI3)                        //Schedule events for the other
  {                                          //infected classes (qI1, qI3).
    if(wd<wdis && wd<wr && wd<wm && wd<we)
    { A[n].pending = pDeath;                 //If death is earliest
      EventSchedule(n, wd);                  //event, schedule the death and
//-printf("About to schedule death2 from Mutate()\n"); fflush(stdout);
      return 3; }                            //ignore everything else.

    if(wr<wdis && wr<wm && wr<we)            //Otherwise, if transition to
    { A[n].pending = pRemote;                //remote infection would occur
//-printf("About to schedule remote2 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wr);                  //before disease and mutation,
      A[n].tMutate = wm;                     //schedule latency and save
      return 1; }                            //mutation time.

    if(wm<wdis && wm<we)                     //Otherwise, if mutation should
    { A[n].pending = pMutate;                //occur before disease, schedule
//-printf("About to schedule mutation2 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wm);                  //mutation.
      return 4; }

    if(wdis<we)
    { A[n].pending = pDisease;               //Otherwise, if disease occurs
//-printf("About to schedule disease2 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wdis);                //before emigration, schedule
      return 2; }                            //progression to disease.

    { A[n].pending = pEmigrate;              //Otherwise, schedule emigration.
//-printf("About to schedule emigration2 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, we);
      return 5; }
  }


  {                                          //Schedule events for diseased.
    wrep = A[n].tRep;                        //Get time of case report.
//-printf("A[n].tRep is %f\n", A[n].tRep); fflush(stdout);

    if(A[n].state<qD4)                       //If this is pulmonary disease,
    { wt = A[n].tTransm;                     //retrieve time for transmission
      if(wt<wd && wt<wr && wt<wm && wt<we && wt<wrep)
      { A[n].pending = pTransm;              //and if it occurs before
//-printf("About to schedule tranms3 from Mutate()\n"); fflush(stdout);
        EventSchedule(n, wt);                //anything else, schedule it
        A[n].tMutate = wm;                   //and save mutation time.
        return 1; } }

    if(wrep<wd && wrep<wr && wrep<wm && wrep<we)
    { A[n].pending = pRep;                  //If case report should occur
//-printf("About to schedule case report3 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, wrep);               //before anything else, schedule
      A[n].tMutate = wm;                    //it and save mutation time.
      return 6; }

    if(wr<wd && wr<wm && wr<we)              //If recovery will occur before
    { A[n].pending = pRemote;                //death, emigration and mutation,
//-printf("About to schedule emote3 from Mutate()\n"); fflush(stdout);
    EventSchedule(n, wr);                    //schedule recovery.
    return 2; }

    if(wm<wd && wm<we)                       //If mutation will occur before
    { A[n].pending = pMutate;                //death and emigration, schedule
//-printf("About to schedule mutate3 from Mutate()\n"); fflush(stdout);
    EventSchedule(n, wm);                    //mutation.
    return 4; }

    if(wd<we)
    { A[n].pending = pDeath;                 //If death will occur before
//-printf("About to schedule death3 from Mutate()\n"); fflush(stdout);
    EventSchedule(n, wd);                    //emigration, schedule death.
    return 3; }

    { A[n].pending = pEmigrate;              //Otherwise, schedule emigration.
//-printf("About to schedule emigration3 from Mutate()\n"); fflush(stdout);
      EventSchedule(n, we);
      return 5; }
  }

}



/*----------------------------------------------------------------------------*
DEATH

This routine is dispatched when an individual dies, leaving the population.

ENTRY: 'n' indexes an individual who has just died.
       't' contains the current time.
       'A[n].state' contains the present state (can be any compartment).
       'A[n].tBirth' contains the time of birth.
       'ukbid' contains the next available index number for UK-born individuals.
       No event is scheduled for individual 'n'.

EXIT:  Either entry 'n' is sent to the Birth() function, to be initialized as a
        susceptible newborn and function returns '0' (DTYPE==0) or index number
       is recycled with 'Transfer' (DTYPE==1), no birth is generated and
       function returns '1'.
       N[A[n].state] is decremented.
       'deaths' is incremented.
       Counters 'age1', 'age2', and 'agec' are updated.
       Counters in 'N' are updated.
*/

#define DTYPE 1                              //Allows for non-constant population
                                             //size.
Death(int n)
{ int n2; dec age;

//- printf("Starting Death routine...\n"); fflush(stdout);

  deaths += 1;                               //Increment the number of deaths.
  N[A[n].state]-=1;                          //Decrement N[A[n].state].
  age = t-A[n].tBirth;                       //Compute the age at death.

  { age1[0] += age; age2[0] += age*age;      //Accumulate statistics for mean
    agec[0] += 1; }                          //age and its variance.

//-This may not make sense anymore?
//-/*  if(A[n].bstate<=1)                    //If this individual has never had
//-  { age1[1] += age; age2[1] += age*age;   //disease, also accumulate with
//-    agec[1] += 1; }                       //other non-diseased individuals.
//-*/
//-
  if(DTYPE==0)                               //If population size is to be held
  { Birth(n, t);                             //constant, initiate a birth.
    return 0; }

  if(A[n].rob)                               //Avoid unoccupied index numbers
  { n2 = ukbid-1; ukbid--; }                 //in array 'A' by transferring
                                             //highest-numbered individual, 'n2',
  else                                       //to index number 'n'.
  { n2 = immid-1; immid--; }

//-printf("About to call Transfer() from Death() with n=%d  n2=%d\n",n, n2);
  Transfer(n, n2);
  return 1;
}



/*----------------------------------------------------------------------------*
EMIGRATION

This routine logs individuals out of compartments as they leave the study
population, maintaining numbers in each compartment so that the array of
individuals never has to be scanned for that information. Also, the
individual's index number is recycled with 'Transfer' so that array 'A'
is always continuous.

ENTRY: 'n' indexes the individual.
       'tli' contains the time of last immigration.
       'A[n].state' contains the disease state.
       'A[n].rob' contains the region of birth.
       'immig[y]' contains the total number of immigrants each year.

EXIT: N[A[n].state] is decremented.
      'n' is recycled such that the highest-numbered individual within the
       same region of birth as 'n' takes the index number 'n' and array 'A'
       remains continuous.

*/

Emigrate(int n)
{ int n2;

//-printf("Starting Emigration routine...\n"); fflush(stdout);

  N[A[n].state] -= 1;                        //Decrement N[A[n].state].

  if(A[n].rob)                               //Use emigrant's region of birth
  { n2 = ukbid-1; ukbid--; }                 //to find highest index number of
                                             //individual who will take over
  else                                       //emigrant's index number, to
  { n2 = immid-1; immid--; }                 //prevent array from having gaps
                                             //of unoccupied index numbers.
//-printf("About to call Transfer() from Emigrate(), n=%d n2=%d\n", n,n2);
  Transfer(n, n2);
}



/*----------------------------------------------------------------------------*
IMMIGRATION GENERATOR

This routine coordinates an 'external' immigration event generator. The routine
uses a pseudo-individual (index number 'IMM') to schedule the external event.

ENTRY: 't' contains the current time.
       't0' contains the end time of model.
       'pimm' contains the proportion of immigrants who are non-UK born.
       'immid' contains the next available index number for non-UK born.
       'ukbid' contains the next available index number for UK born.
       'IMM' contains the index number for the pseudo-individual used to
        schedule external immigration events handled here.
       'ypi' contains the years per immigration, re-calculated each year
        from data on immigrants per year ('immig').

EXIT:  An immigrant is brought into the population and the next immigrant
         due is scheduled.

*/

ImmigrateG()
{ int y, n;

//-printf("Starting ImmigrateG()....\n"); fflush(stdout);
  y = (int)(t-t0);                           //Get integer year array index.

  if(Rand()<pimm[y]) { n = immid; immid++; } //Determine whether immigrant will
  else               { n = ukbid; ukbid++; } //be UK or non-UK born.

  Immigrate(n);                              //Create immigrant.

  A[IMM].pending = pImmig;                   //Schedule next immigration.
//-printf("About to schedule external immigration from ImmigrateG()\n");
  EventSchedule(IMM, t+ypi);
}



/*----------------------------------------------------------------------------*
BIRTH GENERATOR

This routine initiates a birth and schedules the next birth, at regularly spaced
intervals, acting as the peripheral event generator for births.

ENTRY: 't' contains the current time.
       'ukbid' contains the next available ID number for UK-born.
       'A[BIRTH]' is the individual designated for scheduling births.
       'ypb' contains years per birth.

EXIT:   A new individual is born.
        The next birth is scheduled.
        'ukbid' is incremented.
*/

BirthG()
{
//-printf("Starting BirthG()....\n"); fflush(stdout);
//-note Could check if t==t0 and not birth someone upon initialization at t0.
//-Produces one extra birth at initialization
  Birth(ukbid,t); ukbid += 1;            //Produce a birth and increment the next available index number for UK-born.
  A[BIRTH].pending = pBirth;             //Schedule the next birth for 'ypb'
//-printf("About to schedule external birth from BirthG()\n"); fflush(stdout);
  EventSchedule(BIRTH,t+ypb);            //years into the future.
}



/*----------------------------------------------------------------------------*
CHANGE STATES

This routine logs individuals out of compartments as they leave them and into
new compartments as they enter. It maintains counters of the numbers in each
compartment so that the array of individuals never has to be scanned for that
information.

ENTRY: 'n' indexes the individual.
       'q' contains the new state. This is a number greater than zero, in the
         range 'q0' to 'q1'.
       'A[n].state' contains the old state, either 0 or in the range 'q0' to
         'q1'. If 0, this record is not in use.
       'A[n].bstate' includes the non-susceptible states visited thus far.

EXIT:  'A[n].state' contains the number of the new state ('q' on entry).
       'A[n].bstate' incorporates the new state ('q' on entry).
       'A[n].tEntry' contains the time of entry to the new state.
       'N[u]' is decremented, where 'u' represents the old state.
       'N[v]' is incremented, where 'v' represents the new state.
*/

NewState(int n, int q)
{
//-printf("Starting NewState()...\n"); fflush(stdout);
  if(q>qU)                            //Reduce the number in the old state
    N[A[n].state] -= 1;               //unless individual is entering Uninfected, which only happens at birth or immigration.

  if(N[A[n].state]<0)                 //Make sure the state has not become
    Error2(609.0, "q=",(dec)q,        //negative.
                 " n=",(dec)n);

  A[n].state  = q;                    //Change state.
  //-A[n].tEntry = t;                 //Record the time of entry to this state

//-/*  if(q>1 && q<=8) A[n].bstate |= 1<<(q-2); //Accumulate the states this
//-  else A[n].bstate |= 1<<(q-5);    //individual has visited.
//-*/
  N[A[n].state] += 1;                 //Increase the number in the new state.
}



/*----------------------------------------------------------------------------*
TRANSFER

This routine transfer all information about an individual (including saved event
times) to a new identification number. The routine then cancels the pending
event for that index number and re-schedules it using the new index number.

ENTRY: 'n' is the new index number to be assigned, which has no event scheduled
       'n2' is the current index number of the individual.
       There is an event scheduled for 'n2'.

EXIT:  'n' is the new index number of the individual.
       The event scheduled for n2 is now re-scheduled under 'n' and all other
         data from 'n2' are transferred to 'n'. 'n2' no longer has an event
         scheduled and the index number is free to be used again.
*/

Transfer(int n, int n2)
{
//-printf("Starting Transfer()...\n"); fflush(stdout);
  if(n!=n2)
  { A[n] = A[n2]; EventRenumber(n, n2); }    //Copy data and reschedule as 'n'.
}



/*----------------------------------------------------------------------------*
ADD CUMULATIVE CASE

This routine transfer appropriate information about an individual to the next
available space in the array of cumulative cases. Note, this routine is not
actually needed until the genetic applications of the model are run.

ENTRY: 'n' is the new index number to be assigned.
       't' is the time the individual is added to the array.
       'cc' is a pointer to the file allcases.txt.
       'A[n].state' contains the disease state of the individual.
       'A[n].tBirth' contains the time of birth.
       'A[n].tImm' contains the time of immigration to UK.
       'A[n].tInfected' contains the time infection was acquired.
       'A[n].inf' contains the place infection was acquired, 0=non-UK, 1=UK.
       'A[n].sex' contains the sex of the individual.
       'A[n].rob' contains the region of birth of the individual, 0=non-UK,
        1=UK born.
       'A[n].strain' contains the strain identification number of the
        strain of infection or disease.

EXIT:  'n' is the new index number of the individual.
*/

Cumul(int n, dec t)
{
//-printf("Starting Cumul() routine...\n"); fflush(stdout);

/* Writing to file for full version of model:
* fprintf(cc, "%f\t%d\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n",
*   t, n, t-A[n].tBirth, A[n].tImm, A[n].tInfected, A[n].strain,
*   A[n].state, A[n].sex, A[n].rob, A[n].inf);
*/

/* For reduced version of model to save array space (memory), write out
without tImm & tInfected & strain of model: */
/*
* fprintf(cc, "%f\t%d\t%f\t%d\t%d\t%d\n",
*   t, n, t-A[n].tBirth,
*   A[n].state, A[n].sex, A[n].rob);
*/
}



/*----------------------------------------------------------------------------*
ADD REPORTED CASE

This routine keeps track of the number of reported cases by age, sex, place of
birth, disease site and calendar year. After a case is reported, they are
scheduled for their next event.

ENTRY: 'n' is the new index number to be assigned.
       't' is the time the individual is added to the array.
       't0' is model start time.
       't1' is the model end time.
       'A[n].tBirth' contains the time of birth.
       'A[n].sex' contains the sex of the individual.
       'A[n].rob' contains the region of birth of the individual, 0=non-UK,
        1=UK born.
       'SSAV' is the version of the model running, 0=non-SSA, 1=SSA.
       'repc' holds the numbers of reported cases.
       'A[n].tDeath contains the individual's time of death.
       'A[n].tEmigrate' contains the individual's time of emigration.
       'A[n].tExit' contains the individual's time to remote infection.
       'A[n].tMutate' contains the individual's strain mutation time.
       'A[n].ssa' holds country of birth in SSA version of model,
        0=UK & non-UK other, 1=SSA, 2=SSA,HIV+
       'A[n].tImm' contains the time of immigration to UK.
       'A[n].tInfected' contains the time infection was acquired.
       'A[n].inf' contains the place infection was acquired, 0=non-UK, 1=UK.
       'A[n].strain' contains the strain identification number of the
        strain of infection or disease.
       'A[n].state' contains the disease state of the individual.

EXIT:  'n' is the new index number of the individual.
       'Rep' contains a status code.
        1 A transmission is scheduled.
        2 Recovery is scheduled.
        3 Death is scheduled.
        4 Strain mutation is scheduled.
        5 Emigration is scheduled.
*/

Rep(int n)
{
  int s,r,y,acl,d; dec age,wt, wd, we, wr, wm, wrep;

//-printf("Starting Rep() routine...\n"); fflush(stdout);
//-/*
//-  fprintf(rc, "%f\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\n",
//-  t, n, t-A[n].tBirth, t-A[n].tDisease, A[n].tImm, A[n].tInfected,
//-  A[n].strain, A[n].state, A[n].sex,   A[n].rob, A[n].inf);
//-*/
//-/*
//-//w/o tImm & tInfected & strain:
//-fprintf(rc, "%f\t%d\t%f\t%f\t%d\t%d\t%d\n",
//-  t, n, t-A[n].tBirth, t-A[n].tDisease,
//-  A[n].state, A[n].sex,   A[n].rob);
//-*/
//-
//****More efficient reporting for E&W version of model****//
  age = t-A[n].tBirth;                       //Get age.
  if(age<15) acl=0;                          //Find age class (classes which match
  else if(age<45) acl=1;                     //notification rates).
  else if(age<65) acl=2;
  else acl=3;
  s = A[n].sex;                              //Get sex.
  r = A[n].rob;                              //Get region of birth, 0=non-UK, 1=UK.
  if(SSAV)                                    //If running 'SSA' version of model
    if(A[n].ssa)                             //check for ssa and change 'rob'
      r = 2;                                 //to '2' if SubSah African.
  y = (int)t - (int)t0;                      //Get year for array index.
  if(A[n].state>=qD4) d=0;                   //Get disease site (pulm/non-pulm)
  else d=1;                                  //for arrray index.
  repc[acl][s][r][d][y] += 1;                //Increment cases in appropriate
                                             //compartment.
  A[n].tRep = t1*2+Rand();                   //Set reporting time to time beyond
                                             //model run time so it cannot be
                                             //scheduled again, in another routine.
//-printf("Rescheduled reporting time (future) is: %f\n",A[n].tRep);
//-printf("A[n].tRep = %f\n",A[n].tRep); fflush(stdout);

  wd = A[n].tDeath;                          //Get time of death.
  we = A[n].tEmigrate;                       //Get time of emigration.
  wr = A[n].tExit;                           //Get time to remote infection.
  wm = A[n].tMutate;                         //Get strain mutation time.

  if(A[n].state<qD4)                         //If this is pulmonary disease,
  { wt = A[n].tTransm;                       //get time for transmission
    if(wt<wd && wt<we && wt<wr && wt<wm)     //and if it occurs before recovery,
    { A[n].pending = pTransm;                //mutation, emigration, and death,
//-printf("About to schedule transm from Rep()\n"); fflush(stdout);
      EventSchedule(n, wt);                  //schedule it.
      return 1; }
  }

  if(wr<wd && wr<we && wr<wm)                //If recovery will occur before
  { A[n].pending = pRemote;                  //death, emigration, and mutation,
//-printf("About to schedule remote from Rep()\n"); fflush(stdout);
    EventSchedule(n, wr);                    //schedule recovery.
    return 2; }

  if(wm<wd && wm<we)                         //If mutation will occur before
  { A[n].pending = pMutate;                  //death and emigration, schedule
//-printf("About to schedule mutate from Rep()\n"); fflush(stdout);
    EventSchedule(n, wm);                    //mutation.
    return 4; }

  if(we<wd)                                  //If emigration will occur before
  { A[n].pending = pEmigrate;                //death, schedule emigration and
//-printf("About to schedule emigration from Rep()\n"); fflush(stdout);
    EventSchedule(n,we);                     //ignore all other events.
    return 5; }

  { A[n].pending = pDeath;                   //Otherwise schedule the
//-printf("About to schedule death from Rep()\n"); fflush(stdout);
    EventSchedule(n, wd);                    //individual's death and ignore
    return 3; }                              //all other events.
}



/*----------------------------------------------------------------------------*
LIFESPAN DISTRIBUTION

This routine assigns a lifetime to an individual based on the present year, the
individual's sex and age, and other factors in the condition of the individual.
Various probability distributions may be selected. Exponentially distributed
ages, with a constant chance of death in any year, are included for calibration
with ordinary differential equation versions of this program.

ENTRY: 'sex' contains the individual's sex, 0=male, 1=female.
       'age' contains the individual's present age, years and fractions thereof.
       'mort' contains a mortality factor. For testing, this is the proportion
        of individuals who would die per year it deaths were strictly random
         (i.e., Poisson distributed in time).
       'lifedsn' defines the lifespan distribution computation:
         0 Exponential
         1 Gompertz
         2 Empirical life tables
       't' contains the present time.

EXIT:  'LifeDsn' contains the *remaining* life time (years until death) for the
        individual.

*/

static int lifedsn = 2;             //Type of longevity distribution to be used.

dec LifeDsn(int sex, dec age, dec mort)
{ int yb, y, n; dec w;

  switch(lifedsn)
  {
    case 0: return Expon(mort);              //Constant probability of death.
//  case 1: return Gompertz(sex,age,mort);   //Gompertz-Makeham death.
    case 2:                                  //Empirical life tables.
    {  yb = (int)(t-age);                    //Get year of birth.
       y  = yb-1870; if(y<0) y=0;            //Get array index for birth year.
//-printf("year of birth = %d\n",yb);
//-printf("year of birth array index = %d\n",y);
//-printf("age = %f\n",age);
//-printf("Calling RandF() from LifeDsn(), line 2056\n"); fflush(stdout);
//-printf("M1[%d][%d][0] = %f\n",y,sex,M1[y][sex][0]);
       w = RandF(A1, M1[y][sex], 122, age);
//-printf("Age is %f and number of years left is is %f\n",age, w);
       return w;
    }
    default: Error(922.0);                   //Incorrect life span selection.
  }
  return 0;                                  //(Will never reach this.)
}



/*----------------------------------------------------------------------------*
EMIGRATION TIME DISTRIBUTION

This routine assigns a time to emigration for an individual based on the present
year, the individual's sex and age, and other factors in the condition of the
individual. Exponentially distributed times, with a constant chance of death in
any year, are included for calibration with ordinary differential equation
versions of this program.

ENTRY: 'rob' contains the individual's region of birth, 0=non-UK, 1=UK
       'sex' contains the individual's sex, 0=male, 1=female.
       'age' contains the individual's present age, years.
       'em' contains the emigration rate.
       'emdsn' defines the emigration time distribution computation:
         0 Exponential
         2 Empirical migrant flow data.                                               *Improve documentation
       't' contains the present time.

EXIT:  'EmDsn' contains the *remaining* time in the UK for the                       *Improve documentation
        individual in the updated version of the program. Note that many
        individuals will have dates of emigration beyond the running time of
        the model or beyond their own death date; these individuals will never
        emigrate.
*/

static int emdsn = 0;             //Type of longevity distribution to be used.

dec EmDsn(int rob, int sex, dec age, dec em)
{ int yb, y, a, n; dec w;

  switch(emdsn)
  {
    case 0: return Expon(em);                //Constant probability of
                                             //emigration.
    case 1:                                  //(Not in use!):Use migrant flow data to find
    { //yb = (int)(t-age);                   //time to emigration. Get year
      //y  = yb-1870; if(y<0) y=0;           //of birth. Get array index
      //w = RandF(A1,em[sex][rob][y],age,t); //for birth year.
      //return w;
      return 0;
    }
    default: Error(922.0);                   //Incorrect life span selection.
  }
  return 0;                                  //(Will never reach this.)
}



/*----------------------------------------------------------------------------*
RECOVERY DISTRIBUTION

This routine assigns a time to remote infection based on the present year, the
individual's sex and age, and other factors in the condition of the individual.
Various probability distributions may be selected.

ENTRY: 'sex' contains the individual's sex, 0=male, 1=female.
       'age' contains the individual's present age, years.
       'r' contains a recovery parameter describing the individual. For
        testing this represents the proportion that would recover per
        year if recovery were strictly random (i.e., Poisson distributed in
        time).
       't' contains the present time.

EXIT:  'RecovDsn' contains the time until recovery, in years.
*/

static dec recovdsn = 0;     //Type of recovery distribution to be used.
static dec rmu      = 0.0;   //Centering of recovery distribution, years.
static dec rsigma   = 0.1;   //Half-width of recovery distribution, years.

dec RecovDsn(int s, dec age, dec r)
{ dec w;

  switch((int)recovdsn)                         //Select the type of recovery.
  { case 0: return Expon(r);                    //(completely random)
    case 1: w = 0;                       break; //(completely fixed)
    case 2: w = Uniform(-rsigma,rsigma); break; //(uniform variation)
    case 3: w = LogNormal(rmu, rsigma);  break; //(log-normal variation)
    case 4: w = Gauss(.0, rsigma);       break; //(truncated Gaussian variation)
    case 5: w = Cauchy(.0, rsigma);      break; //(truncated Cauchy variation)
    default: Error(922.);
  }

  return max(1e-9, w+1./r);                     //Return the time for recovery,
}                                               //always after a slight delay.



/*----------------------------------------------------------------------------*
TIME TO DISEASE

This routine assigns a time to disease based on the individual's age, sex,
region of birth, infection status (Recent Infection, Remote Infection,
Reinfection), and -- in the SSAV version of the model -- HIV status. Note,
Recent Infection and Reinfection states are handled similarly, whereas Remote
Infection is handled somewhat differently.

Notes on SSA version of model: Would like to get disease rates correct so things
are comparable between the SSA and non-SSA models. I think it is best to leave
them as a ratio ('ehiv') and then when fitting model, the non-SSA model should
fit higher disease risk for non-UK born than in the SSA model. In the SSA model
they should be lower since a portion of non-UK born individuals will have
elevated risk due to HIV.

ENTRY: 'n' contains the individual's identifier.
       's' contains the individual's sex, 0=male, 1=female.
       'a' contains the individual's integer age.
       'rob' contains the individual's region of birth, 0=non-UK, 1=UK, 2=SSA
       'tinf' contains the time since infection (years).
       'd1' and 'd3' contain the probability of progressing to disease for 'qI1'
         and 'qI3' over the first five years of infection.
       'd2' contains the probability per year of progressing to disease for
         'qI2'.
       'A[n].state' contains the present state of the individual, 'qI1', 'qI2',
         'qI3'.
       't' contains the present time.
       'drr' gives the relative risk of disease over the first five years
         of infection.
       'B1' contains the year of infection, for use with 'drr'.
       'SSAV' contains 1 if SSA version of model is to be run.

EXIT:  'Tdis' contains the number of years until disease development.
*/

dec Tdis(int n, int a, int s, int rob, dec tinf)
{ dec d,age,w;
//-printf("Starting Tdis() routine with calls to RandF(), line 2179 starts\n");

  if(SSAV && A[n].ssa==2)              //If running SSA version of model and
    rob=2;                             //individual is HIV+, adjust 'rob'.

  switch(A[n].state)
  { case qI1:                          //Process Recently Infection.
    { d = d1[s][rob][a]                //First calculate overall disease risk
        *(1-Val(1,tinf,B1,drr,0,5));   //for the first five years, correcting
                                       //for 'tinf'>0 if applicable.

      if(Rand()>d)                     //If disease should NOT occur, schedule
      { w = 2*RT+Rand();               //disease past the running time of
        //-printf("Disease (primary) should NOT happen -- time until is %f years\n",w);
        return w; }                    //model so that it never occurs.
      else
      { w = RandF(B1,drr,6,tinf);      //If disease should occur, randomly
        //-printf("Disease (primary) should happen in %f years\n",w);
        return w; } }                  //choose year, based on relative risk
                                       //over five years.

    case qI3:                          //Process Reinfection.
    { d = d3[s][rob][a]                //First calculate overall disease risk
        *(1-Val(1,tinf,B1,drr,0,5));   //for the first five years, correcting
                                       //for 'tinf'>0 if applicable.

      if(Rand()>d)                     //If disease should NOT occur, schedule
      { w = 2*RT+Rand();               //disease past the running time of
        return w; }                    //model so that it never occurrs.
      else
      { w = RandF(B1,drr,6,tinf);      //If disease should occur, get year
        //-printf("Disease (reinf) should happen in %f years\n",w);
        return w; } }                  //from relative risk over five years.

    case qI2:                          //Process Remote Infection.
    { age = t-A[n].tBirth;
      w = RandF(A2,d2[s][rob],AC+2,age);
//-printf("Time until react disease is %f years, s=%d, rob=%d, age=%f\n",w,s,rob,age);
      return w;      }

    default: Error(922.);   }
  return 0;                            //(Will never reach this.)
}


/*----------------------------------------------------------------------------*
GET RANDOM AGE FOR IMMIGRANT

This routine assigns an age to an individual who is immigrating into the population.
Age is randomly assigned based on probabilities of age classes from data.
Probabilities of age classes are conditional on sex and region
of birth.

ENTRY: 'n' contains the individual's identifier.
       's' contains the individual's sex, 0=male 1=female.
       'r' contains the individual's region of birth, 0=non-UK, 1=UK
       't' contains the present time.
       'immage' contains cumulative probabilities of the age classes; to
         be compatible with future calls to RandF(), the first cumulative
         probability is 0.
       'SSAV' contains 1 if SSA version of model is to be run.

EXIT:  'GetAge' contains the age (in years) of the individual.
*/

dec GetAge(int n, int s, int r)
{ dec rn,age; int y;

  rn = Rand();                     //Get random number.

  if(SSAV)                         //If running SSA version of model,
    if(A[n].ssa) r=2;              //check if SSA and adjust 'r'(rob) if so.

  y = (int)t - (int)t0;            //Get array index for calendar year.


  if(rn<immage[y][s][r][1])        //Assign a random age class within the
    return Rand()*15;              //correct age class, depending on the
  if(rn<immage[y][s][r][2])        //random number draw and cumulative
    return Rand()*10+15;           //probabilities of each age class specified
  if(rn<immage[y][s][r][3])        //by 'immage'.
    return Rand()*10+25;
  if(rn<immage[y][s][r][4])
    return Rand()*10+35;
  if(rn<immage[y][s][r][5])
    return Rand()*15+45;
  age=Expon(0.10)+60;              //For the age class 60+, add age of 60 plus
  if(age>=121) age=120+Rand();     //draw from exponential distribution with
                                   //mean of 10 years.
  //-printf("age from 60+ = %f\n",age);

  return age;
}



/*----------------------------------------------------------------------------*
FILE INITIALIZATION

This function writes headers and stores documentation for output files.
Note: These headers need updating, as output have/do change periodically.

ENTRY: 'cc' contains a pointer to the file "allcases.txt".
       'rc' contains a pointer to the file "repcases.txt".

EXIT:  Headers are written to above output files.
*/

Init()
{
//-printf("Starting file initialization routine...\n"); fflush(stdout);
/*
  fprintf(cc,"Time disease progression\tIndividual\tAge\tTime of immigration\t\
  Time of most recent infection\tStrain\tDisease state\tSex\t\
  Region of birth\tRegion infection acquired\n");

  fprintf(rc, "Time of report\tIndividual\tAge\tDisease duration\t\
  Time of immigration\tTime of infection\tStrain\tDisease state\t\
  Sex\tRegion of birth\tRegion infection acquired\n");
*/

/* Documentation for cumulative cases:
  t                                  Time of progression to disease.
  n                                  Index number for individual
  t-A[n].tBirth                      Age
  A[n].tImm                          Time of immigration to UK.
  A[n].tInfected                     Time individual most recently infected.
  A[n].strain                        Strain type identification number
  A[n].state                         Number of present state
  A[n].sex                           Sex of this individual
  A[n].rob                           Region of birth
  A[n].inf                           Region infection acquired
*/

/* Documentation for reported cases:
  t                //Time of report.
  n                //Index number for individual
  t-A[n].tBirth    //Age
  t-A[n].tDisease  //Disease duration
  A[n].tImm        //Time of immigration to UK.
  A[n].tInfected   //Time individual was infected.
  A[n].strain      //Strain type identification number
  A[n].state       //Disease state
  A[n].sex         //Sex
  A[n].rob         //Region of birth
  A[n].inf         //Region infection acquired
*/

}



/*----------------------------------------------------------------------------*
DATA PROCESSING AND ARRAY INITIALIZATION

This function reads data files into appropriate arrays and also copies smaller
amounts of data, coded here instead of read from files, into appropriate
arrays.

ENTRY: Files are set up for 'FileIO', including:
       'births.txt' contains the number of births.
       'mort.txt' is a file of mortality data (life tables).
       'casefat.txt' holds case fatality rates.

EXIT:  'A1' contains ages which accompany 'M1' for calls to 'RandF'
       'A2' contains ages which accompany 'd2' for calls to 'RandF'
       'B1' contains times (years) which accompany 'drr' for calls to 'RandF'.
       'Ax' contains disease state variables which accompany 'infimm'.
       'drr' is filled with cumulative relative risks of developing disease
         over first five years since infection or reinfection.
       'em' contains emigration rates.
       'p1', 'p2', and 'p3' contain the proportion of disease pulmonary for
         Primary Disease, Reactivation Disease and Reinfection Disease.
       'sdf1', 'sdf2', 'sdf3' are filled with sex factors of disease
       risk for primary, reactivation, and reinfection disease respectively,
       if applicable.
       'duk1p', 'duk2p', 'duk3p' are filled with values for pulmonary disease
       progression rates (by age and sex) for UK-born (if applicable).
       'd1p', 'd2p', 'd3p' are filled with values for pulmonary disease
       progression by age, sex, and rob (non-UK born still not filled in).
       All data read in through FileIO() is stored properly, including:
        arrays 'bcy', 'mort', 'cft' contain data for births, mortality and case
        fatality due to TB.
*/

Data()
{ int i,j,k,temp,a,s,y,r,ac,st; dec temp2;
//-printf("Entering Data() routine....\n");

  for(i=0; i<AC; i++)                //Create array of ages which correspond
  { A1[i]=i;                         //to cumulative probabilities of death
    A2[i]=i; }                       //in M1 (A1) and cumulative probabilities
  A2[122] = 3000;                    //of progression for Remote Infection
  A2[123] = 3001;                    //in d2[][][] (A2) for calls to 'RandF'.
                                     //In A2, 3000,3001 refer
                                     //to years/ages far enough into the future
                                     //that they will not happen.

  for(i=0; i<9; i++)                 //Holds disease state variable, which
    Ax[i] = i+1;                     //which accompanies 'infimm'.

  for(i=0; i<6; i++)                 //Create array of times which correspond
    B1[i]= i;                        //to cumulative times to disease for first
                                     //five years of infection (drr[]).

  drr[0] = 0.0;                      //Fill in values for cumulative,
  drr[1] = 0.604594921;              //relative risk of disease progression
  drr[2] = 0.852478839;              //by year since infection, for the first
  drr[3] = 0.931076179;              //five years. These represent cumulative
  drr[4] = 0.983071342;              //risk before index year. Relative risks
  drr[5] = 1.0;                      //from Vynnycky & Fine 1997, though
                                     //modified to be cumulative risks.

  em[F][UK] =0.00225;                //Set emigration rates from ONS data,
  em[M][UK] =0.00280;                //averaged over 1981-2009.
  em[M][NUK]=0.02888;
  em[F][NUK]=0.02656;
  em[M][SSA]=0.02009;
  em[F][SSA]=0.01528;

  for(a=0; a<121; a++)               //Fill arrays for the proportion
  { p1[a][M][NUK] = 0.528236447;     //of disease which is pulmonary
    p1[a][M][UK]  = 0.468333833;     //by age, sex & rob, with sex and
    p1[a][F][NUK] = 0.740686033;     //rob dependent values taken from
    p1[a][F][UK]  = 0.672633119;     //notification data.
    p2[a][M][NUK] = 0.528236447;
    p2[a][M][UK]  = 0.468333833;
    p2[a][F][NUK] = 0.740686033;
    p2[a][F][UK]  = 0.672633119;
    p3[a][M][NUK] = 0.528236447;
    p3[a][M][UK]  = 0.468333833;
    p3[a][F][NUK] = 0.740686033;
    p3[a][F][UK]  = 0.672633119; }

  if(DPARAM)                         //For new model version with expanded set
  { sdf1[0]   = 1.0;                 //of variable disease progression param-
    sdf2[0]   = 1.0;                 //eters, set ratios of female:male disease
    sdf3[0]   = 1.0;                 //progression (male rates are multiplied
    sdf1[1]   = 1.0;                 //by these to get female disease
    sdf2[1]   = 0.000048/0.000299;   //progression risk/rates).
    sdf3[1]   = 0.0001/0.0825;

    d1uk10[M] = 0.0406;              //Also for new version of model, set
    d2uk10[M] = 0.000000000982;      //default values for disease progression
    d3uk10[M] = 0.0689;              //in males to match values from old
    d1uk20[M] = 0.138;               //version, taken from Vynnycky & Fine.
    d2uk20[M] = 0.000299;
    d3uk20[M] = 0.0825;   }

  if(!DPARAM)                        //For old model version:
  { duk1p[0][M] = 0.0406;            //These are estimated risks of disease
    duk2p[0][M] = 0.000000000982;    //progression for respiratory disease,
    duk3p[0][M] = 0.0689;            //by age (0-10, 20+) and sex (M/F), from
    duk1p[1][M] = 0.138;             //Vynnycky & Fine 1997. For duk1 & duk3,
    duk2p[1][M] = 0.000299;          //they are cumulative risks over the
    duk3p[1][M] = 0.0825;            //first five years since infection/
                                     //reinfection. For duk2, they are annual
    duk1p[0][F] = duk1p[0][M];       //risks of infection for that age. d1 =
    duk2p[0][F] = duk2p[0][M];       //Recent, d2=Remote, and d3=Reinfection.
    duk3p[0][F] = duk3p[0][M];       //Note, primary disease risks and risks
    duk1p[1][F] = duk1p[1][M];       //for ages 0-10 are equal for male and
    duk2p[1][F] = 0.000048;          //female.
    duk3p[1][F] = 0.0001;

    for(a=0; a<10; a++)              //For UK-born aged 0-10, 10-20, & 20+,
    for(s=0; s<2;  s++)              //fill array with risks/rates of disease
    { d1p[a][s][UK] = duk1p[0][s];   //progression by age (single year class),
      d2p[a][s][UK] = duk2p[0][s];   //sex, & rob, still RESPIRATORY only,and
      d3p[a][s][UK] = duk3p[0][s]; } //cumulative risk over first 5 yrs for
                                     //d1p,d3p (infected at age 'a'); annual
    for(a=10; a<20; a++)             //risks by age for d2p.
    for(s=0;  s<2;  s++)
    { d1p[a][s][UK] = duk1p[0][s] + (a-10)*((duk1p[1][s]-duk1p[0][s])/10);
      d2p[a][s][UK] = duk2p[0][s] + (a-10)*((duk2p[1][s]-duk2p[0][s])/10);
      d3p[a][s][UK] = duk3p[0][s] + (a-10)*((duk3p[1][s]-duk3p[0][s])/10); }

    for(a=20; a<121; a++)
    for(s=0;  s<2;   s++)
    { d1p[a][s][UK] = duk1p[1][s];
      d2p[a][s][UK] = duk2p[1][s];
      d3p[a][s][UK] = duk3p[1][s]; } }

  #define RSCALE (SUPER? "r|": "r|=n/5")     //Scaling factor 1/5 for laptops.
  FileIO("births.txt",   fmt[0],  RSCALE);   //Read birth data.
  FileIO("immigs.txt",   fmt[1],  RSCALE);   //Read immigration data.
  FileIO("pimm.txt",     fmt[2],  "r|");     //Read immigrants non-UK born.
  FileIO("ssaim.txt",    fmt[3],  "r|");     //Read immigrants SSA-born.
  FileIO("propmale.txt", fmt[4],  "r|");     //Fraction of births that are male.
  FileIO("hivp.txt",     fmt[5],  "r|");     //Read HIV prevalence.
  FileIO(SSAV? "infimm1.txt":"infimm0.txt",  //Read probabilities of immigrants
         fmt[6],"r|");                       //entering the different disease
                                             //states, 'infimm'.
  FileIO("inf1981.txt", fmt[19], "r|");      //Read disease state probabilities
                                             //for 1981, 'inf1981'.
  FileIO("ssa1981.txt", fmt[9],  "r|");      //Read the proportion of non-UK
                                             //born who are SSA, 'ssa1981'.

  FileIO("n1981.txt",    fmt[10], RSCALE);   //Read 'n1981' as integers.
  FileIO("mort.txt",     fmt[14], "r|");     //Mortality data.
  FileIO("casefat.txt",  fmt[15], "r|");     //Case fatality data.
  FileIO("smear.txt",    fmt[20], "r|");     //Read proportions smear-positive.
  FileIO("N3.txt",       fmt[21], "r|");     //Read actual population sizes.
  FileIO(SSAV? "immsex1.txt":"immsex0.txt",  //Read immigrants who are male.
         fmt[11], "r|");
  FileIO(SSAV? "immage1.txt":"immage0.txt",  //Read immigrants by age class.
         fmt[12], "r|");

  for(i=0; i<BY; i++)                        //Audit the cumulative mortality
  for(s=0; s<2;  s++)                        //data to make sure each table
    monotone(M1[i][s], AC, 1, i, s);         //increases from 0 to 1.

  for(i=0; i<RT; i++)                        //Convert 'immageX' to 'immage',
  for(s=0; s<2;  s++)                        //the array of cumulative
  for(r=0; r<(2+SSAV); r++)                  //probabilities of being in an
    immage[i][s][r][1]=immageX[i][s][r][0];  //age class.

  for(i=0; i<RT; i++)
  for(s=0; s<2;  s++)
  for(r=0; r<(2+SSAV); r++)
  for(a=2; a<7;  a++)
    immage[i][s][r][a]=immageX[i][s][r][a-1]+immage[i][s][r][a-1];

  for(i=0; i<RT; i++)                        //Make sure the last age class
  for(s=0; s<2;  s++)                        //has probability 1, first is 0
  for(r=0; r<(2+SSAV); r++)                  //so this is compatible with 'RandF'
  { immage[i][s][r][6]=1.0;                  //format in case that is employed.
    immage[i][s][r][0]=0.0; }

//-for(a=0; a<121; a++)                          //Set upper end of disease
//-for(r=0; r<2;  r++)                           //state probabilities to 1
//- inf1981[a][r][8] = 1.0;                      //for old way of reading.
}



/*----------------------------------------------------------------------------*
PARAMETER CHANGING

This function updates variables associated with parameters that change with
each model run. This routine must come after any change to parameters, e.g.
through the 'gparam' function. Currently four disease risk parameters are
varied during model fitting: 'd1uk20[M]', 'd2uk20[M]', 'd3uk20[M]' and 'df'.
'd1uk20[M]' is the risk of developing Primary Disease for UK-born males aged
20 years above. 'd2uk20[M]' is the annual risk of developing Reactivation
Disease for UK-born males aged 20 years and above. 'd3uk20[M]' is the risk
of developing Reinfection Disease for UK-born males aged 20 years and above.
'df' is the factor by which UK-born disease risks are multiplied to get
non-UK born risks.

ENTRY: 'drr' contains relative cumulative rates of disease progression by
         time since infection, up to five years.
       'd1uk20[M]' contains the risk of developing Primary Disease for UK-born
         males aged 20 years above
       'd2uk20[M]' contains the annual risk of developing Reactivation Disease
         for UK-born males aged 20 years and above.
       'd3uk20[M]' contains the risk of developing Reinfection Disease for
         UK-born males aged 20 years and above.
       'df' contains the factor by which UK disease risk are multiplied to
         get non-UK born disease risks.
       'sdf1' contains the sex ratios of female:male disease risk (Primary Disease).
       'sdf2' contains the sex ratios of female:male disease risk (Reactivation Disease).
       'sdf3' contains the sex ratios of female:male disease risk (Reinfection Disease).
       'presp' contains the proportion of disease respiratory for children.

EXIT:  'd1'contains Primary Disease risk by sex, rob (UK/nonUK/HIV) and age.
       'd2'contains Reactivation Disease rate by sex, rob (UK/nonUK/HIV) and age.
       'd3'contains Primary Disease risk by sex, rob (UK/nonUK/HIV) and age.

*/

/*
Notes: Disease progression risks.

Disease progression risks are specified separately for Recent, Remote, and
Reinfected individuals (d1,d2,d3), and also allowed to vary by sex, rob and age.
Following Vynnycky & Fine, the age-dependency of risk is specified by two
parameters -- one for risk ages 0-10 (constant) and one for 20+ (constant). Risk
between ages 10-20 is assumed to increase linearly from the rate at 10 to the
rate at age 20.For those over 10 and under 20
overall disease progression rate is: A0+ (age-10)*((A20-A10)/10)

For d1 and d3, these represent
overall, cumulative risks of disease progression in the first five years of
infection or reinfection for infection at a given age. The array 'drr' specifies
the cumulative relative risk over these five years and is used along with d1/d2
to generate a time to disease, if applicable. For d2 the disease progression
risks are annual rates at a given age (and sex/rob) for disease progression. d2
is converted to cumulative risk by age, and the cumulative distribution is used,
along with current age of individual infected, to assign a time to disease.

Disease progression risks estimated by Vynnycky & Fine 1997 for white ethnic
males are used to set UK-born males and UK-born female risk for those under ten
years of age. For all ages, female risks are calculated by multiplying male
risks by the risk ratios for sex, 'sdf1', 'sdf2', and 'sdf3'.  Risks for 20
year-olds are allowed to vary in the model (one for each disease types, so
three parameters). For non-UK born risks, 'df' is multiplied
by the UK-born risks to get non-UK born risks. 'df' is allowed to vary in the
model. HIV-positive risks are further multiplied by 'ehiv'. This means as the
three UK-born risks and 'df' are varied during model fitting, non-UK born
disease progression risks would need to be re-generated each time the model
is run.

**Note, one complication regards pulmonary versus overall disease risk. In
Vynnycky and Fine 1997, risks are for respiratory disease. However, disease risk
needed for the model is overall disease risk, pulmonary plus non-pulmonary.
For those fixed in the model (UK-born under 10), the respiratory risk is
corrected to equal overall (pulmonary + non-pulmonary) risk.

*/

Param()
{ int a,s,r;

//-printf("Entering Param() and df is: %f\n",df);

  dec ep = 0.00000000000001;       //Check that 'ehiv' and 'df' are not
  if(ehiv<ep) ehiv= ep;            //negative, making them 'ep', a small
  if(df<ep)   df  = ep;            //positive number if not.

  if(DPARAM)                       //For new way of varying disease risks:
  { if(d1uk10[M]<ep) d1uk10[M]=ep; //First check that all risks/rates are
    if(d1uk20[M]<ep) d1uk20[M]=ep; //positive (greater than a small positive
    if(d2uk10[M]<ep) d2uk10[M]=ep; //number).
    if(d2uk20[M]<ep) d2uk20[M]=ep;
    if(d3uk10[M]<ep) d3uk10[M]=ep;
    if(d3uk20[M]<ep) d3uk20[M]=ep;

    d1uk10[F] = d1uk10[M]*sdf1[0];  //Also for new version of model, set
    d2uk10[F] = d2uk10[M]*sdf2[0];  //females' values for disease progression,
    d3uk10[F] = d3uk10[M]*sdf3[0];  //relative to males, with sex ratios
    d1uk20[F] = d1uk20[M]*sdf1[1];  //calculated from rates in Vynnycky &
    d2uk20[F] = d2uk20[M]*sdf2[1];  //Fine.
    d3uk20[F] = d3uk20[M]*sdf3[1];

    for(s=0; s<2; s++)
    { d1uk10[s] = d1uk10[s]/presp; //Also for new version of model, set
      d2uk10[s] = d2uk10[s]/presp; //divide by 'presp' to get overall (not
      d3uk10[s] = d3uk10[s]/presp;  //respiratory) progression rates/risks.
      //d1uk20[s] = d1uk20[s]/presp; //Only doing this for the fixed rates,
      //d2uk20[s] = d2uk20[s]/presp; //in those aged 0-10, which are based
      //d3uk20[s] = d3uk20[s]/presp; //on Vynn & Fine (respiratory) rates.
    }                               //Others are for TESTING ONLY!

    for(a=0; a<10; a++)             //Now expand rates to all ages for
    for(s=0; s<2;  s++)             //by assuming constant risk/rate from age
    { d1[s][UK][a] = d1uk10[s];     //0-10, linear increase from 10-20, and
      d2[s][UK][a] = d2uk10[s];     //constant risk/rate from age 20+. Note,
      d3[s][UK][a] = d3uk10[s]; }   //for d1 & d3 these are cumulative risks
                                    //over first 5 yrs for (infected at age 'a')
    for(a=10; a<20; a++)            //while for d2 these are annual rates of
    for(s=0;  s<2;  s++)            //progression.
    { d1[s][UK][a] = d1uk10[s] + (a-10)*((d1uk20[s]-d1uk10[s])/10);
      d2[s][UK][a] = d2uk10[s] + (a-10)*((d2uk20[s]-d2uk10[s])/10);
      d3[s][UK][a] = d3uk10[s] + (a-10)*((d3uk20[s]-d3uk10[s])/10); }

    for(a=20; a<121; a++)
    for(s=0;  s<2;   s++)
    { d1[s][UK][a] = d1uk20[s];
      d2[s][UK][a] = d2uk20[s];
      d3[s][UK][a] = d3uk20[s]; }
  }

  else                                    //For old version of disease risk--
  { for(a=0; a<121; a++)                  //For UK-born get overall disease risks
    for(s=0; s<2;   s++)                  //(d1&d3) & annual risks (d2) from
    { d1[s][UK][a] = d1p[a][s][UK]/presp; //respiratory risks, taken from
      d2[s][UK][a] = d2p[a][s][UK]/presp; //Vynn & Fine 1997 (divided by prop.
      d3[s][UK][a] = d3p[a][s][UK]/presp; } //of all TB which is respiratory
  }                                       //(see Data()). Indexed by sex, rob,
                                          //and age.

  for(a=0; a<121; a++)                 //For non-UK born: multiply UK-born
  for(s=0; s<2;   s++)                 //risks by 'df' to get non-UK born
  { d1[s][NUK][a] = df*d1[s][UK][a];   //disease risks by sex, rob & age.
    d2[s][NUK][a] = df*d2[s][UK][a];
    d3[s][NUK][a] = df*d3[s][UK][a]; }

  for(a=0; a<121; a++)                    //Check that overall disease risks
  for(s=0; s<2;   s++)                    //(d1, d3) are not above 1; also check
  { if(d1[s][NUK][a]>1) d1[s][NUK][a] = 1; //annual rates (d2) are not above 1
    if(d2[s][NUK][a]>1) d2[s][NUK][a] = 1; //for non-UK born.
    if(d3[s][NUK][a]>1) d3[s][NUK][a] = 1; }

  if(SSAV)
  { for(a=0; a<121; a++)                  //Get disease risks for HIV-positive
    for(s=0; s<2;   s++)                  //SSAs by multiplying non-UK born
    { d1[s][HIV][a] = ehiv*d1[s][NUK][a]; //risks/rates by factor 'ehiv'.
      d2[s][HIV][a] = ehiv*d2[s][NUK][a];
      d3[s][HIV][a] = ehiv*d3[s][NUK][a]; }

    for(a=0; a<121; a++)                  //Make sure disease risks for HIV+
    for(s=0; s<2;   s++)                  //are not above 1.
    { if(d1[s][HIV][a]>1) d1[s][HIV][a]=1;
      if(d2[s][HIV][a]>1) d2[s][HIV][a]=1;
      if(d3[s][HIV][a]>1) d3[s][HIV][a]=1; }
  }

  for(s=0; s<2; s++)                      //Fix end of array for d2 (longer).
  for(r=0; r<3; r++)
    d2[s][r][121] = d2[s][r][120];

  for(s=0; s<2; s++)               //Redefine 'd2' as cumulative probability of
  for(r=0; r<3; r++)               //disease progression; first translate risk
    d2[s][r][1] = d2[s][r][0];     //in first year of life (d2[s][r][0]) as
                                   //cumulative risk experienced before age 1
  for(a=2; a<=121; a++)            //(d2[s][r][1]). Then convert rest of array
  for(s=0; s<2; s++)               //to cumulative risks. Indexed by sex, rob
  for(r=0; r<3; r++)               //and age.
    d2[s][r][a] = d2[s][r][a-1]+(1-d2[s][r][a-1])*d2[s][r][a];

  for(s=0; s<2; s++)               //Make sure cumulative probability of disease
  for(r=0; r<3; r++)               //has not gone beyond 1.
  { if(d2[s][r][121]>1) Error(754.1); }

  for(s=0; s<2; s++)               //Finish cumulative distribution so that
  for(r=0; r<3; r++)               //cumulative risk before age 0 is 0 and
  { d2[s][r][0] = 0.0;             //those who should never progress to disease
    d2[s][r][122] = d2[s][r][121]; //are assigned times far into the future so
    d2[s][r][123] = 1.0; }         //that disease progression does not happen.

/*
*int i,j;
*printf("Check for zero in d2\n");
* for(i=0; i<2; i++)
* for(j=0; j<3; j++)
*   printf("%lf\n",d2[i][j][0]);
*
*printf("Check for 1 in (123) d2\n");
* for(i=0; i<2; i++)
* for(j=0; j<3; j++)
*   printf("d2[sex=%d][rob=%d]: %lf\n",i,j,d2[i][j][123]);
*/

}



/*----------------------------------------------------------------------------*
INITIALIZE STARTING POPULATION

This function sets up the intial population by looping through matrix 'n1981'
which holds numbers in the initial population by age, sex, and rob. For the SSA
version of model, 'ssa1981' is used for the proportions of SubSaharan Africans
among non-UK born in 1981 by age and sex.

Notes: In Param(), could multiply 1981 (if they are proportions) by initial pop
size, e.g. 'initpop'. Could make assignments deterministic since the numbers are
so large -- this would take some planning for dealing with remainders etc.

ENTRY: 'n1981' contains numbers by age, sex, and rob for individuals in the
         population at initialization in 1981.
       'ssa1981' contains the proportion of SSAs among non-UK born
         at population initialization.

EXIT:  The intial population is set up; each individual is assigned attributes
        and scheduled for exactly one event (other event times may be stored
        for an individual.
*/

InitPop()
{ int a,s,i,n,st,rob; dec age,wd,we,wv,tinf;
  ukbid = maximm+1;                          //Initialize ID numbers to
  immid = 1;                                 //correct values.

//-printf("About to start InitPop() with calls to RandF()\n");
/* NOTE : Could make function so that UK-born, non-UK born and SSA-born are
initialized with one function, called 3 times......also function so that
states are assigned and processed in a function */

  for(a=0; a<121; a++)                       //First, initialize UK-born
  for(s=0; s<2;   s++)                       //(rob=1) population for all age
  for(i=0; i<n1981[a][s][UK]; i++)            //and sex categories.
  { n = ukbid; ukbid++;                      //Take the next available ID.
    age = a+Rand();                          //Assign age plus random bit.
    A[n].tBirth = t-age;                     //Assign birth time from age.
    A[n].sex = s;                            //Assign sex.
    A[n].rob = rob = UK;                     //Set to UK-born.

    BasicInd(n,UK,age,s);                    //Set up basic individual.
    DisState(n,UK,a); }                      //Assign disease state and
                                             //process accordingly.

  if(SSAV)                                   //Process non-UK born in SSA
    for(a=0; a<121; a++)                     //version of model, taking into
    for(s=0; s<2;   s++)                     //account the proportion of
    for(i=0; i<n1981[a][s][NUK]; i++)        //SSAs and their HIV
    { n = immid; immid++;                    //status.
      age = a+Rand();
      A[n].tBirth = t-age;
      A[n].sex = s;
      A[n].rob = rob = NUK;
      if(Rand()<ssa1981[a][s])               //If SubSaharan African, indicate
      { A[n].ssa = 1;                        //this and assign HIV status.
        rob=SSA;
        if(Rand()<hivp[s][0]) A[n].ssa=2; }

      BasicInd(n,NUK,age,s);                 //Set up basic individual.
      DisState(n,rob,a); }                   //Assign disease state and
                                             //process accordingly.

  else                                       //Process non-UK born for non-SSA
    for(a=0; a<121; a++)                     //version of model.
    for(s=0; s<2; s++)
    for(i=0; i<n1981[a][s][NUK]; i++)
    { n=immid; immid++;
      age=a+Rand();
      A[n].tBirth = t-age;
      A[n].sex = s;
      A[n].rob = rob = NUK;

      BasicInd(n,NUK,age,s);                  //Set up basic individual.
      DisState(n,NUK,a); }                    //Assign disease state and
                                             //process accordingly.
}



/*----------------------------------------------------------------------------*
SET UP BASIC INDIVIDUAL FOR POPULATION INITIALIZATION

This function sets up a basic individual assigned to a death, emigration, or
vaccination time. This is similar to birth but processes individuals of any
age, sex, or region of birth (anyone being initialized when the model starts).
The scheduled event may change if a state other than 'Uninfected' is
assigned when diseaes state is assigned. This function is merely used to
get the individual initialized.

ENTRY:  'A[n]' individual is not scheduled for any event.
        'rob' is 0 for non-UK born (ONUK and SSA), 1 for UK-born.
        'age' is the age of the individual in years.
        'sex' is 0 for males and 1 for females.

EXIT:   'A[n]' is in the Uninfected state and scheduled for its earliest
          event.

*/

BasicInd (int n, int rob, dec age, int s)
{ dec wd, we, wv;

  NewState(n,qU);                               //Assign to Uninfected state.
  A[n].tDeath = wd = t+LifeDsn(s,age,m1[0][0]); //Assign time of death.
  if(wd<A[n].tBirth+age) Error(612.2);          //Check death time.
//- A[n].tEmigrate = we = t+Expon(em[s][1]);   //Assign time of emigration (old).
  A[n].tEmigrate = we                           //Assign time of emigration.
                 = t+EmDsn(rob,s,age,em[s][rob]);
  if(age<v3[rob] && Rand()<v1[rob]*v2[rob])     //Calculate time to vaccination,
    wv = t+(v3[rob]-age)+Rand();                //set to time which never
  else wv = t+2*RT+Rand();                      //happens if it should not occur.

  if(wv<wd && wv<we)                            //If vaccination is the earliest
  { A[n].pending = pVaccin;                     //event, schedule it.
    EventSchedule(n, wv); }

  else if(wd<we)                                //If death is the earliest
  { A[n].tExit = wd;                            //event, schedule it.
    A[n].pending = pDeath;
//-printf("About to schedule death from InitPop()\n"); fflush(stdout);
    EventSchedule(n, wd); }
  else                                          //Or, if emigration is the
  { A[n].tExit = we;                            //earliest event, schedule that.
    A[n].pending = pEmigrate;
//-printf("About to schedule emig from InitPop()\n"); fflush(stdout);
    EventSchedule(n, we); }
}



/*----------------------------------------------------------------------------*
ASSIGN DISEASE STATE FOR INITIAL MEMBER OF POPULATION
This function assigns one of the 8 disease states (8 instead of 11 because
those with disease are not assigned to pulmonary or non-pulmonary until after
being processed).

ENTRY:  'A[n]' is set up as Uninfected.
        'Ax' contains the disease state, used for RandF() along with 'inf1981'.
        'inf1981' contains the probabilities of the different disease states,
          which are numbers in 'Ax'.

EXIT:   'A[n]' has been assigned disease state (may remain unchanged) and
          processed accordingly.

*/

DisState(int n, int rob, int a)
{ int st,str; dec tinf;

//-printf("Calling RandF() from line 2870, InitPop()\n"); fflush(stdout);
//  st = 1+(int)RandF(Ax,inf1981[a][A[n].sex][rob],9,1.0); //Get disease state.
    st = 1+(int)RandF(Ax,inf1981[a][A[n].sex][rob],9,1.0); //Get disease state.
//-printf("Disease state is %d\n",st); fflush(stdout);
    switch(st)
    {
    case 2:                                     //Send to vaccination.
      EventCancel(n);
      Vaccination(n);
      break;

    case 3:                                  //Get the time of infection before
      tinf = Rand()*5;                       //sending to 'Infect'.
      str=StrainNum(rob);
      Infect(n,tinf,str);
      break;

    case 4:                                  //Before sending to 'Remote', put
      EventCancel(n);                        //in temporary disease state to
      NewState(n,qD1);                       //facilitate re-scheduling of
      //A[n].strain = StrainNum(rob);        //events in 'Remote' function.
      Remote(n);
      break;

    case 5:                                  //Before sending to 'Infect', put
      NewState(n,qI2);                       //in "Remote infection" state so
      tinf = Rand()*5;                       //that this is correctly treated
      str=StrainNum(rob);                    //as a reinfection.
      Infect(n,tinf,str);
      break;

    case 6: case 7: case 8:                  //Process disease states by first
      EventCancel(n);                        //putting in correct infection
      NewState(n,st-3);                      //state and then sending to
      //A[n].strain=StrainNum(rob);          //'Disease'.
      Disease(n);
      break;

    default:
      if(st!=1) Error(618.2);
    }
}



/*----------------------------------------------------------------------------*
CHOOSE STRAIN IDENTIFICATION NUMBER

Future versions of the model will track strain type profile identifiers for
all infections. This is a stub for a future routine to choose the strain
number at model initialization or for immigrants at immigration. It is a
placeholder and entry/exit conditions have not been determined yet.

ENTRY:

EXIT:

*/

int StrainNum(int rob)
{
  return 0;
}


/*----------------------------------------------------------------------------*
REPORTING

This function is called periodically to display the cumulative number of
infections and other statistics. THIS FUNCTION SHOULD BE CALLED NO LESS THAN
ONCE PER YEAR, in current set-up, so that births per year ('bpy') and immigrants
per year ('ipy') can be updated at least each year. For generating mid-year
population estimates, this function should be called at least TWICE per year,
once at the beginning of the year and once mid-way through the year.

ENTRY: 'prog' contains the name of the program presently running.
       't' contains the current time.
       't0' contains the starting time.
       'N' contains the value of each dynamical variable, index by the state
         number ('qU', 'qV', 'qI1', ...).
       'startsec' contains the wall-clock time of the start of the run.
       'kernel' defines the contagion kernel.
       'deaths' and 'events' contain the number of deaths and events since
         the counters were cleared.
       //-'relativetime' is set if the simulated times are to be reported relative
       //- to the simulated starting time.
       'rand0' defines the random number sequence used.

EXIT:  The next line of the report has been printed.
       'deaths' and 'events' are cleared.
*/

static int ReportFirst;
ReportInit() { ReportFirst = 0; }

Report(char *prog)
{
  int i, ac,r,y,yr; dec z,age;

//-printf("Starting Report() function \n");
  if(ReportFirst==0)
  { ReportFirst = 1;
    printf("Dataset:     Simulation output of program '%s'\n", prog);
    printf("Kernel:      %s\n",
      kernel==0? "Mean field":
      kernel==1? "Cauchy":
      kernel==2? "Gaussian":
                 "Unspecified");
    printf("Sequence:    %lu\n\n", rand0);

    EventProfile("Initial");

    printf("Label t:       Time, in years and fractions thereof.\n");
    printf("Label N:       Total population size.\n");
    printf("Label Up:      Prevalence of susceptible individuals.\n");
    printf("Label Vp:      Prevalence of immune individuals.\n");
    printf("Label I1p:     Prevalence of new infections.\n");
    printf("Label I2p:     Prevalence of latent infections.\n");
    printf("Label I3p:     Prevalence of reinfections.\n");
    printf("Label D1p:     Prevalence of new/primary disease.\n");
    printf("Label D2p:     Prevalence of reactivation disease.\n");
    printf("Label D3p:     Prevalence of reinfection disease.\n");
    printf("Label D4p:     Prevalence of primary non-pulmonary disease.\n");
    printf("Label D5p:     Prevalence of react non-pulmonary disease.\n");
    printf("Label D6p:     Prevalence of reinf non-pulmonary disease.\n");
    printf("Label U:       Number of susceptible individuals.\n");
    printf("Label V:       Number of immune individuals.\n");
    printf("Label I1:      Number of new infections.\n");
    printf("Label I2:      Number of latent infections.\n");
    printf("Label I3:      Number of reinfections.\n");
    printf("Label D1:      Number of new/primary disease.\n");
    printf("Label D2:      Number of reactivation disease.\n");
    printf("Label D3:      Number of reinfection disease.\n");
    printf("Label D4:      Number of primary non-pulmonary disease.\n");
    printf("Label D5:      Number of react non-pulmonary disease.\n");
    printf("Label D6:      Number of reinf non-pulmonary disease.\n");

    printf("Label Deaths:  Number of deaths since last report.\n");
    printf("Label Events:  Number of events dispatched since last report.\n");
    printf("Label Elapsed: Seconds of elapsed wall-clock time to this point.\n");

    printf("\n|t     |N       |Up      |Vp      "
                    "|I1      |I2      |I3      "
                    "|D1      |D2      |D3      "
                    "|D4      |D5      |D6      "
                    "|U       |V       "
                    "|I1      |I2      |I3      "
                    "|D1      |D2      |D3      "
                    "|D4      |D5      |D6      "
                    "|Deaths  |Events  |Elapsed \n");
  }

  for(z=0,i=q0; i<=q1; i++) z += N[i];       //Get population size.

  printf("|%6.1f|%8.0f|%f|%f|%f|%f|%f|%f|%f|%f|%f|%f|%f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8.0f|%8d|%8d|%5d\n",
    t, z,
    N[qU]/z, N[qV]/z, N[qI1]/z, N[qI2]/z, N[qI3]/z,
                      N[qD1]/z, N[qD2]/z, N[qD3]/z,
                      N[qD4]/z, N[qD5]/z, N[qD6]/z,
                      N[qU],    N[qV],
                      N[qI1],   N[qI2],   N[qI3],
                      N[qD1],   N[qD2],   N[qD3],
                      N[qD4],   N[qD5],   N[qD6],
                      deaths, events, (int)(time(NULL)-startsec));

  fprintf(stderr, "  %.1f\r", t);            //Update status indicator.
  fflush(stdout); fflush(stderr);            //Make sure everything shows.
  deaths = events = 0;                       //Clear time-step counters.

  y = (int)t;                                //Get calendar (integer) year.

  if(y>lup)                                  //Check to see if parameters
  { ypb = 1./bcy[y-(int)t0];                 //sensitive to calendar year
    ypi = 1./immig[y-(int)t0];               //need updating.
    lup = y; }

  if((t-y)>0.3 && (t-y)<0.7 && y>1998)       //Since computation is expensive
  {                                          //only get population sizes when
                                             //necessary (mid-year).

    yr = y-(int)t0;                          //Get year array index.
    for(i=1; i<immid; i++)                   //Loop through immigrants.
    {
      r=0;                                   //If running SSA version of
      if(SSAV && A[i].ssa) r=2;               //model, find out if SSA.

      age = t-A[i].tBirth;                   //Get age and find age class.
      ac = age<15?0: age<45?1: age<65?2: 3;
      N2[ac][A[i].sex][r][yr] += 1; }        //Increment correct compartment
                                             //for this individual.

    r=1;                                     //Loop through UK-born.
    for(i=maximm+1; i<ukbid; i++)
    { age = t-A[i].tBirth;                   //Get age and find age class.
      ac = age<15?0: age<45?1: age<65?2: 3;
      N2[ac][A[i].sex][r][yr] += 1; }
  }
}

/* NOTES:

Pop. sizes by age, sex, region of birth, year: Best way is probably to schedule
in advance an 'ageing event' which would depend on their age and the
relevant age class divisions so that only the next one is scheduled of course.
This may also have implications for STD modelling e.g. Alternatively, this could
be done by having an array of individuals stored in order of birth year, as
Clarence and I worked out with regard to age-dependent mixing. This array would
be accompanied by information on where the first and last member of each age
class (e.g. year of birth) is. One could quickly get the population sizes
anytime during model run that way.

However, here it will be done by simply looping through the entire array, A, of
individuals whenever population sizes are needed. This way things are up and
running in a short time -- this can be improved in the future.
*/



/*----------------------------------------------------------------------------*
CLOSURE.

When processing is finished, this routine writes final statistics, prints
time trajectory of disease cases and returns an array with case notification
rates (or numbers of notifications) observed over the simulation.

ENTRY: 'startsec' contains the starting time, in the format of function
        'time'.
       'trho' and 'nrho' contain the total dispersal distance and the number
        of dispersal events.
       'tinfections' and 'linfections' contain the total number of infections
        targetted and the number that fell within the geographic area.
       'tstep' has accumulated statistics throughout the run.

EXIT:   Final statistics have been displayed.
        'out' contains the notification rates observed over the simulation.
*/

static dec nstep  =  0;                      //Number of steps
static dec tstep1 =  0;                      //Mean time step
static dec tstep2 =  0;                      //Standard deviation in time step
static dec tsmin  =  1E10;                   //Smallest time step
static dec tsmax  = -1E10;                   //Largest time step

static dec trho, nrho;                       //Statistics for local dispersal.
static dec tinfections, linfections;


FinalInit()
{
  nstep = tstep1 = tstep2 = 0;
  tsmin  =  1E10; tsmax  = -1E10;

  trho = nrho = 0;
  tinfections = linfections = 0;
}

Final()
{ int a,s,r,y,d; dec size,w;
  FILE *cases, *pop;

  printf("\n");
  size  = (indiv+3) * sizeof(struct Indiv);
  size += EventProfile("Final");

  tstepfin();
  { printf("Time steps:      Mean %s, Min %s, Max %s, SD %s, N %.0f\n",
            Tval(tstep1), Tval(tsmin), Tval(tsmax), Tval(tstep2), nstep); }

  if(nrho)
    printf("Dispersal:       Mean distance %.1f grid units.\n", trho/nrho);

  { printf("Infections:      Targeted %.0f, out of area %.0f, ratio %.2f%%\n",
      tinfections, tinfections-linfections,
      100.*(tinfections-linfections)/tinfections); }

  if(agec[0])
  { age1[0] /= agec[0]; age2[0] = sqrt(age2[0]/agec[0] - age1[0]*age1[0]);
    printf("All individuals: Mean age %.1f, SD %.1f, N %.0f\n",
      age1[0], age2[0], agec[0]); }

  if(agec[1])
  { age1[1] /= agec[1]; age2[1] = sqrt(age2[1]/agec[1] - age1[1]*age1[1]);
    printf("Disease-free:    Mean age %.1f, SD %.1f, N %.0f\n",
      age1[1], age2[1], agec[1]); }

  printf("\n");
  printf("Memory usage:    %.2f gigabytes\n", size/(1024*1024*1024));

  printf("Elapsed time:    %s\n",
    Tval((dec)(time(NULL)-startsec)/60/60/24/365.25));

  fprintf(stderr, "          \n");
  fflush(stdout); fflush(stderr);


  dec tot,tot2;                               //Population and case totals
                                              //for printing aggregated rates.

  dec ratio;                                 //Ratio of population sizes
                                             //observed in England and Wales to
                                             //population sizes produced by the
                                             //simulation.

/*
* tot=0; tot2=0
* for(y=(1999-(int)t0); y<RT; y++)           //Print notification rates
* { for(r=(1+SSAV); r>=0; r--)               //by rob and year, no age classes,
*   { tot=0; tot2=0;                         //to standard out.
*     for(s=0; s<2; s++)
*     for(a=0; a<4; a++)
*     { tot+=(repc[a][s][r][0][y]+repc[a][s][r][1][y]);
*       tot2+=N2[a][s][r][y]; }
*     printf("%f\t",100000*tot/tot2);
*   }
*   printf("\n");
* }
*/

 outi=0;
 printf("Printing all notification rates by age, sex, and rob\n");
 printf("M,0-14\tM,15-44\tM,45-64\tM,65+\tF,0-14\tF,15-44\tF,45-64\tF,65+\n");
 printf("\n");                               //Print all notification rates
 for(r=0; r<=(1+SSAV); r++)                  //by region of birth,year, sex &
 { for(y=(1999-(int)t0); y<RT; y++)          //age to 'out' & stdout, which can
    { for(s=0; s<2; s++)                     //be captured with grep and '|'.
      for(a=0; a<4; a++)
      { w=100000*(repc[a][s][r][0][y]+repc[a][s][r][1][y])/N2[a][s][r][y];
        printf("|%f ",w);
        out[outi++]=w; }
      printf("\n"); }
    printf("\n"); }

/*  UNADJUSTED
* printf("Printing POPULATION NUMBERS \n");
* for(r=(1+SSAV); r>=0; r--)                 //Print model population sizes.
* { printf("For rob=%d: M,0-14\tM,15-44\tM,45-64\tM,65+\tF,0-14\tF,15-44\tF,45-64\tF,65+\n",r);
*   for(y=(1999-(int)t0); y<RT; y++)
*   { printf("-");
*     for(a=0; a<4; a++)
*     for(s=0; s<2; s++)
*       printf("\t%f",N2[a][s][r][y]);
*     printf("\n"); }
*   printf("\n");
* }
* printf("\n");
*/



/* ADJUSTMENT OF 'repc'(REPORTED CASES IN MODEL) BASED ON POPULATION SIZES
   POPULATION SIZES OBSERVED IN ENGLAND AND WALES VERSUS POPULATION SIZES
   PRODUCED BY THE SIMULATION  */
/*
  printf("Now printing N3 for testing \n");
  FileIO("N3test.txt", fmt[21], "w|");
*/
// FileIO("repc1.txt", fmt[22], "w|");       //Print reported cases before
                                             //adjustment.
//printf("Now printing ratios for testing \n");
  for(y=1999-(int)t0; y<RT; y++)             //Loop over reported cases
  for(r=0; r<(2+SSAV); r++)                  //by year,region of birth, sex,
  for(s=0; s<2; s++)                         //and age class to correct for
  for(a=0; a<4; a++)                         //observed versus actual
  { if(N2[a][s][r][y]==0)                    //population sizes, taking care
      N2[a][s][r][y]=1;                      //to avoid dividing by zero.
    ratio = N3[a][s][r][y]/N2[a][s][r][y];
    //printf("%f\t",ratio);
    for(d=0; d<2; d++)
      repc[a][s][r][d][y]*=ratio; }

// FileIO("repc2.txt", fmt[22], "w|");       //Write out reported cases after
                                             //adjustment.

  outni=0;
  printf("Printing all case notifications by age, sex, and rob\n");
  printf("M,0-14\tM,15-44\tM,45-64\tM,65+\tF,0-14\tF,15-44\tF,45-64\tF,65+\n");
  printf("\n");                              //Print all notifications
  for(r=0; r<=(1+SSAV); r++)                 //by region of birth,year, sex &
  { for(y=(1999-(int)t0); y<RT; y++)         //age to 'outn' & stdout, which can
    { for(s=0; s<2; s++)                     //be captured with grep and '|'.
      for(a=0; a<4; a++)
      { w=repc[a][s][r][0][y]+repc[a][s][r][1][y];
        printf("|%f ",w);
        outn[outni++]=w; }
      printf("\n"); }
    printf("\n"); }


  #ifndef main
  #include "plotting.c"
  #endif


/*
* printf("Printing NUMBERS of notifictions\n");
* for(r=(1+SSAV); r>=0; r--)                 //Print -numbers- of notifications
* { printf("For rob=%d: M,0-14\tM,15-44\tM,45-64\tM,65+\tF,0-14\tF,15-44\tF,45-64\tF,65+\n",r);
*   for(y=(1999-(int)t0); y<RT; y++)         //by year, rob, and sex, to be
*   { printf("-");                           //captured with grep and '-'.
*     for(a=0; a<4; a++)
*     for(s=0; s<2; s++)
*       printf("\t%f",repc[a][s][r][0][y]+repc[a][s][r][1][y]);
*     printf("\n"); }
*   printf("\n");
* }
* printf("\n");
*/

}


/*----------------------------------------------------------------------------*
TIMING STATISTICS

The system can take very small time steps or very large, depending on the number
and frequency of events. This routine keeps track of individual time steps for
statistical tracking, to be reported at the end of the run.

ENTRY: 't' contains the current time.
       'tn' contains the next instant of time.

EXIT:  Statistics have been accumulated for the time step.
*/

tstep(dec t, dec tn)
{ dec dt;

  dt = tn-t;                                 //Compute the length of the step.

  tstep1 += dt;                              //Accumulate the sum and sum of
  tstep2 += dt*dt;                           //squares.

  if(tsmin>dt) tsmin = dt;                   //Accumulate the minimum and
  if(tsmax<dt) tsmax = dt;                   //maximum.

  nstep  += 1;                               //Accumulate the count.
}

/*
TIMING RESULTS

ENTRY:  'tstep' has accumulated statistics throughout the run.

EXIT:   'nstep'  contains the number of time steps.
        'tstep1' contains the length of the average time step, in years.
        'tstep2' contains the root-variance in length of the time steps.
        'tsmin'  contains the shortest time step.
        'tsmax'  contains the longest time step.

NOTE: Here it is called "root variance" rather than "standard deviation"
because the division is by 'n', not 'n-1'.
*/

tstepfin()
{
  if(nstep==0) return;                       //Avoid division by zero.
  tstep1 /= nstep;                           //Compute the mean.
  tstep2 = tstep2/nstep - tstep1*tstep1;     //Compute the variance.
  tstep2 = sqrt(tstep2);                     //Convert to root-variance.
}



/*----------------------------------------------------------------------------*
CHECK FOR MONOTONICITY

This routine checks whether a table of cumulative probabilities is monotonically
increasing and optionally whether it is bracketed by 0 and 1.

ENTRY: 'p' is the table of cumulative pobabilities.
       'n' is the number of entries in the table.
       'b' is set if the table should begin with 0 and end with 1.
       'r1' and 'r2' contain two numbers that may help to identify the location
         of the error. If such numbers will not help, then either or both
         contain zero.

EXIT:  The routine returns if the table appears to be correct. If not, an error
         message is issued and the routine never returns.
*/

int monotone(dec p[], int n, int b, int r1, int r2)
{ int i;

  for(i=1; i<n; i++)                         //Make sure the sequence never
    if(p[i-1]>p[i])                          //decreases.
      Error3(621., "`",r1, "` ",r2, "` ",i);

  if(b && (p[0]!=0||p[n-1]!=1))              //If requested, make sure it begins
    Error2(622., "`",r1, "` ",r2);           //with 0 and ends with 1.

  return 0;                                  //(Will never reach this).
}



/*----------------------------------------------------------------------------*
SERVICE ROUTINES
*/

char *pntab[] =                              //Table of parameter names.
{ "s2[0]","s2[1]","c[0][0]","c[0][1]","c[1][0]","c[1][1]",
  "v1[0]","v1[1]", "v2[0]", "v2[1]","v3[0]","v3[1]","ehiv",
  "r1[0]","r1[1]", "r2[0]", "r2[1]","r3[0]","r3[1]",
  "r4[0]","r4[1]", "r5[0]", "r5[1]","r6[0]","r6[1]",
  "r7[0]","r7[1]", "r8[0]", "r8[1]", "df",
  "d1uk20", "d2uk20", "d3uk20",
  "pmale[0]", "randseq", 0 };

dec *patab[] =                               //Table of parameter addresses.
{ &s2[0], &s2[1], &c[0][0],&c[0][1],&c[1][0],&c[1][1],
  &v1[0], &v1[1], &v2[0], &v2[1], &v3[0], &v3[1], &ehiv,
  &r1[0], &r1[1], &r2[0], &r2[1], &r3[0], &r3[1],
  &r4[0], &r4[1], &r5[0], &r5[1], &r6[0], &r6[1],
  &r7[0], &r7[1], &r8[0], &r8[1], &df,
  &d1uk20[0], &d2uk20[0], &d3uk20[0],
  &pmale[0], &randseq, 0 };

#include "service.c"



/*============================================================================*
NOTES 
_____________________________________________________
Before running on supercomputer:
1) Change 'SUPER' to '1'
2) In common.h change #define indiv to higher value
_____________________________________________________


*/
