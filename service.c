/* SERVICE ROUTINES

Collected here are general service routines common to different versions of
the main program. These routines are attached with a "#include" statement
rather than being compiled as a separate module, which allows some of them
access to data structures in the main program.

/*----------------------------------------------------------------------------*
CONTAGION KERNEL

Two contagion kernels are in this version. The method is to generate a
probability, then run it backwards, through the inverse cumulative contagion
kernel, to obtain a spatial displacement from the source of the infection. The
individual at that displacement is the target.

This program uses continuous time but not continuous space. It is not yet clear
how to make continuous space efficient, but a regular spatial lattice is
extremely efficient. This version makes space into a torus.

For actual applications, arbitrary two-dimensional kernels of any complexity
can be added without loss of speed, using the cumulative inverse method
demonstrated here. That requires a direct-access, high-speed table lookup, to
be added in actual versions.

ENTRY: 'n' indexes the individual about to infect another.
       'kernel' selects the contagion kernel.
        0 Mean field - all individuals in the population are equally accessible.
        1 Cauchy kernel - nearby individuals are more accessible.
        2 Gaussian kernel - similar to above but tighter (not implemented).
       'indiv' contains the number of individuals in the simulation.
       'trho' and 'nrho' contain the total dispersal distance and the number
        of dispersal events, for later calculation of the mean dispersal
        distance.
       'tinfections' and 'linfections' contain the total number of infections
        targetted and the number that fell within the geographic area.

EXIT:  'Kernel' contains the number of the individual targeted for infection.
       'trho', 'nrho', 'tinfections', and 'linfections' are updated as
        necessary.
*/

#define PI 3.141592653589793238462643

//dec trho, nrho, tinfections, linfections;    //Statistics for local dispersal.

/*

int Kernel(int n)
{ int i, j, k, di, dj;
  dec dx, dy, rho, theta;

  switch((int)kernel)
  { case 0:                                  //MEAN-FIELD KERNEL.
      tinfections += 1; linfections += 1;    //Advance the number of infections.
      do k = 1+Rand()*indiv; while(k==n);    //Select equally from all
      return k;                              //individuals in the system.
//-tinfections could be incremented here? Which would mean infections attempted
//-rather than accomplished.
    case 1:                                  //CAUCHY KERNEL.
      while(1)
      { tinfections += 1;                    //Advance the number of infections.

        rho = abs(Cauchy(0., sigma));        //Generate a Cauchy radius and
        trho += rho; nrho += 1;              //maintain dispersal statistics.

        theta = 2*PI*Rand();                 //Generate a random angle and
        dx = rho*cos(theta); di = round(dx); //convert angle and radius to a
        dy = rho*sin(theta); dj = round(dy); //lattice displacement.

        i = (di+A[n].i);                     //Add the displacement to the
        j = (dj+A[n].j);                     //location of this individual and
        if(i<0 || i>=gi) continue;           //make sure it is not out of the
        if(j<0 || j>=gj) continue;           //area.

        linfections += 1;                    //Advance the number of infections
        k = G[i][j]; if(k!=n) break;         //that fall within the geographic
      }                                      //area and return the individual
      return k;                              //occupying that place.

    default: Error1(714., "",kernel);        //(Invalid kernel number.)
  }
}

*/


/*----------------------------------------------------------------------------*
UNIFORM DISTRIBUTION

This routine generates random numbers within a specified interval from a
uniform distribution --- that is, were all values are equally likely.

ENTRY: 'a' defines the lower bound of a uniform distribution.
       'b' defines the upper bound.

EXIT:  'Uniform' contains a random number uniformly distributed in the
        range 'a' to 'b'.
*/

dec Uniform(dec a, dec b)
{
  return Rand()*(b-a)+a;
}

/*----------------------------------------------------------------------------*
EXPONENTIAL DISTRIBUTION

This function generates random time intervals between independent events such
that the events obey a Poisson distribution with mean and variance 'lambda*dt'
in any time unit 'dt'. The resulting intervals obey an exponential
distribution with mean '1/lambda' and variance '1/lambda**2'.

ENTRY: 'lambda' contains the average number of events per unit time, which
        must be greater than zero.
       'LIMITG' defines a multiple of the average time interval. No time step
        will be greater than 'LIMITG/lambda'.

EXIT:  'Expon' contains the next Poisson time interval.
*/

#define LIMITG 10

dec Expon(dec lambda)
{ dec r, expdt;

  while(1)                                   // Generate a uniformly distributed
  { r = Rand(); if(r==0) continue;           // random number.

    expdt = -log(r);                         // Convert it to an exponential
    if(expdt>LIMITG) continue;               // distribution and control its
    if(expdt==0)     continue;               // range.

    return(expdt/lambda);                    // Scale and return the next time
  }                                          // interval.
}

/* NOTE: This routine must make practical concessions to the reality of
fixed-point arithmetic. First, the random number generator can return a value
that is precisely equal to zero. In infinite-precision arithmetic, that would
occur with probability zero, but here it occurs with finite probability, so it
must be ignored to avoid singularities that would otherwise halt execution.
Second, and related, the finite granularity of the random number generator
near zero seems able to generate excessively large time intervals, so those
are ignored also. The entry parameter 'LIMITG' helps accomplish that. Finally,
intervals of zero are prevented to avoid violating the principle that only one
event can occur at a time. An interval of zero would only happen if 'Rand()'
returned a value of unity. It should not, but this routine provides an
additional safeguard against floating-point rounding anomalies.
*/

/*----------------------------------------------------------------------------*
NORMAL DISTRIBUTION

This routine generates random numbers from the Gaussian distribution,
f(x) = (1/sqrt(2*PI))*exp(-x*x/2), then adjusts them to a specified mean and
standard deviation.

\ENTRY:  'mu' contains the mean of the distribution.
         'sigma' contains the standard deviation.
         The random sequence of 'Rand' is ready to use.

\EXIT:   'Gauss' contains a Gaussian-distributed random number of the
          specified mean and standard deviation.
*/

dec Gauss(dec mu, dec sigma)
{ dec v1, v2, w;

  do { v1 = 2.*Rand() - 1.;                  //Generate a point within the
       v2 = 2.*Rand() - 1.;                  //unit circle.
       w = v1*v1 + v2*v2; }
  while(w>1 || w==0);

  w = v2 * sqrt(-2.*log(w)/w);               //Return the normal deviate.
  return mu + w*sigma;
}

/* NOTE: The routine uses the polar method due to G.E.P.Box, M.E.Muller, and
G.Marsaglia, described in D.E.Knuth, The Art of Computer Programming, Volume
2 (Addison Wesley, 1969).  This method involves, roughly, picking a point
uniformly distributed within the unit circle and then projecting it.

NOTE: It would be possible to save random number 'v2' in a static variable and
use it as 'v1' on the next call. That would make the routine slightly faster at
the cost of a little complexity.
*/

/*----------------------------------------------------------------------------*
LOGNORMAL DISTRIBUTION

This routine generates random numbers from the lognormal distribution, given
the mean and standard deviation of the underlying normal distribution.

\ENTRY:  'mu' contains the mean of the underlying normal distribution.
         'sigma' contains the standard deviation.
         The random sequence of 'Rand' is ready to use.

\EXIT:   'LogNormal' contains a lognormally-distributed random number of the
          specified parameters.
*/

dec LogNormal(dec mu, dec sigma)
{
  return exp(mu + sigma*Gauss(0., 1.));
}

/*----------------------------------------------------------------------------*
CAUCHY DISTRIBUTION

This routine generates random numbers from the Cauchy distribution,
f(x) = 1/(1+x*x)/PI, adjusted to specified median and width parameters.

\ENTRY:  'mu' contains the median of the distribution. Note that this is not
          the mean, since the Cauchy distribution has no mean.
         'sigma' contains the width of the distribution, at which the height
          falls to half its maximum value. Note that this is not a standard
          deviation, since the Cauchy distribution has no standard deviation.
         The random sequence of 'Rand' is ready to use.

\EXIT:   'Cauchy' contains a Cauchy-distributed random number of the
          specified median and width.

NOTE: Values are generated from the inverse cumulative Cauchy distribution.
I am not sure I have the width parameter is precisely right. To be checked out
further.
*/

dec Cauchy(dec mu, dec sigma)
{
  return mu + sigma*tan(PI*(Rand()-0.5));
}

/*----------------------------------------------------------------------------*
PARAMETERS.

Parameters have default values that can be overridden by command-line arguments.
Each entry on the command line begins with the parameter name, followed by an
equal sign (=), followed by the parameter value. For example, to set a
parameter named 'ralpha0' to 0.21 and another named 'beta1' to 0.718, and also
to clear three parameters named 'mu0', 'mu1', and 'mu2', the following
command-line string could be used:

       'ralpha0=.21 beta1=.718 mu0=mu1=mu2=0'

ENTRY: 'argc' contains the number of arguments plus one.
       'argv' points to a list of arguments in its second and higher positions.
       'pntab' contains a list of valid parameter names.
       'pstab' contains pointers to the corresponding parameter values.

EXIT:  All parameters correctly specified by 'argv' on entry are set.
        Error messages are issued for any others.
*/

#define STR 100

char msg001[] = "Parameter:   %s=%s\n";
char msg101[] = "E101. Parameter %d (%s) does not have the correct format"
                " (name=value).\n";
char msg102[] = "E102. Parameter %d does not contain a simple numeric"
                " value (contains \"%s\").\n";
char msg103[] = "E103. Parameter %d (%s) is not a recognized name.\n";

gparam(int argc, char *argv[])
{ static int lines;
  int i, j, k, m, n;
  char str[STR+1], *cval;
  dec val, atof();

  for(i=1; i<argc; i++)                      //Advance to the next parameter.
  { n = strlen(argv[i]); if(n<=0) continue;

    for(k=n-1; k>=0; k--)                    //Scan for the beginning of the
      if(argv[i][k]=='=') break;             //parameter value and issue an
    if(k<0)                                  //error message if no equal signs
    { printf(msg101, i, argv[i]);            //are to be found.
      continue; }
    k += 1;

    for(m=n=0, j=k; argv[i][j]; j++)         //Make sure the parameter value is
    { if(argv[i][j]=='-' && j==k) continue;  //a simple decimal number.
      if(argv[i][j]>='0' && argv[i][j]<='9') { m += 1; continue; }
      if(argv[i][j]=='.')                    { n += 1; continue; }
      break; }

    if(argv[i][j] || m<1 || n>1)             //If it is not, issue an error
    { printf(msg102, i, &argv[i][k]);        //message and skip processing it.
      continue; }

    cval = &argv[i][k];                      //Otherwise convert it to internal
    val  = atof(cval);                       //form.

    for(k=0; ; k+=j+1)                       //Begin looping through all names.
    {
      for(j=0; j<STR && argv[i][k+j]; j++)   //Copy the next name until an
      { if(argv[i][k+j]=='=') break;         //equal sign is encountered.
        str[j] = argv[i][k+j]; }
      if(argv[i][k+j]==0) break;
      str[j] = 0;

      for(n=0; pntab[n]; n++)                //Search the table of parameter
        if(strcmp(pntab[n], str)==0) break;  //names for a matching entry.

      if(pntab[n]==0)                        //Issue an error message if the
      { printf(msg103, i, str);              //name is not in the table.
        continue; }

      *patab[n] = val;                       //Set the parameter, display its
      printf(msg001, str, cval);             //name and value, and repeat.
      lines += 1;
    }
  }
  if(lines) printf("\n");                    //Leave a blank line at at the
}                                            //end if anything was printed.



/*----------------------------------------------------------------------------*
TIME FORMATTING

This routine formats a time value, presented as years and fractions thereof,
so that it is labeled as years, months, days, etc.

ENTRY: 'tval' contains the time value, in years.

EXIT:  'Tval' points to a formatted version of 'tval'. This is valid for the
        the next ten calls to this routine.
*/

#define SNUM 10

static int  snum;
static char sval[SNUM][100];
static dec  tunit[] = { 365.25, 24, 60, 60, 1000, 1000, 1000, 1000, 0 };
static char *tstr[] = { "year", "day", "hour", "minute", "second",
  "millisecond", "microsecond", "nanosecond", "femptosecond" };

char *Tval(dec tval)
{ int i;

  if(tval==0)
    i = 4;
  else for(i=0; tval<1.0 && tunit[i]; i++)
    tval *= tunit[i];

  sprintf(sval[snum], "%.2g %s%s", tval, tstr[i], &"s"[tval==1]);
  i = snum; snum += 1; if(snum==SNUM) snum = 0;

  return sval[i];
}

/*----------------------------------------------------------------------------*
SELECT EARLIEST EVENT

This routine searches for the earliest time in a subset of a table of times.

ENTRY: 'tab' contains a table of times of future events.
       'subset' contains a vector of indexes into 'tab' that are to be included
         in the search. The vector ends with a negative value.

EXIT:  'Earliest' contains the index of the earliest time, selected from
         vector 'subset'.
*/

int Earliest(dec tab[], int subset[])
{ int i, m; dec x, w;

  for(x=100000000,m=i=0; subset[i]>=0; i++)
  { w = tab[subset[i]]; if(x>w) { x = w; m = i; } }

  return subset[m];
}

/*-----------------------------------------------------------------------------
DISPLAY PARAMETERS

To provide as complete a record of the run as possible, this routine displays
all parameters that can be entered on the command line. Together with the
program source and the data files, runs should be able to be replicated with
this information.

ENTRY: 'pf' points to an output stream ready to receive text. This is either
         'stdout', 'stderr', or a stream opened by 'fopen' or 'popen'.
       'pntab' points to a list of parameter names.
       'patab' points to a corresponding list of parameter addresses.

EXIT:  A single line has been written to the output stream.
*/

int DisplayParam(FILE *pf)
{
  int i;

  fprintf(pf, "Parameters:");

  for(i=0; patab[i]; i++)
    fprintf(pf, " %s=%g", pntab[i], *patab[i]);

  fprintf(pf, "\n");
  return 0;
}

