/*----------------------------------------------------------------------------*
RANDOM NUMBERS FROM ANY DISTRIBUTION, WITH RESAMPLING

Given an arbitrary cumulative probability distribution and a given value within
that distribution, this routine generates random numbers that represent that
part of the distribution where numbers equal or exceed the given value.

That is confusing, but translated to an application---age specific
mortality---it is less so. For that application, given an arbitrary distribution
defining the probability, measured at birth, of a newborn dying on or before a
given age, and given an age already achieved by an individual in this
population, this routine returns a random number drawn from the distribution
telling probabilistically how much longer the individual will live.

Mathematically, 'P(x)' is the cumulative distribution, with 'P(0)=0' and
'P(infinity)=1', The given value is 'g'. Random numbers are selected from a
transformed distribution, 'F(x) = (P(x+g)-P(g)) / (1-P(g))'. Note that all
"memoryless" distributions are invariant under this transformation, including
the cumulative exponential distribution, 'P(x)=1-e^(-r x)' and the trivial
distribution, 'P(x)=0', where everything lives forever. The trivial distribution
is a special case of the exponential.

Operationally, the cumulative distribution is supplied as a table of values 'V'
matched one-to-one with a table of probabilities 'P'. Below is an example of a
matched pair of tables corresponding to the projected probability of mortality
by age 'V' among males born in the United Kingdom in years 2003-2005. The table
represents 121 years but needs only 63 entries because of near linearities over
certain ranges of the table. In practice, the tables might, for simplicity and
uniformity, be built with one entry per year of age regardless of near
linearities in the data.

   V  P         V  P         V  P         V  P         V  P         V  P
 --- -------  --- -------  --- -------  --- -------  --- -------  --- -------
   0 0.00000   48 0.04596   61 0.11653   73 0.30497   86 0.73269   99 0.99042
   1 0.00564   50 0.05246   62 0.12622   74 0.32979   88 0.80094  100 0.99373
  16 0.00810   51 0.05619   63 0.13716   75 0.35663   89 0.83261  121 1
  21 0.01086   52 0.06025   64 0.14872   76 0.38500   90 0.86191
  27 0.01530   53 0.06469   65 0.16142   77 0.41529   91 0.88735
  31 0.01857   54 0.06939   66 0.17486   78 0.44749   92 0.90970
  35 0.02262   55 0.07448   67 0.18937   79 0.48086   93 0.92918
  39 0.02764   56 0.08005   68 0.20518   80 0.51585   94 0.94588
  41 0.03066   57 0.08597   69 0.22214   81 0.55166   95 0.95917
  43 0.03413   58 0.09259   70 0.24052   83 0.62546   96 0.97040
  45 0.03829   59 0.09967   71 0.26015   84 0.66231   97 0.97908
  47 0.04314   60 0.10752   72 0.28180   85 0.69795   98 0.98562

For a different example, consider a hypothetical manufactured device exposed to
periodic dangers as it moves through a hazardous environment. Suppose the
devices have high probability of failing immediately upon being deployed,
precisely at the one-day boundary, and also precisely at the two-day boundary.
After two days the hypothetical devices have passed through the region of
environmental danger and do not fail. One-fourth of the devices fail at each
stage, leaving one-fourth intact at the end. That would be represented by the
following probability tables passed to this routine.

                           V        P
                 -----------        ----
                           0        0
                           0        0.25
                           1        0.25
                           1        0.50
                           2        0.50
                           2        0.75
                 10000000000        0.75
                 10000000000        1

The value "10000000000" above (tenth power of ten) represents infinity. If the
computer arithmetic is able to encode transfinite values, that encoding can be
used. Otherwise a very large finite value, as above, can be substituted.

The table does not really represent a function, since it does not assign a
single value 'P' to each value 'V'. It is what is called a "relation", which
assigns a set of values 'P' to each value 'V"' In this particular case, the
way to think of it is as the path a pen traces on graph paper, from (0,0) ->
(0,1/4) -> (1,1/4) -> (1,1/2) -> (2,1/2) -> (2,3/4) -> (infinity,3/4) ->
(infinity,1). Of course, the pen wouldn't make it all the way to infinity.

For a final example, although there are faster ways to generate such numbers,
here is a two-entry table to generate random numbers uniformly distributed
between -1 and 1:
                           V    P
                          ---  ---
                          -1    0
                           1    1

[The faster way is '2*Rand()-1'.]

ENTRY: 'V' is a table of strictly increasing values in the set of random
         numbers to be generated.
       'P' is a table of probabilities, each being the probability that a random
         value will be less than or equal to the corresponding value in 'V'.
       'n' is the number of entries in tables 'V' and 'P'.
       'g' is given value in the range 'V[0]' to 'V[n-1]', inclusive.

EXIT:  'RandF' contains a random value drawn from the given distribution,
         starting at value 'g'.
*/

typedef double dec;
typedef double decs;
//typedef float  decs;

dec Val(), Rand();

dec RandF(decs V[], decs P[], int n, dec g)
{ int i; dec r, p, w;

  if(V[0]>g  || V[n-1]<g)  Error(753.1);     //Check the bounds of both tables.
  if(P[0]!=0 || P[n-1]!=1) Error(753.2);

  r = Rand();                                //Generate a uniform random value.
  //printf("%f\t",r);   //FOR TESTING

  if(g!=V[0])                                //Rescale the random value if only
  { p = Val(1, g, V,P, 0,n-1);               //part of the distribution is to be
    r = p + r*(1-p); }                       //sampled.

  i = Loc(P, 0, n, r);                       //Pick a value from the portion of
  w = P[i+1]-P[i];                           //the inverse cumulative distribution
  if(w) w = (r-P[i])/w; else w = 1;          //that applies to the value of 'g'.
  return V[i]-g + w*(V[i+1]-V[i]);
}

/* Note: This routine and the subroutines it calls are simple in their structure
but involve some subtleties in their design. Modest improvements in speed could
be obtained at the cost of complexity in several ways: (1) Trade space for time
by storing the cumulative distribution function as a direct-access table, with
one entry per lattice point in the space of 'V'. That would eliminate the binary
search involved in the function 'Val' and retain it in the function 'Loc'. It
would work for the human-life-table example but not for the manufactured-device
example. (2) Trade much more space for time by storing the inverse cumulative
function as a direct-access table, with one entry per lattice point in the
probability space 'P'. This would eliminate the binary search entirely. It would
work if the slope of the cumulative function never gets very close to zero. (3)
Rescale the random number differently so that the entire table need not be
searched, but only the part covering value 'g' and above. That could eliminate a
few calls of binary-search recursion if 'g' was large.
*/

/*----------------------------------------------------------------------------*
FUNCTION EVALUATION

This routine determines the value 'y' of a function at a specific point 'x',
where the function is defined by a table of values for corresponding points
'(x,y)'. Table 'X' is strictly increasing.

For example, suppose the following tables were passed with parameters 'k=1',
'x=0.5', 'i0=0', and 'i1=4'.

       X       Y
      -1       3
       0       0
       2       2
      10       0

The routine would first locate the value 'x=0.5' as halfway along the
independent variable range -1 and 0. It would then compute and return 1.5, the
value halfway between the corresponding dependent variable range 3 and 0.

ENTRY: 'k' defines the type of interpolation.
        0 None      (not implemented yet)
        1 Linear
        2 Quadratic (not implemented yet)
        3 Cubic     (not implemented yet)
        5 Quintic   (not implemented yet)
       'x' contains the value of the independent variable.
       'X' and 'Y' define the independent and dependent variables, respectively.
       'i0' and 'i1' define the indexes of the first and last entries,
         respectively, in both tables 'X' and 'Y'.

EXIT:  'Val' contains the value of the function at point 'x'. If 'x' is
        below or above the range defined in table 'X', the minimum or
        maximum value, respectively, in table 'Y' is returned.
*/

dec Val(int k, dec x, decs X[], decs Y[], int i0, int i1)
{ int i; dec w;

  if(x<=X[i0]) return Y[i0];                 //Handle variables outside of the
  if(x>=X[i1]) return Y[i1];                 //normal range.

  i = Loc(X, i0, i1-i0+1, x);                //Bracket the independent variable.

  w = X[i+1]-X[i];                           //Interpolate linearly within the
  if(w) w = (x-X[i])/w; else w = 1;          //bracketed range, watching for.
  return Y[i] + w*(Y[i+1]-Y[i]);             //discontinuities.
}

/*----------------------------------------------------------------------------*
TABLE LOOK-UP

This is a recursive binary search to process an ordered table of 'n' entries in
time proportional to 'log2(n)'. It requires the entries to be strictly
increasing, which may require a small increment (e.g., '1E-10') added in the
case of equal entries.

ENTRY: 'T' addresses a strictly increasing table of two or more values.
       'b' contains the beginning entry to be examined in 'T'.
       'n' contains the number of entries to be examined in 'T', at least 2.
       'v' contains the value to be located in 'T', with 'T[b] <= v <= T[b+n-1]'.

EXIT:  'Loc' indexes the local pair of table entries containing 'v', such that
        'T[loc] <= v <= T[loc+1]'.
*/

int Loc(decs T[], int b, int n, dec v)
{
  int m = n/2 + n%2;
  return m<=1? b: v<T[b+m-1]? Loc(T,b,m,v): Loc(T,b+m-1,n-m+1,v);
}

// CLARENCE LEHMAN AND ADRIENNE KEEN, AUGUST 2010.

