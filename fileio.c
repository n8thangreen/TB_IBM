/*----------------------------------------------------------------------------*
FILE INPUT/OUTPUT

Scientific programs may be complex but they often have simple input/output
requirements---merely reading and writing multi-dimensional files of numeric
data. This module is to make that easy. It is inspired by the simple read/
write/ format statements of old-time Fortran.

The structure of the file in main memory (eg ram) and the file in secondary
storage (eg disk) are described in two strings of integers. For example, a
three-dimensional 2000 x 121 x 2 array defined as

    dec D[2000][121][2];

could be described by the integer string

    {-'i',2000,-'a',121,-'s',2} .

Index 'i' could represent a series of conditions, 'a' a sequence of ages, and
's' two sexes. The string defining the format in secondary storage could be

    {-'s',0,1,1,-'i',0,1999,1,-'a',0,99,1} .

This means that in reading or writing the file, index 's' would be run from 0 to
1 in steps of 1, within that 'i' would be run from 0 to 1999 in steps of 1, and
within that 'a' would be run from 0 to 99 in steps of 1. Therefore, of the
484,000 integers in the three-dimensional array, 200,000 representing sex 0
would appear first, followed by 200,000 representing sex 1. Each half would be
consist of 2000 "lines" of 100 age classes each. The remaining 21 ages classes
would remain unread or unwritten.


DEFAULT VALUES

The definition of the file in secondary storage can be simplified by default
values. The default beginning index is 0, the default ending index is the
maximum index for that dimension, and the default step size is 1. Thus the above
definition is equivalent to

    {-'s',-'i',-'a',0,99} .

For simplicity, if the array is to be processed in its natural order, as defined
in 'mm', the string 'sm' defining the format in secondary storage can be left
null.

When writing a file as defined above, it will appear as a single line with a
new-line code at the very end. Any index that is capitalized will cause a
new-line code at the end of each set of values for that index. For example, the
definition

    {-'s',-'i',-'A',0,99}

will write lines of 100 values each. The definition

    {-'s',-'I',-'a',0,99}

will write lines of 2000 x 200 = 400,000 values each, and

    {-'s',-'I',-'A',0,99}

will write lines of 99 values each, with a blank line every 2000 lines.


OUTPUT FORMATTING

Data written to secondary storage can be formatted in several simple ways. By
default, an entire multi-dimensional array is written as a single line, with
values separated by spaces and terminated with a new-line code. However, any
capitals in the formatting string force new-line codes at the end of the
corresponding loop, as described above. In addition, the first character after
the 'w' in the parameter specifies the field separator. For example, suppose an
array is defined as

    dec Z[5][3] =
    { {1,   2,   3,   4,   5   },
      { .1,  .2,  .3,  .4,  .5 },
      {1.1, 2.2, 3.3, 4.4, 5.5 } };

and the indexing specification is

    { {-'i',5,-'j',3}, {-'i',-'J'} }

Then this will be the output for various specifications:

1. "w" or "w ". Write data with spaces and newlines.

    1 2 3 4 5
    0.1 0.2 0.3 0.4 0.5
    1.1 2.2 3.3 4.4 5.5


2. "w\t" Write data with tabs and newlines.

    1	2	3	4	5
    0.1	0.2	0.3	0.4	0.5
    1.1	2.2	3.3	4.4	5.5


3. "w," Write data with commas and newlines.

    1,2,3,4,5
    0.1,0.2,0.3,0.4,0.5
    1.1,2.2,3.3,4.4,5.5


4. "w\n" Write every data element on a new line.

    1
    2
    3
    :  <- (10 lines omitted)
    4.4
    5.5


5. "w|" Write data in self-defining vertical-bar format (called Centinel
    Format). This includes column headings and a basic database title.

    Written by 'FileIO' as 'abc.txt', 2011/03/20 12:54 (Sunday)
    |i |j0 |j1 |j2 |j3 |j4
    |0 |1  |2  |3  |4  |5
    |1 |0.1|0.2|0.3|0.4|0.5
    |2 |1.1|2.2|3.3|4.4|5.5


6. Finally, following the separator character is an optional formatting code, as
   used in 'printf', preceded by an equal sign. For example, the parameter
   "w|=%4.2f" will produce

    Written by 'FileIO' as 'abc.txt', 2011/03/20 12:54 (Sunday)
    |i |j0  |j1  |j2  |j3  |j4
    |0 |1.00|2.00|3.00|4.00|5.00
    |1 |0.10|0.20|0.30|0.40|0.50
    |2 |1.10|2.20|3.30|4.40|5.50


INPUT SCALING

Data can be scaled by a linear transformation as it is placed in main memory.
For example, if each data element 'x' is to be reduced to one-fifth its value
and then increased by one as it is read, the string "r=x/5+1" can be specified
on call to 'FileIO'. If it is to be first reduced by 1 and then doubled, that
would be specified with the string "r=x*2-2". The format is always one of four
forms, "r=x*m+b", "r=x*m-b", "r=x/m+b", or "r=x/m-b", where 'm' and 'b' are
numbers with optional decimal precision.

The results can be truncated to integers. To do that, 'n' is specified instead
of 'x', as in "r=n/5+1". That will perform the linear transformation as before,
then perform a nonlinear truncation.

If the parameter "r" is specified alone, that represents the identity
transformation, "r=x*1+0".


INTERNAL DATA STRUCTURES

Two related input/output packages are combined in this module, sharing many of
the same data structures and much of the same source code. The first is a simple
"stream of numbers" interface, where values appear one after the other separated
by non-numeric codes, with spaces, tabs, commas, line-feeds, or other
non-numeric codes delimiting the individual values.

The second package has a self-defining format called "Centinel". It was part of
a system developed by Adrienne Keen and Clarence Lehman in 2005, for preserving
small and moderate data files over a century-long time scale by printing them
with the same kinds of redundancy codes as are used on magnetic media, so the
data could be scanned and verified for accuracy in the distant future. Here it
is used only for the simplicity and flexibility of its data format.

Index labels in the arrays are restricted to 'a' through 'z', with 'z' reserved
for special use. Capital 'A' through 'Z' represent the same labels but carry
additional meaning when files are written (where they indicate where lines are
to end). Three vectors of 26 elements each are addressed by such labels. The
maximum index value for each label being processed is stored in 'clim', the
cumulative width of that dimension is stored in 'cwth', and the present value of
the index being processed is temporarily stored in 'cval'.

For example, consider an array 'dec W[5][10]', indexed as 'W[c][d]'. The index
'c' is the third letter in the alphabet and the third entry in the vectors, and
'd' is the fourth. Also, in the array index 'c' must be less that 5 and 'd' less
nthan 10. Therefore, 'clim[2]' will contain 5 and 'clim[3]' will contain 10.
Each row of the array contains ten integers, so 'cwth[2]' will contain 10 and
'cwth[3]' will contain 1. The linear index into the array is computed as
'cwth[2]*c+cwth[3]*d'.

Other important data structures include 'clabel[k]', where 'k' is a column
number. This contains the label corresponding to each column in the data file,
as a character 'a' to 'z'. Any index value associated with that column is in
'cindex[n]'. If no index is associated, '-1' is stored. The operation of
'cindex' is explained further below, and in the source code.

Index ranges are kept in 'drange[u][v]' as adjacent pairs of integers, to be
processed as the actual data are stored in the array. For example, if the value
3.14 is to be stored in rows 0 and 4, columns 0, 2, 3, 4, 7, 8, and 9, the range
specification of the file for 'c' could be '0,4' and for 'd' could be
'0,2~4,7~9'. Then 'drange' would contain two rows of data terminated with a
negative integer, as follows.

    drange[0]:  0 0   4 4  -1 0   0 0   0 0 ...
    drange[1]:  0 0   2 4   7 9  -1 0   0 0 ...

Finally, end-of-line indicators are stored in eol[u]. If set, the corresponding
dimension is to begin on a new line in the output file, or to be introduced by a
blank line for more-than-two-dimensional arrays.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fileio.h"

static struct IO io;                         //Input/output control structure.
static FILE *fp;                             //File control structure.
static char *filename;                       //Name of input/output file.
static char *rw;                             //Full read/write entry parameter.
static char *fstr;                           //Entry parameter beyond equal sign.
static int  imax;                            //Number of labels to be processed.
static dec  line;                            //Line number in file.
static dec  xm, xb;                          //Scaling slope and intercept.

static int  init, trunc, count, z;           //Initialization and work areas.
static char eol[MDIM+1];                     //Dimensions with line separators.
static int  clim[MLAB];                      //Number of entries for each label.
static int  cwth[MLAB];                      //Cumulative width for each label.
static int  cval[MLAB];                      //Current value of each label.

static char *extfstr(char*,char*);           //Function prototypes.
static int CentinelRead();
static int CentinelWrite(int,int,int,int,int,int);
static int StoreElement(dec,int);
static int BypLine(), StoreColumns(), StoreData();
static int advance(int), retract(int), preview(int), byp(int), posint(int);
static nest(int,int,int), mbparse(char*);
static field(char*,int), ranges(int,int*,int);


/*----------------------------------------------------------------------------*
MAIN ROUTINE

ENTRY: 'rw1' defines the kind of operation.
         "r"*  Read  (see above for the meaning of the asterisk).
         "w"*  Write (see above for the meaning of the asterisk).
       'fn1' is the name of the file. If null, 'stdin' or 'stdout' applies.
       'io1' defines the input/output structure.
       'io1.data' points to the data structure in memory.
       'io1.mm' points to the definition of the main memory structure, as
         described above.
       'io1.sm' points to the definition of the secondary memory structure, as
         described above.

EXIT:  Data from the file 'filename' has been read into the array specified
        in 'io1'.
       'FileIO' contains the number of array elements populated.
*/

int FileIO(char *fn1, struct IO io1, char *rw1)
{ int c, d, i, n, w; char op[2];

  filename = fn1; io = io1; rw = rw1;        //Save and check the calling
  if(rw[0]!='r' && rw[0]!='w') Error(525.);  //parameters.

  init = trunc = count = 0;                  //Initialize work variables.

  fstr = extfstr(rw, rw[0]=='w'? "%g":"");   //Extract any specification in the
  if(rw[0]=='r')                             //parameter and parse it if it
    if(fstr[0]) mbparse(fstr);               //represents a transformation of
    else        fstr = 0;                    //input data.

  for(i=0; i<MLAB; i++)                      //Clear the arrays addressed by the
    clim[i] = cwth[i] = cval[i] = 0;         //index label.

  for(i=0; i<MDIM; i++) eol[i] = 0;          //Clear the end-of-line array.

  for(i=0; i<MDIM*2 && io.mm[i]; i+=2)       //Count the number of dimensions
  { if(-io.mm[i]-'a'>=MLAB)                  //and audit the index labels and
      Error1(515.1, filename,0);             //sizes.
    if(io.mm[i+1]<1)
      Error1(516., filename,0); }
  n = imax = i/2;

  for(w=1; i>=2; i-=2,w*=d)                  //Compute and store the increasing
  { c = abs(io.mm[i-2]); d = io.mm[i-1];     //width of each dimension for use
    if(c>='A' && c<='Z') c += 'a'-'A';       //in the nested loops.
    clim[c-'a'] = d;
    cwth[c-'a'] = w; }

  for(z=MDIM*4; z>0; z--)                    //Determine the length of the
    if(io.sm[z-1]) break;                    //secondary memory string.

  if(filename && filename[0])                //Open the file for reading or
  { op[0] = rw[0]; op[1] = 0;                //writing, using the standard input
    fp = fopen(filename, op);                //or output stream if no file is
    if(fp==0) Error1(510., filename,0); }    //specified.
  else fp = rw[0]=='r'? stdin: stdout;

  if(rw[0]=='r' && rw[1]=='|')               //Read all input, either processing
    CentinelRead();                          //a self-defining file or looping
  else nest(0, 0, n-1);                      //through all dimensions.

  if(fp==stdin||fp==stdout) return 0;        //Close the file and return.
  fclose(fp); return count;
}

/*----------------------------------------------------------------------------*
NESTED LOOPS

Recursion is used here to achieve an arbitrary degree of nesting. This is
suitable for four reasons: (1) Recursion is almost as efficient as looping in
properly written modern compilers, (2) the fastest innermost loop is run
directly without recursion, (3) variables placed on the stack are kept to a
minimum by passing some as global static variables, and (4) this is an
infrequent operation connected with input/output, so the achieving the highest
possible speed is not of consequence.

ENTRY: 'io' contains the input/output file structure controlling the input/
         output operations.
       'fp' points to the file structure of the operating system.
       'cwth' is the array of widths of each dimension.
       'z' contains the number of elements in 'sm'.
       'o' indexes the current location in 'mm' or 'sm'.
       'base' is the base displacement into the array.
       'level' contains the nesting level. The innermost loop is level 0.
       'trunc' is set if data are to be truncated to integers on input.
       'fstr' is set if data are to be transformed on input.
       'xm' and 'xb' contain the slope and intercept, respectively, for any
         transformation.

EXIT:  The present loop and all nested loops at lower levels are completed.
*/

static nest(int o, int base, int level)
{
  dec w; int c, i, i0, i1, di, iw, sep[2];

  if(o>=z && io.sm[0]) Error(920.);          //Check the range of the offset.

  if(io.sm[0]) c = -io.sm[o++];              //Retrieve the present index label
  else       { c = -io.mm[o++]; o++; }       //and make sure it is within range.
  if(c<0||c-'a'>=MLAB) Error(515.2);

  if(c>='A' && c<='Z')                       //If the label is capitalized, mark
  { c += 'a'-'A'; eol[level] = 1; }          //a new-line and make sure the
  if(clim[c-'a']==0) Error(515.3);           //the label was defined.

  i0 = 0; i1 = clim[c-'a']-1; di = 1;        //Develop the indexing range and
  if(io.sm[0])                               //increments, using default values
  { if(o<z && io.sm[o]>=0) i0 = io.sm[o++];  //if 'mm' is being used or if they
    if(o<z && io.sm[o]>=0) i1 = io.sm[o++];  //are not specified in 'sm'.
    if(o<z && io.sm[o]> 0) di = io.sm[o++]; }
  iw = di*cwth[c-'a']; if(i0>i1) di = -di;

  if(i0>=clim[c-'a'])        Error(517.1);   //Make sure the index values are
  if(i1>=clim[c-'a'])        Error(517.2);   //within range and the range is
  if((abs(i1-i0)+1)%abs(di)) Error(518.0);   //divisible by the increment.

  if(level>0)                                //If this is not the lowest level,
  { for(i=i0; di>0?i<=i1:i>=i1; i+=di)       //loop through this level and all
    { cval[c-'a'] = i;                       //lower levels, adding new-line
      nest(o, base+iw*i, level-1); }         //codes at the end of specified
    if(rw[0]=='w' && eol[level])             //levels.
      fprintf(fp, "\n");
    return; }

  if(rw[0]=='r')                             //At the lowest level for reading,
  { for(i=i0; di>0?i<=i1:i>=i1; i+=di)       //read all data at this level into
    { if(fscanf(fp, "%lf", &w)<1)            //the array, optionally scaling it
        Error1(511., filename,0);            //by a linear transformation and
      if(fstr)  w = w*xm+xb;                 //optionally truncating to an
      if(trunc) w = (int)w;                  //integer.
      io.data[base+i*iw] = w; }
    return; }

  sep[0] = ' '; sep[1] = '\n';               //At the lowest level for writing,
  if(rw[1]) sep[0] = rw[1];                  //retrieve the field separator.

  if(sep[0]!='|')                            //If this is an ordinary format,
  { for(i=i0; di>0?i<=i1:i>=i1; i+=di)       //write just data and separators.
    { fprintf(fp, fstr, io.data[base+i*iw]);
      fprintf(fp, "%c", sep[i==i1 && eol[level]]); }
    return; }

  CentinelWrite(c, i0, i1, di, base, iw);    //Otherwise write a self-defining
}                                            //vertical-bar line.

/*----------------------------------------------------------------------------*
EXTRACT PARAMETER STRING

This routine extracts any string following an equal sign ('=') from a parameter.
For example, if the parameter contains "w\t=%12.4f", then this routine will
return "%12.4f".

ENTRY: 'f' points to a parameter.
       'g' points to a default string, or is null.

EXIT:  'extfstr' points to any string following the left-most equal sign (=). If
         no such string exists, 'extfstr' points to 'g' as passed on entry.
*/

static char *extfstr(char *f, char *g)
{ int i;

  for(i=0; f && f[i] && f[i]!='='; i++);
  if(f==0||f[i]==0||f[i+1]==0) return g;
  return &f[i+1];
}


/*----------------------------------------------------------------------------*
PARSE TRANSFORMATION

This routine parses a specification of the form "x*m+b" to obtain 'm' and 'b'.
A division operator '/' may appear instead of '*' and '-' may appear instead of
'+'. Either the '*m' or the '+b' fields, or both, may be omitted. For integer
truncation, 'n' may be used in the place of 'x'.

ENTRY: 'p' points to a parameter.

EXIT:  'xm' contains the multiplier (the slope).
       'xb' contains the bias (the intercept).
       'trunc' is set if results are to be truncated to integers.
*/

static mbparse(char *p)
{ int i; char *np;

  i = *p++; trunc = i=='n';                  //Check for integer truncation.
  if(i!='n' && i!='x') Error(520.);          //Check the start of the string.
  xm = 1; xb = 0;                            //Assume an identity transformation.

  switch(*p)
  { case 0: return;                          //No transformation specified.

    case '*':                                //Starting with multiplication.
      xm = strtod(++p, &np); goto casemd;

    case '/':                                //Starting with division.
      xm = strtod(++p, &np);
      if(xm==0) Error(521.);
      xm = 1./xm;

    casemd:                                  //Continuing multiplication and
      p = np; switch(*p)                     //division.
      { default:  Error(522.1);
        case 0:   return;
        case '+': case '-': break; }

    case '+': case '-':                      //Addition or subtraction.
      xb = strtod(p, &np);
      if(*np==0) return;

    default: Error(522.2); }                 //Invalid transformation.
}


/*----------------------------------------------------------------------------*
WRITE CENTINEL FORMAT

The Centinel format is for multidimensional numeric arrays that are too
complicated to be stored easily as mere streams of numbers. It allows comments
in the files and has column headings so that each numeric value can be placed
correctly in the array regardless of the ordering of the file. The full format
has redundancy codes on each line and each file so it may be printed and
reliably scanned, but those are not part of this implementation.

For example, consider a four-dimensional file defined with the string

    {-'s',2, -'r',3, -'y',29, -'q',9}

There are 2 values for sex ('s'), 3 for region of birth ('r'), 29 for year
('y'), and 9 for state ('q'). An excerpt of this file might look like the
following. (The data are simply random, not part of anything real.)

    Sample Centinel file

    Label s: Sex, 0=female, 1=male.
    Label r: Region of birth, 0=UK, 1=SSA, 2=Other.
    Label y: Year, relative to 1981.
    Label q: State transition parameter.

    |s |r |y |q0  |q1  |q2  |q3  |q4  |q5  |q6  |q7  |q8
    |0 |0 |0 |4.48|1.12|9.21|1.75|3.36|9.38|9.73|1.22|6.12
    |1 |0 |0 |3.95|8.54|1.61|3.98|6.31|7.43|4.51|5.30|5.03
    |0 |0 |1 |7.11|9.43|6.72|2.33|6.65|3.74|2.97|8.49|5.90
    |1 |0 |1 |4.41|1.92|2.86|8.74|6.55|4.35|2.37|7.75|6.01
    |0 |2 |1 |9.90|1.25|3.29|5.85|3.63|8.05|2.24|9.15|9.74
    |1 |2 |1 |8.45|4.67|2.80|1.83|5.17|8.71|9.92|1.24|8.27
     :  :  :   :    :    :    :    :    :    :    :    :

For purposes of this module, any lines beginning with a vertical bar ('|') are
data lines and other lines are comments. Notice that the lines can appear in any
order, for each line defines all the indexes of the array. This "self defining"
format is to help undetected data errors from infiltrating the program.

ENTRY: 'p' points to a parameter.

EXIT:  'xm' contains the multiplier (the slope).
       'xb' contains the bias (the intercept).
       'trunc' is set if results are to be truncated to integers.
*/

static char fm1[] = "Written by 'FileIO' as file '%s'\n";

static int CentinelWrite(int c, int i0, int i1, int di, int base, int iw)
{ int i, d;

  if(init==0)                                //If this is the first call, write
  { fprintf(fp, fm1, filename);              //a main heading line.

    if(c<'a') c += 'a'-'A';                  //Make the current label lower-case.

    for(i=0; i<MDIM*4; i++)                  //Write the higher-level index
    { if(io.sm[0]) d = -io.sm[i];            //names.
      else         d = -io.mm[i];
      if(d<=0) continue;
      if(d<'a') d += 'a'-'A';
      if(d<'a'||d>'z') Error(515.3);
      if(d==c) break;
      fprintf(fp, "|%c", d); }

    for(i=i0; di>0?i<=i1:i>=i1; i+=di)       //Write column labels for the
      fprintf(fp, "|%c%d",c,i);              //innermost index.
    fprintf(fp, "\n"); init = 1; }

  for(i=0; i<MDIM*4; i++)                    //On every call, write the current
  { if(io.sm[0]) d = -io.sm[i];              //index values of the higher
    else         d = -io.mm[i];              //dimensions.
    if(d<=0) continue;
    if(d<'a') d += 'a'-'A';
    if(d==c) break;
    fprintf(fp, "|%d", cval[d-'a']); }

  for(i=i0; di>0?i<=i1:i>=i1; i+=di)         //Write all data for the innermost
  { fprintf(fp, "|");                        //index.
    fprintf(fp, fstr, io.data[base+i*iw]); }
  fprintf(fp, "\n");

  return 0;                                  //Return to caller.
}

// NOTE: Must write one-dimensional arrays as two columns of numbers.


/*----------------------------------------------------------------------------*
READ CENTINEL FORMAT

This routine processes index fields that are a series of ranges separated by
commas. A range can be a single integer or a pair of integers separated by a
tilde ('~'). Here are some examples:

    Field            Indexes represented
      1                1
      0,1              0 1
      0~1              0 1
      0,3~9,40~38,2    0 3 4 5 6 7 8 9 40 39 38 2

ENTRY: 'fp' points to the file structure.

EXIT:  'count' is incremented by the number of array elements populated.
*/

#define EOC (c=='|' || c=='\n' || c=='\r')   //End of field.

static int  cmax;                            //Maximum column number.
static char clabel[MCOL+1];                  //Label for each column.
static int  cindex[MCOL+1];                  //Index for each column, or -1.
static int  drange[MDIM][CRANGE];            //Ranges for each dimension.
static int  dwth[MDIM];                      //Widths for each dimension.

cerr(dec v) { Error3(v, "`Line ",line,       //Write an error message with the
            " in file ",0., filename,0.); }  //location in the file appended.

static int CentinelRead()
{ int c;

  for(line=1;; line+=1)
  { c = preview(1); if(c==EOF) break;
    if(c!='|')  { BypLine();                continue; }
    if(init==0) { StoreColumns(); init = 1; continue; }
    count += StoreData(); }

  if(count==0) cerr(514.);
  return count;
}


/*----------------------------------------------------------------------------*
STORE COLUMN SPECIFICATIONS

All column names are either single-character lower-case index labels, or that
followed by a positive integer. This routine extracts those labels and integers.

ENTRY: The input file is positioned at the heading line of the data file.
       'imax' contains the number of dimensions.

EXIT:  'clabel' contains the index labels for all columns. followed by a
         negative terminator.
       'cindex' contains the numeric parts of the index values that are followed
         by integers, followed by a negative terminator.
*/

static int StoreColumns()
{ int c, n, f, i;

  for(cmax=0; cmax<MCOL; cmax++)             //Advance to the beginning of the
  { c = advance(0); if(c!='|') cerr(840.1);  //next column.

    byp(0); c = advance(0);                  //Extract the label and store it
    if(c<'a'||c>'z') cerr(524.1);            //for this column.
    clabel[cmax] = c;

    n = posint(0);                           //Extract the following integer,
    if(imax==1 && n<0 && c=='z')             //if any, and store it too, with
    { n = 0; clim[c-'a'] = cwth[c-'a'] =1; } //special handling for variable 'z'
    cindex[cmax] = n;                        //in one-dimensional arrays.

    byp(0); c = advance(0);                  //If another field follows
    if(c=='|') { retract(c); continue; }     //immediately, repeat.

    if(c!='\n' && c!='\r') cerr(523.);       //If the line has ended normally,
    for(f=i=0; i<=cmax; i++)                 //scan to make sure any index
    { if(cindex>=0) f += 1;                  //ranges precede all data values,
      else if(f) cerr(530.); }               //and that data values exist.
    if(f==0)     cerr(531.);

    clabel[cmax+1] = 0;                      //Terminate the string, advance to
    retract(c);                              //the next line, and return.
    BypLine(); return 0; }

  cerr(532.1); return 0;                     //(Too many columns.)
}


/*----------------------------------------------------------------------------*
STORE DATA

This routine processes each line of data, following the column specifications
recognized by 'StoreColumns'. Array indexes can be carried explicitly in
individual columns or implicitly in multiple columns. For example, consider the
following two-dimensional 3x3 array.

    0.11 0.16 0.23
    0.35 0.87 0.99
    0.39 0.26 0.99

With the rows labelled 'i' and the columns 'j', this could be represented
two-dimensionally as

    |i |j0   |j1   |j2
    |0 |0.11 |0.16 |0.23
    |1 |0.35 |0.87 |0.99
    |2 |0.39 |0.26 |0.99

or one-dimensionally as

    |i |j |z
    |0 |0 |0.11
    |0 |1 |0.16
    |0 |2 |0.23
    |1 |0 |0.35
    |1 |1 |0.87
    |1 |2 |0.99
    |2 |0 |0.39
    |2 |1 |0.26
    |2 |2 |0.99

When there is only one data value per line and all array indexes are in columns,
then the data column is labelled 'z', as above.  Ordering of the rows and
columns is not significant. For example, in this case the following random
ordering would be equivalent.

    |j |i   |z
    |1 |0   |0.16
    |0 |2   |0.39
    |0 |1   |0.35
    |2 |0   |0.23
    |0 |0   |0.11
    |1 |1   |0.87
    |1 |2   |0.26
    |2 |1,2 |0.99

The specification '1,2' for 'i' in the last line means that the value 'z' is to
be stored both at 'i=1,j=2' and 'i=2,j-2', as described earlier in this module.

ENTRY: The file is positioned at the beginning of a data line.

EXIT:  'StoreData' contains the number of array elements populated on this call.

WORK:  'drange' and 'dwth' contain the ranges and width for the current line.
*/

#define FLTH 100

static int StoreData()
{ int c, i, j, k, w, cnt, use; char f[FLTH], *nf; dec x;

  for(cnt=k=i=0; i<=cmax; i++)               //Advance to the beginning of the
  { c = advance(0); if(c!='|') cerr(840.2);  //next column and determine whether
    byp(0); use = clim[clabel[i]-'a']>0;     //it is to be used in the array.

    if(use && cindex[i]<0)                   //If this is a set of indexes only,
    { if(k>=MDIM) Error(526.);               //extract all ranges.
      ranges(clabel[i], drange[k],CRANGE);
      dwth[k++] = cwth[clabel[i]-'a']; }

    else if(use)                             //If it is a number to be spread
    { field(f, FLTH); if(f[0]==' ') break;   //about the array, extract its
      x = strtod(f, &nf);                    //value and rescale it as needed,
      if(fstr)  x = x*xm+xb;                 //optionally truncating to an
      if(trunc) x = (int)x;                  //integer.

      for(j=nf-f; f[j]==' '; j++);           //Make sure the field does not
      if(f[j]) cerr(533.1);                  //carry spurious characters.

      w  = cwth[clabel[i]-'a'];              //Determine the base offset in the
      w *= cindex[i];                        //array and distribute the element
      cnt += StoreElement(x, w); }           //through the array.

    else                                     //If it is not used in the array,
    { field(f, FLTH); if(f[0]==' ') break; } //skip over it.

    c = preview(0);                          //If another field follows
    if(c=='|') continue;                     //immediately, repeat.

    if(c!='\n' && c!='\r') cerr(533.2);      //If the line has ended normally,
    BypLine(); return cnt; }                 //advance to the next and return.

  cerr(532.2); return 0;                     //(Too many columns.)
}


/*----------------------------------------------------------------------------*
STORE RANGES

This routine processes the ranges in a single field, delimited by vertical bars,
or by an opening vertical bar and closed by the end of the line.

ENTRY: 'd' contains the column label, 'a'-'z'.
       's' points to an area to receive the ranges.
       'k' contains the length of the area, in integers.

EXIT:  's' contains the ranges from the next field of the input file, up to but
         not including the next vertical bar or up to the end of the line.
*/

static ranges(int d, int s[], int k)
{ int c, i, m, r0, r1;
  static int c1;

  if(d<'a'||d>'z') cerr(524.2);              //Check the label for validity.

  for(i=0; i<k; i+=2)                        //Retrieve the first value of a
  { r0 = posint(0); if(r0<0) cerr(534.1);    //pair and save it as if it is a
    s[i] = s[i+1] = r0;                      //single value.

    c = advance(0);                          //If it is the first of a pair,
    if(c=='~')                               //retrieve the second value of the
    { r1 = posint(0); if(r1<0) cerr(534.2);  //pair and store a range.
      if(r0<=r1) s[i+1] = r1;
      else       s[i]   = r1;
      c = advance(0); }

    m = clim[d-'a']-1;                       //Check the range and ignore any
    if(s[i+1]>m)                             //values that are outside the
    { if(!c1++) cerr(387.); s[i+1] = m; }    //array, giving a warning message
    if(s[i]>m) i -= 2;                       //just once.

    if(c==',') continue;                     //Continue if another range follows
    if(c==' ') { byp(0); c = advance(0); }   //or return if it does not.
    if(EOC) { s[i+2] = -1; retract(c); return; }

    cerr(533.3); }                           //(Spurious characters in field.)
  cerr(534.3);                               //(Too many fields.)
}


/*----------------------------------------------------------------------------*
STORE DATA ELEMENT IN ARRAY LOCATIONS

This routine spreads a single data value across all array locations that apply,
as defined by the ranges of indexes supplied on the input line.

ENTRY: 'x' contains the value to be stored.
       'k0' contains the offset to be added to the offsets of all other indexes.
       'drange' contains the sequence of index ranges.
       'cwth' contains the cumulative width of each index value in 'drange'.
       'imax' contains the maximum index value recorded in 'drange'.

EXIT:  'x' has been stored as specified by 'drange' and 'W'.
       'StoreElement' contains the number of array elements populated on this
         call.

Note: During processing, 'I[i]' defines which range being processed for
dimension 'i' and 'V' contains the current index for each range. The routine
starts at the innermost index and works outward, acting somewhat like a
multi-base counter with carries out of one digit into the higer-order digit. The
order of operations does not matter for correctness, but working from the inner
index out tends to move through arrays in the order they are arranged in memory.
With local memory caches, that can be faster than jumping around the array.
*/

static int StoreElement(dec x, int k0)
{ int i, k, m, w, cnt, V[MDIM], I[MDIM];

  m = imax>1? imax-2: 0;                     //Load the starting index for each
  for(i=m; i>=0; i--)                        //higher dimension.
  { V[i] = drange[i][0]; I[i] = 0; }

  for(cnt=0;;)
  { for(k=k0,i=m; i>=0; i--)                 //Compute the location of the next
      k += V[i]*dwth[i];                     //value and store it in the array.
    if(k<0) return 0;
    io.data[k] = x; cnt += 1;

    for(i=m; i>=0; i--)                      //Advance to the next index,
    { V[i] += 1;                             //carring to higher levels as each
      if(V[i]<=drange[i][I[i]+1]) break;     //lower level wraps around.
      I[i] += 2; V[i] = drange[i][I[i]];
      if(V[i]>=0) break;
      I[i] = 0;  V[i] = drange[i][I[i]]; }

    if(i<0) return cnt; }                    //Return when all values have been
}                                            //stored.


/*----------------------------------------------------------------------------*
EXTRACT FIELD

This routine merely collects characters from the front of the remaining input
file and stores them in a character area as a normal string. It stops collecting
characters when a vertical bar or end of line code appears.

ENTRY: 's' points to an area to receive the field.
       'k' contains the length of the area, in characters.
       'fp' points to the control structure for an open file.

EXIT:  's' contains the next field from the input file, up to but not including
         the next vertical bar or to the end of the line.
*/

static field(char *s, int k)
{ int i, c;

  for(i=0; i<k-1; i++)                       //Retreive the next character.
  { c = advance(0);

    if(EOC)                                  //If the field is exhasuted, close
    { retract(c); s[i] = 0; return; }        //the string and return.

    s[i] = c; }                              //Store the character and repeat.

  cerr(535.);                                //(Field too long.)
}


/*----------------------------------------------------------------------------*
EXTRACT POSITIVE INTEGER

This routine reads the next positive integer from the file. It first bypasses
any leading blanks. Then if the character is not a digit 0-9, it returns a
negative value, indicating that there is no positive integer at the beginning of
the file to be extracted. Otherwise it incorporates all consecutive digits into
the integer and bypasses any blanks trailing the integer.

ENTRY: 'fp' points to the control structure for an open file.

EXIT:  'posint' contains the integer that was at the current front of the file.
         If there was no such integer, 'posint' is negative.
*/

static int posint(int rtneof)
{ int n, d;

  byp(0); n = advance(0);                    //If the file is not positioned
  if(n<'0'||n>'9')                           //at a positive integer, return
  { retract(n); return -1; }                 //a negative value.
  n -= '0';

  while(1)                                   //If it is, read any remaining
  { d = advance(0); if(d<'0'||d>'9') break;  //digits to complete the integer.
    n = n*10+(d-'0'); }

  retract(d); byp(0);                        //Bypass any trailing blanks and
  return n;                                  //return the integer.
}


/*----------------------------------------------------------------------------*
BYPASS LINE

ENTRY: The file is positioned somewhere within a line.

EXIT:  'BypLine' contains a result code.
         0 The file is advanced so that the next character read will begin the
             next new line.
         1 The file is exhausted.
*/

static int BypLine()
{ int c;

  while(1)                                   //Advance to the next characters
  { c = getc(fp); if(c==EOC) return 1;       //and return at end of file or end
    if(c=='\n') return 0;                    //of line.

    if(c=='\r')                              //Keep skipping line returns until
    { c = getc(fp); if(c==EOC) return 1;     //end of file or end of line.
      if(c=='\r') continue;
      if(c=='\n') return 0;
      retract(c); return 0; }
  }
}


/*----------------------------------------------------------------------------*
CHARACTER LEVEL INPUT

The three routines below process the input file, which is treated as a string of
characters separated into lines by end-of-line codes. Lines are terminated
either by zero or more line returns ('\r'), a single line feed ('\n'), or line
returns followed by a line feed. Three different conventions are accommodated:

  \r    Old Apple and Old CPT convention
  \r\n  MsDos convention
  \n    Unix convention

ENTRY: 'rtneof' is set if the program should return to the caller if the end of
         the file is encountered during the operation. Otherwise a fatal error
         is processed.
       'pf' points to the structure of a file opened for input.
       'c' is the character to be retracted, where that is applicable.

EXIT:  The routines return the character at the front of the file.
*/

static int advance(int rtneof)               //Advance to the next character,
{                                            //removing it from the front of
  int c = getc(fp);                          //the remaining file.
  if(c==EOF && rtneof==0) cerr(536.1);
  return c;
}

static int retract(int c)                    //Retract the specified character,
{                                            //adding it back at the front of
  ungetc(c, fp); return c;                   //the remaining file.
}

static int preview(int rtneof)               //Preview what the next character
{                                            //will be without changing the
  int c = advance(rtneof); retract(c);       //state of the file.
  return c;
}

static int byp(int rtneof)                   //Bypass blanks at the front of the
{ int c;                                     //file.
  do c = advance(rtneof); while(c==' ');
  retract(c); return c;
}


/* CLARENCE LEHMAN, MARCH 2011.

This module was developed piecemeal over the course of week busy with other
things, learning what it should be as it developed. Its documentation should now
be unified. The documentation is in a state where it can be understood with a
little study, and is complete enough for us to remember what it does and write a
better version as time permits.
*/

#define TESTPROGRAM0
#ifdef  TESTPROGRAM

/*============================================================================*
TEST PROGRAM FOR FILEIO.C
*/

#include <stdlib.h>

dec K[4][7][9][4];
dec L[9];

struct IO Z[] =
{ { (dec*)L, {-'I',9}                      },
  { (dec*)K, {-'r',4,-'y',7,-'q',9,-'c',4} } };

main()
{ int i;

  ErrorInit();

  FileIO("fileio1.in",  Z[0], "r|");

  for(i=0; i<9; i++)
    printf("%d %f\n", i, L[i]);

  FileIO("fileio1.out", Z[0], "w|=%6.4g");
  system("cat fileio1.out");

  FileIO("fileio2.in",  Z[1], "r|");
  FileIO("fileio2.out", Z[1], "w|=%6.4g");
  system("cat fileio2.out");
}

#endif

// ONE-DIMENSIONAL ARRAY TESTING
// |i   |z      |v
// |0   |0      |-5
// |1   |0.4000 |-5
// |2   |0.8800 |-5
// |3,5 |0.8890 |-5
// |4   |0.9750 |-5
// |6   |0.9954 |-5
// |7   |0.9985 |-5
// |8   |1      |-5
// |9~10|1.1    |-5

// MULTI-DIMENSIONAL ARRAY TESTING
// |r   |y       |g    |q   |c0     |c1     |c2     |c3     |v
// |0~3 |0,2~4,6 | 1   |0   |0      |0      |0      |0      |-5
// |0,3 |0,2~4,6 | 2   |1   |0.8000 |0.6000 |0.4900 |0.4000 |-5
// |0,3 |0,2~4,6 | 3   |2   |0.9793 |0.9294 |0.8800 |0.8800 |-5
// |0,3 |0,2~4,6 | 4   |3,5 |0.9800 |0.9374 |0.8900 |0.8890 |-5
// |0,3 |0,2~4,6 | 5   |4   |0.9979 |0.9947 |0.9818 |0.9750 |-5
// |0,3 |0,2~4,6 | 6   |6   |1      |0.9998 |0.9979 |0.9954 |-5
// |0,3 |0,2~4,6 | 7   |7   |1      |1      |0.9996 |0.9985 |-5
// |0,3 |0,2~4,6 | 8   |8   |1      |1      |1      |1      |-5
//
// |1,2 |0,2~4,6 | 10  |1   |0.8500 |0.3500 |0.3500 |0.3500 |-5
// |1,2 |0,2~4,6 | 11  |2   |0.9750 |0.7000 |0.7000 |0.7000 |-5
// |1,2 |0,2~4,6 | 12  |3   |0.9850 |0.7040 |0.7040 |0.7040 |-5
// |1,2 |0,2~4,6 | 13  |4   |0.9900 |0.9900 |0.9900 |0.9900 |-5
// |1,2 |0,2~4,6 | 14  |5   |0.9985 |0.9980 |0.9980 |0.9980 |-5
// |1,2 |0,2~4,6 | 15  |6~9 |1      |1      |1      |1      |-5

