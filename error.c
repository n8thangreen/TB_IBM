/*----------------------------------------------------------------------------*
MESSAGES

All messages are identified by a three-digit number and may also be described in
words. The wording is collected here so it can be unified in style and structure
across the system and changed to languages other than English. The numbers are
assigned as follows:

   000-099  Not to be used
   100-199  Status messages (e.g., periodically as the program runs)
   200-299  Information messages (no operator action needed)
   300-499  Warning messages (operator action may be needed)
   500-799  Fatal error messages (likely induced by data)
   800-999  Fatal error messages (likely induced by program)

This version of the routine displays the messages on the standard output stream.
Other versions can use message boxes on the screen, notification of other
processes running in the system, and so forth.
*/

static char *message[] =
{
  "F387%s  Warning: An index value is out of range and has been ignored",

  "F501%s  This feature is not yet supported",
  "F510%s  The file cannot be opened",
  "F511%s  The file cannot be completely read",
  "F512%s  The file cannot be completely written",
  "F513%s  The file is not of valid format",
  "F514%s  The file appears to be without data",
  "F515%s  A file I/O index is incorrect",
  "F516%s  A file I/O dimension is not positive",
  "F517%s  A file I/O index is too large",
  "F518%s  A file I/O index range is not divisible by the increment",
  "F520%s  The file I/O transformation does not begin with 'x' or 'n'",
  "F521%s  The file I/0 transformation would divide by zero",
  "F522%s  The file I/0 transformation is syntactically incorrect",
  "F523%s  A file I/0 column heading is incorrect",
  "F524%s  A file I/0 index label is incorrect [must be 'a'-'z']",
  "E525%s  The parameter is incorrect",
  "F526%s  A file I/0 dimension is too large",
  "F530%s  The file I/0 structure must have all indices first",
  "F531%s  The file I/0 structure must have at least one data column",
  "F532%s  The file has too many columns",
  "F533%s  A file I/O field contains spurious characters",
  "F534%s  A file I/0 index field is incorrect",
  "F535%s  A file I/0 field is too large",
  "F536%s  The file ended prematurely",

  "E609%s  The state is out of range",
  "E610%s  The number of individuals is incorrect",
  "E611%s  The vaccination type is incorrect",
  "E612%s  Time of death is in error",
  "E616%s  The strain ID is incorrect",
  "E617%s  The time since infection is incorrect",
  "E618%s  Sorting of times is in error",
  "E619%s  Time of reporting is zero",
  "E620%s  Time to disease is in error",
  "E621%s  A cumulative table is not monotonically increasing",
  "E622%s  A cumulative table is not bounded by 0 and 1",

  "E714%s  The kernel designator is incorrect",
  "E734%s  The event number is out of range",
  "E735%s  An event to be scheduled is already scheduled",
  "E736%s  An event to be cancelled is not yet scheduled",
  "E737%s  A new event would be scheduled in the past",
  "E742%s  Attempt to initialize when the time bins are not empty",
  "E753%s  A binary search table is invalid",
  "E754%s  A cumulative table has gone beyond 1",

  "F818%s  An existing event cannot be found in the time bins",
  "F819%s  The event counter has fallen negative",
  "F820%s  The event list has a broken link",
  "F840%s  An internal inconsistency during file I/O has been detected",
  "F850%s  A birth occurred before the present",

  "F911%s  Not enough memory is available",
  "F920%s  An index is out of range",
  "F921%s  A pointer is null",
  "F922%s  A switch index is incorrect",

  "F996%s  System bus error",
  "F997%s  System segmentation error",
  "F998%s  Unsupported error number",
   0 };

static char *fixed[] =
{ "%c%s  %s defined in the source code",
  "F999%s  Processing cannot continue",
   0 };

/*----------------------------------------------------------------------------*
ISSUE MESSAGE

When an error arises in the system, or when an information message is needed,
this routine is called. Given a message number, it scans tables in this module
for a message that matches the number, then displays the message. If no
corresponding message exists, a standard message is given. For example, if
message number 555 is not in the message table, the following pair of messages
will be issued.

  E555  Error defined in the source code.
  E999  Processing cannot continue.

The second is because the message number is greater than 500, which indicates a
fatal condition. Were the message number smaller, say 333, processing would
continue after the following warning message.

  W333  Warning defined in the source code.

Parameters may be displayed in parentheses following the message to help
determine what is wrong. These parameters arrive in string-number pairs. One of
the functions 'Error', 'Error1', 'Error2', or 'Error3' is called, depending on
how many string-number pairs are being supplied. The pairs follow these rules:

1. If the string ends with a continuation code ('<', '=', '>', or ':"), both the
string and the number are displayed next to each other.

2. If the string begins with a special marker ('`'), both the string and the
number are both displayed as in 1 above, except the special marker is bypassed.

3. If the string does not match 1 or 2 above, only the string is displayed, the
number is bypassed.

For example, here are some sample calls and the message they produce.

  Error(511.);                        F501  The file cannot be read.
  Error1(511., "k=",3.14);            F501  The file cannot be read (k=3.14).
  Error2(511., "w.txt",.0, " k=",3.); F501  The file cannot be read (w.txt k=3).
  Error1(511., "`Space ",1.5);        F501  The file cannot be read (Space 1.5).

A single decimal point may be appended to differentiate places in the code that
the same message occurs. For example:

  Error(501.2);                       F501.2  The file cannot be read.

ENTRY: 'n' contains the message number, assigned as described above.
       'p1' and 'v1', if applicable, hold a character string and numeric value
         to be appended to the message. If 'p1' is null, both the string and the
         number are considered null.
       'p2' and 'v2', if applicable, hold a second character string and numeric
         value to be appended, under the same conditions as 'p1' and 'v1'.
       'p3' and 'v3', if applicable, hold a second character string and numeric
         value to be appended, under the same conditions as 'p1' and 'v1'.

EXIT:  The message has been issued.  If the message number on entry was greater
         than or equal to 500, the function never returns, but calls a routine
         named 'Failure' instead. This failure routine must not return but
         rather must abort the program. (For example, in command-line mode,
         'Failure' can use the command 'exit(3)' to return to the operating
         system.)
*/

typedef double dec;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LINE 1024
#define LNUM   32
#define NPAR    3

static int dropnum(char*);
static int matched(char*, char*);
static msgout(int, char*, char *pc[NPAR], char pn[NPAR][LNUM]);
static nformat(char[], dec);

int Error3(dec n, char *p1,dec v1, char *p2,dec v2, char *p3,dec v3)
{
  int i; char *s, *pc[NPAR], pn[NPAR][LNUM], num[LNUM], line[LINE];

  if(n<100||n>=999.5)                        //Make sure the message number is
  { p1 = "`"; v1 = n; n = 998; p2 = 0; }     //within bounds.

  for(i=0; i<NPAR; i++)                      //Make the array of parameters
  { pc[i] = 0; pn[i][0] = 0; }               //null.

  sprintf(num, "%03.1f", n);                 //Format the message number and
  if(num[4]=='0') num[3] = 0;                //clear any trailing '.0'.

  pc[0] = p1; pc[1] = p2; pc[2] = p3;        //Format the parameter strings and
  if(pc[0]) nformat(pn[0], v1);              //numbers in the array.
  if(pc[1]) nformat(pn[1], v2);
  if(pc[2]) nformat(pn[2], v3);

  for(i=0; message[i]; i++)                  //Scan the table for a message
    if(matched(message[i], num)) break;      //which has a matching number.

  if(message[i])                             //If there is one, format it with
  { sprintf(line, message[i], &num[3]);      //any parameters supplied.
    msgout(2, line, pc,pn); }

  else                                       //If none can be found, issue a
  { s = n<200? "Status":                     //general purpose message.
        n<300? "Information":
        n<500? "Warning": "Failure";
    sprintf(line, fixed[0], s[0], num, s);
    msgout(2, line, pc,pn); }

  if(n>=500)                                 //Either return to the caller or
  { s = num[3]? "  ": "";                    //display the stack and abort the
    sprintf(line, fixed[1], s);              //entire program, depending on the
    pc[0] = 0; msgout(0, line, pc,pn);       //error code.
    fprintf(stderr,"\n"); ErrorTrace();
    while(1) Failure(3); }
  return 0;
}

int Error2(dec n, char *p1,dec v1,
                  char *p2,dec v2) { return Error3(n, p1,v1, p2,v2, 0,0); }
int Error1(dec n, char *p1,dec v1) { return Error3(n, p1,v1, 0,0,   0,0); }
int Error (dec n)                  { return Error3(n, 0,0,   0,0,   0,0); }

/*----------------------------------------------------------------------------*
INTERNAL SERVICE ROUTINES

/* ---- CHECK FOR MATCHING MESSAGE

ENTRY: 'p' points to a message.
       'q' points to a message number, formatted as a three-digit decimal
         character string, optionally with a modifier after a decimal point.

EXIT:  'matched' defines the result.
         0  The message does not match the message number.
         1  The message matches.
*/

static int matched(char *p, char *q)
{
  return !strncmp(p+1, q, 3);
}

/* ---- DISPLAY TEXT

ENTRY: 'r' defines the number of new-line codes to precede the message.
       's' points to the main message to be displayed.
       'pc' and 'pn' contain the character and numeric parameters, if any, to be
         displayed in pairs. They are coded as described at the front of this
         module.

EXIT:  The message has been displayed.
*/

static msgout(int r, char *s, char *pc[NPAR], char pn[NPAR][LNUM])
{
  char i, *s1, *s2;

  for(; r>0; r--) fprintf(stderr,"\n");      //Display any separator lines.

  for(i=0; i<NPAR && pc[i]; i++)             //If a text string is not marked,
    if(pc[i][0]=='`')  pc[i] += 1;           //eliminate the number following.
    else if(dropnum(pc[i])) pn[i][0] = 0;

  s1 = s2 = "";                              //Setup parentheses to enclose any
  if(pc[0]) { s1 = " (";  s2 = ")"; }        //parameters.

  fprintf(stderr,"%s%s", s, s1);             //Display the parameters in
  for(i=0; i<NPAR && pc[i]; i++)             //parentheses, if any parameters
    fprintf(stderr,"%s%s", pc[i], pn[i]);    //were supplied.
  fprintf(stderr,"%s.\n", s2);

  fflush(stdout);                            //Make sure the message is visible.
}

/* ---- CHECK FOR SPECIAL MARKS

ENTRY:  's' points to a character string. It is coded as described at the front
          of this module for the first of a pair of character-number pairs.

EXIT:   'dropnum' is set if the string does not end with a continuation symbol
          ('<', "=", '>', or ':').
*/

static int dropnum(char *s)
{ int c;

  if(s==0 || s[0]==0) return 1;              //Check for null strings.

  c = s[strlen(s)-1];                        //If the string ends with one of
  if(c=='<'||c=='='||c=='>'||c==':')         //the continuation symbols, allow
    return 0;                                //the number to be displayed.

  return 1;                                  //Otherwise suppress the number.
}

/* ---- FORMAT NUMBER

ENTRY: 's' is an area to receive the formatted number.
       'v' contains the number.

EXIT:  's' contains the formatted number.
*/

#define EPSILON 0.000000001
#define OMEGA   10000000000

static nformat(char s[LNUM], dec v)
{
  if(fabs(v-(long int)v)<fabs(v)*EPSILON     //If this is a relatively small
  && fabs(v)<OMEGA)                          //integer, format it directly.
    sprintf(s, "%d", (int)v);

  else                                       //Otherwise allow the use of
    sprintf(s, "%g", v);                     //exponential notation.
}

/* ---- FAILURE OF THE PROGRAM

ENTRY: 'n' contains an exit code, to be returned to the operating system.

EXIT:  The routine never exits, but rather terminates the program.
*/

Failure(int n)
{
  exit(n);
}

/*----------------------------------------------------------------------------*
STACK TRACE

Routines in this section provide a rudimentary trace of the stack when an error
occurs. For full debugging information, the code may be compiled with options
'-gdwarf-2 -g3 -rdynamic'.

Function 'ErrorInit' should be called at the beginning of processing to capture
errors such as invalid pointers. Function 'ErrorTrace' is called internally to
display the stack, showing the path that led to the error. It is only called for
fatal errors.
*/

#include <signal.h>
#include <execinfo.h>
#define  PSTACK 50

static __sighandler_t ErrBus() { Error(996.); return 0; }
static __sighandler_t ErrSeg() { Error(997.); return 0; }

int ErrorInit()
{ static int init;
  if(init) return 0;

  signal(SIGBUS,  (__sighandler_t)ErrBus);   //Capture bus errors.
  signal(SIGSEGV, (__sighandler_t)ErrSeg);   //Capture segmentation errors.

  init = 1; return 0;
}

int ErrorTrace()
{ int n; void *pstack[PSTACK];

  n = backtrace(pstack, PSTACK);
  backtrace_symbols_fd(pstack, n, 2);
  return 0;
}

/*
CLARENCE LEHMAN, APRIL 1988.

1. Adapted from the hardware compiler and expanded to floating point message
   numbers, November 2004 [CLL].

2. Converted to ANSI C format and adapted to new parameter format,
   February 2011 [CLL].

3. Error messages to 'stderr' rather than 'stdout', April 2011 [CLL].
*/

