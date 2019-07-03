/*----------------------------------------------------------------------------*
EVENT SCHEDULING

This module implements an event scheduler/dispatcher, efficient even for queues
containing hundreds of millions of events.

This algorithm may be used for scheduling actual events, as in a real-time
operating system, or for simulated events. In the former case, it is
synchronized with a real-time clock. In the latter, the next event is dispatched
immediately and the simulated clock is set to the time of that event. The
simulation system thus jumps forward from one event to the next, skipping
intervening times and never wasting computer cycles when nothing is happening.

It is easy to understand the scheduling algorithm used here with a physical
analogy. Suppose we have half a million events to be scheduled over the next
five years, and that they appear more-or-less randomly throughout that period.
We want to be able to (1) schedule new events, (2) cancel existing events, and
(3) notify some dispatcher as the time for each event arrives.

Image a series of bins, one per minute for an entire year. The first bin
represents the first minute after midnight on New Year's Day, the second bin
represent the second minute, and so forth to the last bin, which represents the
last minute on December 31st of the year. That is 366 days x 24 hours/day x 60
minutes/hour = 527,040 bins total, each labeled with the month, day, hour, and
minute that it represents. With half a million events to be scheduled, that is
less than one per bin on average.

Also imagine that each event has a ticket with its event number and its
scheduled time, represented at least to the nearest second. Then here is the
procedure for handling events:

1. Scheduling a new event.
   a. Go to the bin representing the month, day, hour, and minute for the
      event, ignoring the year and the second.
   b. Drop the event's ticket into the bin.
   c. That took only a single operation.

2. Cancelling an existing event.
   a. Go to the bin representing the month, day, hour, and minute for the
      event, ignoring the year and the second.
   b. Remove all tickets from the bin and flip through them to find the ticket
      for the event in question.
   c. Destroy that ticket and place the others back in the bin.
   d. That took one operation per ticket in the bin. But on average there is
      only one ticket in the bin.

3. Dispatching the next event.
   a. Go to the bin representing the current time of day.
   b. Remove the tickets from the bin and sort them into chronological order.
   c. Place any tickets for future years back in the bin.
   d. Hand each ticket for this year to a dispatcher as its scheduled time
      arrives.
   e. If any tickets for the current bin arrive while the bin is being
      dispatched, put them in their proper position among the other tickets
      for dispatching now or in the future.
   f. That took one operation per ticket in the bin. But on average there is
      only one ticket in the bin.

The basic algorithm underlying these procedures was discovered in the early days
of computer science. Today it would be called a "hash coded search with modulus
key." It allows near-instantaneous access to any entry and keeps the entries in
order, modulo the number of entries in the list. The first paper clearly
describing this kind of search is by Arnold Dumey (1956), entitled "Indexing for
rapid random access memory systems," in Computers and Automation 6:6-9. It is
interesting that "rapid" in that day meant something like 1/3 second per access!

The running time of the above three procedures, remarkably, is not a function of
the number of pending events! It is only a function of the ratio between the
number of pending events and the number of bins. Therefore, with sufficient
memory, the procedures are maximally efficient. In computer parlance, they are
called "Order 1", or "O(1)".

If there are a hundred million events pending randomly distributed into a
hundred million bins, for example, the average number in each bin is 1 and the
average in each occupied bin is $1/(1-1/e) = 1.58$. Therefore, essentially no
searching/sorting occurs, making the method exceedingly fast.

It also has modest space requirements, with only two additional integers needed
per event (one for a head-of-list index and another for a forward list index).
Therefore, assuming four bytes per integer, managing 100 million events requires
less than 3/4 gigabyte of memory, well within the reach even of laptop
computers. ($8\times10^8/2^{30}$ GB).

Assuming the events occur at random times, then a portion $1/e = 0.36787944$ of
the bins will be empty, the same portion will contain a single entry, $1/2e^2 =
0.18393972$ will contain two entries, and so on according to the Poisson
distribution. Below is an actual example with one million bins and one million
pending events. The expected numbers are Poisson values. The observed numbers
are from the end of an SIR disease run, after many millions of events had been
added, cancelled, and dispatched. There is no significant clustering to slow
things down. Here practice corresponds well with theory.

      N  Observed   Expected

      0    371047     367879
      1    364903     367879
      2    181992     183940
      3     62032      61313
      4     15998      15328
      5      3256       3066
      6       675        511
      7        82         73
      8        14          9
      9         1          1

As implemented in this module, The bins need not correspond to standard time
units such as a minutes. The bin representing the current time is kept sorted so
that each event to be dispatched need only be picked from the front of the list.
That means that any time a new event is added to the current bin, that bin must
be resorted.

In this module the three functions that correspond to the above are
'EventSchedule', 'EventCancel,' and 'EventNext'.  They work essentially as
described above.
*/


/*----------------------------------------------------------------------------*
DATA STRUCTURES

The module has two main data structures, a circular series of time bins 'Q[h]',
which index one event scheduled for the time slots 'h' represented by each bin,
and a series of singly linked lists 'P[i]', which define the earliest event for
each individual, linked to their proper time slots. Note that each bin 'Q[h]'
represents many related time slots, all equal modulo the width of the series of
time bins, 'Qw', described further below.

The number of links 'P[i]' must be equal to the number of individuals, and is
indexed by individual number. but the number of time bins 'Q[h]' may be smaller
or larger than the number of individuals. The size of 'Q' is a matter of
optimization. It is typical to have one time bin for each event that could be
scheduled, meaning each bin will represent a single event on average. That is
optimal if all bin operations are of equal speed.

The width 'Qw' of all bins combined is also a matter of optimization. If it is
much too large, events will tend to cluster near the bin being dispatched. If it
is much too small, events will tend to spread out, with most bins containing
events that are for the more distant future, not immediate events. A suitable
value for 'Qw' can be found by knowledge of the system being analysed, or by
experimental trials to find a good speed of operation.

The linked lists of events in 'P' are not maintained in any particular order,
for speed of addition, but before a bin is dispatched, it is sorted into order
by event times, for speed of dispatching.
*/

#include "common.h"

#define PINIT  if(run1==0) EventInit();

#define PEMPTY -1      //Marker for bins containing no linkages.
#define TN (INDIV+0)   //Maximum number of time bins.
#define PN (INDIV+3)   //Maximum number of time bin forward indexes.
#define TW  20         //Time width of all bins combined (for optimization).

dec t;                 //Current time, last dispatched event.

static int run1;       //Flag to detect if the routine is being reused.

static dec T[PN];      //Time for each scheduled event.
static int P[PN];      //Forward indexes within bins, ending with zero.
static int Q[TN];      //First index for the bin, with zero for empty bins.

static int Qn  = TN;   //Number of elements in 'Q'.
static dec Qw  = TW;   //Interval of time represented for each cycle in 'Q'.
static int Qi  = 0;    //Index of the immediate time bin.
static int Qo  = 1;    //Flag set if the immediate bin is in order.
static int Qe  = 0;    //Number of events in all bins.

static dec Qt0 = 0;    //Earliest time representable this cycle in 'Q'.
static dec Qt1 = TW;   //Earliest time beyond this cycle in 'Q'.

/*----------------------------------------------------------------------------*
INITIALIZE STATIC DATA STRUCTURES

This routine must be called if an entirely new simulation is to be started in
the absence of the program being restarted from the beginning. It makes the
module serially reusable, and is included to cover a deficiency in MPI whereby
'system' and 'popen' corrupt the system.

ENTRY: It is time to start an entirely new run.

EXIT:  Any prior data have been wiped clean.
*/

EventInit()
{ int i;

  for(i=0; i<PN; i++) P[i] = PEMPTY;
  if(run1==0) { run1 = 1; return; }

  for(i=0; i<PN; i++) T[i] = 0;
  for(i=0; i<TN; i++) Q[i] = 0;

  Qn  = TN; Qw  = TW; Qi  = 0;
  Qo  = 1;  Qe  = 0;
  Qt0 = 0;  Qt1 = TW;

  t = 0;
}


/*----------------------------------------------------------------------------*
SET STARTING TIME

If the starting time is much greater than zero, this routine may be called at
the outset so the list of events will start at that time. It is not necessary
that it be called but it can save a good deal of time working forward to the
present if the present is some real-year value, such as 1976.

ENTRY: 't0' contains the time at or before which the first event will occur.
       No other routines in this module have been called.

EXIT:  't' is set to the 't0'.
        Data structures have been set to start at time 't0'.
*/

EventStartTime(dec t0)
{ PINIT;                                     //Initialize if necessary.

  if(Qe) Error(742.);                        //Make sure the bins are empty.

  Qt0 = t0-(Qw/Qn)/2;                        //Set the time boundaries, leaving
  Qt1 = Qt0+Qw;                              //room for rounding errors.

  t = t0;                                    //Set the global time.
}

/*
Note: The starting time is positioned in the middle of the first bin because
of the limits of finite precision arithmetic. Suppose it were not and the
starting time was set precisely to the year 1960, but an independent
calculation somewhere generated the first event as 1959.9999999999, or
equivalent in some other radix. Then the event would not be slotted in the
first bin, but rather in the last, and would not come out of the list in
proper order. Positioning it in the middle prevents that.
*/


/*----------------------------------------------------------------------------*
SCHEDULE NEW EVENT

ENTRY: 'n' contains the number (starting with 1) of a new event.
       'te' contains the time at which the event will occur.
       'P[n]' indicates that the event is unscheduled (equal to 'PEMPTY').

EXIT:  The event has been scheduled, to occur when the proper time arrives.
       'T[n]' records the time 'te' of the event.
       'P[n]' links the event with others in its time bin.

WORK:  The scheduling data structures are prepared as described above.
*/

EventSchedule(int n, dec te)
{ int i; dec tr;

  PINIT;                                     //Initialize if necessary.

  if(n<1||n>=PN)   Error1(734.1, "n=",n);    //Check the index and make sure an
  if(P[n]!=PEMPTY) Error1(735.1, "n=",n);    //event is not already scheduled
  if(te<t) Error2(737., "t=",t, ">",te);     //and is not in the past.

  T[n] = te;                                 //Record the time of the new event.

  tr = (te-Qt0)/Qw; tr -= (int)tr;           //Convert the time to a bin number
  i  = tr*Qn; if(i==Qi) Qo = 0;              //and mark for sorting if needed.

  P[n] = Q[i]; Q[i] = n;                     //Add the event to the list for
  Qe += 1;                                   //that bin and increment the number
}                                            //of events.


/*----------------------------------------------------------------------------*
CANCEL EXISTING EVENT

ENTRY: 'n' contains the number (starting with 1) of the event to be cancelled.
       'T[n]' contains the scheduled time of the event.

EXIT:  The event has been removed from the list.

WORK:  The scheduling data structures are prepared as described above.
*/

EventCancel(int n)
{ int i, j, jp; dec tr;

  PINIT;                                     //Initialize if necessary.

  if(n<1||n>=PN)   Error1(734.2, "n=",n);    //Check the index and make sure an
  if(P[n]==PEMPTY) Error1(736.2, "n=",n);    //event is scheduled.

  tr = (T[n]-Qt0)/Qw; tr -= (int)tr;         //Convert the time to a bin number,
  i  = tr*Qn;                                //modulo the duration of the cycle.

  if(cancel1(n, i)) return;                  //Remove it from its normal bin.

  i = (i-1+Qn) % Qn;                         //If it is not there, check the
  if(cancel1(n, i)) return;                  //bin below (from rounding error).

  i = (i+2+Qn) % Qn;                         //Finally, check the bin above.
  if(cancel1(n, i)) return;

  Error1(818., "n=",n);                      //If the specified event was not in
}                                            //the list, signal an error.

int cancel1(int n, int i)
{ int j, jp;

  for(jp=0,j=Q[i]; j>0; jp=j,j=P[j])         //Scan the list of pending events
    if(j==n)                                 //in this bin and remove the
    { if(jp>0) P[jp] = P[j];                 //specified event. (The average
      else     Q[i]  = P[j];                 //number events in a non-empty bin
      P[j] = PEMPTY;                         //is 1.5)
      Qe -= 1; if(Qe<0)
        Error2(819., "n=",n, " bin=",i);
      return 1; }
  return 0;
}

/*
Note: It is necessary to check not just the expected bin when an event is to be
deleted, but the two adjacent bins as well if the event does not appear. This is
because of the vagaries of finite-precision computer arithmetic. Two bins are
separated by a knife-edge, and what is calculated one time as 1.0 may be
calculated in a slightly different manner as 0.999999999999999, dropping it in
an adjacent bin. In simulations testing this condition, it happened only once
per ten million times or less, so there is no speed issue here, only a
correctness issue.
*/


/*----------------------------------------------------------------------------*
RENUMBER EXISTING EVENT

This routine renumbers events. It can be called, for example, to reuse an entry
when it becomes available, for example from a simulated death, shifting an
existing individual from the end of the array keep the array compact.

ENTRY: 'n' contains the new index number, which has no event scheduled.
       'm' contains the current index number of the event.
       There is an event scheduled for 'm'.

EXIT:  'n' is the new index number of the individual.
       The event originally scheduled as 'm' is re-scheduled as 'n'. Event 'm'
         no longer has an event scheduled and the index is free to be reused.
*/

EventRenumber(int n, int m)
{
  if(n<1||n>=PN) Error1(734.3, "n=",n);      //Check the indexes and make sure
  if(m<1||m>=PN) Error1(734.4, "n=",m);      //they are in range.

  if(n!=m)
  { T[n] = T[m];                             //Transfer the time.
    EventCancel(m);                          //Cancel the old number.
    EventSchedule(n, T[n]); }                //Reschedule as the new number.
}

/*----------------------------------------------------------------------------*
LOCATE NEXT EVENT

ENTRY: 'T' contains the time of the next event for every event scheduled.
       The bin structure is properly initialized.

EXIT:  'EventNext' contains the number of the next event. If zero, no events
         are scheduled.
       't' contains the time of the next event, if 'NextEvent' is not zero.

WORK:  The scheduling data structures are prepared as described above.
*/

int EventNext()
{ int j, n;

  PINIT;                                     //Initialize if necessary.

  while(Qe>0)
  { for(; Qi<Qn; Qo=0,Qi++)                  //Advance to the next non-empty
    { j = Q[Qi]; if(j==0) continue;          //bin.

      if(Qo==0)                              //Sort the bin if it may be
      { j = Q[Qi] = sort(P,j,0); Qo = 1; }   //necessary.

      if(T[j]<Qt1)                           //If the event belongs to this
      { if(P[j]==PEMPTY) Error(820.1);       //pass, remove it from the list,
        Q[Qi] = P[j];                        //decrement the number of events,
        P[j] = PEMPTY; Qe -= 1;              //advance the global time, and
        t = T[j]; return j; } }              //return the event's index.

    Qi = 0; Qt0 += Qw; Qt1 = Qt0+Qw; }       //Circle back to the first bin.

  return 0;                                  //Signal completion of all events.
}

/*
NOTE: It is important to sort the active bin, as is done here. At first thought,
it might simpler and faster just to exhaustively searched the bin each time,
since there are typically only a few items in each bin. But there need not be
only a few. Suppose a million individuals ended up in one bin, say due to some
common start-up condition. Sorting is $n\log n$, exhaustively searching would be
$n^2/2$, and the latter would be prohibitive.
*/


/*----------------------------------------------------------------------------*
DISPLAY PROFILE

This routine displays the distribution of schedule events---how frequently 0, 1,
2, 3, ... events occur in bins, and compares it with a Poisson distribution,
$exp(-x) x^n / n!$

ENTRY: 'label' contains an optional label to be displayed with the results,
         such as "Initial" or "Final".
       The bins have been initialized.

EXIT;  'EventProfile' contains the amount of memory allocated for the main data
        structures.
*/

#define PROF 1001

int EventProfile(char *label)
{ int i, j, n, imax, prof[PROF];
  dec nexp, lambda, eml, ln, nf;

  if(label==0||label[0]==0) label = "Bin";   //Establish a default label.

  for(i=0; i<PROF; i++) prof[i] = 0;         //Clear the profile array.

  for(i=0; i<Qn; i++)                        //Count the number of bins that
  { for(j=Q[i],n=0; j>0; j=P[j],n++)         //have no entries, one entry, two
      if(j<1||j>=PN||n>PN)                   //entries, etc.
        Error1(820.2, "j=",j);
    if(n>PROF-1) n = PROF-1;
    prof[n] += 1; }

  for(i=imax=0; i<PROF; i++)                 //Determine the greatest number
    if(prof[i]) imax = i;                    //in any bin.

  lambda = (dec)Qe/Qn;                       //Compute the parameters for the
  eml = exp(-lambda); ln = nf = 1;           //Poisson distribution.

  printf("%s distribution of %d events:\n",  //Display heading lines.
    label, Qe);
  printf("   N   Observed   Expected\n");

  for(i=0; i<=imax; i++)                     //Display the distribution of bin
  { nexp = Qn*eml*ln/nf;                     //sizes, compared with the complete
    if(prof[i] || nexp>0.5)                  //randomness of the Poisson.
      printf("%4d%c%9d %10.0f\n",
        i, i<PROF-1?' ':'+', prof[i], nexp);
    ln *= lambda; nf *= i+1; }

  printf("\n");                              //Leave a blank line and return
  return sizeof Q + sizeof T + sizeof P;     //with the size of the main data
}                                            //structure.


/*----------------------------------------------------------------------------*
DETERMINE SORTING ORDER

This routine is called from 'sort' when it needs to have two elements compared.

ENTRY: 'p' indexes list element one.
       'q' indexes list element two.

EXIT:  'order' contains a status code.
        '<0' Element one is less than element two.
        '=0' Element one is equal to element two.
        '>0' Element one is greater than element two.
*/

int order(int p, int q)
{ dec w;

  w = T[p]-T[q]; return w<0? -1: w>0? 1: 0;
}

/* (C) CLARENCE LEHMAN, UNIVERSITY OF MINNESOTA, AUGUST 2009.

[To become open-source software as part of the modelling package.]

Modifications

1. Data structures 'P' and 'T' removed from the caller's list of individuals and
   placed under the management of this module. February 2011 [CLL].

2. Safeguards added so the same event number cannot be scheduled in duplicate.
   February 2011 [CLL].

3. Comments and names updated for general distribution, April 2011 [CLL].

4. 'EventInit' added for serial reusability, May 2011 [CLL].
*/

