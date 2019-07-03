/*----------------------------------------------------------------------------*
SORTING

Lists of data in main storage can be lexicographically sorted by this function,
the precise ordering being determined by another function whose address is
supplied as a parameter.  At the beginning of each list element is a forward
index, which defines the next element in the list.  A zero index terminates the
list.

On entry, the indexes may link the data in any order.  On exit, they have been
modified so that they link the data in proper lexicographic order. The sort is
stable, which means that elements which are lexicographically equal retain the
order they had originally.  No data in main storage are actually moved.

The routine is quite fast, requiring at most about 'n log2(n)' comparisons,
where 'n' is the number of elements in the list.  This is near the theoretical
maximum speed for any sort algorithm that must use pairwise comparisons.
Presequencing reduces the number of comparisons required.  Even if the list is
in random order, there is some presequencing.  For example, if we select two
adjacent entries in a random list, the chance is better than even that they will
already be in order (equal elements are considered to be in order).  If the list
is fully ordered on entry, only 'n-1' comparisons are required. Detailed
mathematical aspects of the function's behavior have not been worked out for
general cases.

The main routine handles a few special cases itself---simple ones that occur
frequently in certain applications---then passes the remainder through the full
sorting algorithm.

ENTRY: 'list' points to an array of forward indexes. 'list[0]' is unused.
       'p' indexes the first element of the list, which ends with a zero.
       'n' contains the number of items in the list, if known. If zero, the
         number of items is not known and 'sort' should count.
       'order' is a function to compare two list elements.  Its entry and
         exit conditions are given below.

EXIT:  'sort' indexes the first element in the sorted list, which ends with a
        zero.
*/

#include "common.h"

static int *P, pcurr, pprev, count;

int sort(int list[], int p, int n)
{ int i;

  P = list;                                  //Record the location of the list.

  if(n==0) for(i=p; i; i=P[i],n++);          //Count the number of elements and
  if(n==0||p==0) return 0;                   //return empty lists immediately.

  if(n==1) return p;                         //Do the same for single elements.

  if(n==2)                                   //If the list contains only two
  { if(order(p, P[p])<=0) return p;          //elements, sort it by inspection
    i = P[p]; P[i] = p; P[p] = 0;            //(One and two element lists are
    return i; }                              //the most common in applications
                                             //like hashing.)

  pcurr = p;                                 //Otherwise sort the list with the
  return isort(n);                           //full algorithm.
}

/*----------------------------------------------------------------------------*
ENTRY/EXIT CONDITIONS FOR THE ORDERING SUBROUTINE

The lexicographic ordering routine must be compatible with the following
definition:

ENTRY: 'p' indexes list element one.
       'q' indexes list element two.

EXIT:  'order' contains a status code.
        '<0' Element one is less than element two.
        '=0' Element one is equal to element two.
        '>0' Element one is greater than element two.

int order(int p, int q)
*/

/*----------------------------------------------------------------------------*
RECURSIVE MERGE-SORT

This is a recursive merge-sort.  To understand intuitively how it works, think
about writing a new sort routine in terms of an existing one.  Postulate that we
already have two routines, one which accepts a list of data and sorts it into
order, another which accepts two lists already in order and merges them into a
single list, also in order.  Now ask the question, given these two postulated
routines, how could we write a new sort routine?  An appropriate algorithm
requires only two major steps:  (1) If the number of elements in the list to be
sorted is precisely one, then return immediately, since a list containing a
single element is already sorted.  (2) If the number of elements in the list to
be sorted is greater than one, then divide the list in half (as closely as
integer arithmetic permits), call the existing sort routine to sort the first
half, do the same to sort the second half, then call the merge routine to
combine the two sorted halves.  Provided that the postulated sort and merge
routines work, it is clear that the new sort routine also works.

So using an existing sort routine, we have easily defined another.  This is not
surprising.  But now we note that the new sort routine satisfies all the
requirements of the postulated one, and therefore could be called in its place.
That is, instead of calling the postulated sort routine, we substitute recursive
calls to our newly defined routine.  Through the magic of recursion, we obtain a
sort routine merely by postulating the existence of one!  It remains to write
the postulated merge routine, but this task is elementary.

The appearance of recursion in a non-recursive data structure may be surprising,
but what is more surprising is that such a simple sort algorithm runs as fast as
the best known sort algorithms, and far faster than many algorithms in popular
use.  The algorithm used below is somewhat more advanced than the one just
described, but all the steps described can be found in the code that follows.
The departure consists in the method used to take advantage of presequencing.

A different approach to recursive merge-sorting is given by C. Bron
(Communications of the ACM, 5.15, May, 1972.  His method can be thought of as
sorting from the bottom up, whereas this one sorts from the top down. Moreover,
that method, at least in the form presented, destroys the ordering of elements
with equal keys and cannot take advantage of presequencing.

ENTRY: 'n' defines the minimum number of elements to be sorted.
       'P' points to the list of forward indexes.
       'pcurr' indexes the first element of the list.
       'order' is a function to compare two list elements.

EXIT:  'isort' indexes the first element in the sorted list.  The list ends
         with a zero index.
       'count' defines the number of elements which were actually sorted.  This
         can be greater than or equal to its value on entry.
       'pcurr' indexes the element following the last element sorted. If the
         entire list has been sorted, 'pcurr' is null.
*/

int isort(int n)
{
  int wp1, wp2, count1;

// SORT A SINGLE ELEMENT

  if(n<=1)                                   //If a single element has been
  { if(pcurr==0) return 0;                   //requested, initialize variables
    wp1 = pcurr; count = 0;                  //and check for error in count.

    do                                       //Scan forward in the list to
    { pprev = pcurr; count += 1;             //find the longest string that
      if ((pcurr=P[pcurr])==0) return wp1; } //is already in order.
    while(order(pprev, pcurr)<=0);

    P[pprev] = 0;                            //Break this string from the
    return wp1; }                            //rest of the list and return.

// SORT MORE THAN ONE ELEMENT

  wp1 = isort(n/2);                          //Sort the first part of the list
  if(n<=count) return wp1;                   //and return if the required
  count1 = count;                            //amount was fortuitously sorted.

  wp2 = isort(n-count);                      //If it was not, then sort what
  count += count1;                           //remains to be sorted on this
  return imerge(wp1, wp2);                   //call and merge the two sublists.
}

/*----------------------------------------------------------------------------*
MERGE TWO SORTED LISTS

This function accepts two lists already sorted and merges them together to form
a single sorted list.  One of the lists is designated as primary.  When a
merging ambiguity occurs (that is, when the same key appears in both lists),
elements from the primary list are merged first.

Both lists end with zero index values.  Merging is done by manipulation of index
values; no data are actually moved.

ENTRY: 'P' points to the list of forward indexes.
       'p' indexes the first element of a sorted primary list.
       'q' indexes the first element of a sorted secondary list.
       'order' is a function to compare two list elements.

EXIT:  'imerge' indexes the first element of the merged list.
*/

int imerge(int p, int q)
 {
 int pbeg;

 if(p==0) return q;                          //Handle cases where one or more
 if(q==0) return p;                          //of the lists is null.

 pbeg = p;                                   //Save an index to the beginning
 if(order(p,q)<=0) goto scanp;               //of the list and enter the
 pbeg = q;                                   //proper scanning routine.

 scans:                                      //Scan for a secondary element
   do { pprev = q;                           //greater than or equal to the
        if ((q=P[q])==0) goto ends; }        //current primary element and mend
   while(order(p,q)>0);                      //the secondary list.
   P[pprev] = p;

 scanp:                                      //Scan for a primary element
   do { pprev = p;                           //greater than the current
        if ((p=P[p])==0) goto endp; }        //secondary element, mend the
   while(order(p,q)<=0);                     //primary list, and repeat.
   P[pprev] = q; goto scans;

 endp: p = q;                                //Attach any remaining elements,
 ends: P[pprev] = p;                         //which need not be scanned and
   return pbeg;                              //return the merged list.
 }


/*
CLARENCE LEHMAN, MARCH 1978.

Modifications:

1. Converted to C, February 1988 [CLL].

2. Internal functions marked static, April 1988 [CLL].

3. Function 'merge' made externally callable, May 1988 [CLL].

4. Pointers converted to array indices for simplicity/comprehensibility
   in scientific programs, August 2009 [CLL].

5. Converted from K&R C to Ansi C format, August 2009 [CLL].

6. Separate arrays of forward pointers, February 2011 [CLL].

NOTE:  I created this particular algorithm for the early microcomputers and have
used it ever since. Its present form was my very first C program.  I would code
it a little differently now, but that recoding would not affect its speed.
*/

