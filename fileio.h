/*----------------------------------------------------------------------------*
FILE INPUT/OUTPUT HEADERS
*/

#ifndef TYPEDEF
#define TYPEDEF
typedef double dec;
#endif

#define MDIM     8                //Maximum number of dimensions in array.
#define MLAB    26                //Maximum index label (actually only 'a'-'z').
#define MCOL  1000                //Maximum number of columns in Centinel file.
#define CRANGE 100                //Maximum number of ranges per index label.

struct IO
{ dec    *data;                   //Address of data array.
  int     mm[MDIM*2];             //Structure in main memory.
  int     sm[MDIM*4];             //Structure in secondary memory.
};

int FileIO(char*, struct IO, char*); //Function prototypes.

// CLARENCE LEHMAN, MARCH 2011.

