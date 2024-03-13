/* this header define all the initial values must be changed in future */
/* YES and NO boolean macro */
#define YES        1
#define NO         0

/* SOURCE */ 
#define NOSOURCE 0    // NO correction for source term
#define POWELL 1       // Powel Correction for gradb

/* END  SOURCE */

/* BOUNDARY */
#define PERIODIC 0 
#define OUTFLOW  1

#define NGHOST    2          /* Number of ghost zones */

/* PHYSICS CONDITIONS */
#define HD 0 
#define IDEALMHD 1

/* PROBLEM DIMENSIONS AND TEST */ 
#define DIMX1 0
#define DIMX2 1
#define DIM_2D   2

// PROBLEM TYPE
//#define SHOCKTUBE 0
//#define BLASTWAVE 1
//#define ORSZAG 2
// used in initial_condition.cpp
//#define INITCOND ORSZAG
//used in set_indexes.cpp 


/* BOUNDARY CONDITION *********************************************** */
#define PERIODIC 0 
#define OUTFLOW  1

/* ****************************************************************** */

/* INDEXES DEFINITION *********************************************** */
#define IDIR     0     /*   This sequence (0,1,2) should */
#define JDIR     1     /*   never be changed             */
#define KDIR     2     /*                                */
/* ****************************************************************** */

/* MACRO FUNCTIONS DEFINITION */
#define MAX(a,b)  ((a) > (b) ? (a):(b))
#define MIN(a,b)  ((a) < (b) ? (a):(b))

using namespace std;
