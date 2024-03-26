
/* can be POWELL or NOSOURCE*/
#define SOURCE POWELL
/* can be PERIODIC or OUTFLOW*/
#define BOUND  PERIODIC
/* cna be 1 or 2*/
#define ORDER     2
#define LIMITER LIM_VANLEER

#define GAMMA_EOS 5.0/3.0

/* X GRID */
#define XBEG      0.0
#define XEND      1.0
/* Y GRID */
#define YBEG      0.0
#define YEND      1.0

/* grid dimensiones */
#define NX        1024
#define NY        1024

#define CFL 0.4 //MAX 0.4 for 2D application

/* TIME */
#define TSTART    1e-4
#define TSTOP     1e-3


/* INNER OUTPUT (put 0 if you want only first and last)*/
#define NOUTPUT   0

/* PHYSICS definition. */
/* IDEALMHD or HD */
#define PHYSICS IDEALMHD

/* Dimensions definition. */
/* One can run only DIMX1 to test the problem only in X direction (1D) */
/* ONLY DIMX2 to test the software only in Y direction (1D) */
/* DIM_2D to launch the software with both X and Y */
#define DIMENSIONS DIM_2D
//#define DIMENSIONS DIMX1  

/* VARIABLES DEFINITION ********************************************* */
#if PHYSICS == HD
  #define NVAR      5          /* Number of variables   */
#elif PHYSICS == IDEALMHD
  #define NVAR      8
#endif

/* VARIABLES DEFINITION ********************************************* */
/* COMMON VARS */
#define RHO  0
#define VX1  1
#define VX2  2
#define VX3  3
#define PRS  4
#define ENG  PRS
#define MX1  VX1
#define MX2  VX2
#define MX3  VX3
/* MHD VARS */
#if PHYSICS == IDEALMHD 
  #define BX1  5
  #define BX2  6
  #define BX3  7
#endif
/* ****************************************************************** */

#ifndef _USE_NVTX
#define _USE_NVTX
#endif
