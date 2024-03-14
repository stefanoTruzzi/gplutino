#include "plutino.hpp"


/* ********************************************************************* */
int boundary_conditions(double ***v, int ibeg, int iend, int jbeg, int jend){
//ibeg = 0;  iend = NX+2*NGHOST;
//jbeg = 0; jend = NY+2*NGHOST;
  int i,  j;
  
  mynvtxstart_("BOUNDARY",YELLOW);
  for (j = 0; j < jend + NGHOST; j++){
  for (i = 1; i <= NGHOST; i++){
    for (int nv = 0; nv < NVAR; nv++){
      #if BOUND == PERIODIC
      v[ibeg - i][j][nv] = v[ibeg - i + NX][j][nv];
      v[iend + i][j][nv] = v[iend + i - NX][j][nv];
      #endif
      #if BOUND == OUTFLOW
      v[ibeg - i][j][nv] = v[ibeg][j][nv];
      v[iend + i][j][nv] = v[iend][j][nv];
      #endif
    }
  }}  

  for(int i = 0 ; i < iend + NGHOST; i++){
  for(int j = 1; j <= NGHOST ; j++){
    for (int nv = 0; nv < NVAR; nv++){
      #if BOUND == PERIODIC
      v[i][jbeg - 1][nv] =  v[i][jbeg - j + NY][nv];
      v[i][jend + 1][nv] =  v[i][jend + j - NY][nv];
      #endif
      #if BOUND == OUTFLOW
      v[i][jbeg - j][nv] = v[i][jbeg][nv];
      v[i][jend + j][nv] = v[i][jend][nv];
      #endif
    }
  }}
  mynvtxstop_();

  return(0);
} 
/* ********************************************************************* */
