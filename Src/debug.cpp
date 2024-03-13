#include "plutino.hpp"


void Show1DArray(double *q)
{
  int i;
  for (i = 0; i < NVAR; i++){
    printf ("q[%d] = %12.6e, ", i, q[i]);
  }
  printf ("\n");
}
int simmetry1DCheck(double **q, int ibeg, int iend)
{
  int i; 
  int flag = 1;
  int dim = iend+1; 

  for(i = ibeg; i <= iend/2; i++){
    // +1 if there is the flux instead 
    // int is = iend - (i-ibeg)
    int is = iend+1-(i-ibeg);
    int nv = RHO;
    if( fabs(q[i][nv] + q[is][nv]) > 1.e-12) {
        printf("Symmetry error: F(i = %d) = %f, F(is = %d) = %f DIFF = %f \n", 
        i , q[i][nv], is , q[is][nv],q[i][nv]-q[is][nv]);
        printf ("ibeg, iend = %d, %d\n",ibeg,iend);
        flag = 1; 
      }
    }

  return flag;
}

