#include "plutino.hpp"

/* ********************************************************************* */
double ***Array3D (int nx, int ny, int nz)
/*! 
 * Allocate memory for a 3-D array of any basic data type.
 *
 * \param [in] nx    number of elements in the 3rd dimension
 * \param [in] ny    number of elements in the 2nd dimension
 * \param [in] nz    number of elements in the 1st dimension
 * \param [in] dsize data-type of the array to be allocated
 * 
 * \return A pointer of type (char ***) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1][0...nz-1]
 *          
 *********************************************************************** */
{
  int i, j;
  double ***m;

  m = (double ***) malloc ( nx*sizeof (double **));
  #pragma acc enter data create (m[:nx])

  m[0] = (double **) malloc ( nx*ny* sizeof(double *));  


  m[0][0] = (double *) malloc ( nx*ny*nz* sizeof(double));

/* ---------------------------
       single subscript: i
   --------------------------- */ 
   
    for (i = 1; i < nx; i++) m[i] = m[i - 1] + ny;

/* ---------------------------
     double subscript:
     
      (i,0)  (0,j) 
      (i,j)  
   --------------------------- */
  
  for (j = 1; j < ny; j++) m[0][j] = m[0][j - 1] + nz;
  for (i = 1; i < nx; i++) m[i][0] = m[i - 1][0] + ny*nz;

  for (i = 1; i < nx; i++) {
    for (j = 1; j < ny; j++) {
      m[i][j] = m[i][j - 1] + nz;
      }
    }
  



  return m;
}

/* ********************************************************************* */
double **Array2D (int nx, int ny)
/*! 
 * Allocate memory for a 2-D array of doubles
 *
 * \param [in] nx    number of elements in the 2nd dimension
 * \param [in] ny    number of elements in the 1st dimension
 * 
 * \return A pointer of type (double **) to the allocated memory area 
 *         with index range [0...nx-1][0...ny-1]
 *          
 *********************************************************************** */
{
  int i;
  double **m;
 
  m    = (double **)malloc (nx*   sizeof(double *));
  m[0] = (double *) malloc (nx*ny*sizeof(double));
 
  for (i = 1; i < nx; i++) m[i] = m[i-1] + ny;
 
  return m;
}

int mat_trans(double ***v,double ***vt){
  //copio riga in colonna e viceversa?
  for(int i = 0;  i<NX+2*NGHOST;i++){ //ciclo righe v
    for(int j= 0; j<NY+2*NGHOST; j++ ){ //ciclo colonne v
      // es v è 3 x 4 ---> vt è 4 x 3 
      // NX v = NY vt e NY v = NX t
      for(int nv= 0; nv < NVAR; nv++){
        vt[i][j]=v[j][i];
      }
    }
  }
  return(0);
}

double *Array1D(int nrow){
  double *arr;
  arr =new double [nrow];
  return arr;
}