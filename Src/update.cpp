#include "plutino.hpp"


/* ********************************************************************* */
double update(DataInfo &datainfo, double ***V, double ***R, int ibeg, int iend, int jbeg, int jend, double dx, double dy)
/* 
 *********************************************************************** */
{ 
  int dimensions = 2;
  int j,i,nv;
  int dim_max;
  int debug_flag = 0;
  double dxyzmax,dt, lambda,lambda_max;
  static double **f, **vl,  **vr, **v1D, **lambda_matrix,**source;
  static double *Bn,*divB;

  #ifdef DEBUG_DIVB
  double divB_max = -99999999.9;
  double divB_min = 99999999.9;
  #endif
  
  dim_max = MAX(NX,NY);
  
  if(v1D==NULL){
    lambda_matrix = Array2D(NX+2*NGHOST,NY+2*NGHOST);
    f   = Array2D(dim_max+2*NGHOST,NVAR);
    vl  = Array2D(dim_max+2*NGHOST,NVAR);
    vr  = Array2D(dim_max+2*NGHOST,NVAR);
    v1D = Array2D(dim_max+2*NGHOST,NVAR);
    source = Array2D(dim_max+2*NGHOST,NVAR);
    Bn = Array1D(dim_max+2*NGHOST);
    divB = Array1D(dim_max+2*NGHOST);
  }

  #if (DIMENSIONS == DIMX1) || (DIMENSIONS == DIM_2D) 
  // Start of x-direction R cycle
  SetVectorIndices(IDIR);
  
  #pragma acc parallel loop private (v1D[:dim_max+2*NGHOST][:NVAR]) //takes the next loop (j) and divide it over SM and Threads 
  for(j = jbeg; j <= jend; j++){
    for(i = 0; i < iend+NGHOST; i++){ 
      for(nv = 0; nv < NVAR; nv++) v1D[i][nv] = V[i][j][nv]; // 64megabytes (ABOUT)
      // 2 instructions
      // load V[i][j]...
      // store V1D[i][nv]
      // SM 32 (core) same operations
    }}


  mynvtxstart_("UPDATE",RED);
  mynvtxstart_("X_Cycle UPDATE",BLUE);
  for(j = jbeg; j <= jend; j++){
    for(i = 0; i < iend+NGHOST; i++){ 
      for(nv = 0; nv < NVAR; nv++) v1D[i][nv] = V[i][j][nv]; // 64megabytes (ABOUT)
    }
  //debug_flag = simmetry1DCheck(v1D,ibeg,iend);
    reconstruct ( v1D , vl , vr , ibeg-1 , iend+1 ); 
  for(i=ibeg;i<= iend+1; i++){
    lambda = riemannLF (vl[i] , vr[i] ,  f[i] , IDIR);
    lambda_matrix[i][j] = lambda /dx;
    //printf("lambda_max[%d][%d]  = %f: \n", i,j,lambda_matrix[i][j]);
  }
    
  // riemannHLLC(vl,vr,f,ibeg,iend+1,0); 
  // debug_flag = simmetry1DCheck(f,ibeg,iend);
    // j is fixed outside cycle. dt is in RK
    #if SOURCE == POWELL
      for(int i = ibeg; i<= iend+1 ; i++){
        Bn[i] = 0.5 * (vl[i][BXn] + vr[i][BXn]);
      }
      for(int i = ibeg ; i<= iend; i++){
        divB[i] = (Bn[i+1] - Bn[i]) / dx; 
        powell_source(v1D[i],source[i],divB[i]);

        #ifdef DEBUG_DIVB
        divB_max = MAX(divB_max, divB[i]);
        divB_min = MIN(divB_min, divB[i]);
        #endif
      }
    #endif

    for(i = ibeg; i <= iend; i++){     
      for(nv = 0; nv < NVAR;nv++){
          R[i][j][nv] =  -(f[i+1][nv] - f[i][nv])/dx;
          #if SOURCE == POWELL
          //R[i][j][nv] += dt*0.5*(source[j][nv]);
          //the dt is in the RK
          R[i][j][nv] += source[i][nv];
          #endif
        }
    }
  }
  mynvtxstop_();
  

  time_t stop;

  #endif

  #ifdef DEBUG_DIVB
  printf( "x_divB max = %f --- x_divB min = %f \n",divB_max,divB_min);
  divB_max = -99999999.9;
  divB_min = 99999999.9;
  #endif
  
  #if ((DIMENSIONS == DIMX2 ) || (DIMENSIONS == DIM_2D)) 
  // End of x-direction cycle
  // Start of y-direction R cycle
  // same of x-direction but with -= instead of =
    
  SetVectorIndices(JDIR);
  
  mynvtxstart_("Y_Cycle UPDATE",CYAN);
  for(i = ibeg; i <= iend; i++){
    for(j = 0; j < jend+NGHOST; j++){
      for(nv = 0; nv < NVAR; nv++) v1D[j][nv] = V[i][j][nv];
    }
    reconstruct(v1D, vl, vr, jbeg-1, jend+1);
    for(j=0;j<= jend+1; j++){
      lambda = riemannLF (  vl[j] , vr[j] ,  f[j] , JDIR);
      #if DIMENSION == DIMX2
        lambda_matrix[i][j] = lambda / dy;
      #else
        lambda_matrix[i][j] += lambda / dy;
      #endif
    }
       
//    riemannHLLC(vl, vr, f, jbeg, jend+1, 1);
//    debug_flag = simmetry1DCheck(f,jbeg,jend);
    #if SOURCE == POWELL
      for(int j = jbeg; j<= jend+1 ; j++){
        Bn[j] = 0.5 * (vl[j][BXn] + vr[j][BXn]);
      }
      for(int j = jbeg ; j<= jend; j++){
        divB[j] = (Bn[j+1] - Bn[j]) / dy; 
        powell_source(v1D[j],source[j],divB[j]);

        #ifdef DEBUG_DIVB
        divB_max = MAX(divB_max, divB[i]);
        divB_min = MIN(divB_min, divB[i]);
        #endif
      }
    #endif

    for(j = jbeg; j <= jend; j++){
      for(nv = 0; nv < NVAR; nv++){
        #if DIMENSIONS == DIMX2
          R[i][j][nv] = -(f[j+1][nv] - f[j][nv])/dy; 
        #else
          R[i][j][nv] -= (f[j+1][nv] - f[j][nv])/dy; 
          #if SOURCE == POWELL
            //R[i][j][nv] += dt*(source[j][nv]);
            //the dt is in the RK
            R[i][j][nv] += source[j][nv];
          #endif
        #endif
      }
    }
  }
  #endif
  mynvtxstop_();
  mynvtxstop_();
  #ifdef DEBUG_DIVB
  printf( "y_divB max = %f --- y_divB min = %f \n",divB_max,divB_min);
  #endif

  #if DIMENSIONS == DIMX1
  dy = dx;
  #endif
  #if DIMENSIONS == DIMX2 
  dx = dy;
  #endif

  lambda_max = 0.0;
  for(int i = ibeg;i<=iend;i++){
    for(int j=jbeg; j<= jend; j++)
    lambda_max = MAX(lambda_matrix[i][j],lambda_max);
  }
  dxyzmax = MAX(dx,dy);
  dt = (double)CFL* DIMENSIONS/lambda_max;
  
  datainfo.lambda_max = lambda_max;
  datainfo.dxyzmax = dxyzmax;
  datainfo.dx = dx;
  datainfo.dy = dy;
  datainfo.dt_update = dt;


  return dt;
}
