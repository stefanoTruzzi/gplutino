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
 

  if(!v1D){
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

  double *v1d = &(v1D[0][0]);
  double *vll = &(vl[0][0]);
  double *vrr = &(vr[0][0]);
  double *fnew = &(f[0][0]);
  double *sourcenew = &(source[0][0]);
  double *lambda_matrix_new = &(lambda_matrix[0][0]);

  //double *lambda_mat = &(lambda_matrix[0][0]);

  int max_size=(dim_max+2*NGHOST)*NVAR; 
  
  //#pragma acc parallel loop private (v1d[:(dim_max+2*NGHOST)*NVAR],vll[:(dim_max+2*NGHOST)*NVAR],vrr[:(dim_max+2*NGHOST)*NVAR]) //takes the next loop (j) and divide it over SM and Threads 
  /*
  for(j = jbeg; j <= jend; j++){
    #pragma acc loop vector collapse(2)
    for(i = 0; i < iend+NGHOST; i++){
      //for(nv = 0; nv < NVAR; nv++) v1D[i][nv] = V[i][j][nv]; // 64megabytes (ABOUT)
      for(nv = 0; nv < NVAR; nv++) {
        v1d[i*NVAR+nv] = V[i][j][nv];
        if (i*NVAR+nv >= max_size) printf("bug\n"); 
       }
       // 64megabytes (ABOUT)
       // 2 instructions
       // load V[i][j]...
       // store V1D[i][nv]
      // SM 32 (core) same operations
    }}
  */

  mynvtxstart_("UPDATE",RED);
  mynvtxstart_("X_Cycle UPDATE",BLUE);
  #pragma acc parallel loop private (lambda,v1d[:(dim_max+2*NGHOST)*NVAR],vll[:(dim_max+2*NGHOST)*NVAR],vrr[:(dim_max+2*NGHOST)*NVAR],fnew[:(dim_max+2*NGHOST)*NVAR],sourcenew[:(dim_max+2*NGHOST)*NVAR]) //takes the next loop (j) and divide it over SM and Threads 
  for(j = jbeg; j <= jend; j++){
    #pragma acc loop vector collapse(2)
    for(i = 0; i < iend+NGHOST; i++){
      for(nv = 0; nv < NVAR; nv++) v1d[i*NVAR+nv] = V[i][j][nv];
      //for(nv = 0; nv < NVAR; nv++) v1D[i][nv] = V[i][j][nv]; // 64megabytes (ABOUT)
    }
  //debug_flag = simmetry1DCheck(v1D,ibeg,iend);

    //#pragma acc loop vector
    reconstruct ( v1d , vll , vrr, ibeg-1 , iend+1 ); 
    //reconstruct(v1D , ) 
  #pragma acc loop vector
  for(i=ibeg;i<= iend+1; i++){
   // lambda = riemannLF (vl[i] , vr[i] ,  f[i] , IDIR);
    //#pragma acc loop seq
    
    lambda = riemannLF (&vll[i*NVAR], &vrr[i*NVAR] , &fnew[i*NVAR], IDIR, i);
    lambda_matrix[i][j] = lambda /dx;
    //lambda_matrix_new[(i*dim_max+2*NGHOST)+j] = lambda / dx;
    //printf("lambda_max[%d][%d]  = %f: \n", i,j,lambda_matrix[i][j]);
  }
  #pragma acc wait
  // #pragma acc wait
  // riemannHLLC(vl,vr,f,ibeg,iend+1,0); 
  // debug_flag = simmetry1DCheck(f,ibeg,iend);
    // j is fixed outside cycle. dt is in RK
    #if SOURCE == POWELL
    #pragma acc loop vector
      for(int i = ibeg; i<= iend+1 ; i++)  
        // Bn[i] = 0.5 * (vl[i][BXn] + vr[i][BXn]);
        Bn[i] = 0.5 * (vll[i*NVAR+BXn] + vrr[i*NVAR+BXn]);
    
    #pragma acc loop vector
      for(int i = ibeg ; i<= iend; i++){
        divB[i] = (Bn[i+1] - Bn[i]) / dx; 
        //   powell_source(v1D[i],source[i],divB[i]);
        //#pragma acc loop seq
        powell_source(&v1d[i*NVAR],&sourcenew[i*NVAR],divB[i]);

        #ifdef DEBUG_DIVB
        divB_max = MAX(divB_max, divB[i]);
        divB_min = MIN(divB_min, divB[i]);
        #endif
      }
    #endif
  
    #pragma acc loop vector collapse(2)
    for(i = ibeg; i <= iend; i++){     
      for(nv = 0; nv < NVAR;nv++){
          // R[i][j][nv] =   -(f[i+1][nv]      -      f[i][nv])/dx;
          R[i][j][nv] =  -(fnew[(i+1)*NVAR+nv] - fnew[(i)*NVAR+nv])/dx;
          #if SOURCE == POWELL
          //R[i][j][nv] += dt*0.5*(source[j][nv]);
          //the dt is in the RK
          //    R[i][j][nv] += source[i][nv];
          R[i][j][nv] += sourcenew[(i)*NVAR+nv];
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
  #pragma acc parallel loop private (v1d[:(dim_max+2*NGHOST)*NVAR],vll[:(dim_max+2*NGHOST)*NVAR],vrr[:(dim_max+2*NGHOST)*NVAR],fnew[:(dim_max+2*NGHOST)*NVAR],sourcenew[:(dim_max+2*NGHOST)*NVAR]) //takes the next loop (j) and divide it over SM and Threads 
  for(i = ibeg; i <= iend; i++){
    
    #pragma acc loop vector collapse(2)
    for(j = 0; j < jend+NGHOST; j++){
      for(nv = 0; nv < NVAR; nv++) 
        v1d[j*NVAR + nv] = V[i][j][nv];
        // v1D[j][nv] = V[i][j][nv];
    }
    //#pragma acc loop vector
    reconstruct ( v1d , vll , vrr , jbeg-1 , jend+1 );    
      
    #pragma acc loop vector
    for(j=0;j<= jend+1; j++){
      //#pragma acc loop seq
      //lambda = riemannLF (  vl[j] , vr[j] ,  f[j] , JDIR);
      lambda = riemannLF (&vll[j*NVAR] , &vrr[j*NVAR] , &fnew[j*NVAR], JDIR,j);
      #if DIMENSION == DIMX2
        lambda_matrix[i][j] = lambda / dy;
        //lambda_matrix_new[(j*dim_max+2*NGHOST)+i] = lambda / dy;
      #else
        lambda_matrix[i][j] += lambda / dy;
        //lambda_matrix_new[(j*dim_max+2*NGHOST)+i] += lambda / dy;
      #endif
    }
    #pragma acc wait
       
//    riemannHLLC(vl, vr, f, jbeg, jend+1, 1);
//    debug_flag = simmetry1DCheck(f,jbeg,jend);
    #if SOURCE == POWELL
    
      #pragma acc loop vector
      for(int j = jbeg; j<= jend+1 ; j++){
        //Bn[j] = 0.5 * (vl[j][BXn] + vr[j][BXn]);
        Bn[j] = 0.5 * (vll[j*NVAR+BXn] + vrr[j*NVAR+BXn]);
      }
      
      #pragma acc loop vector
      for(int j = jbeg ; j<= jend; j++){
        divB[j] = (Bn[j+1] - Bn[j]) / dy; 
        //         powell_source(v1D[j],source[j],divB[j]);
        //#pragma acc loop seq
        powell_source(&v1d[j*NVAR],&sourcenew[j*NVAR],divB[j]);

        #ifdef DEBUG_DIVB
        divB_max = MAX(divB_max, divB[j]);
        divB_min = MIN(divB_min, divB[j]);
        #endif
      }
    #endif

    #pragma acc loop vector collapse(2)
    for(j = jbeg; j <= jend; j++){
      for(nv = 0; nv < NVAR; nv++){
        #if DIMENSIONS == DIMX2
        // R[i][j][nv] = -(f[j+1][nv] - f[j][nv])/dy; 
          R[i][j][nv] =  -(fnew[(j+1)*NVAR+nv] - fnew[(j)*NVAR+nv])/dy; 
        #else
        // R[i][j][nv] -=  (f[j+1][nv] - f[j][nv])/dy; 
          R[i][j][nv] -=  (fnew[(j+1)*NVAR+nv] - fnew[(j)*NVAR+nv])/dy; 
          #if SOURCE == POWELL
            //R[i][j][nv] += dt*(source[j][nv]);
            //the dt is in the RK
            // R[i][j][nv] += source[j][nv];
            R[i][j][nv] += sourcenew[(j)*NVAR+nv];
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
  for(int j = jbeg; j<=jend; j++){
    for(int i=ibeg; i<= iend; i++) 
      //lambda_max = MAX(lambda_matrix_new[(i*dim_max+2*NGHOST)+j],lambda_max);
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
