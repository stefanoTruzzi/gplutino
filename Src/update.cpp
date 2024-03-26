#include "plutino.hpp"


/* ********************************************************************* */
double update(DataInfo &datainfo, double ***V, double ***R, int ibeg, int iend, int jbeg, int jend, double dx, double dy, Indices &indices)
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
  indices.SetVectorIndices(IDIR);

  double *v1d = &(v1D[0][0]);
  double *vll = &(vl[0][0]);
  double *vrr = &(vr[0][0]);
  double *fnew = &(f[0][0]);
  double *sourcenew = &(source[0][0]);
  double *lambda_matrix_new = &(lambda_matrix[0][0]);
  double *Vptr = &V[0][0][0];
  //V = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR)
  int sx = (NY+2*NGHOST) * NVAR;
  int max_size=(dim_max+2*NGHOST)*NVAR; 
  
  mynvtxstart_("UPDATE",RED);
  mynvtxstart_("X_Cycle UPDATE",BLUE);

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
 
  //#pragma acc data copyin (V[:(dim_max+2*NGHOST)*NVAR][:(dim_max+2*NGHOST)*NVAR][:NVAR]) copy(R[:(dim_max+2*NGHOST)*NVAR][:(dim_max+2*NGHOST)*NVAR][:NVAR]) copyout(lambda_matrix[:NX+2*NGHOST][:NY+2*NGHOST])
  //{
    
  #pragma acc parallel loop private (divB[:(iend+1)],Bn[:iend+2],v1d[:(dim_max+2*NGHOST)*NVAR], vll[:(dim_max+2*NGHOST)*NVAR] , vrr[:(dim_max+2*NGHOST)*NVAR], fnew[:(dim_max+2*NGHOST)*NVAR], sourcenew[:(dim_max+2*NGHOST)*NVAR]) //takes the next loop (j) and divide it over SM and Threads 
    for(j = jbeg; j <= jend; j++){

    #pragma acc loop vector collapse(2)
    for(i = 0; i < iend+NGHOST; i++){
       for(nv = 0; nv < NVAR; nv++){ 
        // v1d[i*NVAR+nv] = V[i][j][nv]; 
        // changed this to obtain reduce the GPU memory troughput. 
        // before the array was not loaded efficiently 
        // V[][][] is a *** and contains a lot of ** and * 
        // Vptr contains only 1 * !!! 
        v1d[i*NVAR+nv] = Vptr[sx*i+NVAR*j+nv];}
    }
    /*
    for(j=jbeg;j<=jend ; j++){
      for(i=0;i<=iend+NGHOST;i++){
        for(nv = 0;nv<NVAR-3;nv++){printf("%f %f\t",V[i][j][nv],v1d[i*NVAR+nv]);}
        printf("\n");
      }
    }*/

    reconstruct ( v1d , vll , vrr, ibeg-1 , iend+1 ); 
        
    #pragma acc loop vector
    for(i=ibeg;i<= iend+1; i++){
      double lambda = riemannLF (&vll[i*NVAR], &vrr[i*NVAR] , &fnew[i*NVAR], IDIR, i, indices);
      lambda_matrix[i][j] = lambda /dx;
    } 

    #if SOURCE == POWELL

    #pragma acc loop vector
    for(int i = ibeg; i<= iend+1 ; i++)  
      Bn[i] = 0.5 * (vll[i*NVAR+indices.BXn] + vrr[i*NVAR+indices.BXn]);

    #pragma acc loop vector
    for(int i = ibeg ; i<= iend; i++){
      divB[i] = (Bn[i+1] - Bn[i]) / dx; 
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
          R[i][j][nv] =  -(fnew[(i+1)*NVAR+nv] - fnew[(i)*NVAR+nv])/dx;
          #if SOURCE == POWELL
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
    
  indices.SetVectorIndices(JDIR);
  
  mynvtxstart_("Y_Cycle UPDATE",CYAN);
  #pragma acc parallel loop reduction(max:lambda_max) private (divB[:(jend+1)],Bn[:jend+2],v1d[:(dim_max+2*NGHOST)*NVAR], vll[:(dim_max+2*NGHOST)*NVAR] , vrr[:(dim_max+2*NGHOST)*NVAR], fnew[:(dim_max+2*NGHOST)*NVAR], sourcenew[:(dim_max+2*NGHOST)*NVAR])//takes the next loop (j) and divide it over SM and Threads 
  for(i = ibeg; i <= iend; i++){    
    #pragma acc loop vector collapse(2)
    for(j = 0; j < jend+NGHOST; j++){ 
      // this must be broken in chunk of 32 
      for(nv = 0; nv < NVAR; nv++) 
        v1d[j*NVAR + nv] = V[i][j][nv];
        // v1D[j][nv] = V[i][j][nv];
      }
  
    //#pragma acc parallel loop //private (vll[:(dim_max+2*NGHOST)*NVAR],vrr[:(dim_max+2*NGHOST)*NVAR]) //takes the next loop (j) and divide it over SM and Threads 
    //for(i = ibeg; i <= iend; i++){
    //#pragma acc loop vector
    reconstruct ( v1d , vll , vrr , jbeg-1 , jend+1 );    
    
    #pragma acc loop vector
    for(j=0;j<= jend+1; j++){
      double lambda = riemannLF (&vll[j*NVAR] , &vrr[j*NVAR] , &fnew[j*NVAR], JDIR,j, indices);
      #if DIMENSION == DIMX2
        lambda_matrix[i][j] = lambda / dy;
      #else
        //lambda_matrix[i][j] += lambda / dy;
        lambda = lambda_matrix[i][j]+ lambda / dy;
        lambda_max = MAX(lambda,lambda_max);
      #endif
    }
    
    #if SOURCE == POWELL
      #pragma acc loop vector
      for(int j = jbeg; j<= jend+1 ; j++)
        Bn[j] = 0.5 * (vll[j*NVAR+indices.BXn] + vrr[j*NVAR+indices.BXn]);
      #pragma acc loop vector
      for(int j = jbeg ; j<= jend; j++){
        divB[j] = (Bn[j+1] - Bn[j]) / dy; 
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
          R[i][j][nv] =  -(fnew[(j+1)*NVAR+nv] - fnew[(j)*NVAR+nv])/dy; 
        #else
          R[i][j][nv] -=  (fnew[(j+1)*NVAR+nv] - fnew[(j)*NVAR+nv])/dy; 
          #if SOURCE == POWELL
            R[i][j][nv] += sourcenew[(j)*NVAR+nv];
          #endif
        #endif

      }

    }
  }
//}
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
  /*
  lambda_max = 0.0;
  
  mynvtxstart_("CHECK_MAX_SPEED",RAPIDS);

  #pragma acc parallel loop collapse(2) reduction(max:lambda_max)
  for(int j = jbeg; j<=jend; j++){
    for(int i=ibeg; i<= iend; i++) 
      //lambda_max = MAX(lambda_matrix_new[(i*dim_max+2*NGHOST)+j],lambda_max);
      lambda_max = MAX(lambda_matrix[i][j],lambda_max);
  }
  mynvtxstop_();
  */
  dxyzmax = MAX(dx,dy);
  dt = (double)CFL* DIMENSIONS/lambda_max;
  
  datainfo.lambda_max = lambda_max;
  datainfo.dxyzmax = dxyzmax;
  datainfo.dx = dx;
  datainfo.dy = dy;
  datainfo.dt_update = dt;

  return dt;
}
