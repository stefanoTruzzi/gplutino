#include "plutino.hpp"


/* ********************************************************************* */
int reconstruct(double *v,double *vl, double *vr, int ibeg, int iend)
/*
 *
 *********************************************************************** */
{
  double slope_best;
  double slope_next;
  double slope_prev;



  #if ORDER == 1

  #pragma acc loop vector collapse(2)
  for(int i = ibeg; i <= iend;i++){
    for(int nv = 0; nv<NVAR;nv++){   
      //vl[i+1][nv]= v[i][nv]; old
      vl[(i+1)*NVAR+nv] = v[(i+1)*NVAR+nv];
      //vr[i][nv]  = v[i][nv]; old
      vr[(i)*NVAR+nv] = v[(i)*NVAR+nv];
    }
  }
  
  #elif ORDER == 2
  
  #pragma acc loop vector collapse(2)
  for(int i = ibeg; i <=iend;i++){
    for(int nv = 0; nv<NVAR;nv++){
      //slope_next = v[i+1][nv]     - v[i][nv];
      slope_next =  v[(i+1)*NVAR+nv]-v[(i)*NVAR+nv];
      //slope_prev = v[i][nv]     -  v[i-1][nv];
      slope_prev = v[(i)*NVAR+nv] - v[(i-1)*NVAR+nv];;
      
      #if LIMITER == LIM_MINMOD
        if( slope_next * slope_prev > 0.0){ 
          if(fabs(slope_next) < fabs(slope_prev)) slope_best = slope_next;
          else                                    slope_best = slope_prev;
        }
        else slope_best = 0.0;
      #else
        if( slope_next * slope_prev > 0.0) 
          slope_best = 2.0 * slope_next * slope_prev / (slope_next + slope_prev);
        else slope_best = 0.0;
      #endif
      // vl[i+1][nv] =      v[i][nv] + 0.5*slope_best; 
      vl[(i+1)*NVAR+nv] = v[(i)*NVAR+nv] + 0.5*slope_best;   // like pluto left  is +0.5 
      // vr[i][nv] =      v[i][nv] - 0.5*slope_best; 
      vr[(i)*NVAR+nv] = v[(i)*NVAR+nv] - 0.5*slope_best;   // like pluto left  is +0.5 
    }
  }

  /*
  #elif ORDER == 3
  for(int i = ibeg; i <=iend;i++){

    for(int nv = 0; nv<NVAR;nv++){
      slope_a = v[i+1][nv]-v[i][nv];
      slope_b = v[i][nv]-v[i-1][nv];
      
      if( slope_a * slope_b > 0.0){ 
        if(fabs(slope_a) < fabs(slope_b)){
            deltau = slope_a;
        }
        else deltau = slope_b;
      }
      else deltau = 0.0;
      //qL[i]   = q[i] + 0.5*dq;  pluto
      //qR[i-1] = q[i] - 0.5*dq;  pluto
      //vr[i][nv]   = v[i][nv] + 0.5*deltau; MINE OLD
      //vl[i+1][nv] = v[i][nv] - 0.5*deltau; MINE OLD
      vr[i][nv]   = v[i][nv] - 0.5*deltau;   // like pluto right is -0.5
      vl[i+1][nv] = v[i][nv] + 0.5*deltau;   // like pluto left  is +0.5 
    }
  }
  */
  #endif


  return 0;
}
