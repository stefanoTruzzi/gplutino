#include "plutino.hpp"


/* ********************************************************************* */
int reconstruct(double **v,double **vl, double **vr, int ibeg, int iend)
/*
 *
 *********************************************************************** */
{
  double deltau;
  double slope_a;
  double slope_b;

  #if ORDER == 1
  for(int i = ibeg; i <= iend;i++){
    for(int nv = 0; nv<NVAR;nv++){     
      // PLUTO REF
      // qL[i]   = q[i];
      // qR[i-1] = q[i];

      vl[i+1][nv]= v[i][nv];
      vr[i][nv]  = v[i][nv];
      
      
    }
  }
  
  #elif ORDER == 2
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
      vl[i+1][nv] = v[i][nv] + 0.5*deltau;   // like pluto left  is +0.5 
      vr[i][nv] = v[i][nv] - 0.5*deltau;   // like pluto left  is +0.5 

    }
  }

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
  #endif



  return 0;
}
