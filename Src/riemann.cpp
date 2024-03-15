#include "plutino.hpp"


/* ********************************************************************* */
double riemannLF(double *vl, double *vr, double *f, int direction)
/*  f_index e' l indice della riga che viene passato
 *
 *********************************************************************** */
{
  double ul[NVAR], ur[NVAR];
  double sourcel[NVAR], sourcer[NVAR];
  double fl[NVAR], fr[NVAR];
  double lambda_max;

  //calcolo ul che mi servira' per aggiornare i flussi
  // ognu ul ha 4 componenti le calcolo fuori dato che sono tutte diverse
  prim_to_cons(vl,ul);
  prim_to_cons(vr,ur);

  lambda_max = MaxSignalSpeedlr(vl,vr);
  flux_single(vl, ul, fl);
  flux_single(vr, ur, fr);

  for(int nv = 0; nv < NVAR; nv++){
    f[nv] = 0.5*(fl[nv]+fr[nv]) - 0.5*lambda_max*(ur[nv]-ul[nv]);  
    
  }
  //Show1DArray(f[i]);
  
  //printf(" %f \t", lambda_max);

  return(lambda_max);
}

int riemannHLLC(double **vl, double **vr, double **f, int ibeg, int iend, int direction)
{ 
  //int    VXn = VX1 + direction;
  //int    MXn = MX1 + direction;

  double ul[NVAR], ur[NVAR];
  double uls[NVAR], urs[NVAR];
  double fl[NVAR], fr[NVAR];
  double sl, sr, sstar;
  double cl, cr,cmax;
  double lambda_l, lambda_r, lambda_max;
  double coeff_left, coeff_right;
  double scrh;

  for(int i = ibeg ; i<= iend; i++){
    prim_to_cons (vl[i],ul);
    prim_to_cons (vr[i],ur);
    flux_single (vl[i], ul, fl);
    flux_single (vr[i], ur, fr);

    cl = sqrt(GAMMA_EOS*vl[i][PRS]/vl[i][RHO]);
    cr = sqrt(GAMMA_EOS*vr[i][PRS]/vr[i][RHO]);
    lambda_l = vl[i][VXn];
    lambda_r = vr[i][VXn];
    // https://www.reading.ac.uk/maths-and-stats/-/media/project/uor-main/schools-departments/maths/documents/ckongriemann.pdf?la=en&hash=5CA2DCA05ADB8CFA4E542D0E6EEB41D1
    sl = MIN(lambda_l - cl, lambda_r - cr); 
    sr = MAX(lambda_l + cl, lambda_r + cr); 
    
    sstar = vr[i][PRS] - vl[i][PRS] 
          + vl[i][RHO]*vl[i][VXn]*(sl - vl[i][VXn]) 
          - vr[i][RHO]*vr[i][VXn]*(sr - vr[i][VXn]);
    sstar /= vl[i][RHO]*(sl - vl[i][VXn]) - vr[i][RHO]*(sr - vr[i][VXn]);
    
//    ustar_single(ul, vl[i], uls, sl, sstar, direction);
//    ustar_single(ur, vr[i], urs, sr, sstar, direction);

    coeff_left = vl[i][RHO]*(sl - vl[i][VXn])/(sl - sstar);
    uls[RHO] = coeff_left * 1.0; 
    uls[MX1] = coeff_left * vl[i][VX1];
    uls[MX2] = coeff_left * vl[i][VX2];
    uls[MX3] = coeff_left * vl[i][VX3];
    uls[MXn] = coeff_left * sstar;
    scrh     = sstar + vl[i][PRS]/(vl[i][RHO]*(sl - vl[i][VXn]));
    uls[ENG] = coeff_left * (ul[ENG]/vl[i][RHO] + (sstar - vl[i][VXn])*scrh);

    coeff_right = vr[i][RHO]*(sr - vr[i][VXn])/(sr - sstar);
    urs[RHO] = coeff_right * 1.0; 
    urs[MX1] = coeff_right * vr[i][VX1];
    urs[MX2] = coeff_right * vr[i][VX2];
    urs[MX3] = coeff_right * vr[i][VX3];
    urs[MXn] = coeff_right * sstar;
    scrh     = sstar + vr[i][PRS]/(vr[i][RHO]*(sr - vr[i][VXn]));
    urs[ENG] = coeff_right * (ur[ENG]/vr[i][RHO] + (sstar - vr[i][VXn])*scrh);
    
    if (0.0 <= sl) {
      for(int nv = 0; nv < NVAR; nv++) f[i][nv] = fl[nv];
    }else if (sl < 0.0 && sstar >= 0.0) {
      for(int nv = 0; nv < NVAR; nv++) f[i][nv] = fl[nv] + sl * (uls[nv] - ul[nv]);
    } else if (sstar <= 0.0 && sr > 0.0) {
      for(int nv = 0; nv < NVAR; nv++) f[i][nv] = fr[nv] + sr * (urs[nv] - ur[nv]);
    } else if (sr <= 0.0) {
      for(int nv = 0; nv < NVAR; nv++) f[i][nv] = fr[nv];
    }
    //printf ("flux(%d) = ",i); Show1DArray(f[i]);
  }

  return(0);
}