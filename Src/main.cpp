#include "plutino.hpp"

/* TODO */
/* 
  1 - print TIME FOR HACKATHON
  2 - print other info for debug
  3 - set RIEMANN SOLVER MACRO
  4 - test 
*/




/* ********************************************************************* */
int main()
/* 
 *    
 *
 *********************************************************************** */


{
  /* variables initialization
     can be move to another file called init?  
     it is convinient move all in MACRO.h where are stored all the variables
  */  
  DataInfo datainfo;
  int    noutput = NOUTPUT;      /* Save noutput (+1) files */
  int    ibeg = NGHOST, jbeg = NGHOST;
  int    iend = ibeg + NX - 1, jend = jbeg + NY -1 ;
  int    i , j, nv, nstep;
  double cfl   =  double(CFL);               /* courant number initialization (must be max 0.4 for 2D problem)*/
  double xbeg  =  double(XBEG);              /* start x (horizontal) domain */
  double xend  =  double(XEND);              /* end x (horizontal) domain */
  double ybeg  =  double(YBEG);              /* start y (vertical) domain */
  double yend  =  double(YEND);              /* end   y (vertical domain) */
  double tstop =  double(TSTOP);              /* Final time     */
  double tin = double(TSTART);
  double x[NX + 2*NGHOST], dx;
  double y[NY + 2*NGHOST], dy;
  double t, dt, dtdx, a;  
  double ***U;               /* Conservation Solution array */
  double ***V;               /* Primitive solution array    */
  
  U = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR); // Conservative
  V = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR); // Primitive 
  
  mynvtxstart_("GRID_INIT",DARKGREEN);
  /* START GRID INITALIZATION ********************************************* */
  dx = (xend - xbeg)/(double)NX;
  dy = (yend - ybeg)/(double)NY;
  
  // x vector/axis initialization 
  for (int i = 0; i <= iend + NGHOST; i++){
    x[i] = xbeg + (0.5 + i - ibeg)*dx;
  }
  
  // y vector/axis initialization 
  for (int j = 0; j <= jend + NGHOST; j++){
    y[j] = ybeg + (0.5 + j - jbeg)*dy;
  }
  
  // set initial condition using inital_condition.cpp file
  for (i = 0; i < NX+2*NGHOST; i++){
  for (j = 0; j < NY+2*NGHOST; j++){
    Initial_Condition(x[i],y[j],V[i][j]);
  }}

  // Conservative initialization/calculation from Primitive
  for(i = ibeg; i <= iend; i++){
  for(j = jbeg; j <= jend; j++){
    prim_to_cons(V[i][j], U[i][j]); 
  }}
  /* END OF GRID INITIALIZATION ******************************************* */
  mynvtxstop_();

  /* starting time */
  
  t=tin;  
  nstep = 0;
  double lambda_max,lambda;  
  
  mynvtxstart_("SAVE INIT COND",WHITE);
  Output (x,y, V, ibeg, iend, jbeg, jend);
  mynvtxstop_();
  
  mynvtxstart_("Time cycle",GREEN);
  time_t start_time_cycle;
  time(&start_time_cycle);   // get current time.
  while (t <= tstop){
    // calculation of max speed (lambda max)
    // print values in file
    if ((noutput > 0) && ((int)floor(noutput*(t+dt)/tstop) == (int)ceil(noutput*(t)/tstop)))
    {
      Output (x,y, V, ibeg, iend, jbeg, jend);
    }
    /* RK ORDER (3 did not work for now)*/
    #if ORDER == 1
    dt = rk_step (datainfo, U, V, ibeg, iend, jbeg, jend, dx, dy);
    #elif ORDER == 2
    dt = rk2_step(datainfo, U, V, ibeg, iend, jbeg, jend, dx, dy);
    #elif ORDER == 3  
    /* NOT WORK FOR NOW! */
    dt = rk3_step(datainfo, U, V, ibeg, iend, jbeg, jend, dx, dy);
    #else
      #error Invalid order
    #endif
    // update step and t
    printf("Nstep: %d , time: %f , dt %f lambda_max %f, dxyzmax %f \n" ,nstep, t,datainfo.dt_update, datainfo.lambda_max, datainfo.dxyzmax);

    nstep++;
    t=t+dt;
  }
  time_t stop_time_cycle;
  time(&stop_time_cycle);   // get current time after time pass.
  double diff_t = difftime(stop_time_cycle, start_time_cycle);

  printf("Total_Time = %f seconds \n",diff_t);
  mynvtxstop_();

  mynvtxstart_("SAVE LAST COND",WHITE);
  Output (x,y, V, ibeg, iend, jbeg, jend);
  mynvtxstop_();

  return (0);
}
