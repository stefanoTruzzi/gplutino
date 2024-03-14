#include "plutino.hpp"


double rk_step(DataInfo &datainfo, double ***U,double ***V,int ibeg,int iend,int jbeg, int jend, double dx, double dy)
{
	static double ***R;
  double dt;
  if(R == NULL){
    R  = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR); 
  }

  boundary_conditions(V,ibeg,iend,jbeg,jend); //gli passo il 3d tutto intero e ne faccio 4 dentro 
  dt = update(datainfo, V, R, ibeg, iend, jbeg, jend, dx,dy); //gli passo il 3d
  for(int i= ibeg; i <= iend; i++){
  for(int j= jbeg; j <= jend; j++){
		for(int nvar = 0; nvar < NVAR; nvar++){
			U[i][j][nvar] = U[i][j][nvar] + dt*R[i][j][nvar]; // +dt*Rt[j][i][nvar]
    }
    cons_to_prim(U[i][j],V[i][j]); 
  }}
  return (dt);
}

double rk2_step(DataInfo &datainfo, double ***U,double ***V,int ibeg,int iend,int jbeg, int jend, double dx, double dy)
{	
	static double ***R, ***R1, ***U1, ***V1;
  double dt;

  if ( U1 == NULL){
    U1 = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR);
    R  = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR);
    R1 = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR);
	  V1 = Array3D(NX + 2*NGHOST, NY+2*NGHOST, NVAR);
    }
    
	boundary_conditions(V,ibeg,iend,jbeg,jend);
  dt = update(datainfo, V, R, ibeg, iend, jbeg, jend, dx,dy); //gli passo il 3d
	//calcolo R e poi aggiorno U
	// U1 = U+ dt*R1
  
  mynvtxstart_("RK update STAGE 1", YELLOW);
	for(int i= ibeg; i <= iend; i++){
  for(int j= jbeg; j <= jend; j++){
		for(int nvar =0; nvar < NVAR; nvar++){
			U1[i][j][nvar] = U[i][j][nvar] + dt*R[i][j][nvar]; 
    }
    cons_to_prim(U1[i][j],V1[i][j]); 
  }}
  mynvtxstop_();
  boundary_conditions(V1,ibeg,iend,jbeg,jend);
  update(datainfo, V1,R1, ibeg, iend, jbeg, jend, dx, dy);
    
  mynvtxstart_("RK update STAGE 2", ORANGE);
  for(int i= ibeg; i <= iend; i++){
  for(int j= jbeg; j <= jend; j++){
		for(int nvar =0; nvar < NVAR; nvar++){
	  U[i][j][nvar] = 0.5*(U1[i][j][nvar]+U[i][j][nvar]) + 0.5*dt*R1[i][j][nvar]; // +dt*Rt[j][i][nvar]
    //U[i][j][nvar] = U[i][j][nvar] + dt * 0.5 * ( R[i][j][nvar] + R1[i][j][nvar]); // +dt*Rt[j][i][nvar]
    }
    cons_to_prim(U[i][j],V[i][j]); 
  }}
  mynvtxstop_();

	return (dt);
}

double rk3_step(DataInfo &datainfo, double ***U,double ***V,int ibeg,int iend,int jbeg, int jend, double dx, double dy){
	
	static double ***R, ***R1, ***U1, ***V1, ***R2, ***U2;
  double dt; 

  if ( U1 == NULL){
    U1 = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
    U2 = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
    R  = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
    R1 = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
    R2 = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
	  V1 = Array3D(NX + 2*NGHOST, NY + 2*NGHOST, NVAR);
  }

	boundary_conditions(V,ibeg,iend,jbeg,jend);
  dt = update(datainfo, V, R, ibeg, iend, jbeg, jend, dx,dy); //gli passo il 3d
	//calcolo R e poi aggiorno U
	// U1 = U+ dt*R1
	for(int i= ibeg; i <= iend; i++){
  for(int j= ibeg; j <= jend; j++){
		for(int nvar =0; nvar < NVAR; nvar++){
			U1[i][j][nvar] = U[i][j][nvar] + dt*R[i][j][nvar]; // +dt*Rt[j][i][nvar]
    }
    cons_to_prim(U1[i][j],V1[i][j]); 
  }}

  boundary_conditions(V1,ibeg,iend,jbeg,jend);
  update(datainfo, V1,R1, ibeg, iend, jbeg, jend, dx, dy);
  for(int i= ibeg; i <= iend; i++){
  for(int j= ibeg; j <= jend; j++){
		for(int nvar = 0; nvar < NVAR; nvar++){
	  U1[i][j][nvar] = 0.75*U[i][j][nvar] + 0.25*U1[i][j][nvar] + 0.25*dt*R1[i][j][nvar]; // +dt*Rt[j][i][nvar]
    }
    cons_to_prim(U1[i][j],V1[i][j]); 
  }}

  boundary_conditions(V1,ibeg,iend,jbeg,jend);
  update(datainfo, V1,R2, ibeg, iend, jbeg, jend, dx, dy);
  for(int i= ibeg; i <= iend; i++){
  for(int j= ibeg; j <= jend; j++){
		for(int nvar =0; nvar < NVAR; nvar++){
	  U[i][j][nvar]=2.0*U[i][j][nvar] + 0.25*dt*R1[i][j][nvar] + (2.0/3.0)*dt*R2[i][j][nvar]; // +dt*Rt[j][i][nvar]
    }
    cons_to_prim(U[i][j],V[i][j]); 
  }}

	return (dt);
}



