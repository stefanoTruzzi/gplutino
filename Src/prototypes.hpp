/* arrays.cpp ********************************************************** */
double **Array2D (int, int);
double ***Array3D (int, int, int);
double *Array1D (int);
int mat_trans(double ***,double ***);
void Show1DArray(double *q);
/* ********************************************************************* */

/* boundary.cpp ******************************************************** */
int boundary_conditions(double ***, int, int, int, int);
/* ********************************************************************* */

/* cons_and_prim.cpp *************************************************** */
int cons_to_prim(double *, double *);

#pragma acc routine seq
int prim_to_cons(double *, double *);
/* ********************************************************************* */

/* debug.cpp *********************************************************** */
void Show1DArray(double *);
int simmetry1DCheck(double **, int , int);
/* ********************************************************************* */

/* flux.cpp ************************************************************ */
#pragma acc routine seq
int flux_single(double *, double *, double *);
int ustar_single(double *, double *, double *, double, double, int);
/* ********************************************************************* */

/* initial_condition.cpp *********************************************** */
void   Initial_Condition (double , double , double *);
/* ********************************************************************* */

/* input_output.cpp **************************************************** */
int    Output    (double *,double *, double ***,int,int, int, int);
/* ********************************************************************* */

/* reconstruct.cpp ***************************************************** */
#pragma acc routine vector
int reconstruct(double *,double *, double *, int, int);
/* ********************************************************************* */

/* riemann.cpp ********************************************************* */
#pragma acc routine seq
double riemannLF(double *, double *, double *, int, int);
int riemannHLLC(double **, double **, double **, int, int, int, double);
/* ********************************************************************* */

/* rungekutta.cpp ****************************************************** */
double rk_step (DataInfo&, double ***, double ***, int, int, int, int, double, double);
double rk2_step(DataInfo&, double ***, double ***, int, int, int, int, double, double);
double rk3_step(DataInfo&, double ***, double ***, int, int, int, int, double, double);
/* ********************************************************************* */

/* update.cpp ********************************************************** */
double update(DataInfo&, double ***, double ***, int, int, int, int, double, double);
/* ********************************************************************* */

/* wave_speed.cpp ****************************************************** */
#pragma acc routine seq
double MaxSignalSpeed(double *);
#pragma acc routine seq
double MaxSignalSpeedlr(double *, double *);
/* ********************************************************************* */

/* source.cpp ********************************************************** */
#pragma acc routine seq
int powell_source(double *, double *, double);
/* ********************************************************************* */

/* set_indexes.cpp ***************************************************** */
void SetVectorIndices (int);
void SetVectorIndices2D (int);
/* ********************************************************************* */

/* nvtx.cpp ***************************************************** */
void mynvtxstart_(const char *, int);
void mynvtxstop_();
/* ********************************************************************* */
