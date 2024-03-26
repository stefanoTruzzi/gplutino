struct DataInfo{
  double lambda_max;
  
  double dxyzmax;

  double divB_max;
  double divB_min;

  double dt_update;
  double dx;
  double dy;
};

struct Indices
{
  int VXn, MXn;
  int VXb, MXb;
  int VXt, MXt;

  int BXn,BXt,BXb;
  int EXn,EXt,EXb;

  void Indices::SetVectorIndices (int);


};