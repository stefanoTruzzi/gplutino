#include "plutino.hpp"


/* ********************************************************************* */
int Output (double *x, double *y,double ***v, int ibeg, int iend , int jbeg, int jend)
/*! 
 * Write output in ASCII format.
 *
 *********************************************************************** */
{
  int i,j;
  //create file to write
  ofstream file;
  //files counter
  static int nfile = 0;
  //base filename
  string filename="data.";
  std::ostringstream ss;
  // set the number of char to complete the filename number length (setw). 
  // set the filling number in before the filename (setfill)
  ss <<filename <<std::setw(4) << std::setfill('0')<<nfile;
  std::string s = ss.str();
  // set the file extension
  s+= ".out";
  // open the file
  file.open(s);
  std::cout << "WRITING "<<s<<" TO DISK" << endl; 
  #if PHYSICS == IDEALMHD
  for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){
      // write data in the file (like cout)
      file << scientific; 
      file << x[i] << " " << y[j] << " " 
           << v[i][j][RHO] << " " 
           << v[i][j][VX1] << " " << v[i][j][VX2] << " " << v[i][j][VX3] << " " 
           << v[i][j][PRS] << " "
           << v[i][j][BX1] << " " << v[i][j][BX2] << " " << v[i][j][BX3] << " "
           << endl;
    } 
    file << endl;
  }
  #else 
   for (i = ibeg; i <= iend; i++){
    for (j = jbeg; j <= jend; j++){
      // write data in the file (like cout)
      file << scientific; 
      file << x[i] << " " << y[j] << " " 
           << v[i][j][RHO] << " " 
           << v[i][j][VX1] << " " << v[i][j][VX2] << " " << v[i][j][VX3] << " " 
           << v[i][j][PRS] <<endl;
    } 
    file << endl;
  }
  #endif
  // close file after the write process
  file.close();
  // increase the static counter
  nfile++;
  return 0;
}
