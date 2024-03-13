/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Perform index permutation and set domain integration indexes.

  The function SetIndexes() performs a cyclic permutation of the
  array indices corresponding to vector components (velocity,
  momentum, magnetic field, etc...).
  Indices are stored as global variables, see globals.h

  \author A. Mignone (mignone@to.infn.it)\n
  \date   July 09, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "plutino.hpp"


/* ********************************************************************* */
void SetVectorIndices (int dir)
/*!
 * Set vector indices and integration index range.
 *
 * \param [in] dir  the direction index
 *
 *********************************************************************** */
{
  if (dir == IDIR) {   /* -- Order: X-Y-Z  -- */
    VXn = MXn = VX1;
    VXt = MXt = VX2;
    VXb = MXb = VX3;
    #if (PHYSICS == IDEALMHD)
    BXn = BX1; 
    BXt = BX2; 
    BXb = BX3;
    #endif
  }else if (dir == JDIR){ /* -- Order: Y-Z-X  -- */
    VXn = MXn = VX2;
    VXt = MXt = VX3;
    VXb = MXb = VX1;

    #if (PHYSICS == IDEALMHD)
    BXn = BX2;
    BXt = BX3;
    BXb = BX1;
    #endif

  }else if (dir == KDIR){ /* -- Order: Z-X-Y   -- */

    VXn = MXn = VX3;
    VXt = MXt = VX1;
    VXb = MXb = VX2;

    #if (PHYSICS == IDEALMHD)
    BXn = BX3;
    BXt = BX1;
    BXb = BX2;
    #endif

  }
}
