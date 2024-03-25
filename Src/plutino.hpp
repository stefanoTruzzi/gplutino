#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <time.h>
#include "macro.hpp"
#include "definition.hpp"
#include "structs.hpp"
#include "prototypes.hpp"


extern int VXn, VXt, VXb;
#pragma acc declare create (VXn,VXt,VXb)
extern int MXn, MXt, MXb;
#pragma acc declare create (MXn,MXt,MXb)
extern int BXn, BXt, BXb;
#pragma acc declare create (BXn,BXt,BXb)
extern int EXn, EXt, EXb;
#pragma acc declare create (EXn,EXt,BXb)


#define GREEN (0x0000FF00)
#define BLUE (0x000000FF)
#define YELLOW (0x00FFFF00)
#define SILVER (0x00C0C0C0)
#define RED (0x00FF0000)
#define PURPLE (0x800080)
#define ORANGE (0xFFA500)
#define CYAN (0x00FFFF)
#define DARKGREEN (0x006400)
#define WHITE (0xFFFFFF)
#define RAPIDS (0x7400FF)