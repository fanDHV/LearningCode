#include "myCGAL.h"

#include "global-settings.h"

double EPSILON = 1.0e-1;
double ALPHA_EPSILON = CGAL_PI / 180.0; // in radian
int MAX_ITERATIONS = 1000;
double VIBRATE_RATIO = 0.8; // 0.5
int VIBRATE_NEDGES = 0;
// UPDATE_METHOD updateMethod;
bool AGARWAL_ORIGINAL = true;
bool PATH_LOG = false;
bool existHorizontalPath_ = true;

// trang: for debug mode
bool isDebug_ = false;
bool isDebug_2 = false;
bool isDebug_3 = false;
bool isDebug_4 = true;
bool globalComputing_ = true;
bool globalClassicalComputing_ = true;
double SIN_THETA = std::sin(15* CGAL_PI / 180.0); //theta=15 degree ~ radian
//Le: add to calculate SP, SGP_distance
bool SGP = true;
//bool PREPROCES = true;
//bool Z_DESCENDING = true;
