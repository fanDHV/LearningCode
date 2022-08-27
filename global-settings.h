#ifndef __GLOBAL_SETTINGS_H_
#define __GLOBAL_SETTINGS_H_

#define MM                        1

typedef enum UPDATE_METHOD {
  LINE_UPDATE,
  LINE_HALF_UPDATE,
  CIRCULAR_UPDATE
} UPDATE_METHOD;
extern double EPSILON;
extern double ALPHA_EPSILON;
extern int MAX_ITERATIONS;
extern double VIBRATE_RATIO;
extern int VIBRATE_NEDGES;
// extern UPDATE_METHOD updateMethod;
extern bool AGARWAL_ORIGINAL;
extern bool PATH_LOG;
extern bool existHorizontalPath_;
extern bool isDebug_;
extern bool isDebug_2;
extern bool isDebug_3;
extern bool isDebug_4;
extern bool globalComputing_;
extern bool globalClassicalComputing_;
extern double SIN_THETA;
//Le adds to calculate SP, SGP_distance
extern bool SGP;
// extern bool PREPROCES;
//extern bool Z_DESCENDING;

#endif
