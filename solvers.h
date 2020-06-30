#ifndef _SOLVERS_H
#define _SOLVERS_H

#include "globals.h"

void Central_diff(double* dof, double fb_init);
void Steepest_descent(void);  
void Forward_track(double fb_init); 
double CostFunction(int case_optimize, double fb_init); 
double MultiFieldError(void);
double ComplexityPenalty(void);
double SurfaceArea(void); 
double SingleFieldError(void);
#endif
