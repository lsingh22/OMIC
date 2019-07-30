#ifndef _SOLVERS_H
#define _SOLVERS_H

#include "globals.h"

void Central_diff( void ); // Calculate numerical derivatives
void Steepest_descent( void ); // Find descent direction 
void Forward_track( void ); // Start moving in descent direction 
double CostFunction( int case_optimize );
#endif
