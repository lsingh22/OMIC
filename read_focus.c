//This file reads the FOCUS output 

#include "read_focus.h"

int Nseg;
int Ncoils;

void ReadFocus(char* output_file){
   
   int ncid, varid;
   Nseg = 120;
   
   nc_open(output_file, NC_NOWRITE, &ncid);
   nc_inq_varid(ncid, "Ncoils", &varid);


}
