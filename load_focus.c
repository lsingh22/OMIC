//This module loads .focus file information

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <netcdf.h>
#include <stdlib.h>

//Inputs: path_to_focus_file

int main(int argc, char **argv) {
	char ch;
	char* fname;
	double var;
	FILE *fp;
	
	//printf("%lu\n", sizeof(argv[1]));
	fname = argv[1];
	//fname = "example/hsx.multi12_00.focus";
	fp = fopen(fname, "r");

	if (fp == NULL) {
		perror("Error opening file -- exiting. \n");
		exit(1);
	}
	
char* line = NULL;
size_t len = 0;
int count = 0;
char* data;
int* n;
int* m;
int i,j,k,Ncoil;

//Read in total coil number

getline(&line, &len, fp);
//printf("%c\n", *(line+1));
getline(&line, &len, fp);

data = strtok(line," ");
Ncoil = atoi(data); //second line is Ncoil

while(fgets(line, sizeof(line), fp) != NULL) {
	if (line[1] == '#') continue;
	else data = strtok(line," ");
	     printf("%s", data);
}



//For each coil read current, length, FS mode amps, and centroid locations

/*
while( getline(&line, &len, fp) != -1) {
       if( count == 1 ){
          data = strtok(line," ");
          N = atoi(data);
          n = (int *) malloc(N*sizeof(int));
          m = (int *) malloc(N*sizeof(int));
       }
       if( count > 3 && count <= 62 ){
          data = strtok(line," ");
          *(n+count-4) = atoi(data);
          data = strtok(NULL," ");
          *(m+count-4) = atoi(data);
          data = strtok(NULL," ");
          *(rmnc+count-4) = atof(data);
          data = strtok(NULL," ");
          data = strtok(NULL," ");
          data = strtok(NULL," ");
          *(zmns+count-4) = atof(data);
       }
       count++;
   }
   fclose(fp);

*/
return 0;
}

