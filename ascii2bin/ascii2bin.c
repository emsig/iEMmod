#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(DOUBLE)
#define REAL double
#define getparfloat getpardouble
#else
#define REAL float
#endif

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    REAL r,i;
} complex;
#endif/* complex */


int main ( int argc, char *argv[] ) {
	FILE    *fid_ascii, *fp_out;
	int	iel, elements, nwrite;
        REAL	temp1, temp2, temp3;
	complex *data;
	char 	filename[1000], filename_out[1000], tempstring[1000];

	if (argc!=4) {
		printf("This program converts an ASCII file in the format of EMmod into a binary EMmod file.\n");
		printf("Usage: ./ascii2bin file_in number_of_elements file_out\n");
		printf("where\n");
		printf("- file_in is the name of the ASCII input file.\n");
		printf("- number_of_elements is an integer specifying how many elements need to be read.\n");
		printf("  If the dataset contains for example 512 times 512 datapoints, number_of_elements\n"); 
		printf("  is 262144.\n");
		printf("- file_out is the name of the binary output file.\n");
		printf("\n");
	} else {
		// Reading name of the source ASCII file
		strcpy(filename,argv[1]);
		// Reading the amount of elements to read
		elements = strtol(argv[2],NULL,10);
		// Reading the name of the target binary file
		strcpy(filename_out,argv[3]);
		// Open input file
		fid_ascii = fopen(filename,"r");
		if (fid_ascii==NULL) {
			printf("Could not open input file: %s\n",filename);
			return 0;
		}
		// Jumping the header line
		fgets(tempstring,1000,fid_ascii);
		// Reading the data from the ASCII file
		data = (complex *)calloc(elements,sizeof(complex));
		for (iel=0; iel<elements; iel++) {
			fscanf(fid_ascii,"%lf %lf %lf %lf %lf",&temp1,&temp2,&temp3,&data[iel].r,&data[iel].i);
		}
		// Close input file
		fclose(fid_ascii);

		// Open output file
	        fp_out = fopen(filename_out,"w");
		if (fp_out==NULL) {
			printf("Could not open output file: %s\n",filename_out);
			return 0;
		}
		// Write data in binary format to output file
	        nwrite = fwrite( &data[0].r, sizeof(complex), elements, fp_out);
		// Close output file
        	fclose(fp_out);
	}
	return 0;
}

