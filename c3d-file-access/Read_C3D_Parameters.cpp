/**************************************************************/
/***						Read C3D Parameters.c					 	***/
/**************************************************************/

/**************************************************************/
/***	Reads in parameter records from C3D file passed			***/
/***	in "infile" (files should be opened before calling.	***/
/**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Bcmlib/utils.h"
//#include "c3dutil.h"

// All data in header is specified as "16 bit words", or the equivalent of
// one unsigned char

void Read_C3D_Parameters(
		unsigned char	**mlabels,
		unsigned char	**alabels,
		float				*gscale,
		float				*ascale,
		int				*zero_off,
		FILE				*infile)

{ 
	float				ConvertDecToFloat  (char bytes[4]);
	float				ConvertIntelToFloat(char bytes[4]);
	void				RPF(int *offset, char *type, char *dim, FILE *infile);
	unsigned char	nrow, ncol, i, j;
	unsigned short	start_parameters, length_parameters;
	char				cdum, dim, type, gname[25], pname[25], DEC_float[4];
	int				offset, processor_type;
	fpos_t		gbyte;
	size_t 		res;

// Start at the beginning
	res = fseek(infile, 0, SEEK_SET);

// bytes 1 is the position of the Parameter Section 
	res = fread(&cdum, sizeof(char), 1, infile);		// byte 1
	start_parameters = (unsigned short)(cdum-1)*512;

// Start at Parameter Section
	res = fseek(infile, start_parameters, SEEK_SET);

// bytes 1 and 2 are part of ID key (1 and 80)
	res = fread(&cdum, sizeof(char), 1, infile);		// byte 1
	res = fread(&cdum, sizeof(char), 1, infile);		// byte 2

// byte 3 holds # of parameter records to follow
	res = fread(&cdum, sizeof(char), 1, infile);		// byte 3
	length_parameters = (unsigned short)cdum*512;

// byte 4 is processor type, Vicon uses DEC (type 2)
	res = fread(&cdum, sizeof(char), 1, infile);		// byte 3
	processor_type = (int)(cdum-83);

	
// Because it is unknown how many groups and how many
// paremeters/group are stored, must scan for each variable 
// of interest.
// Note: Only data of interest passed back from this routine
// is read. Program can be modified to read other data values.

// 1st scan for group POINT
// Parameters stored in POINT are:
//		1. LABELS
//		2. DESCRIPTIONS
//		3. USED
//		4. UNITS
//		5. SCALE
//		6. RATE
//		7. DATA_START
//		8. FAMES
//		9. INITIAL_COMMAND
//		10.X_SCREEN
//		11.Y_SCREEN
#ifdef ___CONSTRUCTION___
	printf("searching for group POINT ... ");
	while (strncmp(gname, "POINT", 5) != 0) {
		fread(&gname, sizeof(char), 5, infile);
		fseek(infile, -4*(long)sizeof(char), SEEK_CUR);
	}
	fseek(infile, 4*sizeof(char), SEEK_CUR);

// Record this file position to go back to for each parameter in group
	fgetpos(infile, &gbyte);
	printf("found group POINT\n");

// Scan for each parameter of interest in group POINT
	while (strncmp(pname, "LABELS", 6) != 0)	{
		fread(&pname, sizeof(char), 6, infile);
		fseek(infile, -5*(long)sizeof(char), SEEK_CUR);
	}
// reposition to end of LABELS
	fseek(infile, 5*sizeof(char), SEEK_CUR);
	printf("\tfound parameter LABELS\n");

	RPF(&offset, &type, &dim, infile);	// dim should be 2 for a 2D array of labels[np][4]

	// read in array dimensions: should be 4 x np
	fread(&ncol, sizeof(char), 1, infile);
	fread(&nrow, sizeof(char), 1, infile);
	set_matrix(mlabels, nrow, ncol);

	for (i = 0; i < nrow; i++)
	for (j = 0; j < ncol; j++)
		fread(&mlabels[i][j], sizeof(char), 1, infile); 
#endif

// Scan for group ANALOG
// Parameters stored in group ANALOG are 
//		1. LABELS
//		2. DESCRIPTIONS
//		3. GEN_SCALE
//		4. SCALE
//		5. OFFSET
//		6. UNITS
//		7. USED
//		8. RATE
	printf("searching for group ANALOG ... ");
	while (strncmp(gname, "ANALOG", 6) != 0) {
		fread(&gname, sizeof(char), 6, infile);
		fseek(infile, -5*(long)sizeof(char), SEEK_CUR);
	}
	fseek(infile, 5*sizeof(char), SEEK_CUR);

// Record this file position to go back to for each parameter in group
	fgetpos(infile, &gbyte);
	printf("found group ANALOG\n");

// Scan for each parameter of interest in group ANALOG
// 1. GEN_SCALE
	while (strncmp(pname, "GEN_SCALE", 9) != 0)	{
		fread(&pname, sizeof(char), 9, infile);
		fseek(infile, -8*(long)sizeof(char), SEEK_CUR);
	}
// reposition to end of string GEN_SCALE
	fseek(infile, 8*sizeof(char), SEEK_CUR);
	printf("\tfound parameter GEN_SCALE\n");

	RPF(&offset, &type, &dim, infile);	// dim should be 0 for scalar, type=3 (real)
	fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
	fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
	fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3rd byte
	fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4th byte
	*gscale = *(float*)DEC_float; //*gscale = ConvertDecToFloat(DEC_float); 

// 2. SCALE - scale factor for each analog channel
// Don't reposition to start of group! Else will read SCALE in GEN_SCALE
//fseek(infile, gbyte, 0);
	while (strncmp(pname, "SCALE", 5) != 0) {
		fread(&pname, sizeof(char), 5, infile);
		fseek(infile, -4*(long)sizeof(char), SEEK_CUR);
	}
// reposition to end of string SCALE
	fseek(infile, 4*sizeof(char), SEEK_CUR);
	printf("\tfound parameter SCALE\n");

	RPF(&offset, &type, &dim, infile);	// dim should be 1 for 1D array, type=3 (real)
// Read # of components in array
	fread(&ncol, sizeof(char), 1, infile);
	ascale = (float *)new_struct(ncol*sizeof(float));
	for (i = 0; i < ncol; i++)	{
		fread(&DEC_float[0], sizeof DEC_float[0], 1, infile); // 1st byte
		fread(&DEC_float[1], sizeof DEC_float[0], 1, infile); // 2nd byte
		fread(&DEC_float[2], sizeof DEC_float[0], 1, infile); // 3rd byte
		fread(&DEC_float[3], sizeof DEC_float[0], 1, infile); // 4th byte
		ascale[i] = *(float*)DEC_float; //ascale[i] = ConvertDecToFloat(DEC_float); 
	}

// 3. OFFSET - zero offset for each analog channel
// Reposition to start of group
	fseek(infile, (long)gbyte, 0);
	while (strncmp(pname, "OFFSET", 6) != 0) {
		fread(&pname, sizeof(char), 6, infile);
		fseek(infile, -5*(long)sizeof(char), SEEK_CUR);
	}
// reposition to end of string OFFSET
	fseek(infile, 5*sizeof(char), SEEK_CUR);
	printf("\tfound parameter OFFSET\n");

	RPF(&offset, &type, &dim, infile);	// dim should be 1 for 1D array, type=2 (integer)
// Read # of components in array
	fread(&ncol, sizeof(char), 1, infile);
	zero_off = (int *)new_struct(ncol*sizeof(int));
	for (i = 0; i < ncol; i++) {
		fread(&DEC_float[0], sizeof DEC_float[0], 1, infile); // 1st byte
		fread(&DEC_float[1], sizeof DEC_float[0], 1, infile); // 2nd byte
		zero_off[i] = 256*DEC_float[1] + DEC_float[0]; 
	}

// 4. LABELS
// Reposition to start of group
	fseek(infile, (long)gbyte, 0);
	while (strncmp(pname, "LABELS", 6) != 0) {
		fread(&pname, sizeof(char), 6, infile);
		fseek(infile, -5*(long)sizeof(char), SEEK_CUR);
	}
// reposition to end of LABELS
	fseek(infile, 5*sizeof(char), SEEK_CUR);
	printf("\tfound parameter LABELS\n");

	RPF(&offset, &type, &dim, infile);	// dim should be 2 for a 2D array of labels[np][4]
// read in array dimensions: should be 4 x np
	fread(&ncol, sizeof(char), 1, infile);
	fread(&nrow, sizeof(char), 1, infile);
	set_matrix(alabels, nrow, ncol);

	for (i = 0; i < nrow; i++)
	for (j = 0; j < ncol; j++)
		fread(&alabels[i][j], sizeof(char), 1, infile); 

	return;
}

void RPF(int *offset, char *type, char *dim, FILE *infile)
{
	char	offset_low, offset_high;
// Read Parameter Format for the variable following the 
// parameter name
//		offset = number of bytes to start of next parameter (2 Bytes, reversed order: low, high)
//		T = parameter type: -1 = char, 1 = boolean, 2 = int, 3 = float
//		D = total size of array storage (incorrect, ignore)
//		d = number of dimensions for an array parameter (d=0 for single value)
//		dlen = length of data
	fread(&offset_low,  sizeof(char), 1, infile);	// byte 1
	fread(&offset_high, sizeof(char), 1, infile);	// byte 2
	*offset = 256*offset_high + offset_low;
	fread(type, sizeof(char), 1, infile);				// byte 3
	fread(dim,  sizeof(char), 1, infile);				// byte 4
}




