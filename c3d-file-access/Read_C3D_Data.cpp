/**************************************************************/
/***						Read C3D Data.c			 					***/
/**************************************************************/

/**************************************************************/
/***	Reads in both 3D position and analog data from C3D 	***/
/***	file passed in "infile" (files should be opened			***/
/***	before calling.)													***/
/**************************************************************/
/***	4/20/00 - Altered from previous routine Read C3D File ***/
/***	to include reading of residual and camera bits			***/
/**************************************************************/

#include <stdio.h>
#include <stdlib.h>

// All data in header is specified as "16 bit words", or the equivalent of
// one unsigned char

void Read_C3D_Data_int(
	unsigned short	num_markers,				// number of marker trajectories
	unsigned short	num_analog_channels,		// number of analog channels
	unsigned short	first_field,				// first frame to read
	unsigned short	last_field,					// last frame to read
	unsigned short	start_byte,					// starting record number
	unsigned short	analog_frames_per_field,// analog samples/ video frame
	// marker xyz positions[1..num_markers][first_field..last_field]
	short			**x,							
	short			**y,
	short			**z,
	char			**residual,						// marker residual
	char			**num_cam,						// cameras used in construction
	// analog data[1..num_analog_channels][1..(last_field-first_field)*analog_frames_per_field
	short			**analog,
	// input file to read: make sure it is opened as binary
	FILE			*infile)

{
	unsigned short	frame, marker, sample, channel;
	size_t 		res;

// Startbyte is the starting location of tvd/adc data
	start_byte = 512*(start_byte-1);

// Position file pointer to start of data record
	res = fseek(infile, start_byte, 0);

// For each frame
	for (frame = 0; frame <= last_field-first_field; frame++) {
	// Read Video data
	// For each marker
		for (marker = 0; marker < num_markers; marker++)	{
		// Read X,Y,Z positions
			res = fread(&x[marker][frame], sizeof(short), 1, infile);
			res = fread(&y[marker][frame], sizeof(short), 1, infile);
			res = fread(&z[marker][frame], sizeof(short), 1, infile);

		// Read residual value and # cameras (1 byte each)
			res = fread(&residual[marker][frame], sizeof(char), 1, infile);
			res = fread(&num_cam [marker][frame], sizeof(char), 1, infile);
		}
	// Read Analog Data
	// For each analog sample/frame
		for (sample = 0; sample < analog_frames_per_field; sample++)	{
		// For each channel of analog data
			for (channel = 0; channel < num_analog_channels; channel++)
				res = fread(&analog[channel][frame*analog_frames_per_field+sample], sizeof(short), 1, infile);
		}
	}
}

void Read_C3D_Data_float(
	unsigned short	num_markers,				// number of marker trajectories
	unsigned short	num_analog_channels,		// number of analog channels
	unsigned short	first_field,				// first frame to read
	unsigned short	last_field,					// last frame to read
	unsigned short	start_byte,					// starting record number
	unsigned short	analog_frames_per_field,// analog samples/ video frame
	// marker xyz positions[1..num_markers][first_field..last_field]
	float			**x,							
	float			**y,
	float			**z,
	char			**residual,						// marker residual
	char			**num_cam,						// cameras used in construction
	// analog data[1..num_analog_channels][1..(last_field-first_field)*analog_frames_per_field
	short			**analog,
	// input file to read: make sure it is opened as binary
	FILE			*infile)

{
	unsigned short	frame, marker, sample, channel;
	char				DEC_float[4];
	int				temp_int;
	size_t 		res;

// Startbyte is the starting location of tvd/adc data
	start_byte = 512*(start_byte-1);

// Position file pointer to start of data record
	res = fseek(infile, start_byte, 0);

// For each frame
	for (frame = 0; frame <= last_field-first_field; frame++) {
	// Read Video data
	// For each marker
		for (marker = 0; marker < num_markers; marker++)	{
		// Read X,Y,Z positions
			res = fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
			res = fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
			res = fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3st byte
			res = fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4nd byte
			x[marker][frame] = *(float*)DEC_float; //x[marker][frame] = ConvertDecToFloat(DEC_float); 

			res = fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
			res = fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
			res = fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3st byte
			res = fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4nd byte
			y[marker][frame] = *(float*)DEC_float; //y[marker][frame] = ConvertDecToFloat(DEC_float); 

			res = fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
			res = fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
			res = fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3st byte
			res = fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4nd byte
			z[marker][frame] = *(float*)DEC_float; //z[marker][frame] = ConvertDecToFloat(DEC_float); 

		// Read residual value and # cameras (1 byte each)
			res = fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
			res = fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
			res = fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3nd byte
			res = fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4nd byte

			temp_int = (int)*(float*)DEC_float;
			residual[marker][frame] = (char)(temp_int);
			num_cam [marker][frame] = (char)(temp_int>>4);
		}
	// Read Analog Data
	// For each analog sample/frame
		for (sample = 0; sample < analog_frames_per_field; sample++)	{
		// For each channel of analog data
			for (channel = 0; channel < num_analog_channels; channel++)
				res = fread(&analog[channel][frame*analog_frames_per_field+sample], sizeof(short), 1, infile);
		}
	}
}
