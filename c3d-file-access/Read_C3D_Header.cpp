/**************************************************************/
/***						Read C3D Header 3.c					 		***/
/**************************************************************/

/**************************************************************/
/***	Reads in header information from C3D file passed		***/
/***	in "infile" (files should be opened before calling.	***/
/***	DEPENDENCIES: Binconvert.c										***/
/**************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "c3dutil.h"

void Read_C3D_Header(
	unsigned short	*num_markers, 
	unsigned short	*num_channels,
	unsigned short	*first_field,
	unsigned short	*last_field,
	float				*scale_factor,
	unsigned short	*start_record_num,
	unsigned short	*frames_per_field,
	float				*video_rate,
	FILE				*infile)

{
	unsigned short	key1, max_gap;
	char				DEC_float[4];
	size_t res;

// Key1, byte = 1,2; word = 1
	res = fread(&key1, sizeof(key1), 1, infile); 
 
// Number of 3D points per field, byte = 3,4; word = 2
	res = fread(num_markers, sizeof(*num_markers), 1, infile); 

// Number of analog channels per field byte = 5,6; word = 3
	res = fread(num_channels, sizeof(*num_channels), 1, infile); 

// Field number of first field of video data, byte = 7,8; word = 4
	res = fread(first_field, sizeof(*first_field), 1, infile); 

// Field number of last field of video data, byte = 9,10; word = 5	
	res = fread(last_field, sizeof(*last_field), 1, infile); 

// Maximum interpolation gap in fields, byte = 11,12; word = 6
	res = fread(&max_gap, sizeof(max_gap), 1, infile); 

// Scaling Factor, bytes = 13,14,15,16; word = 7,8
	res = fread(&DEC_float[0], sizeof(char), 1, infile); // 1st byte
	res = fread(&DEC_float[1], sizeof(char), 1, infile); // 2nd byte
	res = fread(&DEC_float[2], sizeof(char), 1, infile); // 3rd byte
	res = fread(&DEC_float[3], sizeof(char), 1, infile); // 4th byte
	*scale_factor = *(float*)DEC_float; //*scale_factor = ConvertDecToFloat(DEC_float);

// Starting record number, byte = 17,18; word = 9
	res = fread(start_record_num, sizeof(*start_record_num), 1, infile);

// Number of analog frames per video field, byte = 19,20; word = 10
	res = fread(frames_per_field, sizeof(*frames_per_field), 1, infile);

// Analog channels sampled
	if (*frames_per_field != 0) 
		*num_channels /= *frames_per_field;

// Video rate in Hz, bytes = 21,22; word = 11
	res = fread(&DEC_float[0], sizeof(DEC_float[0]), 1, infile); // 1st byte
	res = fread(&DEC_float[1], sizeof(DEC_float[0]), 1, infile); // 2nd byte
	res = fread(&DEC_float[2], sizeof(DEC_float[0]), 1, infile); // 3st byte
	res = fread(&DEC_float[3], sizeof(DEC_float[0]), 1, infile); // 4nd byte
	*video_rate = *(float*)DEC_float; //*video_rate = ConvertDecToFloat(DEC_float); 

// 370 does not use the rest of the header
// Words 12 - 148 unused; position file pointer at byte 149*2=298

	//for (i=12;i<=148;i++) 
	//	res = fread(&cdum, sizeof cdum, 1, infile);
	//res = fseek(infile,298,0);

// Key2, byte = 297,298; word = 149
	//res = fread(&key2, sizeof key2, 1, infile);

// Number of defined time events, byte = 299,300; word = 150
	//res = fread(num_time_events, sizeof *num_time_events, 1, infile);

// Skip byte 301,302; word = 151
	//res = fread(&cdum, sizeof cdum, 1, infile);
	//res = fseek(infile,302,0);

// Event times, bytes 303-374;words = 152-187: 9 events (0-8) of 4 bytes
	//for (i=0;i<9;i++) {
		//res = fread(&DEC_float[0], sizeof DEC_float[0], 1, infile); // 1st byte
		//res = fread(&DEC_float[1], sizeof DEC_float[0], 1, infile); // 2nd byte
		//res = fread(&DEC_float[2], sizeof DEC_float[0], 1, infile); // 3rd byte
		//res = fread(&DEC_float[3], sizeof DEC_float[0], 1, infile); // 4th byte
		//	event_time[i] = ConvertDecToFloat(DEC_float); 
	//}

// Event switches, bytes = 375-394; words 188-197
	//for (i=0;i<9;i++) {
	//	res = fread(&event_switch[i], sizeof event_switch[i], 1, infile);
	//}

// Event Labels, bytes = 395-466; words 198-233: 9 labels (0-8) of 4 unsigned chars
	//for (i=0;i<9;i++) {
	//	res = fread(&event_label[i], sizeof event_label[i], 1, infile);
	//}

// Bytes 234 - 255 unused
	//for (i=234;i<=255;i++) 	res = fread(&cdum, sizeof cdum, 1, infile);
}

