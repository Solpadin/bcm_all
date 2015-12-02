/**************************************************************/
/***							C3D Util.c				 					***/
/**************************************************************/
#ifndef ___C3D_UTIL___
#define ___C3D_UTIL___

/**************************************************************/
/***	Utilities to convert floating points to and from 		***/
/***	the DEC file format used by Vicon.							***/
/**************************************************************/
// Records are 256 (16 bit) words long.

// All floating point numbers are stored in DEC format. 
// The following C code may be used to convert to PC format and back again:

float ConvertIntelToFloat(char bytes[4])
{
    char p[4];
    p[0] = bytes[0];
    p[1] = bytes[1];
    p[2] = bytes[2];
    p[3] = bytes[3];
 
    return *(float*)p;
}

float ConvertDecToFloat(char bytes[4])
{
    char p[4];
    p[0] = bytes[2];
    p[1] = bytes[3];
    p[2] = bytes[0];
    p[3] = bytes[1];
    if (p[0] || p[1] || p[2] || p[3])
        --p[3];          // adjust exponent

    return *(float*)p;
}

#endif

