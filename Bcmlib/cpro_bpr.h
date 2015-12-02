/*===========================================*/
/*                  PRO_BPR                  */
/*===========================================*/
#ifndef ___PRO_BPR___
#define ___PRO_BPR___
#define ___PRO_BPR_cpp___
#include "cprofile.h"

////////////////////////////////////////////////////
//...��������� ������ � ������� ������������� �����;
enum Num_table_BlocksProject {
		NUMTOPOTAB = 0,	//...��� �������� ������� ��������� ������;
		NUMGRIDTAB,			//...��� �������� ������� ������� ����� ������; 
		NUMBLOCKSPROJECT, //...����� ����� ������ � ������� �������� ������ �����������;
};

///////////////////////////////////////////////////////////
//...������������ ����� ������� � �������� ��������� �����;
const unsigned NUM_OF_TOPO_ELEMENTS = 6;

////////////////////////////
//...������� ������� ������;
Table * GetTopoTable(int N_group, int N_points = NUM_OF_TOPO_ELEMENTS);
Table * GetGridTable(int N_group);

//////////////////////////////////////////////
//...����� ������� � �������� ��������� �����;
int topo_num(char * element, int & N_points);

////////////////////////////////////////////////////
//...��������������� ������� ��� ����������� ������;
char * BlockTopoReading(char * pchar, void * context, int head_reading = 0);
char * BlockGridReading(char * pchar, void * context);

#endif
