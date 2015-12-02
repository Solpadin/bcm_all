/*===========================================*/
/*                  PRO_LPR                  */
/*===========================================*/
#ifndef ___PRO_LPR___
#define ___PRO_LPR___
#define ___PRO_LPR_cpp___
#include "cprofile.h"

////////////////////////////////////////////////////
//...��������� ������ � ������� ������������� �����;
enum Num_table_LineProject {
		NUMANKERTAB = 0,    //...��� �������� ������� ������������� �����; 
		NUMTOWERTAB,        //...��� �������� ������� ���������� ����;
		NUMANKERHANGINGTAB, //...��� �������� ������� ������� �������� ������ � ������� �� �������� �����;
		NUMSTANDHANGINGTAB, //...��� �������� ������� ������� �������� ������ � ������� �� �������������� �����;
		NUMSECTIONTAB,      //...��� �������� ������� ����������� �����;
		NUMPROFILETAB,      //...��� �������� ������� ������� �����;
		NUMCLIMATETAB,      //...��� �������� ������� ������������� ������� �����;
		NUMCABLEMODETAB,    //...��� �������� ������� ������� ���������� ������;     
		NUMGROZOMODETAB,    //...��� �������� ������� ������� ���������� ����������;
		NUMWIRE1MODETAB,    //...��� �������� ������� ������� ���������� �������;
		NUMCABLEMONTAGETAB, //...��� �������� ��������� ������� ������;
		NUMGROZOMONTAGETAB, //...��� �������� ��������� ������� ����������;
		NUMWIRE1MONTAGETAB, //...��� �������� ��������� ������� �������;
		NUMCABLESCHEMETAB,  //...��� �������� ������� ���� ����������� ��� ������;
		NUMGROZOSCHEMETAB,  //...��� �������� ������� ���� ����������� ��� ����������;
		NUMWIRE1SCHEMETAB,  //...��� �������� ������� ���� ����������� ��� �������;
		NUMTENSIONTAB,      //...��� �������� ������� ������� ��� �������, ������ � ����������;
		NUMCABLEDAMPHERTAB, //...��� �������� ������� ������������ ���������� ��� ������;
		NUMGROZODAMPHERTAB, //...��� �������� ������� ������������ ���������� ��� ����������;
		NUMWIRE1DAMPHERTAB, //...��� �������� ������� ������������ ���������� ��� �������;
		NUMLINEPROJECT,     //...����� ����� ������ � ������� ������������� �����;
};

////////////////////////////////////////////////////////
//...������������ ����� �������, �������� � �����������;
const unsigned MAXCABLES  = 1;
const unsigned MAXGROZOS  = 1;
const unsigned MAXCIRCUIT = 6;

////////////////////////////////////////////
//...������������ ����� ��������� �� ������;
const unsigned MAXDAMPHER = 4;

////////////////////////////////
//...������� ������� � ��������;
const unsigned NUMGPS_POS = 4; //...������� ������ ������ GPS ��������� (GPS_N) � ������� �����;
const unsigned NUMCIRCUIT = 13;//...������� ����� �������� �� ����� (N_circuit -- "��������") � ������� �����;
const unsigned NUMMODE_P1 = 5; //...������� ������ ������ ������� ���������� (Mode_p1) � ������� ������� ����������;
const unsigned NUMTENSION = 9; //...������� ������ ������ ������� (T_break) � ������� ������� ����������;
const unsigned NUMEFMODUL = 16;//...������� ������ ������� (EF) � ������� ������� ����������;
const unsigned NUMSCTNAME = 3; //...������� ������������ ����������� (Name) � ������� �����������;
const unsigned NUMSCTTYPE = 11;//...������� ������� ����������� (Type) � ������� �����������;
const unsigned NUMAC70_72 = 9; //...������� ������� ��-70/72 (����������!) � ������� �������� ���� �� (_WIRE4);
const unsigned NUMSTENPOS = 5; //...����� ���������� ��� ������ ������� � ������� �������;
const unsigned NUMMAXSPAN = 5; //...������� ������������ ����� ����������� ������� � ���������� �������;

//////////////////////////////////////////////////////////////////////////////////////////////////
//...������������ ����� ������� � ���� ������ ����-���������� ������� � � �������� ��������� ����;
const unsigned NUM_OF_CABLE_BASE_ELEMENTS = 8;
const unsigned NUM_OF_TOWER_BASE_ELEMENTS = 22;

////////////////////////////
//...������� ������� ������;
Table * GetLineTable (int N_group);
Table * GetTowerTable(int N_group);
Table * GetLineHangTable(int N_group);
Table * GetLineSectTable(int N_group);
Table * GetLineProfileTable(int N_group, int N_points = 0);
Table * GetLineClimateTable(int N_group);
Table * GetLineModeTable(int N_group);
Table * GetLineTensionTable(int N_group);
Table * GetLineDampherTable(int N_group);
Table * GetLineDampherTable_Aeol(int N_group);
Table * GetMontageTable (int N_group, int N, double * temper = NULL);

///////////////////////////////////////////////////
//...������ ������ ��� ������� ���������� ��������;
Table * GetTenLimitTable(int N_group = 1, char * table_names[] = NULL);

//////////////////////////////////////////
//...������� ��� ������������� ��� ������;
Table * GetWireDataBaseTable   (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetCableDataBaseTable  (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetTowerDataBaseTable	 (int N_group, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * GetDampherDataBaseTable(int N_group);
Table * GetDampherSchemeTable	 (int N_group, int id_name = OK_STATE);

////////////////////////////////////////////////////////////////
// ...���������������� ����� ����� � ����������� �������� �����;
char * operating_num(Table * table, int i, char * buff, int num_project = NULL_STATE);
int	 id_anker_num (Table * table, int i, int id_tower = 0);

////////////////////////////////////////
// ...������������� ������ ������������;
void comment_correct(Table * table, int id_static_char = NULL_STATE);
void comment_correct(void * context);

////////////////////////////////////
// ...������������� ������ ��������;
void hunging_correct(Table * table, Table * tab_ank, int * anker_num, int id_anker = 1);

//////////////////////////////////////////////////////////////
// ...������������� �������: ����������� ������������� ������;
void table_double_correct(Table * table, int excl_pos = ERR_STATE, int id_static_char = NULL_STATE, int id_sequent = OK_STATE);
void table_double_correct(Table * table, Table * tab_ank, int id_static_char = NULL_STATE);

//////////////////////////////////////////////////////////////////
//...������������� ����� ��� ������� ���� ����������� �������;
void ResetLineForDampher(void ** anker_list);

////////////////////////////////////////////////////
//...��������������� ������� ��� ����������� ������;
char * LineReading       (char * pchar, void * context, int head_reading = NULL_STATE);
char * LineTowerReading  (char * pchar, void * context);
char * LineHangReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineSectReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineProfileReading(char * pchar, void * context, int convert = NULL_STATE);
char * LineClimateReading(char * pchar, void * context);
char * LineModeReading   (char * pchar, void * context, int convert = NULL_STATE);
char * LineMontageReading(char * pchar, void * context);
char * LineSchemeReading (char * pchar, void * context);
char * LineDampherReading(char * pchar, void * context, void * ankertab = NULL);
char * LineTensionReading(char * pchar, void * context, void * ankertab = NULL);
char * LineDampherReading_Aeol(char * pchar, void * context, void * ankertab = NULL);
char * DampherDataBaseReading (char * pchar, void * context, int head_reading = NULL_STATE);

////////////////////////////////////////////////////////////////////
//...����������� ������ CSV ������ ��������� ��������� � LPR ������;
void TowerCSVConvertToLPR		(void * context, void ** anker_list, int mufta_num = OK_STATE);
void TowerCSVConvertToLPR		(void * context, void ** anker_list, int beg_pos, int set_anker, int beg_pos_comment);
void TowerCSVSetSectionsToLPR (void * context, void ** anker_list, int id_section, int id_inverse = NULL_STATE);

void SimpleCSVConvertToLPR		(void * context, void ** anker_list, int mufta_num = OK_STATE);
void SimpleCSVSetTensionToLPR	(void * context, void ** anker_list, double t_norm = 18.);

#endif
