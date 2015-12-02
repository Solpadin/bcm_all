/*============================================*/
/*                  PRO_TOWER                 */
/*============================================*/
#ifndef ___PRO_TOWER___
#define ___PRO_TOWER___
#define ___PRO_TOWER3D_cpp___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "cgrid.h"
#else
#include "utils.h"
#endif

////////////////////////////////////////////////////
//...����� ����������� �������� ��� ���� � �� �����;
int    GetTOWERSampleCount(void);
char * GetTOWERSampleName (int N_sm);

///////////////////////////////////////////////////////////////////////
//...������� ����������� ��������� ��� ����������� �������� ����������;
void * CreateTOWERContext(int N_sm);

/////////////////////////////////////////////////////////////////
//...�p����p�������� ������������� ������������ ��p���� ��� ����;
#ifndef ___ABRIDGE_PROFILE_MODE___
int tower3D_init (void * context, CGrid * block_nd = NULL);
#endif
enum Tower_Topo {
	_TOPO_ONE_CIRCLE_TREE = 1,  
	_TOPO_TWO_CIRCLE_TREE,  
	_TOPO_SPLIT_CIRCLE_TREE,  
	_TOPO_LEFT_GROZO_TREE,  
	_TOPO_HORIZONTAL_CROSS_TREE,  
	_TOPO_ONE_CIRCLE_RUMKA,

	_TOPO_ORDINARY_DUPLEX,
	_TOPO_ONE_CIRCUIT_DUPLEX,
	_TOPO_ONE_CIRCUIT_DUPLEX_A,
	_TOPO_ONE_CIRCUIT_DUPLEX_B, //...�������� ��������� ���� ��110-3, ��110-5, ��330-2, ��330-6, ��35-4;

	_TOPO_ORDINARY_TRIPLEX,  
	_TOPO_ONE_CIRCUIT_TRIPLEX,

	_TOPO_SIMPLE_TREE,    //...����������� ������;
	_TOPO_SIMPLE_DUPLEX,  //...����������� ������;
	_TOPO_SIMPLE_TRIPLEX, //...����������� ������;
	_NUM_TOWER_TOPOS,
};
 
/////////////////////////////////////////////////////////////////////
// ...���� ������ ����������� ���� (������� �������������� ��������);
void TowerBase(void ** tower_list, char * tower, double * par, void * context = NULL);

/////////////////////////////////////
//...p������� ������� ��p�� ��� ����;
int tower3D_solv(void * context, int & sec, int & hund, int id_solv);

/////////////////////////////////////////////////
//...���������� ��� 2D ������� (������ ��� ����);
char * GetFuncTOWERName(int num);
enum Func_Tower {
    _FTOWER1 = 0,  
    _FTOWER2,  
    _FTOWER3,  
    _FTOWER4,  
    _FTOWER5,  
    _FTOWER6,
    _NUM_TOWER_FUNCTIONS
};
#define NUM_TOWER_FUNCTIONS _NUM_TOWER_FUNCTIONS
inline int GetFuncTOWERCount  (void) { return(NUM_TOWER_FUNCTIONS);}
inline int GetFuncTOWERDefault(void) { return(_FTOWER1);}

/////////////////////////////////////////////////////
//...��������� ����������� �� ������� ����� ��� ����;
#ifndef ___ABRIDGE_PROFILE_MODE___
void solver_Tower(void * context, CGrid * nd, double * F, double * par = NULL, int id_F = GetFuncTOWERDefault());
#endif

#endif
